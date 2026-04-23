#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "SubstitutionMatrix.h"
#include "tantan.h"
#include "Masker.h"
#include "KSeqWrapper.h"
#include "FastSort.h"
#include "Orf.h"
#include "itoa.h"
#include "FileUtil.h"

#include <cmath>
#include <algorithm>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>

#ifdef OPENMP
#include <omp.h>
#endif

struct EntryMeta {
    std::string header;      // "name comment\n" (with trailing newline)
    size_t seqLen;           // sequence length (no sequence bytes stored)
    size_t sequenceOffset;   // byte offset in FASTA where sequence starts
    int newlineCount;        // newlines within sequence (for multi-line FASTA)
};

struct SortEntry {
    size_t entryIdx;            // index into entries[]
    size_t startPos;            // chunk start (0 if not split)
    size_t chunkLen;            // chunk length
    unsigned int shuffledKey;   // key after createdb shuffle
};

struct ChunkWriteInfo {
    size_t startPos;    // chunk start within FASTA sequence
    size_t chunkLen;    // chunk length
    size_t dataOffset;  // pre-computed byte offset in output data file
};

static int makepaddedseqdbFromFasta(Parameters &par) {
    // ==============================
    // Pass 1: Metadata collection (no sequence storage)
    // ==============================
    std::vector<EntryMeta> entries;
    size_t maxEntrySeqLen = 0;
    int dbType = par.dbType;
    bool canMmap = false;
    {
        size_t sampleCount = 0;
        size_t isNuclCnt = 0;
        const size_t testForNucSequence = 100;
        KSeqWrapper *kseq = KSeqFactory(par.db1.c_str());
        while (kseq->ReadEntry()) {
            const KSeqWrapper::KSeqEntry &e = kseq->entry;
            if (e.name.l == 0) {
                Debug(Debug::ERROR) << "Fasta entry " << entries.size() << " is invalid\n";
                EXIT(EXIT_FAILURE);
            }
            EntryMeta meta;
            meta.header.append(e.name.s, e.name.l);
            if (e.comment.l > 0) {
                meta.header.append(" ", 1);
                meta.header.append(e.comment.s, e.comment.l);
            }
            meta.header.push_back('\n');
            meta.seqLen = e.sequence.l;
            meta.sequenceOffset = e.sequenceOffset;
            meta.newlineCount = e.newlineCount;
            maxEntrySeqLen = std::max(maxEntrySeqLen, e.sequence.l);

            // Auto-detect nucleotide/amino acid (same heuristic as createdb)
            if (dbType == 0) {
                if (sampleCount < 10 || (sampleCount % 100) == 0) {
                    if (sampleCount < testForNucSequence) {
                        size_t cnt = 0;
                        for (size_t j = 0; j < e.sequence.l; j++) {
                            switch (toupper(e.sequence.s[j])) {
                                case 'T': case 'A': case 'G': case 'C': case 'U': case 'N':
                                    cnt++;
                                    break;
                            }
                        }
                        const float nuclDNAFraction = static_cast<float>(cnt) / static_cast<float>(e.sequence.l);
                        if (nuclDNAFraction > 0.9f) {
                            isNuclCnt++;
                        }
                    }
                    sampleCount++;
                }
            }

            entries.push_back(std::move(meta));
        }
        canMmap = (kseq->type == KSeqWrapper::KSEQ_FILE);
        delete kseq;

        if (entries.empty()) {
            Debug(Debug::ERROR) << "The input file has no entries\n";
            EXIT(EXIT_FAILURE);
        }

        if (dbType == 0) {
            if (isNuclCnt == sampleCount) {
                dbType = Parameters::DBTYPE_NUCLEOTIDES;
            } else {
                dbType = Parameters::DBTYPE_AMINO_ACIDS;
            }
        }
    }
    Debug(Debug::INFO) << "Database type: " << Parameters::getDbTypeName(dbType) << "\n";

    // createdb always writes DBTYPE_AMINO_ACIDS to the dbtype file,
    // so makepaddedseqdb always sees DBTYPE_AMINO_ACIDS as input type
    int seqType = Parameters::DBTYPE_AMINO_ACIDS;

    // ==============================
    // Compute output layout
    // ==============================

    // Shuffle (group by id % shuffleSplits, flatten in group order)
    const unsigned int shuffleSplits = par.shuffleDatabase ? 32 : 1;
    std::vector<std::vector<size_t>> splits(shuffleSplits);
    for (size_t i = 0; i < entries.size(); i++) {
        splits[i % shuffleSplits].push_back(i);
    }
    std::vector<size_t> shuffledOrder;
    shuffledOrder.reserve(entries.size());
    for (unsigned int s = 0; s < shuffleSplits; s++) {
        for (size_t j = 0; j < splits[s].size(); j++) {
            shuffledOrder.push_back(splits[s][j]);
        }
    }

    // Generate SortEntries in shuffled order
    bool doSplit = (par.maxSeqLen < 65535);
    std::vector<SortEntry> sortEntries;
    for (size_t k = 0; k < shuffledOrder.size(); k++) {
        size_t origId = shuffledOrder[k];
        size_t seqLen = entries[origId].seqLen;

        if (doSplit && seqLen > static_cast<size_t>(par.maxSeqLen)) {
            size_t effectiveStep = static_cast<size_t>(par.maxSeqLen) - par.sequenceOverlap;
            size_t splitCnt = static_cast<size_t>(
                ceilf(static_cast<float>(seqLen) / static_cast<float>(effectiveStep)));
            for (size_t split = 0; split < splitCnt; split++) {
                size_t startPos = split * effectiveStep;
                size_t chunkLen = std::min(static_cast<size_t>(par.maxSeqLen), seqLen - startPos);
                SortEntry se;
                se.entryIdx = origId;
                se.startPos = startPos;
                se.chunkLen = chunkLen;
                se.shuffledKey = static_cast<unsigned int>(k);
                sortEntries.push_back(se);
            }
        } else {
            SortEntry se;
            se.entryIdx = origId;
            se.startPos = 0;
            se.chunkLen = seqLen;
            se.shuffledKey = static_cast<unsigned int>(k);
            sortEntries.push_back(se);
        }
    }

    // Sort by (chunkLen+2 ASC, index_position DESC)
    size_t totalSize = sortEntries.size();
    std::vector<std::pair<size_t, size_t>> sortArray(totalSize);
    for (size_t i = 0; i < totalSize; i++) {
        sortArray[i] = std::make_pair(i, sortEntries[i].chunkLen + 2);
    }
    SORT_PARALLEL(sortArray.begin(), sortArray.end(),
        [](const std::pair<size_t, size_t> &a, const std::pair<size_t, size_t> &b) {
            if (a.second != b.second) return a.second < b.second;
            return a.first > b.first;
        });

    // Compute data file offsets — cumulative sum of aligned sizes
    const int ALIGN = 4;
    std::vector<size_t> dataOffsets(totalSize);
    size_t cumOffset = 0;
    for (size_t i = 0; i < totalSize; i++) {
        dataOffsets[i] = cumOffset;
        size_t chunkLen = sortEntries[sortArray[i].first].chunkLen;
        size_t encodedSize = chunkLen + (ALIGN - chunkLen % ALIGN) % ALIGN;
        cumOffset += encodedSize;
    }
    size_t totalDataSize = cumOffset;

    // Build entryChunkMap[origFastaId] = chunks to write for that entry
    std::vector<std::vector<ChunkWriteInfo>> entryChunkMap(entries.size());
    for (size_t i = 0; i < totalSize; i++) {
        size_t origIdx = sortArray[i].first;
        const SortEntry &se = sortEntries[origIdx];
        ChunkWriteInfo info;
        info.startPos = se.startPos;
        info.chunkLen = se.chunkLen;
        info.dataOffset = dataOffsets[i];
        entryChunkMap[se.entryIdx].push_back(info);
    }

    // ==============================
    // Write headers, lookup, source (no sequence data needed)
    // ==============================

    size_t maxChunkLen = doSplit ? static_cast<size_t>(par.maxSeqLen) : maxEntrySeqLen;

    DBWriter dbhw(par.hdr2.c_str(), par.hdr2Index.c_str(), 1, false, Parameters::DBTYPE_GENERIC_DB);
    dbhw.open();

    struct LookupInfo {
        unsigned int outputKey;
        std::string parsedName;
        unsigned int fileNumber;
    };
    std::vector<LookupInfo> lookupEntries(totalSize);

    char orfBuffer[1024];

    for (size_t i = 0; i < totalSize; i++) {
        size_t origIdx = sortArray[i].first;
        const SortEntry &se = sortEntries[origIdx];
        const EntryMeta &meta = entries[se.entryIdx];

        unsigned int outputKey = static_cast<unsigned int>(i);
        unsigned int dbKey = se.shuffledKey;

        if (doSplit && par.headerSplitMode == 0) {
            size_t fromPos = se.startPos;
            size_t toPos = se.startPos + se.chunkLen - 1;
            size_t bufferLen = Orf::writeOrfHeader(orfBuffer, dbKey, fromPos, toPos, 0, 0);
            dbhw.writeData(orfBuffer, bufferLen, outputKey, 0);
            lookupEntries[i].parsedName = Util::parseFastaHeader(orfBuffer);
            lookupEntries[i].fileNumber = static_cast<unsigned int>(origIdx);
        } else {
            dbhw.writeData(meta.header.c_str(), meta.header.length(), outputKey, 0);
            lookupEntries[i].parsedName = Util::parseFastaHeader(meta.header.c_str());
            lookupEntries[i].fileNumber = doSplit ? static_cast<unsigned int>(origIdx) : dbKey;
        }
        lookupEntries[i].outputKey = outputKey;
    }

    dbhw.close(true, false);

    // Write lookup file
    if (par.writeLookup == true) {
        std::string lookupFile = par.db2 + ".lookup";
        FILE *file = FileUtil::openAndDelete(lookupFile.c_str(), "w");
        std::string buffer;
        buffer.reserve(2048);
        DBReader<unsigned int>::LookupEntry entry;
        for (size_t i = 0; i < totalSize; i++) {
            entry.id = lookupEntries[i].outputKey;
            entry.entryName = lookupEntries[i].parsedName;
            entry.fileNumber = lookupEntries[i].fileNumber;
            DBReader<unsigned int>::lookupEntryToBuffer(buffer, entry);
            size_t written = fwrite(buffer.c_str(), sizeof(char), buffer.size(), file);
            if (written != buffer.size()) {
                Debug(Debug::ERROR) << "Cannot write to lookup file " << lookupFile << "\n";
                EXIT(EXIT_FAILURE);
            }
            buffer.clear();
        }
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << lookupFile << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    // Write source file
    {
        std::string sourceFile = par.db2 + ".source";
        FILE *source = FileUtil::openAndDelete(sourceFile.c_str(), "w");
        std::string sourceName = FileUtil::baseName(par.db1);
        char buffer[4096];
        size_t len = snprintf(buffer, sizeof(buffer), "0\t%s\n", sourceName.c_str());
        size_t written = fwrite(buffer, sizeof(char), len, source);
        if (written != len) {
            Debug(Debug::ERROR) << "Cannot write to source file " << sourceFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        if (fclose(source) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << sourceFile << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    // ==============================
    // Pass 2: Streaming encode with pwrite
    // ==============================

    // Open output data file and pre-allocate
    std::string dataFilePath = par.db2;
    int dataFd = ::open(dataFilePath.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (dataFd < 0) {
        Debug(Debug::ERROR) << "Cannot open data file " << dataFilePath << "\n";
        EXIT(EXIT_FAILURE);
    }
    if (totalDataSize > 0 && ftruncate(dataFd, totalDataSize) != 0) {
        Debug(Debug::ERROR) << "Cannot pre-allocate data file " << dataFilePath << "\n";
        EXIT(EXIT_FAILURE);
    }

    if (canMmap) {
        // mmap the input FASTA for parallel random access
        int mmapFd = ::open(par.db1.c_str(), O_RDONLY);
        if (mmapFd < 0) {
            Debug(Debug::ERROR) << "Cannot open FASTA file for mmap: " << par.db1 << "\n";
            EXIT(EXIT_FAILURE);
        }
        struct stat fileStat;
        if (fstat(mmapFd, &fileStat) != 0) {
            Debug(Debug::ERROR) << "Cannot stat FASTA file: " << par.db1 << "\n";
            EXIT(EXIT_FAILURE);
        }
        size_t fastaFileSize = fileStat.st_size;
        char *fastaData = (char *)mmap(NULL, fastaFileSize, PROT_READ, MAP_PRIVATE, mmapFd, 0);
        if (fastaData == MAP_FAILED) {
            Debug(Debug::ERROR) << "Cannot mmap FASTA file: " << par.db1 << "\n";
            EXIT(EXIT_FAILURE);
        }
#ifdef HAVE_POSIX_MADVISE
        posix_madvise(fastaData, fastaFileSize, POSIX_MADV_RANDOM);
#endif

        SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);

        Debug::Progress progress(entries.size());

#pragma omp parallel
{
        Masker masker(subMat);
        Sequence seq(maxChunkLen, seqType, &subMat, 0, false, false);
        std::string result;
        result.reserve(maxChunkLen + ALIGN);
        char *seqBuf = new char[maxEntrySeqLen + 1];
        size_t charSeqBufferSize = maxChunkLen + 1;
        unsigned char *charSequence = NULL;
        if (par.maskMode) {
            charSequence = (unsigned char *)malloc(charSeqBufferSize * sizeof(char));
        }

#pragma omp for schedule(dynamic, 100)
        for (size_t entryIdx = 0; entryIdx < entries.size(); entryIdx++) {
            progress.updateProgress();
            const EntryMeta &meta = entries[entryIdx];
            const std::vector<ChunkWriteInfo> &chunks = entryChunkMap[entryIdx];

            // Strip newlines from mmap'd data into contiguous buffer
            const char *raw = fastaData + meta.sequenceOffset;
            size_t rawLen = meta.seqLen + meta.newlineCount;
            size_t cleanLen = 0;
            for (size_t r = 0; r < rawLen && cleanLen < meta.seqLen; r++) {
                if (raw[r] != '\n' && raw[r] != '\r') {
                    seqBuf[cleanLen++] = raw[r];
                }
            }

            for (size_t c = 0; c < chunks.size(); c++) {
                const ChunkWriteInfo &chunk = chunks[c];
                const char *chunkData = seqBuf + chunk.startPos;
                size_t chunkLen = chunk.chunkLen;

                seq.mapSequence(0, 0, chunkData, chunkLen);

                if (charSequence != NULL) {
                    if (static_cast<size_t>(seq.L) >= charSeqBufferSize) {
                        charSeqBufferSize = static_cast<size_t>(seq.L * 1.5);
                        charSequence = (unsigned char *)realloc(charSequence, charSeqBufferSize * sizeof(char));
                    }
                    memcpy(charSequence, seq.numSequence, seq.L);
                    masker.maskSequence(seq, par.maskMode, par.maskProb, par.maskLowerCaseMode, par.maskNrepeats);
                    for (int j = 0; j < seq.L; j++) {
                        result.append(1, (seq.numSequence[j] == masker.maskLetterNum) ? charSequence[j] + 32 : charSequence[j]);
                    }
                } else {
                    for (int j = 0; j < seq.L; j++) {
                        char aa = chunkData[j];
                        result.append(1, (islower(aa)) ? seq.numSequence[j] + 32 : seq.numSequence[j]);
                    }
                }

                const size_t seqPadding = (seq.L % ALIGN == 0) ? 0 : ALIGN - seq.L % ALIGN;
                result.append(seqPadding, static_cast<char>(24));

                pwrite(dataFd, result.c_str(), result.size(), chunk.dataOffset);
                result.clear();
            }
        }

        delete[] seqBuf;
        if (charSequence != NULL) {
            free(charSequence);
        }
}

        munmap(fastaData, fastaFileSize);
        ::close(mmapFd);
    } else {
        // Fallback: sequential Pass 2 for compressed/stdin inputs
        SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
        Masker masker(subMat);
        Sequence seq(maxChunkLen, seqType, &subMat, 0, false, false);
        std::string result;
        result.reserve(maxChunkLen + ALIGN);

        size_t charSeqBufferSize = maxChunkLen + 1;
        unsigned char *charSequence = NULL;
        if (par.maskMode) {
            charSequence = (unsigned char *)malloc(charSeqBufferSize * sizeof(char));
        }

        Debug::Progress progress(entries.size());

        KSeqWrapper *kseq = KSeqFactory(par.db1.c_str());
        size_t entryIdx = 0;
        while (kseq->ReadEntry()) {
            progress.updateProgress();

            const KSeqWrapper::KSeqEntry &e = kseq->entry;
            const std::vector<ChunkWriteInfo> &chunks = entryChunkMap[entryIdx];

            for (size_t c = 0; c < chunks.size(); c++) {
                const ChunkWriteInfo &chunk = chunks[c];
                const char *chunkData = e.sequence.s + chunk.startPos;
                size_t chunkLen = chunk.chunkLen;

                seq.mapSequence(0, 0, chunkData, chunkLen);

                if (charSequence != NULL) {
                    if (static_cast<size_t>(seq.L) >= charSeqBufferSize) {
                        charSeqBufferSize = static_cast<size_t>(seq.L * 1.5);
                        charSequence = (unsigned char *)realloc(charSequence, charSeqBufferSize * sizeof(char));
                    }
                    memcpy(charSequence, seq.numSequence, seq.L);
                    masker.maskSequence(seq, par.maskMode, par.maskProb, par.maskLowerCaseMode, par.maskNrepeats);
                    for (int j = 0; j < seq.L; j++) {
                        result.append(1, (seq.numSequence[j] == masker.maskLetterNum) ? charSequence[j] + 32 : charSequence[j]);
                    }
                } else {
                    for (int j = 0; j < seq.L; j++) {
                        char aa = chunkData[j];
                        result.append(1, (islower(aa)) ? seq.numSequence[j] + 32 : seq.numSequence[j]);
                    }
                }

                const size_t seqPadding = (seq.L % ALIGN == 0) ? 0 : ALIGN - seq.L % ALIGN;
                result.append(seqPadding, static_cast<char>(24));

                ssize_t written = pwrite(dataFd, result.c_str(), result.size(), chunk.dataOffset);
                if (written != static_cast<ssize_t>(result.size())) {
                    Debug(Debug::ERROR) << "Cannot write to data file\n";
                    EXIT(EXIT_FAILURE);
                }

                result.clear();
            }

            entryIdx++;
        }
        delete kseq;

        if (charSequence != NULL) {
            free(charSequence);
        }
    }

    ::close(dataFd);

    // Write index file
    {
        FILE *idxFile = FileUtil::openAndDelete(par.db2Index.c_str(), "w");
        char buffer[1024];
        for (size_t i = 0; i < totalSize; i++) {
            size_t chunkLen = sortEntries[sortArray[i].first].chunkLen;
            size_t len = DBWriter::indexToBuffer(buffer, static_cast<unsigned int>(i), dataOffsets[i], chunkLen + 2);
            size_t written = fwrite(buffer, sizeof(char), len, idxFile);
            if (written != len) {
                Debug(Debug::ERROR) << "Cannot write to index file\n";
                EXIT(EXIT_FAILURE);
            }
        }
        if (fclose(idxFile) != 0) {
            Debug(Debug::ERROR) << "Cannot close index file\n";
            EXIT(EXIT_FAILURE);
        }
    }

    // Write dbtype file
    {
        int paddedDbType = DBReader<unsigned int>::setExtendedDbtype(seqType, Parameters::DBTYPE_EXTENDED_GPU);
        DBWriter::writeDbtypeFile(par.db2.c_str(), paddedDbType, false);
    }

    return EXIT_SUCCESS;
}

int makepaddedseqdb(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    // Match splitsequence defaults so fused FASTA path behaves identically
    par.maxSeqLen = 10000;
    par.parseParameters(argc, argv, command, true, 0, 0);

    // Detect FASTA input: if no .dbtype file exists, run fused FASTA codepath
    bool isFastaInput = !FileUtil::fileExists(par.db1dbtype.c_str());
    if (isFastaInput) {
        return makepaddedseqdbFromFasta(par);
    }

    // Existing DB-input codepath (unchanged)
    const int mode = DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA;
    DBReader<unsigned int> dbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, mode);
    dbr.open(DBReader<unsigned int>::SORT_BY_LENGTH);

    DBReader<unsigned int> dbhr(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads, mode);
    dbhr.open(DBReader<unsigned int>::NOSORT);

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);

    int dbType = DBReader<unsigned int>::setExtendedDbtype(dbr.getDbtype(), Parameters::DBTYPE_EXTENDED_GPU);
    DBWriter dbsw(par.db2.c_str(), par.db2Index.c_str(), par.threads, false, dbType);
    dbsw.open();
    DBWriter dbhw(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, false, Parameters::DBTYPE_GENERIC_DB);
    dbhw.open();

    // need to prune low scoring k-mers through masking

    Debug::Progress progress(dbr.getSize());
#pragma omp parallel
{
    unsigned int thread_idx = 0;
#ifdef OPENMP
    thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
    Masker masker(subMat);
    std::string result;
    result.reserve(par.maxSeqLen);

    const int ALIGN = 4;
    Sequence seq(dbr.getMaxSeqLen(), dbr.getDbtype(), &subMat,  0, false, false);

    size_t firstIt = SIZE_MAX;
    unsigned int seqKey = 0;

    size_t charSeqBufferSize = par.maxSeqLen + 1;
    unsigned char *charSequence = NULL;
    if (par.maskMode) {
        charSequence = (unsigned char*)malloc(charSeqBufferSize * sizeof(char));
    }

#pragma omp for schedule(static)
    for (size_t i = 0; i < dbr.getSize(); i++) {
        progress.updateProgress();

        if (firstIt == SIZE_MAX) {
            firstIt = i;
        }

        size_t id = dbr.getSize() - 1 - i;
        unsigned int key = dbr.getDbKey(id);
        char *data = dbr.getData(id, thread_idx);
        size_t seqLen = dbr.getSeqLen(id);
        seq.mapSequence(id, key, data, seqLen);

        if (charSequence != NULL) {
            if ((size_t)seq.L >= charSeqBufferSize) {
                charSeqBufferSize = seq.L * 1.5;
                charSequence = (unsigned char*)realloc(charSequence, charSeqBufferSize * sizeof(char));
            }
            memcpy(charSequence, seq.numSequence, seq.L);
            masker.maskSequence(seq, par.maskMode, par.maskProb, par.maskLowerCaseMode, par.maskNrepeats);
            for (int i = 0; i < seq.L; i++) {
                result.append(1, (seq.numSequence[i] == masker.maskLetterNum) ? charSequence[i] + 32 : charSequence[i]);
            }
        } else {
            for (int i = 0; i < seq.L; i++) {
                char aa = data[i];
                result.append(1, (islower(aa)) ? seq.numSequence[i] + 32 : seq.numSequence[i]);
            }
        }
        const size_t sequencepadding = (seq.L % ALIGN == 0) ? 0 : ALIGN - seq.L % ALIGN;
        result.append(sequencepadding, static_cast<char>(24)); // Now it is 24 in dinucleotide
        dbsw.writeData(result.c_str(), result.size(), key, thread_idx, false, false);

        // + 2 is needed for newline and null character
        size_t start = dbsw.getStart(thread_idx);
        if (start % 4 != 0) {
            Debug(Debug::ERROR) << "Misalligned entry\n";
            EXIT(EXIT_FAILURE);
        }
        dbsw.writeIndexEntry(firstIt + seqKey, start, seq.L + 2, thread_idx);

        unsigned int headerId = dbhr.getId(key);
        dbhw.writeData(dbhr.getData(headerId, thread_idx), dbhr.getEntryLen(headerId), firstIt + seqKey, thread_idx, false);

        seqKey++;
        result.clear();
    }
    if (charSequence != NULL) {
        free(charSequence);
    }
}
    dbsw.close(true, false);
    dbhw.close(true, false);
    dbhr.close();
    if (par.writeLookup == true) {
        DBReader<unsigned int> readerHeader(par.hdr2.c_str(), par.hdr2Index.c_str(), 1, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        readerHeader.open(DBReader<unsigned int>::NOSORT);
        // create lookup file
        std::string lookupFile = par.db2 + ".lookup";
        FILE* file = FileUtil::openAndDelete(lookupFile.c_str(), "w");
        std::string buffer;
        buffer.reserve(2048);
        DBReader<unsigned int>::LookupEntry entry;
        size_t totalSize = dbr.getSize();
        for (unsigned int id = 0; id < readerHeader.getSize(); id++) {
            char *header = readerHeader.getData(id, 0);
            entry.id = id;
            entry.entryName = Util::parseFastaHeader(header);
            entry.fileNumber = dbr.getDbKey(totalSize - 1 - id);
            readerHeader.lookupEntryToBuffer(buffer, entry);
            int written = fwrite(buffer.c_str(), sizeof(char), buffer.size(), file);
            if (written != (int)buffer.size()) {
                Debug(Debug::ERROR) << "Cannot write to lookup file " << lookupFile << "\n";
                EXIT(EXIT_FAILURE);
            }
            buffer.clear();
        }
        if (fclose(file) != 0) {
            Debug(Debug::ERROR) << "Cannot close file " << lookupFile << "\n";
            EXIT(EXIT_FAILURE);
        }
        readerHeader.close();
    }
    dbr.close();
    return EXIT_SUCCESS;
}
