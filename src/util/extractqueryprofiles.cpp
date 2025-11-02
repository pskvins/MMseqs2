#include "Debug.h"
#include "Parameters.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Matcher.h"
#include "Util.h"
#include "Orf.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "Prefiltering.h"
#include "MathUtil.h"
#include "PSSMCalculator.h"
#include "MultipleAlignment.h"

#include <unistd.h>
#include <climits>
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif

int extractqueryprofiles(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> *reader = new DBReader<unsigned int>(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader->open(DBReader<unsigned int>::NOSORT);

    unsigned int maxSeqLength = reader->getMaxSeqLen();
    // for SIMD memory alignment
    maxSeqLength = (maxSeqLength) / (VECSIZE_INT * 4) + 2;
    maxSeqLength *= (VECSIZE_INT * 4);

    int outputDbtype = Parameters::DBTYPE_HMM_PROFILE;
    DBWriter sequenceWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, outputDbtype);
    sequenceWriter.open();

    DBWriter headerWriter(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, false, Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();

    // BaseMatrix *subMat = Prefiltering::getSubstitutionMatrix(par.scoringMatrixFile, par.alphabetSize, 2.0, false, false);
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, -0.2f);

    Debug::Progress progress(reader->getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        Sequence seq(maxSeqLength, Parameters::DBTYPE_AMINO_ACIDS, &subMat, 0, false, true);
        size_t querySize = 0;
        size_t queryFrom = 0;
        reader->decomposeDomainByAminoAcid(thread_idx, par.threads, &queryFrom, &querySize);
        if (querySize == 0) {
            queryFrom = 0;
        }

        size_t aaBufferSize = par.maxSeqLen + 3 + 1;
        char* aa = NULL;
        if (par.translate == true) {
            aa = (char*)malloc(aaBufferSize * sizeof(char));
        }
        char* msaContent = (char*) mem_align(ALIGN_INT, sizeof(char) * (maxSeqLength + 1) * 1);
        char** msaSequences = (char**) mem_align(ALIGN_INT, sizeof(char*) * 1);
        msaSequences[0] = msaContent;
        float * profile = new float[(maxSeqLength + 1) * Sequence::PROFILE_AA_SIZE];

        char buffer[1024];
        std::string result;
        result.reserve((par.maxSeqLen + 1) * Sequence::PROFILE_READIN_SIZE * sizeof(char));
        size_t bufferLen;

        // Currently only think about single sequence input, so maxSetSize is 1
        PSSMCalculator calculator(
            &subMat, maxSeqLength + 1, 1, par.pcmode, par.pca, par.pcb
        );

        float *pNullBuffer = new float[maxSeqLength + 1];

        for (unsigned int i = queryFrom; i < (queryFrom + querySize); ++i){
            progress.updateProgress();

            unsigned int key = reader->getDbKey(i);
            const char* data = reader->getData(i, thread_idx);
            size_t seqLen = reader->getSeqLen(i);

            seq.mapSequence(i, key, data, seqLen);
            
            // fill the msaContent with the numSequence
            memcpy(msaContent, seq.numSequence, seqLen * sizeof(unsigned char));
            // fill up the sequence buffer for the SIMD profile calculation
            size_t rowSize = seqLen / (VECSIZE_INT*4);
            rowSize = (rowSize+1) * (VECSIZE_INT*4);
            size_t msaPos = seqLen;
            while(msaPos < rowSize) {
                msaContent[msaPos++] = MultipleAlignment::GAP;
            }

            PSSMCalculator::Profile pssmRes = calculator.computePSSMFromMSA(1, seqLen, (const char **) msaSequences, par.wg, 0.0);
            if (par.compBiasCorrection == true) {
                SubstitutionMatrix::calcGlobalAaBiasCorrection(&subMat, pssmRes.pssm, pNullBuffer,
                                                               Sequence::PROFILE_AA_SIZE,
                                                               seqLen);
            }
            pssmRes.toBuffer((const unsigned char*)msaSequences[0], seqLen, subMat, result);
            sequenceWriter.writeData(result.c_str(), result.length(), key, thread_idx);
            result.clear();

            bufferLen = Orf::writeOrfHeader(buffer, key, static_cast<size_t>(0), seqLen - 1, 0, 0);
            headerWriter.writeData(buffer, bufferLen, key, thread_idx);

            // Do the same thing with the reversed profile
            // Reverse the msaSequences[0], pssmRes.consensus, pssmRes.pssm
            std::reverse(msaSequences[0], msaSequences[0] + seqLen);
            std::reverse(pssmRes.consensus, pssmRes.consensus + seqLen);
            for (size_t pos = 0; pos < seqLen; pos++) {
                // Swap following position pairs: (1,15), (2,6), (4,12), (5,7), (8,9), (10,11)
                std::swap(pssmRes.pssm[pos * Sequence::PROFILE_AA_SIZE + 1], pssmRes.pssm[pos * Sequence::PROFILE_AA_SIZE + 15]);
                std::swap(pssmRes.pssm[pos * Sequence::PROFILE_AA_SIZE + 2], pssmRes.pssm[pos * Sequence::PROFILE_AA_SIZE + 6]);
                std::swap(pssmRes.pssm[pos * Sequence::PROFILE_AA_SIZE + 4], pssmRes.pssm[pos * Sequence::PROFILE_AA_SIZE + 12]);
                std::swap(pssmRes.pssm[pos * Sequence::PROFILE_AA_SIZE + 5], pssmRes.pssm[pos * Sequence::PROFILE_AA_SIZE + 7]);
                std::swap(pssmRes.pssm[pos * Sequence::PROFILE_AA_SIZE + 8], pssmRes.pssm[pos * Sequence::PROFILE_AA_SIZE + 9]);
                std::swap(pssmRes.pssm[pos * Sequence::PROFILE_AA_SIZE + 10], pssmRes.pssm[pos * Sequence::PROFILE_AA_SIZE + 11]);
            }
            pssmRes.toBuffer((const unsigned char*)msaSequences[0], seqLen, subMat, result);
            sequenceWriter.writeData(result.c_str(), result.length(), key, thread_idx);
            
            bufferLen = Orf::writeOrfHeader(buffer, key, seqLen - 1, static_cast<size_t>(0), 0, 0);
            headerWriter.writeData(buffer, bufferLen, key, thread_idx);
        }
        if (aa != NULL) {
            free(aa);
        }
        delete[] profile;
        free(msaSequences);
        free(msaContent);
        delete[] pNullBuffer;
    }
    headerWriter.close(true);
    sequenceWriter.close(true);
    reader->close();
    delete reader;


    // make identifiers stable
#pragma omp parallel
    {
#pragma omp single
        {
#pragma omp task
            {
                DBWriter::createRenumberedDB(par.hdr2, par.hdr2Index, "", "");
            }

#pragma omp task
            {
                DBWriter::createRenumberedDB(par.db2, par.db2Index, par.createLookup ? par.db1 : "", par.createLookup ? par.db1Index : "");
            }
        }
    }
    DBReader<unsigned int>::softlinkDb(par.db1, par.db2, DBFiles::SOURCE);

    return EXIT_SUCCESS;
}

