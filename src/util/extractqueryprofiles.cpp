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

void toBuffer(const char* pssm, const unsigned char* sequence, const unsigned char* consensus, const float* neffM, size_t seqLen, BaseMatrix& subMat, std::string& result) {
    for (size_t pos = 0; pos < seqLen; pos++) {
        for (size_t aa = 0; aa < Sequence::PROFILE_AA_SIZE; aa++) {
            result.push_back(pssm[pos * Sequence::PROFILE_AA_SIZE + aa]);
        }
        result.push_back(static_cast<unsigned char>(sequence[pos]));
        result.push_back(static_cast<unsigned char>(subMat.aa2num[static_cast<int>(consensus[pos])]));
        result.push_back(static_cast<unsigned char>(MathUtil::convertNeffToChar(neffM[pos])));
        result.push_back(static_cast<unsigned char>(0));
        result.push_back(static_cast<unsigned char>(0));
    }
}

int extractqueryprofiles(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    reader.open(DBReader<unsigned int>::NOSORT);

    const int inputDbtype = reader.getDbtype();

    unsigned int maxSeqLength = reader.getMaxSeqLen();
    // for SIMD memory alignment
    maxSeqLength = (maxSeqLength) / (VECSIZE_INT * 4) + 2;
    maxSeqLength *= (VECSIZE_INT * 4);

    int outputDbtype = Parameters::DBTYPE_HMM_PROFILE;
    if (par.pcmode == Parameters::PCMODE_CONTEXT_SPECIFIC) {
        outputDbtype = DBReader<unsigned int>::setExtendedDbtype(outputDbtype, Parameters::DBTYPE_EXTENDED_CONTEXT_PSEUDO_COUNTS);
    }
    DBWriter sequenceWriter(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, outputDbtype);
    sequenceWriter.open();

    DBWriter headerWriter(par.hdr2.c_str(), par.hdr2Index.c_str(), par.threads, false, Parameters::DBTYPE_GENERIC_DB);
    headerWriter.open();

    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, -0.2f);

    Debug::Progress progress(reader.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        Sequence seq(maxSeqLength + 1, inputDbtype, &subMat, 0, false, par.compBiasCorrection != 0);
        size_t querySize = 0;
        size_t queryFrom = 0;
        reader.decomposeDomainByAminoAcid(thread_idx, par.threads, &queryFrom, &querySize);
        if (querySize == 0) {
            queryFrom = 0;
        }

        size_t aaBufferSize = par.maxSeqLen + 3 + 1;
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
        char* data_cpy = (char*)malloc(sizeof(char) * maxSeqLength * Sequence::PROFILE_AA_SIZE);
        
        for (unsigned int i = queryFrom; i < (queryFrom + querySize); ++i){
            progress.updateProgress();

            unsigned int key = reader.getDbKey(i);
            const char* data = reader.getData(i, thread_idx);
            size_t seqLen = reader.getSeqLen(i);

            seq.mapSequence(i, key, data, seqLen);

            if (Parameters::isEqualDbtype(inputDbtype, Parameters::DBTYPE_HMM_PROFILE)) {
                // Copy the data
                memcpy(data_cpy, data, sizeof(char) * seqLen * Sequence::PROFILE_AA_SIZE);
                // Nothing to do for the forward strand
                toBuffer(data_cpy, seq.numSequence, seq.numConsensusSequence, seq.neffM, seqLen, subMat, result);
                sequenceWriter.writeData(result.c_str(), result.length(), key, thread_idx);
                result.clear();

                bufferLen = Orf::writeOrfHeader(buffer, key, static_cast<size_t>(0), seqLen - 1, 0, 0);
                headerWriter.writeData(buffer, bufferLen, key, thread_idx);

                // Do the same thing with the reversed profile
                // Reverse the numSequence, numConsensusSequence, seq.neffM, data_cpy
                std::reverse(seq.numSequence, seq.numSequence + seqLen);
                std::reverse(seq.numConsensusSequence, seq.numConsensusSequence + seqLen);
                std::reverse(seq.neffM, seq.neffM + seqLen);
                char tmpPssm[Sequence::PROFILE_AA_SIZE];
                int i_curr = 0;
                int j_curr = (seqLen - 1) * Sequence::PROFILE_AA_SIZE;
                for (size_t pos = 0; pos < seqLen/2; pos++) {
                    memcpy(&tmpPssm[0], data_cpy + i_curr, Sequence::PROFILE_AA_SIZE * sizeof(char));
                    memcpy(data_cpy + i_curr, data_cpy + j_curr, Sequence::PROFILE_AA_SIZE * sizeof(char));
                    memcpy(data_cpy + j_curr, &tmpPssm[0], Sequence::PROFILE_AA_SIZE * sizeof(char));
                    i_curr += Sequence::PROFILE_AA_SIZE;
                    j_curr -= Sequence::PROFILE_AA_SIZE;
                }
                for (size_t pos = 0; pos < seqLen; pos++) {
                    // Swap following position pairs: (1,15), (2,6), (4,12), (5,7), (8,9), (10,11)
                    std::swap(data_cpy[pos * Sequence::PROFILE_AA_SIZE + 1], data_cpy[pos * Sequence::PROFILE_AA_SIZE + 15]);
                    std::swap(data_cpy[pos * Sequence::PROFILE_AA_SIZE + 2], data_cpy[pos * Sequence::PROFILE_AA_SIZE + 6]);
                    std::swap(data_cpy[pos * Sequence::PROFILE_AA_SIZE + 4], data_cpy[pos * Sequence::PROFILE_AA_SIZE + 12]);
                    std::swap(data_cpy[pos * Sequence::PROFILE_AA_SIZE + 5], data_cpy[pos * Sequence::PROFILE_AA_SIZE + 7]);
                    std::swap(data_cpy[pos * Sequence::PROFILE_AA_SIZE + 8], data_cpy[pos * Sequence::PROFILE_AA_SIZE + 9]);
                    std::swap(data_cpy[pos * Sequence::PROFILE_AA_SIZE + 10], data_cpy[pos * Sequence::PROFILE_AA_SIZE + 11]);
                }
                toBuffer(data_cpy, seq.numSequence, seq.numConsensusSequence, seq.neffM, seqLen, subMat, result);
                sequenceWriter.writeData(result.c_str(), result.length(), key, thread_idx);
                
                bufferLen = Orf::writeOrfHeader(buffer, key, seqLen - 1, static_cast<size_t>(0), 0, 0);
                headerWriter.writeData(buffer, bufferLen, key, thread_idx);
                result.clear();
            } else {
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
                // Reverse the msaSequences[0], pssmRes.consensus, pssmRes.pssm, pssmRes.neffM
                std::reverse(msaSequences[0], msaSequences[0] + seqLen);
                std::reverse(pssmRes.consensus, pssmRes.consensus + seqLen);
                std::reverse(pssmRes.neffM, pssmRes.neffM + seqLen);
                char tmpPssm[Sequence::PROFILE_AA_SIZE];
                int i_curr = 0;
                int j_curr = (seqLen - 1) * Sequence::PROFILE_AA_SIZE;
                for (size_t pos = 0; pos < seqLen/2; pos++) {
                    memcpy(&tmpPssm[0], pssmRes.pssm + i_curr, Sequence::PROFILE_AA_SIZE * sizeof(char));
                    memcpy(pssmRes.pssm + i_curr, pssmRes.pssm + j_curr, Sequence::PROFILE_AA_SIZE * sizeof(char));
                    memcpy(pssmRes.pssm + j_curr, &tmpPssm[0], Sequence::PROFILE_AA_SIZE * sizeof(char));
                    i_curr += Sequence::PROFILE_AA_SIZE;
                    j_curr -= Sequence::PROFILE_AA_SIZE;
                }
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
                result.clear();
            }
        }
        delete[] profile;
        free(msaSequences);
        free(msaContent);
        delete[] pNullBuffer;
        free(data_cpy);
    }
    headerWriter.close(true);
    sequenceWriter.close(true);
    reader.close();


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

