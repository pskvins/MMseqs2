#include <climits>
#include "BaseMatrix.h"

#include "simd.h"
#include "Debug.h"
#include "Util.h"
#include "MathUtil.h"
#include "Sequence.h"
#include <vector>

const double BaseMatrix::ANY_BACK = 1E-5;

BaseMatrix::BaseMatrix(){
    // init [amino acid <-> int] mappings

    num2aa = new char[255];
    aa2num = new unsigned char[USHRT_MAX];
    for (int i = 0; i < USHRT_MAX; ++i) {
        aa2num[i] = UCHAR_MAX;
    }
    // Reverse complement table for dinucleotides
    revcomp = new unsigned char[21];
    for (int i = 0; i < 21; ++i) {
        revcomp[i] = static_cast<unsigned char>('X');
    }
    // To handle the first and the last incomplete dinucleotide
    head = new unsigned char[21];
    for (int i = 0; i < 21; ++i) {
        head[i] = static_cast<unsigned char>('X');
    }
    tail = new unsigned char[21];
    for (int i = 0; i < 21; ++i) {
        tail[i] = static_cast<unsigned char>('X');
    }
    // Mapping from dinucleotide to the first or the second nucleotide
    dinucToNuc = new unsigned char[21];
}

BaseMatrix::~BaseMatrix(){
    delete[] num2aa;
    delete[] aa2num;
    delete[] revcomp;
    delete[] head;
    delete[] tail;
    delete[] dinucToNuc;
    delete[] pBack;
    for (int i = 0; i < allocatedAlphabetSize; i++){
        delete[] probMatrix[i];
        delete[] subMatrix[i];
        free(subMatrixPseudoCounts[i]);
    }
    delete[] probMatrix;
    delete[] subMatrixPseudoCounts;
    delete[] subMatrix;
}


void BaseMatrix::initMatrixMemory(int alphabetSize) {
    allocatedAlphabetSize = alphabetSize;

    // init the background probabilities, joint probability and scoring matrices with zeros
    pBack = new double[alphabetSize];
    probMatrix = new double*[alphabetSize];
    subMatrix = new short*[alphabetSize];
    subMatrixPseudoCounts = new float*[alphabetSize];
    // temporary memory setter

    for (int i = 0; i < alphabetSize; i++){
        pBack[i] = 0.0;
        probMatrix[i] = new double[alphabetSize];
        subMatrix[i] = new short[alphabetSize];
        subMatrixPseudoCounts[i] =  (float *) malloc_simd_float(alphabetSize * sizeof(float));
        for (int j = 0; j < alphabetSize; j++){
            probMatrix[i][j] = 0.0;
            subMatrix[i][j] = 0;
            subMatrixPseudoCounts[i][j] = 0.0;
        }
    }
}


void BaseMatrix::print(short** matrix, char* num2aa, int size){
    std::cout << "\n";
    short avg = 0;
    printf("     ");
    for (int i = 0; i < size; i++)
        printf("%4c ", num2aa[i]);
    std::cout << "\n";
    for (int i = 0; i < size; i++){
        printf("%4c ", num2aa[i]);
        for (int j = 0; j < size; j++){
            printf("%4d ", matrix[i][j]);
            avg += matrix[i][j];
        }
        std::cout << "\n";
    }
    std::cout << ((double)avg/(double)(size*size)) << "\n";
}

void BaseMatrix::print(double** matrix, char* num2aa, int size){
    std::cout << "\n";
    double avg = 0.0;
    printf("%7c ", ' ');
    for (int i = 0; i < size; i++)
        printf("%7c ", num2aa[i]);
    std::cout << "\n";
    for (int i = 0; i < size; i++){
        printf("%7c ", num2aa[i]);
        for (int j = 0; j < size; j++){
            printf("%7.4f ", matrix[i][j]);
            avg += matrix[i][j];
        }
        std::cout << "\n";
    }
    std::cout << (avg/(double)(size*size)) << "\n";
}

void BaseMatrix::computeBackground(double ** probMat, double * pBack, int alphabetSize, bool containsX){
    for(int i = 0; i < alphabetSize; i++){
        pBack[i] = 0;
        for(int j = 0; j < alphabetSize; j++){
            pBack[i]+=probMat[i][j];
        }
    }
    if(containsX){
        pBack[alphabetSize-1] = ANY_BACK;
    }

}

std::vector<size_t> BaseMatrix::returnCanonicalIndices(size_t index) {
    std::vector<size_t> indices;
    switch (index) {
        case 16: // AX
            indices = {0, 1, 2, 3};
            break;
        case 17: // CX
            indices = {4, 5, 6, 7};
            break;
        case 18: // GX
            indices = {8, 9, 10, 11};
            break;
        case 19: // TX/UX
            indices = {12, 13, 14, 15};
            break;
        case 20: // XA, XC, XG, XU/XT, XX (merged for 21-state layout; was 20–24 in 25-state)
            indices = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
            break;
        default:
            Debug(Debug::ERROR) << "Unexpected index for dinucleotide matrix\n";
            EXIT(EXIT_FAILURE);
    }
    return indices;
}

void BaseMatrix::generateSubMatrix(double ** probMatrix, double ** subMatrix, float ** subMatrixPseudoCounts,
                                   int size, bool containsX, std::string matrixName) {

    // calculate background distribution for the amino acids
    double *pBack = new double[size];
    computeBackground(probMatrix, pBack, size, containsX);

    //Precompute matrix R for amino acid pseudocounts:
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            subMatrixPseudoCounts[i][j] = probMatrix[i][j] / (pBack[j]); //subMatrixPseudoCounts[a][b]=P(a|b)
        }
    }

    // calculate the substitution matrix
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            subMatrix[i][j] = std::log2(probMatrix[i][j] / (pBack[i] * pBack[j]));
            //printf("%3.4f %3.4f %3.4f \n",  probMatrix[i][j], pBack[i], pBack[j] );

        }
    }

    // // If the matrixName is dinuc.out, recalculate the non-canonical dinucleotide scores
    // if ("dinuc.out" == matrixName) {
    //     std::vector<size_t> dinucs_row, dinucs_col;
    //     // Without 'XX'
    //     for (int i = 16; i < size; i++) {
    //         // Get the indexes that has the same canonical nucleotides (i.e. "AA", "AC", "AG", "AU/T" for "AX", etc.)
    //         dinucs_row = returnCanonicalIndices(i);
    //         for (int j = 0; j < size; j++) {
    //             dinucs_col.clear();
    //             if (j > 15) {
    //                 dinucs_col = returnCanonicalIndices(j);
    //             }
    //             // If dinucs_col is empty
    //             if (dinucs_col.size() == 0) {
    //                 // j <= 15
    //                 double score_sum = 0.0;
    //                 for (size_t row_idx : dinucs_row) {
    //                     score_sum += subMatrix[row_idx][j];
    //                 }
    //                 subMatrix[i][j] = score_sum / static_cast<double>(dinucs_row.size());
    //                 subMatrix[j][i] = subMatrix[i][j];
    //             } else {
    //                 double score_sum = 0.0;
    //                 for (size_t row_idx : dinucs_row) {
    //                     for (size_t col_idx : dinucs_col) {
    //                         score_sum += subMatrix[row_idx][col_idx];
    //                     }
    //                 }
    //                 subMatrix[i][j] = score_sum / static_cast<double>(dinucs_row.size() * dinucs_col.size());
    //                 subMatrix[j][i] = subMatrix[i][j];
    //             }
    //         }
    //     }
    // }

    delete[] pBack;
//    subMatrix[size - 1][size - 1] = -.7;
//    for (int i = 0; i < size; i++) {
//        subMatrix[size - 1][i] = -.7;
//        subMatrix[i][size - 1] = -.7;
//    }
}

// made non-static for testing purpose
void BaseMatrix::generateSubMatrix(double ** probMatrix, float ** subMatrixPseudoCounts, short ** subMatrix, int size, bool containsX, double bitFactor, double scoringBias, std::string matrixName){
    double** sm = new double* [size];
    for (int i = 0; i < size; i++)
        sm[i] = new double[size];

    generateSubMatrix(probMatrix, sm, subMatrixPseudoCounts, size, containsX, matrixName);

    // convert to short data type matrix
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            double pValNBitScale = (bitFactor * sm[i][j] + scoringBias);
            subMatrix[i][j] = (pValNBitScale < 0.0) ? pValNBitScale - 0.5 : pValNBitScale + 0.5;
        }
    }
    scoreBias = scoringBias;
    for (int i = 0; i < size; i++)
        delete[] sm[i];
    delete[] sm;
}

std::string BaseMatrix::getMatrixName() {
    return matrixName;
}

double BaseMatrix::getBackgroundProb(size_t)  {
    Debug(Debug::ERROR) << "getBackground is not Impl. for this type of Matrix \n";
    EXIT(EXIT_FAILURE);
}

size_t BaseMatrix::memorySize(std::string & matrixName, std::string & matrixData){
    size_t matrixDataSize = matrixData.size() * sizeof(char);
    size_t matrixNameSize = matrixName.size() * sizeof(char);
    return matrixDataSize + 1 + matrixNameSize;
}

char * BaseMatrix::serialize(std::string &matrixName, std::string &matrixData ) {
    char* data = (char*) malloc(memorySize(matrixName, matrixData) + 1);
    char* p = data;
    memcpy(p, matrixName.c_str(), matrixName.size() * sizeof(char));
    p += (matrixName.size() * sizeof(char));
    memcpy(p, ":", 1);
    p += 1;
    memcpy(p, matrixData.c_str(), matrixData.size() * sizeof(char));
    p += (matrixData.size() * sizeof(char));;
    memcpy(p, "\0", 1);
    return data;
}

std::pair<std::string, std::string> BaseMatrix::unserialize(const char * data){
    size_t len=0;

    while(data[len] != '\0'){
        len++;
    }
    std::string matrixName;
    std::string matrixData;
    bool found = false;
    for(size_t pos = 0; pos < std::max(len, (size_t)4) - 4 && found == false; pos++){
        if(data[pos] == '.'
           && data[pos+1] == 'o'
           && data[pos+2] == 'u'
           && data[pos+3] == 't'
           && data[pos+4] == ':' ){
            matrixName = std::string(data, pos+4);
            matrixData = std::string(&data[pos+5]);
            found = true;
        }
    }

    if(found == false){
        matrixName = std::string(data);
    }
    return std::make_pair(matrixName, matrixData);
}

std::string BaseMatrix::unserializeName(const char * data) {
    size_t len = 0;
    while(data[len] != '\0'){
        len++;
    }
    for (size_t pos = 0; pos < std::max(len, (size_t) 4) - 4; pos++) {
        if (data[pos] == '.'
            && data[pos + 1] == 'o'
            && data[pos + 2] == 'u'
            && data[pos + 3] == 't'
            && data[pos + 4] == ':') {
            return std::string(data, pos + 4);
        }
    }
    return data;
}
