#include "SubstitutionMatrix.h"
#include "Util.h"
#include "Debug.h"
#include "LambdaCalculation.h"

#include <cstring>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <climits>

SubstitutionMatrix::SubstitutionMatrix(const char *filename, float bitFactor, float scoreBias) : bitFactor(bitFactor) {
    std::pair<std::string, std::string> parsedMatrix = BaseMatrix::unserialize(filename);
    if(parsedMatrix.second != "") {
        // the filename can contain the substituion matrix
        // SUBMATNAME.out:DATA
        // this is used for index databases
        matrixName = parsedMatrix.first;
        matrixData = parsedMatrix.second;
    } else {
        // read amino acid substitution matrix from file
        std::string fileName(parsedMatrix.first.c_str());
        matrixName = Util::base_name(fileName, "/\\");
        if (fileName.length() < 4 || fileName.substr(fileName.length() - 4, 4).compare(".out") != 0) {
            Debug(Debug::ERROR) << "Invalid format of the substitution matrix input file! Only .out files are accepted.\n";
            EXIT(EXIT_FAILURE);
        }
        std::ifstream in(fileName);
        if (in.fail()) {
            Debug(Debug::ERROR) << "Cannot read " << filename << "\n";
            EXIT(EXIT_FAILURE);
        }
        matrixData = std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        in.close();
    }

    std::pair<int, bool> alphSizeAndX = setAaMappingDetectAlphSize(matrixData);
    alphabetSize = alphSizeAndX.first;
    if(alphabetSize == -1){
        Debug(Debug::ERROR) << "Could not estimate alphabet size.\n";
        EXIT(EXIT_FAILURE);
    }
    initMatrixMemory(alphabetSize);
    readProbMatrix(matrixData, alphSizeAndX.second);
    // if(mappingHasAminoAcidLetters()){
        // setupLetterMapping();
    setupDinucleotideLetterMapping();
    // }else {
    //     for (int letter = 0; letter < UCHAR_MAX; letter++) {
    //         char upperLetter = toupper(static_cast<char>(letter));
    //         aa2num[letter] = (aa2num[static_cast<unsigned char>(upperLetter)] == UCHAR_MAX)
    //                          ? alphabetSize-1 : aa2num[static_cast<int>(upperLetter)];
    //     }
    // }

    //print(probMatrix, num2aa, alphabetSize);
    generateSubMatrix(probMatrix, subMatrixPseudoCounts, subMatrix, alphabetSize, true, bitFactor, scoreBias);
}


bool SubstitutionMatrix::estimateLambdaAndBackground(
    const double** scoreMatrix,
    int alphabetSize,
    double* pBack,
    double& lambda
) {
    std::vector<double> letterProbs1(alphabetSize, 0);
    std::vector<double> letterProbs2(alphabetSize, 0);
    lambda = calculate_lambda(scoreMatrix, alphabetSize, letterProbs1, letterProbs2);
    for (int i = 0; i < alphabetSize; i++) {
        pBack[i] = letterProbs1[i];
    }
    if (lambda < 0)
        return false; //bad
    else
        return true; //good
}


void SubstitutionMatrix::calcLocalAaBiasCorrection(const BaseMatrix *m,
                                                   const unsigned char *int_sequence,
                                                   const int N,
                                                   float *compositionBias,
                                                   float scale) {
    const int windowSize = 50;
    for (int i = 0; i < N; i++) {
        const int minPos = std::max(0, (i - windowSize / 2));
        const int maxPos = std::min(N, (i + windowSize / 2));
        const int windowLength = maxPos - minPos;

        // negative score for the amino acids in the neighborhood of i
        int sumSubScores = 0;
        short *subMat = m->subMatrix[int_sequence[i]];
        for (int j = minPos; j < maxPos; j++) {
            sumSubScores += subMat[int_sequence[j]];
        }

        // remove own amino acid
        sumSubScores -= subMat[int_sequence[i]];
        float deltaS_i = (float) sumSubScores;
        // negative avg.
        deltaS_i /= -1.0 * static_cast<float>(windowLength);
        // positive score for the background score distribution for i
        for (int a = 0; a < m->alphabetSize; a++) {
            deltaS_i += m->pBack[a] * static_cast<float>(subMat[a]);
        }
        compositionBias[i] = scale * deltaS_i;
//        std::cout << i << " " << compositionBias[i] << std::endl;
    }
}


void SubstitutionMatrix::calcProfileProfileLocalAaBiasCorrection(short *profileScores,
                                                                 const size_t profileAASize,
                                                                 const int N, size_t alphabetSize) {

    const int windowSize = 50;

    float * pnul  = new float[alphabetSize];
    float * aaSum = new float[alphabetSize];

    memset(pnul, 0, sizeof(float) * alphabetSize);


    for (int pos = 0; pos < N; pos++) {
        const short * subMat = profileScores + (pos * profileAASize);
        for(size_t aa = 0; aa < alphabetSize; aa++) {
            pnul[aa] += subMat[aa]  ;
        }
    }
    for(size_t aa = 0; aa < alphabetSize; aa++)
        pnul[aa] /= N;
    for (int i = 0; i < N; i++){
        const int minPos = std::max(0, (i - windowSize/2));
        const int maxPos = std::min(N, (i + windowSize/2));
        const int windowLength = maxPos - minPos;
        // negative score for the amino acids in the neighborhood of i
        memset(aaSum, 0, sizeof(float) * alphabetSize);

        for (int j = minPos; j < maxPos; j++){
            const short * subMat = profileScores + (j * profileAASize);
            if( i == j )
                continue;
            for(size_t aa = 0; aa < alphabetSize; aa++){
                aaSum[aa] += subMat[aa] - pnul[aa];
            }
        }
        for(size_t aa = 0; aa < alphabetSize; aa++) {
            profileScores[i*profileAASize + aa] = static_cast<int>((profileScores + (i * profileAASize))[aa] - aaSum[aa]/windowLength);
        }
    }
    delete [] aaSum;
    delete [] pnul;
}

void SubstitutionMatrix::calcProfileProfileLocalAaBiasCorrectionAln(int8_t *profileScores,
                                                                    unsigned int N, size_t alphabetSize, BaseMatrix *subMat) {

    const int windowSize = 50;

    float * pnul = new float[N]; // expected score of the prof ag. a random (blosum bg dist) seq
    memset(pnul, 0, sizeof(float) * N);
    float * aaSum = new float[alphabetSize];

    ProfileStates ps(alphabetSize,subMat->pBack);

    for (unsigned int pos = 0; pos < N; pos++) {
        for(size_t aa = 0; aa < alphabetSize; aa++) {
            pnul[pos] += *(profileScores + pos + N*aa) * ps.prior[aa];
        }
    }

    for (unsigned int i = 0; i < N; i++){
        const int minPos = std::max(0, ((int)i - windowSize/2));
        const unsigned int maxPos = std::min(N, (i + windowSize/2));
        const int windowLength = maxPos - minPos;
        // negative score for the amino acids in the neighborhood of i
        memset(aaSum, 0, sizeof(float) * alphabetSize);

        for (unsigned int j = minPos; j < maxPos; j++){
            if (i == j) {
                continue;
            }
            for (size_t aa = 0; aa < alphabetSize; aa++) {
                aaSum[aa] += *(profileScores + aa*N + j) - pnul[j];
            }
        }
        for (size_t aa = 0; aa < alphabetSize; aa++) {
            profileScores[i + aa*N] = static_cast<int8_t>(*(profileScores + i + N*aa) - aaSum[aa]/windowLength);
        }
    }
    delete[] aaSum;
    delete[] pnul;
}



/* Compute aa correction
   => p(a) =  ( \prod_{i=1}^L pi(a) )^(1/L)
   => p(a) = 2^[ (1/L) * \log2 ( \prod_{i=1}^L pi(a) )
   => p(a) = 2^[ (1/L) * \sum_{i=1}^L  \log2 pi(a) ]
   => p(a) = 2^[ (1/L) * \sum_{i=1}^L  \log2 ( pi(a) / f(a) )  + log2 f(a) ]
   => p(a) = f(a) * 2^[ (1/L) * \sum_{i=1}^L  S(pi, a) ]
 */

void SubstitutionMatrix::calcGlobalAaBiasCorrection(const BaseMatrix *m,
                                                    char *profileScores,
                                                    float *pNullBuffer,
                                                    const size_t profileAASize,
                                                    const int N) {
    memset(pNullBuffer, 0, sizeof(float) * N);
    const int windowSize = 50;
    for (int pos = 0; pos < N; pos++) {
        const char * subMat = profileScores + (pos * profileAASize);
        for(size_t aa = 0; aa < 20; aa++) {
            pNullBuffer[pos] += m->pBack[aa] * static_cast<float>(subMat[aa]);
        }
    }
//    for(size_t aa = 0; aa < 20; aa++)
//        pNullBuffer[aa] /= N;
    for (int i = 0; i < N; i++) {
        const int minPos = std::max(0, (i - windowSize / 2));
        const int maxPos = std::min(N, (i + windowSize / 2));
        const int windowLength = maxPos - minPos;
        // negative score for the amino acids in the neighborhood of i
        float aaSum[20];
        memset(aaSum, 0, sizeof(float) * 20);

        for (int j = minPos; j < maxPos; j++) {
            const char *subMat = profileScores + (j * profileAASize);
            if (i == j) {
                continue;
            }
            for (size_t aa = 0; aa < 20; aa++) {
                aaSum[aa] += subMat[aa] - pNullBuffer[j];
            }
        }
        for (size_t aa = 0; aa < 20; aa++) {
            profileScores[i * profileAASize + aa] = static_cast<int>((profileScores + (i * profileAASize))[aa] -
                                                                     aaSum[aa] / windowLength);
//            avg += static_cast<int>((profileScores + (i * profileAASize))[aa] -  aaSum[aa]/windowLength);
        }
    }
}


SubstitutionMatrix::~SubstitutionMatrix() {}

bool SubstitutionMatrix::mappingHasAminoAcidLetters(){
    std::string lettersToCheck = "ATGCDEFHIKLMNPQRSVWYX";
    size_t cnt = 0;
    for(size_t i = 0; i < lettersToCheck.size(); i++){
        cnt += (aa2num[static_cast<int>(lettersToCheck[i])] != UCHAR_MAX);
    }
    return (cnt == lettersToCheck.size());
}

void SubstitutionMatrix::setupLetterMapping(){
    for(int letter = 0; letter < UCHAR_MAX; letter++){
        char upperLetter = toupper(static_cast<char>(letter));
        switch(upperLetter){
            case 'A':
            case 'T':
            case 'G':
            case 'C':
            case 'D':
            case 'E':
            case 'F':
            case 'H':
            case 'I':
            case 'K':
            case 'L':
            case 'M':
            case 'N':
            case 'P':
            case 'Q':
            case 'R':
            case 'S':
            case 'V':
            case 'W':
            case 'Y':
            case 'X':
                this->aa2num[static_cast<int>(letter)] = this->aa2num[static_cast<int>(upperLetter)];
                break;
            case 'J':
                this->aa2num[static_cast<int>(letter)] = this->aa2num[(int)'L'];
                break;
            case 'U':
            case 'O':
                this->aa2num[static_cast<int>(letter)] = this->aa2num[(int)'X'];
                break;
            case 'Z': this->aa2num[static_cast<int>(letter)] = this->aa2num[(int)'E']; break;
            case 'B': this->aa2num[static_cast<int>(letter)] = this->aa2num[(int)'D']; break;
            default:
                this->aa2num[static_cast<int>(letter)] = this->aa2num[(int)'X'];
                break;
        }
    }
}

void SubstitutionMatrix::setupDinucleotideLetterMapping(){
    // Map dinucleotides
    for (int i = 1; i < UCHAR_MAX; i++) { // Never map 0
        unsigned short upperLetter1 = (static_cast<unsigned short>(toupper(static_cast<unsigned char>(i))) << 8);
        unsigned short first = (static_cast<unsigned short>(i) << 8);
        for (int j = 0; j < UCHAR_MAX; j++) {
            unsigned short dinucleotide = first | static_cast<unsigned char>(j);
            unsigned char upperLetter2 = toupper(static_cast<unsigned char>(j));
            /* Map dinuc2alph = {'AA': 'C', 'AC': 'G', 'AG': 'L', 'AU': 'Q',
                                'CA': 'D', 'CC': 'F', 'CG': 'R', 'CU': 'K',
                                'GA': 'M', 'GC': 'A', 'GG': 'P', 'GU': 'I',
                                'UA': 'E', 'UC': 'N', 'UG': 'H', 'UU': 'S'} */
            unsigned short upperLetters = upperLetter1 | upperLetter2;
            switch(upperLetters){
                case 0b0100000101000001: // AA
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'C'];
                    break;
                case 0b0100000101000011: // AC
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'G'];
                    break;
                case 0b0100000101000111: // AG
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'L'];
                    break;
                case 0b0100000101010100: // AT
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'Q'];
                    break;
                case 0b0100000101010101: // AU
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'Q'];
                    break;
                case 0b0100001101000001: // CA
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'D'];
                    break;
                case 0b0100001101000011: // CC
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'F'];
                    break;
                case 0b0100001101000111: // CG
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'R'];
                    break;
                case 0b0100001101010100: // CT
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'K'];
                    break;
                case 0b0100001101010101: // CU
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'K'];
                    break;
                case 0b0100011101000001: // GA
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'M'];
                    break;
                case 0b0100011101000011: // GC
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'A'];
                    break;
                case 0b0100011101000111: // GG
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'P'];
                    break;
                case 0b0100011101010100: // GT
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'I'];
                    break;
                case 0b0100011101010101: // GU
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'I'];
                    break;
                case 0b0101010001000001: // TA
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'E'];
                    break;
                case 0b0101010001000011: // TC
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'N'];
                    break;
                case 0b0101010001000111: // TG
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'H'];
                    break;
                case 0b0101010001010100: // TT
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'H'];
                    break;
                case 0b0101010001010101: // TU
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'H'];
                    break;
                case 0b0101010101000001: // UA
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'E'];
                    break;
                case 0b0101010101000011: // UC
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'N'];
                    break;
                case 0b0101010101000111: // UG
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'H'];
                    break;
                case 0b0101010101010100: // UT
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'H'];
                    break;
                case 0b0101010101010101: // UU
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'S'];
                    break; 
                default:
                    this->aa2num[static_cast<int>(dinucleotide)] = this->aa2num[(int)'X'];
                    break;
            }
        }
    }

    // Set up revcomp
    /* {'C': 'S', 'G': 'I', 'L': 'K', 'Q': 'Q', 'D': 'H', 
    'F': 'P', 'R': 'R', 'K': 'L', 'M': 'N', 'A': 'A', 
    'P': 'F', 'I': 'G', 'E': 'E', 'N': 'M', 'H': 'D', 'S': 'C', 'X': 'X'} */
    this->revcomp[this->aa2num[(int)'C']] = this->aa2num[(int)'S'];
    this->revcomp[this->aa2num[(int)'G']] = this->aa2num[(int)'I'];
    this->revcomp[this->aa2num[(int)'L']] = this->aa2num[(int)'K'];
    this->revcomp[this->aa2num[(int)'Q']] = this->aa2num[(int)'Q'];
    this->revcomp[this->aa2num[(int)'D']] = this->aa2num[(int)'H'];
    this->revcomp[this->aa2num[(int)'F']] = this->aa2num[(int)'P'];
    this->revcomp[this->aa2num[(int)'R']] = this->aa2num[(int)'R'];
    this->revcomp[this->aa2num[(int)'K']] = this->aa2num[(int)'L'];
    this->revcomp[this->aa2num[(int)'M']] = this->aa2num[(int)'N'];
    this->revcomp[this->aa2num[(int)'A']] = this->aa2num[(int)'A'];
    this->revcomp[this->aa2num[(int)'P']] = this->aa2num[(int)'F'];
    this->revcomp[this->aa2num[(int)'I']] = this->aa2num[(int)'G'];
    this->revcomp[this->aa2num[(int)'E']] = this->aa2num[(int)'E'];
    this->revcomp[this->aa2num[(int)'N']] = this->aa2num[(int)'M'];
    this->revcomp[this->aa2num[(int)'H']] = this->aa2num[(int)'D'];
    this->revcomp[this->aa2num[(int)'S']] = this->aa2num[(int)'C'];
    this->revcomp[this->aa2num[(int)'X']] = this->aa2num[(int)'X'];
    // For paddings
    this->revcomp[16] = 16;
    this->revcomp[17] = 17;
    this->revcomp[18] = 18;
    this->revcomp[19] = 19;
    
    // Set up tail
    /* {'C': 'G', 'G': 'F', 'L': 'A', 'Q': 'N',
    'D': 'G', 'F': 'F', 'R': 'A', 'K': 'N',
    'M': 'G', 'A': 'F', 'P': 'A', 'I': 'N',
    'E': 'G', 'N': 'F', 'H': 'A', 'S': 'N', 'X': 'X'} */
    this->tail[this->aa2num[(int)'C']] = this->aa2num[(int)'G'];
    this->tail[this->aa2num[(int)'G']] = this->aa2num[(int)'F'];
    this->tail[this->aa2num[(int)'L']] = this->aa2num[(int)'A'];
    this->tail[this->aa2num[(int)'Q']] = this->aa2num[(int)'N'];
    this->tail[this->aa2num[(int)'D']] = this->aa2num[(int)'G'];
    this->tail[this->aa2num[(int)'F']] = this->aa2num[(int)'F'];
    this->tail[this->aa2num[(int)'R']] = this->aa2num[(int)'A'];
    this->tail[this->aa2num[(int)'K']] = this->aa2num[(int)'N'];
    this->tail[this->aa2num[(int)'M']] = this->aa2num[(int)'G'];
    this->tail[this->aa2num[(int)'A']] = this->aa2num[(int)'F'];
    this->tail[this->aa2num[(int)'P']] = this->aa2num[(int)'A'];
    this->tail[this->aa2num[(int)'I']] = this->aa2num[(int)'N'];
    this->tail[this->aa2num[(int)'E']] = this->aa2num[(int)'G'];
    this->tail[this->aa2num[(int)'N']] = this->aa2num[(int)'F'];
    this->tail[this->aa2num[(int)'H']] = this->aa2num[(int)'A'];
    this->tail[this->aa2num[(int)'S']] = this->aa2num[(int)'N'];
    this->tail[this->aa2num[(int)'X']] = this->aa2num[(int)'X'];
}


int SubstitutionMatrix::parseAlphabet(char *word, char *num2aa, int *aa2num) {
    char *charReader = word;
    int minAAInt = INT_MAX;
    // find amino acid with minimal int value
    while (isalpha(*charReader)) {
        const char aa = *charReader;
        const int intAA = aa2num[static_cast<int>(aa)];
        minAAInt = std::max(minAAInt, intAA);
        charReader++;
    }
    if(minAAInt==-1){

    }
    char minAAChar = num2aa[minAAInt];
    // do alphabet reduction
    charReader = word;
    while (isalpha(*charReader)) {
        const char aa = *charReader;
        const int intAA = aa2num[static_cast<int>(aa)];
        aa2num[static_cast<int>(aa)] = minAAInt;
        num2aa[intAA] = minAAChar;
        charReader++;
    }
    return minAAInt;
}

bool SubstitutionMatrix::printLambdaAndBackground = false;
void SubstitutionMatrix::readProbMatrix(const std::string &matrixData, const bool containsX) {
    std::stringstream in(matrixData);
    std::string line;
    bool probMatrixStart = false;
    const char *words[256];
    bool hasLambda = false;
    bool hasBackground = false;
    while (in.good()) {
        getline(in, line);
        size_t wordCnt = Util::getWordsOfLine((char *) line.c_str(), words, 256);
        // skip comments
        if (line[0] == '#') {
            if (line.find("# Background (precomputed optional):") == 0) {
                for (size_t i = 4; i < wordCnt; i++) {
                    double f = strtod(words[i], NULL);
                    pBack[i-4] = f;
                }
                hasBackground = true;
            }
            if (line.find("# Lambda     (precomputed optional):") == 0) {
                double f = strtod(words[4], NULL);
                lambda = f;
                hasLambda = true;
            }
            continue;
        }

        if (wordCnt > 1 && probMatrixStart == false) {
            probMatrixStart = true;
            continue;
        }
        if (wordCnt > 1 && probMatrixStart == true) {
            if (isalpha(words[0][0]) == false) {
                Debug(Debug::ERROR) << "First element in probability line must be an alphabet letter.\n";
                EXIT(EXIT_FAILURE);
            }
            int aa = static_cast<int>(aa2num[toupper(words[0][0])]);
            for (int i = 0; i < alphabetSize; i++) {
                double f = strtod(words[i + 1], NULL);
                probMatrix[aa][i] = f; // divided by 2 because we scale bit/2 ot bit
            }
        }
    }
    bool xIsPositive = false;
    if( containsX == true ){
        for (int j = 0; j < alphabetSize; j++) {
            int xIndex = static_cast<int>(aa2num[static_cast<int>('X')]);
            if ((probMatrix[xIndex][j] > 0) || (probMatrix[j][xIndex] > 0)) {
                xIsPositive = true;
                break;
            }
        }
    }


    if (containsX == false) {
        Debug(Debug::ERROR) << "Please add X to your substitution matrix.\n";
        EXIT(EXIT_FAILURE);
    }

    if(hasLambda == false || hasBackground == false){
        if (estimateLambdaAndBackground(const_cast<const double **>(probMatrix), alphabetSize - ((xIsPositive) ? 0 : 1),
                                        pBack, lambda) == false) {
            Debug(Debug::ERROR) << "Computing inverse of substitution matrix failed\n";
            EXIT(EXIT_FAILURE);
        }

        pBack[static_cast<int>(aa2num[static_cast<int>('X')])] = ANY_BACK;

        if (printLambdaAndBackground) {
            Debug(Debug::INFO) << "# Background (precomputed optional):";
            for (int i = 0; i < alphabetSize; ++i) {
                Debug(Debug::INFO) << " " << SSTR((float)pBack[i], 5);
            }
            Debug(Debug::INFO) << "\n";
            Debug(Debug::INFO) << "# Lambda     (precomputed optional): " << SSTR((float)lambda, 5) << "\n";
        }
    }
    if(xIsPositive == false){
        for (int i = 0; i < alphabetSize - 1; i++) {
            pBack[i] = pBack[i] * (1.0 - pBack[static_cast<int>(aa2num[static_cast<int>('X')])]);
        }
    }
    // Reconstruct Probability Sab=(Pab/Pa*Pb) -> Pab = exp(Sab) * Pa * Pb
    for (int i = 0; i < alphabetSize; i++) {
        //smat[i] = smatData+((subMat.alphabetSize-1)*i);
        for (int j = 0; j < alphabetSize; j++) {
            probMatrix[i][j] = std::exp(lambda * probMatrix[i][j]) * pBack[i] * pBack[j];
        }
    }
}

std::pair<int, bool> SubstitutionMatrix::setAaMappingDetectAlphSize(std::string &matrixData){
    std::stringstream in(matrixData);
    std::string line;
    const char *words[256];
    int alphabetSize = 0;
    bool containsX;
    while (in.good()) {
        getline(in, line);
        size_t wordCnt = Util::getWordsOfLine((char *) line.c_str(), words, 256);

        if (line[0] == '#') {
            continue;
        }
        if (wordCnt > 1) {
            for (size_t i = 0; i < wordCnt; i++) {
                if (isalpha(words[i][0]) == false) {
                    Debug(Debug::ERROR) << "Probability matrix must start with alphabet header.\n";
                    EXIT(EXIT_FAILURE);
                }
                int aa = toupper(words[i][0]);
                aa2num[aa] = static_cast<unsigned char>(i);
                num2aa[i] = aa;
                if (aa == 'X') {
                    containsX = true;
                }
//                column_aa[i] = parseAlphabet(words[i], num2aa, aa2num);
            }
            alphabetSize = wordCnt;
            return std::make_pair(alphabetSize, containsX);
        }
    }
    return std::make_pair(-1, false);
}




