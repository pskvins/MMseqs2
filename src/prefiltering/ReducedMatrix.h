#ifndef ReducedMatrix_H
#define ReducedMatrix_H
#include "BaseMatrix.h"
#include "Debug.h"

class ReducedMatrix : public BaseMatrix {
    public:
        ReducedMatrix(double **probMatrix, float ** rMatrix,
                      unsigned char* aa2num, char* num2aa, size_t orgAlphabetSize,
                      size_t reducedAlphabetSize, float bitFactor);
        virtual ~ReducedMatrix();

        void setupLetterMapping() {
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
                                        this->aa2num[static_cast<int>(letter)] = this->aa2num[static_cast<int>('L')];
                                break;
                                case 'U':
                                case 'O':
                                        this->aa2num[static_cast<int>(letter)] = this->aa2num[static_cast<int>('X')];
                                break;
                                case 'Z': this->aa2num[static_cast<int>(letter)] = this->aa2num[static_cast<int>('E')]; break;
                                case 'B': this->aa2num[static_cast<int>(letter)] = this->aa2num[static_cast<int>('D')]; break;
                                default:
                                        this->aa2num[static_cast<int>(letter)] = this->aa2num[static_cast<int>('X')];
                                break;
                        }
                }
        };

        void setupDinucleotideLetterMapping(){
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
        };
    private:

        /*contains the original matrix before the alphabet reduction*/
        unsigned char* orig_aa2num;
        char*  orig_num2aa;
        /* size of the original alphabet*/
        size_t origAlphabetSize;

        // base class aa2num and num2aa mappings contain now:
        // aa2num: mapping aa (orig. alphabet) -> int code of the representative amino acid
        // num2aa: mapping int code (orig. alphabet) -> the representative amino acid char

        // reducedAlphabet contains only the "representative" amino acids
        std::vector<char> reducedAlphabet;

        /* The function adds two rows of a given m x n input matrix and produces a
         * m-1 x n matrix. row1 and row2 (where row1 < row2) are the rows that are
         * to be added. row1 is replaced by row1 + row2 in the output matrix.
         */
        void addTwoRows(double ** input, double ** output, size_t size, size_t row1, size_t row2 );
        /* The function adds two columns of a given m x n input matrix and produces a
         * m x n-1 matrix. col1 and col2 (where col1 < col2) are the columns that are
         * to be added. col1 is replaced by col1 + col2 in the output matrix
         */
        void addTwoColumns(double ** input, double ** output, size_t size, size_t col1, size_t col2 );

        /* Copy from array to array */
        void copyMatrix(double ** input,double ** output, size_t size);
        /* This function generates the corresponding substitution matrix given the corresponding probability
         * matrix.
         */
        void coupleBases(double ** input, double ** output, size_t size, size_t base1, size_t base2);
        /* This function calculates the mutual information of a given substitution
         * matrix (subMatrix) given that its Probability matrix (pMatrix) is also known.
         * numRows and numCols are the number of rows and columns respectively in the
         * matrices that contain meaningful information.
         *
         * Mutual information measures the information that two random variables X and Y share:
         * it measures how much knowing one of these variables reduces uncertainty about the other.
         * (Knowing the amino acid X - how much we know about the aligned amino acid?).
         */
        double calculateMutualInformation(double ** pMatrix, double ** subMatrix, size_t size);
        /* This function finds the two best bases to couple such that we loose the minimum amount of information.
         * Returns the amount of mutual information in the best pairing.
         */
        std::pair<size_t,size_t> coupleWithBestInfo(double ** pinput, double ** pMatrix, float ** rMatrix, size_t size);

};


#endif
