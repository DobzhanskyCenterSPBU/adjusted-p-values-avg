//
// Created by EL on 12.10.17.
//

#ifndef PLINK_FILE_WORK_INPUTPLINK_H
#define PLINK_FILE_WORK_INPUTPLINK_H

#include "InputData.h"
#include <vector>
using namespace std;

class InputPLINK: public InputData {
    char* filename; // The name of the file without extension (.bed, .bim, .fam files will be read)
public:
    InputPLINK (char* fileName);
    vector<vector<unsigned short>> createGenotypeMatrix(int lowInd, int upInd);
    vector<unsigned short> createPhenotypeVector();
};

#endif //PLINK_FILE_WORK_INPUTPLINK_H