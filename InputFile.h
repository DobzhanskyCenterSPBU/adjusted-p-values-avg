//
// Created by EL on 08.05.17.
//

#ifndef P_VALUE_INPUTFILE_H
#define P_VALUE_INPUTFILE_H

#include "InputData.h"
#include <string>
#include <fstream>
#include <vector>

using namespace std;

class InputFile: public InputData {
public:
    ifstream handle;
    InputFile (string fileName);
    vector<vector<unsigned char>> createGenotypeMatrix(int lowInd, int upInd);
    ~InputFile();
};

#endif //P_VALUE_INPUTFILE_H
