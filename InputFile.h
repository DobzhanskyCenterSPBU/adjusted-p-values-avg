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
    string filename;
public:
    InputFile (string fileName);
    vector<vector<unsigned short>> createGenotypeMatrix(int lowInd, int upInd);
};

#endif //P_VALUE_INPUTFILE_H
