//
// Created by EL on 08.05.17.
//

#ifndef P_VALUE_INPUTDATA_H
#define P_VALUE_INPUTDATA_H

#include <vector>
using namespace std;

class InputData {
public:
    InputData() {};
    ~InputData() {};
    virtual vector<vector<unsigned char>> createGenotypeMatrix(int lowInd, int upInd) = 0;
};

#endif //P_VALUE_INPUTDATA_H
