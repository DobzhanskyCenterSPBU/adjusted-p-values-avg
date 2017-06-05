//
// Created by EL on 08.05.17.
//

#include "InputFile.h"
#include <sstream>
#include <iostream>

using namespace std;

InputFile::InputFile (string fileName){
    this->handle.open(fileName);
}

vector<vector<unsigned char>> InputFile::createGenotypeMatrix(int lowInd, int upInd) {

    vector<vector<unsigned char>> result; //(upInd - lowInd + 1)
    char delim = ' '; // Split lines on a delimiter ' '
    string symbol;
    string curLine;
    int line_num = 0;

    if (lowInd < 0) lowInd = 0;
    if (upInd < 0) upInd = 0;

    // Make sure tellg() won't return -1
    this->handle.clear();

    // Define the length of a line in the given file
    this->handle.seekg(0);
    getline(this->handle, curLine);
    int length = (int)(curLine.length());

    // Read all lines between lowInd and upInd
    for (int i = lowInd; i <= upInd; i++) {
        this->handle.seekg((length+1)*i); // Go to the necessary line
        if (handle.tellg() == -1) return result;
        getline(this->handle, curLine);

        // Add all symbols from current line to the resulting matrix
        if (curLine != ""){
            result.push_back({});
            stringstream ss;
            ss.str(curLine);
            while (getline(ss, symbol, delim)) {
                result[line_num].push_back((unsigned char) symbol[0]);
            }
            line_num++;
        }
    }
    return result;
}

InputFile::~InputFile(){
    this->handle.close();
}