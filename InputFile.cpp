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

vector<char> InputFile::createGenotypeVector(int lowInd, int upInd) {
    vector<char> result;
    char delim = ' '; // Split lines on a delimiter ' '
    string symbol;
    string curLine;

    // Define the length of a line in the given file
    this->handle.seekg(0);
    getline(this->handle, curLine);
    int length = curLine.length();

    // Read all lines between lowInd and upInd
    for (int i = lowInd; i <= upInd; i++) {
        this->handle.seekg((length+1)*i); // Go to the necessary line
        getline(this->handle, curLine);
        //cout << i << ' ' << curLine << endl;

        // Add all symbols from current line to the resulting vector
        stringstream ss;
        ss.str(curLine);
        while (getline(ss, symbol, delim)) {
            result.push_back((char) symbol[0]);
        }
    }
    return result;
}

InputFile::~InputFile(){
    this->handle.close();
}