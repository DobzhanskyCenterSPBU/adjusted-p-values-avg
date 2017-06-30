//
// Created by EL on 08.05.17.
//

#include "InputFile.h"
#include <sstream>
#include <iostream>

using namespace std;

InputFile::InputFile (string fileName){
    this->filename = fileName;
}

vector<vector<unsigned short>> InputFile::createGenotypeMatrix(int lowInd, int upInd) {

    ifstream handle;
    handle.open(this->filename);

    if(!handle.is_open())
    {
        cerr << "No such file!" << endl;
        exit(1); // What exit code should be used? How to skip to the next iteration?
    }

    vector<vector<unsigned short>> result; //(upInd - lowInd + 1)
    char delim = ' '; // Split lines on a delimiter ' '
    string symbol;
    string curLine;
    int line_num = 0;

    if (lowInd < 0) lowInd = 0;
    if (upInd < 0) upInd = 0;

    // Define the length of a line in the given file
    handle.seekg(0);
    getline(handle, curLine);
    int length = (int)(curLine.length());

    // Read all lines between lowInd and upInd
    for (int i = lowInd; i <= upInd; i++) {
        handle.seekg((length+1)*i); // Go to the necessary line
        if (handle.tellg() == -1) return result;
        getline(handle, curLine);

        // Add all symbols from current line to the resulting matrix
        if (curLine != ""){
            result.push_back({});
            stringstream ss;
            ss.str(curLine);
            while (getline(ss, symbol, delim)) {
                result[line_num].push_back((unsigned short)stoi(symbol));
            }
            line_num++;
        }
    }

    handle.close();

    return result;
}
