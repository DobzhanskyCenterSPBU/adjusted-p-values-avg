//
// Created by EL on 06.06.17.
//

#include "PerformanceTests.h"
#include "PValue.h"
#include "InputFile.h"
#include <vector>
#include <string>
#include <iostream>
#include <chrono>

using namespace std;
using namespace std::chrono;



void simplePerformanceTest(int maxRep, string test_type) {

    /*
     * Measures the time it takes to run the adjustPValue function on a G matrix with 1000x100 elements, where
     * the test type is defined by test_type and the maximum number of replications is defined by maxRep.
     */

    // Data initialization

    // Generate genotype matrix and write to file
    string filename = "performancetest.txt";
    int values[] = {0, 1, 2, 3};
    int length = sizeof(values) / sizeof(int);
    int randomNumber;
    ofstream myfile;
    myfile.open (filename);
    string cur_line;
    for (int rowNum = 0; rowNum < 20; ++rowNum) {
        cur_line = "";
        for (int colNum = 0; colNum < 400; ++colNum) {
            randomNumber = values[rand() % length];
            cur_line += to_string(randomNumber);
            cur_line += " ";
        }
        myfile << cur_line + "\n";
    }
    myfile.close();

    // Generate the phenotype vector
    vector<unsigned short> phenotype;
    for (int i = 0; i < 400; ++i) {
        randomNumber = values[rand() % length];
        phenotype.push_back((unsigned short)randomNumber);
    }

    vector<TableEntry> tests;
    tests = {{"","",test_type,0,999,0}};
    ExecutionParameters parameters;
    parameters.maxReplications = maxRep;
    parameters.k = 10;
    parameters.isAdaptive = false;
    InputFile genotype(filename);

    // Measure time
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    // Call adjustPValue
    vector<double> results = PValue::adjustPValue(tests, genotype, phenotype, parameters);

    // Measure time
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    cout << endl << "Duration in microseconds (maxReplications = " << maxRep << ", test type = " << test_type <<
         "): " << duration << endl;

    // Print results
    cout << endl << "Adjusted P-values from performance test:" << endl;
    for (vector<double>::iterator it = results.begin(); it != results.end(); ++it) {
        cout << ' ' << *it;
    }
    cout << '\n';
}


vector<unsigned short> createPhenotypeVector(string filename){
    ifstream handle;
    handle.open(filename);

    if(!handle.is_open())
    {
        cerr << "No such file!" << endl;
        exit(1); // What exit code should be used? How to skip to the next iteration?
    }

    vector<unsigned short> result;
    char delim = ' '; // Split lines on a delimiter ' '
    string symbol;
    string curLine;

    // Define the length of a line in the given file
    handle.seekg(0);
    getline(handle, curLine);

    // Add all symbols from  line to the resulting vector
    stringstream ss;
    ss.str(curLine);
    while (getline(ss, symbol, delim)) {
        result.push_back((unsigned short)stoi(symbol));
    }
    handle.close();
    return result;
}

long long int realDataPerformanceTest(int maxRep, string genotypeFile, string phenotypeFile, vector<TableEntry> tests){

    vector<unsigned short> phenotype = createPhenotypeVector(phenotypeFile);

    ExecutionParameters parameters;
    parameters.maxReplications = maxRep;
    parameters.k = 10;
    parameters.isAdaptive = false;
    InputFile genotype(genotypeFile);

    // Measure time
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    // Call adjustPValue
    vector<double> results = PValue::adjustPValue(tests, genotype, phenotype, parameters);

    // Measure time
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    //return duration;

    cout << endl << "Duration in microseconds (maxReplications = " << maxRep << "): " << duration << endl;


    // Print results
    cout << endl << "Adjusted P-values from real data performance test:" << endl;
    for (vector<double>::iterator it = results.begin(); it != results.end(); ++it) {
        cout << ' ' << *it;
    }
    cout << '\n';

}

void runMultiplePerformanceTests(int maxRep){
    string genotypeName, phenotypeName;


    vector<vector<TableEntry>> tests;
    tests = {{{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}},
             {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 23, 0}}, {{"","","cd", 0, 20, 0}},
             {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}},
             {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 23, 0}}, {{"","","cd", 0, 20, 0}},
             {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}};

    for (int i = 1; i <= 20; ++i) {
        genotypeName = "Tests2/newgenotypedatatest" + to_string(i);
        phenotypeName = "Tests2/newphenotypedatatest" + to_string(i);

        realDataPerformanceTest(maxRep, genotypeName, phenotypeName, tests[i-1]);
    }
}

void runMultiplePerformanceTestsMultipleTimes(int maxRep){
    string genotypeName, phenotypeName;

    vector<vector<TableEntry>> tests2;
    tests = {{{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}},
             {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 23, 0}}, {{"","","cd", 0, 20, 0}},
             {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}},
             {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 23, 0}}, {{"","","cd", 0, 20, 0}},
             {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}, {{"","","cd", 0, 20, 0}}};

    int iter_num = 20;
    long long int avg_time;
    for (int i = 1; i <= 20; ++i) {
        genotypeName = "Tests/newgenotypedatatest" + to_string(i);
        phenotypeName = "Tests/newphenotypedatatest" + to_string(i);

        avg_time = 0;
        for (int j = 0; j < iter_num; ++j) {
            avg_time += realDataPerformanceTest(maxRep, genotypeName, phenotypeName, tests[i-1]);
        }
        cout << "Average time for test " << i << " (in microseconds): " << avg_time/iter_num << endl;
    }
}

