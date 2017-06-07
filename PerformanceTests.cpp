//
// Created by EL on 06.06.17.
//

#include "PerformanceTests.h"
#include "PValue.h"
#include "InputFile.h"
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
    for (int rowNum = 0; rowNum < 1000; ++rowNum) {
        cur_line = "";
        for (int colNum = 0; colNum < 100; ++colNum) {
            randomNumber = values[rand() % length];
            cur_line += to_string(randomNumber);
            cur_line += " ";
        }
        myfile << cur_line + "\n";
    }
    myfile.close();

    // Generate the phenotype vector
    vector<unsigned char> phenotype;
    for (int i = 0; i < 100; ++i) {
        randomNumber = values[rand() % length];
        phenotype.push_back(randomNumber + '0');
    }

    vector<TestsData> tests;
    tests = {{test_type,0,999}};
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
