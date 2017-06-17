#include <iostream>
#include "PValue.h"
#include "InputFile.h"
#include <chrono>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "PerformanceTests.h"

#include <boost/math/special_functions/gamma.hpp>

using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[]) {
    // Run unit tests
    testing::InitGoogleTest(&argc, argv);
    RUN_ALL_TESTS();
    //simplePerformanceTest(10, "a");
    realDataPerformanceTest(1000000);


    /*

    // Data initialization
    vector<TestsData> tests;
    tests = {{"cd",2,5}, {"r",6,8}, {"d",10,11}, {"a",12,14}, {"d",16,19}};

    ExecutionParameters parameters;
    parameters.maxReplications = 10; // Should be around 10^(-8)
    parameters.k = 10;
    parameters.isAdaptive = true;

    InputFile genotype("test.txt");
    /*
    vector<vector<unsigned char>> res = genotype.createGenotypeMatrix(-3,1);
    cout << endl << "results" << endl;
    for (vector<vector<unsigned char>>::size_type i = 0; i < res.size(); i++) {
        for (vector<unsigned char>::size_type j = 0; j < res[i].size(); j++ ) {
            cout << res[i][j] << ' ';
        }
        cout << '\n';
    }
    cout << '\n';
    */

    /*

    vector<unsigned short> phenotype = {1, 3, 0, 2, 1, 0};

    // Call adjustPValue
    vector<double> results = PValue::adjustPValue(tests, genotype, phenotype, parameters);

    // Print results
    cout << endl << "Adjusted P-values:" << endl;
    for (vector<double>::iterator it = results.begin(); it != results.end(); ++it) {
        cout << ' ' << *it;
    }
    cout << '\n';


    */

    return 0;
}