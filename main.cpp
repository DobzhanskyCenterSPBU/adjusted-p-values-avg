#include <iostream>
#include "PValue.h"
#include "InputFile.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace std;

int main(int argc, char* argv[]) {

    testing::InitGoogleTest(&argc, argv);
    RUN_ALL_TESTS();

    // Data initialization
    vector<TestsData> tests;
    tests = {{"cd",2,5}, {"r",6,8}, {"d",10,11}, {"a",12,14}, {"d",16,19}};

    ExecutionParameters parameters;
    parameters.maxReplications = 10; // Should be around 10^(-8)
    parameters.k = 10;
    parameters.isAdaptive = true;

    InputFile genotype("/Users/EL/CLionProjects/AdjustPValue/test.txt");
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

    vector<unsigned char> phenotype = {'1', '3', '0', '2', '1', '0'};

    // Call adjustPValue
    vector<double> results = PValue::adjustPValue(tests, genotype, phenotype, parameters);

    // Print results
    for (vector<double>::iterator it = results.begin(); it != results.end(); ++it) {
        cout << ' ' << *it;
    }
    cout << '\n';

    // Destroy InputFile object

    return 0;
}