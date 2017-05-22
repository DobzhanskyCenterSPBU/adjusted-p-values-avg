#include <iostream>
#include "PValue.h"
#include "InputFile.h"

using namespace std;

int main() {

    // Data initialization
    TestsData tests;
    tests.ID = {"cd","r","d","a","d"};
    tests.lower = {2, 6, 10, 12, 16};
    tests.upper = {5, 8, 11, 14, 19};

    ExecutionParameters parameters;
    parameters.maxReplications = 10; // Should be around 10^(-8)
    parameters.k = 10;
    parameters.isAdaptive = true;

    InputFile genotype("/Users/EL/CLionProjects/P_value/test.txt");
    /*
    vector<char> res = genotype.createGenotypeVector(18,19);
    cout << endl << "results" << endl;
    for (vector<char>::iterator it = res.begin(); it != res.end(); ++it)
        cout << ' ' << *it;
    cout << '\n';
    */

    vector<char> phenotype = {'1', '3', '0', '2', '1', '0'};

    // Call adjustPValue
    vector<double> results = PValue::adjustPValue(tests, genotype, phenotype, parameters);

    // Print results
    for (vector<double >::iterator it = results.begin(); it != results.end(); ++it)
        cout << ' ' << *it;
    cout << '\n';

    // Destroy InputFile object
    return 0;
}