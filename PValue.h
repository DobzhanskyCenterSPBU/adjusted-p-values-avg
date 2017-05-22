//
// Created by EL on 22.05.17.
//

#ifndef ADJUSTPVALUE_PVALUE_H
#define ADJUSTPVALUE_PVALUE_H


#include <vector>
#include <string>
#include "InputData.h"
using namespace std;

// Defines the grouping of the genotype = {"cd", "d", "r", "a"}
enum AlternativeHypothesisType {
    eCD, // add full names
    eR,
    eD,
    eA
};

// Tests table. Contains information about the sliding windows for each test.
// All arrays have the same dimensions.
struct TestsData {
    // ! change to hyp type
    vector<string> ID; // Defines the grouping of the genotype. ID = {"cd", "d", "r", "a"}
    vector<int> lower; // Lower bound of the current window
    vector<int> upper; // Upper bound of the current window
};

// Execution parameters
struct ExecutionParameters {
    int maxReplications; // Maximum amount of replications
    int k; // Maximum amount of P-values (D_cur) that are larger than the initial P-value (D_main)
    bool isAdaptive; // Defines whether K will be used
};

class PValue {
public:
    static vector<double> adjustPValue(TestsData const & tests, InputData & G,
                                       vector<char> const & A, ExecutionParameters const & cont);
private:
    static int prepareData(vector<char> & cur_G, vector<char> & cur_A, InputData & G,
                           vector<char> const & A, int low, int up, string id);
    static double calcPValue(vector<char> const & cur_G, vector<char> const & cur_A);
    static AlternativeHypothesisType hashIt (string const& inString);
    static int calcNumElem(vector<char>::const_iterator first, vector<char>::const_iterator last,
                           char p = ' ');
};


#endif //ADJUSTPVALUE_PVALUE_H
