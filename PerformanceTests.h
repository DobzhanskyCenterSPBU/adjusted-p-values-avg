//
// Created by EL on 06.06.17.
//

#include <string>
#include <vector>
#include "PValue.h"

using namespace std;

#ifndef ADJUSTPVALUE_PERFORMANCETESTS_H
#define ADJUSTPVALUE_PERFORMANCETESTS_H
void simplePerformanceTest(int maxRep, string test_type);
vector<unsigned short> createPhenotypeVector(string filename);
void realDataPerformanceTest(int maxRep, string genotypeFile, string phenotypeFile, vector<TestsData> tests);
void runMultiplePerformanceTests(int maxRep);
#endif //ADJUSTPVALUE_PERFORMANCETESTS_H
