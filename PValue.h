//
// Created by EL on 22.05.17.
//

#ifndef ADJUSTPVALUE_PVALUE_H
#define ADJUSTPVALUE_PVALUE_H


#include <vector>
#include <string>
#include "InputData.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
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
    string ID; // Defines the grouping of the genotype. ID = {"cd", "d", "r", "a"}
    int lower; // Lower bound of the current window
    int upper; // Upper bound of the current window
};

// Execution parameters
struct ExecutionParameters {
    int maxReplications; // Maximum amount of replications
    int k; // Maximum amount of P-values (D_cur) that are larger than the initial P-value (D_main)
    bool isAdaptive; // Defines whether K will be used
};

class PValue {
public:
    static vector<double> adjustPValue(vector<TestsData> const &tests, InputData &G,
                                       vector<unsigned char> const &A, ExecutionParameters const &cont);
private:
    static void prepareData(vector<vector<unsigned char>> & cur_G, vector<unsigned char> & cur_A, InputData & G,
                           vector<unsigned char> const & A, TestsData cur_test);
    static double calcPValue(vector<vector<unsigned char>> const & cur_G, vector<unsigned char> const & cur_A, string ID);
    static AlternativeHypothesisType hashIt (string const& inString);
    static int calcNumElem(vector<unsigned char> const & phenotype);
    static int calcNumElem(vector<vector<unsigned char>> const & genotype);
    static void doubleSizeOfMatrices(vector<vector<unsigned char>> & cur_G, vector<unsigned char> & cur_A);

    FRIEND_TEST(pvalue_check, calculate_number_of_elements_vector_checked);
    FRIEND_TEST(pvalue_check, calculate_number_of_elements_matrix_checked);
    FRIEND_TEST(pvalue_check, calculate_number_of_elements_empty_vector_check);
    FRIEND_TEST(pvalue_check, calculate_number_of_elements_empty_matrix_check);
    FRIEND_TEST(pvalue_check, calculate_number_of_elements_vector_all_zeros_checked);
    FRIEND_TEST(pvalue_check, calculate_number_of_elements_matrix_all_zeros_checked);
    FRIEND_TEST(pvalue_check, calculate_number_of_elements_matrix_all_threes_checked);
    FRIEND_TEST(pvalue_check, prepare_data_cd_check);
    FRIEND_TEST(pvalue_check, prepare_data_r_check);
    FRIEND_TEST(pvalue_check, prepare_data_d_check);
    FRIEND_TEST(pvalue_check, prepare_data_a_check);
};

#endif //ADJUSTPVALUE_PVALUE_H
