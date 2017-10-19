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
    eCD, // Codominant
    eR,  // Dominant
    eD,  // Recessive
    eA   // Allelic
};

// Tests table. Contains information about the sliding windows for each test.
// All arrays have the same dimensions.
struct TestsData {
    // ! change to hyp type
    string ID; // Defines the grouping of the genotype. ID = {"cd", "d", "r", "a"}
    int lower; // Lower bound of the current window
    int upper; // Upper bound of the current window
};

// New format for the Top Hits Table
struct TableEntry {
    string test_index;        // From .csv file, not used for computations
    string chromosome_index;  // From .csv file, not used for computations
    string ID;                // Defines the grouping of the genotype. ID = {"cd", "d", "r", "a"}, comes from user(?)
    int lower;                // Lower bound of the current window, from .csv file
    int upper;                // Upper bound of the current window, from .csv file
    double adjusted_p_value;  // Place for the results of computations
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
                                       vector<unsigned short> const &A, ExecutionParameters const &cont);
    static vector<unsigned short> mapPhenotypeValuesToChar(vector<string> const &phenotype);
private:
    static void prepareData(vector<vector<unsigned short>> & cur_G, vector<unsigned short> & cur_A, InputData & G,
                           vector<unsigned short> const & A, TestsData cur_test);
    static double calcPValue(vector<vector<unsigned short>> const & cur_G, vector<unsigned short> const & cur_A, string ID);
    static AlternativeHypothesisType hashIt (string const& inString);
    static int calcNumElem(vector<unsigned short> const & phenotype);
    static int calcNumElem(vector<unsigned short> const & genotype, unsigned short p);
    static void doubleSizeOfMatrices(vector<vector<unsigned short>> & cur_G, vector<unsigned short> & cur_A);
    static void fillVMatrix(vector<unsigned short> const & cur_G, vector<unsigned short> const & cur_A,
                             vector<vector<int>> & V);
    static double calculateChiSqr(vector<vector<int>> const & V, int V_rows, int V_cols);
    static int myRandom(int i);
    static vector<unsigned short> phenotypeRandomPermutation(vector<unsigned short>& phenotype);
    static bool checkNumElem(vector<unsigned short> const & genotype, unsigned short p);
    static int calcNumElementsInGenotype(vector<vector<int>> V);

    FRIEND_TEST(pvalue_check, calculate_number_of_elements_vector_checked);
    FRIEND_TEST(pvalue_check, calculate_number_of_elements_genotype_vector_checked);
    FRIEND_TEST(pvalue_check, calculate_number_of_elements_vector_all_zeros_checked);
    FRIEND_TEST(pvalue_check, calculate_number_of_elements_genotype_vector_all_zeros_checked);
    FRIEND_TEST(pvalue_check, calculate_number_of_elements_genotype_vector_all_threes_checked);
    FRIEND_TEST(pvalue_check, calculate_number_of_elements_empty_vector_check);
    FRIEND_TEST(pvalue_check, calculate_number_of_elements_empty_genotype_vector_check);
    FRIEND_TEST(pvalue_check, prepare_data_cd_check);
    FRIEND_TEST(pvalue_check, prepare_data_r_check);
    FRIEND_TEST(pvalue_check, prepare_data_d_check);
    FRIEND_TEST(pvalue_check, prepare_data_a_check);
    FRIEND_TEST(pvalue_check, fill_v_matrix_3x4_check);
    FRIEND_TEST(pvalue_check, fill_v_matrix_2x4_check);
    FRIEND_TEST(pvalue_check, fill_v_matrix_with_3_in_genotype_check);
    FRIEND_TEST(pvalue_check, calculate_chi_sqr_V2x4);
    FRIEND_TEST(pvalue_check, calculate_chi_sqr_V3x4);
    FRIEND_TEST(pvalue_check, calculate_chi_sqr_V3x4_proportional);
    FRIEND_TEST(pvalue_check, calculate_chi_sqr_empty_matrix);
    FRIEND_TEST(pvalue_check, calculate_chi_sqr_row_of_zeros);
    FRIEND_TEST(pvalue_check, calculate_chi_sqr_column_of_zeros);
    FRIEND_TEST(pvalue_check, calc_p_value_check);
    FRIEND_TEST(pvalue_check, prepare_data_check_empty_genotype_matrix);
    FRIEND_TEST(pvalue_check, prepare_data_check_empty_phenotype_vector);
    FRIEND_TEST(pvalue_check, check_num_elem);
    FRIEND_TEST(pvalue_check, check_num_elem_2);
    FRIEND_TEST(pvalue_check, check_num_elem_3);
    FRIEND_TEST(pvalue_check, check_num_elem_4);
    FRIEND_TEST(pvalue_check, check_num_elem_5);
};

#endif //ADJUSTPVALUE_PVALUE_H
