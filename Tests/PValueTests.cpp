//
// Created by EL on 23.05.17.
//

#include <string>
#include <vector>
#include <boost/math/special_functions/gamma.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../PValue.h"
#include "../InputFile.h"

// Calculating the number of unique elements in a vector
TEST(pvalue_check, calculate_number_of_elements_vector_checked){
    vector<unsigned char> test_vec = {'0', '2', '1', '1', '0', '2', '1', '0', '1', '0', '2', '1'};
    int result = PValue::calcNumElem(test_vec);
    ASSERT_EQ(3, result);
};

// Calculating the number of unique elements in a matrix
TEST(pvalue_check, calculate_number_of_elements_matrix_checked){
    vector<vector<unsigned char>> test_matrix = {{'1', '0', '1', '0', '2', '1'}, {'0', '0', '0', '1', '1', '0'},
                                                 {'2', '1', '0', '0', '2', '1'}};
    int result = PValue::calcNumElem(test_matrix);
    ASSERT_EQ(3, result);
};

// Calculating the number of unique elements in a vector filled with zeros
TEST(pvalue_check, calculate_number_of_elements_vector_all_zeros_checked){
    vector<unsigned char> test_vec = {'0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0'};
    int result = PValue::calcNumElem(test_vec);
    ASSERT_EQ(1, result);
};

// Calculating the number of unique elements in a matrix filled with zeros
TEST(pvalue_check, calculate_number_of_elements_matrix_all_zeros_checked){
    vector<vector<unsigned char>> test_matrix = {{'0', '0', '0', '0', '0', '0'}, {'0', '0', '0', '0', '0', '0'},
                                                 {'0', '0', '0', '0', '0', '0'}};
    int result = PValue::calcNumElem(test_matrix);
    ASSERT_EQ(1, result);
};

// Calculating the number of unique elements in a matrix where all elements are '3'
// The element '3' should not be considered
TEST(pvalue_check, calculate_number_of_elements_matrix_all_threes_checked){
    vector<vector<unsigned char>> test_matrix = {{'3', '3', '3', '3', '3', '3'}, {'3', '3', '3', '3', '3', '3'},
                                                 {'3', '3', '3', '3', '3', '3'}};
    int result = PValue::calcNumElem(test_matrix);
    ASSERT_EQ(0, result);
};

// Calculating the number of unique elements in an empty vector
TEST(pvalue_check, calculate_number_of_elements_empty_vector_check){
    vector<unsigned char> empty = {};
    int result = PValue::calcNumElem(empty);
    ASSERT_EQ(0, result);
};

// Calculating the number of unique elements in an empty matrix
TEST(pvalue_check, calculate_number_of_elements_empty_matrix_check){
    vector<vector<unsigned char>> empty = {};
    int result = PValue::calcNumElem(empty);
    ASSERT_EQ(0, result);
};

// Checking the lower gamma function from the boost library
// The returned value for the given parameters is compared with the value obtained on this resource:
// http://keisan.casio.com/exec/system/1180573447
TEST(pvalue_check, lower_gamma_function_check){
    int p1 = 6;
    double p2 = 8.0;
    ASSERT_NEAR(97.0516, boost::math::tgamma_lower(p1, p2), 0.001);
};

// Testing the function prepareData for test type = "cd"
TEST(pvalue_check, prepare_data_cd_check){
    vector<vector<unsigned char>> cur_G = {};
    vector<unsigned char> cur_A = {};
    InputFile genotype("/Users/EL/CLionProjects/AdjustPValue/test.txt");
    vector<unsigned char> A = {'1', '3', '0', '2', '1', '0'};

    vector<vector<unsigned char>> result_G = {{'0', '0', '0', '1', '1' ,'0'}, {'2', '1', '0', '0', '2', '1'},
                                              {'1', '0', '0', '2', '2', '1'}, {'0', '0', '1', '1', '1', '0'}};
    vector<unsigned char> result_A = {'1', '3', '0', '2', '1', '0'};

    TestsData cur_test = {"cd",2,5};
    PValue::prepareData(cur_G, cur_A, genotype, A, cur_test);

    ASSERT_EQ(result_A, cur_A);
    ASSERT_EQ(result_G, cur_G);
}

// Testing the function prepareData for test type = "r"
TEST(pvalue_check, prepare_data_r_check){
    vector<vector<unsigned char>> cur_G = {};
    vector<unsigned char> cur_A = {};
    InputFile genotype("/Users/EL/CLionProjects/AdjustPValue/test.txt");
    vector<unsigned char> A = {'1', '3', '0', '2', '1', '0'};

    vector<vector<unsigned char>> result_G = {{'0', '0', '0', '0', '0' ,'0'}, {'1', '0', '0', '0', '1', '0'},
                                              {'0', '0', '0', '1', '1', '0'}, {'0', '0', '0', '0', '0', '0'}};
    vector<unsigned char> result_A = {'1', '3', '0', '2', '1', '0'};

    TestsData cur_test = {"r",2,5};
    PValue::prepareData(cur_G, cur_A, genotype, A, cur_test);

    ASSERT_EQ(result_A, cur_A);
    ASSERT_EQ(result_G, cur_G);
}

// Testing the function prepareData for test type = "d"
TEST(pvalue_check, prepare_data_d_check){
    vector<vector<unsigned char>> cur_G = {};
    vector<unsigned char> cur_A = {};
    InputFile genotype("/Users/EL/CLionProjects/AdjustPValue/test.txt");
    vector<unsigned char> A = {'1', '3', '0', '2', '1', '0'};

    vector<vector<unsigned char>> result_G = {{'0', '0', '0', '1', '1' ,'0'}, {'1', '1', '0', '0', '1', '1'},
                                              {'1', '0', '0', '1', '1', '1'}, {'0', '0', '1', '1', '1', '0'}};
    vector<unsigned char> result_A = {'1', '3', '0', '2', '1', '0'};

    TestsData cur_test = {"d",2,5};
    PValue::prepareData(cur_G, cur_A, genotype, A, cur_test);

    ASSERT_EQ(result_A, cur_A);
    ASSERT_EQ(result_G, cur_G);
}

// Testing the function prepareData for test type = "a"
TEST(pvalue_check, prepare_data_a_check){
    vector<vector<unsigned char>> cur_G = {};
    vector<unsigned char> cur_A = {};
    InputFile genotype("/Users/EL/CLionProjects/AdjustPValue/test.txt");
    vector<unsigned char> A = {'1', '3', '0', '2', '1', '0'};

    vector<vector<unsigned char>> result_G = {{'0', '0', '0', '0', '0' ,'0',   '0', '0', '0', '1', '1' ,'0'},
                                              {'1', '0', '0', '0', '1', '0',   '1', '1', '0', '0', '1', '1'},
                                              {'0', '0', '0', '1', '1', '0',   '1', '0', '0', '1', '1', '1'},
                                              {'0', '0', '0', '0', '0', '0',   '0', '0', '1', '1', '1', '0'}};

    vector<unsigned char> result_A = {'1', '3', '0', '2', '1', '0',  '1', '3', '0', '2', '1', '0'};

    TestsData cur_test = {"a",2,5};
    PValue::prepareData(cur_G, cur_A, genotype, A, cur_test);

    ASSERT_EQ(result_A, cur_A);
    ASSERT_EQ(result_G, cur_G);
}

// Testing how a vector of strings is converted to a vector of unsigned chars
TEST(pvalue_check, map_phenotype_values_to_char_check){
    vector<string> phen = {"blue", "red", "green", "red", "black", "blue", "green"};
    vector<unsigned char> new_phen;
    new_phen = PValue::mapPhenotypeValuesToChar(phen);
    vector<unsigned char> result = {'0', '1', '2', '1', '3', '0', '2'};
    ASSERT_EQ(result, new_phen);
}

// Testing how an empty vector of strings is converted to a vector of unsigned chars
TEST(pvalue_check, map_phenotype_values_to_char_empty_vector_check){
    vector<string> phen = {};
    vector<unsigned char> new_phen;
    new_phen = PValue::mapPhenotypeValuesToChar(phen);
    vector<unsigned char> result = {};
    ASSERT_EQ(result, new_phen);
}


