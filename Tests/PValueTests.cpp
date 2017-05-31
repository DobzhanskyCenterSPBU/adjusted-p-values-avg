//
// Created by EL on 23.05.17.
//

#include <string>
#include <vector>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../PValue.h"
#include "../InputFile.h"


vector<unsigned char> test_vec = {'0', '2', '1', '1', '0', '2', '1', '0', '1', '0', '2', '1'};
vector<vector<unsigned char>> test_matrix = {{'1', '0', '1', '0', '2', '1'}, {'0', '0', '0', '1', '1', '0'}, {'2', '1', '0', '0', '2', '1'}};


TEST(pvalue_check, calculate_number_of_elements_vector_checked){
    int result = PValue::calcNumElem(test_vec);
    ASSERT_EQ(3, result);
};

TEST(pvalue_check, calculate_number_of_elements_vector_all_zeros_checked){
    vector<unsigned char> test_vec = {'0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0'};
    int result = PValue::calcNumElem(test_vec);
    ASSERT_EQ(1, result);
};

TEST(pvalue_check, calculate_number_of_elements_matrix_all_zeros_checked){
    vector<vector<unsigned char>> test_matrix = {{'0', '0', '0', '0', '0', '0'}, {'0', '0', '0', '0', '0', '0'}, {'0', '0', '0', '0', '0', '0'}};
    int result = PValue::calcNumElem(test_matrix);
    ASSERT_EQ(1, result);
};

TEST(pvalue_check, calculate_number_of_elements_matrix_all_threes_checked){
    vector<vector<unsigned char>> test_matrix = {{'3', '3', '3', '3', '3', '3'}, {'3', '3', '3', '3', '3', '3'}, {'3', '3', '3', '3', '3', '3'}};
    int result = PValue::calcNumElem(test_matrix);
    ASSERT_EQ(0, result);
};

TEST(pvalue_check, calculate_number_of_elements_matrix_checked){
    int result = PValue::calcNumElem(test_matrix);
    ASSERT_EQ(3, result);
};

TEST(pvalue_check, calculate_number_of_elements_empty_vector_check){
    vector<unsigned char> empty = {};
    int result = PValue::calcNumElem(empty);
    ASSERT_EQ(0, result);
};

TEST(pvalue_check, calculate_number_of_elements_empty_matrix_check){
    vector<vector<unsigned char>> empty = {};
    int result = PValue::calcNumElem(empty);
    ASSERT_EQ(0, result);
};

TEST(pvalue_check, lower_gamma_function_check){ // Add comment + link
    int p1 = 6;
    double p2 = 8.0;
    ASSERT_LT(boost::math::tgamma_lower(p1, p2) - 97.0516, 0.001); // ASSERT_DOUBLE_EQ!
};

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


TEST(pvalue_check, prepare_data_a_second_check){
    vector<vector<unsigned char>> cur_G = {};
    vector<unsigned char> cur_A = {};
    InputFile genotype("/Users/EL/CLionProjects/AdjustPValue/test.txt");
    vector<unsigned char> A = {'1', '3', '0', '2', '1', '0'};

    vector<vector<unsigned char>> result_G = {{'0', '1', '0', '0', '0' ,'0',   '0', '1', '1', '1', '0' ,'0'},
                                              {'0', '1', '0', '0', '1', '0',   '0', '1', '1', '0', '1', '1'},
                                              {'0', '0', '0', '1', '0', '0',   '0', '1', '1', '1', '1', '1'}};

    vector<unsigned char> result_A = {'1', '3', '0', '2', '1', '0',  '1', '3', '0', '2', '1', '0'};

    TestsData cur_test = {"a",12,14};
    PValue::prepareData(cur_G, cur_A, genotype, A, cur_test);

    ASSERT_EQ(result_A, cur_A);
    ASSERT_EQ(result_G, cur_G);
}


