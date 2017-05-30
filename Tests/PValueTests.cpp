//
// Created by EL on 23.05.17.
//

#include <string>
#include <vector>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../PValue.h"


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



