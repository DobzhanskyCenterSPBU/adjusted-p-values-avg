//
// Created by EL on 23.05.17.
//

#include <string>
#include <vector>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../PValue.h"


vector<unsigned char> test_vec = {'0', '2', '1', '1', '0', '2', '1', '0', '1', '0', '2', '1'};


TEST(pvalue_check, calculate_number_of_elements_no_extra_symbol_checked){
    int result = PValue::calcNumElem(test_vec.begin(), test_vec.end());
    ASSERT_EQ(3, result);
};

TEST(pvalue_check, calculate_number_of_elements_extra_symbol_checked){

    int result = PValue::calcNumElem(test_vec.begin(), test_vec.end(), '2');
    ASSERT_EQ(2, result);

    result = PValue::calcNumElem(test_vec.begin(), test_vec.end(), '3');
    ASSERT_EQ(3, result);

    vector<unsigned char> test_vec2 = {'2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2', '2'};

    result = PValue::calcNumElem(test_vec2.begin(), test_vec2.end(), '2');
    ASSERT_EQ(0, result);

    result = PValue::calcNumElem(test_vec2.begin(), test_vec2.end(), '3');
    ASSERT_EQ(1, result);
};

TEST(pvalue_check, calculate_number_of_elements_empty_vector_check){
    vector<unsigned char> empty = {};
    int result = PValue::calcNumElem(empty.begin(), empty.end());
    ASSERT_EQ(0, result);

    result = PValue::calcNumElem(empty.begin(), empty.end(), '5');
    ASSERT_EQ(0, result);
};



