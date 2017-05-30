//
// Created by EL on 22.05.17.
//

#include <string>
#include <vector>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../InputFile.h"

using testing::Eq;

/*
namespace{
    class InputData: public testing::Test{
    public:
        InputFile testGenotype("/Users/EL/CLionProjects/AdjustPValue/test.txt"); // "/Users/EL/CLionProjects/AdjustPValue/test.txt"
        InputData(){
            testGenotype();
        };
    };
}

TEST_F(InputData, readPartOfGenotypeMatrix){

}
 */



// How to initialize the object once?
InputFile genotype("/Users/EL/CLionProjects/AdjustPValue/test.txt");

// Reads 3 lines from the middle
TEST(read_file_check, read_part_of_genotype_matrix_from_middle){
    vector<vector<unsigned char>> result = {{'1', '0', '1', '0', '2', '1'}, {'0', '0', '0', '1', '1', '0'}, {'2', '1', '0', '0', '2', '1'}};
    ASSERT_EQ(result, genotype.createGenotypeMatrix(1,3));
}

// Reads 2 lines from the beginning
TEST(read_file_check, read_part_of_genotype_matrix_from_beginning){
    vector<vector<unsigned char>> result = {{'0', '2', '1', '1', '0', '2'}, {'1', '0', '1', '0', '2', '1'}};
    ASSERT_EQ(result, genotype.createGenotypeMatrix(0,1));
}

// Reads 1 line from the end
TEST(read_file_check, read_part_of_genotype_matrix_from_end){
    vector<vector<unsigned char>> result = {{'0', '1', '3', '2', '0', '0'}};
    ASSERT_EQ(result, genotype.createGenotypeMatrix(19,19));
}

// Read non-existent lines
TEST(read_file_check, read_part_of_genotype_matrix_not_existing_lines){
    vector<vector<unsigned char>> result = {};
    ASSERT_EQ(result, genotype.createGenotypeMatrix(25,71));
}

// Part of the lines at the end are non-existent
TEST(read_file_check, read_part_of_genotype_matrix_part_of_lines_non_existent_end){ // ?????????????????????
    InputFile genotype("/Users/EL/CLionProjects/AdjustPValue/test.txt");
    vector<vector<unsigned char>> result = {{'0', '1', '2', '0', '2', '1'}, {'0', '1', '3', '2', '0', '0'}};
    ASSERT_EQ(result, genotype.createGenotypeMatrix(18,28));
}

// Part of the lines at the beginning are non-existent
TEST(read_file_check, read_part_of_genotype_matrix_part_of_lines_non_existent_beg){ // ?????????????????????
    InputFile genotype("/Users/EL/CLionProjects/AdjustPValue/test.txt");
    vector<vector<unsigned char>> result = {{'0', '2', '1', '1', '0', '2'}, {'1', '0', '1', '0', '2', '1'}};
    ASSERT_EQ(result, genotype.createGenotypeMatrix(-3,1));
}

// Add tests for incorrect input and overflow