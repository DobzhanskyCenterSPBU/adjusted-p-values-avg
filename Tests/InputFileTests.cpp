//
// Created by EL on 22.05.17.
//

#include <string>
#include <vector>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../InputFile.h"

using testing::Eq;

string filename = "test.txt";

// Reads 3 lines from the middle
TEST(read_file_check, read_part_of_genotype_matrix_from_middle){

    InputFile genotype(filename);
    vector<vector<unsigned short>> result = {{1, 0, 1, 0, 2, 1}, {0, 0, 0, 1, 1, 0}, {2, 1, 0, 0, 2, 1}};
    ASSERT_EQ(result, genotype.createGenotypeMatrix(1,3));
}

// Read non-existent lines
TEST(read_file_check, read_part_of_genotype_matrix_not_existing_lines){
    InputFile genotype(filename);
    vector<vector<unsigned short>> result = {};
    ASSERT_EQ(result, genotype.createGenotypeMatrix(25,71));

    result = {{0, 2, 1, 1, 0, 2}, {1, 0, 1, 0, 2, 1}};
    ASSERT_EQ(result, genotype.createGenotypeMatrix(0,1));
}

// Part of the lines at the end are non-existent
TEST(read_file_check, read_part_of_genotype_matrix_part_of_lines_non_existent_end){
    InputFile genotype(filename);
    vector<vector<unsigned short>> result = {{0, 1, 2, 0, 2, 1}, {0, 1, 3, 2, 0, 0}};
    ASSERT_EQ(result, genotype.createGenotypeMatrix(18,28));
}

// Reads 2 lines from the beginning
TEST(read_file_check, read_part_of_genotype_matrix_from_beginning){
    InputFile genotype(filename);
    vector<vector<unsigned short>> result = {{0, 2, 1, 1, 0, 2}, {1, 0, 1, 0, 2, 1}};
    ASSERT_EQ(result, genotype.createGenotypeMatrix(0,1));
}

// Reads 1 line from the end
TEST(read_file_check, read_part_of_genotype_matrix_from_end){
    InputFile genotype(filename);
    vector<vector<unsigned short>> result = {{0, 1, 3, 2, 0, 0}};
    ASSERT_EQ(result, genotype.createGenotypeMatrix(19,19));
}

// Part of the lines at the beginning are non-existent
TEST(read_file_check, read_part_of_genotype_matrix_part_of_lines_non_existent_beg){
    InputFile genotype(filename);
    vector<vector<unsigned short>> result = {{0, 2, 1, 1, 0, 2}, {1, 0, 1, 0, 2, 1}};
    ASSERT_EQ(result, genotype.createGenotypeMatrix(-3,1));
}