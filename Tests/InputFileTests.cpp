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
    vector<unsigned char> result = {'1', '0', '1', '0', '2', '1', '0', '0', '0', '1', '1', '0', '2', '1', '0', '0', '2', '1'};
    ASSERT_EQ(result, genotype.createGenotypeVector(1,3));
}

// Reads 2 lines from the beginning
TEST(read_file_check, read_part_of_genotype_matrix_from_beginning){
    vector<unsigned char> result = {'0', '2', '1', '1', '0', '2', '1', '0', '1', '0', '2', '1'};
    ASSERT_EQ(result, genotype.createGenotypeVector(0,1));
}

// Reads 1 line from the end
TEST(read_file_check, read_part_of_genotype_matrix_from_end){
    vector<unsigned char> result = {'0', '1', '3', '2', '0', '0'};
    ASSERT_EQ(result, genotype.createGenotypeVector(19,19));
}

// Add tests for incorrect input and overflow