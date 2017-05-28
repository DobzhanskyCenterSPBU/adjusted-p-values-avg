//
// Created by EL on 22.05.17.
//

#include <string>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../InputFile.h"

using testing::Eq;

namespace{
    class ClassDeclaration: public testing::Test{
    public:
        InputFile genotype("/Users/EL/CLionProjects/AdjustPValue/test.txt");
        ClassDeclaration(){
            genotype;
        };
    };
}

TEST_F(ClassDeclaration, readPartOfGenotypeMatrix){

}