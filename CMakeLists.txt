cmake_minimum_required(VERSION 3.7)
project(AdjustPValue)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp InputFile.h InputData.h InputFile.cpp InputFile.h PValue.cpp PValue.h)
add_executable(AdjustPValue ${SOURCE_FILES} InputFile.cpp InputFile.h PValue.cpp PValue.h)