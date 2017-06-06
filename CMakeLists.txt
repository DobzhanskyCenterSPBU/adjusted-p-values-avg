cmake_minimum_required(VERSION 3.7)
project(AdjustPValue)

set(CMAKE_CXX_STANDARD 14)
add_subdirectory(lib/googletest-master)
include_directories(lib/googletest-master/googletest/include)
include_directories(lib/googletest-master/googlemock/include)

set(SOURCE_FILES main.cpp InputFile.h InputData.h InputFile.cpp InputFile.h PValue.cpp PValue.h Tests/InputFileTests.cpp Tests/PValueTests.cpp PerformanceTests.cpp PerformanceTests.h)
add_executable(AdjustPValue ${SOURCE_FILES} InputFile.cpp InputFile.h PValue.cpp PValue.h Tests/InputFileTests.cpp Tests/PValueTests.cpp PerformanceTests.cpp PerformanceTests.h)
target_link_libraries(AdjustPValue gtest gtest_main)