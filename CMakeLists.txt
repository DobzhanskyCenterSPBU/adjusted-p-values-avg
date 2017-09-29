cmake_minimum_required(VERSION 3.7)
project(AdjustPValue)

set(CMAKE_CXX_STANDARD 14)
add_subdirectory(lib/googletest-master)
include_directories(lib/googletest-master/googletest/include)
include_directories(lib/googletest-master/googlemock/include)

# added -fopenmp
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fext-numeric-literals -fopenmp")

# added mysql-connector-c++
set(FULL_PATH_TO_MYSQL_CONNECTOR_CPP_DIR /usr/local/Cellar/mysql-connector-c++/1.1.9)
include_directories(/usr/local/Cellar/mysql-connector-c++/1.1.9/include)
include_directories(/usr/local/Cellar/mysql-connector-c++/1.1.9/include/cppconn)
link_directories(/usr/local/Cellar/mysql-connector-c++/1.1.9/lib)
link_libraries(/usr/local/Cellar/mysql-connector-c++/1.1.9/lib/libmysqlcppconn.dylib)

set(SOURCE_FILES main.cpp InputFile.h InputData.h InputFile.cpp InputFile.h PValue.cpp PValue.h Tests/InputFileTests.cpp Tests/PValueTests.cpp PerformanceTests.cpp PerformanceTests.h InputDataBase.cpp InputDataBase.h Tests/InputDataBaseTests.cpp)
add_executable(AdjustPValue ${SOURCE_FILES} InputFile.cpp InputFile.h PValue.cpp PValue.h Tests/InputFileTests.cpp Tests/PValueTests.cpp PerformanceTests.cpp PerformanceTests.h InputDataBase.cpp InputDataBase.h)
target_link_libraries(AdjustPValue gtest gtest_main)