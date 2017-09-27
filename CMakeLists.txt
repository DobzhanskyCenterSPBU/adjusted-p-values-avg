cmake_minimum_required(VERSION 3.5)
project(AdjustPValue)

set(CMAKE_CXX_STANDARD 14)
add_subdirectory(lib/googletest-master)
include_directories(lib/googletest-master/googletest/include)
include_directories(lib/googletest-master/googlemock/include)

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fext-numeric-literals -fopenmp -lmysqlcppconn -lmysqlclient")

# added mysql-connector-c++
set(FULL_PATH_TO_MYSQL_CONNECTOR_CPP_DIR /usr)
include_directories(/usr/include)
include_directories(/usr/include/cppconn)
include_directories(/usr/share/doc)
link_directories(/usr/lib)
#link_libraries(/usr/lib/libmysqlcppconn.a)
link_libraries(/usr/lib/libmysqlcppconn.so)
link_libraries(/usr/lib/libmysqlcppconn-static.a)

set(SOURCE_FILES main.cpp InputFile.h InputData.h InputFile.cpp InputFile.h PValue.cpp PValue.h Tests/InputFileTests.cpp Tests/PValueTests.cpp PerformanceTests.cpp PerformanceTests.h InputDataBase.cpp InputDataBase.h)
add_executable(AdjustPValue ${SOURCE_FILES} InputFile.cpp InputFile.h PValue.cpp PValue.h Tests/InputFileTests.cpp Tests/PValueTests.cpp PerformanceTests.cpp PerformanceTests.h InputDataBase.cpp InputDataBase.h)
target_link_libraries(AdjustPValue gtest gtest_main mysqlcppconn)
