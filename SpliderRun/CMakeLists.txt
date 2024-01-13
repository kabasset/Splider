CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(SpliderRun)

elements_depends_on_subdirs(Splider)

find_package(GSL)

elements_add_library(SpliderRun src/lib/*.cpp
                     INCLUDE_DIRS Splider
                     LINK_LIBRARIES Splider
                     PUBLIC_HEADERS SpliderRun)

elements_add_executable(SpliderBenchmark src/program/SpliderBenchmark.cpp
                     INCLUDE_DIRS SpliderRun
                     LINK_LIBRARIES SpliderRun)

elements_add_unit_test(SpliderDemo tests/src/SpliderDemo_test.cpp 
                     EXECUTABLE SpliderRun_SpliderDemo_test
                     LINK_LIBRARIES SpliderRun GSL
                     TYPE Boost)
