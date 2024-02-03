CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(SpliderRun)

elements_depends_on_subdirs(Splider)
elements_depends_on_subdirs(ElementsKernel)

find_package(Boost COMPONENTS program_options unit_test_framework REQUIRED)
find_package(GSL)

elements_add_executable(SpliderBenchmark src/program/SpliderBenchmark.cpp
                     INCLUDE_DIRS Splider Boost
                     LINK_LIBRARIES Splider Boost)
elements_add_executable(SpliderBiBenchmark src/program/SpliderBiBenchmark.cpp
                     INCLUDE_DIRS Splider Boost
                     LINK_LIBRARIES Splider Boost)
elements_add_executable(SpliderSin src/program/SpliderSin.cpp
                     INCLUDE_DIRS Splider Boost
                     LINK_LIBRARIES Splider Boost)

elements_add_unit_test(SpliderDemo tests/src/SpliderDemo_test.cpp 
                     EXECUTABLE SpliderRun_SpliderDemo_test
                     LINK_LIBRARIES Splider Boost GSL ElementsKernel
                     TYPE Boost)
