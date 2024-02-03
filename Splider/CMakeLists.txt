CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(Splider)

elements_depends_on_subdirs(Linx)

find_package(GSL)

elements_add_library(
    Splider src/lib/*.cpp
    INCLUDE_DIRS Linx GSL
    LINK_LIBRARIES Linx GSL
    PUBLIC_HEADERS Splider)

elements_add_unit_test(
    BiSpline tests/src/BiSpline_test.cpp 
    EXECUTABLE Splider_BiSpline_test
    LINK_LIBRARIES Splider GSL
    TYPE Boost)
elements_add_unit_test(
    C2 tests/src/C2_test.cpp 
    EXECUTABLE Splider_C2_test
    LINK_LIBRARIES Splider GSL
    TYPE Boost)
elements_add_unit_test(
    Cospline tests/src/Cospline_test.cpp 
    EXECUTABLE Splider_Cospline_test
    LINK_LIBRARIES Splider GSL
    TYPE Boost)
elements_add_unit_test(
    Lagrange tests/src/Lagrange_test.cpp 
    EXECUTABLE Splider_Lagrange_test
    LINK_LIBRARIES Splider GSL
    TYPE Boost)
elements_add_unit_test(
    Spline tests/src/Spline_test.cpp 
    EXECUTABLE Splider_Spline_test
    LINK_LIBRARIES Splider GSL
    TYPE Boost)
elements_add_unit_test(
    Partition tests/src/Partition_test.cpp
    EXECUTABLE Splider_Partition_test
    LINK_LIBRARIES Splider GSL
    TYPE Boost)
                     
                     