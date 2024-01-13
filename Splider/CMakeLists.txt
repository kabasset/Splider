CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(Splider)

elements_depends_on_subdirs(Linx ElementsKernel)

find_package(GSL)

elements_add_library(Splider src/lib/*.cpp
                     INCLUDE_DIRS Linx GSL ElementsKernel
                     LINK_LIBRARIES Linx GSL ElementsKernel
                     PUBLIC_HEADERS Splider)

elements_add_unit_test(Spline tests/src/Spline_test.cpp 
                     EXECUTABLE Splider_Spline_test
                     LINK_LIBRARIES Splider GSL
                     TYPE Boost)
elements_add_unit_test(BiSpline tests/src/BiSpline_test.cpp 
                     EXECUTABLE Splider_BiSpline_test
                     LINK_LIBRARIES Splider GSL
                     TYPE Boost)

