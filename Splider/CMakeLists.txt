CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(Splider)

elements_depends_on_subdirs(LinxCore)

elements_add_library(Splider src/lib/*.cpp
                     INCLUDE_DIRS LinxCore
                     LINK_LIBRARIES LinxCore
                     PUBLIC_HEADERS Splider)

elements_add_unit_test(Spline tests/src/Spline_test.cpp 
                     EXECUTABLE Splider_Spline_test
                     LINK_LIBRARIES Splider
                     TYPE Boost)

