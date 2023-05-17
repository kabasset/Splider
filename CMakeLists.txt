CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

find_package(ElementsProject)

set(USE_SPHINX OFF
    CACHE STRING "Disable Sphinx API documentation."
    FORCE)

set(CMAKE_CXX_FLAGS_RELEASE "-flto -O3"
    CACHE STRING "Enable -o3 for release builds."
    FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O3"
    CACHE STRING "Enable -o3 for release builds."
    FORCE)

elements_project(Splider 1.0 USE Linx 1.0)

