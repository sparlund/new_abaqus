cmake_minimum_required(VERSION 3.13)

project(new_abaqus)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_subdirectory(src/src build)
include_directories(src/include)
add_subdirectory(tests)
set_property(GLOBAL PROPERTY LINKER_LANGUAGE CXX)
