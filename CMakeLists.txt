cmake_minimum_required(VERSION 3.13)

project(new_abaqus)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
add_subdirectory(./src)
include_directories(./src/include)
add_subdirectory(tests)
set_property(GLOBAL PROPERTY LINKER_LANGUAGE CXX)
