# Vertical remapping algorithm.
# Modifying Author: Jason Torchisnky
# Date: Summer 2021

cmake_minimum_required(VERSION 3.10)

enable_language(Fortran)

project(vert_remap_redux
        VERSION 1.0
        DESCRIPTION "Modification of the vertical remapping algorithm for HOMME"
       )

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory(src)
