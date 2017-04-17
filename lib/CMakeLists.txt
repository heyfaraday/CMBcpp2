cmake_minimum_required(VERSION 3.6)

project(lib)

set(SOURCE_LIB aml.cpp constants.cpp distance.cpp fft.cpp functionals.cpp io.cpp parameters.cpp
        pml.cpp utils.cpp)
set(HEADER_LIB aml.hpp constants.hpp distance.hpp fft.hpp functionals.hpp io.hpp parameters.hpp
        pml.hpp utils.hpp)

add_library(lib STATIC ${SOURCE_LIB} ${HEADER_LIB})