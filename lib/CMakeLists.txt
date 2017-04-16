cmake_minimum_required(VERSION 3.6)

project(lib)

set(SOURCE_LIB alm.cpp constants.cpp fft.cpp functionals.cpp io.cpp parameters.cpp utils.cpp)
set(HEADER_LIB alm.hpp constants.hpp fft.hpp functionals.hpp io.hpp parameters.hpp utils.hpp)

add_library(lib STATIC ${SOURCE_LIB} ${HEADER_LIB})