cmake_minimum_required(VERSION 3.6)

project(lib2)

set(SOURCE_LIB_2 Alm_Base.cpp Cl_Base.cpp Map_Base.cpp aml.cpp constants.cpp parameters.cpp)
set(HEADER_LIB_2 Alm_Base.hpp Cl_Base.hpp Map_Base.hpp aml.hpp constants.hpp parameters.hpp)

add_library(lib2 STATIC ${SOURCE_LIB_2} ${HEADER_LIB_2})
