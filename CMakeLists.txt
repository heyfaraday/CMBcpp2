cmake_minimum_required(VERSION 3.7)

project(CMBcpp)

set(CMAKE_CXX_STANDARD 14)

set (CMAKE_CXX_FLAGS "-Weverything -Ofast")

include_directories(/usr/local/cHealpix/include)

set(CMAKE_EXE_LINKER_FLAGS "-lfftw3l")

set(SOURCE_FILES_main main.cpp)
set(SOURCE_FILES_v_0_0_ src/v_0_0_.cpp)
set(SOURCE_FILES_v_0_1_ src/v_0_1_.cpp)

add_executable(main ${SOURCE_FILES_main})
add_executable(v_0_0_ ${SOURCE_FILES_v_0_0_})
add_executable(v_0_1_ ${SOURCE_FILES_v_0_1_})

include_directories(lib)
add_subdirectory(lib)

target_link_libraries(main lib)
target_link_libraries(v_0_0_ lib)
target_link_libraries(v_0_1_ lib)

include_directories(bin)