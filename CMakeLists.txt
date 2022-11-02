cmake_minimum_required(VERSION 3.0)
project(SubgraphSampler)

set(CMAKE_C_COMPILER /opt/homebrew/bin/gcc-12)
set(CMAKE_CXX_COMPILER /opt/homebrew/bin/g++-12)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -Ofast -pthread")
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR})

file(GLOB SOURCES src/*.cc)

add_subdirectory(main)