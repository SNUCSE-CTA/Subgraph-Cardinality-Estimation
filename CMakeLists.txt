cmake_minimum_required(VERSION 3.0)
project(SubgraphSampler)

FIND_PACKAGE( Boost )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -Ofast")
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR})

file(GLOB SOURCES src/*.cc)
add_subdirectory(main)