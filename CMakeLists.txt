cmake_minimum_required(VERSION 3.0)
project(SubgraphSampler)

add_compile_options(-Ofast)
#add_compile_options(-Wall)
add_compile_options(-std=c++2a)

#add_compile_options(-fsanitize=address)
#add_link_options(-fsanitize=address)
FIND_PACKAGE( Boost )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR})

file(GLOB SOURCES src/*.cc)
add_subdirectory(main)

#LD_PRELOAD=/usr/local/lib/libprofiler.dylib CPUPROFILE=SubgraphSampler.prof CPUPROFILESIGNAL=12 ./main/SubgraphSampler
#pprof --lines --pdf ./main/SubgraphSampler SubgraphSampler.prof > SubgraphSampler.pdf