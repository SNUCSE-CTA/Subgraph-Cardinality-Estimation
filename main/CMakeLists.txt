add_compile_options(-Ofast)
#add_compile_options(-Wall)
add_compile_options(-std=c++2a)
#add_compile_options(-DNDEBUG)

# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)
FIND_PACKAGE( Boost )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a -Ofast -pthread")
include_directories(${Boost_INCLUDE_DIRS})

add_executable(SubgraphSampler main.cc ${SOURCES})
target_link_libraries(SubgraphSampler ${Boost_LIBRARIES})
