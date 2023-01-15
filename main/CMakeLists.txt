add_compile_options(-Ofast)
#add_compile_options(-Wall)
add_compile_options(-std=c++2a)
#add_compile_options(-DNDEBUG)

#add_compile_options(-fsanitize=address)
#add_link_options(-fsanitize=address)
FIND_PACKAGE( Boost )
<<<<<<< HEAD
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a -Ofast -pthread")
=======
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a -Ofast -pthread -g")
>>>>>>> 44c064df666959375baa6dbdedb954e4483a9973
include_directories(${Boost_INCLUDE_DIRS})

add_executable(SubgraphSampler main.cc ${SOURCES})
target_link_libraries(SubgraphSampler ${Boost_LIBRARIES})
