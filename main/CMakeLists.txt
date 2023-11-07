add_compile_options(-Wno-deprecated-declarations)
FIND_PACKAGE( Boost )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O2")
include_directories(${Boost_INCLUDE_DIRS})

add_executable(FaSTest main.cc ${SOURCES})
target_link_libraries(FaSTest ${Boost_LIBRARIES})
