add_compile_options(-Ofast)
#add_compile_options(-Wall)
add_compile_options(-std=c++20)
#add_compile_options(-DNDEBUG)

# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)

add_executable(SubgraphSampler main.cc ${SOURCES})
