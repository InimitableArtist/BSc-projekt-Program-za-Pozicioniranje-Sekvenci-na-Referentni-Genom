cmake_minimum_required(VERSION 3.10)

project(sekvence_mapper)

add_executable(mapper mapper.cpp)
add_subdirectory(src)

set(CMAKE_CXX_STANDARD 11)


target_link_libraries(
    mapper
    bioparser
    alignment
    minimizers
)