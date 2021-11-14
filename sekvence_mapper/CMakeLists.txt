cmake_minimum_required(VERSION 3.10)

project(sekvence_mapper)

add_executable(mapper mapper.cpp)
add_subdirectory(src)

target_link_libraries(
    mapper
    bioparser
    alignment
    minimizers
)