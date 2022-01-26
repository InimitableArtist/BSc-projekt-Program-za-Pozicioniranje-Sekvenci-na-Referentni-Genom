cmake_minimum_required(VERSION 3.10)

project(sekvence_mapper)

add_executable(mapper mapper.cpp)
add_subdirectory(src)

set(CMAKE_CXX_STANDARD 11)
include(CTest)

target_link_libraries(
    mapper
    bioparser
    alignment
    minimizers
    thread_pool
)
option(align_tests "Align tests" ON)
if (align_tests)
    enable_testing()
    find_package(GTest REQUIRED)
    add_executable(align_test tests/alignment_tests.cpp)
    target_link_libraries(align_test PUBLIC alignment GTest::Main)
    add_test(NAME align_test
            COMMAND align_test)
endif()

if (NOT TARGET thread_pool)
    add_subdirectory(thread_pool EXCLUDE_FROM_ALL)
endif()