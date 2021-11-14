add_library(alignment alignment.h alignment.cpp)
add_library(minimizers minimizers.h minimizers.cpp)

if (NOT TARGET bioparser)
    add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif ()