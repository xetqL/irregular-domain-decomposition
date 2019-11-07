cmake_minimum_required(VERSION 3.12)

project(adlbirreg)

set(CMAKE_CXX_STANDARD 14)
include_directories(/home/xetql/eigen-eigen-323c052e1731)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        release-1.8.0
)
FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif()

enable_testing()

include_directories(include)
add_subdirectory(extern/cgal)
add_subdirectory(test)
add_subdirectory(src)