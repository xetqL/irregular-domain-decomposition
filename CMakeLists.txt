cmake_minimum_required(VERSION 3.14)
project(adlbirreg)

set(CMAKE_CXX_STANDARD 14)
include_directories(/home/xetql/eigen-eigen-323c052e1731)
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})
add_executable(adlbirreg main.cpp)
target_link_libraries(adlbirreg ${CGAL_LIBS})
