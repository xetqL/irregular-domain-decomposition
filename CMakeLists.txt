cmake_minimum_required(VERSION 3.13)

project(adlbirreg)

set(CMAKE_CXX_STANDARD 14)

include_directories(/home/xetql/eigen-eigen-323c052e1731)

find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} )

find_package(CGAL REQUIRED)

include(${CGAL_USE_FILE})

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(adlbirreg main.cpp Communicator.hpp Partition.cpp Partition.hpp GeometricUtils.hpp Domain.cpp Domain.hpp Cell.cpp Cell.hpp Utils.hpp GeometricUtils.cpp Communicator.cpp)

add_executable(test_communcation comm.cpp Communicator.hpp)

target_link_libraries(adlbirreg ${CGAL_LIBS})

target_link_libraries(adlbirreg ${MPI_C_LIBRARIES})

target_link_libraries(test_communcation ${MPI_C_LIBRARIES})
