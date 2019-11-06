find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} )

find_package(CGAL REQUIRED)

include(${CGAL_USE_FILE})

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(adlbirreg main.cpp  Partition.cpp   Domain.cpp  Cell.cpp  GeometricUtils.cpp Communicator.cpp Utils.cpp DiffusiveOptimizer.cpp)

add_executable(test_communcation comm.cpp Communicator.cpp)

target_link_libraries(adlbirreg ${CGAL_LIBS})

target_link_libraries(adlbirreg ${MPI_C_LIBRARIES})

target_link_libraries(test_communcation ${MPI_C_LIBRARIES})