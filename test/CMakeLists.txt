find_package(MPI)

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/)

set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

find_package(CGAL)

include(${CGAL_USE_FILE})

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(TestLinearHashMapUtils TestLinearHashMapUtils.cpp ../src/Utils.cpp ../src/Communicator.cpp ../src/GeometricUtils.cpp)

add_test(TestLinearHashMapUtils "${PROJECT_SOURCE_DIR}/test/TestLinearHashMapUtils")

target_link_libraries(TestLinearHashMapUtils gtest_main)

target_link_libraries(TestLinearHashMapUtils ${CGAL_LIBS})

target_link_libraries(TestLinearHashMapUtils ${MPI_C_LIBRARIES})
