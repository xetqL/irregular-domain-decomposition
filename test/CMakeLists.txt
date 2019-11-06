add_executable(TestA test1.cpp)

target_link_libraries(TestA gtest_main)

add_test(TestA "${PROJECT_SOURCE_DIR}/test/TestA")