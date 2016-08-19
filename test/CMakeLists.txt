include_directories(${GTEST_INCLUDE_DIRS})

add_executable(unittests main.cc)
target_link_libraries(unittests ${GTEST_BOTH_LIBRARIES})
add_test(NAME main-unit-test COMMAND unittests)
