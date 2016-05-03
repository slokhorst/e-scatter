add_executable(unittests main.cc)
target_link_libraries(unittests gtest gtest_main)
add_test(NAME main-unit-test COMMAND unittests)

