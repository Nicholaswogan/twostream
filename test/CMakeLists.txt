add_library(twostream_test_utils twostream_test_utils.f90)

add_executable(benchmark_twostream benchmark_twostream.f90)
target_link_libraries(benchmark_twostream twostream twostream_test_utils)

add_executable(test_twostream test_twostream.f90)
target_link_libraries(test_twostream twostream twostream_test_utils)

add_executable(test_twostream_reproduce test_twostream_reproduce.f90)
target_link_libraries(test_twostream_reproduce twostream twostream_test_utils)