add_library(rtt_math src/isometry.cpp)
add_library(rtt::math ALIAS rtt_math)
target_include_directories(rtt_math PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_math PUBLIC Eigen3::Eigen)
rtt_target_defaults(rtt_math LIBRARY)

if(RTT_BUILD_TESTS)
  rtt_add_test(rtt_math_tests
    SOURCES tests/test_units.cpp tests/test_isometry.cpp
    LIBS rtt::math)
endif()
