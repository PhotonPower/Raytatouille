add_library(rtt_trace src/ray_batch.cpp)
add_library(rtt::trace ALIAS rtt_trace)
target_include_directories(rtt_trace PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_trace PUBLIC rtt::compile rtt::math)
rtt_target_defaults(rtt_trace LIBRARY)

if(RTT_BUILD_TESTS)
  rtt_add_test(rtt_trace_tests
    SOURCES tests/test_ray_batch.cpp
    LIBS rtt::trace)
endif()
