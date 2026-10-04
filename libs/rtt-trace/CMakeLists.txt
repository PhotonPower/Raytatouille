add_library(rtt_trace src/ray_batch.cpp src/apply_event.cpp src/sequential.cpp)
add_library(rtt::trace ALIAS rtt_trace)
target_include_directories(rtt_trace PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_trace PUBLIC rtt::compile rtt::math PRIVATE TBB::tbb)
rtt_target_defaults(rtt_trace LIBRARY)

if(RTT_BUILD_TESTS)
  # rtt::io only loads the reference systems in the tests; TBB limits the threads in the
  # determinism test.
  rtt_add_test(rtt_trace_tests
    SOURCES tests/test_ray_batch.cpp tests/test_apply_event.cpp tests/test_sequential.cpp
    LIBS rtt::trace rtt::io TBB::tbb)
endif()
