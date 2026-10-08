add_library(rtt_trace src/ray_batch.cpp src/apply_event.cpp src/sequential.cpp src/sources.cpp
  src/run_control.cpp)
add_library(rtt::trace ALIAS rtt_trace)
target_include_directories(rtt_trace PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
# rtt-paraxial provides pupils and the paraxial chief ray for aiming (ADR 0016).
# rtt-coating: layer type in EventMedia; rtt-polar: interface PRT matrices (ADR 0021).
target_link_libraries(rtt_trace PUBLIC rtt::compile rtt::paraxial rtt::coating rtt::math
  PRIVATE rtt::polar TBB::tbb)
rtt_target_defaults(rtt_trace LIBRARY)

if(RTT_BUILD_TESTS)
  # rtt::io only loads the reference systems in the tests; TBB limits the threads in the
  # determinism test.
  rtt_add_test(rtt_trace_tests
    SOURCES tests/test_ray_batch.cpp tests/test_apply_event.cpp tests/test_sequential.cpp
            tests/test_sources.cpp tests/test_interactions.cpp tests/test_ray_paths.cpp
            tests/test_run_control.cpp
    LIBS rtt::trace rtt::io rtt::polar rtt::coating TBB::tbb)
  target_compile_definitions(rtt_trace_tests PRIVATE
    RTT_CATALOG_DIR="${PROJECT_SOURCE_DIR}/tests/catalogs")
endif()
