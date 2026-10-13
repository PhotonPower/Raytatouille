add_library(rtt_tolerance src/random.cpp)
add_library(rtt::tolerance ALIAS rtt_tolerance)
target_include_directories(rtt_tolerance PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
# Layer "workflows" (docs/architecture.md, ADR 0036). For now only the random perturbations
# (#259), which need the tolerance items of rtt-model; the trial evaluation (#257) adds
# rtt-optim, rtt-analysis, rtt-trace and rtt-compile.
target_link_libraries(rtt_tolerance PUBLIC rtt::model)
rtt_target_defaults(rtt_tolerance LIBRARY)

if(RTT_BUILD_TESTS)
  # TBB runs the draws on several threads in the determinism test.
  rtt_add_test(rtt_tolerance_tests
    SOURCES tests/test_random.cpp
    LIBS rtt::tolerance TBB::tbb)
  # src/: the kernels of the normal distribution (random_detail.hpp) are tested directly.
  target_include_directories(rtt_tolerance_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
endif()
