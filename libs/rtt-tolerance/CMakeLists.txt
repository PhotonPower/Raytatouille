add_library(rtt_tolerance src/classes.cpp src/expand.cpp src/random.cpp)
add_library(rtt::tolerance ALIAS rtt_tolerance)
target_include_directories(rtt_tolerance PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
# Layer "workflows" (docs/architecture.md, ADR 0036). The random perturbations (#259) need the
# tolerance items of rtt-model; expand (#260) validates the system and throws
# compile::CompileError like optimize. The trial evaluation (#257) adds rtt-optim, rtt-analysis
# and rtt-trace.
target_link_libraries(rtt_tolerance PUBLIC rtt::model rtt::compile)
rtt_target_defaults(rtt_tolerance LIBRARY)

if(RTT_BUILD_TESTS)
  # TBB runs the draws on several threads in the determinism tests; rtt-io reads the reference
  # files of the class tests.
  rtt_add_test(rtt_tolerance_tests
    SOURCES tests/test_classes.cpp tests/test_random.cpp
    LIBS rtt::tolerance rtt::io TBB::tbb)
  # src/: the kernels of the normal distribution (random_detail.hpp) are tested directly.
  target_include_directories(rtt_tolerance_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
endif()
