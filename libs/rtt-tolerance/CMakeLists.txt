add_library(rtt_tolerance src/classes.cpp src/expand.cpp src/monte_carlo.cpp src/random.cpp
  src/sensitivity.cpp src/trial.cpp)
add_library(rtt::tolerance ALIAS rtt_tolerance)
target_include_directories(rtt_tolerance PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
# Layer "workflows" (docs/architecture.md, ADR 0036). The random perturbations (#259) need the
# tolerance items of rtt-model; expand (#260) validates the system and throws
# compile::CompileError like optimize. The trial evaluation (#257) needs rtt-optim (merit
# function and optimize with the trial overlay, #256) and rtt-trace (RunControl) in its header,
# rtt-diagnostics for the warning codes. The sensitivity (#258) and the Monte Carlo (#259) run
# their trials in parallel with TBB.
target_link_libraries(rtt_tolerance PUBLIC rtt::model rtt::compile rtt::optim rtt::trace
  PRIVATE rtt::diagnostics TBB::tbb)
rtt_target_defaults(rtt_tolerance LIBRARY)

if(RTT_BUILD_TESTS)
  # TBB runs the draws on several threads in the determinism tests; rtt-io reads the reference
  # files of the class tests.
  rtt_add_test(rtt_tolerance_tests
    SOURCES tests/test_classes.cpp tests/test_m7_acceptance.cpp tests/test_monte_carlo.cpp
      tests/test_object_reference.cpp tests/test_radius_sign.cpp tests/test_random.cpp
      tests/test_sensitivity.cpp tests/test_trial.cpp
    LIBS rtt::tolerance rtt::io TBB::tbb)
  # src/: the kernels of the normal distribution (random_detail.hpp), the RSS reduction
  # (sensitivity_detail.hpp) and the statistics (monte_carlo_detail.hpp) are tested directly.
  target_include_directories(rtt_tolerance_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  # Glass catalogues of the reference systems: nu_d for the class limits of abbe (#260).
  target_compile_definitions(rtt_tolerance_tests PRIVATE
    RTT_CATALOG_DIR="${PROJECT_SOURCE_DIR}/tests/catalogs")
endif()
