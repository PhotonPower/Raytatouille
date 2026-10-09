add_library(rtt_optim src/bounds.cpp src/levenberg_marquardt.cpp src/merit.cpp src/optimize.cpp
  src/variables.cpp)
add_library(rtt::optim ALIAS rtt_optim)
target_include_directories(rtt_optim PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
# Layer "workflows" (ADR 0030, point 13): rtt-trace for RunControl and RunMonitor (ADR 0004
# addendum #83), rtt-model and rtt-compile (libraries, diagnostics) in the public headers;
# the operands call rtt-paraxial and rtt-analysis. Not rtt-io: the patch is written here.
# No Eigen: the n x n systems are solved by a single-threaded Cholesky factorization in index
# order (MNT Algorithm A.4), so results do not depend on threads (ADR 0030, point 11).
target_link_libraries(rtt_optim PUBLIC rtt::trace rtt::model rtt::compile
  PRIVATE rtt::analysis rtt::paraxial rtt::diagnostics TBB::tbb)
rtt_target_defaults(rtt_optim LIBRARY)

if(RTT_BUILD_TESTS)
  # TBB limits the threads in the determinism test; rtt::io and nlohmann load the reference
  # systems, read the edit form and apply the result patch in the tests only.
  rtt_add_test(rtt_optim_tests
    SOURCES tests/test_bounds.cpp tests/test_levenberg_marquardt.cpp tests/test_optimize.cpp
      tests/test_variables.cpp
    LIBS rtt::optim rtt::io rtt::analysis rtt::paraxial TBB::tbb nlohmann_json::nlohmann_json)
endif()
