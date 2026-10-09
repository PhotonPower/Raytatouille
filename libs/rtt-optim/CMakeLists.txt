add_library(rtt_optim src/bounds.cpp src/levenberg_marquardt.cpp src/variables.cpp)
add_library(rtt::optim ALIAS rtt_optim)
target_include_directories(rtt_optim PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
# rtt-trace for RunControl and RunMonitor (cancellation and progress, ADR 0004 addendum #83);
# rtt-model for the variables of a run (ADR 0030, point 5).
# No Eigen: the n x n systems are solved by a single-threaded Cholesky factorization in index
# order (MNT Algorithm A.4), so results do not depend on threads (ADR 0030, point 11).
target_link_libraries(rtt_optim PUBLIC rtt::trace rtt::model PRIVATE TBB::tbb)
rtt_target_defaults(rtt_optim LIBRARY)

if(RTT_BUILD_TESTS)
  # TBB limits the threads in the determinism test; rtt::io and nlohmann load the reference
  # systems and read the edit form in the tests only.
  rtt_add_test(rtt_optim_tests
    SOURCES tests/test_bounds.cpp tests/test_levenberg_marquardt.cpp tests/test_variables.cpp
    LIBS rtt::optim rtt::io TBB::tbb nlohmann_json::nlohmann_json)
endif()
