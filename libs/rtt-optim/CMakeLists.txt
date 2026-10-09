add_library(rtt_optim src/bounds.cpp src/levenberg_marquardt.cpp)
add_library(rtt::optim ALIAS rtt_optim)
target_include_directories(rtt_optim PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
# rtt-trace only for RunControl and RunMonitor (cancellation and progress, ADR 0004 addendum #83).
# No Eigen: the n x n systems are solved by a single-threaded Cholesky factorization in index
# order (MNT Algorithm A.4), so results do not depend on threads (ADR 0030, point 11).
target_link_libraries(rtt_optim PUBLIC rtt::trace PRIVATE TBB::tbb)
rtt_target_defaults(rtt_optim LIBRARY)

if(RTT_BUILD_TESTS)
  # TBB limits the threads in the determinism test.
  rtt_add_test(rtt_optim_tests
    SOURCES tests/test_bounds.cpp tests/test_levenberg_marquardt.cpp
    LIBS rtt::optim TBB::tbb)
endif()
