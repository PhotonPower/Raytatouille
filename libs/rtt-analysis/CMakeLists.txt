add_library(rtt_analysis src/spot.cpp src/opd.cpp)
add_library(rtt::analysis ALIAS rtt_analysis)
target_include_directories(rtt_analysis PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_analysis PUBLIC rtt::trace rtt::paraxial rtt::compile rtt::math)
rtt_target_defaults(rtt_analysis LIBRARY)

if(RTT_BUILD_TESTS)
  # rtt::io only loads the reference systems in the tests; TBB limits the threads in the
  # determinism test.
  rtt_add_test(rtt_analysis_tests
    SOURCES tests/test_spot.cpp tests/test_fans.cpp tests/test_opd.cpp
    LIBS rtt::analysis rtt::io TBB::tbb)
endif()
