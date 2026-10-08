add_library(rtt_geom src/instantiate.cpp)
add_library(rtt::geom ALIAS rtt_geom)
target_include_directories(rtt_geom PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_geom PUBLIC rtt::math)
rtt_target_defaults(rtt_geom LIBRARY)

if(RTT_BUILD_TESTS)
  rtt_add_test(rtt_geom_tests
    SOURCES tests/test_conic.cpp tests/test_intersect.cpp tests/test_asphere.cpp
            tests/test_newton.cpp tests/test_phase.cpp
    LIBS rtt::geom)
endif()
