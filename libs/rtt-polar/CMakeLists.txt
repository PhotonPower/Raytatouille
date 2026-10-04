add_library(rtt_polar src/fresnel.cpp src/prt_analysis.cpp)
add_library(rtt::polar ALIAS rtt_polar)
target_include_directories(rtt_polar PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_polar PUBLIC rtt::math)
rtt_target_defaults(rtt_polar LIBRARY)

if(RTT_BUILD_TESTS)
  rtt_add_test(rtt_polar_tests
    SOURCES tests/test_fresnel.cpp tests/test_prt.cpp
    LIBS rtt::polar)
endif()
