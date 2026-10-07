# Header-only registry of the stable diagnostic codes (ADR 0022), layer Basis.
add_library(rtt_diagnostics INTERFACE)
add_library(rtt::diagnostics ALIAS rtt_diagnostics)
target_include_directories(rtt_diagnostics INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_compile_features(rtt_diagnostics INTERFACE cxx_std_20)

if(RTT_BUILD_TESTS)
  rtt_add_test(rtt_diagnostics_tests SOURCES tests/test_codes.cpp LIBS rtt::diagnostics)
endif()
