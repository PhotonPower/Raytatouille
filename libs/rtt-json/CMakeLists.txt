# Header-only strict JSON parsing for the file readers (ADR 0020), layer Basis.
add_library(rtt_json INTERFACE)
add_library(rtt::json ALIAS rtt_json)
target_include_directories(rtt_json INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_json INTERFACE nlohmann_json::nlohmann_json)
target_compile_features(rtt_json INTERFACE cxx_std_20)

if(RTT_BUILD_TESTS)
  rtt_add_test(rtt_json_tests SOURCES tests/test_strict.cpp LIBS rtt::json)
endif()
