add_library(rtt_io src/json_io.cpp src/json_format.cpp)
add_library(rtt::io ALIAS rtt_io)
target_include_directories(rtt_io PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(rtt_io PUBLIC rtt::model PRIVATE nlohmann_json::nlohmann_json)
rtt_target_defaults(rtt_io LIBRARY)

if(RTT_BUILD_TESTS)
  rtt_add_test(rtt_io_tests
    SOURCES tests/test_roundtrip.cpp tests/test_parse_errors.cpp
    LIBS rtt::io)
endif()

if(RTT_BUILD_TESTS)
  target_include_directories(rtt_io_tests PRIVATE ${PROJECT_SOURCE_DIR}/libs/rtt-model/tests)
endif()
