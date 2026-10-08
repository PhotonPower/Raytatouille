add_library(rtt_io src/json_io.cpp src/json_format.cpp src/json_patch.cpp src/edit.cpp)
add_library(rtt::io ALIAS rtt_io)
target_include_directories(rtt_io PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(rtt_io PUBLIC rtt::model PRIVATE rtt::json rtt::diagnostics
  nlohmann_json::nlohmann_json)
rtt_target_defaults(rtt_io LIBRARY)

if(RTT_BUILD_TESTS)
  rtt_add_test(rtt_io_tests
    SOURCES tests/test_roundtrip.cpp tests/test_parse_errors.cpp tests/test_segment_materials.cpp
            tests/test_json_patch.cpp tests/test_patch_fuzz.cpp tests/test_edit.cpp
    LIBS rtt::io rtt::json rtt::diagnostics nlohmann_json::nlohmann_json)
endif()

if(RTT_BUILD_TESTS)
  # src/: the JSON Patch engine (json_patch.hpp) is tested directly on generic JSON.
  target_include_directories(rtt_io_tests PRIVATE ${PROJECT_SOURCE_DIR}/libs/rtt-model/tests
    ${CMAKE_CURRENT_SOURCE_DIR}/src)
endif()
