add_library(rtt_io src/json_io.cpp src/json_format.cpp src/json_patch.cpp src/edit.cpp src/zmx.cpp
  src/export_path.cpp src/zmx_export.cpp src/surface_table.cpp)
add_library(rtt::io ALIAS rtt_io)
target_include_directories(rtt_io PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(rtt_io PUBLIC rtt::model PRIVATE rtt::json rtt::diagnostics rtt::material
  rtt::compile rtt::paraxial rtt::trace
  nlohmann_json::nlohmann_json)
rtt_target_defaults(rtt_io LIBRARY)

if(RTT_BUILD_TESTS)
  rtt_add_test(rtt_io_tests
    SOURCES tests/test_roundtrip.cpp tests/test_parse_errors.cpp tests/test_segment_materials.cpp
            tests/test_json_patch.cpp tests/test_patch_fuzz.cpp tests/test_edit.cpp
            tests/test_format_0_3.cpp tests/test_for_each_param.cpp tests/test_format_0_4.cpp
            tests/test_format_optimization.cpp tests/test_format_0_5.cpp
            tests/test_format_tolerances.cpp tests/test_format_0_7.cpp tests/test_format_0_8.cpp
            tests/test_format_0_9.cpp tests/test_format_rotational.cpp
            tests/test_format_forbes.cpp
            tests/test_aperture_stop.cpp tests/test_aperture_reading.cpp
    LIBS rtt::io rtt::json rtt::diagnostics rtt::model_test_support nlohmann_json::nlohmann_json)
endif()

if(RTT_BUILD_TESTS)
  # ZMX import (ADR 0044), ZMX export and surface table (ADR 0045): the cases of tests/zmx
  # against the expected systems and files, the first order and a spot per field, so the test
  # needs the tracing and analysis libraries; the reference systems with coatings compile with the
  # test coatings.
  rtt_add_test(rtt_io_zmx_tests
    SOURCES tests/test_zmx.cpp tests/test_zmx_export.cpp tests/test_surface_table.cpp
    LIBS rtt::io rtt::material rtt::compile rtt::paraxial rtt::analysis rtt::trace
         rtt::diagnostics rtt::coating
         nlohmann_json::nlohmann_json)
  target_compile_definitions(rtt_io_zmx_tests PRIVATE
    RTT_CATALOG_DIR="${PROJECT_SOURCE_DIR}/tests/catalogs"
    RTT_ZMX_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/zmx")
endif()

if(RTT_BUILD_TESTS)
  # src/: the JSON Patch engine (json_patch.hpp) is tested directly on generic JSON.
  target_include_directories(rtt_io_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
endif()
