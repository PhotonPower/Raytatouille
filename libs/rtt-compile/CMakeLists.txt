add_library(rtt_compile src/compile.cpp src/layout.cpp src/ghosts.cpp)
add_library(rtt::compile ALIAS rtt_compile)
target_include_directories(rtt_compile PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_compile PUBLIC rtt::diagnostics rtt::model rtt::geom rtt::material
  rtt::coating rtt::math)
rtt_target_defaults(rtt_compile LIBRARY)

if(RTT_BUILD_TESTS)
  # rtt::io only loads the reference systems in the tests; the library does not depend on it.
  rtt_add_test(rtt_compile_tests
    SOURCES tests/test_compile.cpp tests/test_media.cpp tests/test_coatings.cpp tests/test_layout.cpp
            tests/test_diagnostic_codes.cpp tests/test_errors.cpp tests/test_ghosts.cpp
    LIBS rtt::compile rtt::io)
  target_include_directories(rtt_compile_tests PRIVATE ${PROJECT_SOURCE_DIR}/libs/rtt-model/tests)
  target_compile_definitions(rtt_compile_tests PRIVATE
    RTT_CATALOG_DIR="${PROJECT_SOURCE_DIR}/tests/catalogs")
endif()
