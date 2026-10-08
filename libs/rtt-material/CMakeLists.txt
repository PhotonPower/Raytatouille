add_library(rtt_material src/material.cpp src/dispersion.cpp src/tabulated.cpp src/agf.cpp)
add_library(rtt::material ALIAS rtt_material)
target_include_directories(rtt_material PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_material PUBLIC rtt::math PRIVATE rtt::diagnostics)
rtt_target_defaults(rtt_material LIBRARY)

if(RTT_BUILD_TESTS)
  find_package(Threads REQUIRED)
  rtt_add_test(rtt_material_tests
    SOURCES tests/test_material.cpp tests/test_dispersion.cpp tests/test_tabulated.cpp
            tests/test_agf.cpp tests/test_air.cpp tests/test_library.cpp
            tests/test_load_warnings.cpp
    LIBS rtt::material rtt::diagnostics Threads::Threads)
  # Test catalogues (tests/catalogs) next to the reference systems.
  target_compile_definitions(rtt_material_tests PRIVATE
    RTT_CATALOG_DIR="${PROJECT_SOURCE_DIR}/tests/catalogs")
endif()
