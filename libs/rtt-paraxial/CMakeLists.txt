add_library(rtt_paraxial src/paraxial.cpp src/seidel.cpp src/paraxial_rays.cpp src/prescription.cpp)
add_library(rtt::paraxial ALIAS rtt_paraxial)
target_include_directories(rtt_paraxial PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
# Only rtt-compile and lower layers (ADR 0016: rtt-compile < rtt-paraxial < rtt-trace).
target_link_libraries(rtt_paraxial PUBLIC rtt::compile rtt::math)
rtt_target_defaults(rtt_paraxial LIBRARY)

if(RTT_BUILD_TESTS)
  # rtt::io only loads the reference systems in the tests; the library does not depend on it.
  rtt_add_test(rtt_paraxial_tests
    SOURCES tests/test_paraxial.cpp tests/test_seidel.cpp tests/test_errors.cpp
            tests/test_prescription.cpp tests/test_ideal_lens_interim.cpp
    LIBS rtt::paraxial rtt::io)
  # SCHOTT catalogue for the dispersive field-conversion tests (#35).
  target_compile_definitions(rtt_paraxial_tests PRIVATE
    RTT_CATALOG_DIR="${PROJECT_SOURCE_DIR}/tests/catalogs")
endif()
