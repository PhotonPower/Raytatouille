add_library(rtt_analysis src/spot.cpp src/opd.cpp src/field.cpp src/chromatic.cpp)
add_library(rtt::analysis ALIAS rtt_analysis)
target_include_directories(rtt_analysis PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_analysis PUBLIC rtt::trace rtt::paraxial rtt::compile rtt::math
  PRIVATE rtt::diagnostics)
rtt_target_defaults(rtt_analysis LIBRARY)

if(RTT_BUILD_TESTS)
  # rtt::io only loads the reference systems in the tests; TBB limits the threads in the
  # determinism test; rtt::polar (Stokes, diattenuation) and rtt::coating (catalogues) only for
  # the M3 acceptance (#63).
  rtt_add_test(rtt_analysis_tests
    SOURCES tests/test_spot.cpp tests/test_fans.cpp tests/test_opd.cpp tests/test_field.cpp
            tests/test_chromatic.cpp tests/test_m2_acceptance.cpp tests/test_m3_acceptance.cpp
            tests/test_no_stop.cpp tests/test_error_context.cpp tests/test_losses.cpp
            tests/test_warnings.cpp tests/test_loss_counter.cpp
    LIBS rtt::analysis rtt::io rtt::polar rtt::coating TBB::tbb)
  target_compile_definitions(rtt_analysis_tests PRIVATE
    RTT_CATALOG_DIR="${PROJECT_SOURCE_DIR}/tests/catalogs")
endif()
