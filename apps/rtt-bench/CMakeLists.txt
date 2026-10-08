# Run-time measurement of the interactive performance (#83); see main.cpp.
add_executable(rtt_bench main.cpp)
target_link_libraries(rtt_bench PRIVATE rtt::analysis rtt::trace rtt::io rtt::material TBB::tbb)
target_compile_definitions(rtt_bench PRIVATE
  RTT_REFERENCE_DIR="${RTT_REFERENCE_DIR}"
  RTT_CATALOG_DIR="${PROJECT_SOURCE_DIR}/tests/catalogs")
rtt_target_defaults(rtt_bench)

if(RTT_BUILD_TESTS)
  # Smoke test only (one run, no timing check): the measurement itself runs in CI.
  add_test(NAME bench.smoke COMMAND rtt_bench --runs 1 --quick)
endif()
