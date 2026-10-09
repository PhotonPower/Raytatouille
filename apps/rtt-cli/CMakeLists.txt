add_executable(rtt_cli main.cpp)
set_target_properties(rtt_cli PROPERTIES OUTPUT_NAME rtt)
target_link_libraries(rtt_cli PRIVATE rtt::io rtt::model)
target_compile_definitions(rtt_cli PRIVATE RTT_VERSION="${PROJECT_VERSION}")
rtt_target_defaults(rtt_cli)

if(RTT_BUILD_TESTS)
  add_test(NAME cli.version COMMAND rtt_cli --version)
  add_test(NAME cli.validate.singlet
    COMMAND rtt_cli validate ${RTT_REFERENCE_DIR}/m0/singlet.rtt.json)
  add_test(NAME cli.validate.missing COMMAND rtt_cli validate does-not-exist.rtt.json)
  set_tests_properties(cli.validate.missing PROPERTIES WILL_FAIL TRUE)
  # Diagnostics are printed with their stable code (ADR 0022); the exit code stays 1.
  add_test(NAME cli.validate.code
    COMMAND rtt_cli validate ${CMAKE_CURRENT_SOURCE_DIR}/tests/two_references.rtt.json)
  set_tests_properties(cli.validate.code PROPERTIES
    PASS_REGULAR_EXPRESSION
    "two_references\\.rtt\\.json: error \\[wavelengths\\.reference_count\\] /wavelengths: ")
  add_test(NAME cli.validate.code.exit
    COMMAND rtt_cli validate ${CMAKE_CURRENT_SOURCE_DIR}/tests/two_references.rtt.json)
  set_tests_properties(cli.validate.code.exit PROPERTIES WILL_FAIL TRUE)
  add_test(NAME cli.format.check
    COMMAND rtt_cli format --check ${RTT_REFERENCE_DIR}/m0/michelson.rtt.json)
  # A pickup dropped by the migration 0.3 -> 0.4 is reported, not lost silently (#162); the
  # file stays valid (warnings only).
  add_test(NAME cli.validate.pickup_dropped
    COMMAND rtt_cli validate ${CMAKE_CURRENT_SOURCE_DIR}/tests/pickup_0_3.rtt.json)
  set_tests_properties(cli.validate.pickup_dropped PROPERTIES
    PASS_REGULAR_EXPRESSION
    "pickup_0_3\\.rtt\\.json: warning \\[io\\.pickup_dropped\\] /root/children/1/pose/position/2/pickup: pickup '2 \\* 3'")
endif()
