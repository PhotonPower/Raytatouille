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
  add_test(NAME cli.format.check
    COMMAND rtt_cli format --check ${RTT_REFERENCE_DIR}/m0/michelson.rtt.json)
endif()
