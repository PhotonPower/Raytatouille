add_executable(rtt_cli main.cpp)
set_target_properties(rtt_cli PROPERTIES OUTPUT_NAME rtt)
target_link_libraries(rtt_cli PRIVATE rtt::io rtt::model rtt::material)
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
  # rtt import (ADR 0044): a good case, a file that cannot be imported (exit 1, with its code and
  # line) and a wrong call (exit 2). The output goes to the build tree.
  set(_zmx ${PROJECT_SOURCE_DIR}/libs/rtt-io/tests/zmx)
  set(_agf ${PROJECT_SOURCE_DIR}/tests/catalogs/schott.agf)
  add_test(NAME cli.import.singlet
    COMMAND rtt_cli import ${_zmx}/singlet.zmx ${CMAKE_CURRENT_BINARY_DIR}/singlet.rtt.json
            --catalog ${_agf})
  set_tests_properties(cli.import.singlet PROPERTIES PASS_REGULAR_EXPRESSION "imported to")
  add_test(NAME cli.import.no_catalog
    COMMAND rtt_cli import ${_zmx}/singlet.zmx ${CMAKE_CURRENT_BINARY_DIR}/none.rtt.json)
  set_tests_properties(cli.import.no_catalog PROPERTIES WILL_FAIL TRUE)
  add_test(NAME cli.import.no_catalog.code
    COMMAND rtt_cli import ${_zmx}/singlet.zmx ${CMAKE_CURRENT_BINARY_DIR}/none.rtt.json)
  set_tests_properties(cli.import.no_catalog.code PROPERTIES
    PASS_REGULAR_EXPRESSION "singlet\\.zmx: error line [0-9]+: .*\\[zmx\\.glass_unknown\\]")
  add_test(NAME cli.import.usage COMMAND rtt_cli import ${_zmx}/singlet.zmx)
  set_tests_properties(cli.import.usage PROPERTIES WILL_FAIL TRUE)
  # Exit 1 for a file that cannot be imported, 2 for a wrong call (ADR 0044, point 1).
  add_test(NAME cli.import.exit_1
    COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:rtt_cli> -DEXPECT=1
            "-DARGS=import;${_zmx}/singlet.zmx;${CMAKE_CURRENT_BINARY_DIR}/none.rtt.json"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/expect_exit.cmake)
  add_test(NAME cli.import.exit_2
    COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:rtt_cli> -DEXPECT=2
            "-DARGS=import;${_zmx}/singlet.zmx"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/expect_exit.cmake)
  # The written file is canonical: rtt format --check accepts it on every platform.
  set_tests_properties(cli.import.singlet PROPERTIES FIXTURES_SETUP zmx_import)
  add_test(NAME cli.import.canonical
    COMMAND rtt_cli format --check ${CMAKE_CURRENT_BINARY_DIR}/singlet.rtt.json)
  set_tests_properties(cli.import.canonical PROPERTIES FIXTURES_REQUIRED zmx_import)
endif()
