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
  # rtt export (ADR 0045, point 5): a good case whose file imports again, the options, a system
  # that cannot be exported (exit 1, with its code) and wrong calls (exit 2).
  set(_agf_m2 ${PROJECT_SOURCE_DIR}/tests/catalogs/m2/schott.agf)
  add_test(NAME cli.export.singlet
    COMMAND rtt_cli export ${_zmx}/export_singlet.rtt.json
            ${CMAKE_CURRENT_BINARY_DIR}/export_singlet.zmx --catalog ${_agf})
  set_tests_properties(cli.export.singlet PROPERTIES PASS_REGULAR_EXPRESSION "exported to"
    FIXTURES_SETUP zmx_export)
  add_test(NAME cli.export.imports
    COMMAND rtt_cli import ${CMAKE_CURRENT_BINARY_DIR}/export_singlet.zmx
            ${CMAKE_CURRENT_BINARY_DIR}/export_singlet.rtt.json --catalog ${_agf})
  set_tests_properties(cli.export.imports PROPERTIES PASS_REGULAR_EXPRESSION "imported to"
    FIXTURES_REQUIRED zmx_export)
  add_test(NAME cli.export.options
    COMMAND rtt_cli export ${_zmx}/export_singlet.rtt.json ${CMAKE_CURRENT_BINARY_DIR}/main.ZMX
            --catalog ${_agf} --path main --config 0)
  set_tests_properties(cli.export.options PROPERTIES PASS_REGULAR_EXPRESSION "exported to")
  add_test(NAME cli.export.code
    COMMAND rtt_cli export ${RTT_REFERENCE_DIR}/m2/cooke_triplet.rtt.json
            ${CMAKE_CURRENT_BINARY_DIR}/none.zmx --catalog ${_agf_m2})
  set_tests_properties(cli.export.code PROPERTIES
    PASS_REGULAR_EXPRESSION "cooke_triplet\\.rtt\\.json: error .*\\[export\\.element_in_glass\\]")
  add_test(NAME cli.export.exit_1
    COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:rtt_cli> -DEXPECT=1
            "-DARGS=export;${RTT_REFERENCE_DIR}/m2/cooke_triplet.rtt.json;${CMAKE_CURRENT_BINARY_DIR}/none.zmx;--catalog;${_agf_m2}"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/expect_exit.cmake)
  add_test(NAME cli.export.exit_2
    COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:rtt_cli> -DEXPECT=2
            "-DARGS=export;${_zmx}/export_singlet.rtt.json;${CMAKE_CURRENT_BINARY_DIR}/out.txt;--catalog;${_agf}"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/expect_exit.cmake)
  # The surface table (ADR 0045, point 7): written, read into a base system, the result is
  # canonical; a table with an error exits 1 with its code; a table without a base system and a
  # ZMX file with one are wrong calls (exit 2).
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/semicolons.csv "row;surface;element\n")
  add_test(NAME cli.table.export
    COMMAND rtt_cli export ${_zmx}/export_singlet.rtt.json
            ${CMAKE_CURRENT_BINARY_DIR}/export_singlet.csv --catalog ${_agf})
  set_tests_properties(cli.table.export PROPERTIES PASS_REGULAR_EXPRESSION "exported to"
    FIXTURES_SETUP table_export)
  add_test(NAME cli.table.import
    COMMAND rtt_cli import ${CMAKE_CURRENT_BINARY_DIR}/export_singlet.csv
            ${CMAKE_CURRENT_BINARY_DIR}/table_singlet.rtt.json
            --base ${_zmx}/export_singlet.rtt.json --catalog ${_agf})
  set_tests_properties(cli.table.import PROPERTIES PASS_REGULAR_EXPRESSION "imported to"
    FIXTURES_REQUIRED table_export FIXTURES_SETUP table_import)
  add_test(NAME cli.table.canonical
    COMMAND rtt_cli format --check ${CMAKE_CURRENT_BINARY_DIR}/table_singlet.rtt.json)
  set_tests_properties(cli.table.canonical PROPERTIES FIXTURES_REQUIRED table_import)
  add_test(NAME cli.table.code
    COMMAND rtt_cli import ${CMAKE_CURRENT_BINARY_DIR}/semicolons.csv
            ${CMAKE_CURRENT_BINARY_DIR}/none.rtt.json --base ${_zmx}/export_singlet.rtt.json --catalog ${_agf})
  set_tests_properties(cli.table.code PROPERTIES
    PASS_REGULAR_EXPRESSION "semicolons\\.csv: error .*\\[table\\.")
  add_test(NAME cli.table.exit_1
    COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:rtt_cli> -DEXPECT=1
            "-DARGS=import;${CMAKE_CURRENT_BINARY_DIR}/semicolons.csv;${CMAKE_CURRENT_BINARY_DIR}/none.rtt.json;--base;${_zmx}/export_singlet.rtt.json;--catalog;${_agf}"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/expect_exit.cmake)
  add_test(NAME cli.table.exit_2_no_base
    COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:rtt_cli> -DEXPECT=2
            "-DARGS=import;${_zmx}/export_singlet.csv;${CMAKE_CURRENT_BINARY_DIR}/none.rtt.json;--catalog;${_agf}"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/expect_exit.cmake)
  add_test(NAME cli.table.exit_2_zmx_base
    COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:rtt_cli> -DEXPECT=2
            "-DARGS=import;${_zmx}/singlet.zmx;${CMAKE_CURRENT_BINARY_DIR}/none.rtt.json;--base;${_zmx}/export_singlet.rtt.json"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/expect_exit.cmake)
  add_test(NAME cli.export.exit_2_config
    COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:rtt_cli> -DEXPECT=2
            "-DARGS=export;${_zmx}/export_singlet.rtt.json;${CMAKE_CURRENT_BINARY_DIR}/out.zmx;--config;x"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/expect_exit.cmake)
endif()
