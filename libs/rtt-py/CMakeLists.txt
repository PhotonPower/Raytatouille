# Python package raytatouille (ADR 0002, ADR 0018): extension module raytatouille._core built
# with nanobind, the pure Python modules from python/raytatouille and generated type stubs.
#
# nanobind is a Python build dependency (pyproject.toml, ADR 0018); its CMake package is found
# through the interpreter that builds the package.

find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_ROOT
  RESULT_VARIABLE rtt_nanobind_result)
if(NOT rtt_nanobind_result EQUAL 0)
  message(FATAL_ERROR "nanobind not found for ${Python_EXECUTABLE}; "
                      "pip install -r libs/rtt-py/tests/requirements.txt")
endif()
find_package(nanobind CONFIG REQUIRED)

# Package directory in the build tree: _core, the Python modules and the stubs side by side, so
# that pytest can import raytatouille from here (PYTHONPATH = RTT_PY_BUILD_DIR).
set(RTT_PY_BUILD_DIR "${CMAKE_BINARY_DIR}/python")
set(RTT_PY_PACKAGE_DIR "${RTT_PY_BUILD_DIR}/raytatouille")

nanobind_add_module(_core NB_STATIC
  src/module.cpp
  src/errors.cpp
  src/bind_model.cpp
  src/bind_compile.cpp
  src/bind_paraxial.cpp
  src/bind_trace.cpp)
target_link_libraries(_core PRIVATE
  rtt::io rtt::model rtt::material rtt::compile rtt::paraxial rtt::trace TBB::tbb)
rtt_target_defaults(_core)
# The $<1:...> keeps multi-config generators from appending a configuration directory.
set_target_properties(_core PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY "$<1:${RTT_PY_PACKAGE_DIR}>"
  RUNTIME_OUTPUT_DIRECTORY "$<1:${RTT_PY_PACKAGE_DIR}>")

set(RTT_PY_SOURCES
  python/raytatouille/__init__.py
  python/raytatouille/errors.py
  python/raytatouille/paraxial.py
  python/raytatouille/trace.py
  python/raytatouille/py.typed)
set(RTT_PY_COPIES "")
foreach(file IN LISTS RTT_PY_SOURCES)
  get_filename_component(name "${file}" NAME)
  add_custom_command(
    OUTPUT "${RTT_PY_PACKAGE_DIR}/${name}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/${file}" "${RTT_PY_PACKAGE_DIR}/${name}"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${file}")
  list(APPEND RTT_PY_COPIES "${RTT_PY_PACKAGE_DIR}/${name}")
endforeach()
add_custom_target(rtt_py_package ALL DEPENDS ${RTT_PY_COPIES})

# TBB is a shared library; on Windows its DLL has to sit next to _core. Only TBB: the runtime
# DLLs of _core also include python3*.dll, which must not be copied into the package.
if(WIN32)
  add_custom_command(TARGET _core POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$<TARGET_FILE:TBB::tbb>" "${RTT_PY_PACKAGE_DIR}")
endif()

nanobind_add_stub(rtt_py_stub
  MODULE raytatouille._core
  OUTPUT "${RTT_PY_PACKAGE_DIR}/_core.pyi"
  PYTHON_PATH "${RTT_PY_BUILD_DIR}"
  DEPENDS _core rtt_py_package)

# Installation into the wheel (scikit-build-core): the pure Python modules come from
# wheel.packages in pyproject.toml.
install(TARGETS _core LIBRARY DESTINATION raytatouille RUNTIME DESTINATION raytatouille)
install(FILES "${RTT_PY_PACKAGE_DIR}/_core.pyi" DESTINATION raytatouille)
if(WIN32)
  install(FILES "$<TARGET_FILE:TBB::tbb>" DESTINATION raytatouille)
endif()

if(RTT_BUILD_TESTS)
  # Writes the C++ results of the cases in tests/test_bitwise.py as .npy files.
  add_executable(rtt_py_reference tests/rtt_py_reference.cpp)
  target_link_libraries(rtt_py_reference PRIVATE
    rtt::io rtt::material rtt::compile rtt::paraxial rtt::trace TBB::tbb)
  rtt_target_defaults(rtt_py_reference)

  add_test(NAME rtt_py_pytest
    COMMAND "${Python_EXECUTABLE}" -m pytest -q "${CMAKE_CURRENT_SOURCE_DIR}/tests"
            -p no:cacheprovider)
  set_tests_properties(rtt_py_pytest PROPERTIES
    ENVIRONMENT
      "PYTHONPATH=${RTT_PY_BUILD_DIR};RTT_PY_REFERENCE_EXE=$<TARGET_FILE:rtt_py_reference>;RTT_REFERENCE_DIR=${RTT_REFERENCE_DIR};RTT_CATALOG_DIR=${PROJECT_SOURCE_DIR}/tests/catalogs")
  add_test(NAME rtt_py_mypy
    COMMAND "${Python_EXECUTABLE}" -m mypy --config-file "${PROJECT_SOURCE_DIR}/pyproject.toml"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests")
  set_tests_properties(rtt_py_mypy PROPERTIES
    ENVIRONMENT "MYPYPATH=${RTT_PY_BUILD_DIR}")
endif()
