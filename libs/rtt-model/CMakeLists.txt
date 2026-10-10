add_library(rtt_model src/element.cpp src/expression.cpp src/expression_eval.cpp
  src/parameters.cpp src/pose.cpp src/validate.cpp)
add_library(rtt::model ALIAS rtt_model)
target_include_directories(rtt_model PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_model PUBLIC rtt::diagnostics rtt::math)
rtt_target_defaults(rtt_model LIBRARY)
# Expressions of the parameter table are evaluated without FMA contraction, so that every platform
# gives the same bits (ADR 0029, point 2); MSVC does not contract without /fp:contract.
set_source_files_properties(src/expression_eval.cpp PROPERTIES COMPILE_OPTIONS
  "$<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-ffp-contract=off>")

if(RTT_BUILD_TESTS)
  # Test helpers for model systems (tests/test_support.hpp), shared with the tests of the
  # libraries above rtt-model.
  add_library(rtt_model_test_support INTERFACE)
  add_library(rtt::model_test_support ALIAS rtt_model_test_support)
  target_include_directories(rtt_model_test_support INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  target_link_libraries(rtt_model_test_support INTERFACE rtt::model)

  rtt_add_test(rtt_model_tests
    SOURCES tests/test_validate.cpp tests/test_pose.cpp tests/test_segment_materials.cpp
            tests/test_validate_0_4.cpp tests/test_validate_optimization.cpp
            tests/test_expression.cpp tests/test_parameters.cpp tests/test_validate_ideal_lens.cpp
    LIBS rtt::model rtt::model_test_support)
  # src/: the expression parser (expression.hpp) is tested directly.
  target_include_directories(rtt_model_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
endif()
