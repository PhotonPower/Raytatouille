add_library(rtt_model src/element.cpp src/pose.cpp src/validate.cpp)
add_library(rtt::model ALIAS rtt_model)
target_include_directories(rtt_model PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_model PUBLIC rtt::diagnostics rtt::math)
rtt_target_defaults(rtt_model LIBRARY)

if(RTT_BUILD_TESTS)
  # Test helpers for model systems (tests/test_support.hpp), shared with the tests of the
  # libraries above rtt-model.
  add_library(rtt_model_test_support INTERFACE)
  add_library(rtt::model_test_support ALIAS rtt_model_test_support)
  target_include_directories(rtt_model_test_support INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  target_link_libraries(rtt_model_test_support INTERFACE rtt::model)

  rtt_add_test(rtt_model_tests
    SOURCES tests/test_validate.cpp tests/test_pose.cpp tests/test_segment_materials.cpp
    LIBS rtt::model rtt::model_test_support)
endif()
