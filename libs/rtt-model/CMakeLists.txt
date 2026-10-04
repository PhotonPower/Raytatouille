add_library(rtt_model src/element.cpp src/pose.cpp src/validate.cpp)
add_library(rtt::model ALIAS rtt_model)
target_include_directories(rtt_model PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_model PUBLIC rtt::math)
rtt_target_defaults(rtt_model LIBRARY)

if(RTT_BUILD_TESTS)
  rtt_add_test(rtt_model_tests
    SOURCES tests/test_validate.cpp tests/test_pose.cpp tests/test_segment_materials.cpp
    LIBS rtt::model)
endif()
