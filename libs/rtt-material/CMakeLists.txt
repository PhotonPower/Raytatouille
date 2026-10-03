add_library(rtt_material src/material.cpp)
add_library(rtt::material ALIAS rtt_material)
target_include_directories(rtt_material PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_material PUBLIC rtt::math)
rtt_target_defaults(rtt_material LIBRARY)

if(RTT_BUILD_TESTS)
  find_package(Threads REQUIRED)
  rtt_add_test(rtt_material_tests
    SOURCES tests/test_material.cpp
    LIBS rtt::material Threads::Threads)
endif()
