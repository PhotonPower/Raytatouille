add_library(rtt_compile src/compile.cpp)
add_library(rtt::compile ALIAS rtt_compile)
target_include_directories(rtt_compile PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_compile PUBLIC rtt::model rtt::geom rtt::material rtt::math)
rtt_target_defaults(rtt_compile LIBRARY)

if(RTT_BUILD_TESTS)
  # rtt::io only loads the reference systems in the tests; the library does not depend on it.
  rtt_add_test(rtt_compile_tests
    SOURCES tests/test_compile.cpp tests/test_media.cpp
    LIBS rtt::compile rtt::io)
  target_include_directories(rtt_compile_tests PRIVATE ${PROJECT_SOURCE_DIR}/libs/rtt-model/tests)
endif()
