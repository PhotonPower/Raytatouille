add_library(rtt_coating src/transfer_matrix.cpp src/thickness.cpp src/ideal.cpp src/tabulated.cpp)
add_library(rtt::coating ALIAS rtt_coating)
target_include_directories(rtt_coating PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
# Layer physics (docs/architecture.md): only rtt-math; indices are passed in, no rtt-material.
target_link_libraries(rtt_coating PUBLIC rtt::math)
rtt_target_defaults(rtt_coating LIBRARY)

if(RTT_BUILD_TESTS)
  rtt_add_test(rtt_coating_tests
    SOURCES tests/test_transfer_matrix.cpp tests/test_coatings.cpp
    LIBS rtt::coating)
endif()
