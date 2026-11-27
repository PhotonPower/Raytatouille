# Non-sequential ray tracing (ADR 0057): scene with bodies and generated mantles, BVH.
add_library(rtt_nonseq src/bvh.cpp src/scene.cpp src/trace.cpp)
add_library(rtt::nonseq ALIAS rtt_nonseq)
target_include_directories(rtt_nonseq PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(rtt_nonseq PUBLIC rtt::trace rtt::compile rtt::geom rtt::math rtt::model
  rtt::diagnostics PRIVATE TBB::tbb)
rtt_target_defaults(rtt_nonseq LIBRARY)

if(RTT_BUILD_TESTS)
  rtt_add_test(rtt_nonseq_tests
    SOURCES tests/test_bvh.cpp tests/test_scene.cpp tests/test_trace.cpp
    LIBS rtt::nonseq rtt::model_test_support TBB::tbb)
endif()
