# Shared target settings for all Raytatouille libraries, apps and tests.

if(RTT_ENABLE_SANITIZERS AND NOT MSVC)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all)
  add_link_options(-fsanitize=address,undefined)
endif()

if(RTT_ENABLE_CLANG_TIDY)
  find_program(RTT_CLANG_TIDY_EXE NAMES clang-tidy clang-tidy-18 REQUIRED)
endif()

# rtt_target_defaults(<target> [LIBRARY])
#   Applies warnings and (for libraries) clang-tidy.
function(rtt_target_defaults target)
  cmake_parse_arguments(ARG "LIBRARY" "" "" ${ARGN})
  target_compile_features(${target} PUBLIC cxx_std_20)

  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /permissive- /utf-8 /EHsc /Zc:__cplusplus)
    if(RTT_WARNINGS_AS_ERRORS)
      target_compile_options(${target} PRIVATE /WX)
    endif()
  else()
    target_compile_options(${target} PRIVATE
      -Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor -Wold-style-cast
      -Wcast-align -Woverloaded-virtual -Wnull-dereference -Wdouble-promotion
      -Wformat=2 -Wimplicit-fallthrough)
    if(RTT_WARNINGS_AS_ERRORS)
      target_compile_options(${target} PRIVATE -Werror)
    endif()
  endif()

  if(ARG_LIBRARY AND RTT_ENABLE_CLANG_TIDY)
    set_target_properties(${target} PROPERTIES CXX_CLANG_TIDY "${RTT_CLANG_TIDY_EXE}")
  endif()
endfunction()

# rtt_add_test(<name> SOURCES ... LIBS ...)
function(rtt_add_test name)
  cmake_parse_arguments(ARG "" "" "SOURCES;LIBS" ${ARGN})
  add_executable(${name} ${ARG_SOURCES})
  target_link_libraries(${name} PRIVATE ${ARG_LIBS} Catch2::Catch2WithMain)
  target_compile_definitions(${name} PRIVATE RTT_REFERENCE_DIR="${RTT_REFERENCE_DIR}")
  rtt_target_defaults(${name})
  catch_discover_tests(${name})
endfunction()
