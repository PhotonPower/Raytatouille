cmake_minimum_required(VERSION 3.25)

project(Raytatouille
  VERSION 0.1.0
  DESCRIPTION "Raytracing engine for optical design"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(RTT_BUILD_TESTS "Build unit tests" ON)
option(RTT_BUILD_APPS "Build command line tools" ON)
option(RTT_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(RTT_ENABLE_SANITIZERS "Build with AddressSanitizer and UBSan (GCC/Clang)" OFF)
option(RTT_ENABLE_CLANG_TIDY "Run clang-tidy on library sources during the build" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(RttTargets)

find_package(Eigen3 CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)

if(RTT_BUILD_TESTS)
  find_package(Catch2 3 CONFIG REQUIRED)
  include(CTest)
  include(Catch)
  enable_testing()
endif()

set(RTT_REFERENCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests/reference")

add_subdirectory(libs/rtt-math)
add_subdirectory(libs/rtt-geom)
add_subdirectory(libs/rtt-material)
add_subdirectory(libs/rtt-model)
add_subdirectory(libs/rtt-compile)
add_subdirectory(libs/rtt-trace)
add_subdirectory(libs/rtt-io)

if(RTT_BUILD_APPS)
  add_subdirectory(apps/rtt-cli)
endif()
