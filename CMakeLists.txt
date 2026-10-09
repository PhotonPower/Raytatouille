cmake_minimum_required(VERSION 3.25)

project(Raytatouille
  VERSION 0.6.0
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
option(RTT_BUILD_PYTHON "Build the Python package raytatouille (libs/rtt-py, ADR 0002)" OFF)

# The static libraries are linked into the Python extension module, a shared library.
if(RTT_BUILD_PYTHON)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(RttTargets)

find_package(Eigen3 CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED)

if(RTT_BUILD_TESTS)
  find_package(Catch2 3 CONFIG REQUIRED)
  include(CTest)
  include(Catch)
  enable_testing()
endif()

set(RTT_REFERENCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests/reference")

add_subdirectory(libs/rtt-math)
add_subdirectory(libs/rtt-json)
add_subdirectory(libs/rtt-diagnostics)
add_subdirectory(libs/rtt-geom)
add_subdirectory(libs/rtt-polar)
add_subdirectory(libs/rtt-material)
add_subdirectory(libs/rtt-coating)
add_subdirectory(libs/rtt-model)
add_subdirectory(libs/rtt-compile)
add_subdirectory(libs/rtt-paraxial)
add_subdirectory(libs/rtt-trace)
add_subdirectory(libs/rtt-analysis)
add_subdirectory(libs/rtt-io)

if(RTT_BUILD_PYTHON)
  add_subdirectory(libs/rtt-py)
endif()

if(RTT_BUILD_APPS)
  add_subdirectory(apps/rtt-cli)
  add_subdirectory(apps/rtt-bench)
endif()
