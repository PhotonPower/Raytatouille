cmake_minimum_required(VERSION 3.25)

# macOS 26 is the minimum (ADR 0047): Apple's libc++ has std::from_chars for double only from
# there. The default must be set before project(); MACOSX_DEPLOYMENT_TARGET in the environment
# (e.g. from scikit-build-core) or -DCMAKE_OSX_DEPLOYMENT_TARGET take precedence. FORCE, because
# a build directory configured before this default holds the variable as an empty cache entry,
# which a plain set(... CACHE) would keep (#373, follow-up).
if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin" AND NOT CMAKE_OSX_DEPLOYMENT_TARGET
   AND NOT DEFINED ENV{MACOSX_DEPLOYMENT_TARGET})
  set(CMAKE_OSX_DEPLOYMENT_TARGET "26.0" CACHE STRING "Minimum macOS version (ADR 0047)" FORCE)
endif()

project(Raytatouille
  VERSION 0.9.0
  DESCRIPTION "Raytracing engine for optical design"
  LANGUAGES CXX)

if(APPLE AND CMAKE_OSX_DEPLOYMENT_TARGET VERSION_LESS 26.0)
  message(FATAL_ERROR "Raytatouille needs macOS 26 or later (ADR 0047), but "
    "CMAKE_OSX_DEPLOYMENT_TARGET is '${CMAKE_OSX_DEPLOYMENT_TARGET}'. Set it to 26.0, or "
    "MACOSX_DEPLOYMENT_TARGET=26.0 in the environment.")
endif()

# Apple's linker (ld-prime) warns when CMake repeats a static library on the link line, which it
# does on purpose for the link order between the rtt libraries; the warning is noise.
if(APPLE)
  add_link_options(LINKER:-no_warn_duplicate_libraries)
endif()

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
include(RttPocketfft)  # rtt::pocketfft (ADR 0033)

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
add_subdirectory(libs/rtt-optim)
add_subdirectory(libs/rtt-tolerance)
add_subdirectory(libs/rtt-io)

if(RTT_BUILD_PYTHON)
  add_subdirectory(libs/rtt-py)
endif()

if(RTT_BUILD_APPS)
  add_subdirectory(apps/rtt-cli)
  add_subdirectory(apps/rtt-bench)
endif()
