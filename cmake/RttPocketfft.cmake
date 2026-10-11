# pocketfft (header-only, BSD-3-Clause; ADR 0010, ADR 0033 point 1) as the interface target
# rtt::pocketfft.
#
# The dependency comes from vcpkg (vcpkg.json). Presets with system packages find the header with
# find_path; if it is not there (no package for Ubuntu 24.04 or MSYS2), FetchContent loads exactly
# the commit of the vcpkg port in the builtin-baseline, checked by the same SHA512. There is no
# copy in the repository. tools/check_pins.py compares the two pins below with the vcpkg port in CI.
#
# The include directory is SYSTEM on purpose: GCC 16 reports -Wmaybe-uninitialized inside the
# header at -O2 (pocketfft_hdronly.h:3307), and we build with -Werror; clang-tidy skips SYSTEM
# headers through its HeaderFilter. Only libs/rtt-analysis/src/fft.cpp includes the header.

# Pins of the vcpkg port "pocketfft" 2024-11-30 (ports/pocketfft/portfile.cmake at the
# builtin-baseline). Keep both lines in this form: tools/check_pins.py reads them.
set(RTT_POCKETFFT_REF "0fa0ef591e38c2758e3184c6c23e497b9f732ffa")
set(RTT_POCKETFFT_SHA512 "b606512d452ba87aa0024a464c34df52be2e808f8cca5abbcd63a37830a1a4e4f45ca971f5da2ee5456694ab7c6d840c0b0af82be0db6af3df0c6aa4c00336bf")

option(RTT_FETCH_POCKETFFT
       "Download pocketfft at the vcpkg pin if find_path does not find pocketfft_hdronly.h" ON)

find_path(RTT_POCKETFFT_INCLUDE_DIR pocketfft_hdronly.h)
if(NOT RTT_POCKETFFT_INCLUDE_DIR)
  if(RTT_FETCH_POCKETFFT)
    include(FetchContent)
    # The archive has no CMakeLists.txt: MakeAvailable only downloads and extracts it.
    FetchContent_Declare(pocketfft
      URL "https://github.com/mreineck/pocketfft/archive/${RTT_POCKETFFT_REF}.tar.gz"
      URL_HASH "SHA512=${RTT_POCKETFFT_SHA512}"
      DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    FetchContent_MakeAvailable(pocketfft)
    set(RTT_POCKETFFT_INCLUDE_DIR "${pocketfft_SOURCE_DIR}" CACHE PATH
        "Directory of pocketfft_hdronly.h" FORCE)
  else()
    message(FATAL_ERROR
      "pocketfft_hdronly.h not found. Use vcpkg (vcpkg.json), set RTT_POCKETFFT_INCLUDE_DIR, or "
      "enable RTT_FETCH_POCKETFFT to download the pinned version.")
  endif()
endif()
if(NOT EXISTS "${RTT_POCKETFFT_INCLUDE_DIR}/pocketfft_hdronly.h")
  message(FATAL_ERROR "pocketfft_hdronly.h is not in RTT_POCKETFFT_INCLUDE_DIR "
                      "('${RTT_POCKETFFT_INCLUDE_DIR}')")
endif()

add_library(rtt_pocketfft INTERFACE)
add_library(rtt::pocketfft ALIAS rtt_pocketfft)
target_include_directories(rtt_pocketfft SYSTEM INTERFACE "${RTT_POCKETFFT_INCLUDE_DIR}")
