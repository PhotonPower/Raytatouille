#pragma once

// Jones pupil cases of rtt_py_reference (#208 B2): the same calls as CASES in
// test_bitwise_jones.py (test program of rtt-py, not part of the library).

#include <filesystem>

namespace rtt::py::reference {

/// Cases of test_bitwise_jones.py; writes <out>/<case>.<name>.npy.
void run_jones_cases(const std::filesystem::path& reference_dir,
                     const std::filesystem::path& out,
                     int threads);

}  // namespace rtt::py::reference
