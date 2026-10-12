#pragma once

// M6 acceptance cases of rtt_py_reference (#212): the same calls as CASES in test_bitwise_m6.py
// on the reference systems of the M6 acceptance (test program of rtt-py, not part of the library).

#include <filesystem>

namespace rtt::py::reference {

/// Cases of test_bitwise_m6.py; writes <out>/<case>.<name>.npy.
void run_m6_cases(const std::filesystem::path& reference_dir,
                  const std::filesystem::path& out,
                  int threads);

}  // namespace rtt::py::reference
