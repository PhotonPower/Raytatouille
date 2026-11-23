#pragma once

// M8 cases of rtt_py_reference (#324): the same calls as test_bitwise_m8.py (test program of
// rtt-py, not part of the library).

#include <filesystem>

namespace rtt::py::reference {

/// Cases of test_bitwise_m8.py; writes <out>/<case>.<name>.npy. The systems are the freeform
/// mirrors (S3) and the aperture systems (F1) under <reference_dir>/m8.
void run_m8_cases(const std::filesystem::path& reference_dir,
                  const std::filesystem::path& out,
                  int threads);

}  // namespace rtt::py::reference
