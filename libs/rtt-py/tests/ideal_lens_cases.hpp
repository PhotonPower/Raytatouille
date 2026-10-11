#pragma once

// Ideal lens cases of rtt_py_reference (#178): the same calls as test_bitwise_ideal_lens.py
// (test program of rtt-py, not part of the library).

#include <filesystem>

namespace rtt::py::reference {

/// Cases of test_bitwise_ideal_lens.py; writes <out>/<case>.<name>.npy.
void run_ideal_lens_cases(const std::filesystem::path& reference_dir,
                          const std::filesystem::path& out,
                          int threads);

}  // namespace rtt::py::reference
