#pragma once

// Pupil cases of rtt_py_reference (#209, P1a): the same calls as CASES in test_bitwise_pupil.py
// (test program of rtt-py, not part of the library).

#include <filesystem>

namespace rtt::py::reference {

/// Cases of test_bitwise_pupil.py; writes <out>/<case>.<name>.npy.
void run_pupil_cases(const std::filesystem::path& reference_dir,
                     const std::filesystem::path& catalog_dir,
                     const std::filesystem::path& out,
                     int threads);

}  // namespace rtt::py::reference
