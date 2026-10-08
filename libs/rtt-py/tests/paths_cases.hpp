#pragma once

// Path evaluation and ghost cases of rtt_py_reference (#133): the same calls as CASES in
// test_bitwise_paths.py (test program of rtt-py, not part of the library).

#include <filesystem>

namespace rtt::py::reference {

/// Cases of test_bitwise_paths.py; writes <out>/<case>.<name>.npy.
void run_paths_cases(const std::filesystem::path& reference_dir,
                     const std::filesystem::path& catalog_dir,
                     const std::filesystem::path& out,
                     int threads);

}  // namespace rtt::py::reference
