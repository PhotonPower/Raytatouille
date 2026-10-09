#pragma once

// Crystal cases of rtt_py_reference (#134): the same calls as CASES in test_bitwise_crystal.py
// (test program of rtt-py, not part of the library).

#include <filesystem>

namespace rtt::py::reference {

/// Cases of test_bitwise_crystal.py; writes <out>/<case>.<name>.npy.
void run_crystal_cases(const std::filesystem::path& reference_dir,
                       const std::filesystem::path& catalog_dir,
                       const std::filesystem::path& out,
                       int threads);

}  // namespace rtt::py::reference
