#pragma once

// Configuration cases of rtt_py_reference (#169): the same calls as CASES in
// test_bitwise_configurations.py (test program of rtt-py, not part of the library).

#include <filesystem>

namespace rtt::py::reference {

/// Cases of test_bitwise_configurations.py; writes <out>/<case>.<name>.npy.
void run_configuration_cases(const std::filesystem::path& reference_dir,
                             const std::filesystem::path& catalog_dir,
                             const std::filesystem::path& out,
                             int threads);

}  // namespace rtt::py::reference
