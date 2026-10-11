#pragma once

// Q and P per slot cases of rtt_py_reference (#207): the same calls as test_bitwise_q.py (test
// program of rtt-py, not part of the library).

#include <filesystem>

namespace rtt::py::reference {

/// Cases of test_bitwise_q.py; writes <out>/<case>.<name>.npy.
void run_q_cases(const std::filesystem::path& reference_dir,
                 const std::filesystem::path& out,
                 int threads);

}  // namespace rtt::py::reference
