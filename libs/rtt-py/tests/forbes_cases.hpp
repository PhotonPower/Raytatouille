#pragma once

// Q-con and Q-bfs cases of rtt_py_reference (#314, ADR 0043, acceptance 12): the same
// calls as test_bitwise_forbes.py (test program of rtt-py, not part of the library).

#include <filesystem>

namespace rtt::py::reference {

/// Cases of test_bitwise_forbes.py; writes <out>/<case>.<name>.npy. The systems are the
/// examples under <reference_dir>/../../examples/forbes.
void run_forbes_cases(const std::filesystem::path& reference_dir,
                      const std::filesystem::path& out,
                      int threads);

}  // namespace rtt::py::reference
