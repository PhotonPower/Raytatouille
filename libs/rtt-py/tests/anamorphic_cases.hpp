#pragma once

// Anamorphic cases of rtt_py_reference (#315, ADR 0042, acceptance 9): the same calls as
// test_bitwise_anamorphic.py (test program of rtt-py, not part of the library).

#include <filesystem>

namespace rtt::py::reference {

/// Cases of test_bitwise_anamorphic.py; writes <out>/<case>.<name>.npy. The systems are the
/// examples under <reference_dir>/../../examples/anamorphic.
void run_anamorphic_cases(const std::filesystem::path& reference_dir,
                          const std::filesystem::path& out,
                          int threads);

}  // namespace rtt::py::reference
