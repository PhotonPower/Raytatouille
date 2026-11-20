#pragma once

// Cases of the C++ result writer (ADR 0052, point 1) for rtt_py_reference: one or more results of
// each of the 11 types, written with rtt::io::to_result_json as <out>/json_<case>.result.json;
// test_result_json.py makes the same calls in Python and compares to_json byte for byte (test
// program of rtt-py, not part of the library).

#include <filesystem>

namespace rtt::py::reference {

/// Writes <out>/json_<case>.result.json for the cases of test_result_json.py.
void run_result_json_cases(const std::filesystem::path& reference_dir,
                           const std::filesystem::path& catalog_dir,
                           const std::filesystem::path& out,
                           int threads);

}  // namespace rtt::py::reference
