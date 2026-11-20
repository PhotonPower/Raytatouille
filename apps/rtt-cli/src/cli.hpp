#pragma once

// The command line rtt as a function (ADR 0052): main() calls it with the standard streams, the
// tests with string streams. Exit codes: 0 success, 1 invalid system or failed computation, 2
// wrong call (ADR 0052, point 4).

#include <ostream>
#include <string>
#include <vector>

namespace rtt::cli {

/// Runs the command line with the arguments after the program name; data go to `out` (binary),
/// diagnostics and messages to `err` (and, for validate, format, import and export, to `out` as
/// before). Returns the exit code.
[[nodiscard]] int run(const std::vector<std::string>& args, std::ostream& out, std::ostream& err);

}  // namespace rtt::cli
