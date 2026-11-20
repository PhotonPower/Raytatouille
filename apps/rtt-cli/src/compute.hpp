#pragma once

// rtt trace and rtt analyze (ADR 0052, points 2 to 5), part of rtt::cli::run.

#include <ostream>
#include <string>
#include <vector>

namespace rtt::cli {

/// rtt trace ... or rtt analyze <kind> ... with the arguments including the command; the exit
/// code of ADR 0052, point 4.
[[nodiscard]] int run_compute(const std::vector<std::string>& args,
                              std::ostream& out,
                              std::ostream& err);

}  // namespace rtt::cli
