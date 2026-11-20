#pragma once

// rtt trace and rtt analyze (ADR 0052, points 2 to 5), part of rtt::cli::run.

#include <ostream>
#include <string>
#include <vector>

#include "rtt/trace/run_control.hpp"

namespace rtt::cli {

/// rtt trace ..., rtt analyze <kind> ... or rtt optimize ... with the arguments including the
/// command; `cancel` stops the computation (exit 130). The exit code of ADR 0052, point 4.
[[nodiscard]] int run_compute(const std::vector<std::string>& args,
                              std::ostream& out,
                              std::ostream& err,
                              const trace::CancelToken& cancel);

}  // namespace rtt::cli
