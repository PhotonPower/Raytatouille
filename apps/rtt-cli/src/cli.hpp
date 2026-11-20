#pragma once

// The command line rtt as a function (ADR 0052): main() calls it with the standard streams and
// the interrupt token, the tests with string streams. Exit codes: 0 success, 1 invalid system or
// failed computation, 2 wrong call, 130 cancelled (ADR 0052, point 4).

#include <ostream>
#include <string>
#include <vector>

#include "rtt/trace/run_control.hpp"

namespace rtt::cli {

/// Exit code after a cancellation (Ctrl-C), on every platform (ADR 0052, point 4).
inline constexpr int kCancelled = 130;

/// Runs the command line with the arguments after the program name; data go to `out` (binary),
/// diagnostics and messages to `err` (and, for validate, format, import and export, to `out` as
/// before). `cancel` stops trace, analyze and optimize (ADR 0052, point 6). Returns the exit code.
[[nodiscard]] int run(const std::vector<std::string>& args,
                      std::ostream& out,
                      std::ostream& err,
                      const trace::CancelToken& cancel);

/// run() with a token nobody cancels.
[[nodiscard]] int run(const std::vector<std::string>& args, std::ostream& out, std::ostream& err);

/// The token of the process that the interrupt handler cancels (ADR 0052, point 6). The one
/// piece of mutable global state of the command line, an exception from rule 4 for apps/rtt-cli
/// that ADR 0052 records.
[[nodiscard]] const trace::CancelToken& interrupt_token();

/// What the handler does on Ctrl-C (SIGINT; CTRL_C_EVENT and CTRL_BREAK_EVENT under Windows):
/// cancels interrupt_token() with one lock-free atomic store. Returns true the first time and
/// false after that, when the default action (end the process at once) follows.
bool on_interrupt() noexcept;

/// Installs the handler for Ctrl-C (main() of rtt).
void install_interrupt_handler();

}  // namespace rtt::cli
