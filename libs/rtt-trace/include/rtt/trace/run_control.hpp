#pragma once

/// @file run_control.hpp
/// Cancellation and progress reporting for long computations (#83, G4; docs/architecture.md,
/// "Abbruch und Fortschritt").
///
/// A RunControl is passed to the overloads of SequentialTracer::trace(), make_rays() and the
/// analyses that take one. Without it (or with an empty one) the computation runs exactly as
/// before. With it the work is split into blocks of `block_size` rays:
/// - before a block starts, the CancelToken is checked; after a cancellation request no new
///   block starts, and the call throws Cancelled once the running blocks are done (at the API
///   boundary, never from inside the parallel part);
/// - after each block the progress callback may be called (serialised, from any thread, at most
///   every `min_interval`; the last call of a completed stage reports done == total);
/// - an exception thrown by the callback stops the run like a cancellation and is rethrown
///   after the parallel part.
/// Results do not depend on the block size or the number of threads: every ray is computed
/// independently and exactly once (ADR 0004, addendum #83).

#include <atomic>
#include <chrono>
#include <cstddef>
#include <exception>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string_view>

namespace rtt::trace {

/// Shared cancellation flag. Copies refer to the same flag, so a GUI keeps one copy and passes
/// another to the computation; request_cancel() may be called from any thread at any time.
class CancelToken {
 public:
  CancelToken() : flag_(std::make_shared<std::atomic<bool>>(false)) {}
  // Copy only: a moved-from token would hold no flag. Declaring the copy operations suppresses
  // the implicit moves, so a "move" copies the shared pointer and both stay valid.
  CancelToken(const CancelToken&) = default;
  CancelToken& operator=(const CancelToken&) = default;
  ~CancelToken() = default;

  /// Asks the computation that holds this token to stop.
  void request_cancel() noexcept { flag_->store(true, std::memory_order_relaxed); }

  /// True once request_cancel() was called on any copy.
  [[nodiscard]] bool cancelled() const noexcept { return flag_->load(std::memory_order_relaxed); }

 private:
  std::shared_ptr<std::atomic<bool>> flag_;
};

/// Progress of one stage of a computation.
struct Progress {
  std::size_t done = 0;    ///< work items (rays) finished in this stage
  std::size_t total = 0;   ///< work items of this stage
  std::string_view stage;  ///< name of the stage, e.g. "aim" or "trace"
};

/// Default number of rays per block when a RunControl is active: short enough for a quick
/// reaction to a cancellation (one block per worker: well below a millisecond for a trace of
/// typical lenses, a few milliseconds for real aiming), long enough that the per-block overhead
/// (an atomic add, a flag check) is negligible (ADR 0004, addendum #83).
inline constexpr std::size_t kRunBlockSize = 256;

/// Cancellation and progress for one call. Empty (no token and no callback) means: run as
/// without a RunControl.
struct RunControl {
  std::optional<CancelToken> cancel;  ///< checked before every block
  /// Called with the progress of the current stage; serialised (never concurrently), possibly
  /// from a worker thread, at most every `min_interval`, and once with done == total when a
  /// stage completes. An exception stops the run and is rethrown by the call.
  std::function<void(const Progress&)> progress;
  std::chrono::milliseconds min_interval{50};  ///< minimum time between two progress calls
  std::size_t block_size = kRunBlockSize;  ///< rays per block (>= 1); results do not depend on it

  /// True if the run is controlled (a token or a callback is set).
  [[nodiscard]] bool active() const noexcept {
    return cancel.has_value() || static_cast<bool>(progress);
  }
};

/// Thrown at the API boundary when a computation stopped after a cancellation request. The
/// output of the call (e.g. the RayBatch of a trace) is then unspecified.
class Cancelled : public std::runtime_error {
 public:
  Cancelled() : std::runtime_error("cancelled") {}
};

/// Bookkeeping of one controlled stage, for the implementations of rtt-trace and rtt-analysis.
/// stop() and add() may be called from worker threads and never throw; finish() is called by
/// the calling thread after the parallel part and throws there (rule 3: no exceptions inside
/// the loops).
class RunMonitor {
 public:
  /// @param control the run control of the call (must outlive the monitor)
  /// @param total   work items of the stage
  /// @param stage   name of the stage (must outlive the monitor; a string literal)
  RunMonitor(const RunControl& control, std::size_t total, std::string_view stage);

  /// True if no new block should start: cancellation requested or the callback failed.
  [[nodiscard]] bool stop() const noexcept;

  /// Counts `n` finished work items and may report progress (serialised, throttled). An
  /// exception of the callback is stored and makes stop() true.
  void add(std::size_t n) noexcept;

  /// After the parallel part: rethrows a stored callback exception, throws Cancelled after a
  /// cancellation request, otherwise reports done == total once (an exception of that call
  /// propagates).
  void finish();

 private:
  void report() noexcept;

  const RunControl& control_;
  std::size_t total_;
  std::string_view stage_;
  std::atomic<std::size_t> done_{0};
  std::atomic<bool> failed_{false};
  std::mutex mutex_;          ///< serialises the callback; guards last_ and error_
  std::exception_ptr error_;  ///< first exception of the callback
  std::chrono::steady_clock::time_point last_;  ///< time of the last report
  bool reported_ = false;      ///< a report was made (the first one is not throttled)
  std::size_t last_done_ = 0;  ///< done of the last report
};

}  // namespace rtt::trace
