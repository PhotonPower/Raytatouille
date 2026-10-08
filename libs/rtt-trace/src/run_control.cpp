#include "rtt/trace/run_control.hpp"

#include <chrono>
#include <cstddef>
#include <exception>
#include <mutex>
#include <string_view>

namespace rtt::trace {

RunMonitor::RunMonitor(const RunControl& control, std::size_t total, std::string_view stage)
    : control_(control), total_(total), stage_(stage) {}

bool RunMonitor::stop() const noexcept {
  return failed_.load(std::memory_order_relaxed) ||
         (control_.cancel && control_.cancel->cancelled());
}

void RunMonitor::add(std::size_t n) noexcept {
  done_.fetch_add(n, std::memory_order_relaxed);
  report();
}

void RunMonitor::report() noexcept {
  if (!control_.progress || failed_.load(std::memory_order_relaxed)) return;
  // Another worker reporting right now: skip this report (throttling allows it). The values
  // reported under the lock are read from the counter there, so they never decrease.
  const std::unique_lock lock(mutex_, std::try_to_lock);
  if (!lock.owns_lock()) return;
  // Re-check under the lock: another worker may have failed while this one waited.
  if (failed_.load(std::memory_order_relaxed)) return;
  const auto now = std::chrono::steady_clock::now();
  if (reported_ && now - last_ < control_.min_interval) return;
  const std::size_t current = done_.load(std::memory_order_relaxed);
  try {
    control_.progress(Progress{current, total_, stage_});
  } catch (...) {
    if (!error_) error_ = std::current_exception();
    failed_.store(true, std::memory_order_relaxed);
  }
  last_ = now;
  reported_ = true;
  last_done_ = current;
}

void RunMonitor::finish() {
  {
    const std::lock_guard lock(mutex_);
    if (error_) std::rethrow_exception(error_);
  }
  if (control_.cancel && control_.cancel->cancelled()) throw Cancelled();
  if (!control_.progress) return;
  // The final report of a completed stage, unless the last report already said done == total;
  // an exception here propagates directly.
  const std::lock_guard lock(mutex_);
  if (reported_ && last_done_ == total_) return;
  control_.progress(Progress{total_, total_, stage_});
}

}  // namespace rtt::trace
