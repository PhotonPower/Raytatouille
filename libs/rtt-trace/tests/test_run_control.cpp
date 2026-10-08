#include <oneapi/tbb/task_arena.h>

#include <atomic>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

// Cancellation and progress (#83, run_control.hpp): results never depend on the control, the
// block size or the thread count; a cancellation stops after the running blocks and throws
// Cancelled at the API boundary; exceptions of the progress callback are rethrown there.

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::trace::Aiming;
using rtt::trace::Cancelled;
using rtt::trace::CancelToken;
using rtt::trace::Progress;
using rtt::trace::RayBatch;
using rtt::trace::RunControl;
using rtt::trace::SequentialTracer;

namespace {

CompiledSystem singlet() {
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json"), lib);
}

/// Hexapolar rays of all three fields of the singlet (1 + 3 r (r + 1) per field).
RayBatch rays_of(const CompiledSystem& cs, int rings, Aiming aiming) {
  const std::vector<std::uint16_t> fields{0, 1, 2};
  return rtt::trace::make_rays(cs, PathId{0}, fields, 0, rtt::trace::HexapolarPupil{rings}, aiming);
}

template <typename T>
bool same_bits(std::span<const T> a, std::span<const T> b) {
  return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size_bytes()) == 0;
}

/// Every column of two batches, bit for bit.
bool identical(const RayBatch& a, const RayBatch& b) {
  if (a.size() != b.size()) return false;
  bool ok = same_bits(a.pos_x(), b.pos_x()) && same_bits(a.pos_y(), b.pos_y()) &&
            same_bits(a.pos_z(), b.pos_z()) && same_bits(a.dir_x(), b.dir_x()) &&
            same_bits(a.dir_y(), b.dir_y()) && same_bits(a.dir_z(), b.dir_z()) &&
            same_bits(a.opl(), b.opl()) && same_bits(a.weight(), b.weight()) &&
            same_bits(a.status(), b.status()) && same_bits(a.last_surface(), b.last_surface());
  for (std::size_t i = 0; ok && i < a.size(); ++i) ok = a.prt_matrix(i) == b.prt_matrix(i);
  return ok;
}

}  // namespace

TEST_CASE("run control: results do not depend on the control, block size or threads (#83)",
          "[run_control]") {
  const CompiledSystem cs = singlet();
  const RayBatch start = rays_of(cs, 20, Aiming::Real);  // 3 * 1261 rays
  RayBatch reference = start;
  [[maybe_unused]] const auto s0 = SequentialTracer().trace(cs, PathId{0}, reference);

  for (const std::size_t block : {std::size_t{256}, std::size_t{7}}) {
    for (const int threads : {1, 4}) {
      INFO("block " << block << ", threads " << threads);
      RunControl control;
      control.cancel = CancelToken();
      control.progress = [](const Progress&) {};
      control.block_size = block;
      RayBatch rays = start;
      oneapi::tbb::task_arena arena(threads);
      const auto stats =
          arena.execute([&] { return SequentialTracer().trace(cs, PathId{0}, rays, control); });
      REQUIRE(identical(rays, reference));
      REQUIRE(stats.rays == s0.rays);
    }
  }

  SECTION("make_rays with a control gives the same rays") {
    RunControl control;
    control.progress = [](const Progress&) {};
    control.block_size = 5;
    const std::vector<std::uint16_t> fields{0, 1, 2};
    const RayBatch aimed = rtt::trace::make_rays(
        cs, PathId{0}, fields, 0, rtt::trace::HexapolarPupil{20}, Aiming::Real, control);
    REQUIRE(identical(aimed, start));
  }
}

TEST_CASE("run control: progress is monotone and ends with done == total (#83)", "[run_control]") {
  const CompiledSystem cs = singlet();
  RayBatch rays = rays_of(cs, 20, Aiming::Paraxial);
  std::vector<Progress> calls;
  std::vector<std::string> stages;
  RunControl control;
  control.min_interval = std::chrono::milliseconds{0};
  control.block_size = 100;
  control.progress = [&](const Progress& p) {
    calls.push_back(p);
    stages.emplace_back(p.stage);
  };
  [[maybe_unused]] const auto stats = SequentialTracer().trace(cs, PathId{0}, rays, control);
  REQUIRE(calls.size() >= 2);
  for (std::size_t k = 0; k < calls.size(); ++k) {
    REQUIRE(stages[k] == "trace");
    REQUIRE(calls[k].total == rays.size());
    if (k > 0) REQUIRE(calls[k].done >= calls[k - 1].done);
  }
  REQUIRE(calls.back().done == rays.size());

  SECTION("throttled: with a long minimum interval only the first and the final call") {
    calls.clear();
    stages.clear();
    control.min_interval = std::chrono::hours{1};
    RayBatch again = rays_of(cs, 20, Aiming::Paraxial);
    [[maybe_unused]] const auto s = SequentialTracer().trace(cs, PathId{0}, again, control);
    REQUIRE(calls.size() <= 2);
    REQUIRE(calls.back().done == again.size());
  }

  SECTION("make_rays reports the stage aim") {
    std::vector<Progress> aim_calls;
    std::vector<std::string> aim_stages;
    RunControl aim_control;
    aim_control.min_interval = std::chrono::milliseconds{0};
    aim_control.block_size = 50;
    aim_control.progress = [&](const Progress& p) {
      aim_calls.push_back(p);
      aim_stages.emplace_back(p.stage);
    };
    const std::vector<std::uint16_t> fields{0, 1};
    const RayBatch aimed = rtt::trace::make_rays(
        cs, PathId{0}, fields, 0, rtt::trace::HexapolarPupil{10}, Aiming::Real, aim_control);
    REQUIRE(aim_calls.size() >= 2);
    for (const auto& s : aim_stages) REQUIRE(s == "aim");
    REQUIRE(aim_calls.back().done == aimed.size());
    REQUIRE(aim_calls.back().total == aimed.size());
  }
}

TEST_CASE("run control: a cancelled call throws Cancelled (#83)", "[run_control]") {
  const CompiledSystem cs = singlet();

  SECTION("cancelled before the start") {
    RunControl control;
    control.cancel = CancelToken();
    control.cancel->request_cancel();
    RayBatch rays = rays_of(cs, 10, Aiming::Paraxial);
    REQUIRE_THROWS_AS(SequentialTracer().trace(cs, PathId{0}, rays, control), Cancelled);
    const std::vector<std::uint16_t> fields{0};
    REQUIRE_THROWS_AS(rtt::trace::make_rays(cs, PathId{0}, fields, 0,
                                            rtt::trace::HexapolarPupil{10}, Aiming::Real, control),
                      Cancelled);
  }

  SECTION("cancelled from a second thread returns within the documented time") {
    // At most one block per worker runs after the request: 256 rays of the singlet take well
    // below a millisecond, so 1 s leaves room for a loaded CI runner. The batch is large
    // (3 * 109 621 rays) so that cancelling saves a measurable amount of work.
    RayBatch rays = rays_of(cs, 191, Aiming::Paraxial);
    RunControl control;
    control.cancel = CancelToken();
    control.min_interval = std::chrono::milliseconds{0};
    // The request comes from another thread, started and joined inside the first progress
    // report: so it certainly falls into the parallel part (the trace cannot finish first).
    std::atomic<bool> started{false};
    std::chrono::steady_clock::time_point requested;
    control.progress = [&](const Progress&) {
      if (started.exchange(true)) return;
      std::thread canceller([&] {
        requested = std::chrono::steady_clock::now();
        control.cancel->request_cancel();
      });
      canceller.join();
    };
    bool cancelled = false;
    try {
      [[maybe_unused]] const auto s = SequentialTracer().trace(cs, PathId{0}, rays, control);
    } catch (const Cancelled&) {
      cancelled = true;
    }
    const auto returned = std::chrono::steady_clock::now();
    REQUIRE(started.load());
    REQUIRE(cancelled);
    REQUIRE(returned - requested < std::chrono::seconds{1});
  }
}

TEST_CASE("run control: an exception of the progress callback is rethrown after the run (#83)",
          "[run_control]") {
  const CompiledSystem cs = singlet();
  RayBatch rays = rays_of(cs, 20, Aiming::Paraxial);
  RunControl control;
  control.min_interval = std::chrono::milliseconds{0};
  control.block_size = 64;
  int calls = 0;
  control.progress = [&](const Progress&) {
    ++calls;
    throw std::domain_error("stop from the callback");
  };
  try {
    [[maybe_unused]] const auto s = SequentialTracer().trace(cs, PathId{0}, rays, control);
    FAIL("no exception");
  } catch (const std::domain_error& e) {
    REQUIRE(std::string(e.what()) == "stop from the callback");
  }
  REQUIRE(calls == 1);  // the run stops at the first failing call
}
