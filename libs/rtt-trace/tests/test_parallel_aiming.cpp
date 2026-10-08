#include <oneapi/tbb/info.h>
#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <atomic>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sources.hpp"

// Parallel ray aiming in make_rays (#119): every ray is aimed independently from its own
// paraxial start value, so the rays do not depend on the number of threads, the control or the
// block size (ADR 0004, addendum #119); a cancellation stops the parallel aiming after the
// running blocks and throws Cancelled at the API boundary.

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::trace::Aiming;
using rtt::trace::Cancelled;
using rtt::trace::CancelToken;
using rtt::trace::Progress;
using rtt::trace::RayBatch;
using rtt::trace::RunControl;

namespace {

template <typename T>
bool same_bits(std::span<const T> a, std::span<const T> b) {
  return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size_bytes()) == 0;
}

/// Every column make_rays() sets, bit for bit (PRT and weight as well: make_rays sets them to
/// the identity and 1).
bool identical(const RayBatch& a, const RayBatch& b) {
  if (a.size() != b.size()) return false;
  bool ok = same_bits(a.pos_x(), b.pos_x()) && same_bits(a.pos_y(), b.pos_y()) &&
            same_bits(a.pos_z(), b.pos_z()) && same_bits(a.dir_x(), b.dir_x()) &&
            same_bits(a.dir_y(), b.dir_y()) && same_bits(a.dir_z(), b.dir_z()) &&
            same_bits(a.wl(), b.wl()) && same_bits(a.field(), b.field()) &&
            same_bits(a.pupil_x(), b.pupil_x()) && same_bits(a.pupil_y(), b.pupil_y()) &&
            same_bits(a.opl(), b.opl()) && same_bits(a.weight(), b.weight()) &&
            same_bits(a.status(), b.status()) && same_bits(a.last_surface(), b.last_surface());
  for (std::size_t i = 0; ok && i < a.size(); ++i) ok = a.prt_matrix(i) == b.prt_matrix(i);
  return ok;
}

/// make_rays() in an arena of `threads` workers, with or without a control.
RayBatch aim_in_arena(const CompiledSystem& cs,
                      PathId path,
                      std::span<const std::uint16_t> fields,
                      std::uint16_t wavelength,
                      int threads,
                      const RunControl* control) {
  oneapi::tbb::task_arena arena(threads);
  return arena.execute([&] {
    const rtt::trace::HexapolarPupil sampling{4};
    return control != nullptr
               ? rtt::trace::make_rays(cs, path, fields, wavelength, sampling, Aiming::Real,
                                       *control)
               : rtt::trace::make_rays(cs, path, fields, wavelength, sampling, Aiming::Real);
  });
}

/// All reference systems that compile with the test catalogues, sorted by path.
std::vector<std::filesystem::path> reference_files() {
  std::vector<std::filesystem::path> files;
  for (const auto& entry : std::filesystem::recursive_directory_iterator(RTT_REFERENCE_DIR)) {
    const std::string name = entry.path().filename().string();
    if (entry.is_regular_file() && name.size() > 9 &&
        name.compare(name.size() - 9, 9, ".rtt.json") == 0) {
      files.push_back(entry.path());
    }
  }
  std::sort(files.begin(), files.end());
  return files;
}

}  // namespace

TEST_CASE("make_rays: bitwise the same for every thread count, control and block size (#119)",
          "[sources][aiming][parallel]") {
  rtt::material::MaterialLibrary lib1;
  lib1.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  rtt::material::MaterialLibrary lib2;
  lib2.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf");
  rtt::coating::CoatingLibrary coatings;
  coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/demo.json");
  coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/m3.json");
  const int all = oneapi::tbb::info::default_concurrency();
  std::size_t compared = 0;
  for (const auto& file : reference_files()) {
    std::optional<CompiledSystem> cs;
    for (const rtt::material::MaterialLibrary* lib : {&lib1, &lib2}) {
      try {
        cs.emplace(rtt::compile::compile(rtt::io::load_system(file), *lib, coatings));
        break;
      } catch (const std::exception&) {
      }
    }
    if (!cs) continue;  // e.g. a coating without catalogue prefix (#59)
    std::vector<std::uint16_t> fields(cs->fields().points.size());
    std::iota(fields.begin(), fields.end(), std::uint16_t{0});
    for (std::uint32_t p = 0; p < cs->paths().size(); ++p) {
      for (std::uint16_t w = 0; w < cs->wavelengths_um().size(); ++w) {
        INFO(file.filename().string() << ", path " << p << ", wavelength " << w);
        std::optional<RayBatch> serial;
        try {
          serial.emplace(aim_in_arena(*cs, PathId{p}, fields, w, 1, nullptr));
        } catch (const std::exception&) {
          // No aiming on this path (no stop, ...): the same exception for every variant.
          REQUIRE_THROWS(aim_in_arena(*cs, PathId{p}, fields, w, 4, nullptr));
          continue;
        }
        for (const int threads : {4, all}) {
          REQUIRE(identical(aim_in_arena(*cs, PathId{p}, fields, w, threads, nullptr), *serial));
        }
        for (const std::size_t block : {std::size_t{7}, std::size_t{256}}) {
          RunControl control;
          control.cancel = CancelToken();
          control.progress = [](const Progress&) {};
          control.block_size = block;
          for (const int threads : {1, 4}) {
            REQUIRE(identical(aim_in_arena(*cs, PathId{p}, fields, w, threads, &control), *serial));
          }
        }
        compared += serial->size();
      }
    }
  }
  // The reference systems were actually aimed: 4758 rays per variant when this was written.
  REQUIRE(compared > 4000);
}

TEST_CASE("make_rays: a cancellation from a second thread stops the parallel aiming (#119)",
          "[sources][aiming][parallel]") {
  // As in the trace (#83): the request comes from another thread, started and joined inside the
  // first progress report, so it certainly falls into the running aiming. At most one block
  // per worker runs after it: 256 real aimings of about 20 us each, about 5 ms (documented
  // latency of the stage aim); 1 s leaves room for a loaded CI runner.
  const rtt::material::MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json"), lib);
  const std::vector<std::uint16_t> fields{0, 1, 2};
  RunControl control;
  control.cancel = CancelToken();
  control.min_interval = std::chrono::milliseconds{0};
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
    // 3 * 10 981 rays: about 0.6 s of serial real aiming, far more than the reported latency.
    [[maybe_unused]] const RayBatch rays = rtt::trace::make_rays(
        cs, PathId{0}, fields, 0, rtt::trace::HexapolarPupil{60}, Aiming::Real, control);
  } catch (const Cancelled&) {
    cancelled = true;
  }
  const auto returned = std::chrono::steady_clock::now();
  REQUIRE(started.load());
  REQUIRE(cancelled);
  REQUIRE(returned - requested < std::chrono::seconds{1});
}
