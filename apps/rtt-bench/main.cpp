// rtt-bench: run-time measurement of the interactive performance (#83, G4).
//
// Spot of the Cooke triplet (tests/reference/m2/cooke_triplet.rtt.json) with 1027 rays
// (hexapolar, 18 rings) at the reference wavelength, on axis and at the outer field, with one
// thread and with all threads. Times the aiming (make_rays), the trace and the whole spot
// analysis; median and minimum over --runs runs after a warm-up. --quick (smoke test): one
// thread, axis only.
//
// Output: one JSON line on stdout. With --summary FILE a Markdown table is appended to FILE
// (GitHub job summary); with --github a GitHub warning is printed when the median of a spot
// exceeds warn-factor * target. The program never fails because of the times: the target
// ("spot with about 1000 rays of the Cooke triplet below 10 ms") is a rough guide of the
// maintainer, and CI runners vary.

#include <oneapi/tbb/info.h>
#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace {

constexpr double kTargetMs = 10.0;   // rough guide of the maintainer (#83)
constexpr double kWarnFactor = 3.0;  // warn only on a clear excess
constexpr int kRings = 18;           // 1 + 3 * 18 * 19 = 1027 rays

struct Timing {
  double median_ms = 0.0;
  double min_ms = 0.0;
};

Timing measure(int runs, const std::function<void()>& f) {
  for (int i = 0; i < (runs > 1 ? 3 : 0); ++i) f();  // warm-up, none for a single run
  std::vector<double> ms;
  ms.reserve(static_cast<std::size_t>(runs));
  for (int i = 0; i < runs; ++i) {
    const auto t0 = std::chrono::steady_clock::now();
    f();
    const auto t1 = std::chrono::steady_clock::now();
    ms.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
  }
  std::sort(ms.begin(), ms.end());
  return {ms[ms.size() / 2], ms.front()};
}

struct Result {
  int threads = 1;
  std::uint16_t field = 0;
  std::size_t rays = 0;
  Timing aim;
  Timing trace;
  Timing spot;
};

}  // namespace

int main(int argc, char** argv) {
  int runs = 21;
  std::string summary;
  bool github = false;
  bool quick = false;
  std::string reference_dir = RTT_REFERENCE_DIR;
  std::string catalog_dir = RTT_CATALOG_DIR;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--runs" && i + 1 < argc) {
      runs = std::max(1, std::atoi(argv[++i]));
    } else if (arg == "--summary" && i + 1 < argc) {
      summary = argv[++i];
    } else if (arg == "--github") {
      github = true;
    } else if (arg == "--quick") {
      quick = true;  // smoke test: one thread, axis only
    } else if (arg == "--reference-dir" && i + 1 < argc) {
      reference_dir = argv[++i];
    } else if (arg == "--catalog-dir" && i + 1 < argc) {
      catalog_dir = argv[++i];
    } else {
      std::fprintf(stderr,
                   "usage: rtt_bench [--runs N] [--quick] [--summary FILE] [--github] "
                   "[--reference-dir DIR] [--catalog-dir DIR]\n");
      return 2;
    }
  }

  try {
    rtt::material::MaterialLibrary lib;
    lib.add_catalog(catalog_dir + "/m2/schott.agf");
    const rtt::compile::CompiledSystem cs = rtt::compile::compile(
        rtt::io::load_system(reference_dir + "/m2/cooke_triplet.rtt.json"), lib);
    const rtt::compile::PathId path{0};
    const std::uint16_t wl = cs.reference_wavelength();
    const auto outer = static_cast<std::uint16_t>(cs.fields().points.size() - 1);
    const int all = oneapi::tbb::info::default_concurrency();
    std::vector<int> thread_counts{1};
    if (all > 1 && !quick) thread_counts.push_back(all);

    std::vector<Result> results;
    for (const int threads : thread_counts) {
      oneapi::tbb::task_arena arena(threads);
      std::vector<std::uint16_t> field_list{0};
      if (!quick) field_list.push_back(outer);
      for (const std::uint16_t field : field_list) {
        Result r;
        r.threads = threads;
        r.field = field;
        const std::vector<std::uint16_t> fields{field};
        arena.execute([&] {
          rtt::trace::RayBatch rays = rtt::trace::make_rays(
              cs, path, fields, wl, rtt::trace::HexapolarPupil{kRings}, rtt::trace::Aiming::Real);
          r.rays = rays.size();
          r.aim = measure(runs, [&] {
            rays = rtt::trace::make_rays(cs, path, fields, wl, rtt::trace::HexapolarPupil{kRings},
                                         rtt::trace::Aiming::Real);
          });
          const rtt::trace::RayBatch start = rays;
          r.trace = measure(runs, [&] {
            rtt::trace::RayBatch batch = start;
            [[maybe_unused]] const auto stats =
                rtt::trace::SequentialTracer().trace(cs, path, batch);
          });
          const rtt::analysis::SpotOptions options{rtt::trace::HexapolarPupil{kRings},
                                                   rtt::trace::Aiming::Real};
          r.spot = measure(runs, [&] {
            [[maybe_unused]] const auto d = rtt::analysis::spot(cs, path, field, wl, options);
          });
        });
        results.push_back(r);
      }
    }

    // JSON line.
    std::printf(
        "{\"benchmark\": \"cooke_triplet_spot\", \"runs\": %d, \"cpus\": %d, "
        "\"target_ms\": %.1f, \"warn_factor\": %.1f, \"results\": [",
        runs, all, kTargetMs, kWarnFactor);
    for (std::size_t i = 0; i < results.size(); ++i) {
      const Result& r = results[i];
      std::printf(
          "%s{\"threads\": %d, \"field\": %u, \"rays\": %zu, \"aim_ms\": %.3f, "
          "\"trace_ms\": %.3f, \"spot_ms\": %.3f, \"spot_min_ms\": %.3f}",
          i == 0 ? "" : ", ", r.threads, static_cast<unsigned>(r.field), r.rays, r.aim.median_ms,
          r.trace.median_ms, r.spot.median_ms, r.spot.min_ms);
    }
    std::printf("]}\n");

    if (!summary.empty()) {
      std::ofstream md(summary, std::ios::app);
      md << "### Benchmark: spot of the Cooke triplet (" << results.front().rays
         << " rays, median of " << runs << " runs)\n\n"
         << "Rough guide: " << kTargetMs << " ms; warning above " << kWarnFactor * kTargetMs
         << " ms (never an error).\n\n"
         << "| threads | field | aim (ms) | trace (ms) | spot (ms) | spot min (ms) |\n"
         << "| --- | --- | --- | --- | --- | --- |\n";
      for (const Result& r : results) {
        char line[200];
        std::snprintf(line, sizeof line, "| %d | %u | %.3f | %.3f | %.3f | %.3f |\n", r.threads,
                      static_cast<unsigned>(r.field), r.aim.median_ms, r.trace.median_ms,
                      r.spot.median_ms, r.spot.min_ms);
        md << line;
      }
      md << "\n";
    }
    if (github) {
      for (const Result& r : results) {
        if (r.spot.median_ms > kWarnFactor * kTargetMs) {
          std::printf(
              "::warning title=Benchmark::spot of the Cooke triplet (%zu rays, field %u, "
              "%d threads) took %.1f ms (median), more than %.0f x the rough guide of "
              "%.0f ms\n",
              r.rays, static_cast<unsigned>(r.field), r.threads, r.spot.median_ms, kWarnFactor,
              kTargetMs);
        }
      }
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "rtt_bench: %s\n", e.what());
    return 1;
  }
  return 0;
}
