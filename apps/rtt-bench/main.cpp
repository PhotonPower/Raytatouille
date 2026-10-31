// rtt-bench: run-time measurement of the interactive performance (#83, G4).
//
// Spot of the Cooke triplet (tests/reference/m2/cooke_triplet.rtt.json) with 1027 rays
// (hexapolar, 18 rings) at the reference wavelength, on axis and at the outer field, with one
// thread and with all threads. Times the aiming (make_rays), the trace and the whole spot
// analysis; median and minimum over --runs runs after a warm-up. --quick (smoke test): one
// thread, axis only.
//
// A second system, "cooke_triplet_asphere", is the same triplet with every sphere turned into an
// even asphere with the same radius and A4 = kAsphereA4: the spheres use the analytic conic
// intersection, the aspheres the Newton intersection (rtt-geom), whose hot path ADR 0041 (M8 S0)
// changes. The rough guide applies to the Cooke triplet only.
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
#include <variant>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace {

constexpr double kTargetMs = 10.0;   // rough guide of the maintainer (#83)
constexpr double kWarnFactor = 3.0;  // warn only on a clear excess
constexpr int kRings = 18;           // 1 + 3 * 18 * 19 = 1027 rays
// A4 of the asphere variant, mm^-3: about 1e-4 mm of sag at the rim of the lenses (r ~ 6 mm),
// enough for a few Newton steps from the conic start value, small against the design.
constexpr double kAsphereA4 = 1e-7;

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
  const char* system = "";
  int threads = 1;
  std::uint16_t field = 0;
  std::size_t rays = 0;
  Timing aim;
  Timing trace;
  Timing spot;
  std::size_t arrived = 0;  ///< rays of the spot that reached the image (content guard)
};

/// The system with every Conic base shape replaced by an EvenAsphere with the same radius and
/// conic constant and A4 = kAsphereA4 (the Cooke triplet has only top-level elements).
rtt::model::System with_aspheres(rtt::model::System s) {
  for (rtt::model::Node& node : s.root.children) {
    auto* element = std::get_if<rtt::model::Element>(&node.value);
    if (element == nullptr) continue;
    for (rtt::model::Surface& surface : element->surfaces) {
      if (const auto* conic = std::get_if<rtt::model::Conic>(&surface.shape.base)) {
        surface.shape.base =
            rtt::model::EvenAsphere{conic->radius, conic->conic, {rtt::model::Param(kAsphereA4)}};
      }
    }
  }
  return s;
}

/// Times aiming, trace and spot of `cs` for every thread count and field.
std::vector<Result> bench(const char* name,
                          const rtt::compile::CompiledSystem& cs,
                          int runs,
                          bool quick,
                          const std::vector<int>& thread_counts) {
  const rtt::compile::PathId path{0};
  const std::uint16_t wl = cs.reference_wavelength();
  const auto outer = static_cast<std::uint16_t>(cs.fields().points.size() - 1);
  std::vector<Result> results;
  for (const int threads : thread_counts) {
    oneapi::tbb::task_arena arena(threads);
    std::vector<std::uint16_t> field_list{0};
    if (!quick) field_list.push_back(outer);
    for (const std::uint16_t field : field_list) {
      Result r;
      r.system = name;
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
          [[maybe_unused]] const auto stats = rtt::trace::SequentialTracer().trace(cs, path, batch);
        });
        const rtt::analysis::SpotOptions options{rtt::trace::HexapolarPupil{kRings},
                                                 rtt::trace::Aiming::Real};
        r.spot = measure(runs, [&] {
          [[maybe_unused]] const auto d = rtt::analysis::spot(cs, path, field, wl, options);
        });
        r.arrived = rtt::analysis::spot(cs, path, field, wl, options).rays_arrived;
      });
      results.push_back(r);
    }
  }
  return results;
}

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
    const rtt::model::System cooke =
        rtt::io::load_system(reference_dir + "/m2/cooke_triplet.rtt.json");
    const rtt::compile::CompiledSystem cs = rtt::compile::compile(cooke, lib);
    const rtt::compile::CompiledSystem aspheres = rtt::compile::compile(with_aspheres(cooke), lib);
    const int all = oneapi::tbb::info::default_concurrency();
    std::vector<int> thread_counts{1};
    if (all > 1 && !quick) thread_counts.push_back(all);

    std::vector<Result> results = bench("cooke_triplet", cs, runs, quick, thread_counts);
    for (const Result& r : bench("cooke_triplet_asphere", aspheres, runs, quick, thread_counts)) {
      results.push_back(r);
    }

    // JSON line.
    std::printf(
        "{\"benchmark\": \"cooke_triplet_spot\", \"runs\": %d, \"cpus\": %d, "
        "\"target_ms\": %.1f, \"warn_factor\": %.1f, \"results\": [",
        runs, all, kTargetMs, kWarnFactor);
    for (std::size_t i = 0; i < results.size(); ++i) {
      const Result& r = results[i];
      std::printf(
          "%s{\"system\": \"%s\", \"threads\": %d, \"field\": %u, \"rays\": %zu, "
          "\"arrived\": %zu, \"aim_ms\": %.3f, \"trace_ms\": %.3f, \"spot_ms\": %.3f, "
          "\"spot_min_ms\": %.3f}",
          i == 0 ? "" : ", ", r.system, r.threads, static_cast<unsigned>(r.field), r.rays,
          r.arrived, r.aim.median_ms, r.trace.median_ms, r.spot.median_ms, r.spot.min_ms);
    }
    std::printf("]}\n");

    if (!summary.empty()) {
      std::ofstream md(summary, std::ios::app);
      md << "### Benchmark: spot of the Cooke triplet (" << results.front().rays
         << " rays, median of " << runs << " runs)\n\n"
         << "Rough guide: " << kTargetMs << " ms; warning above " << kWarnFactor * kTargetMs
         << " ms (never an error).\n\n"
         << "The variant with even aspheres (Newton intersection) has no guide.\n\n"
         << "| system | threads | field | aim (ms) | trace (ms) | spot (ms) | spot min (ms) |\n"
         << "| --- | --- | --- | --- | --- | --- | --- |\n";
      for (const Result& r : results) {
        char line[240];
        std::snprintf(line, sizeof line, "| %s | %d | %u | %.3f | %.3f | %.3f | %.3f |\n", r.system,
                      r.threads, static_cast<unsigned>(r.field), r.aim.median_ms, r.trace.median_ms,
                      r.spot.median_ms, r.spot.min_ms);
        md << line;
      }
      md << "\n";
    }
    if (github) {
      for (const Result& r : results) {
        if (std::string(r.system) == "cooke_triplet" &&
            r.spot.median_ms > kWarnFactor * kTargetMs) {
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
