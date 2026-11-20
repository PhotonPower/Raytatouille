// rtt trace and rtt analyze against the C++ API (ADR 0052, point 8): every command writes byte for
// byte rtt::io::to_result_json of the same API call; CSV reads back to the same values; wrong
// calls exit 2, failed computations 1, each with a message. The command line runs in this
// process through rtt::cli::run, so the test runs in every CI job, with or without Python.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "cli.hpp"
#include "rtt/analysis/nominal_reference.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/psf.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/analysis/zernike.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/io/result_json.hpp"
#include "rtt/io/result_names.hpp"
#include "rtt/material/material.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace fs = std::filesystem;
namespace trace = rtt::trace;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;

namespace {

const std::string kRef = RTT_REFERENCE_DIR;
const std::string kSchott = std::string(RTT_CATALOG_DIR) + "/schott.agf";
const std::string kSinglet = kRef + "/m1/singlet_const.rtt.json";
const std::string kAchromat = kRef + "/m2/achromat.rtt.json";
const std::string kParaboloid = kRef + "/m2/paraboloid_stop.rtt.json";

struct Run {
  int code = 0;
  std::string out;
  std::string err;
};

Run run(const std::vector<std::string>& args) {
  std::ostringstream out;
  std::ostringstream err;
  Run r;
  r.code = rtt::cli::run(args, out, err);
  r.out = out.str();
  r.err = err.str();
  return r;
}

CompiledSystem compiled(const std::string& file, bool schott = false) {
  rtt::material::MaterialLibrary lib;
  if (schott) lib.add_catalog(kSchott);
  return rtt::compile::compile(rtt::io::load_system(file), lib);
}

/// The CLI output equals `expected` (exit 0, JSON with content).
void same(const Run& r, const std::string& expected) {
  INFO(r.err);
  REQUIRE(r.code == 0);
  CHECK(r.out == expected);
  CHECK(r.out.size() > 100);  // content guard
}

}  // namespace

TEST_CASE("rtt analyze writes the result of the API call (ADR 0052, point 8)", "[cli]") {
  const CompiledSystem singlet = compiled(kSinglet);
  const CompiledSystem achromat = compiled(kAchromat, true);
  const CompiledSystem paraboloid = compiled(kParaboloid);
  const PathId path{0};
  SECTION("first-order") {
    same(run({"analyze", "first-order", kAchromat, "--catalog", kSchott}),
         rtt::io::to_result_json(
             rtt::paraxial::first_order(achromat, path, achromat.reference_wavelength())));
  }
  SECTION("seidel with a pair") {
    same(run({"analyze", "seidel", kAchromat, "--catalog", kSchott, "--pair", "0,2"}),
         rtt::io::to_result_json(rtt::paraxial::seidel(
             achromat, path, achromat.reference_wavelength(), rtt::paraxial::ChromaticPair{0, 2})));
  }
  SECTION("spot, polychromatic by default, with a sampling") {
    rtt::analysis::SpotOptions o;
    o.sampling = rtt::trace::RandomPupil{300, 7};
    same(run({"analyze", "spot", kAchromat, "--catalog", kSchott, "--field", "2", "--rays",
              "random:300:7"}),
         rtt::io::to_result_json(rtt::analysis::spot(achromat, path, 2, std::nullopt, o)));
  }
  SECTION("opd with grid and focus shift") {
    rtt::analysis::OpdOptions o;
    o.grid = 17;
    o.focus_shift = 0.25;
    same(run({"analyze", "opd", kSinglet, "--field", "2", "--grid", "17", "--focus-shift", "0.25"}),
         rtt::io::to_result_json(
             rtt::analysis::opd_map(singlet, path, 2, singlet.reference_wavelength(), o)));
  }
  SECTION("opd and spot with the reference of a nominal system (ADR 0041, point 3)") {
    // The grid sag mirror has no paraxial model; it is aimed with the reference of the Zernike
    // mirror (a different system with a paraxial model). Without --reference it fails.
    const std::string trial_file = kRef + "/m8/grid_sag_mirror.rtt.json";
    const std::string nominal_file = kRef + "/m8/zernike_mirror.rtt.json";
    const CompiledSystem trial = compiled(trial_file);
    const CompiledSystem nominal = compiled(nominal_file);
    const rtt::analysis::NominalReference ref = rtt::analysis::nominal_reference(nominal, path);
    rtt::analysis::OpdOptions o;
    o.grid = 9;
    o.reference = ref;
    same(run({"analyze", "opd", trial_file, "--grid", "9", "--reference", nominal_file}),
         rtt::io::to_result_json(
             rtt::analysis::opd_map(trial, path, 0, trial.reference_wavelength(), o)));
    rtt::analysis::SpotOptions so;
    so.reference = ref;
    same(run({"analyze", "spot", trial_file, "--reference", nominal_file}),
         rtt::io::to_result_json(rtt::analysis::spot(trial, path, 0, std::nullopt, so)));
    CHECK(run({"analyze", "opd", trial_file, "--grid", "9"}).code == 1);  // content guard
  }
  SECTION("zernike: pupil function, then the fit") {
    rtt::analysis::PupilFunctionOptions o;
    o.grid = 24;
    rtt::analysis::ZernikeOptions z;
    z.terms = 15;
    same(run({"analyze", "zernike", kSinglet, "--field", "2", "--grid", "24", "--terms", "15"}),
         rtt::io::to_result_json(rtt::analysis::zernike_fit(
             rtt::analysis::pupil_function(singlet, path, 2, singlet.reference_wavelength(), o),
             z)));
  }
  SECTION("psf") {
    rtt::analysis::PsfOptions o;
    o.pupil.grid = 32;
    o.padding = 2;
    same(run({"analyze", "psf", kParaboloid, "--grid", "32", "--padding", "2"}),
         rtt::io::to_result_json(rtt::analysis::psf(paraboloid, path, 0, o)));
  }
  SECTION("mtf with map and complex OTF") {
    rtt::analysis::MtfOptions o;
    o.pupil.grid = 32;
    o.outputs.map = true;
    o.outputs.complex_otf = true;
    same(run({"analyze", "mtf", kParaboloid, "--grid", "32", "--map", "--complex-otf"}),
         rtt::io::to_result_json(rtt::analysis::mtf(paraboloid, path, 0, o)));
  }
}

TEST_CASE("rtt trace writes the traced batch, the paths and the stats", "[cli]") {
  const CompiledSystem singlet = compiled(kSinglet);
  const PathId path{0};
  std::vector<std::uint16_t> fields(singlet.fields().points.size());
  std::iota(fields.begin(), fields.end(), std::uint16_t{0});
  trace::RayBatch rays =
      trace::make_rays(singlet, path, fields, singlet.reference_wavelength(), trace::GridPupil{5});
  trace::RayBatch again = rays;
  trace::RayPaths paths;
  const trace::TraceStats stats = trace::SequentialTracer().trace(singlet, path, rays);
  static_cast<void>(trace::SequentialTracer().trace(singlet, path, again, paths));
  REQUIRE(stats.count(trace::RayStatus::Alive) == rays.size());  // content guard
  REQUIRE(rays.size() > 10);
  same(run({"trace", kSinglet, "--rays", "grid:5"}), rtt::io::to_result_json(rays));
  same(run({"trace", kSinglet, "--rays", "grid:5", "--paths"}), rtt::io::to_result_json(paths));
  const fs::path dir = fs::temp_directory_path() / "rtt_cli_tests";
  fs::create_directories(dir);
  const fs::path out = dir / "rays.result.json";
  const fs::path st = dir / "stats.result.json";
  const Run r =
      run({"trace", kSinglet, "--rays", "grid:5", "-o", out.string(), "--stats", st.string()});
  REQUIRE(r.code == 0);
  CHECK(r.out.empty());
  const auto text = [](const fs::path& f) {
    // rdbuf, not istreambuf_iterator: that trips a false -Wnull-dereference in GCC 13 at -O3.
    const std::ifstream in(f, std::ios::binary);
    std::ostringstream s;
    s << in.rdbuf();
    return s.str();
  };
  CHECK(text(out) == rtt::io::to_result_json(rays));
  CHECK(text(st) == rtt::io::to_result_json(stats));
  CHECK_FALSE(fs::exists(fs::path(out.string() + ".rtt-partial")));
}

/// The rows of a CSV text: the header, then the fields of each row.
std::vector<std::vector<std::string>> rows(const std::string& text) {
  std::vector<std::vector<std::string>> out;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line)) {
    std::vector<std::string> fields;
    std::string field;
    std::istringstream l(line);
    while (std::getline(l, field, ',')) fields.push_back(field);
    if (!line.empty() && line.back() == ',') fields.emplace_back();
    out.push_back(std::move(fields));
  }
  return out;
}

/// A CSV number back as a double: empty is NaN; compared bit for bit with `expected`.
bool same_value(const std::string& field, double expected) {
  if (field.empty()) return std::isnan(expected);
  return std::stod(field) == expected;
}

TEST_CASE("CSV reads back to the values of the API (ADR 0052, point 5)", "[cli]") {
  const CompiledSystem singlet = compiled(kSinglet);
  const CompiledSystem achromat = compiled(kAchromat, true);
  const CompiledSystem paraboloid = compiled(kParaboloid);
  const PathId path{0};
  SECTION("spot, monochromatic and polychromatic") {
    rtt::analysis::SpotOptions o;
    o.sampling = trace::HexapolarPupil{3};
    const rtt::analysis::SpotDiagram mono =
        rtt::analysis::spot(singlet, path, 2, std::uint16_t{0}, o);
    const auto m = rows(run({"analyze", "spot", kSinglet, "--field", "2", "--wavelength", "0",
                             "--rays", "hexapolar:3", "--format", "csv"})
                            .out);
    REQUIRE(m.size() == mono.points.size() + 1);
    CHECK(m[0] == std::vector<std::string>{"x", "y"});
    for (std::size_t i = 0; i < mono.points.size(); ++i) {
      CHECK(same_value(m[i + 1][0], mono.points[i].x));
      CHECK(same_value(m[i + 1][1], mono.points[i].y));
    }
    const rtt::analysis::SpotDiagram poly = rtt::analysis::spot(achromat, path, 2, std::nullopt, o);
    const auto q = rows(run({"analyze", "spot", kAchromat, "--catalog", kSchott, "--field", "2",
                             "--rays", "hexapolar:3", "--format", "csv"})
                            .out);
    REQUIRE(q.size() == poly.points.size() + 1);
    CHECK(q[0] == std::vector<std::string>{"x", "y", "wavelength", "weight"});
    std::set<std::string> wavelengths;
    for (std::size_t i = 0; i < poly.points.size(); ++i) {
      CHECK(same_value(q[i + 1][0], poly.points[i].x));
      CHECK(same_value(q[i + 1][1], poly.points[i].y));
      CHECK(q[i + 1][2] == std::to_string(poly.points[i].wavelength));
      CHECK(same_value(q[i + 1][3], poly.points[i].weight));
      wavelengths.insert(q[i + 1][2]);
    }
    CHECK(wavelengths.size() > 1);  // content guard: several wavelengths
  }
  SECTION("mtf with the phase of the OTF") {
    rtt::analysis::MtfOptions o;
    o.pupil.grid = 16;
    o.outputs.complex_otf = true;
    const rtt::analysis::Mtf mtf = rtt::analysis::mtf(paraboloid, path, 0, o);
    const auto r = rows(
        run({"analyze", "mtf", kParaboloid, "--grid", "16", "--complex-otf", "--format", "csv"})
            .out);
    CHECK(r[0] == std::vector<std::string>{"direction", "frequency", "modulation", "phase"});
    std::size_t row = 1;
    for (const auto& [name, curve] :
         {std::pair{"tangential", &mtf.tangential}, std::pair{"sagittal", &mtf.sagittal}}) {
      for (std::size_t i = 0; i < curve->frequency.size(); ++i, ++row) {
        REQUIRE(row < r.size());
        CHECK(r[row][0] == name);
        CHECK(same_value(r[row][1], curve->frequency[i]));
        CHECK(same_value(r[row][2], curve->modulation[i]));
        CHECK(same_value(r[row][3], std::arg(curve->otf[i])));
      }
    }
    CHECK(row == r.size());
    CHECK(row > 10);  // content guard
  }
  SECTION("trace and trace --paths") {
    const std::uint16_t field0[] = {0};
    trace::RayBatch rays = trace::make_rays(paraboloid, path, field0,
                                            paraboloid.reference_wavelength(), trace::GridPupil{5});
    trace::RayBatch again = rays;
    trace::RayPaths paths;
    static_cast<void>(trace::SequentialTracer().trace(paraboloid, path, rays));
    static_cast<void>(trace::SequentialTracer().trace(paraboloid, path, again, paths));
    const auto b = rows(
        run({"trace", kParaboloid, "--field", "0", "--rays", "grid:5", "--format", "csv"}).out);
    REQUIRE(b.size() == rays.size() + 1);
    CHECK(b[0].size() == 18);
    for (std::size_t i = 0; i < rays.size(); ++i) {
      const std::vector<std::string>& f = b[i + 1];
      CHECK(same_value(f[0], rays.pos_x()[i]));
      CHECK(same_value(f[2], rays.pos_z()[i]));
      CHECK(same_value(f[5], rays.dir_z()[i]));
      CHECK(f[6] == std::to_string(rays.wl()[i]));
      CHECK(same_value(f[7], rays.opl()[i]));
      CHECK(f[12] == std::to_string(rays.last_surface()[i]));
      CHECK(f[13] == rtt::io::result_name(rays.status()[i]));
      CHECK(same_value(f[17], rays.mode_index()[i]));
    }
    const auto r = rows(run({"trace", kParaboloid, "--field", "0", "--rays", "grid:5", "--paths",
                             "--format", "csv"})
                            .out);
    REQUIRE(r.size() == paths.ray_count() * paths.slots + 1);
    CHECK(r[0] == std::vector<std::string>{"ray", "slot", "status", "x", "y", "z", "dx", "dy", "dz",
                                           "opl"});
    for (std::size_t k = 0; k + 1 < r.size(); ++k) {
      const std::vector<std::string>& f = r[k + 1];
      CHECK(f[0] == std::to_string(paths.ray_indices[k / paths.slots]));
      CHECK(f[2] == rtt::io::result_name(paths.status[k]));
      for (std::size_t c = 0; c < 3; ++c) {
        CHECK(same_value(f[3 + c], paths.position[3 * k + c]));
        CHECK(same_value(f[6 + c], paths.direction[3 * k + c]));
      }
      CHECK(same_value(f[9], paths.opl[k]));
    }
  }
}

TEST_CASE("exit codes and messages (ADR 0052, point 4)", "[cli]") {
  const auto usage = [](const std::vector<std::string>& args) {
    const Run r = run(args);
    INFO(r.err);
    CHECK(r.code == 2);
    CHECK_FALSE(r.err.empty());
    CHECK(r.out.empty());
  };
  usage({"analyze"});
  usage({"analyze", "bogus", kSinglet});
  usage({"analyze", "spot", kSinglet, "--grid", "3"});  // spot has no --grid
  usage({"analyze", "opd", kSinglet, "--grid"});        // a missing value
  usage({"analyze", "opd", kSinglet, "--grid", "x"});
  usage({"analyze", "opd", kSinglet, "--format", "csv"});  // opd has no CSV
  usage({"analyze", "opd", kSinglet, "--field", "99"});
  usage({"trace", kSinglet, "--rays", "spiral:3"});
  usage({"trace", kSinglet, "--polarization"});  // only with --paths
  usage({"trace"});
  const auto failed = [](const std::vector<std::string>& args, const std::string& text) {
    const Run r = run(args);
    INFO(r.err);
    CHECK(r.code == 1);
    CHECK(r.err.find(text) != std::string::npos);
    CHECK(r.out.empty());
  };
  failed({"analyze", "first-order", kAchromat}, "[material.unknown]");  // no catalogue
  failed({"analyze", "spot", kRef + "/m0/does_not_exist.rtt.json"}, "error");
  failed({"analyze", "spot",
          std::string(RTT_SOURCE_DIR) + "/apps/rtt-cli/tests/two_references.rtt.json"},
         "[wavelengths.reference_count]");
  // A system without a paraxial model: the grid sag mirror of #316.
  failed({"analyze", "first-order", kRef + "/m8/grid_sag_mirror.rtt.json"}, "error");
  // An output that cannot be written: exit 1, no partial file.
  const fs::path missing = fs::temp_directory_path() / "rtt_cli_no_such_dir" / "x.json";
  failed({"analyze", "first-order", kSinglet, "-o", missing.string()}, "cannot write");
  CHECK_FALSE(fs::exists(missing.parent_path()));
  // The warnings of an analysis go to stderr with the file, also with CSV (ADR 0052, point 4).
  const Run w = run({"analyze", "mtf", kParaboloid, "--grid", "32", "--format", "csv"});
  CHECK(w.code == 0);
  CHECK(w.err.find(kParaboloid + ": warning [") != std::string::npos);
  // --help of every kind and of trace names every option of the parser and exits 0.
  const auto every_option = [](const std::vector<std::string>& args,
                               const std::vector<std::string>& options) {
    const Run h = run(args);
    CHECK(h.code == 0);
    for (const std::string& option : options) {
      INFO(option);
      CHECK(h.out.find(option) != std::string::npos);
    }
  };
  const std::vector<std::string> common = {"--catalog", "--coatings",   "--config", "--path",
                                           "--threads", "-o, --output", "--format", "--help"};
  const auto with = [&](std::vector<std::string> own) {
    own.insert(own.end(), common.begin(), common.end());
    return own;
  };
  every_option({"trace", "--help"},
               with({"--field", "--wavelength", "--rays", "--aiming", "--paths",
                     "--max-recorded-rays", "--polarization", "--stats"}));
  every_option({"trace", "-h"}, with({"--rays"}));
  every_option({"analyze", "first-order", "--help"}, with({"--wavelength"}));
  every_option({"analyze", "seidel", "--help"}, with({"--wavelength", "--pair"}));
  every_option({"analyze", "spot", "--help"},
               with({"--field", "--wavelength", "--aiming", "--reference", "--rays"}));
  every_option({"analyze", "opd", "--help"}, with({"--field", "--wavelength", "--aiming",
                                                   "--reference", "--focus-shift", "--grid"}));
  every_option({"analyze", "zernike", "--help"},
               with({"--field", "--wavelength", "--aiming", "--focus-shift", "--apodization",
                     "--grid", "--terms"}));
  every_option({"analyze", "psf", "--help"},
               with({"--field", "--wavelengths", "--aiming", "--focus-shift", "--apodization",
                     "--grid", "--padding"}));
  every_option({"analyze", "mtf", "--help"},
               with({"--field", "--wavelengths", "--aiming", "--focus-shift", "--apodization",
                     "--grid", "--map", "--complex-otf"}));
  CHECK(run({"analyze", "--help"}).out.find("zernike") != std::string::npos);
}
