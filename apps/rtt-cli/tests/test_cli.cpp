// rtt trace and rtt analyze against the C++ API (ADR 0052, point 8): every command writes byte for
// byte rtt::io::to_result_json of the same API call; CSV reads back to the same values; wrong
// calls exit 2, failed computations 1, each with a message. The command line runs in this
// process through rtt::cli::run, so the test runs in every CI job, with or without Python.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <numeric>
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
  SECTION("opd with the reference of a nominal system") {
    rtt::analysis::OpdOptions o;
    o.grid = 9;
    o.reference = rtt::analysis::nominal_reference(singlet, path);
    same(run({"analyze", "opd", kSinglet, "--grid", "9", "--reference", kSinglet}),
         rtt::io::to_result_json(
             rtt::analysis::opd_map(singlet, path, 0, singlet.reference_wavelength(), o)));
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

TEST_CASE("CSV reads back to the values of the API (ADR 0052, point 5)", "[cli]") {
  const CompiledSystem singlet = compiled(kSinglet);
  rtt::analysis::SpotOptions o;
  o.sampling = trace::HexapolarPupil{3};
  const rtt::analysis::SpotDiagram d =
      rtt::analysis::spot(singlet, PathId{0}, 2, std::uint16_t{0}, o);
  const Run r = run({"analyze", "spot", kSinglet, "--field", "2", "--wavelength", "0", "--rays",
                     "hexapolar:3", "--format", "csv"});
  REQUIRE(r.code == 0);
  std::istringstream in(r.out);
  std::string line;
  std::getline(in, line);
  CHECK(line == "x,y");  // monochromatic
  std::size_t i = 0;
  while (std::getline(in, line)) {
    REQUIRE(i < d.points.size());
    const std::size_t comma = line.find(',');
    CHECK(std::stod(line.substr(0, comma)) == d.points[i].x);
    CHECK(std::stod(line.substr(comma + 1)) == d.points[i].y);
    ++i;
  }
  CHECK(i == d.points.size());
  CHECK(i > 10);  // content guard
  const Run m = run({"analyze", "mtf", kParaboloid, "--grid", "16", "--format", "csv"});
  REQUIRE(m.code == 0);
  CHECK(m.out.starts_with("direction,frequency,modulation\ntangential,0.0,1.0\n"));
  CHECK(m.out.find("\nsagittal,") != std::string::npos);
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
  // A system without a paraxial model: the Zernike mirror of #316 with an astigmatic term.
  failed({"analyze", "first-order", kRef + "/m8/grid_sag_mirror.rtt.json"}, "error");
  const Run h = run({"analyze", "psf", "--help"});
  CHECK(h.code == 0);
  for (const char* option : {"--padding", "--grid", "--wavelengths", "--catalog", "-o, --output"}) {
    CHECK(h.out.find(option) != std::string::npos);
  }
  CHECK(run({"trace", "--help"}).code == 0);
  CHECK(run({"analyze", "--help"}).out.find("zernike") != std::string::npos);
}
