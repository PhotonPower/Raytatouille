// Cases of the C++ result writer (ADR 0052, point 1): the same calls as CASES in
// test_result_json.py, each result written with rtt::io::to_result_json.

#include "result_json_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/psf.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/analysis/zernike.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/io/result_json.hpp"
#include "rtt/material/material.hpp"
#include "rtt/optim/optimize.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace fs = std::filesystem;

namespace rtt::py::reference {

namespace {

void write(const fs::path& out, const std::string& name, const std::string& text) {
  std::ofstream file(out / ("json_" + name + ".result.json"), std::ios::binary);
  if (!file) throw std::runtime_error("cannot write the case " + name);
  file << text;
}

}  // namespace

void run_result_json_cases(const fs::path& reference_dir,
                           const fs::path& catalog_dir,
                           const fs::path& out,
                           int threads) {
  using compile::CompiledSystem;
  using compile::PathId;
  const material::MaterialLibrary plain;
  material::MaterialLibrary schott;
  schott.add_catalog(catalog_dir / "schott.agf");
  const CompiledSystem singlet =
      compile::compile(io::load_system(reference_dir / "m1/singlet_const.rtt.json"), plain);
  const CompiledSystem achromat =
      compile::compile(io::load_system(reference_dir / "m2/achromat.rtt.json"), schott);
  const CompiledSystem paraboloid =
      compile::compile(io::load_system(reference_dir / "m2/paraboloid_stop.rtt.json"), plain);
  const model::System optim_singlet = io::load_system(reference_dir / "m5/singlet_optim.rtt.json");
  const PathId path{0};

  oneapi::tbb::task_arena arena(threads);
  arena.execute([&] {
    // rays: all fields, hexapolar 2 rings, reference wavelength; traced with all paths recorded.
    std::vector<std::uint16_t> fields(singlet.fields().points.size());
    std::iota(fields.begin(), fields.end(), std::uint16_t{0});
    trace::RayBatch rays = trace::make_rays(singlet, path, fields, singlet.reference_wavelength(),
                                            trace::HexapolarPupil{2});
    trace::RayPaths paths;
    const trace::TraceStats stats = trace::SequentialTracer().trace(singlet, path, rays, paths);
    write(out, "rays", io::to_result_json(rays));
    write(out, "paths", io::to_result_json(paths));
    write(out, "stats", io::to_result_json(stats));

    write(
        out, "first_order",
        io::to_result_json(paraxial::first_order(achromat, path, achromat.reference_wavelength())));
    write(out, "seidel",
          io::to_result_json(paraxial::seidel(achromat, path, achromat.reference_wavelength(),
                                              paraxial::ChromaticPair{0, 2})));

    analysis::SpotOptions so;
    so.sampling = trace::RandomPupil{300, 7};
    write(out, "spot", io::to_result_json(analysis::spot(achromat, path, 2, std::nullopt, so)));

    analysis::OpdOptions oo;
    oo.grid = 17;
    write(out, "opd",
          io::to_result_json(
              analysis::opd_map(singlet, path, 2, singlet.reference_wavelength(), oo)));

    analysis::PupilFunctionOptions po;
    po.grid = 24;
    write(out, "zernike",
          io::to_result_json(analysis::zernike_fit(
              analysis::pupil_function(singlet, path, 2, singlet.reference_wavelength(), po))));

    analysis::PsfOptions pso;
    pso.pupil.grid = 32;
    write(out, "psf", io::to_result_json(analysis::psf(paraboloid, path, 0, pso)));

    analysis::MtfOptions mo;
    mo.pupil.grid = 32;
    mo.outputs.map = true;
    mo.outputs.complex_otf = true;
    write(out, "mtf", io::to_result_json(analysis::mtf(paraboloid, path, 0, mo)));

    write(out, "optimize", io::to_result_json(optim::optimize(optim_singlet, plain)));
  });
}

}  // namespace rtt::py::reference
