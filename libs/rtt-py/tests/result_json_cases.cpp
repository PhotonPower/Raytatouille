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
#include "rtt/model/model.hpp"
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
  const model::System optim_configurations =
      io::load_system(reference_dir / "m5/singlet_optim_configurations.rtt.json");
  // The paraboloid with its mirror aperture reduced to 21 mm (beam radius 30 mm): vignetting, as
  // in analysis_cases.cpp.
  model::System clipped = io::load_system(reference_dir / "m2/paraboloid_stop.rtt.json");
  std::get<model::Element>(clipped.root.children[1].value).surfaces[0].aperture =
      model::CircularAperture{21.0, 0.0};
  const CompiledSystem paraboloid21 = compile::compile(clipped, plain);
  // Without a stop: no pupils in the first order.
  const CompiledSystem no_stop =
      compile::compile(io::load_system(reference_dir / "m1/paraboloid_mirror.rtt.json"), plain);
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

    // Branches the cases above do not reach (review of PR A):
    // Q in the batch and P and Q per slot.
    trace::RayBatch q_rays = trace::make_rays(singlet, path, fields, singlet.reference_wavelength(),
                                              trace::HexapolarPupil{2});
    q_rays.enable_q();
    trace::RayPaths q_paths;
    q_paths.polarization = true;
    static_cast<void>(trace::SequentialTracer().trace(singlet, path, q_rays, q_paths));
    write(out, "rays_q", io::to_result_json(q_rays));
    write(out, "paths_q", io::to_result_json(q_paths));
    // Lost rays: status values, lost_at, NaN slots, losses with a worst surface.
    const std::uint16_t field0[] = {0};
    trace::RayBatch lost = trace::make_rays(
        paraboloid21, path, field0, paraboloid21.reference_wavelength(), trace::GridPupil{21});
    trace::RayPaths lost_paths;
    const trace::TraceStats lost_stats =
        trace::SequentialTracer().trace(paraboloid21, path, lost, lost_paths);
    write(out, "rays_lost", io::to_result_json(lost));
    write(out, "paths_lost", io::to_result_json(lost_paths));
    write(out, "stats_lost", io::to_result_json(lost_stats));
    analysis::SpotOptions grid21;
    grid21.sampling = trace::GridPupil{21};
    write(out, "spot_lost",
          io::to_result_json(analysis::spot(paraboloid21, path, 0, std::nullopt, grid21)));
    // Seidel without a pair, a monochromatic spot, the MTF with the default outputs.
    write(out, "seidel_plain",
          io::to_result_json(paraxial::seidel(achromat, path, achromat.reference_wavelength())));
    write(out, "spot_mono",
          io::to_result_json(
              analysis::spot(singlet, path, 2, std::uint16_t{0}, analysis::SpotOptions{})));
    analysis::MtfOptions plain_mtf;
    plain_mtf.pupil.grid = 32;
    write(out, "mtf_plain", io::to_result_json(analysis::mtf(paraboloid, path, 0, plain_mtf)));
    // First order without a stop (no pupils).
    write(out, "first_order_no_stop",
          io::to_result_json(paraxial::first_order(no_stop, path, no_stop.reference_wavelength())));
    // A variable row per configuration (VariableValue.row and configuration set).
    write(out, "optimize_configurations",
          io::to_result_json(optim::optimize(optim_configurations, plain)));
  });
}

}  // namespace rtt::py::reference
