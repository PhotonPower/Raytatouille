// The aiming of the optimization as a run option (#291; ADR 0030, addendum #291):
// OptimizeOptions::aiming and the MeritFunction parameter, real by default (bit for bit as
// before, the other tests), paraxial on request. The merit function sees the same bundles as the
// analyses with the same aiming: a spot_rms operand equals analysis::spot, a ray operand the ray
// of make_rays, and the rms_spot generator generator_residuals() with that aiming.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <span>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/perturbation.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/optimization.hpp"
#include "rtt/model/system.hpp"
#include "rtt/optim/generators.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/optimize.hpp"
#include "rtt/optim/overlay.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::PathId;
using rtt::model::System;
using rtt::trace::Aiming;

namespace {

/// m1/two_lenses_stop_between (the stop between two lenses, so real and paraxial aiming differ
/// by the pupil aberration of the first lens) with a spot_rms operand of field 2 and the
/// rms_spot generator.
System singlet() {
  System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/two_lenses_stop_between.rtt.json");
  rtt::model::SpotRmsOperand spot;
  spot.path = "main";
  spot.field = 2;
  s.optimization.operands.push_back(spot);
  rtt::model::SpotGenerator g;
  g.path = "main";
  s.optimization.generators.push_back(g);
  return s;
}

const rtt::material::MaterialLibrary& library() {
  static const rtt::material::MaterialLibrary lib;
  return lib;
}

}  // namespace

TEST_CASE("the merit function aims as the analyses with the same aiming", "[optim][291]") {
  const System s = singlet();
  const rtt::compile::CompiledSystem cs = rtt::compile::compile(s, library());
  double paraxial_value = 0.0;
  double real_value = 0.0;
  for (const Aiming aiming : {Aiming::Real, Aiming::Paraxial}) {
    INFO("aiming " << static_cast<int>(aiming));
    const rtt::optim::MeritFunction merit(s, library(), nullptr, nullptr, aiming);
    REQUIRE(merit.aiming() == aiming);
    const rtt::optim::MeritEvaluation e = merit.evaluate(merit.start());
    // the spot operand is analysis::spot with that aiming (hexapolar, 6 rings, the default)
    rtt::analysis::SpotOptions o;
    o.sampling = rtt::trace::HexapolarPupil{6};
    o.aiming = aiming;
    const auto spot = rtt::analysis::spot(cs, PathId{0}, 2, cs.reference_wavelength(), o);
    CHECK(e.values[0] == spot.stats.rms_centroid);
    // the generator residuals are those of generator_residuals() with that aiming
    std::vector<double> residuals(rtt::optim::generator_size(s.optimization.generators[0], cs));
    rtt::optim::GeneratorStats stats;
    rtt::optim::generator_residuals(s.optimization.generators[0], cs, residuals, stats, aiming);
    REQUIRE(e.residuals.size() == 1 + residuals.size());
    for (std::size_t i = 0; i < residuals.size(); ++i) CHECK(e.residuals[1 + i] == residuals[i]);
    (aiming == Aiming::Real ? real_value : paraxial_value) = e.values[0];
  }
  // content guard: the two aimings aim different bundles at field 2 (pupil aberration)
  CHECK(std::abs(real_value - paraxial_value) > 1e-9);
  // the default is real, bit for bit
  const rtt::optim::MeritFunction plain(s, library(), nullptr);
  CHECK(plain.aiming() == Aiming::Real);
  CHECK(plain.evaluate(plain.start()).values[0] == real_value);
}

TEST_CASE("optimize with paraxial aiming", "[optim][291]") {
  // the two lenses with the image distance variable and the rms_spot generator: the run ends
  // with either aiming, and the image distance it finds depends on the aiming of the spot
  System s = singlet();
  s.optimization.operands.clear();
  std::get<rtt::model::Element>(s.root.children[3].value).pose.position[2].variable = true;
  rtt::optim::OptimizeOptions real;
  rtt::optim::OptimizeOptions paraxial;
  paraxial.aiming = Aiming::Paraxial;
  const rtt::optim::OptimResult a = rtt::optim::optimize(s, library(), nullptr, real);
  const rtt::optim::OptimResult b = rtt::optim::optimize(s, library(), nullptr, paraxial);
  CHECK(a.status != rtt::optim::LmStatus::Failed);
  CHECK(b.status != rtt::optim::LmStatus::Failed);
  REQUIRE(a.variables.size() == b.variables.size());
  bool differs = false;
  for (std::size_t i = 0; i < a.variables.size(); ++i) {
    differs = differs || a.variables[i].end != b.variables[i].end;
  }
  CHECK(differs);
  // the default run is the run with aiming real
  const rtt::optim::OptimResult c = rtt::optim::optimize(s, library());
  CHECK(c.patch == a.patch);
}

TEST_CASE("the parabasal operands of a trial follow the aiming", "[optim][291]") {
  // efl and image_fnumber of a trial (ADR 0030, addendum #256) whose first lens, before the
  // stop, is decentred by 0.5 mm: the chief and parabasal rays aimed real and paraxially are
  // different rays (in the nominal system they coincide, the real aiming starts at the paraxial
  // solution), and both values stay close to the nominal ones.
  System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/two_lenses_stop_between.rtt.json");
  for (const auto q :
       {rtt::model::FirstOrderQuantity::Efl, rtt::model::FirstOrderQuantity::ImageFNumber}) {
    rtt::model::FirstOrderOperand op;
    op.path = "main";
    op.quantity = q;
    s.optimization.operands.push_back(op);
  }
  rtt::compile::Perturbation decentred;
  decentred.nodes["L1"].translation_mm = rtt::math::Vec3(0.0, 0.5, 0.0);
  const rtt::optim::TrialOverlay overlay{decentred,
                                         rtt::optim::nominal_basis(s, library(), nullptr)};
  const auto values = [&](Aiming aiming) {
    const rtt::optim::MeritFunction merit(s, library(), nullptr, &overlay, aiming);
    return merit.evaluate(merit.start()).values;
  };
  const std::vector<double> real = values(Aiming::Real);
  const std::vector<double> paraxial = values(Aiming::Paraxial);
  REQUIRE(real.size() == 2);
  for (std::size_t i = 0; i < 2; ++i) {
    INFO("operand " << i);
    CHECK(real[i] != paraxial[i]);                                      // content guard
    CHECK(std::abs(real[i] - paraxial[i]) < 1e-3 * std::abs(real[i]));  // nearly the same
  }
  CHECK(std::abs(real[0] - 38.2) < 0.05);  // the EFL of the two lenses
}
