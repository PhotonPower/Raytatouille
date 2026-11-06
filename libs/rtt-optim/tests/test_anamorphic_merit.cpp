// Optimization and the overlay of a tolerance trial on anamorphic paths (#315, ADR 0042,
// point 4; acceptance 7): the path check of the merit function (ADR 0030, point 3 and addendum
// #256) lets operands from real rays run on an anamorphic path whose aiming context per section
// builds (ADR 0042, point 3); the paraxial operands and everything that needs an exit pupil stay
// merit.operand_unsupported. OPD of an anamorphic system is a ParaxialError, also with the
// reference of the (anamorphic) nominal system.
//
// The system is the plano-convex cylinder lens of acceptance 1 (f_y = 100 mm, afocal in x).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/nominal_reference.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/errors.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/overlay.hpp"
#include "rtt/paraxial/paraxial.hpp"

using Catch::Matchers::ContainsSubstring;
using rtt::model::System;
using rtt::optim::MeritEvaluation;
using rtt::optim::MeritFunction;
using rtt::optim::OptimError;
using rtt::optim::TrialOverlay;

namespace {

const rtt::material::MaterialLibrary& library() {
  static const rtt::material::MaterialLibrary materials;
  return materials;
}

/// The cylinder lens of ADR 0042, acceptance 1: stop r = 10 at z = 0, cylinder R = 51.68 (power
/// in y) at z = 5, plane rear side at z = 9, n = 1.5168, VACUUM, detector at z = 106.363 (the
/// paraxial y focus); object at infinity, fields 0 and 1 degree in y, aperture `aperture`.
System cylinder_lens(const std::string& aperture = R"({"type": "epd", "value": 10.0})") {
  return rtt::io::parse_system(
      R"({"schema_version": "0.9.0", "name": "cylinder lens",)"
      R"( "units": {"length": "mm", "wavelength": "um"}, "environment": {"medium": "VACUUM"},)"
      R"( "wavelengths": [{"um": 0.5876, "reference": true}], "aperture": )" +
      aperture +
      R"(, "fields": {"points": [{}, {"y": 1.0}]},)"
      R"( "root": {"type": "assembly", "name": "system", "children": [)"
      R"({"type": "stop", "name": "stop", "surfaces": [{"id": "STO",)"
      R"( "aperture": {"type": "circular", "radius": 10.0}}]},)"
      R"({"type": "lens", "name": "L1", "pose": {"position": [0.0, 0.0, 5.0]},)"
      R"( "material": "CONST:1.5168", "surfaces": [)"
      R"({"id": "L1.S1", "shape": {"base": {"type": "cylinder", "radius": 51.68}},)"
      R"( "aperture": {"type": "rectangular", "half_width_x": 12.7, "half_width_y": 12.7}},)"
      R"({"id": "L1.S2", "pose": {"position": [0.0, 0.0, 4.0]},)"
      R"( "aperture": {"type": "rectangular", "half_width_x": 12.7, "half_width_y": 12.7}}]},)"
      R"({"type": "detector", "name": "image", "pose": {"position": [0.0, 0.0, 106.363]},)"
      R"( "surfaces": [{"id": "IMG"}]}]}, "paths": [{"name": "main", "events": "auto"}]})");
}

rtt::model::SpotRmsOperand spot(std::uint16_t field) {
  rtt::model::SpotRmsOperand s;
  s.path = "main";
  s.field = field;
  s.rings = 3;
  return s;
}

rtt::model::RayOperand ray() {
  rtt::model::RayOperand r;
  r.path = "main";
  r.field = 1;
  r.surface = rtt::model::SurfaceId("IMG");
  r.py = 0.7;
  return r;
}

rtt::model::FirstOrderOperand efl() {
  rtt::model::FirstOrderOperand f;
  f.path = "main";
  f.quantity = rtt::model::FirstOrderQuantity::Efl;
  return f;
}

rtt::model::OpdRmsOperand opd() {
  rtt::model::OpdRmsOperand o;
  o.path = "main";
  o.grid = 9;
  return o;
}

TrialOverlay overlay_of(const System& nominal) {
  return TrialOverlay{{}, rtt::optim::nominal_basis(nominal, library(), nullptr), 1e-3};
}

/// The merit.operand_unsupported diagnostics of constructing the merit function of `s`.
std::vector<rtt::model::Diagnostic> unsupported(const System& s, const TrialOverlay* overlay) {
  try {
    const MeritFunction merit(s, library(), nullptr, overlay);
  } catch (const OptimError& e) {
    return e.diagnostics();
  }
  FAIL("the merit function was constructed");
  return {};
}

/// Spot and ray operands run: finite values, every ray of the spots arrives.
void check_real_ray_operands(const System& s, const TrialOverlay* overlay) {
  const MeritFunction merit(s, library(), nullptr, overlay);
  const MeritEvaluation e = merit.evaluate(merit.start());
  REQUIRE(e.values.size() == 3);
  for (std::size_t i = 0; i < 3; ++i) {
    INFO("operand " << i << ": " << e.undefined[i]);
    CHECK(std::isfinite(e.values[i]));
    CHECK(e.undefined[i].empty());
  }
  for (std::size_t i = 0; i < 2; ++i) {
    INFO("spot operand " << i);
    REQUIRE(e.operand_arrived[i].size() == 37);  // hexapolar, 3 rings
    for (const std::uint8_t a : e.operand_arrived[i]) CHECK(a == 1);
  }
  // Content guard: the line focus in y, the x extent of the bundle stays (afocal in x), so the
  // spot is several mm large and differs between the fields.
  CHECK(e.values[0] > 1.0);
  CHECK(e.values[0] != e.values[1]);
  // The ray operand: py = 0.7 of field 1 lands near the paraxial image height f_y tan(1 deg)
  // in y (f_y = 100 mm, detector in the paraxial y focus).
  CHECK(std::abs(e.values[2] - 100.0 * std::tan(std::numbers::pi / 180.0)) < 0.05);
  REQUIRE(e.generators.size() == 1);
  CHECK(e.valid());
}

}  // namespace

TEST_CASE("OPD of an anamorphic system is a ParaxialError, also with the nominal reference",
          "[optim][anamorphic]") {
  // ADR 0042, point 4: the exit pupil per section is not part of S2; nominal_reference takes
  // the exit pupil from first_order of the nominal system, which is anamorphic itself.
  const rtt::compile::CompiledSystem cs = rtt::compile::compile(cylinder_lens(), library());
  const rtt::compile::PathId path{0};
  REQUIRE(rtt::paraxial::is_anamorphic(cs, path));
  rtt::analysis::OpdOptions options;
  options.grid = 9;
  CHECK_THROWS_AS(rtt::analysis::opd_map(cs, path, 0, 0, options), rtt::paraxial::ParaxialError);
  options.reference = rtt::analysis::nominal_reference(cs, path);
  // The aiming data of the reference exist now (per section); the exit pupil does not.
  CHECK_NOTHROW(options.reference->aim(0));
  CHECK_THROWS_AS(options.reference->exit_pupil(), rtt::paraxial::ParaxialError);
  CHECK_THROWS_AS(rtt::analysis::opd_map(cs, path, 0, 0, options), rtt::paraxial::ParaxialError);
}

TEST_CASE("merit on an anamorphic path: real-ray operands run, paraxial and OPD unsupported",
          "[optim][anamorphic]") {
  for (const bool with_overlay : {false, true}) {
    INFO((with_overlay ? "overlay of a tolerance trial" : "no overlay"));
    System s = cylinder_lens();
    s.optimization.operands = {spot(0), spot(1), ray()};
    rtt::model::SpotGenerator generator;
    generator.path = "main";
    s.optimization.generators = {generator};
    const TrialOverlay overlay = overlay_of(s);
    const TrialOverlay* const o = with_overlay ? &overlay : nullptr;
    check_real_ray_operands(s, o);

    // efl and OPD on the same path, and the wavefront generator: unsupported, with the reason.
    s.optimization.operands = {spot(0), efl(), opd()};
    s.optimization.generators = {generator, rtt::model::WavefrontGenerator{}};
    std::get<rtt::model::WavefrontGenerator>(s.optimization.generators[1]).path = "main";
    const std::vector<rtt::model::Diagnostic> d = unsupported(s, o);
    REQUIRE(d.size() == 3);
    CHECK(d[0].code == "merit.operand_unsupported");
    CHECK(d[0].location == "/optimization/operands/1");
    CHECK_THAT(d[0].message, ContainsSubstring("anamorphic"));
    CHECK(d[1].code == "merit.operand_unsupported");
    CHECK(d[1].location == "/optimization/operands/2");
    CHECK_THAT(d[1].message, ContainsSubstring("anamorphic"));
    CHECK_THAT(d[1].message, ContainsSubstring("exit pupil"));
    CHECK(d[2].code == "merit.operand_unsupported");
    CHECK(d[2].location == "/optimization/generators/1");
    CHECK_THAT(d[2].message, ContainsSubstring("exit pupil"));
  }
}

TEST_CASE("merit on an anamorphic path: the error of the aiming context is the reason",
          "[optim][anamorphic]") {
  // ADR 0042, point 4: an anamorphic path counts as aimable if the aiming context per section
  // builds; an image-space F-number is ambiguous there (point 3, table).
  for (const bool with_overlay : {false, true}) {
    INFO((with_overlay ? "overlay of a tolerance trial" : "no overlay"));
    System s = cylinder_lens(R"({"type": "image_fnumber", "value": 10.0})");
    s.optimization.operands = {spot(0)};
    const TrialOverlay overlay = overlay_of(s);
    const std::vector<rtt::model::Diagnostic> d = unsupported(s, with_overlay ? &overlay : nullptr);
    REQUIRE(d.size() == 1);
    CHECK(d[0].code == "merit.operand_unsupported");
    CHECK(d[0].location == "/optimization/operands/0");
    CHECK_THAT(d[0].message, ContainsSubstring("the path is not supported"));
    CHECK_THAT(d[0].message, ContainsSubstring("F-number"));
  }
}

TEST_CASE("merit on an anamorphic path without a stop: NoStopError with and without overlay",
          "[optim][anamorphic]") {
  // ADR 0030, addendum #315: a path without a stop is not a reason of the path check; the
  // NoStopError of the aiming is thrown on, with the overlay as without it.
  System s = cylinder_lens();
  s.root.children.erase(s.root.children.begin());  // the stop element
  s.optimization.operands = {spot(0)};
  const TrialOverlay overlay = overlay_of(s);
  CHECK_THROWS_AS(MeritFunction(s, library(), nullptr, nullptr), rtt::compile::NoStopError);
  CHECK_THROWS_AS(MeritFunction(s, library(), nullptr, &overlay), rtt::compile::NoStopError);
}
