// The limit of the affine pupil mapping for the FFT analyses (ADR 0033, addendum #286): above
// PupilFunctionOptions::mapping_limit psf, through_focus, mtf, encircled_energy and psf_at reject
// with AnalysisError, code pupil.mapping_too_nonlinear at /paths/i; pupil_function itself does
// not reject; mapping_limit = 1 switches the check off.
//
// Ideal lens (ADR 0031): tan theta = h / f. With the entrance pupil radius a = f tan theta_max
// and t = tan theta_max, a ray at pupil radius rho leaves at tan theta = rho t, so its direction
// sine is s(rho) = rho t / sqrt(1 + rho^2 t^2) (vacuum, n' = 1), and the mapping at p = 0 is
// |M| = t. With delta = |s - s0 - M p| / ||M||_2 (ADR 0032, point 4):
//   delta(rho) = rho (1 - 1 / sqrt(1 + rho^2 t^2)),
// increasing in rho, 1 - cos theta_max at rho = 1: 0.13397 for NA = sin theta_max = 0.5.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/analysis/psf.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"

using rtt::analysis::AnalysisError;
using rtt::analysis::PupilFunction;
using rtt::analysis::PupilFunctionOptions;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;

namespace {

constexpr double kF = 100.0;  // mm
constexpr int kGrid = 16;

Surface plane(const std::string& id) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(0.0);
  return s;
}

/// Ideal lens f = 100 in vacuum at z = 0, the stop just before it (z = -5), the image surface in
/// the rear focal plane; object at infinity, one field on the axis, NA' = sin theta_max with
/// EPD = 2 f tan theta_max.
CompiledSystem ideal_lens(double tan_max) {
  rtt::model::System s;
  s.name = "ideal lens";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(2.0 * kF * tan_max)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  Surface stop = plane("STO");
  stop.aperture = rtt::model::CircularAperture{1.1 * kF * tan_max, 0.0};
  Surface lens = plane("IL");
  lens.aperture = rtt::model::CircularAperture{1.2 * kF * tan_max, 0.0};
  lens.interaction = rtt::model::IdealLens{Param(kF), std::nullopt};
  s.root.children = {
      {Element{"stop", ElementKind::Stop, Pose::along_z(-5.0), std::nullopt, {stop}}},
      {Element{"lens", ElementKind::ThinElement, Pose::along_z(0.0), std::nullopt, {lens}}},
      {Element{"D", ElementKind::Detector, Pose::along_z(kF), std::nullopt, {plane("IMG")}}}};
  s.paths = {{"main", true, {}}};
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

/// The paraboloid at F/1.7 (#202): delta_max about 0.0225 (ADR 0032, point 4).
CompiledSystem paraboloid() {
  return rtt::compile::compile(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m2/paraboloid_stop.rtt.json"),
      rtt::material::MaterialLibrary{});
}

PupilFunctionOptions options(double limit) {
  PupilFunctionOptions o;
  o.grid = kGrid;
  o.mapping_limit = limit;
  return o;
}

double delta_max(const PupilFunction& p) {
  REQUIRE(p.mapping.has_value());
  return p.mapping->delta_max;  // NOLINT(bugprone-unchecked-optional-access): REQUIRE above
}

bool has(const std::vector<rtt::model::Diagnostic>& w, const std::string& code) {
  return std::any_of(w.begin(), w.end(), [&](const auto& d) { return d.code == code; });
}

/// `run` throws AnalysisError with the code and location of the rejection, naming delta_max.
void require_rejected(const std::function<void()>& run, const std::string& what) {
  INFO(what);
  bool thrown = false;
  try {
    run();
  } catch (const AnalysisError& e) {
    thrown = true;
    CHECK(e.code() == std::optional<std::string>("pupil.mapping_too_nonlinear"));
    CHECK(e.location() == std::optional<std::string>("/paths/0"));
    CHECK(std::string(e.what()).find("delta_max") != std::string::npos);
  }
  CHECK(thrown);
}

/// All five FFT analyses of `cs` with `o`, and those on the pupil function.
void require_all_rejected(const CompiledSystem& cs, const PupilFunctionOptions& o) {
  const PathId path{0};
  rtt::analysis::PsfOptions po;
  po.pupil = o;
  require_rejected([&] { static_cast<void>(rtt::analysis::psf(cs, path, 0, po)); }, "psf");
  const std::vector<double> shifts{-0.01, 0.0, 0.01};
  require_rejected(
      [&] { static_cast<void>(rtt::analysis::through_focus(cs, path, 0, 0, shifts, o)); },
      "through_focus");
  rtt::analysis::MtfOptions mo;
  mo.pupil = o;
  require_rejected([&] { static_cast<void>(rtt::analysis::mtf(cs, path, 0, mo)); }, "mtf");
  rtt::analysis::EncircledEnergyOptions eo;
  eo.pupil = o;
  const std::vector<double> radii{1e-3};
  require_rejected(
      [&] { static_cast<void>(rtt::analysis::encircled_energy(cs, path, 0, radii, eo)); },
      "encircled_energy");
  // pupil_function itself does not reject; the analyses of a pupil function do.
  const PupilFunction p = rtt::analysis::pupil_function(cs, path, 0, 0, o);
  CHECK(p.mapping_limit == o.mapping_limit);
  const std::vector<rtt::analysis::Point2> x{{0.0, 0.0}};
  require_rejected([&] { static_cast<void>(rtt::analysis::psf_at(p, x)); }, "psf_at");
  require_rejected([&] { static_cast<void>(rtt::analysis::mtf(p)); }, "mtf of a pupil function");
  require_rejected([&] { static_cast<void>(rtt::analysis::encircled_energy(p, radii)); },
                   "encircled_energy of a pupil function");
  // refocus() keeps the limit.
  CHECK(rtt::analysis::refocus(cs, p, 0.01).mapping_limit == o.mapping_limit);
}

}  // namespace

TEST_CASE("mapping limit: the ideal lens at NA 0.5 is rejected", "[analysis][psf][mapping]") {
  const double t = std::tan(std::asin(0.5));
  const CompiledSystem cs = ideal_lens(t);
  const PupilFunctionOptions o = options(PupilFunctionOptions{}.mapping_limit);
  CHECK(o.mapping_limit == 0.05);  // the default (maintainer decision on #286)
  const PupilFunction p = rtt::analysis::pupil_function(cs, PathId{0}, 0, 0, o);
  // delta_max at the largest |p| of the arrived cells; tolerance 1e-5 relative: M comes from
  // central differences at +-1e-3 (relative error t^2 1e-6 / 2), the rest is rounding.
  double rho = 0.0;
  for (std::size_t j = 0; j < p.fraction.size(); ++j) {
    if (p.fraction[j] > 0.0 && p.status[j] == rtt::trace::RayStatus::Alive) {
      rho = std::max(rho, std::hypot(p.px[j], p.py[j]));
    }
  }
  REQUIRE(rho > 0.9);  // content guard: the rim cells arrived
  const double expected = rho * (1.0 - 1.0 / std::sqrt(1.0 + rho * rho * t * t));
  INFO("rho_max " << rho << ", expected " << expected << ", delta_max " << delta_max(p));
  CHECK(std::abs(delta_max(p) - expected) <= 1e-5 * expected);
  CHECK(delta_max(p) < 1.0 - std::cos(std::asin(0.5)));  // 0.13397 at rho = 1
  CHECK(delta_max(p) > 0.05);
  CHECK(has(p.warnings, "pupil.mapping_nonlinear"));
  require_all_rejected(cs, o);
}

TEST_CASE("mapping limit: 1 switches the check off", "[analysis][psf][mapping]") {
  // ADR 0033, addendum #286: mapping_limit = 1 never rejects, whatever delta_max; the result is
  // that of the analyses before the limit, with the warning.
  const CompiledSystem cs = ideal_lens(std::tan(std::asin(0.5)));
  rtt::analysis::PsfOptions po;
  po.pupil = options(1.0);
  const rtt::analysis::Psf psf = rtt::analysis::psf(cs, PathId{0}, 0, po);
  CHECK(has(psf.warnings, "pupil.mapping_nonlinear"));
  CHECK(psf.strehl > 0.99);  // W = 0 for the ideal lens
  // The limit changes no number: the paraboloid (delta_max < 0.05) gives bitwise the same PSF
  // with 1 and with the default.
  const CompiledSystem para = paraboloid();
  rtt::analysis::PsfOptions a;
  a.pupil = options(1.0);
  rtt::analysis::PsfOptions b;
  b.pupil = options(0.05);
  const rtt::analysis::Psf pa = rtt::analysis::psf(para, PathId{0}, 0, a);
  const rtt::analysis::Psf pb = rtt::analysis::psf(para, PathId{0}, 0, b);
  CHECK(pa.intensity == pb.intensity);
  CHECK(pa.strehl == pb.strehl);
}

TEST_CASE("mapping limit: the paraboloid at F/1.7 warns, and the bound is inclusive",
          "[analysis][psf][mapping]") {
  // delta_max about 0.0225: above mapping_warning 1e-2, below mapping_limit 0.05.
  const CompiledSystem cs = paraboloid();
  const PupilFunctionOptions o = options(0.05);
  rtt::analysis::PsfOptions po;
  po.pupil = o;
  const rtt::analysis::Psf psf = rtt::analysis::psf(cs, PathId{0}, 0, po);
  CHECK(has(psf.warnings, "pupil.mapping_nonlinear"));
  const PupilFunction p = rtt::analysis::pupil_function(cs, PathId{0}, 0, 0, o);
  const double d = delta_max(p);
  CHECK(d > 0.02);
  CHECK(d < 0.03);
  // Content guard: 0.02 rejects. The bound itself is allowed (delta_max > limit rejects), the
  // next double below it rejects.
  require_all_rejected(cs, options(0.02));
  po.pupil = options(d);
  CHECK_NOTHROW(rtt::analysis::psf(cs, PathId{0}, 0, po));
  po.pupil = options(std::nextafter(d, 0.0));
  require_rejected([&] { static_cast<void>(rtt::analysis::psf(cs, PathId{0}, 0, po)); },
                   "psf just below delta_max");
}

TEST_CASE("mapping limit: values outside [0, 1] are invalid", "[analysis][psf][mapping]") {
  const CompiledSystem cs = paraboloid();
  for (const double bad : {-0.1, 1.5, std::numeric_limits<double>::quiet_NaN()}) {
    INFO(bad);
    CHECK_THROWS_AS(rtt::analysis::pupil_function(cs, PathId{0}, 0, 0, options(bad)),
                    std::invalid_argument);
    rtt::analysis::PsfOptions po;
    po.pupil = options(bad);
    CHECK_THROWS_AS(rtt::analysis::psf(cs, PathId{0}, 0, po), std::invalid_argument);
  }
  // Below mapping_warning is allowed (no ordering rule, #286 F3): 0 rejects every nonlinearity.
  PupilFunctionOptions o = options(0.0);
  o.mapping_warning = 0.5;
  CHECK_NOTHROW(rtt::analysis::pupil_function(cs, PathId{0}, 0, 0, o));
  require_all_rejected(cs, o);
}
