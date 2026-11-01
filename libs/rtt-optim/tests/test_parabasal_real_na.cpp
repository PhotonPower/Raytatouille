// The F-number from parabasal rays with object_na read real (#296, ADR 0040, point 5): the trial's
// image_fnumber divides the measured stop hits by the paraxial stop radius, so in the nominal
// system it agrees with rtt::paraxial (u = NA / n) within the Richardson bound of ADR 0030,
// addendum #256, although the rays aim at the real cone.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <string>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/overlay.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/sources.hpp"

using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::FieldType;
using rtt::model::FirstOrderQuantity;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::optim::MeritFunction;
using rtt::optim::TrialOverlay;

namespace {

constexpr double kEps = std::numeric_limits<double>::epsilon();
constexpr double kDelta = 1e-3;  // TrialOverlay::parabasal_step (ADR 0036, point 2)

const rtt::material::MaterialLibrary& library() {
  static const rtt::material::MaterialLibrary materials;
  return materials;
}

Surface plane(const std::string& id, double z_mm = 0.0, double radius = 0.0) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  if (radius != 0.0) s.shape.base = rtt::model::Conic{Param(radius), Param(0.0)};
  return s;
}

/// The two lenses of test_parabasal.cpp (CONST:1.6, CONST:1.5, stop r = 8 between them, finite
/// object at 200 mm, vacuum) with object-space NA 0.04 read `reading`.
System doublet(rtt::model::ApertureReading reading) {
  System s;
  s.name = "two lenses, finite object";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  s.object = {false, Param(200.0)};
  s.aperture = {rtt::model::SystemApertureType::ObjectSpaceNA, Param(0.04)};
  s.aperture.reading = reading;
  s.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}}};
  Surface stop = plane("STO");
  stop.aperture = rtt::model::CircularAperture{8.0, 0.0};
  s.root.children = {
      {Element{"L1",
               ElementKind::Lens,
               Pose::along_z(0.0),
               "CONST:1.6",
               {plane("L1.S1", 0.0, 50.0), plane("L1.S2", 5.0, -200.0)}}},
      {Element{"stop", ElementKind::Stop, Pose::along_z(10.0), std::nullopt, {stop}}},
      {Element{"L2",
               ElementKind::Lens,
               Pose::along_z(15.0),
               "CONST:1.5",
               {plane("L2.S1", 0.0, 80.0), plane("L2.S2", 4.0, -80.0)}}},
      {Element{"D", ElementKind::Detector, Pose::along_z(120.0), std::nullopt, {plane("IMG")}}}};
  return s;
}

/// The immersion front of test_real_na.cpp (object in oil n = 1.52 at z = -0.1, index-matched
/// lens ending in a sphere of R = 1 mm around the object point, plane stop r = 10 in vacuum at
/// z = 2) with object-space NA 1.45 read `reading`, and a positive lens (CONST:1.6, R = +-15 mm)
/// at z = 20 that forms a real image.
System immersion_objective(rtt::model::ApertureReading reading) {
  System s;
  s.name = "immersion front and focusing lens";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  s.object = {false, Param(0.1)};
  s.object.medium = "CONST:1.52";
  s.aperture = {rtt::model::SystemApertureType::ObjectSpaceNA, Param(1.45)};
  s.aperture.reading = reading;
  s.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}}};
  Surface stop = plane("STO");
  stop.aperture = rtt::model::CircularAperture{10.0, 0.0};
  s.root.children = {
      {Element{"F",
               ElementKind::Lens,
               Pose::along_z(0.0),
               "CONST:1.52",
               {plane("F.S1", 0.0), plane("F.S2", 0.9, -1.0)}}},
      {Element{"stop", ElementKind::Stop, Pose::along_z(2.0), std::nullopt, {stop}}},
      {Element{"L",
               ElementKind::Lens,
               Pose::along_z(20.0),
               "CONST:1.6",
               {plane("L.S1", 0.0, 15.0), plane("L.S2", 4.0, -15.0)}}},
      {Element{"D", ElementKind::Detector, Pose::along_z(60.0), std::nullopt, {plane("IMG")}}}};
  return s;
}

double f_number(const System& system, const TrialOverlay* overlay = nullptr) {
  System s = system;
  rtt::model::FirstOrderOperand op;
  op.quantity = FirstOrderQuantity::ImageFNumber;
  op.path = "main";
  s.optimization.operands = {op};
  const MeritFunction merit(s, library(), nullptr, overlay);
  const rtt::optim::MeritEvaluation e = merit.evaluate(merit.start());
  INFO(e.undefined.at(0));
  return e.values.at(0);
}

TrialOverlay overlay_of(const System& nominal, double step) {
  return TrialOverlay{{}, rtt::optim::nominal_basis(nominal, library(), nullptr), step};
}

}  // namespace

/// The Richardson check of test_parabasal.cpp for image_fnumber of `s` (see the test cases).
void check_richardson(const System& s) {
  const double x0 = f_number(s);
  const TrialOverlay one = overlay_of(s, kDelta);
  const TrialOverlay two = overlay_of(s, 2.0 * kDelta);
  const double x1 = f_number(s, &one);
  const double x2 = f_number(s, &two);
  REQUIRE(std::isfinite(x0));
  REQUIRE(std::isfinite(x1));
  REQUIRE(std::isfinite(x2));
  const rtt::compile::CompiledSystem cs = rtt::compile::compile(s, library());
  const double n_image = rtt::paraxial::first_order(cs, {0}, 0).image_index;
  const double u_m = 1.0 / (2.0 * n_image * x0);
  const double r = 20.0 * kEps / (kDelta * u_m) * std::abs(x0);
  INFO("X0 " << x0 << ", X(delta) - X0 " << x1 - x0 << ", X(2 delta) - X(delta) " << x2 - x1);
  CHECK(std::abs(x1 - x0) <= 2.0 / 3.0 * std::abs(x2 - x1) + r);
}

TEST_CASE("parabasal: image_fnumber with object_na read real", "[optim][parabasal][real-na]") {
  // Bound as in test_parabasal.cpp ("nominal systems within the Richardson bound"):
  // |X(delta) - X_0| <= (2/3) |X(2 delta) - X(delta)| + 20 eps / (delta u_m) |X_0|, u_m the image
  // marginal slope 1 / (2 n' F#). X_0 is rtt::paraxial (no overlay), the same for both readings.
  // The parabasal offsets are delta |R_s,real| on the stop with reading real (ADR 0040, point 5),
  // larger than with paraxial. A division by the real R_s instead would be off by
  // tan U / sin U - 1 = 8e-4 at NA 0.04, far above the bound (checked in the red run).
  for (const auto reading :
       {rtt::model::ApertureReading::Paraxial, rtt::model::ApertureReading::Real}) {
    INFO("reading " << static_cast<int>(reading));
    check_richardson(doublet(reading));
  }
  // Content guard: the readings differ in the aim target (R_s), not in rtt::paraxial.
  const rtt::compile::CompiledSystem real =
      rtt::compile::compile(doublet(rtt::model::ApertureReading::Real), library());
  CHECK(rtt::trace::stop_radius(real, {0}, 0) != rtt::trace::paraxial_stop_radius(real, {0}, 0));
  CHECK(f_number(doublet(rtt::model::ApertureReading::Real)) ==
        f_number(doublet(rtt::model::ApertureReading::Paraxial)));
}

TEST_CASE("parabasal: image_fnumber at NA 1.45 read real", "[optim][parabasal][real-na]") {
  // ADR 0040, point 5, side effect: with reading real the parabasal offsets on the stop grow to
  // delta |R_s,real| = delta 2.1 tan U, 1 / cos U = 3.3 times delta 2.1 NA / n at NA 1.45 in
  // n = 1.52, so the O(delta^2) error grows about 11-fold. The Richardson bound of ADR 0036, point
  // 2, still holds for both readings (same check as above).
  for (const auto reading :
       {rtt::model::ApertureReading::Paraxial, rtt::model::ApertureReading::Real}) {
    INFO("reading " << static_cast<int>(reading));
    check_richardson(immersion_objective(reading));
  }
  // Content guard: the real stop target is 1 / cos U times the paraxial one (both 2.1 mm from the
  // object, slopes tan U and NA / n).
  const rtt::compile::CompiledSystem real =
      rtt::compile::compile(immersion_objective(rtt::model::ApertureReading::Real), library());
  const double ratio =
      rtt::trace::stop_radius(real, {0}, 0) / rtt::trace::paraxial_stop_radius(real, {0}, 0);
  CHECK(std::abs(ratio - 1.0 / std::cos(std::asin(1.45 / 1.52))) <= 1e-9);
}
