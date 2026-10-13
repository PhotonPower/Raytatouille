// Paraxial operands of a trial from parabasal rays (ADR 0030, addendum #256; ADR 0036, point 2):
// efl, bfl, image_fnumber and magnification from five rays around the chief ray, with the nominal
// basis of the trial's overlay. On the nominal system they agree with rtt::paraxial within the
// Richardson bound derived in the addendum; on a tilted or decentered ideal lens they meet hand
// values of the collinear mapping (ADR 0031).

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/perturbation.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/overlay.hpp"
#include "rtt/paraxial/paraxial.hpp"

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

Surface plane(const std::string& id, double z_mm = 0.0, std::optional<double> radius = {}) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  if (radius) s.shape.base = rtt::model::Conic{Param(*radius), Param(0.0)};
  return s;
}

Element stop_at(double z, double radius) {
  Surface s = plane("STO");
  s.aperture = rtt::model::CircularAperture{radius, 0.0};
  return Element{"stop", ElementKind::Stop, Pose::along_z(z), std::nullopt, {s}};
}

Element detector_at(double z) {
  return Element{"D", ElementKind::Detector, Pose::along_z(z), std::nullopt, {plane("IMG")}};
}

System base(const std::string& name) {
  System s;
  s.name = name;
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  return s;
}

/// Ideal lens f in vacuum (ADR 0031): stop r = 10 at z = -5, lens at z = 0 (aperture 30), detector
/// in the rear focal plane; object at infinity, EPD 10 (r_EP = 5), field 0.
System ideal_lens(double f) {
  System s = base("ideal lens");
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  Surface lens = plane("IL");
  lens.aperture = rtt::model::CircularAperture{30.0, 0.0};
  lens.interaction = rtt::model::IdealLens{Param(f), std::nullopt};
  s.root.children = {
      {stop_at(-5.0, 10.0)},
      {Element{"lens", ElementKind::ThinElement, Pose::along_z(0.0), std::nullopt, {lens}}},
      {detector_at(f)}};
  return s;
}

/// Biconvex singlet CONST:1.5, R = +-60, 6 mm, stop r = 12.5 at z = -2; object at infinity, EPD
/// 20 (F/3): enough spherical aberration that the O(delta^2) term is visible.
System singlet(double r1 = 60.0, double r2 = -60.0) {
  System s = base("singlet");
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(20.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.children = {{stop_at(-2.0, 12.5)},
                     {Element{"L",
                              ElementKind::Lens,
                              Pose::along_z(0.0),
                              "CONST:1.5",
                              {plane("L.S1", 0.0, r1), plane("L.S2", 6.0, r2)}}},
                     {detector_at(60.0)}};
  return s;
}

/// Two lenses (CONST:1.6, CONST:1.5) with the stop between them, finite object at 200 mm,
/// object height fields, object-space NA 0.04.
System doublet_finite() {
  System s = base("two lenses, finite object");
  s.object = {false, Param(200.0)};
  s.aperture = {rtt::model::SystemApertureType::ObjectSpaceNA, Param(0.04)};
  s.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}, {0.0, 3.0, 1.0}}};
  s.root.children = {{Element{"L1",
                              ElementKind::Lens,
                              Pose::along_z(0.0),
                              "CONST:1.6",
                              {plane("L1.S1", 0.0, 50.0), plane("L1.S2", 5.0, -200.0)}}},
                     {stop_at(10.0, 8.0)},
                     {Element{"L2",
                              ElementKind::Lens,
                              Pose::along_z(15.0),
                              "CONST:1.5",
                              {plane("L2.S1", 0.0, 80.0), plane("L2.S2", 4.0, -80.0)}}},
                     {detector_at(120.0)}};
  return s;
}

/// The concave mirror of tests/reference/m2/paraboloid_stop.rtt.json (R = -200, EPD 60, light
/// leaves towards -z) with the conic `k`: k = -1 the paraboloid (no spherical aberration on
/// axis), k = 0 the sphere. The reflection reverses the image space (rtt::paraxial: image
/// direction -1), the folded case of ADR 0036, point 2; a 45 degree fold is not rotationally
/// symmetric and has no paraxial model.
System mirror(double k) {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m2/paraboloid_stop.rtt.json");
  auto& m = std::get<Element>(s.root.children.at(1).value);
  m.surfaces.at(0).shape.base = rtt::model::Conic{Param(-200.0), Param(k)};
  return s;
}

/// The value of one first-order operand on `path` "main" (or the system's first path), with the
/// overlay if given.
double value(const System& system,
             FirstOrderQuantity quantity,
             const TrialOverlay* overlay = nullptr) {
  System s = system;
  rtt::model::FirstOrderOperand op;
  op.quantity = quantity;
  op.path = s.paths.at(0).name;
  s.optimization.operands = {op};
  s.optimization.generators.clear();
  const MeritFunction merit(s, library(), nullptr, overlay);
  const rtt::optim::MeritEvaluation e = merit.evaluate(merit.start());
  INFO(e.undefined.at(0));
  return e.values.at(0);
}

TrialOverlay overlay_of(const System& nominal,
                        rtt::compile::Perturbation perturbation = {},
                        double step = kDelta) {
  return TrialOverlay{std::move(perturbation),
                      rtt::optim::nominal_basis(nominal, library(), nullptr), step};
}

/// The variant of the addendum for the axis-parallel bundle: object at infinity, field 0 degree,
/// system aperture StopSize.
System infinite(const System& s) {
  System v = s;
  v.object.at_infinity = true;
  v.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  v.aperture = {rtt::model::SystemApertureType::StopSize, Param(0.0)};
  return v;
}

/// Largest |coordinate| of the vertices of the system's surfaces, mm.
double size_of(const rtt::compile::CompiledSystem& cs) {
  double l = 1.0;
  for (const auto& surface : cs.surfaces()) {
    l = std::max(l, surface.to_global.translation().cwiseAbs().maxCoeff());
  }
  return l;
}

}  // namespace

TEST_CASE("parabasal: nominal systems within the Richardson bound", "[optim][parabasal]") {
  // Addendum #256: X(delta) = X_0 + c delta^2 + O(delta^4) for central differences, so
  // X(2 delta) - X(delta) = 3 c delta^2 + O(delta^4) and |X(delta) - X_0| <= (2/3) |X(2 delta) -
  // X(delta)| + r (factor 2 against the O(delta^4) rest). Rounding r: the differences of the slopes
  // are delta |u_m| large (u_m the marginal slope per pupil unit) with a few eps of rounding, r =
  // 20 eps / (delta |u_m|) |X_0|; for bfl the coordinates of size L against the offset delta r_EP,
  // r = 20 eps L / (delta r_EP) |X_0|. X_0 is the value of rtt::paraxial (the merit function
  // without overlay); image_fnumber is the paraxial working F-number 1/(2 n' |u'|) of the axial
  // bundle (prescription().paraxial_working_f_number).
  struct Case {
    std::string name;
    System system;
    std::vector<FirstOrderQuantity> quantities;
    bool aberrated;  // content guard: X(2 delta) != X(delta) beyond the rounding
  };
  using Q = FirstOrderQuantity;
  const std::vector<Case> cases = {
      {"singlet", singlet(), {Q::Efl, Q::Bfl, Q::ImageFNumber}, true},
      {"diverging singlet", singlet(-60.0, 60.0), {Q::Efl, Q::Bfl}, true},
      {"two lenses, finite object",
       doublet_finite(),
       {Q::Efl, Q::Bfl, Q::ImageFNumber, Q::Magnification},
       true},
      {"spherical mirror", mirror(0.0), {Q::Efl, Q::Bfl, Q::ImageFNumber}, true},
      {"paraboloid mirror", mirror(-1.0), {Q::Efl, Q::Bfl, Q::ImageFNumber}, false},
      {"ideal lens", ideal_lens(100.0), {Q::Efl, Q::Bfl, Q::ImageFNumber}, false},
      {"object-space telecentric singlet",
       rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/telecentric_singlet.rtt.json"),
       {Q::ImageFNumber, Q::Magnification},
       false},
  };
  for (const Case& c : cases) {
    INFO(c.name);
    const TrialOverlay one = overlay_of(c.system);
    const TrialOverlay two = overlay_of(c.system, {}, 2.0 * kDelta);
    // u_m and r_EP of the bundles: the axial bundle of the system (image side, from the working
    // F-number) and the axis-parallel bundle of the variant (entrance pupil of StopSize).
    const rtt::compile::CompiledSystem cs = rtt::compile::compile(c.system, library());
    const double n_image = rtt::paraxial::first_order(cs, {0}, 0).image_index;
    for (const Q q : c.quantities) {
      INFO("quantity " << static_cast<int>(q));
      const double x0 = value(c.system, q);
      const double x1 = value(c.system, q, &one);
      const double x2 = value(c.system, q, &two);
      REQUIRE(std::isfinite(x0));
      REQUIRE(std::isfinite(x1));
      REQUIRE(std::isfinite(x2));
      double r = 0.0;
      if (q == Q::Efl || q == Q::Bfl) {
        const rtt::compile::CompiledSystem ci =
            rtt::compile::compile(infinite(c.system), library());
        const rtt::paraxial::FirstOrder fo = rtt::paraxial::first_order(ci, {0}, 0);
        REQUIRE(fo.entrance_pupil);
        REQUIRE(fo.entrance_pupil->diameter);
        REQUIRE(fo.efl);
        const double r_ep = *fo.entrance_pupil->diameter / 2.0;
        r = q == Q::Efl ? 20.0 * kEps / (kDelta * r_ep / std::abs(*fo.efl)) * std::abs(x0)
                        : 20.0 * kEps * size_of(ci) / (kDelta * r_ep) * std::abs(x0);
      } else {
        const double f_number = value(c.system, Q::ImageFNumber);
        const double u_m = 1.0 / (2.0 * n_image * f_number);
        r = 20.0 * kEps / (kDelta * u_m) * std::abs(x0);
      }
      INFO("X0 " << x0 << ", X(delta) - X0 " << x1 - x0 << ", X(2 delta) - X(delta) " << x2 - x1
                 << ", r " << r);
      CHECK(std::abs(x1 - x0) <= 2.0 / 3.0 * std::abs(x2 - x1) + r);
      // Sign: the same as rtt::paraxial (efl < 0 for the diverging lens; the mirrors).
      CHECK(std::signbit(x1) == std::signbit(x0));
      if (c.aberrated && q == Q::Efl) CHECK(std::abs(x2 - x1) > r);
    }
  }
}

TEST_CASE("parabasal: a tilted ideal lens, hand values", "[optim][parabasal]") {
  // The ideal lens maps collinearly (ADR 0031): an axis-parallel bundle meets in one point F at
  // the distance L = f / n_z from the lens centre along the undeviated chief ray, n the unit lens
  // axis after the tilt. A parallel ray at the transverse offset o (in x, y) crosses the lens plane
  // n . p = 0 at the axial position -g . o with g = (n_x, n_y) / n_z, so its slope against the
  // chief ray in image space is o / (L + g . o). The pair at +-y e_j therefore gives
  // ((y e_j)/(L + y g_j) + (y e_j)/(L - y g_j)) / (2 y) = e_j L / (L^2 - y^2 g_j^2) =: e_j k_j, and
  // the mean principal value of diag(k_x, k_y) is (k_x + k_y)/2: efl = 2 / (k_x + k_y), the limit
  // f / n_z for y -> 0, with y = delta r_EP = 5e-3 mm. bfl: every ray of the bundle runs through F,
  // so h + s w = 0 holds exactly at s = L for every ray, bfl = L.
  // The tilts about x and about x and y (an oblique axis: the map is not diagonal in the axes of
  // the start offsets once g has two components).
  const double f = 100.0;
  const double y = kDelta * 5.0;
  const double r = 20.0 * kEps / (kDelta * 5.0 / f) * f;  // rounding as in the Richardson case
  const System nominal = ideal_lens(f);
  for (const auto& [a_deg, b_deg] :
       std::vector<std::pair<double, double>>{{10.0, 0.0}, {7.0, 7.0}}) {
    INFO("tilt x " << a_deg << ", y " << b_deg << " degree");
    const double a = a_deg * std::numbers::pi / 180.0;
    const double b = b_deg * std::numbers::pi / 180.0;
    rtt::compile::Perturbation p;
    p.nodes["lens"].rotation_rad = rtt::math::Vec3(a, b, 0.0);
    const TrialOverlay o = overlay_of(nominal, p);
    // n = Rx(a) Ry(b) e_z = (sin b, -sin a cos b, cos a cos b) (R = Rx Ry Rz, perturbation.hpp).
    const double nx = std::sin(b);
    const double ny = -std::sin(a) * std::cos(b);
    const double nz = std::cos(a) * std::cos(b);
    const double l = f / nz;
    const double kx = l / (l * l - y * y * (nx / nz) * (nx / nz));
    const double ky = l / (l * l - y * y * (ny / nz) * (ny / nz));
    const double efl = value(nominal, FirstOrderQuantity::Efl, &o);
    const double bfl = value(nominal, FirstOrderQuantity::Bfl, &o);
    INFO("efl " << efl << " against " << 2.0 / (kx + ky) << ", bfl " << bfl << " against " << l);
    CHECK(std::abs(efl - 2.0 / (kx + ky)) <= r);
    CHECK(std::abs(bfl - l) <= r);
    // Content guard: the tilt is seen (f / n_z differs from f by far more than r).
    CHECK(std::abs(efl - f) > 1e3 * r);
  }
}

TEST_CASE("parabasal: a decentered ideal lens, hand values", "[optim][parabasal]") {
  // Lens centre at (d, 0, 0), d = 0.1 mm; the chief ray (stop centre, axis-parallel) crosses the
  // lens at the origin and runs to the focus F = (d, 0, f), L = |F| = sqrt(f^2 + d^2). Slopes
  // against the chief ray (tan of the angle, |c x r| / (c . r) with c = (d, 0, f)):
  // - x pair, ray from (y, 0, 0) to F, r = (d - y, 0, f): tan = f y / (L^2 - d y); the pair gives
  //   k_x = f L^2 / (L^4 - d^2 y^2);
  // - y pair, ray from (0, y, 0), r = (d, -y, f): |c x r| = y L, c . r = L^2, tan = y / L, k_y =
  // 1/L. efl = 2 / (k_x + k_y) (about f + d^2/(2f) + ...: the decenter tilts the image-side chief
  // ray); bfl = L (all rays through F).
  const double f = 100.0;
  const double d = 0.1;
  const double y = kDelta * 5.0;
  const double r = 20.0 * kEps / (kDelta * 5.0 / f) * f;
  const System nominal = ideal_lens(f);
  rtt::compile::Perturbation p;
  p.nodes["lens"].translation_mm = rtt::math::Vec3(d, 0.0, 0.0);
  const TrialOverlay o = overlay_of(nominal, p);
  const double l = std::hypot(f, d);
  const double kx = f * l * l / (l * l * l * l - d * d * y * y);
  const double ky = 1.0 / l;
  const double efl = value(nominal, FirstOrderQuantity::Efl, &o);
  const double bfl = value(nominal, FirstOrderQuantity::Bfl, &o);
  INFO("efl " << efl << " against " << 2.0 / (kx + ky) << ", bfl " << bfl << " against " << l);
  CHECK(std::abs(efl - 2.0 / (kx + ky)) <= r);
  CHECK(std::abs(bfl - l) <= r);
  // Without the overlay the decentered system is no paraxial system at all.
  System moved = nominal;
  std::get<Element>(moved.root.children[1].value).pose.position[0] = Param(d);
  CHECK_THROWS_AS(value(moved, FirstOrderQuantity::Efl), rtt::optim::OptimError);
}

TEST_CASE("parabasal: object-space telecentric, efl and bfl are not supported in a trial",
          "[optim][parabasal]") {
  // Addendum #256: with the entrance pupil at infinity an axis-parallel bundle cannot be aimed at
  // the stop (every such bundle converges on the stop centre); efl and bfl of a trial are a start
  // error with the reason, never a silent value. image_fnumber and magnification (the bundle of the
  // axis point, aimed as object-space telecentric, #96) are defined.
  const System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/telecentric_singlet.rtt.json");
  const TrialOverlay o = overlay_of(s);
  for (const FirstOrderQuantity q : {FirstOrderQuantity::Efl, FirstOrderQuantity::Bfl}) {
    try {
      static_cast<void>(value(s, q, &o));
      FAIL("no start error");
    } catch (const rtt::optim::OptimError& e) {
      REQUIRE(e.diagnostics().size() == 1);
      CHECK(e.diagnostics()[0].code == "merit.operand_unsupported");
      CHECK(e.diagnostics()[0].location == "/optimization/operands/0");
      CHECK(e.diagnostics()[0].message.find("axis-parallel") != std::string::npos);
    }
    // Without overlay the nominal values exist (content guard: the gap is the trial's only).
    CHECK(std::isfinite(value(s, q)));
  }
  // Their values: the Richardson case "object-space telecentric singlet" above.
  for (const FirstOrderQuantity q :
       {FirstOrderQuantity::ImageFNumber, FirstOrderQuantity::Magnification}) {
    CHECK(std::isfinite(value(s, q, &o)));
  }
}

TEST_CASE("parabasal: magnification with the object at infinity stays a start error",
          "[optim][parabasal]") {
  const System s = singlet();
  const TrialOverlay o = overlay_of(s);
  CHECK_THROWS_AS(value(s, FirstOrderQuantity::Magnification, &o), rtt::optim::OptimError);
}
