// Sensitivity, RSS and inverse sensitivity (ADR 0036, point 6; #258). The systems are built
// here with ideal lenses (ADR 0031), plates and plane mirrors in VACUUM; every expected value is
// derived in the comment of its test case before the run.
//
// Bound for traced coordinates: the coordinates of these systems are at most about 200 mm, and a
// trace of a few events rounds each coordinate a few times by u = 1.1e-16 relative, so a ray
// coordinate on the image surface is exact to about 10 * 1.1e-16 * 200 mm = 2.2e-13 mm. The
// bound kEps = 1e-12 mm is about five times that; a contribution (Phi(+T) - Phi(-T)) / 2 has the
// same bound.

#include <oneapi/tbb/task_arena.h>

#include <atomic>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/tolerance/sensitivity.hpp"
#include "rtt/tolerance/trial.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sources.hpp"
#include "sensitivity_detail.hpp"

using rtt::model::Criterion;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::EventKind;
using rtt::model::FieldType;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::RayCoordinate;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::model::ToleranceDistribution;
using rtt::model::ToleranceItem;
using rtt::model::ToleranceKind;
using rtt::tolerance::CriterionKind;
using rtt::tolerance::CriterionSensitivity;
using rtt::tolerance::ItemContribution;
using rtt::tolerance::SensitivityStatus;
using rtt::tolerance::ToleranceSensitivity;
using rtt::tolerance::TrialContext;
using rtt::tolerance::TrialStatus;

namespace {

constexpr double kEps = 1e-12;  // mm, see the head of the file

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

bool within(const std::optional<double>& value, double expected, double bound = kEps) {
  return value.has_value() && std::abs(*value - expected) <= bound;
}

Surface plane(const std::string& id, double z_mm, double aperture_radius) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  s.aperture = rtt::model::CircularAperture{aperture_radius, 0.0};
  return s;
}

/// In VACUUM at 0.5876 um, EPD 10 mm: the stop "stop" (STO, r = 10 mm) at z = -1, the ideal lens
/// "lens" (IL, f, r = 30 mm) at z = 0. Object at infinity with the fields 0 and 0.5 degree, or at
/// `object_distance` in front of the lens (also the lens's design conjugate) with the object
/// heights 0 and 1 mm. Without a detector and without paths.
System lens_system(double f, std::optional<double> object_distance = std::nullopt) {
  System s;
  s.name = "sensitivity";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  if (object_distance) {
    s.object.at_infinity = false;
    s.object.distance = Param(*object_distance);
    s.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}, {0.0, 1.0, 1.0}}};
  } else {
    s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 0.5, 1.0}}};
  }
  s.root.name = "root";
  Surface lens = plane("IL", 0.0, 30.0);
  std::optional<Param> conjugate;
  if (object_distance) conjugate = Param(*object_distance);
  lens.interaction = rtt::model::IdealLens{Param(f), conjugate};
  s.root.children = {
      {Element{
          "stop", ElementKind::Stop, Pose::along_z(-1.0), std::nullopt, {plane("STO", 0.0, 10.0)}}},
      {Element{"lens", ElementKind::ThinElement, Pose::along_z(0.0), std::nullopt, {lens}}}};
  return s;
}

/// Adds the detector "detector" (IMG) at z and the automatic path "main".
void add_detector(System& s, double z, double aperture_radius = 50.0) {
  s.root.children.push_back({Element{"detector",
                                     ElementKind::Detector,
                                     Pose::along_z(z),
                                     std::nullopt,
                                     {plane("IMG", 0.0, aperture_radius)}}});
  s.paths = {{"main", true, {}}};
}

ToleranceItem item(const std::string& id,
                   const std::string& node,
                   ToleranceKind kind,
                   double limit,
                   ToleranceDistribution distribution = ToleranceDistribution::Uniform) {
  ToleranceItem i;
  i.id = id;
  i.target.anchor = rtt::model::ToleranceAnchor::Node;
  i.target.name = node;
  i.kind = kind;
  i.limit = limit;
  i.distribution = distribution;
  return i;
}

/// A ray coordinate on `surface` of `path` (ADR 0030, point 3) with |Phi - Phi_0| <= max_change.
Criterion ray(const std::string& id,
              RayCoordinate coordinate,
              std::uint16_t field,
              double py,
              double max_change,
              const std::string& path = "main",
              const std::string& surface = "IMG") {
  rtt::model::RayOperand op;
  op.coordinate = coordinate;
  op.path = path;
  op.surface = SurfaceId(surface);
  op.field = field;
  op.py = py;
  Criterion c;
  c.id = id;
  c.operand = op;
  c.max_change = max_change;
  return c;
}

const rtt::material::MaterialLibrary& materials() {
  static const rtt::material::MaterialLibrary library;
  return library;
}

TrialContext context(const System& s) {
  return rtt::tolerance::prepare(s, materials(), nullptr, s.tolerances.items);
}

/// The value of `x`; fails the test case if there is none. clang-tidy does not take REQUIRE as a
/// check of an optional (#277), so the tests read optionals only through this.
template <typename T>
const T& required(const std::optional<T>& x) {
  const T* p = x ? &*x : nullptr;
  REQUIRE(p != nullptr);
  return *p;
}

/// The contribution of item i to criterion k.
const ItemContribution& at(const ToleranceSensitivity& r, std::size_t k, std::size_t i) {
  REQUIRE(k < r.criteria.size());
  REQUIRE(i < r.criteria[k].items.size());
  return r.criteria[k].items[i];
}

/// A linear contribution Phi = Phi_0 + slope x: plus, minus, delta and worst as derived, shift 0.
void check_linear(const ItemContribution& c, double delta) {
  INFO("delta " << c.delta.value_or(-1e300) << " against " << delta);
  CHECK(c.excluded.empty());
  CHECK(within(c.delta, delta));
  CHECK(within(c.shift, 0.0));
  CHECK(within(c.plus, delta));
  CHECK(within(c.minus, -delta));
  REQUIRE(c.worst.has_value());
  CHECK(std::abs(std::abs(required(c.worst)) - std::abs(delta)) <= 2.0 * kEps);
  CHECK(!c.q.has_value());
}

/// A linear contribution about Phi_0: Phi(+-T) = Phi_0 +- delta, shift 0.
void check_linear_about(const ItemContribution& c, double phi0, double delta) {
  INFO("delta " << c.delta.value_or(-1e300) << " against " << delta);
  CHECK(c.excluded.empty());
  CHECK(within(c.delta, delta));
  CHECK(within(c.shift, 0.0));
  CHECK(within(c.plus, phi0 + delta));
  CHECK(within(c.minus, phi0 - delta));
}

}  // namespace

TEST_CASE("sensitivity: c = sigma / T of the three distributions", "[tolerance][sensitivity]") {
  // Truncated normal (ADR 0036, point 5): sigma = T / 2 truncated at +-2 sigma, variance
  // (1 - 2 a phi(a) / Z) sigma^2 with a = 2, phi(a) = exp(-a^2 / 2) / sqrt(2 pi) and
  // Z = erf(a / sqrt(2)) (NIST DLMF, Eqs. (7.2.1), (7.2.2)), so c = sqrt(1 - 2 a phi(a) / Z) / 2.
  // Computed here at run time with the same libm (AGENTS.md, #218); 0.43981 only as a check.
  const volatile double a_in = 2.0;
  const double a = a_in;
  const double phi = std::exp(-a * a / 2.0) / std::sqrt(2.0 * std::numbers::pi);
  const double z = std::erf(a / std::numbers::sqrt2);
  const double c = std::sqrt(1.0 - 2.0 * a * phi / z) / 2.0;
  REQUIRE(std::abs(c - 0.43981) <= 1e-5);
  CHECK(std::abs(rtt::tolerance::sigma_ratio(ToleranceDistribution::Normal) - c) <= 1e-15);
  // Uniform on [-T, T]: sigma = T / sqrt(3) (NIST e-Handbook, Sec. 1.3.6.6.2); endpoint: |x| = T.
  CHECK(same(rtt::tolerance::sigma_ratio(ToleranceDistribution::Uniform), 1.0 / std::sqrt(3.0)));
  CHECK(rtt::tolerance::sigma_ratio(ToleranceDistribution::Endpoint) == 1.0);
  CHECK_THROWS_AS(rtt::tolerance::sigma_ratio(static_cast<ToleranceDistribution>(7)),
                  std::invalid_argument);
}

TEST_CASE("sensitivity: spot_rms and opd_rms are quadratic, all other operands linear",
          "[tolerance][sensitivity]") {
  using rtt::tolerance::criterion_kind;
  CHECK(criterion_kind(rtt::model::SpotRmsOperand{}) == CriterionKind::Quadratic);
  CHECK(criterion_kind(rtt::model::OpdRmsOperand{}) == CriterionKind::Quadratic);
  CHECK(criterion_kind(rtt::model::RayOperand{}) == CriterionKind::Linear);
  CHECK(criterion_kind(rtt::model::FirstOrderOperand{}) == CriterionKind::Linear);
  CHECK(criterion_kind(rtt::model::ParamValueOperand{}) == CriterionKind::Linear);
}

TEST_CASE("sensitivity: d(image position)/d(decenter) = 1 - m", "[tolerance][sensitivity]") {
  // The ideal lens maps collinearly (ADR 0031): in its own frame the object point y_o goes to
  // m y_o exactly. Decentered by delta along y, the lens sees the axial object point at -delta
  // and images it to -m delta in its frame, delta (1 - m) globally. The criterion is ray_y of the
  // chief ray on the detector in the image plane; every ray of the object point meets the image
  // point, so the aiming does not matter. Phi is linear in delta and odd: delta_i = (1 - m) T,
  // shift 0, Phi_0 = 0, all within kEps.
  SECTION("object at 2 f: m = -1, slope 2") {
    // f = 50, object 100 mm in front of the lens, image 100 mm behind it (Gauss: 1/s' = 1/f -
    // 1/s with s = 100, s' = 100), m = -s'/s = -1.
    System s = lens_system(50.0, 100.0);
    add_detector(s, 100.0);
    s.tolerances.items = {item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1)};
    s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 0, 0.0, 1.0)};
    const TrialContext ctx = context(s);
    const ToleranceSensitivity r = rtt::tolerance::sensitivity(ctx);
    CHECK(r.status == SensitivityStatus::Completed);
    REQUIRE(r.criteria.size() == 1);
    CHECK(r.criteria[0].kind == CriterionKind::Linear);
    CHECK(within(r.criteria[0].nominal, 0.0));
    check_linear(at(r, 0, 0), 2.0 * 0.1);
    // No compensators: the trials carry no LM data.
    REQUIRE(r.items.size() == 1);
    const rtt::tolerance::TrialResult& plus = required(r.items[0].plus);
    CHECK(plus.status == TrialStatus::Valid);
    CHECK(plus.compensators.empty());
    CHECK(!plus.lm_status.has_value());
  }
  SECTION("object at infinity: m = 0, slope 1") {
    System s = lens_system(100.0);
    add_detector(s, 100.0);
    s.tolerances.items = {item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1)};
    s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 0, 0.0, 1.0)};
    const ToleranceSensitivity r = rtt::tolerance::sensitivity(context(s));
    REQUIRE(r.criteria.size() == s.tolerances.criteria.size());
    CHECK(within(r.criteria[0].nominal, 0.0));
    check_linear(at(r, 0, 0), 0.1);
  }
}

TEST_CASE("sensitivity: a tilted plane-parallel plate shifts the chief ray",
          "[tolerance][sensitivity]") {
  // Plate t = 10 mm in CONST:1.5 behind the ideal lens (f = 100, object at infinity), tilted
  // about x by theta. The lateral shift of a ray at incidence theta is t sin(theta - theta') /
  // cos(theta') (docs/quellen.md, "Versatz der planparallelen Platte", from de Greve, Eq. (22), law
  // of refraction). With sin theta' = sin theta / n, cos theta' = sqrt(n^2 - sin^2 theta) / n:
  //   t (sin theta cos theta' - cos theta sin theta') / cos theta'
  //   = t sin theta [1 - cos theta / sqrt(n^2 - sin^2 theta)] =: s(theta).
  // Sign: a positive tilt about x turns +y towards +z (ADR 0036, point 1b), so the plate normal
  // becomes (0, -sin theta, cos theta). Inside the plate the ray direction is cos theta' n +
  // sin theta' t with the tangent t = (0, cos theta, sin theta); its y component sin(theta' -
  // theta) is negative, so the chief ray (on the axis, through the centre of the stop in front
  // of the lens and undeviated by it) leaves the plate at -s(theta): Phi(theta) = -s(theta), odd.
  // The shift of a parallel ray does not depend on where the plate sits, so the pivot does not
  // matter. delta_i = -s(T), shift 0.
  System s = lens_system(100.0);
  const double t = 10.0;
  const double n = 1.5;
  Surface front = plane("P.S1", 0.0, 20.0);
  Surface back = plane("P.S2", t, 20.0);
  s.root.children.push_back(
      {Element{"plate", ElementKind::Plate, Pose::along_z(40.0), "CONST:1.5", {front, back}}});
  add_detector(s, 100.0);
  const double limit_deg = 2.0;
  s.tolerances.items = {item("plate.tx", "plate", ToleranceKind::TiltXDeg, limit_deg)};
  s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 0, 0.0, 1.0)};
  const ToleranceSensitivity r = rtt::tolerance::sensitivity(context(s));
  REQUIRE(r.criteria.size() == s.tolerances.criteria.size());
  const double theta = limit_deg * std::numbers::pi / 180.0;
  const double shift =
      t * std::sin(theta) *
      (1.0 - std::cos(theta) / std::sqrt(n * n - std::sin(theta) * std::sin(theta)));
  REQUIRE(shift > 0.05);  // content guard: about t theta (1 - 1/n) = 0.058 mm
  CHECK(within(r.criteria[0].nominal, 0.0));
  check_linear(at(r, 0, 0), -shift);
}

TEST_CASE("sensitivity: double pass", "[tolerance][sensitivity]") {
  SECTION("a mirror tilted by theta deflects by 2 theta") {
    // Ideal lens f = 100 at z = 0, plane mirror M at z = 50, detector BIMG at z = 25 on the way
    // back; path STO, IL, M, BIMG. The chief ray of field 0 runs along the axis to the mirror
    // vertex. Mirror tilted about x through its vertex (the pivot of its pose, ADR 0035, point 2)
    // by theta: normal n = (0, -sin theta, cos theta) (ADR 0036, point 1b). Reflection of
    // d = (0, 0, 1), d' = d - 2 (d . n) n (de Greve, Eq. (13); docs/quellen.md), gives
    // d' = (0, sin 2 theta, -cos 2 theta): deflected by 2 theta. After L = 25 mm back along -z
    // the ray is at y = L tan 2 theta (BIMG is not rotated, so its local y is the global y).
    // Odd in theta: delta_i = L tan 2T, shift 0.
    System s = lens_system(100.0);
    Surface mirror = plane("M", 0.0, 20.0);
    mirror.interaction = rtt::model::IdealMirror{};
    s.root.children.push_back(
        {Element{"mirror", ElementKind::Mirror, Pose::along_z(50.0), std::nullopt, {mirror}}});
    s.root.children.push_back({Element{"back",
                                       ElementKind::Detector,
                                       Pose::along_z(25.0),
                                       std::nullopt,
                                       {plane("BIMG", 0.0, 50.0)}}});
    s.paths = {{"double pass",
                false,
                {{SurfaceId{"STO"}, EventKind::Transmit, 0},
                 {SurfaceId{"IL"}, EventKind::Transmit, 0},
                 {SurfaceId{"M"}, EventKind::Reflect, 0},
                 {SurfaceId{"BIMG"}, EventKind::Transmit, 0}}}};
    const double limit_deg = 0.5;
    s.tolerances.items = {item("mirror.tx", "mirror", ToleranceKind::TiltXDeg, limit_deg)};
    s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 0, 0.0, 1.0, "double pass", "BIMG")};
    const ToleranceSensitivity r = rtt::tolerance::sensitivity(context(s));
    REQUIRE(r.criteria.size() == s.tolerances.criteria.size());
    const double theta = limit_deg * std::numbers::pi / 180.0;
    CHECK(within(r.criteria[0].nominal, 0.0));
    check_linear(at(r, 0, 0), 25.0 * std::tan(2.0 * theta));
  }
  SECTION("the lens passed twice sees the same decenter") {
    // Ideal lens f = 100 at z = 0, plane mirror at d = 50, detector BIMG at z = -50; path STO,
    // IL, M, IL, BIMG (unfolded: two equal lenses 2 d apart, as in rtt-paraxial
    // test_ideal_lens_paraxial.cpp). The ideal lens maps collinearly (ADR 0031): a ray at
    // height y with slope t along the propagation leaves it with slope t - (y - delta) / f, delta
    // the decenter of the lens. Chief ray of field 0, slopes along the propagation, back pass
    // with u = -z (y unchanged):
    //   pass 1: y = 0, t = 0 -> t1 = delta / f; at the mirror y = d delta / f;
    //   plane mirror: the slope along the propagation stays; back at the lens y = 2 d delta / f;
    //   pass 2: t2 = delta / f - (2 d delta / f - delta) / f = 2 delta / f - 2 d delta / f^2;
    //   at BIMG, e = 50 behind the lens: y = 2 d delta / f + e t2.
    // With f = 100, d = 50, e = 50: y = delta (1 + 50 (0.02 - 0.01)) = 1.5 delta. Linear in delta:
    // delta_i = 1.5 T. Only the first pass decentered would give delta (slope 1), only the second
    // e delta / f = 0.5 delta: slope 1.5 shows that both passes see the same perturbation.
    System s = lens_system(100.0);
    Surface mirror = plane("M", 0.0, 20.0);
    mirror.interaction = rtt::model::IdealMirror{};
    s.root.children.push_back(
        {Element{"mirror", ElementKind::Mirror, Pose::along_z(50.0), std::nullopt, {mirror}}});
    s.root.children.push_back({Element{"back",
                                       ElementKind::Detector,
                                       Pose::along_z(-50.0),
                                       std::nullopt,
                                       {plane("BIMG", 0.0, 50.0)}}});
    s.paths = {{"double pass",
                false,
                {{SurfaceId{"STO"}, EventKind::Transmit, 0},
                 {SurfaceId{"IL"}, EventKind::Transmit, 0},
                 {SurfaceId{"M"}, EventKind::Reflect, 0},
                 {SurfaceId{"IL"}, EventKind::Transmit, 0},
                 {SurfaceId{"BIMG"}, EventKind::Transmit, 0}}}};
    s.tolerances.items = {item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1)};
    s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 0, 0.0, 1.0, "double pass", "BIMG")};
    const ToleranceSensitivity r = rtt::tolerance::sensitivity(context(s));
    REQUIRE(r.criteria.size() == s.tolerances.criteria.size());
    CHECK(within(r.criteria[0].nominal, 0.0));
    check_linear(at(r, 0, 0), 1.5 * 0.1);
  }
}

TEST_CASE("sensitivity: RSS, prediction and the inverse of a linear criterion are exact",
          "[tolerance][sensitivity]") {
  // The 2 f system of the first test (m = -1): lens decenter T1 = 0.1 (slope 2, uniform) and
  // detector decenter T2 = 0.05 (slope -1: the image point stays, the detector's local y moves
  // by -delta; endpoint). Criterion ray_y of the chief ray with max_change A = 0.03.
  //   delta_1 = 0.2, delta_2 = -0.05; RSS sqrt(0.2^2 + 0.05^2);
  //   prediction sqrt((c_1 delta_1)^2 + (c_2 delta_2)^2) with c_1 = 1/sqrt(3), c_2 = 1;
  //   inverse: A = max_change (Phi_0 = 0 lies inside), share A / sqrt(2) (n = 2),
  //   T'_i = T_i share / |delta_i|: T'_1 = share / 2, T'_2 = share;
  //   confirmation at +-T'_i: |delta_i| = share each, the RSS sqrt(2) share = A.
  System s = lens_system(50.0, 100.0);
  add_detector(s, 100.0);
  s.tolerances.items = {
      item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1, ToleranceDistribution::Uniform),
      item("det.dy", "detector", ToleranceKind::DecenterY, 0.05, ToleranceDistribution::Endpoint)};
  const double allowance = 0.03;
  s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 0, 0.0, allowance)};
  const TrialContext ctx = context(s);
  const ToleranceSensitivity r = rtt::tolerance::sensitivity(ctx);
  REQUIRE(r.status == SensitivityStatus::Completed);
  REQUIRE(r.criteria.size() == 1);
  const CriterionSensitivity& c = r.criteria[0];
  check_linear(at(r, 0, 0), 0.2);
  check_linear(at(r, 0, 1), -0.05);
  CHECK(c.complete);
  CHECK(c.excluded.empty());
  CHECK(within(c.rss, std::sqrt(0.2 * 0.2 + 0.05 * 0.05), 2.0 * kEps));
  REQUIRE(r.items.size() == 2);
  CHECK(same(r.items[0].sigma_ratio, 1.0 / std::sqrt(3.0)));
  CHECK(r.items[1].sigma_ratio == 1.0);
  const double c1 = 0.2 / std::sqrt(3.0);
  CHECK(within(c.predicted_sigma, std::sqrt(c1 * c1 + 0.05 * 0.05), 2.0 * kEps));
  CHECK(!c.predicted_mean_square.has_value());

  const double share = allowance / std::sqrt(2.0);
  CHECK(c.inverse.allowance == allowance);
  CHECK(c.inverse.note.empty());
  CHECK(within(c.inverse.share, share, 1e-15));
  REQUIRE(c.inverse.limits.size() == 2);
  // T' = T share / |delta| with |delta| exact to kEps: relative error kEps / |delta| <= 2e-11.
  CHECK(within(c.inverse.limits[0], share / 2.0, 1e-12));
  CHECK(within(c.inverse.limits[1], share, 1e-12));
  CHECK(c.inverse.unlimited.empty());
  REQUIRE(r.limits.size() == 2);
  CHECK(r.limits[0] == c.inverse.limits[0]);
  CHECK(r.limits[1] == c.inverse.limits[1]);

  const rtt::tolerance::Confirmation& k = required(r.confirmation);
  REQUIRE(k.limits.size() == 2);
  CHECK(k.limits[0] == required(r.limits[0]));
  CHECK(k.limits[1] == required(r.limits[1]));
  REQUIRE(k.criteria.size() == 1);
  REQUIRE(k.criteria[0].items.size() == 2);
  // |delta'_i| = slope T'_i: kEps from the trace plus |slope| 1e-12 from T'.
  CHECK(within(k.criteria[0].items[0].delta, share, 3e-12));
  CHECK(within(k.criteria[0].items[1].delta, -share, 3e-12));
  CHECK(within(k.criteria[0].rss, allowance, 5e-12));
  CHECK(k.criteria[0].inverse.limits.empty());

  SECTION("without confirmation") {
    rtt::tolerance::SensitivityOptions options;
    options.confirm = false;
    const ToleranceSensitivity without = rtt::tolerance::sensitivity(ctx, options);
    CHECK(!without.confirmation.has_value());
    REQUIRE(without.criteria.size() == 1);
    CHECK(without.criteria[0].inverse.limits == c.inverse.limits);
  }
}

TEST_CASE("sensitivity: a quadratic criterion, q = k^2 T^2 exactly", "[tolerance][sensitivity]") {
  // Ideal lens f = 100, object at infinity, detector in the focal plane, tolerance shift_z of
  // the detector by x; no compensators (a focus compensator would undo it). Every ray leaves the
  // lens at height h towards the focus (0, 0, f): at z = f + x it is at -h x / f. The spot
  // (centroid 0 by symmetry) has RMS radius |x| rho / f with rho^2 the mean of |h|^2 over the
  // hexapolar pupil (6 rings, equal weights: one wavelength, no losses in VACUUM). With real
  // aiming the stop (the entrance pupil, nothing in front of it) is met at h = r_EP (px, py),
  // r_EP = 5 mm, and the parallel bundle keeps this height up to the lens.
  //   So Phi(x) = k |x| with k = rho / f, Phi_0 = 0, and q = (Phi(T)^2 + Phi(-T)^2) / 2 = k^2 T^2.
  // Bound (derived, not measured): Phi_{+-} = k T + e_{+-}, |e| <= kEps, |Phi_0| <= kEps (REQUIRE),
  // so |q - k^2 T^2| <= 2 k T kEps + 2 kEps^2.
  System s = lens_system(100.0);
  add_detector(s, 100.0);
  const double limit = 0.5;
  s.tolerances.items = {item("det.dz", "detector", ToleranceKind::ShiftZ, limit)};
  const double upper = 0.01;
  rtt::model::SpotRmsOperand spot;
  spot.path = "main";
  spot.field = 0;
  spot.rings = 6;
  spot.reference = rtt::model::SpotReference::Centroid;
  Criterion criterion;
  criterion.id = "spot";
  criterion.operand = spot;
  criterion.max = upper;
  s.tolerances.criteria = {criterion};

  // rho from the pupil coordinates of the hexapolar sampling.
  const rtt::compile::CompiledSystem cs = rtt::compile::compile(s, materials());
  const std::uint16_t fields[] = {0};
  const auto path = cs.find_path("main");
  const rtt::trace::RayBatch pupil =
      rtt::trace::make_rays(cs, required(path), fields, 0, rtt::trace::HexapolarPupil{6});
  REQUIRE(pupil.size() > 100);
  double sum = 0.0;
  for (std::size_t j = 0; j < pupil.size(); ++j) {
    sum += pupil.pupil_x()[j] * pupil.pupil_x()[j] + pupil.pupil_y()[j] * pupil.pupil_y()[j];
  }
  const double rho = 5.0 * std::sqrt(sum / static_cast<double>(pupil.size()));
  const double k = rho / 100.0;

  const ToleranceSensitivity r = rtt::tolerance::sensitivity(context(s));
  REQUIRE(r.criteria.size() == s.tolerances.criteria.size());
  const CriterionSensitivity& c = r.criteria[0];
  CHECK(c.kind == CriterionKind::Quadratic);
  REQUIRE(c.nominal.has_value());
  REQUIRE(std::abs(required(c.nominal)) <= kEps);
  const ItemContribution& q = at(r, 0, 0);
  CHECK(within(q.plus, k * limit));
  CHECK(within(q.minus, k * limit));
  const double bound = 2.0 * k * limit * kEps + 2.0 * kEps * kEps;
  CHECK(within(q.q, k * k * limit * limit, bound));
  CHECK(!q.delta.has_value());
  CHECK(!q.shift.has_value());
  // RSS: Phi_RSS^2 = Phi_0^2 + q; prediction Phi_0^2 + c^2 q with c^2 = 1/3 (uniform).
  CHECK(within(c.rss, k * limit, 2.0 * kEps));
  CHECK(within(c.predicted_mean_square, k * k * limit * limit / 3.0, bound));
  CHECK(!c.predicted_sigma.has_value());
  // Inverse: U = max, A^2 = U^2 - Phi_0^2, share A^2 (n = 1), T' = T sqrt(share / q) = U / k up
  // to the relative error of q, bound / (k T)^2, and Phi_0^2 / U^2 <= 1e-20.
  const double relative = bound / (k * k * limit * limit) + 1e-15;
  REQUIRE(c.inverse.allowance.has_value());
  CHECK(std::abs(required(c.inverse.allowance) - upper) <= 1e-15);
  REQUIRE(c.inverse.limits.size() == 1);
  CHECK(within(c.inverse.limits[0], upper / k, relative * upper / k));
  // Confirmation at T' = U / k: Phi = k T' = U.
  const rtt::tolerance::Confirmation& confirmation = required(r.confirmation);
  REQUIRE(confirmation.criteria.size() == 1);
  CHECK(within(confirmation.criteria[0].rss, upper, relative * upper + 2.0 * kEps));
}

TEST_CASE("sensitivity: items that limit nothing, and the smallest limit over the criteria",
          "[tolerance][sensitivity]") {
  // Ideal lens f = 100, object at infinity, detector at the focus. Items: lens decenter along
  // y (T = 0.1) and along x (T = 0.1); both have slope 1 on the coordinate they move.
  // Criteria: chief ray_y (max_change 0.04), chief ray_x (0.01), ray_y of the ray at py = 0.5
  // (0.02; it meets the focus, so the same slope 1).
  // A decenter along y moves nothing out of the meridional plane, and a decenter along x leaves
  // the y coordinates alone: for those pairs Phi(+T) = Phi(-T) bitwise (mirror symmetry) and
  // delta = 0,
  // so the item limits nothing for that criterion (ADR 0036, point 6).
  //   y item: T' = share = A / sqrt(2) for both y criteria, the smallest is 0.02 / sqrt(2);
  //   x item: limited only by ray_x, T' = 0.01 / sqrt(2).
  System s = lens_system(100.0);
  add_detector(s, 100.0);
  s.tolerances.items = {item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1),
                        item("lens.dx", "lens", ToleranceKind::DecenterX, 0.1)};
  s.tolerances.criteria = {ray("chief_y", RayCoordinate::Y, 0, 0.0, 0.04),
                           ray("chief_x", RayCoordinate::X, 0, 0.0, 0.01),
                           ray("zone_y", RayCoordinate::Y, 0, 0.5, 0.02)};
  const ToleranceSensitivity r = rtt::tolerance::sensitivity(context(s));
  REQUIRE(r.criteria.size() == s.tolerances.criteria.size());
  REQUIRE(r.criteria.size() == 3);
  check_linear(at(r, 0, 0), 0.1);
  check_linear(at(r, 1, 1), 0.1);
  check_linear(at(r, 2, 0), 0.1);
  CHECK(at(r, 1, 0).delta == 0.0);
  CHECK(at(r, 0, 1).delta == 0.0);
  CHECK(at(r, 2, 1).delta == 0.0);
  CHECK(r.criteria[0].inverse.unlimited == std::vector<std::size_t>{1});
  CHECK(r.criteria[1].inverse.unlimited == std::vector<std::size_t>{0});
  CHECK(r.criteria[2].inverse.unlimited == std::vector<std::size_t>{1});
  CHECK(!r.criteria[1].inverse.limits[0].has_value());
  CHECK(!r.criteria[0].inverse.limits[1].has_value());
  REQUIRE(r.limits.size() == 2);
  CHECK(within(r.limits[0], 0.02 / std::sqrt(2.0), 1e-12));
  CHECK(within(r.limits[1], 0.01 / std::sqrt(2.0), 1e-12));
  const rtt::tolerance::Confirmation& confirmation = required(r.confirmation);
  REQUIRE(confirmation.limits.size() == 2);
  CHECK(confirmation.limits[0] == required(r.limits[0]));
  CHECK(confirmation.limits[1] == required(r.limits[1]));
}

TEST_CASE("sensitivity: an item that acts only in second order limits nothing",
          "[tolerance][sensitivity]") {
  // ADR 0036, addendum #258 (review 13, P2-1). The plate of the plate test (t = 10 mm,
  // CONST:1.5, behind the ideal lens f = 100) tilted about y. The ray at py = 0.5 lies in the
  // y-z plane; the system with the plate at -theta is the mirror image (x -> -x) of the one at
  // +theta, and the mirror maps this ray onto itself, so Phi(theta) = Phi(-theta): ray_y is even
  // in theta, dPhi/dtheta = 0 at the nominal point. After rounding delta_i is noise, at most a
  // few ulp of the coordinate (|y| < 1 mm here), far below eps_M max(|Phi|, A) = 1e-9 * 0.01;
  // the item must be reported as limiting nothing, not get a T' of kilometres.
  // Content guard: the second-order change shift_i is real. The ray crosses the plate at the
  // slope h / f = 0.025; tilted by theta about y its path in the glass changes by a relative
  // O(theta^2) = 1.2e-3 (theta = 2 deg), so the lateral offset of about t (1 - 1/n) 0.025 =
  // 0.08 mm changes by about 1e-4 mm. Guard |shift_i| > 1e-7 mm, far above the rounding.
  System s = lens_system(100.0);
  Surface front = plane("P.S1", 0.0, 20.0);
  Surface back = plane("P.S2", 10.0, 20.0);
  s.root.children.push_back(
      {Element{"plate", ElementKind::Plate, Pose::along_z(40.0), "CONST:1.5", {front, back}}});
  add_detector(s, 100.0);
  s.tolerances.items = {item("plate.ty", "plate", ToleranceKind::TiltYDeg, 2.0)};
  s.tolerances.criteria = {ray("zone", RayCoordinate::Y, 0, 0.5, 0.01)};
  const ToleranceSensitivity r = rtt::tolerance::sensitivity(context(s));
  REQUIRE(r.criteria.size() == s.tolerances.criteria.size());
  const ItemContribution& c = at(r, 0, 0);
  REQUIRE(c.delta.has_value());
  REQUIRE(c.shift.has_value());
  CHECK(std::abs(c.shift.value_or(0.0)) > 1e-7);
  CHECK(std::abs(c.delta.value_or(1.0)) <= 1e-9 * 0.01);
  CHECK(r.criteria[0].inverse.unlimited == std::vector<std::size_t>{0});
  CHECK(!r.criteria[0].inverse.limits[0].has_value());
  REQUIRE(r.limits.size() == 1);
  CHECK(!r.limits[0].has_value());
  // Not limited: the confirmation keeps T_i.
  CHECK(required(r.confirmation).limits.at(0) == 2.0);
}

TEST_CASE("sensitivity: a side that is not defined excludes the item for that criterion",
          "[tolerance][sensitivity]") {
  // Addendum #258, F1. Ideal lens f = 100, object at infinity, detector at the focus with a
  // clear aperture of radius 1 mm. The chief ray of field 1 (0.5 degree) meets it at
  // y_1 = 100 tan(0.5 deg) = 0.873 mm. Detector decenter T = 0.3: its local y is y_1 - delta,
  // 0.573 at +T, 1.173 > 1 at -T (vignetted, not defined). Lens decenter T = 0.1: slope 1,
  // y_1 +- 0.1 stays inside. Criteria: chief ray_y of field 1 (A) and of field 0 (B).
  //   A: detector item excluded with the reason of the -T side, RSS incomplete = |delta_lens|;
  //   B: both items, RSS sqrt(0.3^2 + 0.1^2).
  System s = lens_system(100.0);
  add_detector(s, 100.0, 1.0);
  s.tolerances.items = {item("det.dy", "detector", ToleranceKind::DecenterY, 0.3),
                        item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1)};
  s.tolerances.criteria = {ray("field1", RayCoordinate::Y, 1, 0.0, 0.5),
                           ray("field0", RayCoordinate::Y, 0, 0.0, 0.5)};
  const ToleranceSensitivity r = rtt::tolerance::sensitivity(context(s));
  REQUIRE(r.criteria.size() == s.tolerances.criteria.size());
  REQUIRE(r.status == SensitivityStatus::Completed);
  REQUIRE(r.items.size() == 2);
  const rtt::tolerance::TrialResult* minus = &required(r.items[0].minus);
  CHECK(minus->status == TrialStatus::Invalid);
  REQUIRE(minus->criteria.size() == 2);
  CHECK(!minus->criteria[0].value.has_value());
  CHECK(minus->criteria[1].value.has_value());

  const CriterionSensitivity& a = r.criteria[0];
  CHECK(!at(r, 0, 0).excluded.empty());
  CHECK(!at(r, 0, 0).delta.has_value());
  CHECK(a.excluded == std::vector<std::size_t>{0});
  CHECK(!a.complete);
  REQUIRE(a.nominal.has_value());  // y_1 = 0.873, not 0
  check_linear_about(at(r, 0, 1), a.nominal.value_or(0.0), 0.1);
  CHECK(within(a.rss, 0.1));
  // The excluded item gets no limit, and still counts in n = 2 (addendum #258).
  CHECK(within(a.inverse.share, 0.5 / std::sqrt(2.0), 1e-15));
  CHECK(!a.inverse.limits[0].has_value());
  CHECK(a.inverse.limits[1].has_value());

  const CriterionSensitivity& b = r.criteria[1];
  CHECK(b.complete);
  check_linear(at(r, 1, 0), -0.3);
  check_linear(at(r, 1, 1), 0.1);
  CHECK(within(b.rss, std::sqrt(0.3 * 0.3 + 0.1 * 0.1), 2.0 * kEps));
}

TEST_CASE("sensitivity: a side not defined at T' excludes the item in the confirmation",
          "[tolerance][sensitivity]") {
  // Addendum #258, F1 with the T' of the inverse (ADR 0036, addendum #257: values beyond T_i are
  // legitimate). Ideal lens f = 100, object at infinity, detector at the focus with a clear
  // aperture of radius 1 mm; the chief ray of field 1 meets it at y_1 = 100 tan(0.5 deg) =
  // 0.873 mm. Detector decenter T = 0.05 (slope -1, both sides inside: 0.823, 0.923). Criterion
  // chief ray_y of field 1 with max_change A = 0.5: n = 1, share A, T' = T A / 0.05 = 0.5. In
  // the confirmation the -T' side puts the ray at 1.373 > 1: vignetted, not defined. The item
  // is excluded there with the reason of that side; the trial is Invalid, nothing throws.
  System s = lens_system(100.0);
  add_detector(s, 100.0, 1.0);
  s.tolerances.items = {item("det.dy", "detector", ToleranceKind::DecenterY, 0.05)};
  s.tolerances.criteria = {ray("field1", RayCoordinate::Y, 1, 0.0, 0.5)};
  const ToleranceSensitivity r = rtt::tolerance::sensitivity(context(s));
  REQUIRE(r.criteria.size() == s.tolerances.criteria.size());
  REQUIRE(r.status == SensitivityStatus::Completed);
  REQUIRE(r.criteria[0].nominal.has_value());
  check_linear_about(at(r, 0, 0), r.criteria[0].nominal.value_or(0.0), -0.05);
  REQUIRE(r.limits.size() == 1);
  CHECK(within(r.limits[0], 0.5, 1e-10));
  const rtt::tolerance::Confirmation& k = required(r.confirmation);
  REQUIRE(k.minus.size() == 1);
  CHECK(k.plus[0].status == TrialStatus::Valid);
  CHECK(k.minus[0].status == TrialStatus::Invalid);
  REQUIRE(k.criteria.size() == 1);
  REQUIRE(k.criteria[0].items.size() == 1);
  const ItemContribution& c = k.criteria[0].items[0];
  CHECK(c.excluded.rfind("-T: ", 0) == 0);
  CHECK(c.plus.has_value());
  CHECK(!c.delta.has_value());
  CHECK(k.criteria[0].excluded == std::vector<std::size_t>{0});
  CHECK(!k.criteria[0].complete);
}

TEST_CASE("sensitivity: the quadratic RSS of a negative sum is 0 with a warning",
          "[tolerance][sensitivity]") {
  // ADR 0036, point 6: Phi_RSS^2 = Phi_0^2 + sum q_i; q_i < 0 if an item improves on average.
  // Binary fractions, so the sums below are exact.
  using rtt::tolerance::detail::quadratic_rss;
  const std::vector<double> better = {-0.75, -0.5};
  const auto negative = quadratic_rss(1.0, better);  // 1 - 1.25 = -0.25
  CHECK(negative.negative);
  CHECK(negative.value == 0.0);
  const std::vector<double> worse = {0.25, 0.5};
  const auto positive = quadratic_rss(0.5, worse);  // 0.25 + 0.75 = 1
  CHECK(!positive.negative);
  CHECK(positive.value == 1.0);
  const std::vector<double> even = {-0.25};
  const auto zero = quadratic_rss(0.5, even);  // exactly 0: no warning
  CHECK(!zero.negative);
  CHECK(zero.value == 0.0);
}

TEST_CASE("sensitivity: the overload with a System expands and keeps the expansion",
          "[tolerance][sensitivity]") {
  // ADR 0036, addendum #260: expand, prepare, the run; bitwise the same as with the context of
  // the same items (no class: expand gives the explicit items in file order).
  System s = lens_system(50.0, 100.0);
  add_detector(s, 100.0);
  s.tolerances.items = {item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1),
                        item("det.dy", "detector", ToleranceKind::DecenterY, 0.05)};
  s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 0, 0.0, 0.03)};
  const ToleranceSensitivity a = rtt::tolerance::sensitivity(s, materials(), nullptr);
  const ToleranceSensitivity b = rtt::tolerance::sensitivity(context(s));
  const rtt::tolerance::ToleranceItems* expansion = &required(a.expansion);
  REQUIRE(expansion->items.size() == 2);
  CHECK(expansion->items[0].id == "lens.dy");
  CHECK(expansion->items[1].id == "det.dy");
  CHECK(!b.expansion.has_value());
  REQUIRE(a.items.size() == 2);
  CHECK(a.items[1].item.id == "det.dy");
  REQUIRE(a.criteria.size() == 1);
  REQUIRE(a.criteria[0].rss.has_value());
  REQUIRE(b.criteria.size() == 1);
  REQUIRE(b.criteria[0].rss.has_value());
  CHECK(same(a.criteria[0].rss.value_or(0.0), b.criteria[0].rss.value_or(1.0)));
  CHECK(within(a.criteria[0].rss, std::sqrt(0.2 * 0.2 + 0.05 * 0.05), 2.0 * kEps));
  // An invalid set is an error of expand, before any trial.
  System invalid = s;
  invalid.tolerances.items[1].limit = -1.0;
  CHECK_THROWS_AS(rtt::tolerance::sensitivity(invalid, materials(), nullptr),
                  rtt::compile::CompileError);
}

TEST_CASE("sensitivity: a criterion without Phi_0 gets no sensitivity and no inverse",
          "[tolerance][sensitivity]") {
  // Addendum #258, F1: the ray at py = 2.5 meets the plane of the stop at 12.5 mm, outside its
  // radius of 10 mm, already in the nominal system. The other criterion is not affected.
  System s = lens_system(100.0);
  add_detector(s, 100.0);
  s.tolerances.items = {item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1)};
  s.tolerances.criteria = {ray("outside", RayCoordinate::Y, 0, 2.5, 0.5),
                           ray("chief", RayCoordinate::Y, 0, 0.0, 0.5)};
  const ToleranceSensitivity r = rtt::tolerance::sensitivity(context(s));
  REQUIRE(r.criteria.size() == s.tolerances.criteria.size());
  CHECK(r.status == SensitivityStatus::Completed);
  CHECK(r.nominal.status == TrialStatus::Invalid);
  const CriterionSensitivity& a = r.criteria[0];
  CHECK(!a.nominal.has_value());
  CHECK(!a.nominal_reason.empty());
  CHECK(!at(r, 0, 0).excluded.empty());
  CHECK(!at(r, 0, 0).delta.has_value());
  CHECK(!a.rss.has_value());
  CHECK(!a.predicted_sigma.has_value());
  CHECK(!a.inverse.allowance.has_value());
  CHECK(a.inverse.note == "nominal value not defined");
  check_linear(at(r, 1, 0), 0.1);
  CHECK(r.criteria[1].rss.has_value());
}

TEST_CASE("sensitivity: cancelled runs give a partial result", "[tolerance][sensitivity]") {
  // Addendum #258, F2: only items with both sides, no RSS, no inverse, no confirmation.
  System s = lens_system(50.0, 100.0);
  add_detector(s, 100.0);
  s.tolerances.items = {item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1),
                        item("det.dy", "detector", ToleranceKind::DecenterY, 0.05),
                        item("lens.dx", "lens", ToleranceKind::DecenterX, 0.1)};
  s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 0, 0.0, 0.03)};
  const TrialContext ctx = context(s);
  const auto check_partial = [](const ToleranceSensitivity& r) {
    REQUIRE(r.criteria.size() == 1);
    REQUIRE(r.limits.size() == r.items.size());
    for (const auto& limit : r.limits) CHECK(!limit.has_value());
    CHECK(r.status == SensitivityStatus::Cancelled);
    CHECK(!r.confirmation.has_value());
    for (std::size_t i = 0; i < r.items.size(); ++i) {
      INFO("item " << i);
      const bool both = r.items[i].plus.has_value() && r.items[i].minus.has_value();
      CHECK(at(r, 0, i).excluded.empty() == both);
      CHECK(at(r, 0, i).delta.has_value() == both);
    }
    CHECK(!r.criteria[0].rss.has_value());
    CHECK(!r.criteria[0].predicted_sigma.has_value());
    CHECK(!r.criteria[0].inverse.allowance.has_value());
    CHECK(r.criteria[0].inverse.limits.empty());
  };
  SECTION("before the first evaluation") {
    rtt::trace::RunControl control;
    control.cancel = rtt::trace::CancelToken{};
    control.cancel->request_cancel();
    const ToleranceSensitivity r = rtt::tolerance::sensitivity(ctx, {}, &control);
    check_partial(r);
    for (const auto& i : r.items) {
      CHECK(!i.plus.has_value());
      CHECK(!i.minus.has_value());
    }
  }
  SECTION("already in prepare (nominal() cancelled, #257)") {
    rtt::trace::RunControl control;
    control.cancel = rtt::trace::CancelToken{};
    control.cancel->request_cancel();
    const TrialContext cancelled =
        rtt::tolerance::prepare(s, materials(), nullptr, s.tolerances.items, {}, &control);
    REQUIRE(cancelled.nominal().status == TrialStatus::Cancelled);
    // No control for the run itself: it must stop because of nominal() alone.
    const ToleranceSensitivity r = rtt::tolerance::sensitivity(cancelled);
    check_partial(r);
    // The cause is the cancellation, not an undefined nominal value.
    CHECK(at(r, 0, 0).excluded == "cancelled");
    CHECK(r.criteria[0].inverse.note.empty());
    for (const auto& i : r.items) {
      CHECK(!i.plus.has_value());
      CHECK(!i.minus.has_value());
    }
  }
  SECTION("in the confirmation: the first part stays complete") {
    // Addendum #258, F2: only the confirmation is dropped. lens.dx limits nothing for ray_y
    // (decenter along x), so only the first two items get a T'.
    rtt::trace::RunControl control;
    control.cancel = rtt::trace::CancelToken{};
    control.min_interval = std::chrono::milliseconds{0};
    std::atomic<int> confirmation_reports{0};
    control.progress = [&](const rtt::trace::Progress& p) {
      if (p.stage != "confirmation") return;
      ++confirmation_reports;
      control.cancel->request_cancel();
    };
    const ToleranceSensitivity r = rtt::tolerance::sensitivity(ctx, {}, &control);
    CHECK(confirmation_reports >= 1);
    CHECK(r.status == SensitivityStatus::Cancelled);
    CHECK(!r.confirmation.has_value());
    REQUIRE(r.criteria.size() == 1);
    const CriterionSensitivity& c = r.criteria[0];
    CHECK(c.complete);
    CHECK(c.rss.has_value());
    CHECK(c.predicted_sigma.has_value());
    CHECK(c.inverse.share.has_value());
    REQUIRE(c.inverse.limits.size() == 3);
    REQUIRE(r.limits.size() == 3);
    CHECK(r.limits[0].has_value());
    CHECK(r.limits[1].has_value());
    CHECK(!r.limits[2].has_value());
    for (const auto& i : r.items) {
      CHECK(i.plus.has_value());
      CHECK(i.minus.has_value());
    }
  }
  SECTION("after the first progress report") {
    rtt::trace::RunControl control;
    control.cancel = rtt::trace::CancelToken{};
    control.min_interval = std::chrono::milliseconds{0};
    std::atomic<int> reports{0};
    control.progress = [&](const rtt::trace::Progress&) {
      if (++reports == 1) control.cancel->request_cancel();
    };
    const ToleranceSensitivity r = rtt::tolerance::sensitivity(ctx, {}, &control);
    check_partial(r);
    CHECK(reports >= 1);
  }
}

TEST_CASE("sensitivity: the compensator takes back a lens shift", "[tolerance][sensitivity]") {
  // #258: +-T with compensators. Ideal lens f = 100, object at infinity, stop at z = -1,
  // detector in the focal plane. Item: lens shift_z by x (T = 0.3). Criterion: ray_y of the chief
  // ray of field 1 (u = tan 0.5 deg) on the detector. The ideal lens maps the slope t of a ray
  // at height h to t - h / f (ADR 0031, collinear); the chief ray passes the stop centre, so it
  // meets the lens at h = u (1 + x).
  // - Without compensator the detector stays at 100:
  //     Phi(x) = h + (100 - x)(u - h / f) = u [1 + x + (100 - x)(99 - x) / 100]
  //            = 100 u - 0.99 u x + u x^2 / 100,
  //   so delta = (Phi(T) - Phi(-T)) / 2 = -0.99 u T and shift = (Phi(T) + Phi(-T)) / 2 - Phi_0 =
  //   u T^2 / 100, both exact up to the trace rounding kEps.
  // - With the detector position as compensator (travel +-2, merit: rms spot of field 0 about
  //   the chief ray) the detector moves to the focus 100 + x, and Phi = h + 100 (u - h / f) =
  //   100 u for every x: delta = 0 up to the compensator error e_p times the slope of the chief
  //   ray after the lens, below 0.01 (#257: |e_p| about 3e-8, bound there; here 1e-8 mm in
  //   total is generous). The compensator moves by +-T within 1e-6.
  System s = lens_system(100.0);
  add_detector(s, 100.0);
  s.tolerances.items = {item("lens.dz", "lens", ToleranceKind::ShiftZ, 0.3)};
  s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 1, 0.0, 0.01)};
  const double u = std::tan(0.5 * std::numbers::pi / 180.0);

  const ToleranceSensitivity open = rtt::tolerance::sensitivity(context(s));
  REQUIRE(open.criteria.size() == 1);
  const ItemContribution& o = at(open, 0, 0);
  CHECK(o.excluded.empty());
  CHECK(within(open.criteria[0].nominal, 100.0 * u));
  CHECK(within(o.delta, -0.99 * u * 0.3));
  CHECK(within(o.shift, u * 0.3 * 0.3 / 100.0, 2.0 * kEps));

  rtt::model::SpotGenerator focus;
  focus.path = "main";
  focus.fields = std::vector<std::uint16_t>{0};
  focus.reference = rtt::model::SpotReference::Chief;
  s.tolerances.compensator_merit = rtt::model::Optimization{{}, {focus}};
  rtt::model::ToleranceTarget detector;
  detector.anchor = rtt::model::ToleranceAnchor::Node;
  detector.name = "detector";
  detector.param = "/pose/position/2";
  s.tolerances.compensators = {{detector, -2.0, 2.0}};
  const ToleranceSensitivity r = rtt::tolerance::sensitivity(context(s));
  REQUIRE(r.status == SensitivityStatus::Completed);
  REQUIRE(r.criteria.size() == 1);
  const ItemContribution& c = at(r, 0, 0);
  CHECK(c.excluded.empty());
  CHECK(std::abs(c.delta.value_or(1.0)) <= 1e-8);
  CHECK(within(r.criteria[0].nominal, 100.0 * u, 1e-8));
  REQUIRE(r.items.size() == 1);
  const rtt::tolerance::TrialResult& plus = required(r.items[0].plus);
  const rtt::tolerance::TrialResult& minus = required(r.items[0].minus);
  CHECK(plus.lm_status.has_value());
  REQUIRE(plus.compensators.size() == 1);
  REQUIRE(minus.compensators.size() == 1);
  CHECK(std::abs(plus.compensators[0].end - plus.compensators[0].start - 0.3) <= 1e-6);
  CHECK(std::abs(minus.compensators[0].end - minus.compensators[0].start + 0.3) <= 1e-6);
  // The compensated item limits nothing for this criterion (delta is noise).
  CHECK(r.criteria[0].inverse.unlimited == std::vector<std::size_t>{0});
}

TEST_CASE("sensitivity: bitwise the same for 1, 4 and all threads", "[tolerance][sensitivity]") {
  System s = lens_system(50.0, 100.0);
  add_detector(s, 100.0);
  s.tolerances.items = {item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1),
                        item("det.dy", "detector", ToleranceKind::DecenterY, 0.05),
                        item("lens.dz", "lens", ToleranceKind::ShiftZ, 0.2)};
  rtt::model::SpotRmsOperand spot;
  spot.path = "main";
  spot.field = 1;
  Criterion quadratic;
  quadratic.id = "spot";
  quadratic.operand = spot;
  quadratic.max = 0.05;
  s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 1, 0.0, 0.03), quadratic};
  const TrialContext ctx = context(s);
  const auto run = [&](int threads) {
    oneapi::tbb::task_arena arena(threads);
    return arena.execute([&] { return rtt::tolerance::sensitivity(ctx); });
  };
  const auto values = [](const ToleranceSensitivity& r) {
    std::vector<double> v;
    const auto add = [&](const std::optional<double>& x) { v.push_back(x.value_or(-1e300)); };
    for (const CriterionSensitivity& c : r.criteria) {
      for (const ItemContribution& i : c.items) {
        for (const auto* x : {&i.plus, &i.minus, &i.worst, &i.delta, &i.shift, &i.q}) add(*x);
      }
      add(c.rss);
      add(c.predicted_sigma);
      add(c.predicted_mean_square);
      add(c.inverse.share);
      for (const auto& x : c.inverse.limits) add(x);
    }
    for (const auto& x : r.limits) add(x);
    if (r.confirmation) {
      for (const double x : r.confirmation->limits) v.push_back(x);
      for (const CriterionSensitivity& c : r.confirmation->criteria) add(c.rss);
    }
    return v;
  };
  const ToleranceSensitivity first = run(1);
  // Content guard: a complete run with values, not only placeholders.
  REQUIRE(first.status == SensitivityStatus::Completed);
  REQUIRE(first.criteria.size() == 2);
  CHECK(first.criteria[0].rss.has_value());
  CHECK(first.criteria[1].rss.has_value());
  CHECK(first.criteria[1].nominal.value_or(0.0) > 0.0);  // the spot of field 1 is not a point
  CHECK(first.confirmation.has_value());
  const std::vector<double> one = values(first);
  REQUIRE(one.size() > 20);
  for (const int threads : {4, oneapi::tbb::task_arena::automatic}) {
    INFO(threads << " threads");
    const std::vector<double> other = values(run(threads));
    REQUIRE(other.size() == one.size());
    bool all_same = true;
    for (std::size_t i = 0; i < one.size(); ++i) all_same = all_same && same(one[i], other[i]);
    CHECK(all_same);
  }
}
