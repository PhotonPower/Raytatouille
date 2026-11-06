// Aiming per section for anamorphic paths (#315, ADR 0042, point 3; acceptance 1 (real rays),
// 5 (start plane), 7 (errors of the aiming) and 8 (aiming)).
//
// The pupil data per section (z_ep,s, r_ep,s, stop scale) come from first_order_section and
// trace_ray_section, which test_paraxial_section.cpp checks against an independent matrix
// computation; here they only place the targets of point 3. The real rays are checked against
// the geometry: the hit on the stop surface against the target of point 3, and the axis
// crossing of the cylinder lens against an independent two-dimensional trace in this file
// (refraction in vector form, B. de Greve, "Reflections and Refractions in Ray Tracing" (2006),
// Eq. (22) with (23) and (28); docs/quellen.md, row "Brechung in Vektorform").
//
// Start rays are aimed from the object side; in the hand-made batches they start at z = -1, in
// front of the first surface (AGENTS.md).

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using Catch::Matchers::ContainsSubstring;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::math::Vec3;
using rtt::paraxial::FirstOrder;
using rtt::paraxial::Section;
using rtt::trace::Aiming;
using rtt::trace::PointListPupil;
using rtt::trace::PupilPoint;
using rtt::trace::RayBatch;
using rtt::trace::RayPaths;
using rtt::trace::RayStatus;

namespace {

constexpr double kN = 1.5168;  // CONST:1.5168
constexpr PathId kMain{0};

// ---------------------------------------------------------------- test systems (JSON) -----

/// "%.17g" reads back as the same double.
std::string num(double v) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.17g", v);
  return buf;
}

std::string pose(double z, double rz = 0.0) {
  return R"({"position": [0.0, 0.0, )" + num(z) + R"(], "rotation_deg": [0.0, 0.0, )" + num(rz) +
         "]}";
}

std::string stop(double z, double radius) {
  return R"({"type": "stop", "name": "stop", "pose": )" + pose(z) +
         R"(, "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": )" +
         num(radius) + "}}]}";
}

/// Lens "L1" at global z with thickness t: front surface L1.S1 with `shape` (JSON of
/// shape.base), rear surface plane; `aperture` JSON for both surfaces or empty for none.
std::string lens(double z,
                 const std::string& shape,
                 double t,
                 const std::string& material = "CONST:1.5168",
                 const std::string& aperture =
                     R"({"type": "rectangular", "half_width_x": 12.7, "half_width_y": 12.7})") {
  const std::string ap = aperture.empty() ? "" : R"(, "aperture": )" + aperture;
  return R"({"type": "lens", "name": "L1", "pose": )" + pose(z) + R"(, "material": ")" + material +
         R"(", "surfaces": [{"id": "L1.S1", "shape": {"base": )" + shape + "}" + ap +
         R"(}, {"id": "L1.S2", "pose": )" + pose(t) + ap + "}]}";
}

/// Thin element "IL" at global z with `interaction` (ideal lens or ideal cylinder lens).
std::string thin(double z, const std::string& interaction, const std::string& name = "lens") {
  return R"({"type": "thin_element", "name": ")" + name + R"(", "pose": )" + pose(z) +
         R"(, "surfaces": [{"id": ")" + name + R"(.IL", "aperture": {"type": "circular", )" +
         R"("radius": 30.0}, "interaction": )" + interaction + "}]}";
}

std::string cylinder_lens(double f, double psi_deg = 0.0) {
  return R"({"type": "ideal_cylinder_lens", "focal_length": )" + num(f) + R"(, "axis_deg": )" +
         num(psi_deg) + "}";
}

std::string detector(double z) {
  return R"({"type": "detector", "name": "image", "pose": )" + pose(z) +
         R"(, "surfaces": [{"id": "IMG"}]})";
}

struct Options {
  std::string aperture = R"({"type": "epd", "value": 10.0})";
  std::string fields = R"({"points": [{}]})";
  std::string object;  ///< JSON of "object", empty: at infinity
  std::string environment = "VACUUM";
  std::string wavelengths = R"([{"um": 0.5876, "reference": true}])";
};

std::string system_json(const std::vector<std::string>& children, const Options& o) {
  std::string list;
  for (std::size_t i = 0; i < children.size(); ++i) {
    if (i > 0) list += ", ";
    list += children[i];
  }
  const std::string object = o.object.empty() ? "" : R"( "object": )" + o.object + ",";
  return R"({"schema_version": "0.9.0", "name": "aiming test",)"
         R"( "units": {"length": "mm", "wavelength": "um"}, "environment": {"medium": ")" +
         o.environment + R"("},)" + object + R"( "wavelengths": )" + o.wavelengths +
         R"(, "aperture": )" + o.aperture + R"(, "fields": )" + o.fields +
         R"(, "root": {"type": "assembly", "name": "system", "children": [)" + list +
         R"(]}, "paths": [{"name": "main", "events": "auto"}]})";
}

const rtt::material::MaterialLibrary& schott() {
  struct Schott {  // MaterialLibrary is not movable
    rtt::material::MaterialLibrary lib;
    Schott() { lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf"); }
  };
  static const Schott s;
  return s.lib;
}

CompiledSystem build(const std::vector<std::string>& children, const Options& o = {}) {
  return rtt::compile::compile(rtt::io::parse_system(system_json(children, o)), schott());
}

const std::string kCylinder = R"({"type": "cylinder", "radius": 51.68})";
const std::string kToroid = R"({"type": "toroid", "radius": 51.68, "radius_of_rotation": 77.52})";

/// Plano-convex cylinder lens (ADR 0042, acceptance 1): stop r = 10 at z = 0, curved side
/// R = 51.68 in the y section at z = 5, plane rear side at z = 9, n = 1.5168, VACUUM; detector at
/// z = 110. f_y = R/(n - 1) = 100 mm, x section afocal.
CompiledSystem cylinder_system(const Options& o = {}) {
  return build({stop(0.0, 10.0), lens(5.0, kCylinder, 4.0), detector(110.0)}, o);
}

/// Toroid in front of the stop (ADR 0042, acceptance 6 and 8; the system of
/// test_paraxial_section.cpp): toroid lens at z = 5 (t = 4), stop r = 4 at z = 20, detector at
/// z = 110; the entrance pupil differs per section in position and size.
CompiledSystem toroid_system(Options o = {}) {
  if (o.aperture == Options{}.aperture) o.aperture = R"({"type": "stop_size"})";
  return build({lens(5.0, kToroid, 4.0), stop(20.0, 4.0), detector(110.0)}, o);
}

// ------------------------------------------------------------------- pupil per section -----

/// Paraxial entrance pupil of one section: global z, radius and the stop height of a paraxial
/// ray through the EP point at unit height (ADR 0042, point 3).
struct SectionPupil {
  double z = 0.0;
  double r = 0.0;
  double stop_scale = 0.0;
};

SectionPupil pupil(const CompiledSystem& cs, Section s, std::uint16_t wl = 0) {
  const FirstOrder fo = rtt::paraxial::first_order_section(cs, kMain, wl, s);
  REQUIRE(fo.entrance_pupil.has_value());
  const rtt::paraxial::Pupil& ep = fo.entrance_pupil.value();
  REQUIRE(ep.z.has_value());
  REQUIRE(ep.diameter.has_value());
  SectionPupil p;
  p.z = ep.z.value();
  p.r = ep.diameter.value() / 2.0;
  const std::size_t stop = rtt::compile::stop_event(cs, kMain).value();
  p.stop_scale = rtt::paraxial::trace_ray_section(cs, kMain, wl, p.z, 1.0, 0.0, s)[stop].y;
  return p;
}

/// Index of the event at surface `id` on the main path.
std::size_t event_of(const CompiledSystem& cs, const std::string& id) {
  const auto& events = cs.path(kMain).events;
  for (std::size_t e = 0; e < events.size(); ++e) {
    if (cs.surfaces()[events[e].surface].id.str() == id) return e;
  }
  FAIL("no event at " << id);
  return 0;
}

/// x, y of the ray `ray` at global z (straight line through its position).
std::array<double, 2> at_z(const rtt::trace::RayState& ray, double z) {
  const double t = (z - ray.pos.z()) / ray.dir.z();
  return {ray.pos.x() + t * ray.dir.x(), ray.pos.y() + t * ray.dir.y()};
}

PointListPupil points(const std::vector<PupilPoint>& p) {
  return PointListPupil{p};
}

// ------------------------------------------------ independent trace of acceptance 1 -----

/// Refraction d' = mu d + (mu cos_i - cos_t) N in the y-z plane (de Greve, Eq. (22), (23),
/// (28)), N towards the incident medium (N.d < 0), mu = n / n'.
std::array<double, 2> refract(const std::array<double, 2>& d,
                              const std::array<double, 2>& n_hat,
                              double mu) {
  const double cos_i = -(n_hat[0] * d[0] + n_hat[1] * d[1]);
  const double sin2_t = mu * mu * (1.0 - cos_i * cos_i);
  const double k = mu * cos_i - std::sqrt(1.0 - sin2_t);
  return {mu * d[0] + k * n_hat[0], mu * d[1] + k * n_hat[1]};
}

/// Axis crossing (global z) of the ray parallel to the axis at height h in the y-z plane through
/// the plano-convex cylinder lens of cylinder_system(): the y section is a circle of radius R
/// around (y, z) = (0, 5 + R), then the plane z = 9 into vacuum.
double axis_crossing_by_hand(double h) {
  constexpr double kR = 51.68;
  const double zc = 5.0 + kR;
  const double z1 = zc - std::sqrt(kR * kR - h * h);
  const std::array<double, 2> n1{h / kR, (z1 - zc) / kR};  // towards -z, N.d < 0
  const std::array<double, 2> d1 = refract({0.0, 1.0}, n1, 1.0 / kN);
  const double y2 = h + d1[0] / d1[1] * (9.0 - z1);
  const std::array<double, 2> d2 = refract(d1, {0.0, -1.0}, kN);
  return 9.0 - y2 * d2[1] / d2[0];
}

/// Aimed rays (one field, one wavelength) traced to the end of the path with recording.
struct Traced {
  RayBatch start;
  RayBatch end;
  RayPaths paths;
};

Traced trace(const CompiledSystem& cs,
             std::uint16_t field,
             const rtt::trace::PupilSampling& sampling,
             Aiming aiming = Aiming::Real) {
  Traced t;
  const std::uint16_t fields[] = {field};
  t.start = rtt::trace::make_rays(cs, kMain, fields, 0, sampling, aiming);
  t.end = t.start;
  static_cast<void>(rtt::trace::SequentialTracer().trace(cs, kMain, t.end, t.paths));
  return t;
}

}  // namespace

// ------------------------------------------------------------------- acceptance 1 (real) -----

TEST_CASE("cylinder lens, real rays: x-z ray stays parallel, y rays cross at bfl + O(y^2)",
          "[trace][anamorphic][aiming]") {
  // ADR 0042, acceptance 1 (real rays). Paraxially f_y = R/(n - 1) and bfl = f_y - t/n (the
  // paraxial part is in test_paraxial_section.cpp); the rear focal point lies at
  // z_F = 9 + R/(n - 1) - 4/n. The rays are aimed for real with the aiming per section (the
  // EP is the stop at z = 0 in both sections, r_ep = 5): pupil (0.6, 0) is the height x = 3,
  // pupil (0, 0.1) and (0, 0.2) the heights y = 0.5 and 1.
  const CompiledSystem cs = cylinder_system();
  const Traced t = trace(cs, 0, points({{0.6, 0.0}, {0.0, 0.1}, {0.0, 0.2}}));
  for (std::size_t i = 0; i < t.end.size(); ++i) {
    INFO("ray " << i);
    REQUIRE(t.start.status()[i] == RayStatus::Alive);
    REQUIRE(t.end.status()[i] == RayStatus::Alive);
  }
  // The ray in the x-z plane: the cylinder is flat in x, so the ray is not deviated.
  CHECK(std::abs(t.end.pos_x()[0] - 3.0) <= 1e-12);
  CHECK(std::abs(t.end.pos_y()[0]) <= 1e-12);
  CHECK(std::abs(t.end.dir_x()[0]) <= 1e-15);
  CHECK(std::abs(t.end.dir_y()[0]) <= 1e-15);
  // The rays in the y-z plane: axis crossing z_c(h), against the hand trace and against
  // z_F + O(h^2) (third-order longitudinal aberration, (z_c - z_F)/h^2 about constant).
  const double z_focus = 9.0 + 51.68 / (kN - 1.0) - 4.0 / kN;
  std::array<double, 2> per_h2{};
  const std::array<double, 2> heights{0.5, 1.0};
  for (std::size_t k = 0; k < 2; ++k) {
    const std::size_t i = k + 1;
    const double h = heights[k];
    INFO("h = " << h);
    // Aimed through the stop at y = h, parallel to the axis (object at infinity, field 0).
    CHECK(std::abs(t.start.pos_y()[i] - h) <= 1e-12);
    CHECK(std::abs(t.end.pos_x()[i]) <= 1e-12);
    CHECK(std::abs(t.end.dir_x()[i]) <= 1e-15);
    const double z_c = t.end.pos_z()[i] - t.end.pos_y()[i] * t.end.dir_z()[i] / t.end.dir_y()[i];
    CHECK(std::abs(z_c - axis_crossing_by_hand(h)) <= 1e-9);
    per_h2[k] = (z_c - z_focus) / (h * h);
  }
  // Content guard: the longitudinal aberration is not zero (about -0.0109 mm per mm^2).
  CHECK(per_h2[1] < -1e-3);
  CHECK(std::abs(per_h2[0] - per_h2[1]) <= 1e-3 * std::abs(per_h2[1]));
}

// ------------------------------------------------------------------- acceptance 8 -----

TEST_CASE("toroid in front of the stop: chief rays hit the stop centre, marginal rays the target",
          "[trace][anamorphic][aiming]") {
  // ADR 0042, acceptance 8. stop_size: the target of pupil point (px, py) is
  // (px r_ep,x stop_scale_x, py r_ep,y stop_scale_y) = r_stop (px, py) (point 3, table); the EP
  // is the virtual, upright image of the stop through the lens, so stop_scale_s > 0.
  Options o;
  o.fields = R"({"points": [{}, {"x": 2.0}, {"y": 2.0}, {"x": -1.5, "y": 1.0}]})";
  const CompiledSystem cs = toroid_system(o);
  const SectionPupil px = pupil(cs, Section::x);
  const SectionPupil py = pupil(cs, Section::y);
  // Content guard: the pupils differ per section (acceptance 6).
  REQUIRE(std::abs(px.z - py.z) > 0.1);
  REQUIRE(std::abs(px.r - py.r) > 0.01);
  REQUIRE(px.stop_scale > 0.0);
  REQUIRE(py.stop_scale > 0.0);
  CHECK(std::abs(px.r * px.stop_scale - 4.0) <= 1e-12);
  CHECK(std::abs(py.r * py.stop_scale - 4.0) <= 1e-12);

  const std::vector<PupilPoint> pupil_points = {{0.0, 0.0},  {1.0, 0.0},  {0.0, 1.0}, {-1.0, 0.0},
                                                {0.0, -1.0}, {-0.6, 0.8}, {0.5, 0.5}};
  const std::size_t stop_slot = event_of(cs, "STO") + 1;
  double paraxial_miss = 0.0;
  for (std::uint16_t field = 0; field < 4; ++field) {
    const Traced t = trace(cs, field, points(pupil_points));
    const Traced p = trace(cs, field, points(pupil_points), Aiming::Paraxial);
    for (std::size_t i = 0; i < pupil_points.size(); ++i) {
      INFO("field " << field << ", pupil (" << pupil_points[i].px << ", " << pupil_points[i].py
                    << ")");
      REQUIRE(t.start.status()[i] == RayStatus::Alive);
      REQUIRE(t.end.status()[i] == RayStatus::Alive);
      REQUIRE(t.paths.count[i] == t.paths.slots);
      const Vec3 hit = t.paths.position_at(i, stop_slot);
      CHECK(std::abs(hit.z() - 20.0) <= 1e-12);
      const double tx = 4.0 * pupil_points[i].px;
      const double ty = 4.0 * pupil_points[i].py;
      CHECK(std::hypot(hit.x() - tx, hit.y() - ty) <= rtt::trace::kAimTolerance);
      // A paraxial start ray may be vignetted at the stop (outside r = 4); its hit there is
      // recorded all the same (count = stop event + 2). Any other loss must not feed the guard.
      const bool passed = p.end.status()[i] == RayStatus::Alive;
      const bool vignetted_at_stop = p.end.status()[i] == RayStatus::Vignetted &&
                                     p.paths.lost_at[i] == static_cast<std::int32_t>(stop_slot - 1);
      REQUIRE((passed || vignetted_at_stop));
      REQUIRE(p.paths.count[i] > stop_slot);
      const Vec3 q = p.paths.position_at(i, stop_slot);
      paraxial_miss = std::max(paraxial_miss, std::hypot(q.x() - tx, q.y() - ty));
    }
  }
  // Content guard: the paraxial start rays miss the targets by more than the tolerance, so the
  // real aiming has work to do.
  CHECK(paraxial_miss > 1e-6);
}

TEST_CASE("paraxial start ray per section: through (z_ep,x; a) and (z_ep,y; b)",
          "[trace][anamorphic][aiming]") {
  // ADR 0042, point 3, "Paraxialer Startstrahl", a = s_x r_ep,x, b = s_y r_ep,y.
  constexpr double kDeg = std::numbers::pi / 180.0;
  SECTION("object at infinity: direction from the field") {
    Options o;
    o.fields = R"({"points": [{"x": 1.0, "y": 2.0}]})";
    const CompiledSystem cs = toroid_system(o);
    const SectionPupil px = pupil(cs, Section::x);
    const SectionPupil py = pupil(cs, Section::y);
    const Vec3 d = Vec3(std::tan(1.0 * kDeg), std::tan(2.0 * kDeg), 1.0).normalized();
    for (const PupilPoint& p : std::vector<PupilPoint>{{0.0, 0.0}, {0.3, -0.7}, {-1.0, 0.0}}) {
      INFO("pupil (" << p.px << ", " << p.py << ")");
      const rtt::trace::AimedRay r =
          rtt::trace::aim_ray(cs, kMain, std::uint16_t{0}, 0, p.px, p.py, Aiming::Paraxial);
      REQUIRE(r.ray.status == RayStatus::Alive);
      CHECK((r.ray.dir - d).norm() <= 1e-15);
      CHECK(std::abs(at_z(r.ray, px.z)[0] - p.px * px.r) <= 1e-12);
      CHECK(std::abs(at_z(r.ray, py.z)[1] - p.py * py.r) <= 1e-12);
      // The start lies in front of the toroid (z = 5).
      CHECK(r.ray.pos.z() < 5.0);
    }
  }
  SECTION("finite object: slopes (a - h_x)/(z_ep,x - z_obj), (b - h_y)/(z_ep,y - z_obj)") {
    Options o;
    o.object = R"({"at_infinity": false, "distance": 200.0})";
    o.fields = R"({"type": "object_height", "points": [{"x": 1.0, "y": 2.0}]})";
    const CompiledSystem cs = toroid_system(o);
    const SectionPupil px = pupil(cs, Section::x);
    const SectionPupil py = pupil(cs, Section::y);
    for (const PupilPoint& p : std::vector<PupilPoint>{{0.0, 0.0}, {0.3, -0.7}, {0.0, 1.0}}) {
      INFO("pupil (" << p.px << ", " << p.py << ")");
      const rtt::trace::AimedRay r =
          rtt::trace::aim_ray(cs, kMain, std::uint16_t{0}, 0, p.px, p.py, Aiming::Paraxial);
      REQUIRE(r.ray.status == RayStatus::Alive);
      CHECK(r.ray.pos == Vec3(1.0, 2.0, -200.0));
      const double mx = (p.px * px.r - 1.0) / (px.z + 200.0);
      const double my = (p.py * py.r - 2.0) / (py.z + 200.0);
      CHECK((r.ray.dir - Vec3(mx, my, 1.0).normalized()).norm() <= 1e-15);
    }
  }
  SECTION("field angle with a finite object: object point per section from z_ep,s") {
    Options o;
    o.object = R"({"at_infinity": false, "distance": 200.0})";
    o.fields = R"({"points": [{"x": 1.0, "y": 2.0}]})";
    const CompiledSystem cs = toroid_system(o);
    const SectionPupil px = pupil(cs, Section::x);
    const SectionPupil py = pupil(cs, Section::y);
    const rtt::trace::AimedRay r =
        rtt::trace::aim_ray(cs, kMain, std::uint16_t{0}, 0, 0.0, 0.0, Aiming::Paraxial);
    REQUIRE(r.ray.status == RayStatus::Alive);
    CHECK(std::abs(r.ray.pos.x() - (-200.0 - px.z) * std::tan(1.0 * kDeg)) <= 1e-12);
    CHECK(std::abs(r.ray.pos.y() - (-200.0 - py.z) * std::tan(2.0 * kDeg)) <= 1e-12);
    CHECK(r.ray.pos.z() == -200.0);
    // The chief ray passes the EP centre of each section.
    CHECK(std::abs(at_z(r.ray, px.z)[0]) <= 1e-12);
    CHECK(std::abs(at_z(r.ray, py.z)[1]) <= 1e-12);
  }
  SECTION("paraxial image height: unit image height per section") {
    Options o;
    o.fields = R"({"type": "paraxial_image_height", "points": [{"x": 1.0, "y": 2.0}]})";
    const CompiledSystem cs = toroid_system(o);
    // Unit image height per section: the chief ray of unit slope through the EP centre, at the
    // paraxial image plane of that section (ADR 0042, point 3, "Paraxiale Bildhöhe"). It uses
    // first_order_section and trace_ray_section, which test_paraxial_section.cpp verifies
    // independently against a 2x2 matrix computation.
    const auto unit = [&](Section s) {
      const SectionPupil ep = pupil(cs, s);
      const FirstOrder fo = rtt::paraxial::first_order_section(cs, kMain, 0, s);
      const auto ray = rtt::paraxial::trace_ray_section(cs, kMain, 0, ep.z, 0.0, 1.0, s);
      return ray.back().y + (fo.image_z.value() - ray.back().z) * ray.back().u;
    };
    const double ux = unit(Section::x);
    const double uy = unit(Section::y);
    REQUIRE(std::abs(ux - uy) > 1e-3 * std::abs(uy));  // content guard: they differ
    const rtt::trace::AimedRay r =
        rtt::trace::aim_ray(cs, kMain, std::uint16_t{0}, 0, 0.0, 0.0, Aiming::Paraxial);
    REQUIRE(r.ray.status == RayStatus::Alive);
    CHECK(std::abs(r.ray.dir.x() / r.ray.dir.z() - 1.0 / ux) <= 1e-14);
    CHECK(std::abs(r.ray.dir.y() / r.ray.dir.z() - 2.0 / uy) <= 1e-14);
  }
}

TEST_CASE("object-space telecentric anamorphic path: slopes in both sections",
          "[trace][anamorphic][aiming]") {
  // Ideal lens f = 50 at z = 0 with the stop in its rear focal plane (z = 50): the EP is at
  // infinity in both sections. An ideal cylinder lens behind the stop makes the path anamorphic
  // without touching the EP. object_na 0.05 read paraxially: r_ep,s = u_m = 0.05 in both
  // sections (point 3, table); pupil coordinates are slopes (#96).
  Options o;
  o.object = R"({"at_infinity": false, "distance": 100.0})";
  o.aperture = R"({"type": "object_na", "value": 0.05})";
  o.fields = R"({"type": "object_height", "points": [{"x": 1.0, "y": -2.0}]})";
  const CompiledSystem cs =
      build({thin(0.0, R"({"type": "ideal_lens", "focal_length": 50.0})", "front"), stop(50.0, 5.0),
             thin(60.0, cylinder_lens(40.0), "rear"), detector(150.0)},
            o);
  REQUIRE(rtt::paraxial::is_anamorphic(cs, kMain));
  const rtt::trace::AimedRay chief =
      rtt::trace::aim_ray(cs, kMain, std::uint16_t{0}, 0, 0.0, 0.0, Aiming::Paraxial);
  CHECK(chief.ray.pos == Vec3(1.0, -2.0, -100.0));
  CHECK(chief.ray.dir == Vec3::UnitZ());
  const rtt::trace::AimedRay marginal =
      rtt::trace::aim_ray(cs, kMain, std::uint16_t{0}, 0, 0.6, 0.8, Aiming::Paraxial);
  CHECK((marginal.ray.dir - Vec3(0.03, 0.04, 1.0).normalized()).norm() <= 1e-15);
  // Real aiming: the chief ray hits the stop centre (the ideal lens is exact).
  const rtt::trace::AimedRay real =
      rtt::trace::aim_ray(cs, kMain, std::uint16_t{0}, 0, 0.0, 0.0, Aiming::Real);
  REQUIRE(real.ray.status == RayStatus::Alive);
  CHECK(std::abs(at_z(real.ray, 0.0)[0] - 1.0) <= 1e-12);
}

TEST_CASE("aim reference of an anamorphic system aims like the system itself",
          "[trace][anamorphic][aiming]") {
  // ADR 0042, point 3, "AimReference": the reference keeps the data per section.
  Options o;
  o.fields = R"({"type": "paraxial_image_height", "points": [{}, {"x": 1.0}, {"y": -2.0}]})";
  const CompiledSystem cs = toroid_system(o);
  const rtt::trace::AimReference ref = rtt::trace::aim_reference(cs, kMain, 0);
  const std::uint16_t fields[] = {0, 1, 2};
  const RayBatch a = rtt::trace::make_rays(cs, kMain, fields, 0, rtt::trace::HexapolarPupil{3});
  const RayBatch b =
      rtt::trace::make_rays(cs, kMain, fields, 0, rtt::trace::HexapolarPupil{3}, Aiming::Real, ref);
  REQUIRE(a.size() == b.size());
  double spread = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    INFO("ray " << i);
    CHECK(a.status()[i] == RayStatus::Alive);
    CHECK(b.status()[i] == RayStatus::Alive);
    CHECK(a.pos_x()[i] == b.pos_x()[i]);
    CHECK(a.pos_y()[i] == b.pos_y()[i]);
    CHECK(a.pos_z()[i] == b.pos_z()[i]);
    CHECK(a.dir_x()[i] == b.dir_x()[i]);
    CHECK(a.dir_y()[i] == b.dir_y()[i]);
    CHECK(a.dir_z()[i] == b.dir_z()[i]);
    spread = std::max(spread, std::abs(a.dir_x()[i]) + std::abs(a.dir_y()[i]));
  }
  CHECK(spread > 1e-3);  // content guard: the rays differ
  CHECK(rtt::trace::stop_radius(ref) == 4.0);
  CHECK(rtt::trace::paraxial_stop_radius(ref) == 4.0);
}

TEST_CASE("stop radius of an anamorphic path", "[trace][anamorphic][aiming]") {
  // ADR 0042, point 3, "stop_radius und paraxial_stop_radius".
  SECTION("stop_size: r_stop directly") {
    const CompiledSystem cs = toroid_system();
    CHECK(rtt::trace::stop_radius(cs, kMain, 0) == 4.0);
    CHECK(rtt::trace::paraxial_stop_radius(cs, kMain, 0) == 4.0);
  }
  SECTION("EPD with equal sections at the stop: the value of section y") {
    // The cylinder lens lies behind the stop: r_ep,s stop_scale_s = 5 in both sections.
    const CompiledSystem cs = cylinder_system();
    CHECK(rtt::trace::stop_radius(cs, kMain, 0) == 5.0);
    CHECK(rtt::trace::paraxial_stop_radius(cs, kMain, 0) == 5.0);
  }
  SECTION("EPD with an anamorphic part in front of the stop: the sections differ, error") {
    Options o;
    o.aperture = R"({"type": "epd", "value": 4.0})";
    const CompiledSystem cs = toroid_system(o);
    const SectionPupil px = pupil(cs, Section::x);
    const SectionPupil py = pupil(cs, Section::y);
    REQUIRE(std::abs(px.stop_scale - py.stop_scale) > 1e-3);  // content guard
    CHECK_THROWS_WITH(rtt::trace::stop_radius(cs, kMain, 0),
                      ContainsSubstring("differs between the sections"));
    CHECK_THROWS_WITH(rtt::trace::paraxial_stop_radius(cs, kMain, 0),
                      ContainsSubstring("differs between the sections"));
    // Aiming itself works: the target is an ellipse.
    const rtt::trace::AimedRay r =
        rtt::trace::aim_ray(cs, kMain, std::uint16_t{0}, 0, 1.0, 0.0, Aiming::Real);
    CHECK(r.ray.status == RayStatus::Alive);
  }
}

// ------------------------------------------------------------------- acceptance 5 -----

TEST_CASE("start plane: strongly anamorphic biconic without aperture, bundle beyond the incircle",
          "[trace][anamorphic][aiming]") {
  // ADR 0042, acceptance 5 (start plane). Concave biconic front R_x = -200, R_y = -10
  // (R_x/R_y = 20) at z = 0 without aperture, plane rear side at z = 2, n = 1.5; stop r = 2 at
  // z = 40, stop_size. The incircle of the domain has radius |R_y| = 10; the field of 20 degree
  // in x puts the bundle at |x| of about 14 on the biconic, inside the ellipse of the domain but
  // outside the incircle. The surface curves back against the light (sag < 0), and the EP lies
  // behind it in both sections (z_ep,x != z_ep,y): the start plane must lie before every hit.
  Options o;
  o.aperture = R"({"type": "stop_size"})";
  o.fields = R"({"points": [{}, {"x": 20.0}, {"x": -20.0, "y": 3.0}]})";
  const CompiledSystem cs =
      build({lens(0.0, R"({"type": "biconic", "radius_x": -200.0, "radius_y": -10.0})", 2.0,
                  "CONST:1.5", ""),
             stop(40.0, 2.0), detector(80.0)},
            o);
  const SectionPupil px = pupil(cs, Section::x);
  const SectionPupil py = pupil(cs, Section::y);
  REQUIRE(px.z > 2.0);
  REQUIRE(py.z > 2.0);
  REQUIRE(std::abs(px.z - py.z) > 1.0);
  const std::size_t front = event_of(cs, "L1.S1") + 1;
  double beyond = 0.0;
  for (std::uint16_t field = 0; field < 3; ++field) {
    const Traced t = trace(cs, field, rtt::trace::HexapolarPupil{4});
    for (std::size_t i = 0; i < t.end.size(); ++i) {
      INFO("field " << field << ", ray " << i);
      REQUIRE(t.start.status()[i] == RayStatus::Alive);
      REQUIRE(t.end.status()[i] == RayStatus::Alive);
      const Vec3 hit = t.paths.position_at(i, front);
      CHECK(t.start.pos_z()[i] < hit.z());
      CHECK(std::isfinite(t.start.pos_z()[i]));
      beyond = std::max(beyond, std::hypot(hit.x(), hit.y()));
    }
  }
  CHECK(beyond > 10.0);  // content guard: hits outside the incircle
}

// ------------------------------------------------------------------- acceptance 7 -----

TEST_CASE("aiming errors of anamorphic paths (ADR 0042, point 3)", "[trace][anamorphic][aiming]") {
  const std::uint16_t field0[] = {0};
  const auto rays = [&](const CompiledSystem& cs, std::uint16_t wl = 0) {
    return rtt::trace::make_rays(cs, kMain, field0, wl, rtt::trace::HexapolarPupil{2});
  };
  SECTION("entrance pupil at infinity in one section only") {
    // Ideal cylinder lens f = 50 (power in y) at z = 0, stop in its rear focal plane: EP at
    // infinity in y, the stop itself in x.
    const CompiledSystem cs =
        build({thin(0.0, cylinder_lens(50.0)), stop(50.0, 5.0), detector(100.0)});
    REQUIRE_FALSE(rtt::paraxial::first_order_section(cs, kMain, 0, Section::y)
                      .entrance_pupil.value()
                      .z.has_value());
    REQUIRE(rtt::paraxial::first_order_section(cs, kMain, 0, Section::x)
                .entrance_pupil.value()
                .z.has_value());
    CHECK_THROWS_AS(rays(cs), std::invalid_argument);
    CHECK_THROWS_WITH(rays(cs), ContainsSubstring("one section"));
    CHECK_THROWS_WITH(rtt::trace::aim_reference(cs, kMain, 0), ContainsSubstring("one section"));
  }
  SECTION("entrance pupil at infinity in one section only at the reference wavelength") {
    // N-BK7 cylinder lens, the stop in the rear focal plane of the y section at the reference
    // wavelength 0.5876 um; at 0.4861 um the EP is finite in both sections. Paraxial image
    // heights convert at the reference wavelength (#31), so the context needs it.
    Options o;
    o.wavelengths = R"([{"um": 0.4861}, {"um": 0.5876, "reference": true}])";
    o.fields = R"({"type": "paraxial_image_height", "points": [{}]})";
    const auto with_stop = [&](double z) {
      return build({lens(5.0, kCylinder, 4.0, "SCHOTT:N-BK7"), stop(z, 5.0), detector(300.0)}, o);
    };
    const double z_focus =
        rtt::paraxial::first_order_section(with_stop(200.0), kMain, 1, Section::y)
            .rear_focal_z.value();
    const CompiledSystem cs = with_stop(z_focus);
    REQUIRE_FALSE(rtt::paraxial::first_order_section(cs, kMain, 1, Section::y)
                      .entrance_pupil.value()
                      .z.has_value());
    REQUIRE(rtt::paraxial::first_order_section(cs, kMain, 0, Section::y)
                .entrance_pupil.value()
                .z.has_value());
    REQUIRE(rtt::paraxial::first_order_section(cs, kMain, 0, Section::x)
                .entrance_pupil.value()
                .z.has_value());
    CHECK_THROWS_WITH(
        rays(cs, 0), ContainsSubstring("one section") && ContainsSubstring("reference wavelength"));
  }
  SECTION("pupil labels differ between the sections") {
    // Finite object at z = -20; ideal cylinder lens f = 10 (power in y) at z = 0, stop at
    // z = 12: z_ep,y = 12/(1 - 12/10) = -60 lies beyond the object (label_sign -1), z_ep,x = 12
    // does not (+1). ADR 0007, addendum #35 B.
    Options o;
    o.object = R"({"at_infinity": false, "distance": 20.0})";
    o.aperture = R"({"type": "epd", "value": 2.0})";
    o.fields = R"({"type": "object_height", "points": [{}]})";
    const CompiledSystem cs =
        build({thin(0.0, cylinder_lens(10.0)), stop(12.0, 1.0), detector(40.0)}, o);
    REQUIRE(std::abs(pupil(cs, Section::y).z + 60.0) <= 1e-9);
    REQUIRE(pupil(cs, Section::x).z == 12.0);
    CHECK_THROWS_WITH(rays(cs), ContainsSubstring("label"));
  }
  SECTION("object_na read real") {
    Options o;
    o.object = R"({"at_infinity": false, "distance": 100.0})";
    o.aperture = R"({"type": "object_na", "value": 0.05, "reading": "real"})";
    o.fields = R"({"type": "object_height", "points": [{}]})";
    const CompiledSystem cs = cylinder_system(o);
    CHECK_THROWS_WITH(rays(cs), ContainsSubstring("not for anamorphic systems (ADR 0042)"));
    // Read paraxially: r_ep,s = u_m |z_ep,s - z_obj| = 0.05 * 100 in both sections (EP = stop).
    o.aperture = R"({"type": "object_na", "value": 0.05})";
    const CompiledSystem paraxial = cylinder_system(o);
    const rtt::trace::AimedRay r =
        rtt::trace::aim_ray(paraxial, kMain, std::uint16_t{0}, 0, 0.0, 1.0, Aiming::Paraxial);
    CHECK(std::abs(at_z(r.ray, 0.0)[1] - 5.0) <= 1e-12);
  }
  SECTION("image-space F-number") {
    Options o;
    o.aperture = R"({"type": "image_fnumber", "value": 10.0})";
    const CompiledSystem cs = cylinder_system(o);
    CHECK_THROWS_AS(rays(cs), std::invalid_argument);
    CHECK_THROWS_WITH(rays(cs), ContainsSubstring("F-number"));
  }
  SECTION("paraxial image height with an afocal section") {
    // The cylinder lens is afocal in x with the object at infinity; the error holds also for a
    // field point with x = 0, because the field type holds for the whole system.
    Options o;
    o.fields = R"({"type": "paraxial_image_height", "points": [{"y": 1.0}]})";
    const CompiledSystem cs = cylinder_system(o);
    const std::string message =
        "paraxial image height needs a finite paraxial image in both sections";
    CHECK_THROWS_AS(rays(cs), std::invalid_argument);
    CHECK_THROWS_WITH(rays(cs), ContainsSubstring(message));
    // The reference keeps the error for the field conversion (#255).
    const rtt::trace::AimReference ref = rtt::trace::aim_reference(cs, kMain, 0);
    CHECK_THROWS_WITH(rtt::trace::make_rays(cs, kMain, field0, 0, rtt::trace::HexapolarPupil{2},
                                            Aiming::Real, ref),
                      ContainsSubstring(message));
  }
}
