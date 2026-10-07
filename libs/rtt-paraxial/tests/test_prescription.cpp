#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/paraxial/prescription.hpp"
#include "rtt/paraxial/seidel.hpp"

// Sources of the reference formulas (docs/quellen.md):
// - y-nu trace: Greivenkamp, OPTI-201/202 lecture notes, Sec. 9, p. 9-2: n'u' = n u - y phi,
//   phi = (n' - n) c, y' = y + u t.
// - Angle of incidence and Lagrange invariant: Sasian, OPTI 517 L4, p. 24: A = n i = n u +
//   n y c, H = n u_bar y - n u y_bar; Greivenkamp, OPTI-502 Sec. 9, p. 9-41 (invariant on
//   refraction and transfer; H = -n u y_bar in an object or image plane).
// - NA = n sin U ~ n u (p. 9-34), F/# = f / D_EP = 1 / (2 n' |u'|) for an object at infinity
//   (p. 9-35), working F/#_W = 1 / (2 n |u|) (p. 9-36) and F/#_W = (1 - m) F/# for a thin lens
//   with the stop at the lens (p. 9-37, 9-38).
// - Thick lens: lensmaker formula as in test_paraxial.cpp (Hecht, Optics, Ch. 6).
//
// Tolerances (a priori): the reference values are a few floating-point operations on the
// inputs and the y-nu trace adds about one rounding per operation and event, so the relative
// error of a ray value is a small multiple of 1e-16 for these few events; 1e-12 relative (and
// 1e-12 mm or rad absolute for values that are 0) leaves a margin of more than 1e3. For the
// Lagrange invariant after an event the rounding errors of y, u, y_bar, u_bar enter with the
// size of the two products, so it is compared relative to |n u_bar y| + |n u y_bar| there.

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::paraxial::ParaxialError;
using rtt::paraxial::Prescription;
using rtt::paraxial::prescription;
using rtt::paraxial::PrescriptionRay;
using rtt::paraxial::PrescriptionSurface;

namespace {

constexpr double kRel = 1e-12;
constexpr double kAbs = 1e-12;  // mm or rad, for values that are 0

// ------------------------------------------------------------------ builders -----

/// Reference values are formulas for a surrounding index of exactly 1 (VACUUM; AIR is Ciddor
/// air since #25).
System base_system() {
  System s;
  s.name = "prescription test";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(20.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 5.0, 1.0}}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  return s;
}

/// Surface at local z with radius R in mm (std::nullopt = plane).
Surface surface(const std::string& id, double z, std::optional<double> radius) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z);
  if (radius) s.shape.base = rtt::model::Conic{Param(*radius), Param(0.0)};
  return s;
}

/// Lens with vertex at global z, index n ("CONST:1.500000"), radii R1, R2, centre thickness d.
Element lens(const std::string& name,
             double z,
             double n,
             std::optional<double> r1,
             std::optional<double> r2,
             double d) {
  return Element{name,
                 ElementKind::Lens,
                 Pose::along_z(z),
                 "CONST:" + std::to_string(n),
                 {surface(name + ".S1", 0.0, r1), surface(name + ".S2", d, r2)}};
}

Element stop(const std::string& name, double z, double radius) {
  Surface s = surface(name + ".S", 0.0, std::nullopt);
  s.aperture = rtt::model::CircularAperture{radius, 0.0};
  return Element{name, ElementKind::Stop, Pose::along_z(z), std::nullopt, {s}};
}

Element mirror(const std::string& name, double z, std::optional<double> radius) {
  Surface s = surface(name + ".S", 0.0, radius);
  s.interaction = rtt::model::IdealMirror{};
  return Element{name, ElementKind::Mirror, Pose::along_z(z), std::nullopt, {s}};
}

Element detector(const std::string& name, double z) {
  return Element{name,
                 ElementKind::Detector,
                 Pose::along_z(z),
                 std::nullopt,
                 {surface(name + ".S", 0.0, std::nullopt)}};
}

void add(System& s, Element e) {
  s.root.children.push_back({std::move(e)});
}

CompiledSystem compile(const System& s) {
  const MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

Prescription prescription_of(const System& s) {
  return prescription(compile(s), PathId{0}, 0);
}

System load_in_vacuum(const std::string& name) {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + name);
  s.environment.medium = "VACUUM";
  return s;
}

double tan_deg(double deg) {
  return std::tan(deg * std::numbers::pi / 180.0);
}

/// Compares with relative tolerance, or absolute for an expected 0.
void require_near(double actual, double expected) {
  if (expected == 0.0) {
    REQUIRE_THAT(actual, WithinAbs(0.0, kAbs));
  } else {
    REQUIRE_THAT(actual, WithinRel(expected, kRel));
  }
}

void require_near(const std::optional<double>& actual, double expected) {
  REQUIRE(actual.has_value());
  require_near(*actual, expected);
}

/// Expected ray at one event: height, slope after, angle of incidence.
struct Expected {
  double y, u, i;
};

void require_ray(const std::optional<PrescriptionRay>& ray, const Expected& e) {
  REQUIRE(ray.has_value());
  INFO("expected y = " << e.y << ", u = " << e.u << ", i = " << e.i);
  require_near(ray->y, e.y);
  require_near(ray->u, e.u);
  require_near(ray->i, e.i);
}

/// The Lagrange invariant after every event equals the object-space value (see tolerances).
void require_lagrange_constant(const Prescription& p) {
  REQUIRE(p.lagrange_invariant.has_value());
  const double h = *p.lagrange_invariant;
  for (const PrescriptionSurface& s : p.surfaces) {
    REQUIRE(s.marginal.has_value());
    REQUIRE(s.chief.has_value());
    REQUIRE(s.lagrange.has_value());
    const double scale =
        std::abs(s.n * s.chief->u * s.marginal->y) + std::abs(s.n * s.marginal->u * s.chief->y);
    INFO("surface " << s.surface << ": H = " << *s.lagrange << ", object space " << h);
    REQUIRE(std::abs(*s.lagrange - h) <= kRel * scale);
  }
}

}  // namespace

// ------------------------------------------------------------------- tests -----

TEST_CASE("thin lens at infinity: rays, angles of incidence, F/#, NA and H by hand",
          "[prescription]") {
  // Stop (r = 10) and a plano-convex lens of zero thickness at z = 0 (n = 1.5, R1 = 50:
  // f = R1 / (n - 1) = 100 mm, exact for d = 0), detector at z = 100. EPD 20, field 0 and 5 deg.
  System s = base_system();
  add(s, stop("STO", 0.0, 10.0));
  add(s, lens("L", 0.0, 1.5, 50.0, std::nullopt, 0.0));
  add(s, detector("IMG", 100.0));
  const Prescription p = prescription_of(s);
  const double t = tan_deg(5.0);

  REQUIRE(p.surfaces.size() == 4);  // stop, L.S1, L.S2, IMG
  const std::vector<double> z = {0.0, 0.0, 0.0, 100.0};
  const std::vector<double> n = {1.0, 1.5, 1.0, 1.0};
  for (std::size_t k = 0; k < 4; ++k) {
    REQUIRE(p.surfaces[k].z == z[k]);
    REQUIRE(p.surfaces[k].n == n[k]);
  }

  // Marginal ray y = 10, u = 0. S1: i = u + y c = 10 / 50, n'u' = -y (n' - n) c = -0.1.
  // S2 (plane): i = u = -0.1 / 1.5, n'u' = n u = -0.1. Detector: y = 10 - 100 * 0.1 = 0.
  require_ray(p.surfaces[0].marginal, {10.0, 0.0, 0.0});
  require_ray(p.surfaces[1].marginal, {10.0, -0.1 / 1.5, 0.2});
  require_ray(p.surfaces[2].marginal, {10.0, -0.1, -0.1 / 1.5});
  require_ray(p.surfaces[3].marginal, {0.0, -0.1, -0.1});
  // Chief ray through the stop centre with u = tan 5 deg: y = 0 at the lens, i = u there.
  require_ray(p.surfaces[0].chief, {0.0, t, t});
  require_ray(p.surfaces[1].chief, {0.0, t / 1.5, t});
  require_ray(p.surfaces[2].chief, {0.0, t, t / 1.5});
  require_ray(p.surfaces[3].chief, {100.0 * t, t, t});

  REQUIRE(p.marginal_start.has_value());
  REQUIRE(p.marginal_start->y == 10.0);
  REQUIRE(p.marginal_start->u == 0.0);
  REQUIRE(p.chief_start.has_value());
  require_near(p.chief_start->u, t);

  require_near(p.total_track, 100.0);
  REQUIRE_FALSE(p.object_distance.has_value());
  // F/# = f / D_EP = 100 / 20 = 1 / (2 n' |u'|) = 1 / (2 * 0.1) (Greivenkamp p. 9-35).
  require_near(p.paraxial_working_f_number, 5.0);
  require_near(p.paraxial_image_na, 0.1);
  // H = n (u_bar y - u y_bar) = tan 5 deg * 10 in the stop (y_bar = 0, Greivenkamp p. 9-41).
  require_near(p.lagrange_invariant, 10.0 * t);
  require_lagrange_constant(p);
  REQUIRE(p.first_order.efl.has_value());
  require_near(p.first_order.efl, 100.0);
}

TEST_CASE("thin lens at finite conjugates: F/#_W = (1 - m) F/#, H in object and image plane",
          "[prescription]") {
  // Same lens, object at z = -200 (2 f), object height 5 mm, detector in the image at z = 200.
  System s = base_system();
  s.object = {false, Param(200.0)};
  s.fields = {rtt::model::FieldType::ObjectHeight, {{0.0, 0.0, 1.0}, {0.0, 5.0, 1.0}}};
  add(s, stop("STO", 0.0, 10.0));
  add(s, lens("L", 0.0, 1.5, 50.0, std::nullopt, 0.0));
  add(s, detector("IMG", 200.0));
  const Prescription p = prescription_of(s);
  REQUIRE(p.surfaces.size() == 4);

  // Marginal ray from (-200, 0) to the pupil rim: u = 10 / 200 = 0.05; behind the lens
  // u' = u - y / f = 0.05 - 0.1 = -0.05, in the image y = 0.
  require_ray(p.surfaces[0].marginal, {10.0, 0.05, 0.05});
  require_ray(p.surfaces[2].marginal, {10.0, -0.05, -0.05 / 1.5});
  require_ray(p.surfaces[3].marginal, {0.0, -0.05, -0.05});
  // Chief ray from (-200, 5) through the pupil centre: u_bar = -5 / 200 = -0.025, image
  // height -5 = m * 5.
  require_ray(p.surfaces[3].chief, {-5.0, -0.025, -0.025});

  require_near(p.object_distance, 200.0);
  require_near(p.total_track, 200.0);
  REQUIRE(p.first_order.lateral_magnification.has_value());
  require_near(p.first_order.lateral_magnification, -1.0);
  // F/#_W = 1 / (2 n' |u'|) = 10 = (1 - m) f / D_EP with m = -1 (Greivenkamp p. 9-37/9-38,
  // exact here: thin lens with the stop at the lens).
  require_near(p.paraxial_working_f_number, 10.0);
  require_near(p.paraxial_working_f_number,
               (1.0 - *p.first_order.lateral_magnification) * *p.first_order.efl / 20.0);
  require_near(p.paraxial_image_na, 0.05);
  // In the object plane y = 0: H = -n u y_bar = -0.05 * 5 = -0.25; in the image plane
  // H = -n' u' y_bar' = -(-0.05)(-5) (Greivenkamp p. 9-41, 9-42).
  require_near(p.lagrange_invariant, -0.25);
  require_near(p.surfaces[3].lagrange, -0.25);
  require_lagrange_constant(p);
}

TEST_CASE("reference singlet: thick-lens rays, F/# = EFL / EPD, H constant", "[prescription]") {
  // tests/reference/m1/singlet_const.rtt.json in VACUUM: stop (r = 10) at z = 0, plano-convex
  // lens n = 1.5168, R1 = 51.68 at z = 5, d = 4, image at z = 106.363; EPD 20, max field 5 deg.
  const CompiledSystem cs = compile(load_in_vacuum("m1/singlet_const.rtt.json"));
  const Prescription p = prescription(cs, PathId{0}, 0);
  REQUIRE(p.surfaces.size() == 4);  // STO, L1.S1, L1.S2, IMG
  const double n = 1.5168;
  const double c = 1.0 / 51.68;
  const double t = tan_deg(5.0);

  // Marginal ray (y-nu by hand): S1 at y = 10, n u1 = -10 c (n - 1); transfer 4 mm in glass;
  // S2 plane: n'u' = n u1.
  const double u1 = -10.0 * c * (n - 1.0) / n;
  const double y2 = 10.0 + 4.0 * u1;
  require_ray(p.surfaces[1].marginal, {10.0, u1, 10.0 * c});
  require_ray(p.surfaces[2].marginal, {y2, n * u1, u1});
  // Chief ray: y_bar = 5 t at S1, n u_bar1 = t - 5 t c (n - 1).
  const double ub1 = (t - 5.0 * t * c * (n - 1.0)) / n;
  require_ray(p.surfaces[1].chief, {5.0 * t, ub1, t + 5.0 * t * c});
  require_ray(p.surfaces[2].chief, {5.0 * t + 4.0 * ub1, n * ub1, ub1});

  require_near(p.total_track, 106.363);
  // Object at infinity: F/# = EFL / EPD (Greivenkamp p. 9-35); EFL = R1 / (n - 1) = 100.
  REQUIRE(p.first_order.efl.has_value());
  require_near(p.paraxial_working_f_number, *p.first_order.efl / 20.0);
  require_near(p.paraxial_working_f_number, 1.0 / (2.0 * 10.0 * c * (n - 1.0)));
  require_near(p.paraxial_image_na, 10.0 * c * (n - 1.0));
  require_near(p.lagrange_invariant, 10.0 * t);
  require_lagrange_constant(p);
}

TEST_CASE("afocal Kepler telescope: no F/#, angular magnification -f1/f2, H constant",
          "[prescription]") {
  // f1 = 100 mm, f2 = 25 mm (thin, R = 50 and 12.5, n = 1.5), separation f1 + f2, stop
  // (r = 10) at the objective, EPD 20, field 2 deg.
  System s = base_system();
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 2.0, 1.0}}};
  add(s, stop("STO", 0.0, 10.0));
  add(s, lens("L1", 0.0, 1.5, 50.0, std::nullopt, 0.0));
  add(s, lens("L2", 125.0, 1.5, 12.5, std::nullopt, 0.0));
  const Prescription p = prescription_of(s);
  REQUIRE(p.surfaces.size() == 5);
  const double t = tan_deg(2.0);

  // Marginal ray: u = -10 / 100 after L1, y = 10 - 125 * 0.1 = -2.5 at L2, then parallel.
  // L2.S1: i = u + y c = -0.1 - 2.5 / 12.5 = -0.3, n'u' = -0.1 + 2.5 * 0.5 / 12.5 = 0.
  require_ray(p.surfaces[3].marginal, {-2.5, 0.0, -0.3});
  require_ray(p.surfaces[4].marginal, {-2.5, 0.0, 0.0});
  REQUIRE_FALSE(p.first_order.efl.has_value());
  REQUIRE_FALSE(p.paraxial_working_f_number.has_value());
  REQUIRE(p.paraxial_image_na.has_value());
  REQUIRE_THAT(*p.paraxial_image_na, WithinAbs(0.0, kAbs));
  // Angular magnification from first_order: -f1 / f2 = -4; the chief ray leaves with -4 t.
  require_near(p.first_order.angular_magnification, -4.0);
  REQUIRE(p.surfaces[4].chief.has_value());
  require_near(p.surfaces[4].chief->u, -4.0 * t);
  require_near(p.total_track, 125.0);
  require_near(p.lagrange_invariant, 10.0 * t);
  require_lagrange_constant(p);
}

TEST_CASE("Lagrange invariant constant after every event: two lenses, Cooke triplet, mirror",
          "[prescription]") {
  SECTION("two singlets with the stop between") {
    const CompiledSystem cs = compile(load_in_vacuum("m1/two_lenses_stop_between.rtt.json"));
    const Prescription p = prescription(cs, PathId{0}, 0);
    require_lagrange_constant(p);
    // Consistency with the pupils of first_order() (Greivenkamp p. 9-41: the chief ray passes
    // the centre of the stop and of both pupils, the marginal ray the rim of the EP):
    // - the chief ray has y_bar = 0 at the stop event;
    // - the marginal ray, transferred from its start to the EP plane, has height EPD / 2;
    // - after the last event the chief ray crosses the axis in the exit pupil plane.
    // Tolerance as above: 1e-12 relative to the size of the terms (1e-12 mm for y_bar = 0).
    const auto& fo = p.first_order;
    REQUIRE(fo.entrance_pupil.has_value());
    REQUIRE(fo.exit_pupil.has_value());
    REQUIRE(fo.entrance_pupil->z.has_value());
    REQUIRE(fo.entrance_pupil->diameter.has_value());
    REQUIRE(fo.exit_pupil->z.has_value());
    const auto& events = cs.path(PathId{0}).events;
    bool stop_seen = false;
    for (std::size_t k = 0; k < events.size(); ++k) {
      if (cs.surfaces()[events[k].surface].element_kind == ElementKind::Stop) {
        REQUIRE(p.surfaces[k].chief.has_value());
        REQUIRE_THAT(p.surfaces[k].chief->y, WithinAbs(0.0, kAbs));
        stop_seen = true;
        break;
      }
    }
    REQUIRE(stop_seen);
    REQUIRE(p.marginal_start.has_value());
    const auto& m = *p.marginal_start;
    require_near(m.y + (*fo.entrance_pupil->z - m.z) * m.u, 0.5 * *fo.entrance_pupil->diameter);
    const PrescriptionSurface& last = p.surfaces.back();
    REQUIRE(last.chief.has_value());
    const double travel = (*fo.exit_pupil->z - last.z) * last.chief->u;
    REQUIRE(std::abs(last.chief->y + travel) <=
            kRel * (std::abs(last.chief->y) + std::abs(travel)));
  }
  SECTION("Cooke triplet with constant indices") {
    // tests/reference/m2/cooke_triplet.rtt.json; the SCHOTT glasses are replaced by their
    // nd (N-LAK9 1.69100, N-SF5 1.67271): the invariant does not depend on the glass.
    System s = load_in_vacuum("m2/cooke_triplet.rtt.json");
    for (auto& child : s.root.children) {
      auto& e = std::get<Element>(child.value);
      if (e.material == "SCHOTT:N-LAK9") e.material = "CONST:1.69100";
      if (e.material == "SCHOTT:N-SF5") e.material = "CONST:1.67271";
    }
    const Prescription p = prescription(compile(s), PathId{0}, 0);
    REQUIRE(p.surfaces.size() == 8);
    require_lagrange_constant(p);
  }
  SECTION("concave mirror with a stop: signed index after the reflection") {
    System s = base_system();
    add(s, stop("STO", 0.0, 10.0));
    add(s, mirror("M", 100.0, -200.0));
    add(s, detector("IMG", 0.0));
    const Prescription p = prescription_of(s);
    REQUIRE(p.surfaces.size() == 3);
    REQUIRE(p.surfaces[1].n == -1.0);
    // At the mirror (c = -1/200 in global orientation): i = u + y c = -0.05;
    // n'u' = n u - y c (n' - n) = -10 (-1/200)(-2) = -0.1, so u' = -0.1 / n' = 0.1.
    require_ray(p.surfaces[1].marginal, {10.0, 0.1, -0.05});
    require_lagrange_constant(p);
    // Unfolded length: 100 mm to the mirror and 100 mm back (the z span is 0).
    require_near(p.total_track, 200.0);
    // Mirror focal length |R| / 2 = 100: F/# = 100 / 20 = 1 / (2 |n' u'|).
    require_near(p.paraxial_working_f_number, 5.0);
  }
}

TEST_CASE("same rays as seidel(): rays bit for bit, A and H to rounding", "[prescription]") {
  // The rays come from the same functions (paraxial_rays.hpp, trace_ray), so they agree bit for
  // bit. A = n i and H are the same expressions, but evaluated in two translation units
  // (seidel.cpp, prescription.cpp); an optimiser may contract y c + u or the products of H
  // differently, which changes a result by about one rounding of its largest term. So they
  // are compared to 4 ulp of the term sizes |n u| + |n y c| and |n u_bar y| + |n u y_bar|.
  constexpr double kUlp = 4.0 * std::numeric_limits<double>::epsilon();
  const auto check = [&](const CompiledSystem& cs) {
    const Prescription p = prescription(cs, PathId{0}, 0);
    const rtt::paraxial::Seidel sd = rtt::paraxial::seidel(cs, PathId{0}, 0);
    REQUIRE(p.surfaces.size() == sd.surfaces.size());
    REQUIRE(p.marginal_start.has_value());
    REQUIRE(p.chief_start.has_value());
    REQUIRE(p.lagrange_invariant.has_value());
    REQUIRE(p.marginal_start->z == sd.marginal.z);
    REQUIRE(p.marginal_start->y == sd.marginal.y);
    REQUIRE(p.marginal_start->u == sd.marginal.u);
    REQUIRE(p.chief_start->z == sd.chief.z);
    REQUIRE(p.chief_start->y == sd.chief.y);
    REQUIRE(p.chief_start->u == sd.chief.u);
    const double n0 = p.first_order.object_index;
    const PrescriptionSurface& first = p.surfaces.front();
    REQUIRE(first.marginal.has_value());
    REQUIRE(first.chief.has_value());
    REQUIRE(std::abs(*p.lagrange_invariant - sd.lagrange) <=
            kUlp * (std::abs(n0 * sd.chief.u * first.marginal->y) +
                    std::abs(n0 * sd.marginal.u * first.chief->y)));
    double n_before = n0;
    double u = sd.marginal.u;
    double u_bar = sd.chief.u;
    for (std::size_t k = 0; k < p.surfaces.size(); ++k) {
      const PrescriptionSurface& q = p.surfaces[k];
      INFO("event " << k);
      REQUIRE(q.marginal.has_value());
      REQUIRE(q.chief.has_value());
      REQUIRE(q.lagrange.has_value());
      REQUIRE(q.surface == sd.surfaces[k].surface);
      REQUIRE(q.marginal->y == sd.surfaces[k].y);
      REQUIRE(q.chief->y == sd.surfaces[k].y_bar);
      // A = n i with n before the event (Sasian L4, p. 24).
      // y c = i - u up to rounding, so |u| + |i - u| is the size of the terms of i.
      REQUIRE(std::abs(n_before * q.marginal->i - sd.surfaces[k].a) <=
              kUlp * std::abs(n_before) * (std::abs(u) + std::abs(q.marginal->i - u)) +
                  kUlp * std::abs(sd.surfaces[k].a));
      REQUIRE(std::abs(n_before * q.chief->i - sd.surfaces[k].a_bar) <=
              kUlp * std::abs(n_before) * (std::abs(u_bar) + std::abs(q.chief->i - u_bar)) +
                  kUlp * std::abs(sd.surfaces[k].a_bar));
      REQUIRE(std::abs(*q.lagrange - sd.surfaces[k].lagrange) <=
              kUlp * (std::abs(q.n * q.chief->u * q.marginal->y) +
                      std::abs(q.n * q.marginal->u * q.chief->y)));
      n_before = q.n;
      u = q.marginal->u;
      u_bar = q.chief->u;
    }
  };
  SECTION("singlet") {
    check(compile(load_in_vacuum("m1/singlet_const.rtt.json")));
  }
  SECTION("two singlets with the stop between") {
    check(compile(load_in_vacuum("m1/two_lenses_stop_between.rtt.json")));
  }
  SECTION("finite object with an object-height field") {
    System s = base_system();
    s.object = {false, Param(300.0)};
    s.fields = {rtt::model::FieldType::ObjectHeight, {{0.0, 3.0, 1.0}}};
    add(s, stop("STO", 0.0, 10.0));
    add(s, lens("L", 20.0, 1.6, 40.0, -120.0, 5.0));
    check(compile(s));
  }
}

TEST_CASE("without a stop or pupil the affected values are none, without an error",
          "[prescription]") {
  SECTION("no stop, object at infinity, EPD: marginal ray only") {
    System s = base_system();
    add(s, lens("L", 0.0, 1.5, 50.0, std::nullopt, 0.0));
    const Prescription p = prescription_of(s);
    REQUIRE_FALSE(p.first_order.entrance_pupil.has_value());
    REQUIRE(p.surfaces.size() == 2);
    require_ray(p.surfaces[1].marginal, {10.0, -0.1, -0.1 / 1.5});
    REQUIRE_FALSE(p.surfaces[1].chief.has_value());
    REQUIRE_FALSE(p.surfaces[1].lagrange.has_value());
    REQUIRE_FALSE(p.chief_start.has_value());
    REQUIRE_FALSE(p.lagrange_invariant.has_value());
    require_near(p.paraxial_working_f_number, 5.0);
    require_near(p.paraxial_image_na, 0.1);
  }
  SECTION("no stop, object at infinity, image-space F-number: EPD = EFL / F#") {
    System s = base_system();
    s.aperture = {rtt::model::SystemApertureType::ImageSpaceFNumber, Param(4.0)};
    add(s, lens("L", 0.0, 1.5, 50.0, std::nullopt, 0.0));
    const Prescription p = prescription_of(s);
    REQUIRE(p.surfaces.size() == 2);
    require_ray(p.surfaces[1].marginal, {12.5, -0.125, -0.125 / 1.5});
    require_near(p.paraxial_working_f_number, 4.0);
  }
  SECTION("paraboloid mirror without a stop: F/# = EFL / EPD") {
    // tests/reference/m1/paraboloid_mirror.rtt.json: R = -200, EPD 60, focus at z = -100.
    const CompiledSystem cs = compile(load_in_vacuum("m1/paraboloid_mirror.rtt.json"));
    const Prescription p = prescription(cs, PathId{0}, 0);
    REQUIRE(p.surfaces.size() == 2);
    require_near(p.total_track, 100.0);
    require_near(p.paraxial_working_f_number, 100.0 / 60.0);
    require_near(p.paraxial_image_na, 0.3);
    REQUIRE_FALSE(p.lagrange_invariant.has_value());
  }
  SECTION("no stop, finite object: no rays") {
    System s = base_system();
    s.object = {false, Param(200.0)};
    s.fields = {rtt::model::FieldType::ObjectHeight, {{0.0, 1.0, 1.0}}};
    add(s, lens("L", 0.0, 1.5, 50.0, std::nullopt, 0.0));
    const Prescription p = prescription_of(s);
    REQUIRE(p.surfaces.size() == 2);
    REQUIRE(p.surfaces[1].z == 0.0);
    REQUIRE(p.surfaces[1].n == 1.0);
    REQUIRE_FALSE(p.surfaces[1].marginal.has_value());
    REQUIRE_FALSE(p.surfaces[1].chief.has_value());
    REQUIRE_FALSE(p.marginal_start.has_value());
    REQUIRE_FALSE(p.paraxial_working_f_number.has_value());
    REQUIRE_FALSE(p.paraxial_image_na.has_value());
    REQUIRE_FALSE(p.lagrange_invariant.has_value());
    require_near(p.object_distance, 200.0);
  }
  SECTION("entrance pupil at infinity (stop in the rear focal plane): marginal ray only") {
    System s = base_system();
    add(s, lens("L", 0.0, 1.5, 50.0, std::nullopt, 0.0));
    add(s, stop("STO", 100.0, 10.0));
    const Prescription p = prescription_of(s);
    REQUIRE(p.first_order.entrance_pupil.has_value());
    REQUIRE_FALSE(p.first_order.entrance_pupil->z.has_value());
    REQUIRE(p.surfaces.size() == 3);
    require_ray(p.surfaces[2].marginal, {0.0, -0.1, -0.1});
    REQUIRE_FALSE(p.surfaces[2].chief.has_value());
    REQUIRE_FALSE(p.lagrange_invariant.has_value());
    require_near(p.paraxial_working_f_number, 5.0);
  }
}

TEST_CASE("prescription input errors", "[prescription]") {
  SECTION("not rotationally symmetric") {
    const CompiledSystem cs =
        compile(rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m0/michelson.rtt.json"));
    REQUIRE_THROWS_AS(prescription(cs, PathId{0}, 0), ParaxialError);
  }
  SECTION("path and wavelength indices") {
    const CompiledSystem cs = compile(load_in_vacuum("m1/singlet_const.rtt.json"));
    REQUIRE_THROWS_AS(prescription(cs, PathId{1}, 0), ParaxialError);
    REQUIRE_THROWS_AS(prescription(cs, PathId{0}, 3), ParaxialError);
  }
  SECTION("invalid field definition for the chief ray") {
    System s = load_in_vacuum("m1/singlet_const.rtt.json");
    s.fields = {rtt::model::FieldType::ObjectHeight, {{0.0, 1.0, 1.0}}};
    REQUIRE_THROWS_AS(prescription(compile(s), PathId{0}, 0), ParaxialError);
  }
}
