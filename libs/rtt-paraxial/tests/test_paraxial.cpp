#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/air.hpp"
#include "rtt/paraxial/paraxial.hpp"

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
using rtt::paraxial::first_order;
using rtt::paraxial::FirstOrder;
using rtt::paraxial::ParaxialError;
using rtt::paraxial::trace_ray;

namespace {

constexpr double kRel = 1e-12;  // relative tolerance, issue #7
constexpr double kAbs = 1e-12;  // mm, for values that are zero

// ------------------------------------------------------------------ builders -----

/// The reference values below are formulas for a surrounding index of exactly 1, so the tests
/// run in VACUUM; since #25 AIR is Ciddor air (n = 1.00027).
System base_system() {
  System s;
  s.name = "paraxial test";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
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

/// Lens with vertex at global z, index n (CONST), radii R1, R2 and centre thickness d.
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

void add(System& s, Element e) {
  s.root.children.push_back({std::move(e)});
}

CompiledSystem compile(const System& s) {
  const MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

FirstOrder first_order_of(const System& s) {
  return first_order(compile(s), PathId{0}, 0);
}

void require_rel(const std::optional<double>& actual, double expected) {
  REQUIRE(actual.has_value());
  REQUIRE_THAT(*actual, WithinRel(expected, kRel));
}

const rtt::paraxial::Pupil& entrance_pupil(const FirstOrder& fo) {
  REQUIRE(fo.entrance_pupil.has_value());
  return *fo.entrance_pupil;
}

const rtt::paraxial::Pupil& exit_pupil(const FirstOrder& fo) {
  REQUIRE(fo.exit_pupil.has_value());
  return *fo.exit_pupil;
}

// ------------------------------------------------------ reference formulas -----

/// Thick lens in air: lensmaker formula (issue #7)
///   1/f = (n - 1) [1/R1 - 1/R2 + (n - 1) d / (n R1 R2)],
///   h1 = -f (n - 1) d / (R2 n),  h2 = -f (n - 1) d / (R1 n),
/// h1, h2 = positions of the principal points H1, H2 relative to the vertices V1, V2, positive
/// in the direction of the light; ffl = V1 -> F1 = f - h1 (F1 upstream of V1 positive) and
/// bfl = V2 -> F2 = f + h2. 1/R = 0 for a plane surface.
/// Derivation by hand from the y-nu equations (Greivenkamp, OPTI-201/202 lecture notes, Sec. 9,
/// p. 9-2: n'u' = nu - y phi, phi = (n' - n) C, y' = y + u t), phi1 = (n - 1)/R1,
/// phi2 = (1 - n)/R2:
///   ray (y, nu) = (1, 0): S1 -> nu = -phi1; transfer d in glass -> y2 = 1 - d phi1 / n;
///   S2 -> nu' = -phi1 - y2 phi2 = -(phi1 + phi2 - d phi1 phi2 / n) = -1/f, which is the
///   lensmaker formula; bfl = -y2 / u' = f (1 - d phi1 / n) = f (1 - (n - 1) d / (n R1)) = f + h2.
///   Reversed light gives ffl = f (1 - d phi2 / n) = f (1 + (n - 1) d / (n R2)) = f - h1.
/// See also Hecht, Optics, Ch. 6 (thick lenses), the chapter cited in issue #7.
struct ThickLens {
  double f, h1, h2, ffl, bfl;
};

ThickLens hecht_thick_lens(double n, double inv_r1, double inv_r2, double d) {
  const double f = 1.0 / ((n - 1.0) * (inv_r1 - inv_r2 + (n - 1.0) * d * inv_r1 * inv_r2 / n));
  const double h1 = -f * (n - 1.0) * d * inv_r2 / n;
  const double h2 = -f * (n - 1.0) * d * inv_r1 / n;
  return {f, h1, h2, f - h1, f + h2};
}

/// Image distance b from the Gaussian lens equation 1/a + 1/b = 1/f with object distance a,
/// distances positive for a real object in front of and a real image behind the lens
/// (Gaussian lens equation). Derivation from the y-nu equations (Greivenkamp, OPTI-201/202
/// lecture notes, Sec. 9, p. 9-2) in air: an axial ray with slope u from the object reaches the
/// lens at y = a u, leaves with u' = u - y / f, and meets the axis b = -y / u' behind it, so
/// 1/b = 1/f - 1/a. For a mirror the same holds with f = |R|/2 along the reflected direction.
/// See also Hecht, Optics, Ch. 5 (thin lenses).
double gauss_image(double a, double f) {
  return 1.0 / (1.0 / f - 1.0 / a);
}

void check_thick_lens(double n, std::optional<double> r1, std::optional<double> r2, double d) {
  const double z1 = 10.0;  // vertex V1
  System s = base_system();
  add(s, lens("L", z1, n, r1, r2, d));
  const FirstOrder fo = first_order_of(s);

  // n is passed as text "CONST:<n>" with std::to_string (6 decimals): use the same value.
  const double n_used = std::stod(std::to_string(n));
  const ThickLens ref = hecht_thick_lens(n_used, r1 ? 1.0 / *r1 : 0.0, r2 ? 1.0 / *r2 : 0.0, d);
  INFO("f = " << ref.f << ", h1 = " << ref.h1 << ", h2 = " << ref.h2);

  require_rel(fo.efl, ref.f);
  REQUIRE_THAT(fo.power, WithinRel(1.0 / ref.f, kRel));
  require_rel(fo.front_focal_length, ref.f);
  require_rel(fo.rear_focal_length, ref.f);
  require_rel(fo.ffl, ref.ffl);
  require_rel(fo.bfl, ref.bfl);
  require_rel(fo.front_principal_z, z1 + ref.h1);
  require_rel(fo.rear_principal_z, z1 + d + ref.h2);
  require_rel(fo.front_focal_z, z1 - ref.ffl);
  require_rel(fo.rear_focal_z, z1 + d + ref.bfl);
  require_rel(fo.image_z, z1 + d + ref.bfl);  // object at infinity: image in F'
  REQUIRE(fo.image_direction == 1);
  REQUIRE(fo.object_index == 1.0);
  REQUIRE(fo.image_index == 1.0);
  REQUIRE_FALSE(fo.lateral_magnification.has_value());
  REQUIRE_FALSE(fo.entrance_pupil.has_value());  // no stop on the path
}

}  // namespace

// ------------------------------------------------------------------- tests -----

TEST_CASE("thick lens: EFL, BFL, FFL and principal planes from the lensmaker formula",
          "[paraxial]") {
  // Mandatory reference case (docs/architecture.md): thick lens, relative 1e-12.
  SECTION("biconvex") {
    check_thick_lens(1.5168, 51.68, -80.0, 6.0);
  }
  SECTION("positive meniscus") {
    check_thick_lens(1.7, 30.0, 60.0, 8.0);
  }
  SECTION("biconcave") {
    check_thick_lens(1.6, -40.0, 55.0, 3.0);
  }
  SECTION("plano-convex, plane side first") {
    check_thick_lens(1.5, std::nullopt, -50.0, 5.0);
  }
}

TEST_CASE("reference singlet with CONST:1.5168", "[paraxial]") {
  // tests/reference/m1/singlet_const.rtt.json: plano-convex, R1 = 51.68 mm, d = 4 mm, vertex
  // at z = 5, stop (radius 10) at z = 0, EPD 20 mm. f = R1 / (n - 1) = 100 mm.
  // Lensmaker reference for n_outside = 1: the file's AIR is replaced by VACUUM here (#25).
  System singlet =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  singlet.environment.medium = "VACUUM";
  const CompiledSystem cs = compile(singlet);
  const ThickLens ref = hecht_thick_lens(1.5168, 1.0 / 51.68, 0.0, 4.0);
  REQUIRE_THAT(ref.f, WithinRel(100.0, 1e-12));

  for (std::uint16_t wl = 0; wl < 3; ++wl) {
    const FirstOrder fo = first_order(cs, PathId{0}, wl);
    require_rel(fo.efl, ref.f);
    require_rel(fo.bfl, ref.bfl);
    require_rel(fo.rear_focal_z, 9.0 + ref.bfl);
    // Principal planes of the mandatory thick-lens case: H = V1 + h1, H' = V2 + h2.
    require_rel(fo.front_principal_z, 5.0 + ref.h1);
    require_rel(fo.rear_principal_z, 9.0 + ref.h2);

    // Stop in front of the lens: entrance pupil = stop, EPD from the system aperture.
    REQUIRE(fo.entrance_pupil.has_value());
    REQUIRE_THAT(*entrance_pupil(fo).z, WithinAbs(0.0, kAbs));
    require_rel(entrance_pupil(fo).diameter, 20.0);
    // Exit pupil: the stop imaged by the lens. Object distance from H (h1 = 0 for a plane
    // second surface, so H = V1): a = 5 mm; image distance from H' with Gauss.
    const double a = 5.0;
    const double b = gauss_image(a, ref.f);
    REQUIRE(fo.exit_pupil.has_value());
    require_rel(exit_pupil(fo).z, 9.0 + ref.h2 + b);
    require_rel(exit_pupil(fo).diameter, 20.0 * std::abs(b / a));
  }
}

TEST_CASE("thin lens limit d -> 0", "[paraxial]") {
  const double n = 1.5;
  const double r1 = 50.0;
  const double r2 = -75.0;
  // d = 0 is the thin lens: 1/f = (n - 1)(1/R1 - 1/R2), principal planes at the lens.
  const double f_thin = 1.0 / ((n - 1.0) * (1.0 / r1 - 1.0 / r2));
  check_thick_lens(n, r1, r2, 0.0);

  // For d > 0 the deviation from the thin lens shrinks linearly with d.
  const auto efl = [&](double d) {
    System s = base_system();
    add(s, lens("L", 0.0, n, r1, r2, d));
    return *first_order_of(s).efl;
  };
  const double e1 = std::abs(efl(1e-2) - f_thin);
  const double e2 = std::abs(efl(1e-3) - f_thin);
  REQUIRE(e2 < e1);
  REQUIRE_THAT(e1 / e2, WithinRel(10.0, 1e-3));
}

TEST_CASE("two thin lenses: 1/f = 1/f1 + 1/f2 - d/(f1 f2)", "[paraxial]") {
  // Thin lenses as two surfaces at the same z (d = 0), so the reference values are exact.
  // f1 = R / (n - 1) = 100 mm (plano-convex), f2 = 1 / ((n - 1)(2/R)) = 100 mm (biconvex).
  // 1/f = 1/f1 + 1/f2 - d/(f1 f2) is the formula of issue #7. BFL and FFL by hand from the
  // y-nu equations (Greivenkamp, OPTI-201/202 lecture notes, Sec. 9 "Paraxial Raytracing",
  // p. 9-2: n'u' = nu - y phi, y' = y + u t) in air, phi_i = 1/f_i:
  //   ray (y, u) = (1, 0): lens 1 -> u = -phi1; transfer d -> y2 = 1 - d phi1;
  //   lens 2 -> u' = -phi1 - y2 phi2 = -phi;  bfl = -y2 / u' = (1 - d/f1) f.
  // The front side follows by symmetry (light reversed): ffl = (1 - d/f2) f. Hence
  //   bfl = f2 (d - f1) / (d - (f1 + f2)),  ffl = f1 (d - f2) / (d - (f1 + f2)).
  const double f1 = 100.0;
  const double f2 = 100.0;
  const double d = 30.0;
  System s = base_system();
  add(s, lens("L1", 0.0, 1.5, 50.0, std::nullopt, 0.0));
  add(s, lens("L2", d, 1.5, 100.0, -100.0, 0.0));
  const FirstOrder fo = first_order_of(s);

  const double f = 1.0 / (1.0 / f1 + 1.0 / f2 - d / (f1 * f2));
  require_rel(fo.efl, f);
  require_rel(fo.bfl, f2 * (d - f1) / (d - (f1 + f2)));
  require_rel(fo.ffl, f1 * (d - f2) / (d - (f1 + f2)));
}

TEST_CASE("spherical mirror: |f| = |R| / 2", "[paraxial]") {
  // The issue states f = R/2 in magnitude: with our radius convention (centre of curvature on
  // the +z side means R > 0) a concave mirror lit from -z has R < 0, and its EFL is positive
  // (converging), so EFL = -R/2 = |R|/2.
  SECTION("concave, object at infinity") {
    System s = base_system();
    add(s, mirror("M", 0.0, -200.0));
    const FirstOrder fo = first_order_of(s);
    require_rel(fo.efl, 100.0);
    require_rel(fo.rear_focal_length, 100.0);
    require_rel(fo.bfl, 100.0);
    require_rel(fo.rear_focal_z, -100.0);
    // Both focal points lie at R/2 in front of the mirror, both principal points at the vertex.
    require_rel(fo.front_focal_z, -100.0);
    require_rel(fo.ffl, 100.0);
    REQUIRE(fo.front_principal_z.has_value());
    REQUIRE_THAT(*fo.front_principal_z, WithinAbs(0.0, kAbs));
    REQUIRE(fo.rear_principal_z.has_value());
    REQUIRE_THAT(*fo.rear_principal_z, WithinAbs(0.0, kAbs));
    REQUIRE(fo.image_direction == -1);
    REQUIRE(fo.image_index == 1.0);
  }
  SECTION("convex, diverging") {
    System s = base_system();
    add(s, mirror("M", 0.0, 200.0));
    require_rel(first_order_of(s).efl, -100.0);
  }
  SECTION("concave, object at the centre of curvature images onto itself, m = -1") {
    System s = base_system();
    s.object = {false, Param(200.0)};  // object at z = -200 = centre of curvature
    add(s, mirror("M", 0.0, -200.0));
    const FirstOrder fo = first_order_of(s);
    require_rel(fo.image_z, -200.0);
    require_rel(fo.lateral_magnification, -1.0);
  }
  SECTION("concave, finite object: mirror equation 1/a + 1/b = 2/|R|") {
    System s = base_system();
    s.object = {false, Param(300.0)};
    add(s, mirror("M", 0.0, -200.0));
    const FirstOrder fo = first_order_of(s);
    const double b = gauss_image(300.0, 100.0);  // 150 mm in front of the mirror
    require_rel(fo.image_z, -b);
    require_rel(fo.lateral_magnification, -b / 300.0);
  }
  SECTION("mirror turned by 180 deg about x acts with the opposite curvature") {
    // Local z points to -z globally, so R = +200 in local coordinates is concave for light
    // from -z. Rx(180 deg) is not exact in floating point (sin pi ~ 1.2e-16); it still counts
    // as rotationally symmetric.
    System s = base_system();
    Element m = mirror("M", 0.0, 200.0);
    m.pose.rotation_deg[0] = Param(180.0);
    add(s, m);
    require_rel(first_order_of(s).efl, 100.0);
  }
}

TEST_CASE("even asphere is paraxially its base conic", "[paraxial]") {
  // Paraxial refraction uses only the vertex curvature: phi = (n' - n) C (Greivenkamp,
  // OPTI-201/202 lecture notes, Sec. 9, p. 9-2; docs/quellen.md). The even asphere is
  // z = conic(r) + A4 r^4 + A6 r^6 + ... (model::EvenAsphere); the polynomial terms have zero
  // curvature at the vertex (d^2/dr^2 of A_2k r^2k vanishes at r = 0 for k >= 2), and the conic
  // constant does not enter the vertex curvature either. So an asphere must give exactly the
  // first-order data of a sphere with the same radius: bit-identical, because both reach the
  // y-nu trace with the same c = 1/R.
  // Not tested: an asphere on a flat base (c = 0, only A4 != 0) would be paraxially a plane, but
  // the model excludes it: model::validate requires a finite, non-zero radius for EvenAsphere
  // just as for Conic, so such a surface never reaches rtt-paraxial.
  const auto with_s1 = [](const rtt::model::BaseShape& base) {
    System s = base_system();
    Element l = lens("L", 10.0, 1.5168, std::nullopt, -80.0, 6.0);
    l.surfaces[0].shape.base = base;
    add(s, l);
    Element m = mirror("M", 60.0, std::nullopt);
    m.surfaces[0].shape.base = base;  // mirror with the same base shape, image space -z
    add(s, m);
    return first_order_of(s);
  };
  const FirstOrder sphere = with_s1(rtt::model::Conic{Param(51.68), Param(0.0)});
  const FirstOrder conic = with_s1(rtt::model::Conic{Param(51.68), Param(-0.5)});
  const FirstOrder asphere = with_s1(rtt::model::EvenAsphere{
      Param(51.68), Param(-0.5), {Param(1e-6), Param(-2e-9), Param(3e-12)}});

  REQUIRE(sphere.efl.has_value());  // focal system: all cardinal points below are set
  for (const FirstOrder* fo : {&conic, &asphere}) {
    REQUIRE(fo->efl.has_value());
    REQUIRE(*fo->efl == *sphere.efl);
    REQUIRE(*fo->bfl == *sphere.bfl);
    REQUIRE(*fo->ffl == *sphere.ffl);
    REQUIRE(*fo->front_principal_z == *sphere.front_principal_z);
    REQUIRE(*fo->rear_principal_z == *sphere.rear_principal_z);
    REQUIRE(fo->image_direction == sphere.image_direction);
  }
  REQUIRE(sphere.image_direction == -1);  // the mirror is part of the path
}

TEST_CASE("double pass: thin lens on a plane mirror refracts with negative index", "[paraxial]") {
  // Thin lens f = 100 mm and a plane mirror, both at z = 0; explicit path through the lens,
  // onto the mirror and back through the lens. The light passes the lens twice with zero
  // spacing, so Phi = 2/f: EFL = f/2 = 50 mm, F' at 50 mm in front of the lens (z = -50).
  // This checks refraction after a reflection (signed index n < 0).
  System s = base_system();
  add(s, lens("L", 0.0, 1.5, 50.0, std::nullopt, 0.0));
  add(s, mirror("M", 0.0, std::nullopt));
  s.paths = {{"double pass",
              false,
              {{SurfaceId("L.S1"), rtt::model::EventKind::Refract, 0},
               {SurfaceId("L.S2"), rtt::model::EventKind::Refract, 0},
               {SurfaceId("M.S"), rtt::model::EventKind::Reflect, 0},
               {SurfaceId("L.S2"), rtt::model::EventKind::Refract, 0},
               {SurfaceId("L.S1"), rtt::model::EventKind::Refract, 0}}}};
  const FirstOrder fo = first_order_of(s);
  require_rel(fo.efl, 50.0);
  require_rel(fo.bfl, 50.0);
  require_rel(fo.rear_focal_z, -50.0);
  REQUIRE(fo.image_direction == -1);
}

TEST_CASE("thin lens in AIR: absolute indices give EFL = 1/Phi and f' = n_air/Phi", "[paraxial]") {
  // Since #25 AIR is Ciddor air with n_air > 1 (decision D2: EFL = 1/Phi in the absolute sense).
  // Thin plano-convex lens, R = 50 mm, n = 1.5 (CONST, absolute) in AIR at 20 degC, 1 atm:
  // phi = (n - n_air) / R (Greivenkamp, OPTI-201/202, Sec. 9, p. 9-2: phi = (n' - n) C; the plane
  // second surface has no power), EFL = 1/phi, and with y = 1, u = 0 as in the test below:
  // n_air u' = -phi, so f' = n_air / phi (p. 9-12: f'_R = n'/phi) and BFL = f' for the thin lens
  // at z = 0; p. 9-14: f_F = -n/phi, front_focal_length = -f_F = n_air / phi.
  System s = base_system();
  s.environment.medium = "AIR";
  add(s, lens("L", 0.0, 1.5, 50.0, std::nullopt, 0.0));
  const FirstOrder fo = first_order_of(s);
  const double n_air = rtt::material::ciddor_air_index(0.5876, 20.0, 1.0);
  const double phi = (1.5 - n_air) / 50.0;
  REQUIRE(n_air > 1.0002);
  require_rel(fo.efl, 1.0 / phi);
  require_rel(fo.rear_focal_length, n_air / phi);
  require_rel(fo.front_focal_length, n_air / phi);
  require_rel(fo.bfl, n_air / phi);
  REQUIRE_THAT(fo.object_index, WithinRel(n_air, kRel));
  REQUIRE_THAT(fo.image_index, WithinRel(n_air, kRel));
}

TEST_CASE("image space in glass: f' = n' f", "[paraxial]") {
  // Single refracting surface air -> n' = 1.5, R = 50 mm (explicit path ending in the glass).
  // phi = (n' - n) / R = 0.01 / mm, f = n / phi = 100 mm, f' = n' / phi = 150 mm, and F' lies
  // f' behind the vertex (Greivenkamp, OPTI-201/202 lecture notes, Sec. 9, p. 9-2 with
  // y = 1, u = 0: n'u' = -phi, F' at -y/u' = n'/phi).
  System s = base_system();
  add(s, lens("L", 0.0, 1.5, 50.0, std::nullopt, 5.0));
  s.paths = {{"into glass", false, {{SurfaceId("L.S1"), rtt::model::EventKind::Refract, 0}}}};
  const FirstOrder fo = first_order_of(s);
  REQUIRE(fo.image_index == 1.5);
  require_rel(fo.efl, 100.0);
  require_rel(fo.front_focal_length, 100.0);
  require_rel(fo.rear_focal_length, 150.0);
  require_rel(fo.bfl, 150.0);
  require_rel(fo.ffl, 100.0);
  require_rel(fo.rear_focal_z, 150.0);
  REQUIRE(fo.front_principal_z.has_value());
  REQUIRE_THAT(*fo.front_principal_z, WithinAbs(0.0, kAbs));
  REQUIRE(fo.rear_principal_z.has_value());
  REQUIRE_THAT(*fo.rear_principal_z, WithinAbs(0.0, kAbs));
}

TEST_CASE("thin lens with a finite object: Gauss lens equation", "[paraxial]") {
  System s = base_system();
  s.object = {false, Param(300.0)};
  add(s, lens("L", 0.0, 1.5, 50.0, std::nullopt, 0.0));  // f = 100 mm
  const FirstOrder fo = first_order_of(s);
  const double b = gauss_image(300.0, 100.0);  // 150 mm
  require_rel(fo.image_z, b);
  require_rel(fo.lateral_magnification, -b / 300.0);
}

TEST_CASE("pupils with the stop in front of, between and behind the lenses", "[paraxial]") {
  // Thin lenses f = 100 mm (plano-convex, R = 50, n = 1.5). The pupils are the images of the
  // stop; the hand values use the Gauss equation with the stop as object (for the entrance
  // pupil with the light reversed: the stop lies a mm behind the lens, its image b mm in front
  // of the lens in reversed light, i.e. at z_lens - b).
  const double f = 100.0;
  const double r = 5.0;

  SECTION("between two lenses, stop size") {
    System s = base_system();
    s.aperture = {rtt::model::SystemApertureType::StopSize, Param(1.0)};
    add(s, lens("L1", 0.0, 1.5, 50.0, std::nullopt, 0.0));
    add(s, stop("STO", 20.0, r));
    add(s, lens("L2", 50.0, 1.5, 50.0, std::nullopt, 0.0));
    const FirstOrder fo = first_order_of(s);

    const double a1 = 20.0;
    const double b1 = gauss_image(a1, f);  // -25: virtual, behind lens 1
    require_rel(entrance_pupil(fo).z, 0.0 - b1);
    require_rel(entrance_pupil(fo).diameter, 2.0 * r * std::abs(b1 / a1));
    const double a2 = 30.0;
    const double b2 = gauss_image(a2, f);  // -300/7: virtual, in front of lens 2
    require_rel(exit_pupil(fo).z, 50.0 + b2);
    require_rel(exit_pupil(fo).diameter, 2.0 * r * std::abs(b2 / a2));
  }
  SECTION("between two lenses, entrance pupil diameter given") {
    System s = base_system();
    s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
    add(s, lens("L1", 0.0, 1.5, 50.0, std::nullopt, 0.0));
    add(s, stop("STO", 20.0, r));
    add(s, lens("L2", 50.0, 1.5, 50.0, std::nullopt, 0.0));
    const FirstOrder fo = first_order_of(s);
    const double m_ep = std::abs(gauss_image(20.0, f) / 20.0);  // stop -> entrance pupil
    const double m_xp = std::abs(gauss_image(30.0, f) / 30.0);  // stop -> exit pupil
    require_rel(entrance_pupil(fo).diameter, 10.0);
    require_rel(exit_pupil(fo).diameter, 10.0 / m_ep * m_xp);
  }
  SECTION("behind the lens") {
    System s = base_system();
    s.aperture = {rtt::model::SystemApertureType::StopSize, Param(1.0)};
    add(s, lens("L1", 0.0, 1.5, 50.0, std::nullopt, 0.0));
    add(s, stop("STO", 50.0, r));
    const FirstOrder fo = first_order_of(s);
    const double b = gauss_image(50.0, f);  // -100
    require_rel(entrance_pupil(fo).z, 0.0 - b);
    require_rel(entrance_pupil(fo).diameter, 2.0 * r * std::abs(b / 50.0));
    require_rel(exit_pupil(fo).z, 50.0);
    require_rel(exit_pupil(fo).diameter, 2.0 * r);
  }
  SECTION("in front of the lens, image-space F-number") {
    System s = base_system();
    s.aperture = {rtt::model::SystemApertureType::ImageSpaceFNumber, Param(4.0)};
    add(s, stop("STO", -10.0, r));
    add(s, lens("L1", 0.0, 1.5, 50.0, std::nullopt, 0.0));
    const FirstOrder fo = first_order_of(s);
    require_rel(entrance_pupil(fo).z, -10.0);
    require_rel(entrance_pupil(fo).diameter, f / 4.0);  // EPD = EFL / F#
  }
  SECTION("object-space NA needs a finite object") {
    System s = base_system();
    s.aperture = {rtt::model::SystemApertureType::ObjectSpaceNA, Param(0.05)};
    s.object = {false, Param(300.0)};
    add(s, stop("STO", 0.0, r));
    add(s, lens("L1", 0.0, 1.5, 50.0, std::nullopt, 0.0));
    // Paraxial u = NA / n; entrance pupil at z = 0, object at z = -300.
    require_rel(entrance_pupil(first_order_of(s)).diameter, 2.0 * 0.05 * 300.0);

    s.object = {true, Param(0.0)};
    const FirstOrder at_infinity = first_order_of(s);
    REQUIRE(at_infinity.entrance_pupil.has_value());
    REQUIRE_FALSE(at_infinity.entrance_pupil->diameter.has_value());
  }
  SECTION("non-circular stop is an error") {
    System s = base_system();
    Element st = stop("STO", 0.0, r);
    st.surfaces[0].aperture = rtt::model::RectangularAperture{5.0, 3.0};
    add(s, st);
    add(s, lens("L1", 10.0, 1.5, 50.0, std::nullopt, 0.0));
    REQUIRE_THROWS_AS(first_order_of(s), ParaxialError);
  }
}

TEST_CASE("afocal Kepler telescope: no focal length, angular magnification -f1/f2", "[paraxial]") {
  // f1 = 100 mm, f2 = 25 mm (R = 12.5, n = 1.5), separation f1 + f2; stop at the objective.
  System s = base_system();
  add(s, stop("STO", 0.0, 10.0));
  add(s, lens("L1", 0.0, 1.5, 50.0, std::nullopt, 0.0));
  add(s, lens("L2", 125.0, 1.5, 12.5, std::nullopt, 0.0));
  const FirstOrder fo = first_order_of(s);
  REQUIRE_THAT(fo.power, WithinAbs(0.0, 1e-14));
  REQUIRE_FALSE(fo.efl.has_value());
  REQUIRE_FALSE(fo.bfl.has_value());
  REQUIRE_FALSE(fo.rear_focal_z.has_value());
  REQUIRE_FALSE(fo.image_z.has_value());
  require_rel(fo.angular_magnification, -4.0);
}

TEST_CASE("plane mirror: angular magnification +1 along the propagation direction", "[paraxial]") {
  System s = base_system();
  add(s, stop("STO", -50.0, 5.0));
  add(s, mirror("M", 0.0, std::nullopt));
  const FirstOrder fo = first_order_of(s);
  REQUIRE_FALSE(fo.efl.has_value());
  require_rel(fo.angular_magnification, 1.0);
  REQUIRE(fo.image_direction == -1);
}

TEST_CASE("trace_ray follows the y-nu equations", "[paraxial]") {
  // Thin lens f = 100 mm at z = 0, then a concave mirror R = -100 at z = 50.
  System s = base_system();
  add(s, lens("L", 0.0, 1.5, 50.0, std::nullopt, 0.0));
  add(s, mirror("M", 50.0, -100.0));
  const CompiledSystem cs = compile(s);
  const auto rays = trace_ray(cs, PathId{0}, 0, -10.0, 2.0, 0.0);
  REQUIRE(rays.size() == 3);
  // L.S1 (z = 0): y = 2, n'u' = nu - y phi with phi = c (n' - n) = (1/50)(0.5) = 0.01.
  REQUIRE_THAT(rays[0].y, WithinRel(2.0, kRel));
  REQUIRE_THAT(rays[0].n, WithinRel(1.5, kRel));
  REQUIRE_THAT(rays[0].u, WithinRel(-0.02 / 1.5, kRel));
  // L.S2 (plane, z = 0): n'u' = nu = -0.02.
  REQUIRE_THAT(rays[1].u, WithinRel(-0.02, kRel));
  REQUIRE_THAT(rays[1].n, WithinRel(1.0, kRel));
  // Transfer to the mirror: y = 2 - 50 * 0.02 = 1. Reflection: n' = -1,
  // phi = c (n' - n) = (-1/100)(-2) = 0.02, n'u' = nu - y phi = -0.02 - 0.02 = -0.04, u' = 0.04.
  REQUIRE(rays[2].z == 50.0);
  REQUIRE_THAT(rays[2].y, WithinRel(1.0, kRel));
  REQUIRE(rays[2].n == -1.0);
  REQUIRE_THAT(rays[2].u, WithinRel(0.04, kRel));
}

TEST_CASE("paths that are not rotationally symmetric are rejected", "[paraxial]") {
  SECTION("Michelson with a tilted beam splitter") {
    const CompiledSystem cs =
        compile(rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m0/michelson.rtt.json"));
    REQUIRE_THROWS_AS(first_order(cs, PathId{0}, 0), ParaxialError);
    REQUIRE_THROWS_AS(trace_ray(cs, PathId{0}, 0, 0.0, 1.0, 0.0), ParaxialError);
  }
  SECTION("decentred lens") {
    System s = base_system();
    Element l = lens("L", 0.0, 1.5, 50.0, -50.0, 5.0);
    l.pose.position[0] = Param(0.1);
    add(s, l);
    REQUIRE_THROWS_AS(first_order_of(s), ParaxialError);
  }
  SECTION("slightly tilted lens") {
    System s = base_system();
    Element l = lens("L", 0.0, 1.5, 50.0, -50.0, 5.0);
    l.pose.rotation_deg[1] = Param(1e-3);
    add(s, l);
    REQUIRE_THROWS_AS(first_order_of(s), ParaxialError);
  }
  SECTION("rotation about the axis is fine") {
    System s = base_system();
    Element l = lens("L", 0.0, 1.5, 50.0, std::nullopt, 0.0);
    l.pose.rotation_deg[2] = Param(37.0);
    add(s, l);
    require_rel(first_order_of(s).efl, 100.0);
  }
  SECTION("phase layer") {
    // An order != 0 at a phase surface has no paraxial model (ADR 0025, point 4; order 0 is
    // accepted since #127, see "order 0 at a phase surface"). Order -1 on
    // refraction, distinct from the section below (order 1).
    System s = base_system();
    Element l = lens("L", 0.0, 1.5, 50.0, -50.0, 5.0);
    l.surfaces[0].phases.push_back(rtt::model::RadialPhase{Param(10.0), {Param(1.0)}});
    add(s, l);
    s.paths = {{"main",
                false,
                {{SurfaceId("L.S1"), rtt::model::EventKind::Refract, -1},
                 {SurfaceId("L.S2"), rtt::model::EventKind::Refract, 0}}}};
    REQUIRE_THROWS_AS(first_order_of(s), ParaxialError);
  }
  SECTION("diffraction event") {
    System s = base_system();
    Element l = lens("L", 0.0, 1.5, 50.0, -50.0, 5.0);
    l.surfaces[0].phases.push_back(rtt::model::RadialPhase{Param(10.0), {Param(1.0)}});
    add(s, l);
    s.paths = {{"main",
                false,
                {{SurfaceId("L.S1"), rtt::model::EventKind::Refract, 1},
                 {SurfaceId("L.S2"), rtt::model::EventKind::Refract, 0}}}};
    REQUIRE_THROWS_AS(first_order_of(s), ParaxialError);
  }
}

TEST_CASE("order 0 at a phase surface is the surface without phase layer (#127)", "[paraxial]") {
  // ADR 0025, points 2 and 4: first_order accepts a surface with phase layers when the path
  // takes order 0 there; the result is bitwise that of the same lens without phase layer.
  System plain = base_system();
  add(plain, lens("L", 0.0, 1.5, 50.0, -50.0, 5.0));
  System phased = base_system();
  Element l = lens("L", 0.0, 1.5, 50.0, -50.0, 5.0);
  l.surfaces[0].phases.push_back(rtt::model::RadialPhase{Param(10.0), {Param(1.0)}});
  l.surfaces[1].phases.push_back(rtt::model::LinearGrating{Param(300.0), 0.0});
  add(phased, l);
  const FirstOrder a = first_order_of(plain);
  const FirstOrder b = first_order_of(phased);
  REQUIRE(b.power == a.power);
  REQUIRE(b.efl == a.efl);
  REQUIRE(b.ffl == a.ffl);
  REQUIRE(b.bfl == a.bfl);
  REQUIRE(b.front_principal_z == a.front_principal_z);
  REQUIRE(b.rear_principal_z == a.rear_principal_z);
  REQUIRE(b.image_z == a.image_z);
}

TEST_CASE("invalid path id or wavelength index is an error", "[paraxial]") {
  System s = base_system();
  add(s, lens("L", 0.0, 1.5, 50.0, -50.0, 5.0));
  const CompiledSystem cs = compile(s);
  REQUIRE_THROWS_AS(first_order(cs, PathId{1}, 0), ParaxialError);
  REQUIRE_THROWS_AS(first_order(cs, PathId{0}, 1), ParaxialError);
  REQUIRE_THROWS_AS(trace_ray(cs, PathId{0}, 7, 0.0, 1.0, 0.0), ParaxialError);
}

// ------------------------------------------------ afocal threshold (#35, item B9) -----
// A thick biconvex lens R1 = r, R2 = -r, n = 1.5 has phi1 = phi2 = 0.5 / r and the power
// Phi = phi1 + phi2 - (d / n) phi1 phi2 (lensmaker formula, as in check_thick_lens), so it is
// afocal for d = 6 r and has Phi = (6 - d) / (6 r^2) near it. The threshold is relative:
// afocal iff |Phi| <= 16 N u S with S the power of the same y-nu trace in magnitudes, the
// transfers with (|z'| + |z|) / |n| (derivation in paraxial.cpp).

TEST_CASE("afocal threshold: a micro lens is afocal despite rounding of its coordinates",
          "[paraxial][afocal]") {
  // r = 3 um at z = 10 mm: the computed Phi is about -1.3e-11 / mm (the vertex z = 10 + d is
  // only known to u * 10 mm, and dPhi/dd = -phi1 phi2 / n = -1.9e4 / mm^2), far above the old
  // absolute threshold 1e-14 / mm, which made the lens focal with EFL of about -8e10 mm.
  // The relative threshold is about 1.3e-9 / mm here.
  System s = base_system();
  add(s, lens("L", 10.0, 1.5, 0.003, -0.003, 0.018));
  const FirstOrder fo = first_order_of(s);
  REQUIRE(std::abs(fo.power) > 1e-12);  // the rounding is really there
  REQUIRE_FALSE(fo.efl.has_value());
  REQUIRE_FALSE(fo.rear_focal_z.has_value());
}

TEST_CASE("afocal threshold: a very weak lens stays focal", "[paraxial][afocal]") {
  // Plano-convex, R = 5e15 mm, n = 1.5: Phi = 0.5 / R = 1e-16 / mm, EFL 1e16 mm. All terms are
  // that small, so the relative threshold is about 1e-30 / mm; the old absolute 1e-14 / mm
  // made the lens afocal. One division and one product: relative 1e-12.
  System s = base_system();
  add(s, lens("L", 10.0, 1.5, 5e15, std::nullopt, 2.0));
  const FirstOrder fo = first_order_of(s);
  REQUIRE(fo.efl.has_value());
  require_rel(fo.efl, 1e16);
}

TEST_CASE("afocal threshold: a nearly afocal lens just above the threshold stays focal",
          "[paraxial][afocal]") {
  // r = 1 mm at z = 10 mm, d = 6 - 6e-12 mm: Phi = 1e-12 / mm, about 50 times the threshold
  // (1.9e-14 / mm), so the threshold is not too generous. Tolerance fixed beforehand: the
  // computed Phi is off by at most 8 N u S = 9.5e-15 / mm, relative 1e-2 of Phi.
  System s = base_system();
  const double d = 6.0 - 6.0e-12;
  add(s, lens("L", 10.0, 1.5, 1.0, -1.0, d));
  const FirstOrder fo = first_order_of(s);
  REQUIRE(fo.efl.has_value());
  const double phi = (6.0 - d) / 6.0;
  REQUIRE(std::abs(fo.power - phi) <= 1e-2 * phi);
}

TEST_CASE("afocal threshold: depends on the absolute position of the surfaces",
          "[paraxial][afocal]") {
  // S grows with |z|, because the coordinates carry a rounding of u |z| each. The same lens
  // (r = 1 mm, Phi = 5e-13 / mm) is focal at z = 10 mm (threshold 1.9e-14 / mm) and afocal at
  // z = 1000 mm (threshold 1.2e-12 / mm). The ordinary cases do not change: Phi = 1e-10 / mm
  // is focal and the afocal lens r = 10 mm, d = 60 mm is afocal at both positions.
  const auto first_order_at = [](double z, double r, double d) {
    System s = base_system();
    add(s, lens("L", z, 1.5, r, -r, d));
    return first_order_of(s);
  };
  const double weak = 6.0 - 3.0e-12;  // Phi = 5e-13 / mm
  REQUIRE(first_order_at(10.0, 1.0, weak).efl.has_value());
  REQUIRE_FALSE(first_order_at(1000.0, 1.0, weak).efl.has_value());
  const double clear = 6.0 - 6.0e-10;  // Phi = 1e-10 / mm
  for (const double z : {10.0, 1000.0}) {
    INFO("z = " << z);
    REQUIRE(first_order_at(z, 1.0, clear).efl.has_value());
    REQUIRE_FALSE(first_order_at(z, 10.0, 60.0).efl.has_value());
  }
}

// ----------------------------- infinity thresholds of pupils and image (#35, rest of B9) -----
// first_order decides three cases "at infinity" from one matrix element each: the entrance pupil
// (front.a), the exit pupil (back.d) and the image of a finite object (nu_out). An element that
// is zero in exact arithmetic comes out as a rounding remainder of about N u S, which would put
// the pupil or the image some 1e16 mm away instead of at infinity. The thresholds are relative,
// as for the afocal test: |element| <= 16 N u S with S the same y-nu trace in magnitudes.
//
// The test lens: plano-convex, n = 1.5, |R| = 64 mm, d = 3 mm, in VACUUM, f = 128 mm. With the
// curved side first the rear focal point lies 126 mm behind S2 (z = 129, telecentric_singlet);
// with the curved side last the front focal point lies 126 mm before S1 (z = -126): a stop or
// an object there sees the other side at infinity. Moving it by k ulp (|k| <= 16, about 4e-13
// mm) changes the element by at most 16 ulp(129) / 128 = 3.5e-15, below the threshold of about
// 1e-14 (derivation in paraxial.cpp); a shift of 1e-6 mm (element 8e-9) stays finite.

namespace {

/// x moved by k units in the last place (k < 0: towards -infinity).
double ulps_from(double x, int k) {
  const double to = (k < 0 ? -1.0 : 1.0) * std::numeric_limits<double>::infinity();
  for (int i = 0; i < std::abs(k); ++i) x = std::nextafter(x, to);
  return x;
}

}  // namespace

TEST_CASE("infinity thresholds: a stop in the rear focal plane up to rounding is telecentric",
          "[paraxial][b9]") {
  for (int k = -16; k <= 16; ++k) {
    INFO("stop at 129 + " << k << " ulp");
    System s = base_system();
    add(s, lens("L", 0.0, 1.5, 64.0, std::nullopt, 3.0));
    add(s, stop("STO", ulps_from(129.0, k), 5.0));
    const FirstOrder fo = first_order_of(s);
    REQUIRE(fo.entrance_pupil.has_value());
    CHECK_FALSE(fo.entrance_pupil->z.has_value());  // entrance pupil at infinity
    CHECK_FALSE(fo.angular_magnification.has_value());
  }
  // Guard: 1e-6 mm away the entrance pupil is finite (about 1.6e10 mm away).
  System s = base_system();
  add(s, lens("L", 0.0, 1.5, 64.0, std::nullopt, 3.0));
  add(s, stop("STO", 129.0 + 1e-6, 5.0));
  const FirstOrder fo = first_order_of(s);
  REQUIRE(fo.entrance_pupil.has_value());
  CHECK(fo.entrance_pupil->z.has_value());
}

TEST_CASE("infinity thresholds: a stop in the front focal plane up to rounding is telecentric",
          "[paraxial][b9]") {
  for (int k = -16; k <= 16; ++k) {
    INFO("stop at -126 + " << k << " ulp");
    System s = base_system();
    add(s, stop("STO", ulps_from(-126.0, k), 5.0));
    add(s, lens("L", 0.0, 1.5, std::nullopt, -64.0, 3.0));
    const FirstOrder fo = first_order_of(s);
    REQUIRE(fo.exit_pupil.has_value());
    CHECK_FALSE(fo.exit_pupil->z.has_value());  // exit pupil at infinity
    CHECK_FALSE(fo.exit_pupil->diameter.has_value());
  }
  System s = base_system();
  add(s, stop("STO", -126.0 + 1e-6, 5.0));
  add(s, lens("L", 0.0, 1.5, std::nullopt, -64.0, 3.0));
  const FirstOrder fo = first_order_of(s);
  REQUIRE(fo.exit_pupil.has_value());
  CHECK(fo.exit_pupil->z.has_value());
}

TEST_CASE("infinity thresholds: an object in the front focal plane up to rounding has its image "
          "at infinity",
          "[paraxial][b9]") {
  for (int k = -16; k <= 16; ++k) {
    INFO("object distance 126 + " << k << " ulp");
    System s = base_system();
    s.object = {false, Param(ulps_from(126.0, k))};
    add(s, lens("L", 0.0, 1.5, std::nullopt, -64.0, 3.0));
    const FirstOrder fo = first_order_of(s);
    CHECK_FALSE(fo.image_z.has_value());
    CHECK_FALSE(fo.lateral_magnification.has_value());
  }
  System s = base_system();
  s.object = {false, Param(126.0 + 1e-6)};
  add(s, lens("L", 0.0, 1.5, std::nullopt, -64.0, 3.0));
  const FirstOrder fo = first_order_of(s);
  CHECK(fo.image_z.has_value());
  CHECK(fo.lateral_magnification.has_value());
}
