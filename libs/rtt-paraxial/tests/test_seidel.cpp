#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/paraxial/seidel.hpp"

// Sources of the formulas (docs/quellen.md): J. Sasian, OPTI 517 lecture notes, L4 "Seidel
// aberration coefficients", p. 23 (Seidel sums and normalisation) and p. 24 (A, A_bar, H, P);
// OPTI 518 lecture notes, L14 "Aspheric surfaces and stop shifting", p. 17 (aspheric cap);
// OPTI 518 L6 "Chromatic aberrations", p. 4 and p. 14 (chromatic terms). Conventions in
// rtt/paraxial/seidel.hpp.
//
// Reference values marked "Python" come from an independent script (own y-nu trace and the
// formulas above, written for #30 without the C++ code); the script also checked all five sums
// against a real ray trace at small aperture and field (agreement to 1e-3 relative, the size of
// the higher-order terms). The script is described in the PR of #30.

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
using rtt::paraxial::ChromaticPair;
using rtt::paraxial::ParaxialError;
using rtt::paraxial::seidel;
using rtt::paraxial::Seidel;
using rtt::paraxial::SeidelSurfaceInput;
using rtt::paraxial::SeidelTerms;
using rtt::paraxial::surface_seidel;

namespace {

constexpr double kRel = 1e-12;

// ------------------------------------------------------------------ builders -----

System base_system() {
  System s;
  s.name = "seidel test";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(20.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  s.paths = {{"main", true, {}}};
  return s;
}

Surface plane_surface(const std::string& id, double z) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z);
  return s;
}

Element stop(double z, double radius) {
  Surface s = plane_surface("STO", 0.0);
  s.aperture = rtt::model::CircularAperture{radius, 0.0};
  return Element{"stop", ElementKind::Stop, Pose::along_z(z), std::nullopt, {s}};
}

Element mirror(double z, const rtt::model::BaseShape& base) {
  Surface s = plane_surface("M.S", 0.0);
  s.shape.base = base;
  s.interaction = rtt::model::IdealMirror{};
  return Element{"M", ElementKind::Mirror, Pose::along_z(z), std::nullopt, {s}};
}

void add(System& s, Element e) {
  s.root.children.push_back({std::move(e)});
}

CompiledSystem compile(const System& s) {
  const MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

System load(const std::string& name) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + name);
}

double tan_deg(double deg) {
  return std::tan(deg * std::numbers::pi / 180.0);
}

void require_sums(const SeidelTerms& t, const std::vector<double>& expected, double rel = kRel) {
  REQUIRE_THAT(t.s1, WithinRel(expected[0], rel));
  REQUIRE_THAT(t.s2, WithinRel(expected[1], rel));
  REQUIRE_THAT(t.s3, WithinRel(expected[2], rel));
  REQUIRE_THAT(t.s4, WithinRel(expected[3], rel));
  REQUIRE_THAT(t.s5, WithinRel(expected[4], rel));
}

/// The system values are the sums of the surface contributions, and the Lagrange invariant
/// n' (u_bar' y - u' y_bar) after every event equals the object-space value.
void check_sum_and_lagrange(const Seidel& r) {
  SeidelTerms sum;
  for (const auto& s : r.surfaces) {
    sum.s1 += s.terms.s1;
    sum.s2 += s.terms.s2;
    sum.s3 += s.terms.s3;
    sum.s4 += s.terms.s4;
    sum.s5 += s.terms.s5;
    sum.c_l += s.terms.c_l;
    sum.c_t += s.terms.c_t;
    REQUIRE_THAT(s.lagrange, WithinRel(r.lagrange, kRel));
  }
  const double scale = std::abs(r.sum.s1) + std::abs(r.sum.s2) + std::abs(r.sum.s3) +
                       std::abs(r.sum.s4) + std::abs(r.sum.s5);
  REQUIRE_THAT(sum.s1, WithinAbs(r.sum.s1, 1e-15 * scale));
  REQUIRE_THAT(sum.s2, WithinAbs(r.sum.s2, 1e-15 * scale));
  REQUIRE_THAT(sum.s3, WithinAbs(r.sum.s3, 1e-15 * scale));
  REQUIRE_THAT(sum.s4, WithinAbs(r.sum.s4, 1e-15 * scale));
  REQUIRE_THAT(sum.s5, WithinAbs(r.sum.s5, 1e-15 * scale));
  REQUIRE(sum.c_l == r.sum.c_l);
  REQUIRE(sum.c_t == r.sum.c_t);
}

}  // namespace

// ------------------------------------------------------------------- tests -----

TEST_CASE("spherical mirror with the stop at the centre of curvature: only S_I and S_IV",
          "[seidel]") {
  // Concave mirror R = -200 mm (c = -1/200 in global orientation, light from -z), stop at the
  // centre of curvature z = -200, EPD 2 r = 20 mm, maximum field theta = 3 deg, t = tan theta.
  // By hand with the formulas of Sasian, OPTI 517 L4 p. 23/24 (n = 1, n' = -1):
  //   marginal ray y = r, u = 0: A = n (u + y c) = r c; n'u' = nu - y c (n' - n) = 2 r c, so
  //   u' = -2 r c and Delta(u/n) = u'/n' - u/n = 2 r c; S_I = -A^2 y Delta(u/n) = -2 r^4 c^3
  //   = 2 r^4 / |R|^3.
  //   chief ray from the stop centre with slope t: y_bar = -t/c at the mirror, so
  //   A_bar = n (t + y_bar c) = 0 and S_II = S_III = 0 (factor A_bar), S_V = 0 (factor A_bar).
  //   H = n (u_bar y - u y_bar) = r t at the stop; P = c (1/n' - 1/n) = -2 c;
  //   S_IV = -H^2 P = 2 r^2 t^2 c = -2 r^2 t^2 / |R|.
  const double r = 10.0;
  const double radius = 200.0;
  const double t = tan_deg(3.0);
  System s = base_system();
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 3.0, 1.0}}};
  add(s, stop(-radius, r));
  add(s, mirror(0.0, rtt::model::Conic{Param(-radius), Param(0.0)}));
  const Seidel res = seidel(compile(s), PathId{0}, 0);

  const double s1 = 2.0 * std::pow(r, 4) / std::pow(radius, 3);
  const double s4 = -2.0 * r * r * t * t / radius;
  REQUIRE_THAT(res.sum.s1, WithinRel(s1, kRel));
  REQUIRE_THAT(res.sum.s4, WithinRel(s4, kRel));
  REQUIRE_THAT(res.sum.s2, WithinAbs(0.0, 1e-15 * s1));
  REQUIRE_THAT(res.sum.s3, WithinAbs(0.0, 1e-15 * s1));
  REQUIRE_THAT(res.sum.s5, WithinAbs(0.0, 1e-15 * s1));
  REQUIRE_THAT(res.lagrange, WithinRel(r * t, kRel));
  REQUIRE(res.surfaces.size() == 2);
  REQUIRE_THAT(res.surfaces[1].a_bar, WithinAbs(0.0, 1e-15));
  REQUIRE(res.sum.c_l == 0.0);
  REQUIRE(res.sum.c_t == 0.0);
  check_sum_and_lagrange(res);
}

TEST_CASE("paraboloid on axis with the stop at the mirror: S_I = 0", "[seidel]") {
  // The conic term cancels the spherical aberration of the base sphere (Sasian, OPTI 518 L14
  // p. 17: a = -eps^2 c^3 y^4 Delta(n) with eps^2 = -K, i.e. a = K c^3 y^4 Delta(n)):
  // sphere S_I = -2 r^4 c^3 (case above), a = (-1) c^3 r^4 (-2) = 2 r^4 c^3.
  // A4 = -c^3/8 on a sphere has the same fourth-order sag as the paraboloid
  // (z = c r^2/2 + c^3 r^4/8 + ... for the sphere, Sasian L14 p. 6 with K = 0), so it cancels
  // S_I as well; with the mirror turned by 180 deg about x the local and global curvature and
  // A4 change sign together.
  const double r = 30.0;
  const double radius = 200.0;
  const double sphere_s1 = 2.0 * std::pow(r, 4) / std::pow(radius, 3);
  const auto s1_of = [&](const rtt::model::BaseShape& base, bool turned) {
    System s = base_system();
    s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(2.0 * r)};
    add(s, stop(0.0, r));
    Element m = mirror(0.0, base);
    if (turned) m.pose.rotation_deg[0] = Param(180.0);
    add(s, m);
    const Seidel res = seidel(compile(s), PathId{0}, 0);
    check_sum_and_lagrange(res);
    return res.sum.s1;
  };
  const double c_local = -1.0 / radius;  // concave for light from -z
  REQUIRE_THAT(s1_of(rtt::model::Conic{Param(-radius), Param(0.0)}, false),
               WithinRel(sphere_s1, kRel));
  REQUIRE_THAT(s1_of(rtt::model::Conic{Param(-radius), Param(-1.0)}, false),
               WithinAbs(0.0, 1e-12 * sphere_s1));
  REQUIRE_THAT(s1_of(rtt::model::EvenAsphere{Param(-radius), Param(-1.0), {}}, false),
               WithinAbs(0.0, 1e-12 * sphere_s1));
  REQUIRE_THAT(s1_of(
                   rtt::model::EvenAsphere{
                       Param(-radius), Param(0.0), {Param(-c_local * c_local * c_local / 8.0)}},
                   false),
               WithinAbs(0.0, 1e-12 * sphere_s1));
  // Turned mirror: local R = +200 is concave for light from -z.
  const double c_turned = 1.0 / radius;
  REQUIRE_THAT(s1_of(rtt::model::Conic{Param(radius), Param(0.0)}, true),
               WithinRel(sphere_s1, kRel));
  REQUIRE_THAT(s1_of(
                   rtt::model::EvenAsphere{
                       Param(radius), Param(0.0), {Param(-c_turned * c_turned * c_turned / 8.0)}},
                   true),
               WithinAbs(0.0, 1e-12 * sphere_s1));
}

TEST_CASE("reference singlet: Seidel sums by hand, Petzval sum and Lagrange invariant",
          "[seidel]") {
  // tests/reference/m1/singlet_const.rtt.json: stop (r = 10) at z = 0, plano-convex lens
  // CONST:1.5168, R1 = 51.68 at z = 5, plane at z = 9; fields 0, 3.5 and 5 deg (maximum 5 deg).
  const CompiledSystem cs = compile(load("m1/singlet_const.rtt.json"));
  const double n = 1.5168;
  const double t = tan_deg(5.0);
  // Python (see top of file).
  const std::vector<double> expected = {0.02177852531749783, 0.0034952015762708061,
                                        0.0078042652442534685, 0.0050463253201162607,
                                        0.0013265226745504148};
  for (std::uint16_t wl = 0; wl < 3; ++wl) {  // CONST: the same at every wavelength
    const Seidel res = seidel(cs, PathId{0}, wl);
    require_sums(res.sum, expected);
    // Lagrange invariant H = n r t (object space at the stop, y = r, u = 0, y_bar = 0).
    REQUIRE_THAT(res.lagrange, WithinRel(10.0 * t, kRel));
    // Petzval sum of the single powered surface: S_IV = -H^2 c (1/n - 1) = H^2 phi / n
    // (Sasian L4 p. 23/24) with phi = (n - 1) c.
    const double h = 10.0 * t;
    REQUIRE_THAT(res.sum.s4, WithinRel(h * h * (n - 1.0) / 51.68 / n, kRel));
    REQUIRE(res.surfaces.size() == 4);  // stop, two lens surfaces, detector
    REQUIRE(res.surfaces[0].terms.s1 == 0.0);
    REQUIRE(res.surfaces[3].terms.s4 == 0.0);
    REQUIRE(res.marginal.z == 0.0);
    REQUIRE(res.marginal.y == 10.0);
    REQUIRE(res.marginal.u == 0.0);
    REQUIRE_THAT(res.chief.u, WithinRel(t, kRel));
    check_sum_and_lagrange(res);
  }

  SECTION("chromatic terms vanish for constant indices") {
    const Seidel res = seidel(cs, PathId{0}, 1, ChromaticPair{0, 2});
    REQUIRE(res.chromatic.has_value());
    REQUIRE(res.sum.c_l == 0.0);
    REQUIRE(res.sum.c_t == 0.0);
    require_sums(res.sum, expected);
  }
  SECTION("the field given as paraxial image height gives the same sums") {
    // EFL = R1 / (n - 1) = 100 mm, object at infinity: the 5 deg chief ray reaches the paraxial
    // image at h' = EFL tan 5 deg. The sign of a field value does not matter (maximum of the
    // radial value).
    System s = load("m1/singlet_const.rtt.json");
    s.fields = {rtt::model::FieldType::ParaxialImageHeight,
                {{0.0, 0.0, 1.0}, {0.0, -100.0 * t, 1.0}, {0.0, 50.0 * t, 1.0}}};
    const Seidel res = seidel(compile(s), PathId{0}, 1);
    require_sums(res.sum, expected, 1e-11);
  }
}

TEST_CASE("two singlets with the stop between: Seidel sums by hand", "[seidel]") {
  // tests/reference/m1/two_lenses_stop_between.rtt.json: stop size 4 mm between the lenses;
  // fields 0, (0, 5) and (-3, 4) deg. Largest radial value tan theta = tan 5 deg
  // (hypot(tan 3, tan 4) = 0.08739 < tan 5 = 0.08749). Entrance pupil (Python): z = 19.196...,
  // r = 4.511...
  const Seidel res = seidel(compile(load("m1/two_lenses_stop_between.rtt.json")), PathId{0}, 0);
  // Python (see top of file).
  require_sums(res.sum, {0.0093739209101174215, -0.0062266351246948425, 0.004220499333761403,
                         0.0029484609542615144, -0.00061714305623997135});
  REQUIRE_THAT(res.lagrange, WithinRel(0.39467248531034937, kRel));
  REQUIRE_THAT(res.marginal.z, WithinRel(19.196437896201594, kRel));
  REQUIRE_THAT(res.marginal.y, WithinRel(4.5111271495581011, kRel));
  REQUIRE_THAT(res.chief.u, WithinRel(tan_deg(5.0), kRel));
  check_sum_and_lagrange(res);
}

TEST_CASE("aspheric lens with a finite object and an object-height field", "[seidel]") {
  // Stop (EPD 2) at z = 0, object at z = -300, object height 3 mm. Lens n = 1.6 at z = 20,
  // thickness 5: S1 even asphere R = 40, K = -0.7, A4 = -3e-6; S2 even asphere R = -120,
  // K = 0, A4 = 2e-6. Marginal ray from (z, y) = (-300, 0) with u = 1/300, chief ray from
  // (-300, 3) with u = -3/300; H = n (u_bar y - u y_bar) = -0.01 mm.
  System s = base_system();
  s.object = {false, Param(300.0)};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(2.0)};
  s.fields = {rtt::model::FieldType::ObjectHeight, {{0.0, 0.0, 1.0}, {0.0, 3.0, 1.0}}};
  add(s, stop(0.0, 1.0));
  Surface s1 = plane_surface("L.S1", 0.0);
  s1.shape.base = rtt::model::EvenAsphere{Param(40.0), Param(-0.7), {Param(-3e-6)}};
  Surface s2 = plane_surface("L.S2", 5.0);
  s2.shape.base = rtt::model::EvenAsphere{Param(-120.0), Param(0.0), {Param(2e-6)}};
  add(s, Element{"L", ElementKind::Lens, Pose::along_z(20.0), "CONST:1.6", {s1, s2}});
  const Seidel res = seidel(compile(s), PathId{0}, 0);

  // Python (see top of file), per surface and in total.
  REQUIRE(res.surfaces.size() == 3);
  require_sums(res.surfaces[1].terms,
               {-1.9186758518518526e-05, 1.1131422222222229e-06, 1.0334733333333334e-06,
                9.3750000000000024e-07, -1.2836200000000002e-06});
  require_sums(res.surfaces[2].terms,
               {-1.5086315192713673e-06, 3.7144022868240428e-06, -2.8287873592423827e-07,
                3.1250000000000003e-07, 1.8867379661277488e-07});
  require_sums(res.sum, {-2.0695390037789892e-05, 4.8275445090462653e-06, 7.5059459740909515e-07,
                         1.2500000000000003e-06, -1.0949462033872252e-06});
  REQUIRE_THAT(res.lagrange, WithinRel(-0.01, kRel));
  REQUIRE_THAT(res.chief.y, WithinRel(3.0, kRel));
  REQUIRE_THAT(res.chief.u, WithinRel(-0.01, kRel));
  check_sum_and_lagrange(res);
}

TEST_CASE("chromatic terms of a thin lens: C_L = y^2 phi / V, C_T = y y_bar phi / V", "[seidel]") {
  // Thin lens in air from two surfaces at the same height y (Sasian, OPTI 518 L6 p. 4 and
  // p. 14: C_L = sum A y Delta(dn/n), C_T = sum A_bar y Delta(dn/n), i.e. 2 d_lambda W020 and
  // d_lambda W111). By hand: Delta(dn/n) is +dn/n at the first and -dn/n at the second surface,
  // so C_L = y dn/n (A_1 - A_2) with A_1 = u + y c1 (air) and A_2 = u'' + y c2 (air after the
  // lens), u'' = u - y phi: A_1 - A_2 = y (c1 - c2) + y phi = y n (c1 - c2). Hence
  // C_L = y^2 dn (c1 - c2) = y^2 phi / V with V = (n - 1) / dn; likewise
  // C_T = y y_bar phi / V (A_bar_1 - A_bar_2 = y_bar n (c1 - c2)).
  const double n = 1.5;
  const double dn = 0.01;  // V = 50
  const double c1 = 1.0 / 50.0;
  const double c2 = -1.0 / 50.0;
  const double phi = (n - 1.0) * (c1 - c2);
  const double v = (n - 1.0) / dn;
  const double y = 5.0;
  const double u = 0.02;  // finite object; the result does not depend on u
  const double u_bar = 0.1;

  for (const double y_bar : {0.0, 2.0}) {
    SeidelSurfaceInput first;
    first.c = c1;
    first.n = 1.0;
    first.n_after = n;
    first.y = y;
    first.u = u;
    first.y_bar = y_bar;
    first.u_bar = u_bar;
    first.dn = 0.0;
    first.dn_after = dn;
    SeidelSurfaceInput second = first;
    second.c = c2;
    second.n = n;
    second.n_after = 1.0;
    second.u = (u - y * c1 * (n - 1.0)) / n;  // n'u' = n u - y c (n' - n)
    second.u_bar = (u_bar - y_bar * c1 * (n - 1.0)) / n;
    second.dn = dn;
    second.dn_after = 0.0;
    const SeidelTerms a = surface_seidel(first);
    const SeidelTerms b = surface_seidel(second);
    REQUIRE_THAT(a.c_l + b.c_l, WithinRel(y * y * phi / v, kRel));
    if (y_bar == 0.0) {
      REQUIRE_THAT(a.c_t + b.c_t, WithinAbs(0.0, 1e-17));
    } else {
      REQUIRE_THAT(a.c_t + b.c_t, WithinRel(y * y_bar * phi / v, kRel));
    }
  }

  SECTION("a mirror adds no colour") {
    // Mirror in glass: n' = -n and dn' = -dn, so Delta(dn/n) = 0.
    SeidelSurfaceInput m;
    m.c = -0.01;
    m.n = n;
    m.n_after = -n;
    m.y = y;
    m.u = u;
    m.y_bar = 1.0;
    m.u_bar = u_bar;
    m.dn = dn;
    m.dn_after = -dn;
    const SeidelTerms t = surface_seidel(m);
    REQUIRE(t.c_l == 0.0);
    REQUIRE(t.c_t == 0.0);
  }
}

TEST_CASE("chromatic terms of a thin lens in a system with a dispersive test material",
          "[seidel]") {
  // Same hand formulas as above (C_L = y^2 phi / V, C_T = y y_bar phi / V), now through
  // seidel() with the media of the compiled system. Test material with
  // n(lambda) = 1.5 - 0.01 (lambda - 0.5876) / (0.6563 - 0.4861): n_d = 1.5,
  // dn = n_F - n_C = 0.01, V = (n_d - 1) / dn = 50. Thin lens R1 = 50, R2 = -50 (two surfaces
  // at z = 0): phi = (n_d - 1)(c1 - c2) = 0.02 / mm. EPD 10 (y = 5 at the lens), field 5 deg.
  class LinearGlass final : public rtt::material::Material {
   public:
    [[nodiscard]] rtt::math::Complex index(double wavelength_um,
                                           double /*temperature_c*/,
                                           double /*pressure_atm*/) const override {
      return {1.5 - 0.01 * (wavelength_um - 0.5876) / (0.6563 - 0.4861), 0.0};
    }
  };
  const double y = 5.0;
  const double phi = 0.5 * (1.0 / 50.0 + 1.0 / 50.0);
  const double v = 50.0;
  const double t = tan_deg(5.0);
  const auto run = [&](double z_stop, ChromaticPair pair) {
    System s = base_system();
    s.wavelengths = {{0.4861, 1.0, false}, {0.5876, 1.0, true}, {0.6563, 1.0, false}};
    s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(2.0 * y)};
    s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 5.0, 1.0}}};
    add(s, stop(z_stop, y));
    Surface s1 = plane_surface("L.S1", 0.0);
    s1.shape.base = rtt::model::Conic{Param(50.0), Param(0.0)};
    Surface s2 = plane_surface("L.S2", 0.0);
    s2.shape.base = rtt::model::Conic{Param(-50.0), Param(0.0)};
    add(s, Element{"L", ElementKind::Lens, Pose::along_z(0.0), "TEST:LINEAR", {s1, s2}});
    MaterialLibrary lib;
    lib.add("TEST:LINEAR", std::make_shared<const LinearGlass>());
    const CompiledSystem cs = rtt::compile::compile(s, lib);
    const Seidel res = seidel(cs, PathId{0}, 1, pair);
    check_sum_and_lagrange(res);
    return res;
  };
  SECTION("stop at the lens: no lateral colour") {
    const Seidel res = run(0.0, ChromaticPair{0, 2});
    REQUIRE_THAT(res.sum.c_l, WithinRel(y * y * phi / v, kRel));
    REQUIRE_THAT(res.sum.c_t, WithinAbs(0.0, 1e-17));
  }
  SECTION("stop 20 mm in front: y_bar = 20 tan theta at the lens") {
    const Seidel res = run(-20.0, ChromaticPair{0, 2});
    REQUIRE_THAT(res.sum.c_l, WithinRel(y * y * phi / v, kRel));
    REQUIRE_THAT(res.sum.c_t, WithinRel(y * 20.0 * t * phi / v, kRel));
    // Reversed pair: dn = n_C - n_F changes the sign of both terms.
    const Seidel rev = run(-20.0, ChromaticPair{2, 0});
    REQUIRE_THAT(rev.sum.c_l, WithinRel(-y * y * phi / v, kRel));
    REQUIRE_THAT(rev.sum.c_t, WithinRel(-y * 20.0 * t * phi / v, kRel));
  }
}

TEST_CASE("Seidel input errors", "[seidel]") {
  SECTION("not rotationally symmetric") {
    System s = load("m1/singlet_const.rtt.json");
    std::get<Element>(s.root.children[1].value).pose.rotation_deg[0] = Param(1e-3);
    REQUIRE_THROWS_AS(seidel(compile(s), PathId{0}, 0), ParaxialError);
  }
  SECTION("no stop on the path") {
    const CompiledSystem cs = compile(load("m1/paraboloid_mirror.rtt.json"));
    REQUIRE_THROWS_AS(seidel(cs, PathId{0}, 0), ParaxialError);
  }
  SECTION("object height with the object at infinity") {
    System s = load("m1/singlet_const.rtt.json");
    s.fields = {rtt::model::FieldType::ObjectHeight, {{0.0, 1.0, 1.0}}};
    REQUIRE_THROWS_AS(seidel(compile(s), PathId{0}, 0), ParaxialError);
  }
  SECTION("wavelength indices") {
    const CompiledSystem cs = compile(load("m1/singlet_const.rtt.json"));
    REQUIRE_THROWS_AS(seidel(cs, PathId{0}, 3), ParaxialError);
    REQUIRE_THROWS_AS(seidel(cs, PathId{0}, 0, ChromaticPair{0, 3}), ParaxialError);
    REQUIRE_THROWS_AS(seidel(cs, PathId{1}, 0), ParaxialError);
  }
}
