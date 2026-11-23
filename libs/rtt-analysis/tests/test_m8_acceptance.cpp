// M8 acceptance (#324): the surfaces of S1 to S3 and the apertures of F1 end to end, read
// through rtt-io from tests/reference/m8 and examples/, compiled, aimed or started by hand, traced
// and analysed. The cases and expected values are those of plan v1 in #324, agreed with the
// maintainer (2026-11-22) before the code; a red case is a finding, not a reason to change the
// test (AGENTS.md, rule 2).
//
// Every expected value is computed here from the source of its formula, not with the code under
// test: Snell's law, an own y-nu trace (Greivenkamp, OPTI-201/202, Sec. 9, n'u' = nu - y phi,
// phi = (n' - n) c), the radial polynomials of Niu/Tian Eq. (4), the area of a circular segment
// and of an ellipse (MathWorld). Values of libm functions on constants are computed at run time
// from volatile inputs, never folded (AGENTS.md, #218). The abbreviations [F07], [F10b] and [F11]
// for the Forbes papers and [UG] for the OpticStudio User Guide are defined in ADR 0043 and
// docs/quellen.md.
//
// Bounds:
// - traced positions and OPL: the Newton intersection stops at |z - sag| <= kNewtonTolerance =
//   1e-12 mm (rtt/geom/intersect.hpp) or at a step of rounding size; a hit point is thus fixed to
//   about 1e-12 mm along z, and the propagation over <= 200 mm with a direction that carries a
//   few ulp (1e-15) adds 2e-13 mm. Bound kPos = 1e-11 mm (a factor 5 over the sum of two hits);
// - directions: a few ulp of the refraction and reflection formulas, bound kDir = 1e-14;
// - first-order values: about 20 roundings of the y-nu trace (both sides) on a well conditioned
//   singlet or mirror, 20 eps = 4.4e-15 relative; bound kRel = 1e-12 relative, as in the plan;
// - transmitted fractions: the hits are binomial, sigma = sqrt(p (1 - p) / N), bound 5 sigma as
//   in rtt-trace/tests/test_polygon_aperture.cpp (#320).

#include <tbb/global_control.h>

#include <algorithm>
#include <array>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

namespace fs = std::filesystem;

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::Element;
using rtt::model::System;
using rtt::paraxial::ParaxialError;
using rtt::paraxial::Section;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kPos = 1e-11;  // mm, see the head of the file
constexpr double kDir = 1e-14;
constexpr double kRel = 1e-12;

System load_reference(const std::string& name) {
  return rtt::io::load_system(fs::path(RTT_REFERENCE_DIR) / "m8" / name);
}

System load_example(const std::string& relative) {
  return rtt::io::load_system(fs::path(RTT_REFERENCE_DIR) / ".." / ".." / "examples" / relative);
}

CompiledSystem compiled(const System& s) {
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

/// Element `i` of the root assembly.
Element& element(System& s, std::size_t i) {
  return std::get<Element>(s.root.children[i].value);
}

/// Run-time value of a double, not folded by the compiler (AGENTS.md, #218).
double runtime(double v) {
  volatile const double x = v;
  return x;
}

double deg(double d) {
  return runtime(d) * kPi / 180.0;
}

struct Start {
  double x, y, z, dx, dy, dz;
};

/// Rays started by hand at the reference wavelength.
RayBatch rays_from(const CompiledSystem& cs, const std::vector<Start>& starts) {
  RayBatch b(starts.size());
  for (std::size_t i = 0; i < starts.size(); ++i) {
    b.pos_x()[i] = starts[i].x;
    b.pos_y()[i] = starts[i].y;
    b.pos_z()[i] = starts[i].z;
    b.dir_x()[i] = starts[i].dx;
    b.dir_y()[i] = starts[i].dy;
    b.dir_z()[i] = starts[i].dz;
    b.wl()[i] = cs.reference_wavelength();
  }
  return b;
}

RayBatch traced(const CompiledSystem& cs, RayBatch b) {
  (void)rtt::trace::SequentialTracer{}.trace(cs, PathId{0}, b);
  return b;
}

std::uint32_t surface(const CompiledSystem& cs, const std::string& id) {
  const auto i = cs.find_surface(rtt::model::SurfaceId(id));
  REQUIRE(i.has_value());
  return i.value_or(std::numeric_limits<std::uint32_t>::max());
}

/// Paraxial focal data of a singlet in vacuum, curved front (vertex curvature c), plane back,
/// index n, thickness t, by the y-nu trace of Greivenkamp Sec. 9 (S. 9-2, 9-12): y1 = 1, u = 0;
/// n u1' = -y1 (n - 1) c; y2 = y1 + t u1'; u2' = n u1' (plane, n' = 1); Phi = -u2' / y1,
/// efl = 1 / Phi, bfl = -y2 / u2'.
struct Focal {
  double efl, bfl;
};
Focal singlet(double c, double n, double t) {
  const double u1 = -(n - 1.0) * c / n;
  const double y2 = 1.0 + t * u1;
  const double u2 = n * u1;
  return {-1.0 / u2, -y2 / u2};
}

/// Mirror in a medium of index n_m (the default environment is air), light along +z, vertex
/// curvature c: n = n_m, n' = -n_m (signed index after a reflection), phi = (n' - n) c =
/// -2 n_m c; n'u' = -y phi gives Phi = -n'u' / y = phi and efl = 1 / phi (Greivenkamp Sec. 9).
double mirror_efl(double c, double n_m) {
  return 1.0 / (-2.0 * n_m * c);
}

/// |n| of the image medium at the reference wavelength, read from the compiled media (not
/// through rtt-paraxial), as in test_m6_acceptance.cpp.
double image_index(const CompiledSystem& cs) {
  const auto& last = cs.path(PathId{0}).events.back();
  return std::abs(cs.media()[last.medium_after].index[cs.reference_wavelength()].real());
}

/// The value of a first-order result that must be defined (NaN fails every comparison).
double got(const std::optional<double>& v) {
  CHECK(v.has_value());
  return v.value_or(std::numeric_limits<double>::quiet_NaN());
}

bool rel_ok(double actual, double expected) {
  return std::abs(actual - expected) <= kRel * std::abs(expected);
}

/// Coefficient of rho^k in the radial polynomial R_n^m, Niu/Tian Eq. (4):
/// R_n^m = sum_s (-1)^s (n - s)! / (s! ((n + m)/2 - s)! ((n - m)/2 - s)!) rho^(n - 2s).
double radial_coefficient(int n, int m, int k) {
  const auto fact = [](int v) {
    double f = 1.0;
    for (int i = 2; i <= v; ++i) f *= static_cast<double>(i);
    return f;
  };
  const int s = (n - k) / 2;
  if ((n - k) % 2 != 0 || s < 0 || s > (n - m) / 2) return 0.0;
  const double sign = (s % 2 == 0) ? 1.0 : -1.0;
  return sign * fact(n - s) / (fact(s) * fact((n + m) / 2 - s) * fact((n - m) / 2 - s));
}

/// Fraction of rays with status Alive after the trace; checks that every other ray is
/// Vignetted.
double alive_fraction(const RayBatch& r) {
  std::size_t alive = 0;
  for (std::size_t i = 0; i < r.size(); ++i) {
    if (r.status()[i] == RayStatus::Alive) {
      ++alive;
    } else {
      CHECK(r.status()[i] == RayStatus::Vignetted);
    }
  }
  return static_cast<double>(alive) / static_cast<double>(r.size());
}

/// Area of the circular segment cut from a circle of radius R by a chord at distance d from
/// the centre: A = R^2 acos(d / R) - d sqrt(R^2 - d^2) (MathWorld "Circular Segment", Eq. (17)).
double segment_area(double r, double d) {
  return r * r * std::acos(d / r) - d * std::sqrt(r * r - d * d);
}

/// Area of the strip |y| <= w inside a circle of radius R: the circle without the two segments
/// at distance w (Eq. (17) above).
double strip_area(double r, double w) {
  return kPi * r * r - 2.0 * segment_area(r, w);
}

std::vector<double> columns(const RayBatch& r) {
  std::vector<double> v;
  for (std::size_t i = 0; i < r.size(); ++i) {
    for (const double x : {r.pos_x()[i], r.pos_y()[i], r.pos_z()[i], r.dir_x()[i], r.dir_y()[i],
                           r.dir_z()[i], r.opl()[i]}) {
      v.push_back(x);
    }
    v.push_back(static_cast<double>(r.status()[i]));
    v.push_back(static_cast<double>(r.last_surface()[i]));
  }
  return v;
}

bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (std::bit_cast<std::uint64_t>(a[i]) != std::bit_cast<std::uint64_t>(b[i])) return false;
  }
  return true;
}

// ---------------------------------------------------------------- systems of the cases -----

/// Case 1: own axis-parallel rays at heights h (in y and one in x) in front of the stop.
RayBatch axicon_rays(const CompiledSystem& cs) {
  return traced(cs, rays_from(cs, {{0.0, 1.0, -1.0, 0.0, 0.0, 1.0},
                                   {0.0, 2.5, -1.0, 0.0, 0.0, 1.0},
                                   {0.0, 4.0, -1.0, 0.0, 0.0, 1.0},
                                   {3.0, 0.0, -1.0, 0.0, 0.0, 1.0}}));
}

/// Case 2: the Q-bfs singlet with the stop opened to r = 13 mm, so that rays at the normalization
/// radius rho = 12.7 mm reach L1.S1 (the stop only clips; the lens is unchanged).
System qbfs_open_stop() {
  System s = load_example("forbes/qbfs_singlet.rtt.json");
  element(s, 0).surfaces[0].aperture = rtt::model::CircularAperture{13.0, 0.0};
  return s;
}

/// Case 7: the grid-sag mirror, its nominal design (the conic without the term) and the analytic
/// surface: the conic plus an XY polynomial a x^2 / R^2 - a y^2 / R^2, R = 22 mm.
constexpr double kGridA = 2e-5;  // mm, a of a (x^2 - y^2) / 22^2 (name of the reference file)
System grid_variant(bool analytic) {
  System s = load_reference("grid_sag_mirror.rtt.json");
  auto& terms = element(s, 1).surfaces[0].shape.terms;
  terms.clear();
  if (analytic) {
    terms.emplace_back(rtt::model::XyPolynomial{
        22.0, {0.0, 0.0, kGridA, 0.0, -kGridA}});  // x, y, x^2, x y, y^2 ([UG])
  }
  return s;
}

/// Points of the grid case: px, py in {-0.7, -0.35, 0, 0.35, 0.7} of the pupil radius 20 mm, so
/// |x|, |y| <= 14 mm < 18 mm: inside the inner cells (the nodes at +-24 mm are the rim), and
/// r <= 19.8 mm inside the stop (21 mm) and the mirror aperture (22 mm).
rtt::trace::PointListPupil grid_points() {
  rtt::trace::PointListPupil p;
  for (const double px : {-0.7, -0.35, 0.0, 0.35, 0.7}) {
    for (const double py : {-0.7, -0.35, 0.0, 0.35, 0.7}) p.points.push_back({px, py});
  }
  return p;
}

RayBatch aimed_random(const CompiledSystem& cs, std::size_t n) {
  const std::vector<std::uint16_t> fields{0};
  return traced(cs, rtt::trace::make_rays(cs, PathId{0}, fields, cs.reference_wavelength(),
                                          rtt::trace::RandomPupil{n, 324}));
}

}  // namespace

TEST_CASE("m8 case 1: an axicon deviates by asin(n sin a) - a and crosses the axis", "[m8]") {
  // examples/rotational/axicon_planoconvex: plane front at z = 5, cone z = r tan(alpha) with
  // alpha = -5 deg at z = 9 (ADR 0043, point 2), n = 1.5, vacuum, detector at z = 40. A ray
  // parallel to the axis at height h meets the mantle at z_e = 9 - h tan 5 deg; there the normal
  // is tilted by 5 deg, so the ray leaves at theta_t = asin(n sin 5 deg) to the normal (Snell)
  // and is turned towards the axis by delta = theta_t - 5 deg. It crosses the axis at
  // z_e + h / tan(delta) and is at height h - (40 - z_e) tan(delta) on the detector.
  const CompiledSystem cs = compiled(load_example("rotational/axicon_planoconvex.rtt.json"));
  const RayBatch r = axicon_rays(cs);
  const double a = deg(5.0);
  const double n = runtime(1.5);
  const double delta = std::asin(n * std::sin(a)) - a;
  REQUIRE(delta > 0.04);  // content guard: about 2.512 deg
  const std::array<double, 4> h{1.0, 2.5, 4.0, 3.0};
  for (std::size_t i = 0; i < r.size(); ++i) {
    INFO("ray " << i << ", h = " << h[i]);
    REQUIRE(r.status()[i] == RayStatus::Alive);
    CHECK(r.last_surface()[i] == surface(cs, "IMG"));
    const double ze = 9.0 - h[i] * std::tan(a);
    const double height = h[i] - (40.0 - ze) * std::tan(delta);
    const bool along_x = i == 3;
    const double px = along_x ? r.pos_x()[i] : r.pos_y()[i];
    const double dt = along_x ? r.dir_x()[i] : r.dir_y()[i];
    CHECK(std::abs(r.pos_z()[i] - 40.0) <= kPos);
    CHECK(std::abs(px - height) <= kPos);
    CHECK(std::abs(dt + std::sin(delta)) <= kDir);
    CHECK(std::abs(r.dir_z()[i] - std::cos(delta)) <= kDir);
    CHECK((along_x ? r.pos_y()[i] : r.pos_x()[i]) == 0.0);
    // OPL from z = -1: 6 mm in vacuum, n (z_e - 5) in glass, (40 - z_e) / cos(delta) behind.
    const double opl = 6.0 + n * (ze - 5.0) + (40.0 - ze) / std::cos(delta);
    CHECK(std::abs(r.opl()[i] - opl) <= kPos);
    // Crossing of the axis (the case in the name), from the traced ray: z_c = z - p d_z / d_t.
    // Errors kPos in z and p and kDir in d_z and d_t propagate to
    // kPos (1 + |d_z / d_t|) + |p| kDir (1 / |d_t| + |d_z| / d_t^2), about 5e-10 mm here.
    const double z_cross = r.pos_z()[i] - px * r.dir_z()[i] / dt;
    const double cross_bound =
        kPos * (1.0 + std::abs(r.dir_z()[i] / dt)) +
        std::abs(px) * kDir * (1.0 / std::abs(dt) + std::abs(r.dir_z()[i]) / (dt * dt));
    CHECK(std::abs(z_cross - (ze + h[i] / std::tan(delta))) <= cross_bound);
  }
  // No paraxial model at the apex (ADR 0043, F6).
  CHECK_THROWS_AS(rtt::paraxial::first_order(cs, PathId{0}, cs.reference_wavelength()),
                  ParaxialError);
}

TEST_CASE("m8 case 2: Q-bfs meets its best-fit conic at the rim and has c + 2 S0 / rho^2", "[m8]") {
  // [F11] Eq. (2.4), as QBfs in surface.hpp: z = c r^2 / (1 + phi) + u^2 (1 - u^2) / sigma
  // sum a_m Q_m(u^2), u = r / rho. At u = 1 the factor u^2 (1 - u^2) is 0: the surface is the
  // base conic there. A ray parallel to the axis at r = rho keeps x, y exactly, so its hit is
  // checked in z alone; bound 2 kNewtonTolerance plus a few ulp of the conic, 2e-12 mm.
  const System s = qbfs_open_stop();
  const CompiledSystem cs = compiled(s);
  const double c = 1.0 / 60.0;  // as in the file: 0.016666666666666666
  const double k = -0.3;
  const double rho = 12.7;
  const double diag = rho / std::sqrt(runtime(2.0));
  RayBatch b = rays_from(cs, {{0.0, rho, -1.0, 0.0, 0.0, 1.0},
                              {-rho, 0.0, -1.0, 0.0, 0.0, 1.0},
                              {diag, diag, -1.0, 0.0, 0.0, 1.0},
                              {0.0, 0.7 * rho, -1.0, 0.0, 0.0, 1.0}});
  rtt::trace::RayPaths paths;
  (void)rtt::trace::SequentialTracer{}.trace(cs, PathId{0}, b, paths);
  const std::uint32_t s1 = surface(cs, "L1.S1");
  std::size_t slot = 0;
  for (std::size_t e = 0; e < paths.event_surfaces.size(); ++e) {
    if (paths.event_surfaces[e] == s1) slot = e + 1;
  }
  REQUIRE(slot > 0);
  const auto conic_z = [&](double r) {
    return 5.0 + c * r * r / (1.0 + std::sqrt(1.0 - (1.0 + k) * c * c * r * r));
  };
  for (std::size_t i = 0; i < 3; ++i) {
    INFO("rim ray " << i);
    REQUIRE(b.status()[i] == RayStatus::Alive);
    const double x = paths.position[(i * paths.slots + slot) * 3 + 0];
    const double y = paths.position[(i * paths.slots + slot) * 3 + 1];
    const double z = paths.position[(i * paths.slots + slot) * 3 + 2];
    CHECK(std::abs(std::hypot(x, y) - rho) <= 4e-15 * rho);
    CHECK(std::abs(z - conic_z(rho)) <= 2e-12);
  }
  // Content guard: inside the rim the Q-bfs part is there (about 0.019 mm at u = 0.7).
  const double z_in = paths.position[(3 * paths.slots + slot) * 3 + 2];
  CHECK(std::abs(z_in - conic_z(0.7 * rho)) > 1e-3);

  // First order: vertex curvature c_v = c + 2 S_0 / rho^2 ([F10b], Eq. (3.16); ADR 0043, point
  // 3) with S_0 = a_0 Q_0(0) + a_1 Q_1(0) = a_0 + 13 a_1 / sqrt(19) ([F07], Eq. (14)), the
  // singlet with n = 1.5168, t = 4 mm and a plane back by the own y-nu trace.
  const double s0 = 0.1 + 13.0 * -0.02 / std::sqrt(runtime(19.0));
  const double cv = c + 2.0 * s0 / (rho * rho);
  const Focal f = singlet(cv, 1.5168, 4.0);
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, cs.reference_wavelength());
  REQUIRE(fo.efl.has_value());
  REQUIRE(fo.bfl.has_value());
  INFO("efl " << got(fo.efl) << " expected " << f.efl);
  CHECK(rel_ok(got(fo.efl), f.efl));
  CHECK(rel_ok(got(fo.bfl), f.bfl));
  CHECK(std::abs(got(fo.efl) - singlet(c, 1.5168, 4.0).efl) > 1.0);  // S_0 changes the focus
}

TEST_CASE("m8 case 3: odd asphere and Q-con follow the vertex curvature", "[m8]") {
  // Odd asphere ([UG] "Odd Asphere"): z = conic + sum A_i r^i; A_2 r^2 adds 2 A_2 to the vertex
  // curvature, A_3 r^3 and A_4 r^4 leave the first order (ADR 0043, point 3): c_v = 0.02 +
  // 2e-4 / mm. Q-con ([F07] Eqs. (2), (3)): the deviation carries u^4, c_v = c = 1/64 / mm.
  struct Singlet {
    const char* file;
    double cv;
  };
  for (const Singlet& sg : {Singlet{"rotational/odd_asphere_singlet.rtt.json", 0.02 + 2.0 * 1e-4},
                            Singlet{"forbes/qcon_singlet.rtt.json", 0.015625}}) {
    INFO(sg.file);
    const CompiledSystem cs = compiled(load_example(sg.file));
    const auto fo = rtt::paraxial::first_order(cs, PathId{0}, cs.reference_wavelength());
    const Focal f = singlet(sg.cv, 1.5168, 4.0);
    REQUIRE(fo.efl.has_value());
    REQUIRE(fo.bfl.has_value());
    INFO("efl " << got(fo.efl) << " expected " << f.efl);
    CHECK(rel_ok(got(fo.efl), f.efl));
    CHECK(rel_ok(got(fo.bfl), f.bfl));
  }
}

TEST_CASE("m8 case 4: cylinder and biconic focus per section", "[m8]") {
  // A plane back adds no power: efl = 1 / ((n - 1) c) of the section, bfl = efl - t / n (y-nu
  // above). The cylinder (R = 51.68 mm in y) has no power in x: that section is afocal.
  const double n = 1.5168;
  const double t = 4.0;
  {
    const CompiledSystem cs = compiled(load_example("anamorphic/cylinder_planoconvex.rtt.json"));
    const std::uint16_t wl = cs.reference_wavelength();
    CHECK_THROWS_AS(rtt::paraxial::first_order(cs, PathId{0}, wl), ParaxialError);
    const auto fy = rtt::paraxial::first_order_section(cs, PathId{0}, wl, Section::y);
    const auto fx = rtt::paraxial::first_order_section(cs, PathId{0}, wl, Section::x);
    const Focal f = singlet(1.0 / 51.68, n, t);
    REQUIRE(fy.efl.has_value());
    REQUIRE(fy.bfl.has_value());
    CHECK(rel_ok(got(fy.efl), f.efl));
    CHECK(rel_ok(got(fy.bfl), f.bfl));
    CHECK_FALSE(fx.efl.has_value());
    // A ray in the x-z plane meets surfaces whose sag does not depend on x (cylinder) or is
    // plane: every normal has no x and no y component there, so the ray leaves parallel.
    const RayBatch r = traced(cs, rays_from(cs, {{3.0, 0.0, -1.0, 0.0, 0.0, 1.0}}));
    REQUIRE(r.status()[0] == RayStatus::Alive);
    CHECK(r.dir_x()[0] == 0.0);
    CHECK(r.dir_y()[0] == 0.0);
    CHECK(r.pos_x()[0] == 3.0);
  }
  {
    const CompiledSystem cs = compiled(load_example("anamorphic/biconic_singlet.rtt.json"));
    const std::uint16_t wl = cs.reference_wavelength();
    CHECK_THROWS_AS(rtt::paraxial::first_order(cs, PathId{0}, wl), ParaxialError);
    for (const auto& [section, radius] :
         {std::pair{Section::x, 77.52}, std::pair{Section::y, 51.68}}) {
      INFO("section " << (section == Section::x ? "x" : "y"));
      const auto fo = rtt::paraxial::first_order_section(cs, PathId{0}, wl, section);
      const Focal f = singlet(1.0 / radius, n, t);
      REQUIRE(fo.efl.has_value());
      REQUIRE(fo.bfl.has_value());
      CHECK(rel_ok(got(fo.efl), f.efl));
      CHECK(rel_ok(got(fo.bfl), f.bfl));
    }
  }
}

TEST_CASE("m8 case 5: Zernike terms (Noll and Fringe) change the paraxial curvature", "[m8]") {
  // tests/reference/m8/zernike_mirror: sphere R = -200 mm, Noll j = 11 with a = 2e-5 mm and
  // Fringe 4 with a = -1e-5 mm, both normalised to R_n = 22 mm. Noll 11 is (n, m) = (4, 0) with
  // the factor sqrt(n + 1) (Niu/Tian Eq. (3)); Fringe 4 is (2, 0) without a factor (Niu/Tian
  // Table 4, Eq. (22)). A term a N R_n^0(r / R_n) adds a N k_0 at the vertex and
  // 2 a N k_2 / R_n^2 to the curvature, k_i the coefficient of rho^i of Eq. (4).
  const CompiledSystem cs = compiled(load_reference("zernike_mirror.rtt.json"));
  const double rn = 22.0;
  const double noll = 2e-5 * std::sqrt(runtime(5.0));
  const double fringe = -1e-5;
  const double c =
      -1.0 / 200.0 +
      2.0 * (noll * radial_coefficient(4, 0, 2) + fringe * radial_coefficient(2, 0, 2)) / (rn * rn);
  const double zv = noll * radial_coefficient(4, 0, 0) + fringe * radial_coefficient(2, 0, 0);
  REQUIRE(radial_coefficient(4, 0, 2) == -6.0);  // R_4^0 = 6 rho^4 - 6 rho^2 + 1
  REQUIRE(radial_coefficient(2, 0, 0) == -1.0);  // R_2^0 = 2 rho^2 - 1
  const double efl = mirror_efl(c, image_index(cs));
  // Light leaves towards -z; F' lies the rear focal length n' / Phi = 1 / (2 |c|) before the
  // vertex.
  const double focus = zv + 1.0 / (2.0 * c);
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, cs.reference_wavelength());
  REQUIRE(fo.efl.has_value());
  REQUIRE(fo.rear_focal_z.has_value());
  INFO("efl " << got(fo.efl) << " expected " << efl);
  CHECK(rel_ok(got(fo.efl), efl));
  CHECK(std::abs(got(fo.rear_focal_z) - focus) <= kRel * std::abs(focus));
  CHECK(std::abs(focus + 100.0) > 1e-3);  // content guard: the terms change the focus

  // A real ray parallel to the axis at height h crosses the axis at z(h) = z_F + A h^2 + B h^4
  // + ... Richardson: E(h) = (4 z(h/2) - z(h)) / 3 = z_F + O(h^4); bound e = 2 |E(h) - E(h/2)|
  // (the rule of ADR 0033, point 10), plus the rounding of the crossing: a direction of a few
  // ulp (1e-15) over the focal length f, divided by the slope h / f, f^2 1e-15 / h = 1e-10 mm
  // at h = 0.1 mm, times 5/3 of the extrapolation; rounding bound 1e-9 mm.
  const std::array<double, 3> h{0.4, 0.2, 0.1};
  std::vector<Start> starts;
  starts.reserve(h.size());
  for (const double hi : h) starts.push_back({0.0, hi, -60.0, 0.0, 0.0, 1.0});
  const RayBatch r = traced(cs, rays_from(cs, starts));
  std::array<double, 3> z{};
  for (std::size_t i = 0; i < 3; ++i) {
    REQUIRE(r.status()[i] == RayStatus::Alive);
    z[i] = r.pos_z()[i] - r.pos_y()[i] * r.dir_z()[i] / r.dir_y()[i];
  }
  const double e1 = (4.0 * z[1] - z[0]) / 3.0;
  const double e2 = (4.0 * z[2] - z[1]) / 3.0;
  const double e = 2.0 * std::abs(e1 - e2);
  INFO("z(h) " << z[0] << " " << z[1] << " " << z[2] << ", E " << e2 << ", e " << e);
  // Acceptance limit of the extrapolation itself, as the limit 1e-3 on e in ADR 0033, point 10:
  // e must stay far below the effect it resolves (the h^2 term, 1.9e-4 mm at h = 0.4 mm), so
  // that the comparison with z_F tests the paraxial focus and not the extrapolation.
  CHECK(e <= 1e-6);
  CHECK(std::abs(e2 - focus) <= e + 1e-9);
  CHECK(std::abs(z[0] - focus) > 1e-6);  // content guard: the real ray is not paraxial
}

TEST_CASE("m8 case 6: an XY polynomial changes the curvature by 2 c_x2 / R^2", "[m8]") {
  // tests/reference/m8/xy_mirror: z = sum c_k (x/R)^i (y/R)^j ([UG] "Extended Polynomial"),
  // R = 22 mm, c_x2 = c_y2 = 1e-3 mm, fourth order rotationally symmetric, x^6 correction: the
  // vertex curvature is -1/200 + 2 c_x2 / R^2, the vertex stays at 0 (no constant term).
  System s = load_reference("xy_mirror.rtt.json");
  const CompiledSystem cs = compiled(s);
  const double c = -1.0 / 200.0 + 2.0 * 1e-3 / (22.0 * 22.0);
  const double efl = mirror_efl(c, image_index(cs));
  const double focus = 1.0 / (2.0 * c);  // rear focal length 1 / (2 |c|) before the vertex
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, cs.reference_wavelength());
  REQUIRE(fo.efl.has_value());
  REQUIRE(fo.rear_focal_z.has_value());
  INFO("efl " << got(fo.efl) << " expected " << efl);
  CHECK(rel_ok(got(fo.efl), efl));
  CHECK(std::abs(got(fo.rear_focal_z) - focus) <= kRel * std::abs(focus));
  CHECK(std::abs(focus + 100.0) > 1e-2);  // content guard
  CHECK_NOTHROW(rtt::paraxial::seidel(cs, PathId{0}, cs.reference_wavelength()));

  // Variant with a term x^2 y (coefficient 7, index 6; decision of the maintainer, 2026-11-23):
  // at the vertex it has no slope and no Hessian, so the surface keeps slope 0, an isotropic
  // Hessian and order >= 2, and the first order is that of the sphere with c above (ADR 0041,
  // point 2, table: first order). It is not rotationally symmetric to fourth order
  // (rotational4 false, ADR 0051, point 4, all terms of third degree 0), so seidel throws (same
  // table: Seidel only with order 4 and rotational4).
  auto& xy = std::get<rtt::model::XyPolynomial>(element(s, 1).surfaces[0].shape.terms.at(0));
  REQUIRE(xy.coefficients.size() > 6);
  REQUIRE(xy.coefficients[6].value == 0.0);
  xy.coefficients[6] = 1e-5;
  const CompiledSystem bent = compiled(s);
  const auto fb = rtt::paraxial::first_order(bent, PathId{0}, bent.reference_wavelength());
  CHECK(rel_ok(got(fb.efl), efl));
  CHECK(std::abs(got(fb.rear_focal_z) - focus) <= kRel * std::abs(focus));
  CHECK_THROWS_AS(rtt::paraxial::seidel(bent, PathId{0}, bent.reference_wavelength()),
                  ParaxialError);
}

TEST_CASE("m8 case 7: a grid sag of a quadric meets the analytic surface", "[m8]") {
  // tests/reference/m8/grid_sag_mirror: 9 x 9 values of a (x^2 - y^2) / 22^2 at spacing 6 mm,
  // no derivatives given. The central differences (z_{i+1} - z_{i-1}) / (2 h) and the bicubic
  // Hermite interpolation (ADR 0051, point 5; Xu, remainder term with the fourth derivative)
  // are exact for a quadric: the grid is the analytic surface up to rounding in the inner cells.
  // The rays are aimed with the nominal design (the conic without term, ADR 0041, point 3) and
  // the same start rays go through the grid mirror and through the conic with the XY polynomial
  // a x^2 / 22^2 - a y^2 / 22^2. Both hits carry kNewtonTolerance; the normals differ by the
  // rounding of the Hermite sum: the term has values up to 2.4e-5 mm and slopes
  // 2 a |x| / 22^2 <= 2 * 2e-5 * 14 / 484 = 1.2e-6 <= 2e-6 in the inner cells, so its rounding
  // is about 2e-6 eps = 4e-22 in the slope, far below kDir; bound kPos on the detector and for
  // the OPL, kDir for the directions.
  const System grid_model = load_reference("grid_sag_mirror.rtt.json");
  {
    const auto& g = std::get<rtt::model::GridSag>(
        std::get<Element>(grid_model.root.children[1].value).surfaces[0].shape.terms.at(0));
    REQUIRE(g.nx == 9);
    REQUIRE(g.ny == 9);
    REQUIRE(g.dz_dx.empty());
    for (int j = 0; j < 9; ++j) {
      for (int i = 0; i < 9; ++i) {
        const double x = 6.0 * (i - 4);
        const double y = 6.0 * (j - 4);
        const double v = kGridA * (x * x - y * y) / 484.0;
        CHECK(std::abs(g.sag[static_cast<std::size_t>(j * 9 + i)] - v) <= 1e-15 * 2.4e-5);
      }
    }
  }
  const CompiledSystem grid = compiled(grid_model);
  const CompiledSystem nominal = compiled(grid_variant(false));
  const CompiledSystem analytic = compiled(grid_variant(true));
  const std::uint16_t wl = grid.reference_wavelength();
  CHECK_THROWS_AS(rtt::paraxial::first_order(grid, PathId{0}, wl), ParaxialError);
  const auto reference = rtt::trace::aim_reference(nominal, PathId{0}, wl);
  const std::vector<std::uint16_t> fields{0};
  const RayBatch start = rtt::trace::make_rays(grid, PathId{0}, fields, wl, grid_points(),
                                               rtt::trace::Aiming::Real, reference);
  const RayBatch a = traced(grid, start);
  const RayBatch b = traced(analytic, start);
  const RayBatch c = traced(nominal, start);
  REQUIRE(a.size() == 25);
  double largest = 0.0;
  double astig = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    INFO("ray " << i);
    REQUIRE(a.status()[i] == RayStatus::Alive);
    REQUIRE(b.status()[i] == RayStatus::Alive);
    CHECK(a.last_surface()[i] == surface(grid, "F"));
    CHECK(std::abs(start.pos_x()[i]) <= 14.0 + 1e-6);
    CHECK(std::abs(start.pos_y()[i]) <= 14.0 + 1e-6);
    for (const auto& [u, v] :
         {std::pair{a.pos_x()[i], b.pos_x()[i]}, std::pair{a.pos_y()[i], b.pos_y()[i]},
          std::pair{a.pos_z()[i], b.pos_z()[i]}, std::pair{a.opl()[i], b.opl()[i]}}) {
      CHECK(std::abs(u - v) <= kPos);
      largest = std::max(largest, std::abs(u - v));
    }
    CHECK(std::abs(a.dir_x()[i] - b.dir_x()[i]) <= kDir);
    CHECK(std::abs(a.dir_y()[i] - b.dir_y()[i]) <= kDir);
    CHECK(std::abs(a.dir_z()[i] - b.dir_z()[i]) <= kDir);
    astig = std::max(astig, std::hypot(a.pos_x()[i] - c.pos_x()[i], a.pos_y()[i] - c.pos_y()[i]));
  }
  INFO("largest deviation " << largest << " mm, effect of the grid " << astig << " mm");
  CHECK(astig > 1e-5);  // content guard: the grid moves the rays against the bare conic
}

TEST_CASE("m8 case 8: a hexagonal mirror passes 1 - 6 segments of the pupil", "[m8]") {
  // tests/reference/m8/hexagon_mirror: paraboloid behind a stop, collimated beam of radius
  // R = 30 mm (EPD 60), hexagon with circumradius 34 mm on the mirror. Parallel rays meet the
  // mirror at their pupil x, y. The hexagon cuts six segments at the inradius
  // a = 34 cos 30 deg (MathWorld "Regular Polygon", Eq. (4)); the half chord sqrt(R^2 - a^2)
  // (MathWorld "Circular Segment", Eq. (8)) is below the half edge 34 sin 30 deg, so the
  // segments do not overlap: p = 1 - 6 A_seg / (pi R^2) with Eq. (17).
  const CompiledSystem cs = compiled(load_reference("hexagon_mirror.rtt.json"));
  constexpr std::size_t kN = 1000000;
  const double r = 30.0;
  const double a = 34.0 * std::cos(deg(30.0));
  REQUIRE(std::sqrt(r * r - a * a) < 34.0 * std::sin(deg(30.0)));
  const double p = 1.0 - 6.0 * segment_area(r, a) / (kPi * r * r);
  const double bound = 5.0 * std::sqrt(p * (1.0 - p) / static_cast<double>(kN));
  const RayBatch rays = aimed_random(cs, kN);
  const double f = alive_fraction(rays);
  INFO("fraction " << f << ", expected " << p << ", bound " << bound);
  CHECK(std::abs(f - p) <= bound);
  CHECK(1.0 - p > 10.0 * bound);  // content guard: the corners matter
  for (std::size_t i = 0; i < rays.size(); ++i) {
    if (rays.status()[i] != RayStatus::Alive) CHECK(rays.last_surface()[i] == surface(cs, "M"));
  }
}

TEST_CASE("m8 case 9: spider, central and decentred obscurations remove their areas", "[m8]") {
  // tests/reference/m8/spider_telescope: pupil R = 30 mm; at the stop a central disk rho = 10 mm
  // and two vanes |y| <= w, |x| <= w with w = 0.5 mm; on the mirror an ellipse 3 x 1.5 mm at
  // (20, 5) mm, turned by 30 deg. Each vane leaves the strip of the annulus,
  // S(R) - S(rho) with the strip area S of Eq. (17) (MathWorld "Circular Segment"); their
  // crossing lies in the central disk (|x|, |y| <= w < rho), so nothing counts twice. The ellipse
  // has the area pi a b (MathWorld "Ellipse", Eq. (81)) and lies in the annulus, away from the
  // vanes: its bounding box about (20, 5) has the half widths sqrt(a^2 cos^2 t + b^2 sin^2 t)
  // and sqrt(a^2 sin^2 t + b^2 cos^2 t), checked below. Parallel rays meet the mirror at their
  // pupil x, y. p = [pi (R^2 - rho^2) - 2 (S(R) - S(rho)) - pi a b] / (pi R^2).
  const CompiledSystem cs = compiled(load_reference("spider_telescope.rtt.json"));
  constexpr std::size_t kN = 1000000;
  const double r = 30.0;
  const double rho = 10.0;
  const double w = 0.5;
  const double ea = 3.0;
  const double eb = 1.5;
  const double t = deg(30.0);
  const double hx =
      std::sqrt(ea * ea * std::cos(t) * std::cos(t) + eb * eb * std::sin(t) * std::sin(t));
  const double hy =
      std::sqrt(ea * ea * std::sin(t) * std::sin(t) + eb * eb * std::cos(t) * std::cos(t));
  REQUIRE(5.0 - hy > w);
  REQUIRE(20.0 - hx > w);
  REQUIRE(std::hypot(20.0 - hx, 5.0 - hy) > rho);
  REQUIRE(std::hypot(20.0 + hx, 5.0 + hy) < r);
  const double open =
      kPi * (r * r - rho * rho) - 2.0 * (strip_area(r, w) - strip_area(rho, w)) - kPi * ea * eb;
  const double p = open / (kPi * r * r);
  const double bound = 5.0 * std::sqrt(p * (1.0 - p) / static_cast<double>(kN));
  const RayBatch rays = aimed_random(cs, kN);
  const double f = alive_fraction(rays);
  INFO("fraction " << f << ", expected " << p << ", bound " << bound);
  CHECK(std::abs(f - p) <= bound);
  // Content guard: the ellipse alone is resolved (its share 5.0e-3 against the bound).
  CHECK(ea * eb / (r * r) > 2.0 * bound);
  // Where the rays end: at the stop (disk and vanes, share [pi rho^2 + 2 (S(R) - S(rho))] /
  // (pi R^2)) and at the mirror (ellipse, a b / R^2), each binomial with its own 5 sigma bound.
  const std::uint32_t stop = surface(cs, "STO");
  const std::uint32_t mirror = surface(cs, "M");
  std::size_t at_stop = 0;
  std::size_t at_mirror = 0;
  for (std::size_t i = 0; i < rays.size(); ++i) {
    if (rays.status()[i] == RayStatus::Alive) continue;
    if (rays.last_surface()[i] == stop) ++at_stop;
    if (rays.last_surface()[i] == mirror) ++at_mirror;
  }
  const double n = static_cast<double>(kN);
  const double p_stop =
      (kPi * rho * rho + 2.0 * (strip_area(r, w) - strip_area(rho, w))) / (kPi * r * r);
  const double p_mirror = ea * eb / (r * r);
  const double f_stop = static_cast<double>(at_stop) / n;
  const double f_mirror = static_cast<double>(at_mirror) / n;
  INFO("stop " << f_stop << " expected " << p_stop << ", mirror " << f_mirror << " expected "
               << p_mirror);
  CHECK(std::abs(f_stop - p_stop) <= 5.0 * std::sqrt(p_stop * (1.0 - p_stop) / n));
  CHECK(std::abs(f_mirror - p_mirror) <= 5.0 * std::sqrt(p_mirror * (1.0 - p_mirror) / n));
  CHECK(at_stop + at_mirror + static_cast<std::size_t>(std::llround(f * n)) == kN);
}

TEST_CASE("m8 case 10: aiming and tracing are bitwise the same at 1, 4 and all threads", "[m8]") {
  // Rule 7 (ADR 0004): the rays of cases 1 to 9 and the spot of the hexagonal mirror (an
  // analysis on top of aiming and trace) at the default thread count and with 1 and 4.
  const CompiledSystem axicon = compiled(load_example("rotational/axicon_planoconvex.rtt.json"));
  const CompiledSystem qbfs = compiled(qbfs_open_stop());
  const CompiledSystem cylinder =
      compiled(load_example("anamorphic/cylinder_planoconvex.rtt.json"));
  const CompiledSystem biconic = compiled(load_example("anamorphic/biconic_singlet.rtt.json"));
  std::vector<CompiledSystem> aimed;
  for (const char* f : {"zernike_mirror.rtt.json", "xy_mirror.rtt.json", "hexagon_mirror.rtt.json",
                        "spider_telescope.rtt.json"}) {
    aimed.push_back(compiled(load_reference(f)));
  }
  aimed.push_back(compiled(load_example("rotational/odd_asphere_singlet.rtt.json")));
  aimed.push_back(compiled(load_example("forbes/qcon_singlet.rtt.json")));
  const CompiledSystem grid = compiled(load_reference("grid_sag_mirror.rtt.json"));
  const CompiledSystem nominal = compiled(grid_variant(false));
  const auto run = [&] {
    std::vector<double> v = columns(axicon_rays(axicon));
    const auto add = [&v](const RayBatch& r) {
      const std::vector<double> c = columns(r);
      v.insert(v.end(), c.begin(), c.end());
    };
    add(traced(qbfs, rays_from(qbfs, {{0.0, 12.7, -1.0, 0.0, 0.0, 1.0}})));
    const std::vector<std::uint16_t> fields{0};
    for (const CompiledSystem* cs : {&cylinder, &biconic}) {
      add(traced(*cs, rtt::trace::make_rays(*cs, PathId{0}, fields, cs->reference_wavelength(),
                                            rtt::trace::HexapolarPupil{4})));
    }
    for (const CompiledSystem& cs : aimed) add(aimed_random(cs, 2000));
    const std::uint16_t wl = nominal.reference_wavelength();
    const auto reference = rtt::trace::aim_reference(nominal, PathId{0}, wl);
    add(traced(grid, rtt::trace::make_rays(grid, PathId{0}, fields, wl, grid_points(),
                                           rtt::trace::Aiming::Real, reference)));
    const rtt::analysis::SpotDiagram spot =
        rtt::analysis::spot(aimed[2], PathId{0}, 0, std::nullopt, {rtt::trace::GridPupil{21}});
    for (const rtt::analysis::SpotPoint& q : spot.points) {
      v.push_back(q.x);
      v.push_back(q.y);
    }
    v.push_back(spot.stats.rms_centroid);
    v.push_back(spot.vignetted_fraction);
    return v;
  };
  const std::vector<double> all = run();
  REQUIRE(all.size() > 10000);
  for (const int threads : {1, 4}) {
    INFO(threads << " threads");
    const tbb::global_control limit(tbb::global_control::max_allowed_parallelism,
                                    static_cast<std::size_t>(threads));
    CHECK(same_bits(run(), all));
  }
}
