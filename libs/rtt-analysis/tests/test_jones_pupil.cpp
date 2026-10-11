// Jones pupil and its maps (#208, ADR 0034, points 4 to 6 and 9, with the addendum #208).
//
// Sources: W.-S. T. Lam, "Anisotropic Ray Trace" (docs/quellen.md): Appendix B, p. 278 (printed
// dipole formulas) and p. 280 (double-pole formulas), Eq. (1.2) (Jones pupil), Eq. (3.4) (P of a
// reflection), Eqs. (4.5)/(4.6) (Q); S. J. Byrnes, arXiv:1603.02720v5, Eq. (6) (Fresnel
// amplitudes, rtt::polar::fresnel_at_angle); the project's Stokes convention
// (docs/architecture.md, "Händigkeit und Stokes").
//
// The reference values of the coaxial metal mirrors (decision of the maintainer for #208: the
// convenience form is accepted on rotationally symmetric paths, the 45 deg fold follows with the
// main form) are built per pupil point without the P of the tracer and without the bases of the
// implementation: the directions of the ray from a recorded path (RayPaths), the lens as the
// rotations R(k_in -> k_out) of its ideal anti-reflection surfaces (ADR 0034, point 1), the
// mirror from Fresnel at the local angle of incidence with s = k x N (Lam (3.4)), the bases from
// Lam's formulas, the phase from W of opd_points().
//
// Tolerances: J, D and delta are formed from unit vectors and Fresnel amplitudes of modulus
// about 0.9 with a few dozen roundings each, a few 1e-16; the reference values take another
// path through the same quantities (rotations against a product of rotations), so 1e-10 as in
// ADR 0034, point 9 leaves a wide margin and still separates the effects of a few 1e-3 that the
// metal mirrors give at these angles (content guards).

#include <tbb/global_control.h>

#include <Eigen/Geometry>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/jones.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/polar/fresnel.hpp"
#include "rtt/polar/prt.hpp"
#include "rtt/polar/prt_analysis.hpp"
#include "rtt/trace/apply_event.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using Cx = std::complex<double>;
using rtt::analysis::AnalysisError;
using rtt::analysis::JonesBasis;
using rtt::analysis::JonesPupil;
using rtt::analysis::JonesPupilOptions;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::math::CMat3;
using rtt::math::Mat3;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::EventKind;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::trace::RayStatus;

namespace {

constexpr double kPi = std::numbers::pi;
const Cx kMetal(1.2, 7.26);  // invented test metal of M3 (m3/metal_periscope)

System load(const std::string& name) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + name);
}

CompiledSystem compiled(const System& s) {
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

Element& child(System& s, std::size_t i) {
  return std::get<Element>(s.root.children[i].value);
}

/// Lam, Appendix B, p. 280: the double-pole basis for k = (x, y, z) about the chief ray along z
/// with the local x axis (1, 0, 0).
std::array<Vec3, 2> lam_double_pole(const Vec3& k) {
  const double x = k.x();
  const double y = k.y();
  const double z = k.z();
  return {Vec3(z + y * y / (1.0 + z), -x * y / (1.0 + z), -x),
          Vec3(-x * y / (1.0 + z), z + x * x / (1.0 + z), -y)};
}

/// Lam, Appendix B, p. 278, the printed dipole formulas: poles at +-z.
std::array<Vec3, 2> lam_dipole_printed(const Vec3& k) {
  const double rho = std::hypot(k.x(), k.y());
  return {Vec3(-k.y() / rho, k.x() / rho, 0.0),
          Vec3(-k.x() * k.z() / rho, -k.y() * k.z() / rho, rho)};
}

/// The double-pole basis about a chief frame (x_c, y_c, k_c), from Lam's formulas in the frame
/// of the chief ray: T Lam(T^T k) with T = [x_c y_c k_c].
std::array<Vec3, 2> lam_double_pole_about(const Vec3& x_c,
                                          const Vec3& y_c,
                                          const Vec3& k_c,
                                          const Vec3& k) {
  Mat3 t;
  t.col(0) = x_c;
  t.col(1) = y_c;
  t.col(2) = k_c;
  const std::array<Vec3, 2> local = lam_double_pole(t.transpose() * k);
  return {t * local[0], t * local[1]};
}

double max_diff(const std::array<Vec3, 2>& a, const std::array<Vec3, 2>& b) {
  return std::max((a[0] - b[0]).cwiseAbs().maxCoeff(), (a[1] - b[1]).cwiseAbs().maxCoeff());
}

/// P of a reflection with amplitudes r_s, r_p (Lam, Eq. (3.4)): s = k_in x N normalised
/// (any s perpendicular to k_in at normal incidence), p = k x s.
CMat3 reflection(const Vec3& k_in, const Vec3& k_out, const Vec3& n, Cx rs, Cx rp) {
  Vec3 s = k_in.cross(n);
  s = s.norm() > 1e-12 ? Vec3(s.normalized()) : Vec3(k_in.unitOrthogonal());
  const Vec3 p_in = k_in.cross(s);
  const Vec3 p_out = k_out.cross(s);
  const Mat3 ss = s * s.transpose();
  const Mat3 pp = p_out * p_in.transpose();
  const Mat3 kk = k_out * k_in.transpose();
  return rs * ss.cast<Cx>() + rp * pp.cast<Cx>() + kk.cast<Cx>();
}

/// Expected P, Q and directions of one ray (reference values, file comment).
struct Expected {
  CMat3 p = CMat3::Identity();
  Mat3 q = Mat3::Identity();
  Vec3 k0 = Vec3::UnitZ();
  Vec3 k = Vec3::UnitZ();
};

/// J = exp(-2 pi i W lambda_ref / lambda) B_out^T P B_in (ADR 0034, point 5).
std::array<Cx, 4> expected_j(const Expected& e,
                             const std::array<Vec3, 2>& in,
                             const std::array<Vec3, 2>& out,
                             double w,
                             double lambda_ref,
                             double lambda) {
  const Cx factor = std::polar(1.0, -2.0 * kPi * w * lambda_ref / lambda);
  std::array<Cx, 4> j{};
  for (int r = 0; r < 2; ++r) {
    for (int c = 0; c < 2; ++c) {
      j[static_cast<std::size_t>(2 * r + c)] =
          factor * out[static_cast<std::size_t>(r)].cast<Cx>().dot(
                       e.p * in[static_cast<std::size_t>(c)].cast<Cx>());
    }
  }
  return j;
}

/// singlet_const (f = 100 mm, EPD 20: F/5) in vacuum with ideal anti-reflection surfaces (P = Q
/// = R(k_in -> k_out), not polarizing) and a flat mirror of the test metal at z = 60 that
/// reflects the converging beam back to a detector at z = 13.637 (the reflected focus), path
/// STO, L1.S1, L1.S2, M, IMG; rotationally symmetric, so make_rays and the OPD work.
System convergent_mirror() {
  System s = load("m1/singlet_const.rtt.json");
  s.environment.medium = "VACUUM";
  for (Surface& surface : child(s, 1).surfaces) {
    surface.interaction = rtt::model::IdealAntiReflection{};
  }
  child(s, 2).pose = Pose::along_z(13.637);
  Surface m;
  m.id = SurfaceId("M");
  s.root.children.push_back(
      {Element{"M", ElementKind::Mirror, Pose::along_z(60.0), std::string("CONST:1.2,7.26"), {m}}});
  s.paths = {{"fold back",
              false,
              {{SurfaceId("STO"), EventKind::Transmit, 0},
               {SurfaceId("L1.S1"), EventKind::Refract, 0},
               {SurfaceId("L1.S2"), EventKind::Refract, 0},
               {SurfaceId("M"), EventKind::Reflect, 0},
               {SurfaceId("IMG"), EventKind::Transmit, 0}}}};
  return s;
}

/// m2/paraboloid_stop (R = -200, k = -1, EPD 60) in vacuum with the test metal as mirror material:
/// angles of incidence up to about 8.5 deg.
System metal_paraboloid() {
  System s = load("m2/paraboloid_stop.rtt.json");
  s.environment.medium = "VACUUM";
  child(s, 1).material = "CONST:1.2,7.26";
  return s;
}

/// The rays of `sampling` traced with the path recorded (the same rays as jones_pupil()).
rtt::trace::RayPaths recorded(const CompiledSystem& cs,
                              std::uint16_t field,
                              std::uint16_t wavelength,
                              const rtt::trace::PupilSampling& sampling) {
  const std::uint16_t fields[] = {field};
  rtt::trace::RayBatch rays =
      rtt::trace::make_rays(cs, PathId{0}, fields, wavelength, sampling, rtt::trace::Aiming::Real);
  rtt::trace::RayPaths paths;
  [[maybe_unused]] const auto stats =
      rtt::trace::SequentialTracer().trace(cs, PathId{0}, rays, paths);
  return paths;
}

/// The direction of the chief ray of the reference wavelength on the image surface.
Vec3 chief_direction(const CompiledSystem& cs, std::uint16_t field) {
  const std::uint16_t fields[] = {field};
  rtt::trace::RayBatch rays =
      rtt::trace::make_rays(cs, PathId{0}, fields, cs.reference_wavelength(),
                            rtt::trace::SinglePupilPoint{0.0, 0.0}, rtt::trace::Aiming::Real);
  [[maybe_unused]] const auto stats = rtt::trace::SequentialTracer().trace(cs, PathId{0}, rays);
  return {rays.dir_x()[0], rays.dir_y()[0], rays.dir_z()[0]};
}

/// Compares jones_pupil() with the reference values per point (file comment); `expected` gives
/// P, Q and the directions of recorded ray r. Returns the number of compared points.
template <class ExpectedOf>
std::size_t compare(const JonesPupil& jp,
                    const rtt::trace::RayPaths& paths,
                    const rtt::analysis::OpdPupilPoints& opd,
                    const Vec3& k0_c,
                    const Vec3& k_c,
                    const Vec3& x_b,
                    ExpectedOf expected) {
  // Chief frames from Lam's formulas (input) and the PSF plane (output, ADR 0032, point 4).
  const std::array<Vec3, 2> in_c = lam_double_pole(k0_c);
  const Vec3 x_out = (x_b - x_b.dot(k_c) * k_c).normalized();
  const Vec3 y_out = k_c.cross(x_out);
  CHECK((jp.x_in - in_c[0]).norm() <= 1e-14);
  CHECK((jp.y_in - in_c[1]).norm() <= 1e-14);
  CHECK((jp.x_out - x_out).norm() <= 1e-14);
  CHECK((jp.y_out - y_out).norm() <= 1e-14);
  std::size_t compared = 0;
  REQUIRE(jp.status.size() == paths.ray_count());
  for (std::size_t r = 0; r < paths.ray_count(); ++r) {
    INFO("point " << r << " at (" << jp.px[r] << ", " << jp.py[r] << ")");
    CHECK(jp.status[r] == opd.points[r].status);
    CHECK(jp.w[r] == opd.points[r].w);
    if (jp.status[r] != RayStatus::Alive) continue;
    const Expected e = expected(r);
    const auto in = lam_double_pole_about(in_c[0], in_c[1], k0_c, e.k0);
    const auto out = lam_double_pole_about(x_out, y_out, k_c, e.k);
    const auto j =
        expected_j(e, in, out, opd.points[r].w, jp.reference_wavelength_um, jp.wavelength_um);
    for (std::size_t i = 0; i < 4; ++i) CHECK(std::abs(jp.j[4 * r + i] - j[i]) <= 1e-10);
    CHECK(jp.transmission[r] == paths.weight[r * paths.slots + paths.count[r] - 1]);
    const auto d = rtt::polar::diattenuation(e.p, e.k0, e.k);
    CHECK(std::abs(jp.diattenuation[r] - d.value) <= 1e-10);
    const auto ret = rtt::polar::physical_retardance(e.p, e.q, e.k0);
    CHECK(std::abs(jp.retardance[r] - ret.value) <= 1e-10);
    ++compared;
  }
  return compared;
}

bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

bool same_bits(const std::vector<Cx>& a, const std::vector<Cx>& b) {
  return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(Cx)) == 0;
}

}  // namespace

// ------------------------------------------------------------------- bases (point 5) -----

TEST_CASE("double-pole basis: Lam, Appendix B, p. 280, also about a turned chief frame",
          "[jones][basis]") {
  // R(k_c -> k) (x_c, y_c) is Lam's double-pole basis (ADR 0034, point 5, checked there with a
  // script); here against the printed formulas at several k, and about a chief frame turned by an
  // arbitrary rotation T: T Lam(T^T k). 1e-14 (rounding of the rotation against the formula).
  for (const Vec3& raw : {Vec3(0.0, 0.0, 1.0), Vec3(0.3, -0.2, 0.93), Vec3(1e-3, 0.0, 1.0),
                          Vec3(0.6, 0.8, 0.0), Vec3(0.1, 0.1, -0.9)}) {
    const Vec3 k = raw.normalized();
    INFO("k = " << k.transpose());
    const auto basis =
        rtt::analysis::double_pole_basis(Vec3::UnitX(), Vec3::UnitY(), Vec3::UnitZ(), k);
    REQUIRE(basis);
    CHECK(max_diff(*basis, lam_double_pole(k)) <= 1e-14);
    const Mat3 t = Eigen::AngleAxisd(0.7, Vec3(1.0, -2.0, 0.5).normalized()).toRotationMatrix();
    const auto turned = rtt::analysis::double_pole_basis(t * Vec3::UnitX(), t * Vec3::UnitY(),
                                                         t * Vec3::UnitZ(), t * k);
    REQUIRE(turned);
    const auto lam = lam_double_pole(k);
    CHECK(max_diff(*turned, {t * lam[0], t * lam[1]}) <= 1e-14);
  }
  // The pole -k_c: no basis within 1e-6 rad, a basis at 1e-5 rad.
  const auto at = [](double angle) {
    return rtt::analysis::double_pole_basis(Vec3::UnitX(), Vec3::UnitY(), Vec3::UnitZ(),
                                            Vec3(std::sin(angle), 0.0, -std::cos(angle)));
  };
  CHECK_FALSE(at(0.0).has_value());
  CHECK_FALSE(at(5e-7).has_value());
  CHECK(at(1e-5).has_value());
}

TEST_CASE("dipole basis: Lam's construction with the pole axis as a parameter (p. 278)",
          "[jones][basis]") {
  // Lam's printed dipole formulas (p. 278) have their poles at +-z, his text puts them at +-y (the
  // pitfall of ADR 0034). x = pole x k / |pole x k|, y = k x x is the construction of both; with
  // the pole axis z it gives the printed formulas. 1e-14.
  for (const Vec3& raw : {Vec3(0.3, -0.2, 0.93), Vec3(0.6, 0.8, 0.0), Vec3(-0.5, 0.1, -0.7)}) {
    const Vec3 k = raw.normalized();
    INFO("k = " << k.transpose());
    const auto basis = rtt::analysis::dipole_basis(Vec3::UnitZ(), k);
    REQUIRE(basis);
    CHECK(max_diff(*basis, lam_dipole_printed(k)) <= 1e-14);
  }
  // Poles at +-y_c: none within 1e-6 rad.
  CHECK_FALSE(rtt::analysis::dipole_basis(Vec3::UnitY(), Vec3::UnitY()).has_value());
  CHECK_FALSE(rtt::analysis::dipole_basis(Vec3::UnitY(), -Vec3::UnitY()).has_value());
  CHECK_FALSE(
      rtt::analysis::dipole_basis(Vec3::UnitY(), Vec3(5e-7, 1.0, 0.0).normalized()).has_value());
  CHECK(rtt::analysis::dipole_basis(Vec3::UnitY(), Vec3(1e-5, 1.0, 0.0).normalized()).has_value());
  // At the chief ray (k = k_c, y_c = k_c x x_c) the dipole basis is the chief frame.
  const auto chief = rtt::analysis::dipole_basis(Vec3::UnitY(), Vec3::UnitZ());
  REQUIRE(chief);
  CHECK(max_diff(*chief, {Vec3::UnitX(), Vec3::UnitY()}) <= 1e-15);
}

TEST_CASE("input chief frame: R(z -> k0_c) (x, y), and diag(1, -1, -1) for -z", "[jones][basis]") {
  const auto z = rtt::analysis::input_chief_frame(Vec3::UnitZ());
  CHECK(z[0] == Vec3::UnitX());
  CHECK(z[1] == Vec3::UnitY());
  const auto minus_z = rtt::analysis::input_chief_frame(Vec3(0.0, 0.0, -1.0));
  CHECK(minus_z[0] == Vec3::UnitX());
  CHECK(minus_z[1] == Vec3(0.0, -1.0, 0.0));
  const Vec3 k = Vec3(0.2, -0.4, 0.8).normalized();
  const auto f = rtt::analysis::input_chief_frame(k);
  CHECK(max_diff(f, lam_double_pole(k)) <= 1e-14);  // R(z -> k) is Lam's double pole (p. 280)
  CHECK((f[0].cross(f[1]) - k).norm() <= 1e-15);
}

TEST_CASE("polarization ellipse: orientation and ellipticity with the project's handedness",
          "[jones][basis]") {
  // e = a e_1 + b e_2; S1 = |a|^2 - |b|^2, S2 = 2 Re(a* b), S3 = 2 Im(a b*) (docs/architecture.md,
  // "Händigkeit und Stokes": (i e_1 + e_2)/sqrt 2 is right-circular, S3 > 0).
  // psi = arg(S1 + i S2)/2 in [0, pi), chi = asin(S3/S0)/2. 1e-15.
  const double h = 1.0 / std::sqrt(2.0);
  struct Case {
    Cx a, b;
    double psi, chi;
  };
  const double t = 20.0 * kPi / 180.0;
  for (const Case& c :
       {Case{1.0, 0.0, 0.0, 0.0}, Case{0.0, 1.0, kPi / 2.0, 0.0}, Case{h, h, kPi / 4.0, 0.0},
        Case{h, -h, 3.0 * kPi / 4.0, 0.0},
        Case{Cx(0.0, h), h, 0.0, kPi / 4.0},  // psi of a circle: atan2(0, 0) = 0
        Case{h, Cx(0.0, h), 0.0, -kPi / 4.0}, Case{std::cos(t), Cx(0.0, std::sin(t)), 0.0, -t}}) {
    INFO("a = " << c.a << ", b = " << c.b);
    const auto e = rtt::analysis::polarization_ellipse(c.a, c.b);
    CHECK(std::abs(e.orientation - c.psi) <= 1e-15);
    CHECK(std::abs(e.ellipticity - c.chi) <= 1e-15);
  }
  const auto none = rtt::analysis::polarization_ellipse(0.0, 0.0);
  CHECK(std::isnan(none.orientation));
  CHECK(std::isnan(none.ellipticity));
}

// ------------------------------------------- coaxial metal mirrors (points 5, 6 and 9) -----

TEST_CASE("Jones pupil of a converging F/5 beam on a flat metal mirror", "[jones][acceptance]") {
  // convergent_mirror() (file comment): per point J, D and the physical retardance against the
  // reference values from the recorded directions, Fresnel at the local angle and Lam's bases;
  // fields 0 and 2 (5 deg), the reference wavelength and 0.4861 um (phase scale lambda_ref /
  // lambda). 1e-10 (file comment).
  const CompiledSystem cs = compiled(convergent_mirror());
  const Vec3 n_mirror = Vec3::UnitZ();
  const rtt::trace::GridPupil grid{9};
  for (const std::uint16_t field : {std::uint16_t{0}, std::uint16_t{2}}) {
    for (const std::uint16_t wl : {std::uint16_t{1}, std::uint16_t{0}}) {
      INFO("field " << field << ", wavelength " << wl);
      JonesPupilOptions o;
      o.sampling = grid;
      const JonesPupil jp = rtt::analysis::jones_pupil(cs, PathId{0}, field, wl, o);
      rtt::analysis::OpdOptions oo;
      const auto opd = rtt::analysis::opd_points(cs, PathId{0}, field, wl, grid, oo);
      const auto paths = recorded(cs, field, wl, grid);
      // k0 of the chief ray of the reference wavelength: the field direction (object at
      // infinity, the same for all wavelengths).
      const std::uint16_t fields[] = {field};
      const auto chief_start =
          rtt::trace::make_rays(cs, PathId{0}, fields, cs.reference_wavelength(),
                                rtt::trace::SinglePupilPoint{0.0, 0.0}, rtt::trace::Aiming::Real);
      const Vec3 k0c(chief_start.dir_x()[0], chief_start.dir_y()[0], chief_start.dir_z()[0]);
      const auto rs_rp = [&](const Vec3& k_in) {
        const double cos_i = std::abs(k_in.dot(n_mirror));
        const auto a = rtt::polar::fresnel_at_angle<double>(Cx(1.0), kMetal, cos_i);
        return std::pair{a.rs, a.rp};
      };
      const auto expected = [&](std::size_t r) {
        // Slots: 0 start, 1 after STO, 2 after L1.S1, 3 after L1.S2, 4 after M, 5 after IMG.
        Expected e;
        e.k0 = paths.direction_at(r, 0);
        const Vec3 k1 = paths.direction_at(r, 2);
        const Vec3 k2 = paths.direction_at(r, 3);
        e.k = paths.direction_at(r, 4);
        const Mat3 lens =
            rtt::trace::rotation_between(k1, k2) * rtt::trace::rotation_between(e.k0, k1);
        const auto [rs, rp] = rs_rp(k2);
        e.p = reflection(k2, e.k, n_mirror, rs, rp) * lens.cast<Cx>();
        e.q = rtt::polar::geometric_transform<double>(k2, e.k, n_mirror, true) * lens;
        return e;
      };
      const std::size_t compared =
          compare(jp, paths, opd, k0c, chief_direction(cs, field), Vec3::UnitX(), expected);
      CHECK(compared >= 50);  // content guard: the 9 x 9 grid inside the circle
      // Content guard: the metal acts, D and delta differ from 0 at the rim and J varies.
      double d_max = 0.0;
      double delta_max = 0.0;
      for (std::size_t r = 0; r < jp.status.size(); ++r) {
        if (jp.status[r] != RayStatus::Alive) continue;
        d_max = std::max(d_max, jp.diattenuation[r]);
        delta_max = std::max(delta_max, jp.retardance[r]);
      }
      CHECK(d_max > 1e-4);
      CHECK(delta_max > 1e-4);
    }
  }
}

TEST_CASE("Jones pupil of a metal paraboloid", "[jones][acceptance]") {
  // metal_paraboloid(): normal of the conic k = -1 at the hit point (x, y): N ~ (-c x, -c y, 1)
  // with c = 1/R (sag z = c r^2 / 2), angles of incidence up to 8.5 deg; per point as above.
  const CompiledSystem cs = compiled(metal_paraboloid());
  const double c = -1.0 / 200.0;
  const rtt::trace::HexapolarPupil hex{4};
  JonesPupilOptions o;
  o.sampling = hex;
  const JonesPupil jp = rtt::analysis::jones_pupil(cs, PathId{0}, 0, 0, o);
  const auto opd = rtt::analysis::opd_points(cs, PathId{0}, 0, 0, hex, rtt::analysis::OpdOptions{});
  const auto paths = recorded(cs, 0, 0, hex);
  const auto expected = [&](std::size_t r) {
    // Slots: 0 start, 1 after STO, 2 after M, 3 after F.
    Expected e;
    e.k0 = paths.direction_at(r, 0);
    const Vec3 k_in = paths.direction_at(r, 1);
    e.k = paths.direction_at(r, 2);
    const Vec3 hit = paths.position_at(r, 2);  // the mirror vertex is the global origin
    const Vec3 n = Vec3(-c * hit.x(), -c * hit.y(), 1.0).normalized();
    const auto a = rtt::polar::fresnel_at_angle<double>(Cx(1.0), kMetal, std::abs(k_in.dot(n)));
    e.p = reflection(k_in, e.k, n, a.rs, a.rp);
    e.q = rtt::polar::geometric_transform<double>(k_in, e.k, n, true);
    return e;
  };
  const std::size_t compared =
      compare(jp, paths, opd, Vec3::UnitZ(), chief_direction(cs, 0), Vec3::UnitX(), expected);
  CHECK(compared == 61);  // hexapolar ring 4: 1 + 6 (1 + 2 + 3 + 4) points, all arrive
  double d_rim = 0.0;
  for (std::size_t r = 0; r < jp.status.size(); ++r) d_rim = std::max(d_rim, jp.diattenuation[r]);
  CHECK(d_rim > 1e-3);  // content guard
}

// ------------------------------------------------------------- frames and options -----

TEST_CASE("Jones pupil: a detector turned about z turns only the output basis", "[jones]") {
  // The detector turned by 30 deg in its plane (rotationally symmetric still): the output frame
  // turns, x'_new = a x' + b y', y'_new = a y' - b x' (k x x' = y', k x y' = -x'), so
  // J_new = [[a, b], [-b, a]] J; D and delta do not depend on the output basis (bitwise).
  // 1e-14 (one 2 x 2 product).
  System s = metal_paraboloid();
  JonesPupilOptions o;
  o.sampling = rtt::trace::HexapolarPupil{2};
  const JonesPupil straight = rtt::analysis::jones_pupil(compiled(s), PathId{0}, 0, 0, o);
  child(s, 2).pose.rotation_deg[2] = Param(30.0);
  const JonesPupil turned = rtt::analysis::jones_pupil(compiled(s), PathId{0}, 0, 0, o);
  const double a = turned.x_out.dot(straight.x_out);
  const double b = turned.x_out.dot(straight.y_out);
  CHECK(std::abs(a - std::cos(kPi / 6.0)) <= 1e-14);  // content guard: the frame turned
  for (std::size_t r = 0; r < straight.status.size(); ++r) {
    INFO("point " << r);
    const Cx* j = &straight.j[4 * r];
    const Cx* k = &turned.j[4 * r];
    CHECK(std::abs(k[0] - (a * j[0] + b * j[2])) <= 1e-14);
    CHECK(std::abs(k[1] - (a * j[1] + b * j[3])) <= 1e-14);
    CHECK(std::abs(k[2] - (-b * j[0] + a * j[2])) <= 1e-14);
    CHECK(std::abs(k[3] - (-b * j[1] + a * j[3])) <= 1e-14);
  }
  CHECK(same_bits(turned.diattenuation, straight.diattenuation));
  CHECK(same_bits(turned.retardance, straight.retardance));
}

TEST_CASE("Jones pupil: phase 0 at the chief ray, and the option to remove the global phase",
          "[jones]") {
  // ADR 0034, point 5 (F2): W = 0 at the chief ray, so J(0, 0) has no OPL phase; the Fresnel
  // phases stay. remove_global_phase multiplies every J by exp(-i psi), psi = arg(det J(0, 0))/2:
  // then det J(0, 0) is real and positive (1e-15 relative).
  const CompiledSystem cs = compiled(convergent_mirror());
  JonesPupilOptions o;
  o.sampling = rtt::trace::SinglePupilPoint{0.0, 0.0};
  const JonesPupil plain = rtt::analysis::jones_pupil(cs, PathId{0}, 0, 1, o);
  REQUIRE(plain.status[0] == RayStatus::Alive);
  CHECK(plain.w[0] == 0.0);
  const Cx det = plain.j[0] * plain.j[3] - plain.j[1] * plain.j[2];
  CHECK(std::abs(std::arg(det)) > 1e-2);  // content guard: the metal gives a phase
  o.remove_global_phase = true;
  const JonesPupil removed = rtt::analysis::jones_pupil(cs, PathId{0}, 0, 1, o);
  CHECK(removed.global_phase == std::arg(det) / 2.0);
  const Cx det_removed = removed.j[0] * removed.j[3] - removed.j[1] * removed.j[2];
  CHECK(std::abs(det_removed.imag()) <= 1e-15 * std::abs(det_removed));
  CHECK(det_removed.real() > 0.0);
  for (std::size_t i = 0; i < 4; ++i) {
    CHECK(std::abs(removed.j[i] - plain.j[i] * std::polar(1.0, -removed.global_phase)) <= 1e-15);
  }
}

TEST_CASE("Jones pupil: the global phase of a polarizer comes from the largest element",
          "[jones]") {
  // An ideal polarizer at the chief ray makes J(0, 0) singular (|det| <= 1e-12 ||J||_F^2): psi is
  // the argument of its largest element, which is then real and positive.
  System s = metal_paraboloid();
  Surface pol;
  pol.id = SurfaceId("POL");
  pol.interaction = rtt::model::IdealPolarizer{{1.0, 1.0, 0.0}, 0.0};
  s.root.children.insert(
      s.root.children.begin() + 1,
      {Element{"POL", ElementKind::ThinElement, Pose::along_z(-40.0), std::nullopt, {pol}}});
  s.paths = {{"polarized",
              false,
              {{SurfaceId("STO"), EventKind::Transmit, 0},
               {SurfaceId("POL"), EventKind::Transmit, 0},
               {SurfaceId("M"), EventKind::Reflect, 0},
               {SurfaceId("F"), EventKind::Transmit, 0}}}};
  const CompiledSystem cs = compiled(s);
  JonesPupilOptions o;
  o.sampling = rtt::trace::SinglePupilPoint{0.0, 0.0};
  o.remove_global_phase = true;
  const JonesPupil jp = rtt::analysis::jones_pupil(cs, PathId{0}, 0, 0, o);
  REQUIRE(jp.status[0] == RayStatus::Alive);
  const Cx det = jp.j[0] * jp.j[3] - jp.j[1] * jp.j[2];
  double frobenius2 = 0.0;
  std::size_t largest = 0;
  for (std::size_t i = 0; i < 4; ++i) {
    frobenius2 += std::norm(jp.j[i]);
    if (std::abs(jp.j[i]) > std::abs(jp.j[largest])) largest = i;
  }
  CHECK(std::abs(det) <= 1e-12 * frobenius2);  // content guard: the singular case
  CHECK(std::abs(jp.j[largest].imag()) <= 1e-15 * std::abs(jp.j[largest]));
  CHECK(jp.j[largest].real() > 0.0);
}

TEST_CASE("Jones pupil: transmission is the weight, W and status are those of opd_points",
          "[jones]") {
  // The same rays as opd_points() (traced with Q requested, bitwise the same otherwise, #207):
  // W and status bit for bit, the transmission the weight of the ray.
  const CompiledSystem cs = compiled(convergent_mirror());
  const rtt::trace::HexapolarPupil hex{3};
  JonesPupilOptions o;
  o.sampling = hex;
  const JonesPupil jp = rtt::analysis::jones_pupil(cs, PathId{0}, 2, 1, o);
  const auto opd = rtt::analysis::opd_points(cs, PathId{0}, 2, 1, hex, rtt::analysis::OpdOptions{});
  const std::uint16_t fields[] = {2};
  rtt::trace::RayBatch rays =
      rtt::trace::make_rays(cs, PathId{0}, fields, 1, hex, rtt::trace::Aiming::Real);
  [[maybe_unused]] const auto stats = rtt::trace::SequentialTracer().trace(cs, PathId{0}, rays);
  REQUIRE(jp.status.size() == opd.points.size());
  for (std::size_t i = 0; i < jp.status.size(); ++i) {
    CHECK(jp.px[i] == opd.points[i].px);
    CHECK(jp.py[i] == opd.points[i].py);
    CHECK(jp.status[i] == opd.points[i].status);
    CHECK(jp.w[i] == opd.points[i].w);
    if (jp.status[i] == RayStatus::Alive) CHECK(jp.transmission[i] == rays.weight()[i]);
  }
  CHECK(jp.losses.launched == opd.points.size());
}

TEST_CASE("Jones pupil: the dipole basis agrees with the double-pole basis at the chief ray",
          "[jones][basis]") {
  // Both output bases are the chief frame (x'_c, y'_c) at k = k_c (Lam, App. B); J(0, 0) of the
  // reference wavelength therefore agrees to rounding (1e-15), and elsewhere they differ.
  const CompiledSystem cs = compiled(metal_paraboloid());
  JonesPupilOptions o;
  o.sampling = rtt::trace::HexapolarPupil{2};
  const JonesPupil dp = rtt::analysis::jones_pupil(cs, PathId{0}, 0, 0, o);
  o.basis = JonesBasis::Dipole;
  const JonesPupil di = rtt::analysis::jones_pupil(cs, PathId{0}, 0, 0, o);
  CHECK(di.basis == JonesBasis::Dipole);
  for (std::size_t i = 0; i < 4; ++i) CHECK(std::abs(di.j[i] - dp.j[i]) <= 1e-15);
  double differ = 0.0;
  for (std::size_t i = 4; i < dp.j.size(); ++i)
    differ = std::max(differ, std::abs(di.j[i] - dp.j[i]));
  CHECK(differ > 1e-6);                            // content guard
  CHECK(same_bits(di.retardance, dp.retardance));  // the maps use the input basis only
}

TEST_CASE("Jones pupil: a path with crystal modes is an AnalysisError at /paths/i",
          "[jones][errors]") {
  // ADR 0034, H6: rejected before any ray is traced, location of the path, no code.
  const CompiledSystem cs = compiled(load("m4/calcite_walkoff.rtt.json"));
  try {
    static_cast<void>(rtt::analysis::jones_pupil(cs, PathId{0}, 0, 0));
    FAIL("no AnalysisError");
  } catch (const AnalysisError& e) {
    CHECK(e.location() == std::optional<std::string>("/paths/0"));
    CHECK_FALSE(e.surface().has_value());
  }
}

TEST_CASE("Jones pupil: bitwise the same at every thread count", "[jones]") {
  const CompiledSystem cs = compiled(convergent_mirror());
  JonesPupilOptions o;
  o.sampling = rtt::trace::GridPupil{15};
  o.remove_global_phase = true;
  const JonesPupil all = rtt::analysis::jones_pupil(cs, PathId{0}, 2, 1, o);
  for (const int threads : {1, 4}) {
    INFO(threads << " threads");
    const tbb::global_control limit(tbb::global_control::max_allowed_parallelism,
                                    static_cast<std::size_t>(threads));
    const JonesPupil jp = rtt::analysis::jones_pupil(cs, PathId{0}, 2, 1, o);
    CHECK(same_bits(jp.j, all.j));
    CHECK(same_bits(jp.w, all.w));
    CHECK(same_bits(jp.diattenuation, all.diattenuation));
    CHECK(same_bits(jp.retardance, all.retardance));
    CHECK(same_bits(jp.retardance_orientation, all.retardance_orientation));
  }
}
