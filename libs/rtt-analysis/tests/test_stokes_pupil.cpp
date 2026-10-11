// Stokes parameters and the polarization map of the pupil (#211; ADR 0034, point 6), from the
// Jones pupil (#208). Stokes definition of the project (polar::stokes; Lam, Table 2.1, p. 69):
// s0 = |a|^2 + |b|^2, s1 = |a|^2 - |b|^2, s2 = 2 Re(a* b), s3 = 2 Im(a b*) of E_out = J E_in in
// the output basis; power scale S = weight / (||J||_F^2 / 2) Stokes(J E_in) (ADR 0021, point 2,
// coordinator's decision for #211).

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/analysis/jones.hpp"
#include "rtt/analysis/stokes.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/types.hpp"
#include "rtt/model/model.hpp"
#include "rtt/polar/fresnel.hpp"
#include "rtt/trace/ray_batch.hpp"

using Cx = std::complex<double>;
using rtt::analysis::JonesPupil;
using rtt::analysis::JonesPupilOptions;
using rtt::analysis::StokesPupil;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
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
constexpr std::size_t kChief = 12;
const Cx kMetal(1.2, 7.26);  // invented test metal of M3 (m3/metal_periscope)

CompiledSystem compiled(const System& s) {
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

/// A system in vacuum at 0.55 um with the given elements and path.
System bench(std::vector<rtt::model::Node> children, std::vector<rtt::model::Event> events) {
  System s;
  s.name = "stokes bench";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.55, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(6.0)};
  s.fields = {rtt::model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  s.root.children = std::move(children);
  s.paths = {{"path", false, std::move(events)}};
  return s;
}

Surface surface(const char* id) {
  Surface s;
  s.id = SurfaceId(id);
  return s;
}

/// 5 x 5 rays along +z from z = 0 at x, y in {-3, -1.5, 0, 1.5, 3}, labelled (x, y) / 3; the
/// chief ray (0, 0) has index kChief.
rtt::trace::RayBatch parallel_rays() {
  rtt::trace::RayBatch rays(25);
  for (std::size_t i = 0; i < 25; ++i) {
    const double x = -3.0 + 1.5 * static_cast<double>(i % 5);
    const double y = -3.0 + 1.5 * static_cast<double>(i / 5);
    rays.pos_x()[i] = x;
    rays.pos_y()[i] = y;
    rays.pupil_x()[i] = x / 3.0;
    rays.pupil_y()[i] = y / 3.0;
  }
  return rays;
}

/// A quarter-wave retarder with its fast axis at 45 deg (local (1, 1, 0) / sqrt(2)) at z = 10 in
/// the collimated beam, an ideal lens f = 50 at z = 20 (ADR 0031: P = Q = R(k_in -> k_out),
/// lossless) and the detector in the focus at z = 70: the rays converge with different k.
System quarter_wave_lens() {
  Surface qwp = surface("QWP");
  const double h = std::sqrt(0.5);
  qwp.interaction = rtt::model::IdealRetarder{{h, h, 0.0}, 0.25};
  Surface lens = surface("IL");
  lens.interaction = rtt::model::IdealLens{Param(50.0), std::nullopt};
  return bench(
      {{Element{"QWP", ElementKind::ThinElement, Pose::along_z(10.0), std::nullopt, {qwp}}},
       {Element{"IL", ElementKind::ThinElement, Pose::along_z(20.0), std::nullopt, {lens}}},
       {Element{
           "image", ElementKind::Detector, Pose::along_z(70.0), std::nullopt, {surface("IMG")}}}},
      {{SurfaceId("QWP"), EventKind::Transmit, 0},
       {SurfaceId("IL"), EventKind::Transmit, 0},
       {SurfaceId("IMG"), EventKind::Transmit, 0}});
}

/// The 45 deg fold of #208 (test_jones_pupil.cpp): metal mirror at z = 10 tilted 45 deg about x,
/// detector at (0, 20, 10) facing +y. Collimated: J = diag(r_s, r_p) times a phase at every
/// point, in the s/p basis (x, p).
System fold45() {
  Pose mirror = Pose::along_z(10.0);
  mirror.rotation_deg[0] = Param(45.0);
  Pose detector = Pose::along_z(10.0);
  detector.position[1] = Param(20.0);
  detector.rotation_deg[0] = Param(-90.0);
  return bench(
      {{Element{"M", ElementKind::Mirror, mirror, std::string("CONST:1.2,7.26"), {surface("M")}}},
       {Element{"image", ElementKind::Detector, detector, std::nullopt, {surface("IMG")}}}},
      {{SurfaceId("M"), EventKind::Reflect, 0}, {SurfaceId("IMG"), EventKind::Transmit, 0}});
}

JonesPupil jones(const System& s, JonesPupilOptions options = {}) {
  return rtt::analysis::jones_pupil(compiled(s), PathId{0}, parallel_rays(), kChief, options);
}

}  // namespace

TEST_CASE("quarter-wave plate at 45 deg: |S3| = S0 over the whole pupil", "[stokes]") {
  // Acceptance of #211: x-polarized light through a lambda/4 plate at 45 deg is circular, and
  // the lossless ideal lens keeps it so at every pupil point (the basis rotation R(k_c -> k)
  // against R(z -> k) is a rotation about k, which leaves circular states unchanged). Scale (b):
  // S0 = weight |J E|^2 / (||J||_F^2 / 2) = weight = 1 for the lossless path (coordinator).
  // Tolerances 1e-12: rounding of a few products of unit vectors.
  const JonesPupil jp = jones(quarter_wave_lens());
  const StokesPupil st = rtt::analysis::stokes_pupil(jp, Cx(1.0), Cx(0.0));
  REQUIRE(st.s0.size() == 25);
  for (std::size_t r = 0; r < 25; ++r) {
    INFO("point " << r);
    REQUIRE(st.status[r] == RayStatus::Alive);
    CHECK(std::abs(std::abs(st.s3[r]) - st.s0[r]) <= 1e-12);
    CHECK(std::abs(st.s1[r]) <= 1e-12);
    CHECK(std::abs(st.s2[r]) <= 1e-12);
    CHECK(std::abs(st.s0[r] - jp.transmission[r]) <= 1e-12);
    CHECK(std::abs(st.s0[r] - 1.0) <= 1e-12);
    CHECK(std::abs(std::abs(st.ellipticity[r]) - kPi / 4.0) <= 1e-6);  // |chi| = pi/4
  }
  // Content guard: the rays converge (different k), and light along the fast axis stays linear.
  CHECK(std::abs(jp.k_out.z() - 1.0) <= 1e-15);
  const StokesPupil along = rtt::analysis::stokes_pupil(jp, Cx(1.0), Cx(1.0));
  for (std::size_t r = 0; r < 25; ++r) {
    INFO("point " << r);
    CHECK(std::abs(along.s3[r]) <= 1e-12);
    CHECK(along.s2[r] > 0.999);
  }
}

TEST_CASE("45 deg fold: the Stokes map is that of diag(r_s, r_p) E_in", "[stokes]") {
  // Acceptance of #211: J = diag(r_s, r_p) times a phase (#208), so E_out = (r_s e1, r_p e2) up
  // to the phase, which cancels in the Stokes parameters; weight = (|r_s|^2 + |r_p|^2) / 2 with
  // no other factor, so the scale is 1. Several input states; 1e-10 as in #211.
  const JonesPupil jp = jones(fold45());
  const auto f = rtt::polar::fresnel_at_angle<double>(Cx(1.0), kMetal, std::cos(kPi / 4.0));
  const double h = std::sqrt(0.5);
  const std::array<std::array<Cx, 2>, 4> inputs{
      {{Cx(1.0), Cx(0.0)}, {Cx(0.0), Cx(1.0)}, {Cx(h), Cx(h)}, {Cx(h), Cx(0.0, h)}}};
  double s3_max = 0.0;
  for (const auto& e : inputs) {
    const StokesPupil st = rtt::analysis::stokes_pupil(jp, e[0], e[1]);
    const Cx a = f.rs * e[0];
    const Cx b = f.rp * e[1];
    const double s0 = std::norm(a) + std::norm(b);
    const double s1 = std::norm(a) - std::norm(b);
    const double s2 = 2.0 * std::real(std::conj(a) * b);
    const double s3 = 2.0 * std::imag(a * std::conj(b));
    const auto ellipse = rtt::analysis::polarization_ellipse(a, b);
    REQUIRE(st.s0.size() == 25);
    for (std::size_t r = 0; r < 25; ++r) {
      INFO("input (" << e[0] << ", " << e[1] << "), point " << r);
      CHECK(std::abs(st.s0[r] - s0) <= 1e-10);
      CHECK(std::abs(st.s1[r] - s1) <= 1e-10);
      CHECK(std::abs(st.s2[r] - s2) <= 1e-10);
      CHECK(std::abs(st.s3[r] - s3) <= 1e-10);
      if (s0 * s0 - s1 * s1 > 1e-6) {  // the orientation is defined
        CHECK(std::abs(std::remainder(st.orientation[r] - ellipse.orientation, kPi)) <= 1e-8);
      }
      CHECK(std::abs(st.ellipticity[r] - ellipse.ellipticity) <= 1e-8);
    }
    s3_max = std::max(s3_max, std::abs(s3));
  }
  CHECK(s3_max > 0.1);  // content guard: the metal makes the 45 deg input elliptical
}

TEST_CASE("handedness: right circular stays S3 = +S0 in free space", "[stokes]") {
  // (i e1 + e2) / sqrt(2) is right circular (Lam, Table 2.1; docs/architecture.md), s3 = +s0,
  // chi = +pi/4. A detector alone: J = identity. The input is normalised by stokes_pupil.
  const JonesPupil jp = jones(bench(
      {{Element{
          "image", ElementKind::Detector, Pose::along_z(10.0), std::nullopt, {surface("IMG")}}}},
      {{SurfaceId("IMG"), EventKind::Transmit, 0}}));
  const StokesPupil st = rtt::analysis::stokes_pupil(jp, Cx(0.0, 3.0), Cx(3.0));
  CHECK(std::abs(st.e1 - Cx(0.0, std::sqrt(0.5))) <= 1e-15);
  CHECK(std::abs(st.e2 - Cx(std::sqrt(0.5))) <= 1e-15);
  for (std::size_t r = 0; r < 25; ++r) {
    INFO("point " << r);
    CHECK(std::abs(st.s3[r] - 1.0) <= 1e-15);
    CHECK(std::abs(st.s0[r] - 1.0) <= 1e-15);
    CHECK(std::abs(st.ellipticity[r] - kPi / 4.0) <= 1e-12);
  }
}

TEST_CASE("power scale with absorption: S0(x) + S0(y) = 2 weight", "[stokes]") {
  // Scale (b): S = s Stokes(J E) with s = weight / (||J||_F^2 / 2) the polarization-independent
  // factor (here the volume absorption of the glass, ADR 0021). For two orthogonal inputs
  // S0(x) + S0(y) = s ||J||_F^2 = 2 weight, independent of the polarization of the path.
  Surface front = surface("B.S1");
  Surface back = surface("B.S2");
  back.pose = Pose::along_z(1.0);
  const System s = bench(
      {{Element{"block",
                ElementKind::Plate,
                Pose::along_z(10.0),
                std::string("CONST:1.5,1e-5"),
                {front, back}}}},
      {{SurfaceId("B.S1"), EventKind::Refract, 0}, {SurfaceId("B.S2"), EventKind::Transmit, 0}});
  rtt::trace::RayBatch rays = parallel_rays();
  for (std::size_t i = 0; i < 25; ++i) {  // 30 deg in the y-z plane
    rays.dir_y()[i] = 0.5;
    rays.dir_z()[i] = std::sqrt(0.75);
  }
  const JonesPupil jp =
      rtt::analysis::jones_pupil(compiled(s), PathId{0}, rays, kChief, JonesPupilOptions{});
  const StokesPupil x = rtt::analysis::stokes_pupil(jp, Cx(1.0), Cx(0.0));
  const StokesPupil y = rtt::analysis::stokes_pupil(jp, Cx(0.0), Cx(1.0));
  for (std::size_t r = 0; r < 25; ++r) {
    INFO("point " << r);
    REQUIRE(jp.status[r] == RayStatus::Alive);
    CHECK(std::abs(x.s0[r] + y.s0[r] - 2.0 * jp.transmission[r]) <= 1e-12);
    double frobenius2 = 0.0;
    for (std::size_t e = 0; e < 4; ++e) frobenius2 += std::norm(jp.j[4 * r + e]);
    CHECK(jp.transmission[r] < 0.95 * 0.5 * frobenius2);  // content guard: s < 1 (absorption)
    CHECK(std::abs(x.s0[r] - y.s0[r]) > 1e-3);            // and the path polarizes (Fresnel)
  }
}

TEST_CASE("Stokes map: points without J, copied columns, input errors", "[stokes]") {
  JonesPupilOptions o;
  o.lost_warning_fraction = 0.0;
  rtt::trace::RayBatch rays = parallel_rays();
  rays.dir_z()[0] = -1.0;  // runs away from the mirror: lost
  const JonesPupil jp = rtt::analysis::jones_pupil(compiled(fold45()), PathId{0}, rays, kChief, o);
  const StokesPupil st = rtt::analysis::stokes_pupil(jp, Cx(1.0), Cx(0.0));
  REQUIRE(st.s0.size() == 25);
  CHECK(st.status == jp.status);
  CHECK(st.px == jp.px);
  CHECK(st.py == jp.py);
  CHECK(st.status[0] != RayStatus::Alive);
  for (const std::vector<double>* v :
       {&st.s0, &st.s1, &st.s2, &st.s3, &st.orientation, &st.ellipticity}) {
    CHECK(std::isnan((*v)[0]));
    CHECK(std::isfinite((*v)[1]));
  }
  REQUIRE(st.warnings.size() == jp.warnings.size());
  REQUIRE(!st.warnings.empty());
  CHECK(st.warnings[0].code == "rays.lost");
  CHECK((st.x_out - jp.x_out).norm() == 0.0);
  CHECK(st.path.index == jp.path.index);
  CHECK(st.wavelength_um == jp.wavelength_um);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CHECK_THROWS_AS(rtt::analysis::stokes_pupil(jp, Cx(0.0), Cx(0.0)), std::invalid_argument);
  CHECK_THROWS_AS(rtt::analysis::stokes_pupil(jp, Cx(nan), Cx(1.0)), std::invalid_argument);
  JonesPupilOptions field;
  field.field_amplitude = true;
  const JonesPupil amplitudes = jones(fold45(), field);
  REQUIRE(amplitudes.field_amplitude);
  CHECK_THROWS_AS(rtt::analysis::stokes_pupil(amplitudes, Cx(1.0), Cx(0.0)), std::invalid_argument);
}
