// Q-con and Q-bfs in the paraxial model (#314, ADR 0043, point 3 and acceptance 3 and 4): Q-con
// represents an even asphere exactly ([F10a], Eqs. (5.6), (5.7)), with the first order bit for
// bit and A4 within its rounding; the first order of Q-bfs follows the vertex curvature
// c + 2 S_0 / rho^2 ([F10b], Eq. (3.16)).

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/geom/forbes.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/prescription.hpp"
#include "rtt/paraxial/seidel.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::Element;
using rtt::model::Param;
using rtt::model::System;
using rtt::paraxial::ChromaticPair;
using rtt::paraxial::FirstOrder;
using rtt::paraxial::ParaxialError;
using rtt::paraxial::Prescription;
using rtt::paraxial::Seidel;

namespace {

/// m1/singlet_const with `base` at L1.S1 (/root/children/1/surfaces/0).
System singlet(const rtt::model::BaseShape& base) {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  std::get<Element>(s.root.children[1].value).surfaces[0].shape.base = base;
  return s;
}

CompiledSystem compiled(const System& s) {
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

void same(const std::optional<rtt::paraxial::Pupil>& a,
          const std::optional<rtt::paraxial::Pupil>& b) {
  REQUIRE(a.has_value() == b.has_value());
  if (!a || !b) return;
  CHECK(a->z == b->z);
  CHECK(a->diameter == b->diameter);
}

void same(const FirstOrder& a, const FirstOrder& b) {
  CHECK(a.object_index == b.object_index);
  CHECK(a.image_index == b.image_index);
  CHECK(a.image_direction == b.image_direction);
  CHECK(a.power == b.power);
  CHECK(a.efl == b.efl);
  CHECK(a.front_focal_length == b.front_focal_length);
  CHECK(a.rear_focal_length == b.rear_focal_length);
  CHECK(a.ffl == b.ffl);
  CHECK(a.bfl == b.bfl);
  CHECK(a.front_focal_z == b.front_focal_z);
  CHECK(a.rear_focal_z == b.rear_focal_z);
  CHECK(a.front_principal_z == b.front_principal_z);
  CHECK(a.rear_principal_z == b.rear_principal_z);
  CHECK(a.image_z == b.image_z);
  CHECK(a.lateral_magnification == b.lateral_magnification);
  CHECK(a.angular_magnification == b.angular_magnification);
  same(a.entrance_pupil, b.entrance_pupil);
  same(a.exit_pupil, b.exit_pupil);
}

void same(const std::optional<rtt::paraxial::PrescriptionRay>& a,
          const std::optional<rtt::paraxial::PrescriptionRay>& b) {
  REQUIRE(a.has_value() == b.has_value());
  if (!a || !b) return;
  CHECK(a->y == b->y);
  CHECK(a->u == b->u);
  CHECK(a->i == b->i);
}

void same(const Prescription& a, const Prescription& b) {
  REQUIRE(a.surfaces.size() == b.surfaces.size());
  for (std::size_t i = 0; i < a.surfaces.size(); ++i) {
    INFO("event " << i);
    CHECK(a.surfaces[i].surface == b.surfaces[i].surface);
    CHECK(a.surfaces[i].z == b.surfaces[i].z);
    CHECK(a.surfaces[i].n == b.surfaces[i].n);
    same(a.surfaces[i].marginal, b.surfaces[i].marginal);
    same(a.surfaces[i].chief, b.surfaces[i].chief);
    CHECK(a.surfaces[i].lagrange == b.surfaces[i].lagrange);
  }
  CHECK(a.total_track == b.total_track);
  CHECK(a.object_distance == b.object_distance);
  CHECK(a.lagrange_invariant == b.lagrange_invariant);
  same(a.first_order, b.first_order);
}

}  // namespace

TEST_CASE("Q-con represents an even asphere exactly", "[paraxial][forbes]") {
  // ADR 0043, acceptance 3. [F10a], Eqs. (5.6), (5.7): u^4 sum s_m Q_m(u^2) = sum A_(2m+4)
  // r^(2m+4) gives sum s_m Q_m(x) = sum t_m x^m with t_m = A_(2m+4) rho^(2m+4); solved here by
  // back substitution with the monomial coefficients of Q_0 ... Q_3 ([F07], Eq. (9)). R = 64 is
  // dyadic, so curvature = 1.0 / 64 is exactly the double that compile forms for the even asphere.
  const double rho = 12.7;
  const std::vector<double> a_even = {-2e-6, 3e-9, -1e-12, 2e-15};  // A4, A6, A8, A10
  const std::array<std::array<double, 4>, 4> q = {{{1.0, 0.0, 0.0, 0.0},
                                                   {-5.0, 6.0, 0.0, 0.0},
                                                   {15.0, -42.0, 28.0, 0.0},
                                                   {-35.0, 168.0, -252.0, 120.0}}};
  std::array<double, 4> t{};
  for (std::size_t m = 0; m < 4; ++m)
    t[m] = a_even[m] * std::pow(rho, static_cast<double>(2 * m + 4));
  std::array<double, 4> s{};
  for (std::size_t m = 4; m-- > 0;) {
    double rest = t[m];
    for (std::size_t n = m + 1; n < 4; ++n) rest -= s[n] * q[n][m];
    s[m] = rest / q[m][m];
  }
  for (const double k : {0.0, -0.8}) {
    INFO("k = " << k);
    std::vector<Param> even_coefficients;
    for (const double v : a_even) even_coefficients.emplace_back(v);
    std::vector<Param> qcon_coefficients;
    for (const double v : s) qcon_coefficients.emplace_back(v);
    const System even = singlet(rtt::model::EvenAsphere{Param(64.0), Param(k), even_coefficients});
    const System qcon =
        singlet(rtt::model::QCon{Param(1.0 / 64.0), Param(k), Param(rho), qcon_coefficients});
    const CompiledSystem a = compiled(even);
    const CompiledSystem b = compiled(qcon);
    for (std::uint16_t w = 0; w < 3; ++w) {
      INFO("wavelength " << w);
      same(rtt::paraxial::first_order(b, PathId{0}, w),
           rtt::paraxial::first_order(a, PathId{0}, w));
      same(rtt::paraxial::prescription(b, PathId{0}, w),
           rtt::paraxial::prescription(a, PathId{0}, w));
    }
    // A4 within (N + 2) eps sum |s_m| C(m + 4, 4) / rho^4, N = 4 terms.
    const double a4 =
        std::get<rtt::geom::QCon<double>>(b.surfaces()[1].shape).vertex_expansion().a4;
    double scale = 0.0;
    const std::array<double, 4> binomial = {1.0, 5.0, 15.0, 35.0};
    for (std::size_t m = 0; m < 4; ++m) scale += std::abs(s[m]) * binomial[m];
    const double rho4 = rho * rho * rho * rho;
    CHECK(std::abs(a4 - a_even[0]) <= 6.0 * std::numeric_limits<double>::epsilon() * scale / rho4);
    // Content guard: the cap enters the Seidel sums of both.
    const Seidel sa = rtt::paraxial::seidel(a, PathId{0}, 1, ChromaticPair{0, 2});
    const Seidel sb = rtt::paraxial::seidel(b, PathId{0}, 1, ChromaticPair{0, 2});
    CHECK(std::abs(sb.sum.s1 - sa.sum.s1) <= 1e-12 * std::abs(sa.sum.s1));
    CHECK(sa.sum.s1 != 0.0);
  }
}

TEST_CASE("Q-bfs: the first order follows c + 2 S_0 / rho^2", "[paraxial][forbes]") {
  // ADR 0043, acceptance 4 (content guard): a plano-convex singlet with the Q-bfs front and the
  // plane back has the power (n - 1) c_v; S_0 = a_0 Q_0(0) + a_1 Q_1(0) = a_0 + 13 a_1 /
  // sqrt(19) ([F07], Eq. (14)).
  const double c = 1.0 / 60.0;
  const double rho = 12.7;
  const double a0 = 0.1;
  const double a1 = -0.02;
  System s = singlet(rtt::model::QBfs{Param(c), Param(-0.3), Param(rho), {Param(a0), Param(a1)}});
  s.environment.medium = "VACUUM";  // n = 1 outside, so the power is (n - 1) c_v
  std::get<Element>(s.root.children[1].value).surfaces[1].shape.base = rtt::model::Plane{};
  const CompiledSystem cs = compiled(s);
  const double s0 = a0 + 13.0 * a1 / std::sqrt(19.0);
  const double cv = c + 2.0 * s0 / (rho * rho);
  const FirstOrder f = rtt::paraxial::first_order(cs, PathId{0}, 1);
  REQUIRE(f.efl.has_value());
  CHECK(std::abs(1.0 / *f.efl - (1.5168 - 1.0) * cv) <= 1e-12 * (1.5168 - 1.0) * cv);
  CHECK(std::abs(cv - c) > 1e-4);  // S_0 changes the focal length
  CHECK_NOTHROW((void)rtt::paraxial::seidel(cs, PathId{0}, 1));
}
