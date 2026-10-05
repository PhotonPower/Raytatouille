// Ideal polarizers and retarders (#60). Source: W.-S. T. Lam, "Anisotropic ray trace",
// dissertation, University of Arizona (docs/quellen.md); conventions in rtt/polar/ideal.hpp and
// docs/architecture.md. The reference cases test convention-free quantities (Malus, circularity,
// rotation angle) plus one handedness check against Lam's example.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>

#include "rtt/polar/ideal.hpp"
#include "rtt/polar/prt_analysis.hpp"

using Cx = std::complex<double>;
using rtt::math::CMat3;
using rtt::math::Vec3;
using CVec3 = Eigen::Vector3cd;

namespace {

constexpr double kPi = std::numbers::pi;

/// Deterministic uniform numbers in [0, 1) ((x >> 11) * 2^-53, as for the random pupil).
class Uniform {
 public:
  explicit Uniform(std::uint64_t seed) : engine_(seed) {}
  double operator()() { return static_cast<double>(engine_() >> 11) * 0x1.0p-53; }

 private:
  std::mt19937_64 engine_;
};

Vec3 random_unit(Uniform& u) {
  const double z = 2.0 * u() - 1.0;
  const double phi = 2.0 * kPi * u();
  const double r = std::sqrt(1.0 - z * z);
  return {r * std::cos(phi), r * std::sin(phi), z};
}

/// Linear polarization angle from e1 towards e2, in (-pi/2, pi/2].
double linear_angle(const rtt::polar::Stokes<double>& s) {
  return 0.5 * std::atan2(s.s2, s.s1);
}

/// Difference of two polarization angles modulo pi, in [0, pi/2].
double angle_distance(double a, double b) {
  const double d = std::remainder(a - b, kPi);
  return std::abs(d);
}

double max_abs(const CMat3& m) {
  return m.cwiseAbs().maxCoeff();
}

}  // namespace

TEST_CASE("embedded Jones matrix: P e1, P e2 and P k", "[ideal]") {
  // P = sum_ij J_ij e_i e_j^T + k k^T with e1 = axis projected perpendicular to k, e2 = k x e1.
  // Tolerance 1e-14.
  Uniform u(1);
  for (int i = 0; i < 50; ++i) {
    const Vec3 k = random_unit(u);
    const Vec3 axis = random_unit(u);
    rtt::polar::Jones2T<double> j;
    j << Cx(u(), u()), Cx(u(), -u()), Cx(-u(), u()), Cx(u(), u());
    const CMat3 p = rtt::polar::embed_jones(j, axis, k);
    const Vec3 e1 = (axis - axis.dot(k) * k).normalized();
    const Vec3 e2 = k.cross(e1);
    const CVec3 col1 = j(0, 0) * e1.cast<Cx>() + j(1, 0) * e2.cast<Cx>();
    const CVec3 col2 = j(0, 1) * e1.cast<Cx>() + j(1, 1) * e2.cast<Cx>();
    REQUIRE((p * e1.cast<Cx>() - col1).norm() <= 1e-14);
    REQUIRE((p * e2.cast<Cx>() - col2).norm() <= 1e-14);
    REQUIRE((p * k.cast<Cx>() - k.cast<Cx>()).norm() <= 1e-14);
  }
}

TEST_CASE("Malus: two polarizers transmit cos^2 theta", "[ideal]") {
  // Architecture, acceptance table: T = cos^2(theta) for perfect polarizers, tolerance 1e-12.
  // With extinction ratio epsilon (power) of the second polarizer and linear light along the
  // first axis, the field has components cos(theta) on the transmission axis and sin(theta) on
  // the crossed axis with amplitude sqrt(epsilon): T = cos^2 + epsilon sin^2. Random directions
  // (fixed seed); the second axis carries a component along k that the projection must remove.
  Uniform u(2);
  for (int i = 0; i < 100; ++i) {
    const Vec3 k = random_unit(u);
    const Vec3 a1 = random_unit(u);
    const Vec3 e1 = (a1 - a1.dot(k) * k).normalized();
    const Vec3 e2 = k.cross(e1);
    const double theta = kPi * u();
    const Vec3 a2 = std::cos(theta) * e1 + std::sin(theta) * e2 + (u() - 0.5) * k;
    const double epsilon = i % 2 == 0 ? 0.0 : 0.01 * u();
    const CMat3 p1 = rtt::polar::ideal_polarizer(a1, 0.0, k);
    const CMat3 p2 = rtt::polar::ideal_polarizer(a2, epsilon, k);
    CVec3 e(Cx(u(), u()), Cx(u(), u()), Cx(u(), u()));
    e -= k.cast<Cx>() * k.cast<Cx>().dot(e);
    const CVec3 after_first = p1 * e;
    const double t = (p2 * after_first).squaredNorm() / after_first.squaredNorm();
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    INFO("theta " << theta << ", epsilon " << epsilon);
    REQUIRE(std::abs(t - (c * c + epsilon * s * s)) <= 1e-12);
  }
}

TEST_CASE("quarter-wave plate at 45 deg: linear becomes circular", "[ideal]") {
  // Architecture, acceptance table: linear -> circular, tolerance 1e-12, convention-free:
  // |S3| = S0, S1 = S2 = 0 for linear light at 45 deg to the fast axis, random directions.
  Uniform u(3);
  for (int i = 0; i < 100; ++i) {
    const Vec3 k = random_unit(u);
    const Vec3 axis = random_unit(u);
    const Vec3 e1 = (axis - axis.dot(k) * k).normalized();
    const Vec3 e2 = k.cross(e1);
    const Vec3 fast = (e1 + e2) / std::sqrt(2.0);
    const CMat3 p = rtt::polar::linear_retarder(fast, 0.25, k);
    const CVec3 out = p * e1.cast<Cx>();
    const auto s = rtt::polar::stokes(out, e1, k);
    REQUIRE(std::abs(s.s0 - 1.0) <= 1e-12);
    REQUIRE(std::abs(std::abs(s.s3) - 1.0) <= 1e-12);
    REQUIRE(std::abs(s.s1) <= 1e-12);
    REQUIRE(std::abs(s.s2) <= 1e-12);
  }
}

TEST_CASE("handedness: Lam's quarter-wave example gives right circular light", "[ideal]") {
  // Lam, Fig. 1.1 (p. 48): x delayed by a quarter wave against y, input at 45 deg -> right
  // circular; Lam, Table 2.1 (p. 69): right circular is (i, 1, 0). Here slow axis x, fast axis
  // y, k = z: the symmetric phases give (e^(i pi/4), e^(-i pi/4))/sqrt(2) = e^(-i pi/4) (i, 1)/
  // sqrt(2), and S3 > 0 in the basis e1 = x. Tolerance 1e-12.
  const Vec3 k(0.0, 0.0, 1.0);
  const CMat3 p = rtt::polar::linear_retarder(Vec3(0.0, 1.0, 0.0), 0.25, k);
  const CVec3 in = CVec3(Cx(1.0), Cx(1.0), Cx(0.0)) / std::sqrt(2.0);
  const CVec3 out = p * in;
  const CVec3 expected =
      std::polar(1.0, -kPi / 4.0) * CVec3(Cx(0.0, 1.0), Cx(1.0), Cx(0.0)) / std::sqrt(2.0);
  REQUIRE((out - expected).norm() <= 1e-12);
  const auto s = rtt::polar::stokes(out, Vec3(1.0, 0.0, 0.0), k);
  REQUIRE(std::abs(s.s3 - 1.0) <= 1e-12);
}

TEST_CASE("half-wave plate at theta rotates linear polarization by 2 theta", "[ideal]") {
  // Architecture, acceptance table, tolerance 1e-12: linear light along e1 through a half-wave
  // plate with its fast axis at theta from e1 leaves linear at 2 theta (mirror image about the
  // fast axis), random directions. The angle is taken from the Stokes parameters (modulo pi).
  Uniform u(4);
  for (int i = 0; i < 100; ++i) {
    const Vec3 k = random_unit(u);
    const Vec3 axis = random_unit(u);
    const Vec3 e1 = (axis - axis.dot(k) * k).normalized();
    const Vec3 e2 = k.cross(e1);
    const double theta = kPi * (u() - 0.5);
    const Vec3 fast = std::cos(theta) * e1 + std::sin(theta) * e2;
    const CMat3 p = rtt::polar::linear_retarder(fast, 0.5, k);
    const auto s = rtt::polar::stokes(CVec3(p * e1.cast<Cx>()), e1, k);
    INFO("theta " << theta);
    REQUIRE(std::abs(s.s3) <= 1e-12);
    REQUIRE(angle_distance(linear_angle(s), 2.0 * theta) <= 1e-12);
  }
}

TEST_CASE("half-wave plate equals Lam's Jones matrix up to a global phase", "[ideal]") {
  // Lam, p. 123: half-wave plate with fast axis at alpha from x, J = [[cos 2a, sin 2a],
  // [sin 2a, -cos 2a]]. With the symmetric phase e^(-+i pi/2) ours is -i times that matrix.
  // k = z, e1 = x. Tolerance 1e-14.
  for (const double alpha : {0.0, 0.3, 1.1, -0.7}) {
    const CMat3 p = rtt::polar::linear_retarder(Vec3(std::cos(alpha), std::sin(alpha), 0.0), 0.5,
                                                Vec3(0.0, 0.0, 1.0));
    CMat3 lam = CMat3::Zero();
    lam(0, 0) = std::cos(2.0 * alpha);
    lam(0, 1) = std::sin(2.0 * alpha);
    lam(1, 0) = std::sin(2.0 * alpha);
    lam(1, 1) = -std::cos(2.0 * alpha);
    CMat3 block = p;
    block(2, 2) = 0.0;
    REQUIRE(max_abs(block - Cx(0.0, -1.0) * lam) <= 1e-14);
  }
}

TEST_CASE("consistency with #57: diattenuation and retardance of the ideal elements", "[ideal]") {
  // rtt/polar/prt_analysis.hpp: D = (1 - epsilon)/(1 + epsilon) with the axis of maximum
  // transmission along a; linear retarder: delta = 2 pi w in [0, pi] with fast axis f (for
  // w = 0.6, delta = 1.2 pi is reduced to 0.8 pi and the slow axis becomes the fast one);
  // circular retarder: delta with the left circular state as fast axis (right is delayed).
  // Tolerance 1e-12.
  Uniform u(5);
  for (int i = 0; i < 20; ++i) {
    const Vec3 k = random_unit(u);
    const Vec3 axis = random_unit(u);
    const Vec3 e1 = (axis - axis.dot(k) * k).normalized();
    const Vec3 e2 = k.cross(e1);
    for (const double epsilon : {0.0, 0.01, 0.3}) {
      const auto d = rtt::polar::diattenuation(rtt::polar::ideal_polarizer(axis, epsilon, k), k, k);
      REQUIRE(std::abs(d.value - (1.0 - epsilon) / (1.0 + epsilon)) <= 1e-12);
      REQUIRE(std::abs(std::abs(d.axis.dot(e1.cast<Cx>())) - 1.0) <= 1e-12);
    }
    for (const double w : {0.1, 0.25, 0.4, 0.6}) {
      const auto r = rtt::polar::retardance(rtt::polar::linear_retarder(axis, w, k), k);
      const double delta = 2.0 * kPi * w;
      const double reduced = delta <= kPi ? delta : 2.0 * kPi - delta;
      const Vec3 expected_fast = delta <= kPi ? e1 : e2;
      INFO("w " << w);
      REQUIRE(std::abs(r.value - reduced) <= 1e-12);
      REQUIRE(std::abs(std::abs(r.fast_axis.dot(expected_fast.cast<Cx>())) - 1.0) <= 1e-12);
    }
    const auto c = rtt::polar::retardance(rtt::polar::circular_retarder(0.2, k), k);
    const CVec3 left = (e1.cast<Cx>() + Cx(0.0, 1.0) * e2.cast<Cx>()) / std::sqrt(2.0);
    REQUIRE(std::abs(c.value - 0.4 * kPi) <= 1e-12);
    REQUIRE(std::abs(std::abs(c.fast_axis.dot(left)) - 1.0) <= 1e-12);
  }
}

TEST_CASE("circular retarder rotates linear polarization by +delta/2", "[ideal]") {
  // ideal.hpp: right circular delayed by delta = rotation by +delta/2 from e1 towards e2,
  // independent of the choice of e1. Tolerance 1e-12.
  Uniform u(6);
  for (int i = 0; i < 50; ++i) {
    const Vec3 k = random_unit(u);
    const Vec3 axis = random_unit(u);
    const Vec3 e1 = (axis - axis.dot(k) * k).normalized();
    const double w = u() * 0.5;
    const CMat3 p = rtt::polar::circular_retarder(w, k);
    const auto s = rtt::polar::stokes(CVec3(p * e1.cast<Cx>()), e1, k);
    REQUIRE(std::abs(s.s3) <= 1e-12);
    REQUIRE(angle_distance(linear_angle(s), kPi * w) <= 1e-12);
  }
}

TEST_CASE("ideal elements: P k = k, retarders unitary, axis sign and k component irrelevant",
          "[ideal]") {
  // Thin elements do not deflect; retarders conserve energy on the transverse plane; the axis
  // enters only through its projection and its sign does not matter. Tolerance 1e-14.
  Uniform u(7);
  for (int i = 0; i < 50; ++i) {
    const Vec3 k = random_unit(u);
    const Vec3 axis = random_unit(u);
    const CMat3 pol = rtt::polar::ideal_polarizer(axis, 0.04, k);
    const CMat3 ret = rtt::polar::linear_retarder(axis, 0.3, k);
    const CMat3 circ = rtt::polar::circular_retarder(0.3, k);
    for (const CMat3* p : {&pol, &ret, &circ}) {
      REQUIRE((*p * k.cast<Cx>() - k.cast<Cx>()).norm() <= 1e-14);
    }
    for (const CMat3* p : {&ret, &circ}) {
      REQUIRE(max_abs(p->adjoint() * *p - CMat3::Identity()) <= 1e-14);
    }
    const Vec3 shifted = -axis + 0.7 * k;
    REQUIRE(max_abs(rtt::polar::ideal_polarizer(shifted, 0.04, k) - pol) <= 1e-14);
    REQUIRE(max_abs(rtt::polar::linear_retarder(shifted, 0.3, k) - ret) <= 1e-14);
  }
}

TEST_CASE("checked variants reject invalid input", "[ideal]") {
  const Vec3 k(0.0, 0.0, 1.0);
  const Vec3 x(1.0, 0.0, 0.0);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  REQUIRE_NOTHROW(rtt::polar::ideal_polarizer_checked(x, 0.0, k));
  REQUIRE_THROWS_AS(rtt::polar::ideal_polarizer_checked(k, 0.0, k), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::polar::ideal_polarizer_checked(Vec3::Zero(), 0.0, k),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::polar::ideal_polarizer_checked(x, -0.1, k), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::polar::ideal_polarizer_checked(x, 1.5, k), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::polar::ideal_polarizer_checked(x, nan, k), std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::polar::ideal_polarizer_checked(x, 0.0, Vec3(0.0, 0.0, 2.0)),
                    std::invalid_argument);
  REQUIRE_NOTHROW(rtt::polar::linear_retarder_checked(x, 0.25, k));
  REQUIRE_THROWS_AS(rtt::polar::linear_retarder_checked(Vec3(0.0, 0.0, -3.0), 0.25, k),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(rtt::polar::linear_retarder_checked(x, nan, k), std::invalid_argument);
  REQUIRE_NOTHROW(rtt::polar::circular_retarder_checked(0.25, k));
  REQUIRE_THROWS_AS(rtt::polar::circular_retarder_checked(0.25, Vec3(nan, 0.0, 1.0)),
                    std::invalid_argument);
}
