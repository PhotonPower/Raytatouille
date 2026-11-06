// Real rays through Q-con and Q-bfs (#314, ADR 0043, acceptance 3 and 7): a ray through a Q-bfs
// with only a_0 agrees with an independent computation (bisection on the closed form of the sag,
// normal from its derivative by hand, vector form of Snell's law), and a Q-con that represents an
// even asphere ([F10a], Eqs. (5.6), (5.7)) gives the same intersection points.

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <string>
#include <variant>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::math::Vec3;
using rtt::model::Element;
using rtt::model::Param;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayPaths;
using rtt::trace::RayStatus;
using rtt::trace::SequentialTracer;

namespace {

const double kPi = std::acos(-1.0);

/// m1/singlet_const in vacuum (stop at z = 0, L1.S1 at z = 5, L1.S2 at z = 9, image at
/// z = 106.363) with the bases `front` and `back` and the index `material`.
System singlet(const rtt::model::BaseShape& front,
               const rtt::model::BaseShape& back,
               const std::string& material) {
  System s = rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  s.environment.medium = "VACUUM";
  Element& lens = std::get<Element>(s.root.children[1].value);
  lens.material = material;
  lens.surfaces[0].shape.base = front;
  lens.surfaces[1].shape.base = back;
  return s;
}

/// Traces rays from `starts` along `directions` and returns their paths; every ray arrives.
RayPaths traced(const System& s, const std::vector<Vec3>& starts, const std::vector<Vec3>& dirs) {
  const rtt::material::MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(s, lib);
  RayBatch rays(starts.size());
  for (std::size_t k = 0; k < rays.size(); ++k) {
    rays.pos_x()[k] = starts[k].x();
    rays.pos_y()[k] = starts[k].y();
    rays.pos_z()[k] = starts[k].z();
    rays.dir_x()[k] = dirs[k].x();
    rays.dir_y()[k] = dirs[k].y();
    rays.dir_z()[k] = dirs[k].z();
  }
  RayPaths paths;
  const auto stats = SequentialTracer().trace(cs, PathId{0}, rays, paths);
  REQUIRE(stats.count(RayStatus::Alive) == rays.size());
  REQUIRE(paths.slots == 5);  // start, stop, L1.S1, L1.S2, image
  return paths;
}

/// Vector form of Snell's law: d unit, n unit normal with n . d > 0, from index n1 into n2.
Vec3 refract(const Vec3& d, const Vec3& n, double n1, double n2) {
  const double mu = n1 / n2;
  const double ci = n.dot(d);
  const double ct = std::sqrt(1.0 - mu * mu * (1.0 - ci * ci));
  return mu * d + (ct - mu * ci) * n;
}

}  // namespace

TEST_CASE("a real ray through a Q-bfs with only a_0 against the closed form", "[trace][forbes]") {
  // ADR 0043, acceptance 7: with k = 0 and only a_0, Q_0 = 1 and the sag is
  // z = c r^2 / (1 + phi) + u^2 (1 - u^2) a_0 / phi, phi = sqrt(1 - c^2 r^2) ([F10b], Eq. (2.1)),
  // dz/dr = c r / phi + a_0 [(2u - 4u^3) / (rho phi) + u^2 (1 - u^2) c^2 r / phi^3].
  const double c = 1.0 / 30.0;
  const double rho = 12.0;
  const double a0 = 0.004;
  const System s = singlet(rtt::model::QBfs{Param(c), Param(0.0), Param(rho), {Param(a0)}},
                           rtt::model::Plane{}, "CONST:1.5168");
  const Vec3 o(0.0, 3.0, -1.0);
  const Vec3 d = Vec3(0.0, 0.05, 1.0).normalized();
  const RayPaths paths = traced(s, {o}, {d});
  const auto sag = [&](double r) {
    const double phi = std::sqrt(1.0 - c * c * r * r);
    const double u = r / rho;
    return c * r * r / (1.0 + phi) + u * u * (1.0 - u * u) * a0 / phi;
  };
  const auto slope = [&](double r) {
    const double phi = std::sqrt(1.0 - c * c * r * r);
    const double u = r / rho;
    return c * r / phi + a0 * ((2.0 * u - 4.0 * u * u * u) / (rho * phi) +
                               u * u * (1.0 - u * u) * c * c * r / (phi * phi * phi));
  };
  const Vec3 vertex(0.0, 0.0, 5.0);
  const auto f = [&](double t) {
    const Vec3 q = o + t * d - vertex;
    return q.z() - sag(std::hypot(q.x(), q.y()));
  };
  double lo = 0.0;
  double hi = 20.0;
  REQUIRE(f(lo) < 0.0);
  REQUIRE(f(hi) > 0.0);
  for (int i = 0; i < 200; ++i) {
    const double mid = lo + (hi - lo) / 2.0;
    if (mid == lo || mid == hi) break;
    (f(mid) < 0.0 ? lo : hi) = mid;
  }
  const Vec3 hit = o + lo * d;
  const Vec3 q = hit - vertex;
  const double r = std::hypot(q.x(), q.y());
  const double g = slope(r);
  const Vec3 normal = Vec3(-g * q.x() / r, -g * q.y() / r, 1.0).normalized();
  const Vec3 out = refract(d, normal, 1.0, 1.5168);
  CHECK((paths.position_at(0, 2) - hit).cwiseAbs().maxCoeff() <= 1e-12);
  CHECK((paths.direction_at(0, 2) - out).cwiseAbs().maxCoeff() <= 1e-12);
  // Content guard: a_0 moves the hit by more than 1e-4 mm from the sphere's.
  CHECK(std::abs(q.z() - c * r * r / (1.0 + std::sqrt(1.0 - c * c * r * r))) > 1e-4);
}

TEST_CASE("a Q-con that represents an even asphere hits the same points", "[trace][forbes]") {
  // ADR 0043, acceptance 3: the even asphere R = 64, k = -0.8, A4 ... A10 and its Q-con form
  // (back substitution of [F10a], Eqs. (5.6), (5.7) as in test_forbes_paraxial.cpp); both
  // intersections run Newton, the points agree within 1e-11 mm.
  const double rho = 12.7;
  const std::array<double, 4> a_even = {-2e-6, 3e-9, -1e-12, 2e-15};
  const std::array<std::array<double, 4>, 4> qm = {{{1.0, 0.0, 0.0, 0.0},
                                                    {-5.0, 6.0, 0.0, 0.0},
                                                    {15.0, -42.0, 28.0, 0.0},
                                                    {-35.0, 168.0, -252.0, 120.0}}};
  std::array<double, 4> s{};
  for (std::size_t m = 4; m-- > 0;) {
    double rest = a_even[m] * std::pow(rho, static_cast<double>(2 * m + 4));
    for (std::size_t n = m + 1; n < 4; ++n) rest -= s[n] * qm[n][m];
    s[m] = rest / qm[m][m];
  }
  std::vector<Param> even_c;
  std::vector<Param> qcon_c;
  for (std::size_t m = 0; m < 4; ++m) {
    even_c.emplace_back(a_even[m]);
    qcon_c.emplace_back(s[m]);
  }
  const System even = singlet(rtt::model::EvenAsphere{Param(64.0), Param(-0.8), even_c},
                              rtt::model::Plane{}, "CONST:1.5168");
  const System qcon = singlet(rtt::model::QCon{Param(1.0 / 64.0), Param(-0.8), Param(rho), qcon_c},
                              rtt::model::Plane{}, "CONST:1.5168");
  std::vector<Vec3> starts;
  std::vector<Vec3> dirs;
  for (int i = -4; i <= 4; ++i) {
    for (int j = -4; j <= 4; ++j) {
      if (i * i + j * j > 16) continue;
      starts.emplace_back(2.0 * i, 2.0 * j, -1.0);
      dirs.push_back(Vec3(0.01 * j, -0.02 * i, 1.0).normalized());
    }
  }
  const RayPaths a = traced(even, starts, dirs);
  const RayPaths b = traced(qcon, starts, dirs);
  double worst = 0.0;
  double rise = 0.0;
  for (std::size_t r = 0; r < starts.size(); ++r) {
    for (std::size_t slot = 0; slot < a.slots; ++slot) {
      worst =
          std::max(worst, (a.position_at(r, slot) - b.position_at(r, slot)).cwiseAbs().maxCoeff());
    }
    rise = std::max(rise, std::abs(a.position_at(r, 2).z() - 5.0));
  }
  CHECK(worst <= 1e-11);
  CHECK(rise > 0.1);  // content guard: the rays meet the curved front off the vertex plane
}
