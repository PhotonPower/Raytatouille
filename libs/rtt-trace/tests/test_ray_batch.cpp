#include <catch2/catch_test_macros.hpp>
#include <cstddef>

#include "rtt/trace/ray_batch.hpp"

using rtt::math::CMat3;
using rtt::math::Complex;
using rtt::trace::kNoSurface;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;

namespace {

/// Every field array must have exactly n entries.
void require_size(const RayBatch& b, std::size_t n) {
  REQUIRE(b.size() == n);
  REQUIRE(b.pos_x().size() == n);
  REQUIRE(b.pos_y().size() == n);
  REQUIRE(b.pos_z().size() == n);
  REQUIRE(b.dir_x().size() == n);
  REQUIRE(b.dir_y().size() == n);
  REQUIRE(b.dir_z().size() == n);
  REQUIRE(b.wl().size() == n);
  REQUIRE(b.opl().size() == n);
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) REQUIRE(b.prt(r, c).size() == n);
  }
  REQUIRE(b.weight().size() == n);
  REQUIRE(b.status().size() == n);
  REQUIRE(b.field().size() == n);
  REQUIRE(b.pupil_x().size() == n);
  REQUIRE(b.pupil_y().size() == n);
  REQUIRE(b.last_surface().size() == n);
}

/// Ray i must be in the documented default state.
void require_default_ray(const RayBatch& b, std::size_t i) {
  REQUIRE(b.pos_x()[i] == 0.0);
  REQUIRE(b.pos_y()[i] == 0.0);
  REQUIRE(b.pos_z()[i] == 0.0);
  REQUIRE(b.dir_x()[i] == 0.0);
  REQUIRE(b.dir_y()[i] == 0.0);
  REQUIRE(b.dir_z()[i] == 1.0);
  REQUIRE(b.wl()[i] == 0);
  REQUIRE(b.opl()[i] == 0.0);
  REQUIRE(b.prt_matrix(i) == CMat3::Identity());
  REQUIRE(b.weight()[i] == 1.0);
  REQUIRE(b.status()[i] == RayStatus::Alive);
  REQUIRE(b.field()[i] == 0);
  REQUIRE(b.pupil_x()[i] == 0.0);
  REQUIRE(b.pupil_y()[i] == 0.0);
  REQUIRE(b.last_surface()[i] == kNoSurface);
}

}  // namespace

TEST_CASE("an empty batch has no rays", "[ray_batch]") {
  const RayBatch b;
  require_size(b, 0);
}

TEST_CASE("new rays start in the documented default state", "[ray_batch]") {
  const RayBatch b(5);
  require_size(b, 5);
  for (std::size_t i = 0; i < b.size(); ++i) require_default_ray(b, i);
}

TEST_CASE("resize keeps existing rays and initialises new ones", "[ray_batch]") {
  RayBatch b(2);
  b.pos_x()[1] = 3.5;
  b.opl()[1] = 7.0;
  b.status()[0] = RayStatus::Vignetted;
  b.resize(4);
  require_size(b, 4);
  REQUIRE(b.pos_x()[1] == 3.5);
  REQUIRE(b.opl()[1] == 7.0);
  REQUIRE(b.status()[0] == RayStatus::Vignetted);
  require_default_ray(b, 2);
  require_default_ray(b, 3);
  b.resize(1);
  require_size(b, 1);
  REQUIRE(b.status()[0] == RayStatus::Vignetted);
}

TEST_CASE("prt(row, col) addresses the matrix element P(row, col)", "[ray_batch]") {
  RayBatch b(3);
  CMat3 p;
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) p(r, c) = Complex(10.0 * r + c, -1.0 * (r + 1));
  }
  b.set_prt_matrix(1, p);
  REQUIRE(b.prt_matrix(1) == p);
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) REQUIRE(b.prt(r, c)[1] == p(r, c));
  }
  // Neighbouring rays are untouched.
  REQUIRE(b.prt_matrix(0) == CMat3::Identity());
  REQUIRE(b.prt_matrix(2) == CMat3::Identity());

  b.prt(2, 0)[2] = Complex(0.0, 4.0);
  REQUIRE(b.prt_matrix(2)(2, 0) == Complex(0.0, 4.0));
}
