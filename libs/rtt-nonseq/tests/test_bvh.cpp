// Slab test and BVH of the non-sequential scene (ADR 0057, point 3): the slab test against
// analytic entry parameters (PBRT 4ed, §6.1.2, §6.8.2), the tree against a brute-force search,
// and its structure (at most 4 primitives per leaf, boxes contain their children, determinism).
// Entry parameters of axis-parallel rays are differences and products of exact small numbers and
// compare exactly; the diagonal one goes through one division (relative error <= 1e-16).

#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <set>
#include <vector>

#include "rtt/nonseq/aabb.hpp"
#include "rtt/nonseq/bvh.hpp"

using rtt::math::Vec3;
using rtt::nonseq::Aabb;
using rtt::nonseq::Bvh;
using rtt::nonseq::BvhNode;
using rtt::nonseq::slab_test;
using rtt::nonseq::SlabRay;

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

Aabb box(const Vec3& lo, const Vec3& hi) {
  Aabb b;
  b.lo = lo;
  b.hi = hi;
  return b;
}

/// Deterministic boxes on a 5 x 5 x 4 lattice with sizes that vary per index.
std::vector<Aabb> lattice() {
  std::vector<Aabb> boxes;
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      for (int k = 0; k < 4; ++k) {
        const Vec3 c(10.0 * i, 7.0 * j, 13.0 * k);
        const double s = 1.0 + 0.25 * ((i + 2 * j + 3 * k) % 5);
        boxes.push_back(box(c - Vec3(s, 0.5 * s, 2.0), c + Vec3(s, 0.5 * s, 2.0)));
      }
    }
  }
  return boxes;
}

/// Leaf of every primitive and the containment of children, recursively.
void check_node(const Bvh& bvh, std::uint32_t index, std::vector<int>& seen) {
  const BvhNode& node = bvh.nodes()[index];
  if (node.count > 0) {
    CHECK(node.count <= Bvh::kMaxLeafSize);
    for (std::uint32_t i = 0; i < node.count; ++i) ++seen[bvh.order()[node.offset + i]];
    return;
  }
  for (const std::uint32_t child : {index + 1, node.offset}) {
    const Aabb& c = bvh.nodes()[child].bounds;
    CHECK((c.lo.array() >= node.bounds.lo.array()).all());
    CHECK((c.hi.array() <= node.bounds.hi.array()).all());
    check_node(bvh, child, seen);
  }
}

}  // namespace

TEST_CASE("slab test gives the analytic entry parameter", "[nonseq][bvh]") {
  const Aabb b = box(Vec3(1.0, -1.0, -1.0), Vec3(3.0, 1.0, 1.0));
  double t_near = -1.0;
  // Along +x from the origin: enters at x = 1.
  REQUIRE(slab_test(b, SlabRay::from(Vec3::Zero(), Vec3::UnitX()), kInf, t_near));
  CHECK(t_near == 1.0);
  // Along -x from x = 5: enters at x = 3, t = 2.
  REQUIRE(slab_test(b, SlabRay::from(Vec3(5.0, 0.0, 0.0), -Vec3::UnitX()), kInf, t_near));
  CHECK(t_near == 2.0);
  // Origin inside: t_near = 0.
  REQUIRE(slab_test(b, SlabRay::from(Vec3(2.0, 0.0, 0.0), Vec3::UnitY()), kInf, t_near));
  CHECK(t_near == 0.0);
  // Diagonal in the x-y plane: enters through x = 1 at t = sqrt(2) (y = 1, on the edge).
  const Vec3 d = Vec3(1.0, 1.0, 0.0).normalized();
  REQUIRE(slab_test(b, SlabRay::from(Vec3(0.0, 0.0, 0.0), d), kInf, t_near));
  CHECK(t_near == Catch::Approx(std::sqrt(2.0)).epsilon(1e-15));
  // Behind the ray, beyond t_max, and parallel outside a slab: no hit.
  CHECK_FALSE(slab_test(b, SlabRay::from(Vec3(5.0, 0.0, 0.0), Vec3::UnitX()), kInf, t_near));
  CHECK_FALSE(slab_test(b, SlabRay::from(Vec3::Zero(), Vec3::UnitX()), 0.5, t_near));
  CHECK_FALSE(slab_test(b, SlabRay::from(Vec3(0.0, 2.0, 0.0), Vec3::UnitX()), kInf, t_near));
  // Axis-parallel ray with its origin on a slab plane (0 * inf = NaN, PBRT 4ed, §6.1.2): the
  // NaN leaves the interval unchanged, the ray along the face y = 1 still meets the box.
  CHECK(slab_test(b, SlabRay::from(Vec3(0.0, 1.0, 0.0), Vec3::UnitX()), kInf, t_near));
  // The same on the x axis (the first slab): origin on the plane x = 1, direction +z.
  CHECK(slab_test(b, SlabRay::from(Vec3(1.0, 0.0, -5.0), Vec3::UnitZ()), kInf, t_near));
  CHECK(t_near == 4.0);
}

TEST_CASE("slab test keeps a far plane hit despite rounding (2 gamma_3)", "[nonseq][bvh]") {
  // A flat box (thickness 0 in z) hit exactly in its plane: t_near = t_far up to rounding;
  // the enlargement by 2 gamma_3 (PBRT 4ed, §6.8.2) keeps it.
  const Aabb b = box(Vec3(-1.0, -1.0, 0.3), Vec3(1.0, 1.0, 0.3));
  double t_near = 0.0;
  const Vec3 d = Vec3(0.1, 0.2, 1.0).normalized();
  CHECK(slab_test(b, SlabRay::from(Vec3(0.0, 0.0, -0.7), d), kInf, t_near));
  CHECK(rtt::nonseq::gamma_bound(3) == Catch::Approx(3.0 * std::ldexp(1.0, -53)).epsilon(1e-15));
}

TEST_CASE("BVH finds every box the ray meets, with at most 4 per leaf", "[nonseq][bvh]") {
  const std::vector<Aabb> boxes = lattice();
  const Bvh bvh(boxes);
  std::vector<int> seen(boxes.size(), 0);
  check_node(bvh, 0, seen);
  CHECK(std::all_of(seen.begin(), seen.end(), [](int n) { return n == 1; }));

  // Rays in several directions through the lattice: the visited set contains every box the
  // ray meets (brute force) when the bound t_max stays infinite.
  const std::vector<std::pair<Vec3, Vec3>> rays = {
      {Vec3(-5.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0)},
      {Vec3(-5.0, -5.0, -5.0), Vec3(25.0, 19.0, 31.0)},     // through the centre of box (2, 2, 2)
      {Vec3(45.0, 30.0, 45.0), Vec3(-35.0, -23.0, -32.0)},  // through the centre of box (1, 1, 1)
      {Vec3(20.0, 14.0, -10.0), Vec3(0.0, 0.0, 1.0)},
      {Vec3(20.0, -10.0, 26.0), Vec3(0.0, 1.0, 0.0)},
  };
  for (const auto& [o, dir] : rays) {
    const SlabRay ray = SlabRay::from(o, dir.normalized());
    std::set<std::uint32_t> visited;
    bvh.traverse(ray, kInf, [&](std::uint32_t prim, double) {
      visited.insert(prim);
      return kInf;
    });
    std::size_t met = 0;
    for (std::uint32_t i = 0; i < boxes.size(); ++i) {
      double t = 0.0;
      if (slab_test(boxes[i], ray, kInf, t)) {
        ++met;
        CHECK(visited.count(i) == 1);
      }
    }
    CHECK(met > 0);  // guard: the rays do pass through the lattice
  }
}

TEST_CASE("BVH is deterministic and skips nodes beyond t_max", "[nonseq][bvh]") {
  const std::vector<Aabb> boxes = lattice();
  const Bvh a(boxes);
  const Bvh b(boxes);
  REQUIRE(a.nodes().size() == b.nodes().size());
  CHECK(a.order() == b.order());
  for (std::size_t i = 0; i < a.nodes().size(); ++i) {
    CHECK(a.nodes()[i].offset == b.nodes()[i].offset);
    CHECK(a.nodes()[i].count == b.nodes()[i].count);
    CHECK(a.nodes()[i].axis == b.nodes()[i].axis);
  }
  // Along +x through the row j = 0, k = 0 from x = -5, outside every box: returning t_max = 0
  // at the first visit prunes every further node, so only the rest of the first leaf is seen.
  const SlabRay ray = SlabRay::from(Vec3(-5.0, 0.0, 0.0), Vec3::UnitX());
  std::size_t visits = 0;
  a.traverse(ray, kInf, [&](std::uint32_t, double) {
    ++visits;
    return 0.0;
  });
  CHECK(visits >= 1);
  CHECK(visits <= Bvh::kMaxLeafSize);
}

TEST_CASE("BVH keeps coincident centroids in one leaf and rejects infinite boxes",
          "[nonseq][bvh]") {
  std::vector<Aabb> same(6, box(Vec3(-1.0, -1.0, -1.0), Vec3(1.0, 1.0, 1.0)));
  const Bvh bvh(same);
  REQUIRE(bvh.nodes().size() == 1);
  CHECK(bvh.nodes()[0].count == 6);
  std::vector<Aabb> bad{box(Vec3::Zero(), Vec3::Constant(kInf))};
  CHECK_THROWS_AS(Bvh(bad), std::invalid_argument);
  CHECK(Bvh(std::vector<Aabb>{}).empty());
}
