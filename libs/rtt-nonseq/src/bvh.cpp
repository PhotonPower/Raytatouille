#include "rtt/nonseq/bvh.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

#include "rtt/nonseq/aabb.hpp"

namespace rtt::nonseq {

SlabRay SlabRay::from(const math::Vec3& origin, const math::Vec3& direction) noexcept {
  SlabRay r;
  r.origin = origin;
  for (int a = 0; a < 3; ++a) {
    // 1 / 0 = +-inf as in PBRT 4ed, §7.3.5 (IEEE 754); the slab test handles it.
    r.inv_dir[a] = 1.0 / direction[a];
    r.dir_is_neg[static_cast<std::size_t>(a)] = r.inv_dir[a] < 0.0;
  }
  return r;
}

bool slab_test(const Aabb& box, const SlabRay& ray, double t_max, double& t_near) noexcept {
  // PBRT 4ed, §6.8.2, Bounds3::IntersectP with precomputed reciprocal direction (§7.3.5): the
  // near and far planes of each slab follow from the sign of the direction; the far parameter
  // grows by 2 gamma_3 so that the test stays conservative.
  const auto lo_hi = [&](int axis, bool neg) { return neg ? box.hi[axis] : box.lo[axis]; };
  const double grow = 1.0 + 2.0 * gamma_bound(3);
  // Start with the whole line and narrow it slab by slab, every axis alike (§6.1.2), so that a
  // NaN from 0 * inf on any axis leaves the interval unchanged.
  double t0 = -std::numeric_limits<double>::infinity();
  double t1 = std::numeric_limits<double>::infinity();
  for (int a = 0; a < 3; ++a) {
    const bool neg = ray.dir_is_neg[static_cast<std::size_t>(a)];
    const double s0 = (lo_hi(a, neg) - ray.origin[a]) * ray.inv_dir[a];
    double s1 = (lo_hi(a, !neg) - ray.origin[a]) * ray.inv_dir[a];
    s1 *= grow;
    if (t0 > s1 || s0 > t1) return false;
    if (s0 > t0) t0 = s0;
    if (s1 < t1) t1 = s1;
  }
  if (!(t0 < t_max && t1 > 0.0)) return false;
  t_near = t0 > 0.0 ? t0 : 0.0;
  return true;
}

namespace {

struct Builder {
  std::span<const Aabb> boxes;
  std::vector<math::Vec3> centroids;
  std::vector<std::uint32_t> order;
  std::vector<BvhNode> nodes;

  /// Builds the subtree over order[begin, end) and returns its node index.
  std::uint32_t build(std::size_t begin, std::size_t end, std::size_t depth) {
    if (depth >= Bvh::kStackSize) {
      throw std::invalid_argument("nonseq::Bvh: tree deeper than the traversal stack");
    }
    const auto index = static_cast<std::uint32_t>(nodes.size());
    nodes.emplace_back();
    Aabb bounds;
    Aabb centroid_bounds;
    for (std::size_t i = begin; i < end; ++i) {
      bounds.grow(boxes[order[i]]);
      centroid_bounds.grow(centroids[order[i]]);
    }
    nodes[index].bounds = bounds;
    const std::size_t n = end - begin;
    const auto make_leaf = [&] {
      if (n > 0xffff) throw std::invalid_argument("nonseq::Bvh: too many primitives in a leaf");
      nodes[index].offset = static_cast<std::uint32_t>(begin);
      nodes[index].count = static_cast<std::uint16_t>(n);
      return index;
    };
    if (n == 1 || !(bounds.surface_area() > 0.0)) return make_leaf();
    // Axis of the largest centroid extent (PBRT 4ed, §7.3.2, MaxDimension); the smaller axis wins
    // on equal extent.
    const math::Vec3 extent = centroid_bounds.hi - centroid_bounds.lo;
    int axis = 0;
    if (extent.y() > extent[axis]) axis = 1;
    if (extent.z() > extent[axis]) axis = 2;
    const double c_lo = centroid_bounds.lo[axis];
    const double c_hi = centroid_bounds.hi[axis];
    if (!(c_hi > c_lo)) {
      return make_leaf();  // all centroids coincide: no split helps (§7.3.2)
    }
    const auto first = order.begin() + static_cast<std::ptrdiff_t>(begin);
    const auto last = order.begin() + static_cast<std::ptrdiff_t>(end);
    std::ptrdiff_t mid = 0;
    if (n == 2) {
      // Two primitives: one each side (§7.3.2 splits at the median for n <= 2).
      if (centroids[*(first + 1)][axis] < centroids[*first][axis]) std::iter_swap(first, first + 1);
      mid = 1;
    } else {
      const auto bucket_of = [&](std::uint32_t prim) {
        const double offset = (centroids[prim][axis] - c_lo) / (c_hi - c_lo);
        const int b = static_cast<int>(Bvh::kBuckets * offset);
        return std::clamp(b, 0, Bvh::kBuckets - 1);
      };
      std::array<std::size_t, Bvh::kBuckets> counts{};
      std::array<Aabb, Bvh::kBuckets> bucket_bounds{};
      for (auto it = first; it != last; ++it) {
        const auto b = static_cast<std::size_t>(bucket_of(*it));
        ++counts[b];
        bucket_bounds[b].grow(boxes[*it]);
      }
      // Cost n_A S_A + n_B S_B of the split after bucket i, swept from both ends (§7.3.2).
      constexpr std::size_t kSplits = Bvh::kBuckets - 1;
      std::array<double, kSplits> costs{};
      std::size_t count_below = 0;
      Aabb below;
      for (std::size_t i = 0; i < kSplits; ++i) {
        below.grow(bucket_bounds[i]);
        count_below += counts[i];
        costs[i] += static_cast<double>(count_below) * below.surface_area();
      }
      std::size_t count_above = 0;
      Aabb above;
      for (std::size_t i = kSplits; i >= 1; --i) {
        above.grow(bucket_bounds[i]);
        count_above += counts[i];
        costs[i - 1] += static_cast<double>(count_above) * above.surface_area();
      }
      std::size_t best = 0;
      for (std::size_t i = 1; i < kSplits; ++i) {
        if (costs[i] < costs[best]) best = i;  // strict: the smaller bucket wins on equal cost
      }
      // Eq. (7.1) with t_trav = 1/2, t_isect = 1, p = S_child / S.
      const double split_cost = 0.5 + costs[best] / bounds.surface_area();
      const auto leaf_cost = static_cast<double>(n);
      if (n <= Bvh::kMaxLeafSize && !(split_cost < leaf_cost)) return make_leaf();
      const auto split = std::stable_partition(first, last, [&](std::uint32_t prim) {
        return static_cast<std::size_t>(bucket_of(prim)) <= best;
      });
      mid = split - first;
    }
    nodes[index].axis = static_cast<std::uint8_t>(axis);
    build(begin, begin + static_cast<std::size_t>(mid), depth + 1);
    const std::uint32_t second = build(begin + static_cast<std::size_t>(mid), end, depth + 1);
    nodes[index].offset = second;
    nodes[index].count = 0;
    return index;
  }
};

}  // namespace

Bvh::Bvh(std::span<const Aabb> boxes) {
  if (boxes.empty()) return;
  Builder b{boxes, {}, {}, {}};
  b.centroids.reserve(boxes.size());
  b.order.reserve(boxes.size());
  for (std::size_t i = 0; i < boxes.size(); ++i) {
    if (!boxes[i].finite()) throw std::invalid_argument("nonseq::Bvh: box not finite");
    b.centroids.push_back(boxes[i].centroid());
    b.order.push_back(static_cast<std::uint32_t>(i));
  }
  b.nodes.reserve(2 * boxes.size());
  b.build(0, boxes.size(), 0);
  nodes_ = std::move(b.nodes);
  order_ = std::move(b.order);
}

}  // namespace rtt::nonseq
