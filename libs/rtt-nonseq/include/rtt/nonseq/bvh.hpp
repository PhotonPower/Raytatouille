#pragma once

/// @file bvh.hpp
/// Bounding volume hierarchy over the faces of the non-sequential scene (ADR 0057, point 3):
/// surface area heuristic with 12 buckets (PBRT 4ed, §7.3.2, Eq. (7.1)), at most 4 faces per
/// leaf, traversal with a fixed stack (§7.3.5). Immutable after construction; safe to read from
/// many threads.

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "rtt/nonseq/aabb.hpp"

namespace rtt::nonseq {

/// One node of the flattened tree (PBRT 4ed, §7.3.4): an interior node has its first child
/// directly after it and the second at `offset`; a leaf holds the primitives
/// order()[offset, offset + count).
struct BvhNode {
  Aabb bounds;
  std::uint32_t offset = 0;  ///< leaf: first index into Bvh::order(); interior: second child
  std::uint16_t count = 0;   ///< number of primitives, 0 for an interior node
  std::uint8_t axis = 0;     ///< split axis of an interior node (0 = x, 1 = y, 2 = z)
};

/// BVH over boxes given by index.
class Bvh {
 public:
  /// Most primitives in a leaf unless all centroids coincide (ADR 0057, point 3).
  static constexpr std::size_t kMaxLeafSize = 4;
  /// Number of SAH buckets (PBRT 4ed, §7.3.2).
  static constexpr int kBuckets = 12;
  /// Depth bound of the traversal stack (PBRT 4ed, §7.3.5: 64 entries).
  static constexpr std::size_t kStackSize = 64;

  Bvh() = default;

  /// Builds the tree over `boxes` (every box finite). Split rule (PBRT 4ed, §7.3.2): the axis of
  /// the largest extent of the centroids (the smaller axis on equal extent); all centroids equal
  /// on it -> one leaf; two primitives -> one each side; otherwise SAH over 12 buckets with
  /// cost 1/2 + (n_A S_A + n_B S_B) / S (t_trav = 1/2, t_isect = 1, Eq. (7.1)), the smallest
  /// bucket on equal cost, and a leaf if it holds at most 4 primitives and its cost n is not
  /// larger. Partitions are stable, so the tree depends only on the boxes and their order.
  /// @throws std::invalid_argument if a box is not finite, the tree would be deeper than the
  ///         traversal stack, or more than 65535 boxes share one centroid
  explicit Bvh(std::span<const Aabb> boxes);

  [[nodiscard]] const std::vector<BvhNode>& nodes() const noexcept { return nodes_; }
  /// Primitive indices in leaf order.
  [[nodiscard]] const std::vector<std::uint32_t>& order() const noexcept { return order_; }
  [[nodiscard]] bool empty() const noexcept { return nodes_.empty(); }

  /// Visits the primitives whose leaf boxes the ray meets (PBRT 4ed, §7.3.5), near child first
  /// by the sign of the direction on the split axis. `visit(index, t_near)` returns the new
  /// upper bound t_max (the t of the best hit so far, or the old value); nodes entered beyond
  /// it are skipped. t_near is the entry parameter of the leaf box.
  template <class Visit>
  void traverse(const SlabRay& ray, double t_max, Visit&& visit) const {
    if (nodes_.empty()) return;
    std::array<std::uint32_t, kStackSize> stack{};
    std::size_t top = 0;
    std::uint32_t current = 0;
    for (;;) {
      const BvhNode& node = nodes_[current];
      double t_near = 0.0;
      if (slab_test(node.bounds, ray, t_max, t_near)) {
        if (node.count > 0) {
          for (std::uint32_t i = 0; i < node.count; ++i) {
            t_max = visit(order_[node.offset + i], t_near);
          }
          if (top == 0) break;
          current = stack[--top];
        } else if (ray.dir_is_neg[node.axis]) {
          stack[top++] = current + 1;
          current = node.offset;
        } else {
          stack[top++] = node.offset;
          current = current + 1;
        }
      } else {
        if (top == 0) break;
        current = stack[--top];
      }
    }
  }

 private:
  std::vector<BvhNode> nodes_;
  std::vector<std::uint32_t> order_;
};

}  // namespace rtt::nonseq
