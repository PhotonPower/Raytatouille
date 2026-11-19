#pragma once

/// @file stacked.hpp
/// A shape stack: a base shape plus additive terms (ADR 0041, point 1; ADR 0051, point 1). One
/// alternative of compile::CompiledShape for every surface with terms; a surface without terms
/// stays its base shape. The base is visited statically per call, the terms likewise.

#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/geom/anamorphic.hpp"
#include "rtt/geom/asphere.hpp"
#include "rtt/geom/axicon.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/geom/forbes.hpp"
#include "rtt/geom/odd_asphere.hpp"
#include "rtt/geom/plane.hpp"
#include "rtt/geom/shape.hpp"
#include "rtt/geom/terms.hpp"

namespace rtt::geom {

/// Every base shape (the alternatives of compile::CompiledShape except Stacked).
template <rtt::math::Real T>
using BaseShapeVariant = std::variant<Plane<T>,
                                      Conic<T>,
                                      EvenAsphere<T>,
                                      Biconic<T>,
                                      Toroid<T>,
                                      Cylinder<T>,
                                      OddAsphere<T>,
                                      Axicon<T>,
                                      QCon<T>,
                                      QBfs<T>>;

/// Every term type.
template <rtt::math::Real T>
using TermVariant = std::variant<ZernikeTerm<T>, XyPolynomialTerm<T>>;

/// z = base(x, y) + sum of the terms, in mm, local coordinates (ADR 0051, point 1).
/// - sag, grad: sums;
/// - base_conic: that of the base (start value of the intersection);
/// - max_radius: the smallest of base and terms (polynomial terms are unbounded);
/// - in_domain: the base's own domain (ADR 0042, point 1) or the disk of its max_radius, and every
///   term's domain;
/// - vertex_expansion: stack_expansion() of base and terms (ADR 0041, point 1).
/// The intersection is always intersect_newton() (ADR 0041, point 4).
template <rtt::math::Real T>
class Stacked final {
 public:
  Stacked(BaseShapeVariant<T> base, std::vector<TermVariant<T>> terms)
      : base_(std::move(base)), terms_(std::move(terms)) {}

  [[nodiscard]] T sag(T x, T y) const noexcept {
    T z = std::visit([&](const auto& b) noexcept { return b.sag(x, y); }, base_);
    for (const TermVariant<T>& t : terms_) {
      z += std::visit([&](const auto& s) noexcept { return s.sag(x, y); }, t);
    }
    return z;
  }

  [[nodiscard]] std::pair<T, T> grad(T x, T y) const noexcept {
    std::pair<T, T> g = std::visit([&](const auto& b) noexcept { return b.grad(x, y); }, base_);
    for (const TermVariant<T>& t : terms_) {
      const std::pair<T, T> d = std::visit([&](const auto& s) noexcept { return s.grad(x, y); }, t);
      g.first += d.first;
      g.second += d.second;
    }
    return g;
  }

  [[nodiscard]] std::pair<T, T> base_conic() const noexcept {
    return std::visit([](const auto& b) noexcept { return b.base_conic(); }, base_);
  }

  [[nodiscard]] std::optional<T> max_radius() const noexcept {
    std::optional<T> r = std::visit([](const auto& b) noexcept { return b.max_radius(); }, base_);
    for (const TermVariant<T>& t : terms_) {
      const std::optional<T> q =
          std::visit([](const auto& s) noexcept { return s.max_radius(); }, t);
      if (q && (!r || *q < *r)) r = q;
    }
    return r;
  }

  [[nodiscard]] bool in_domain(T x, T y) const noexcept {
    const bool base_in = std::visit(
        [&](const auto& b) noexcept {
          using B = std::decay_t<decltype(b)>;
          if constexpr (DomainShape<B, T>) {
            return b.in_domain(x, y);
          } else {
            const std::optional<T> r = b.max_radius();
            return !r || x * x + y * y <= *r * *r;
          }
        },
        base_);
    if (!base_in) return false;
    for (const TermVariant<T>& t : terms_) {
      const bool in = std::visit(
          [&](const auto& s) noexcept {
            if constexpr (requires { s.in_domain(x, y); }) {
              return s.in_domain(x, y);
            } else {
              return true;
            }
          },
          t);
      if (!in) return false;
    }
    return true;
  }

  [[nodiscard]] VertexExpansion<T> vertex_expansion() const {
    const VertexExpansion<T> base =
        std::visit([](const auto& b) { return b.vertex_expansion(); }, base_);
    std::vector<VertexExpansion<T>> terms;
    terms.reserve(terms_.size());
    for (const TermVariant<T>& t : terms_) {
      terms.push_back(std::visit([](const auto& s) { return s.vertex_expansion(); }, t));
    }
    return stack_expansion<T>(base, std::span<const VertexExpansion<T>>(terms));
  }

  [[nodiscard]] const BaseShapeVariant<T>& base() const noexcept { return base_; }
  [[nodiscard]] std::span<const TermVariant<T>> terms() const noexcept { return terms_; }

 private:
  BaseShapeVariant<T> base_;
  std::vector<TermVariant<T>> terms_;
};

}  // namespace rtt::geom
