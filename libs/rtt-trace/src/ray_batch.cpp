#include "rtt/trace/ray_batch.hpp"

#include <cassert>

namespace rtt::trace {

RayBatch::RayBatch(std::size_t size) {
  resize(size);
}

void RayBatch::resize(std::size_t size) {
  // Default state documented in ray_batch.hpp.
  pos_x_.resize(size, 0.0);
  pos_y_.resize(size, 0.0);
  pos_z_.resize(size, 0.0);
  dir_x_.resize(size, 0.0);
  dir_y_.resize(size, 0.0);
  dir_z_.resize(size, 1.0);
  wl_.resize(size, 0);
  opl_.resize(size, 0.0);
  for (std::size_t k = 0; k < prt_.size(); ++k) {
    const bool diagonal = k % 4 == 0;  // row-major indices 0, 4, 8
    prt_[k].resize(size, diagonal ? math::Complex(1.0) : math::Complex(0.0));
  }
  weight_.resize(size, 1.0);
  status_.resize(size, RayStatus::Alive);
  field_.resize(size, 0);
  pupil_x_.resize(size, 0.0);
  pupil_y_.resize(size, 0.0);
  last_surface_.resize(size, kNoSurface);
}

namespace {

/// Row-major index of P(row, col) in RayBatch::prt_.
std::size_t prt_index(int row, int col) noexcept {
  assert(row >= 0 && row < 3 && col >= 0 && col < 3);
  return 3 * static_cast<std::size_t>(row) + static_cast<std::size_t>(col);
}

}  // namespace

std::span<math::Complex> RayBatch::prt(int row, int col) noexcept {
  return prt_[prt_index(row, col)];
}

std::span<const math::Complex> RayBatch::prt(int row, int col) const noexcept {
  return prt_[prt_index(row, col)];
}

math::CMat3 RayBatch::prt_matrix(std::size_t i) const {
  math::CMat3 p;
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) p(r, c) = prt(r, c)[i];
  }
  return p;
}

void RayBatch::set_prt_matrix(std::size_t i, const math::CMat3& p) {
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) prt(r, c)[i] = p(r, c);
  }
}

}  // namespace rtt::trace
