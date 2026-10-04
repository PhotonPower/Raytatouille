#include "rtt/coating/tabulated.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace rtt::coating {
namespace {

constexpr double kAngleTolerance = 1e-12;  // rad

void check_grid(const std::vector<double>& grid, std::size_t min_size, const char* name) {
  if (grid.size() < min_size) {
    throw std::invalid_argument(std::string("coating table: too few ") + name);
  }
  for (std::size_t i = 0; i < grid.size(); ++i) {
    if (!std::isfinite(grid[i])) {
      throw std::invalid_argument(std::string("coating table: ") + name + " must be finite");
    }
    if (i > 0 && !(grid[i] > grid[i - 1])) {
      throw std::invalid_argument(std::string("coating table: ") + name +
                                  " must be strictly increasing");
    }
  }
}

bool finite(std::complex<double> v) {
  return std::isfinite(v.real()) && std::isfinite(v.imag());
}

/// Index i of the interval [grid[i], grid[i + 1]] containing x (clamped) and the weight of
/// grid[i + 1]; a grid with one point gives (0, 0).
std::pair<std::size_t, double> locate(const std::vector<double>& grid, double x) noexcept {
  if (grid.size() == 1 || !(x > grid.front())) return {0, 0.0};
  if (!(x < grid.back())) return {grid.size() - 2, 1.0};
  const auto upper = std::upper_bound(grid.begin(), grid.end(), x);
  const auto i = static_cast<std::size_t>(upper - grid.begin()) - 1;
  return {i, (x - grid[i]) / (grid[i + 1] - grid[i])};
}

Amplitudes<double> blend(const Amplitudes<double>& a, const Amplitudes<double>& b, double w) {
  const auto mix = [w](std::complex<double> x, std::complex<double> y) { return x + w * (y - x); };
  return {mix(a.rs, b.rs), mix(a.rp, b.rp), mix(a.ts, b.ts), mix(a.tp, b.tp)};
}

}  // namespace

TabulatedCoating::TabulatedCoating(std::vector<double> angles_rad,
                                   std::vector<double> wavelengths_um,
                                   std::vector<Amplitudes<double>> values)
    : angles_rad_(std::move(angles_rad)),
      wavelengths_um_(std::move(wavelengths_um)),
      values_(std::move(values)) {
  check_grid(angles_rad_, 2, "angles");
  if (std::abs(angles_rad_.front()) > kAngleTolerance ||
      std::abs(angles_rad_.back() - std::numbers::pi / 2.0) > kAngleTolerance) {
    throw std::invalid_argument("coating table: the angles must cover [0, pi/2] rad");
  }
  angles_rad_.front() = 0.0;
  angles_rad_.back() = std::numbers::pi / 2.0;
  check_grid(wavelengths_um_, 1, "wavelengths");
  if (!(wavelengths_um_.front() > 0.0)) {
    throw std::invalid_argument("coating table: the wavelengths must be > 0 um");
  }
  if (values_.size() != angles_rad_.size() * wavelengths_um_.size()) {
    throw std::invalid_argument("coating table: need one value per angle and wavelength");
  }
  for (const Amplitudes<double>& v : values_) {
    if (!finite(v.rs) || !finite(v.rp) || !finite(v.ts) || !finite(v.tp)) {
      throw std::invalid_argument("coating table: the amplitudes must be finite");
    }
  }
}

Amplitudes<double> TabulatedCoating::amplitudes(double angle_rad,
                                                double wavelength_um) const noexcept {
  const auto [ia, wa] = locate(angles_rad_, angle_rad);
  const auto [iw, ww] = locate(wavelengths_um_, wavelength_um);
  const std::size_t iw1 = wavelengths_um_.size() == 1 ? iw : iw + 1;
  const Amplitudes<double> low = blend(at(ia, iw), at(ia, iw1), ww);
  const Amplitudes<double> high = blend(at(ia + 1, iw), at(ia + 1, iw1), ww);
  return blend(low, high, wa);
}

}  // namespace rtt::coating
