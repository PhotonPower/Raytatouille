#include "rtt/coating/transfer_matrix.hpp"

#include <cmath>
#include <complex>
#include <span>
#include <stdexcept>
#include <string>

namespace rtt::coating {
namespace {

void check_index(std::complex<double> n, const std::string& what) {
  if (!std::isfinite(n.real()) || !std::isfinite(n.imag()) || n == 0.0) {
    throw std::invalid_argument("coating: " + what + " must be finite and not 0");
  }
  if (n.imag() < 0.0) {
    throw std::invalid_argument("coating: " + what +
                                " has kappa < 0 (gain is not supported, Byrnes App. C)");
  }
}

}  // namespace

void check_stack(std::complex<double> n_in,
                 std::span<const Layer<double>> layers,
                 std::complex<double> n_out,
                 double xi,
                 double wavelength_um) {
  if (!std::isfinite(wavelength_um) || !(wavelength_um > 0.0)) {
    throw std::invalid_argument("coating: the wavelength must be finite and > 0 um");
  }
  check_index(n_in, "the incident index");
  check_index(n_out, "the exit index");
  for (std::size_t j = 0; j < layers.size(); ++j) {
    const std::string name = "layer " + std::to_string(j);
    check_index(layers[j].index, "the index of " + name);
    if (!std::isfinite(layers[j].thickness_um) || layers[j].thickness_um < 0.0) {
      throw std::invalid_argument("coating: the thickness of " + name +
                                  " must be finite and >= 0 um");
    }
  }
  if (!std::isfinite(xi) || xi < 0.0) {
    throw std::invalid_argument("coating: xi = n sin(theta) must be finite and >= 0");
  }
  if (n_in.imag() == 0.0 && !(xi < n_in.real())) {
    throw std::invalid_argument(
        "coating: xi = n sin(theta) must be below the incident index (no grazing incidence)");
  }
}

}  // namespace rtt::coating
