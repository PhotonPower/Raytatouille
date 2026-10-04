#include "rtt/polar/fresnel.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace rtt::polar {

namespace {

void check_index(std::complex<double> n, const char* name) {
  if (!std::isfinite(n.real()) || !std::isfinite(n.imag())) {
    throw std::invalid_argument(std::string("fresnel: ") + name + " is not finite");
  }
  if (!(n.real() > 0.0)) {
    throw std::invalid_argument(std::string("fresnel: Re(") + name + ") must be > 0");
  }
  if (n.imag() < 0.0) {
    throw std::invalid_argument(std::string("fresnel: Im(") + name +
                                ") must be >= 0 (kappa >= 0, no gain)");
  }
}

void check_input(std::complex<double> n_i, std::complex<double> n_t, double xi) {
  check_index(n_i, "n_i");
  check_index(n_t, "n_t");
  if (!std::isfinite(xi)) throw std::invalid_argument("fresnel: xi is not finite");
  if (xi < 0.0 || xi > n_i.real()) {
    throw std::invalid_argument("fresnel: xi = Re(n_i) sin(theta_i) must lie in [0, Re(n_i)]");
  }
  const std::complex<double> q_i = normal_component(n_i, xi);
  const std::complex<double> q_t = normal_component(n_t, xi);
  if (q_i + q_t == 0.0 || n_t * n_t * q_i + n_i * n_i * q_t == 0.0) {
    throw std::invalid_argument(
        "fresnel: a Fresnel denominator vanishes (grazing incidence between equal media)");
  }
}

}  // namespace

FresnelAmplitudes<double> fresnel_checked(std::complex<double> n_i,
                                          std::complex<double> n_t,
                                          double xi) {
  check_input(n_i, n_t, xi);
  return fresnel(n_i, n_t, xi);
}

FresnelPower<double> fresnel_power_checked(std::complex<double> n_i,
                                           std::complex<double> n_t,
                                           double xi) {
  check_input(n_i, n_t, xi);
  if (!(normal_component(n_i, xi).real() > 0.0)) {
    throw std::invalid_argument("fresnel: grazing incidence carries no power to the interface");
  }
  return fresnel_power(n_i, n_t, xi);
}

}  // namespace rtt::polar
