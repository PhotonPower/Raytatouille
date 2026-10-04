#include "rtt/coating/ideal.hpp"

#include <cmath>
#include <complex>
#include <stdexcept>

namespace rtt::coating {

IdealCoating::IdealCoating(double reflectance) : reflectance_(reflectance) {
  if (!std::isfinite(reflectance) || reflectance < 0.0 || reflectance > 1.0) {
    throw std::invalid_argument("coating: an ideal reflectance must be in [0, 1]");
  }
}

Amplitudes<double> IdealCoating::amplitudes(std::complex<double> n_in,
                                            std::complex<double> n_out,
                                            double xi) const noexcept {
  const std::complex<double> q_in = normal_component(n_in, xi);
  const std::complex<double> q_out = normal_component(n_out, xi);
  if (!(q_out.real() > 0.0)) return interface_amplitudes(n_in, n_out, xi);  // TIR
  // Power normalisation of Byrnes Eq. (21)/(22): T = |t|^2 w_out / w_in.
  const std::complex<double> cos_in = q_in / n_in;
  const std::complex<double> cos_out = q_out / n_out;
  const double ws = q_in.real() / q_out.real();
  const double wp = (n_in * std::conj(cos_in)).real() / (n_out * std::conj(cos_out)).real();
  const double r = std::sqrt(reflectance_);
  const double t = std::sqrt(1.0 - reflectance_);
  return {-r, r, t * std::sqrt(ws), t * std::sqrt(wp)};
}

}  // namespace rtt::coating
