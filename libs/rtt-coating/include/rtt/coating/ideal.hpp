#pragma once

/// @file ideal.hpp
/// Ideal coatings without retardance: fixed reflectance, everything else transmitted without
/// loss. Amplitudes in the basis [s, p, k] of rtt-polar, Convention A (transfer_matrix.hpp):
/// a reflection without phase difference between s and p has r_s = -sqrt(R), r_p = +sqrt(R)
/// (the limit of an ideal conductor, n -> infinity, of Byrnes Eq. (6)).

#include <complex>

#include "rtt/coating/transfer_matrix.hpp"

namespace rtt::coating {

/// Ideal coating with reflectance R in [0, 1] for s and p and transmittance 1 - R.
/// R = 0: ideal anti-reflection; R = 1: ideal mirror (r_s = -1, r_p = +1, no transmission);
/// in between: ideal beam splitter. The transmitted amplitudes are real and positive, chosen so
/// that T = 1 - R by Byrnes Eq. (21)/(22) for the given media. Under total internal reflection
/// (no propagating wave in the exit medium) nothing can be transmitted: the coating then
/// reflects totally with the amplitudes of the bare interface (Eq. (6)).
class IdealCoating {
 public:
  /// @param reflectance R for s and p, dimensionless
  /// @throws std::invalid_argument if R is not finite or outside [0, 1]
  explicit IdealCoating(double reflectance);

  [[nodiscard]] static IdealCoating anti_reflection() { return IdealCoating(0.0); }
  [[nodiscard]] static IdealCoating mirror() { return IdealCoating(1.0); }

  [[nodiscard]] double reflectance() const noexcept { return reflectance_; }

  /// Amplitudes between the incident medium n_in and the exit medium n_out at the tangential
  /// invariant xi = n sin(theta) (transfer_matrix.hpp).
  /// @pre the incident medium does not absorb; not at grazing incidence
  [[nodiscard]] Amplitudes<double> amplitudes(std::complex<double> n_in,
                                              std::complex<double> n_out,
                                              double xi) const noexcept;

 private:
  double reflectance_;
};

}  // namespace rtt::coating
