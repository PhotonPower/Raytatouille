#pragma once

/// @file polar_batch.hpp
/// Polarization quantities of a traced RayBatch (#62), used by the extension module
/// (raytatouille.polar) and by the test program rtt_py_reference, so that both compute the
/// same values.
///
/// Semantics of P and weight: ADR 0021. P is power-normalised: |P E|^2 is the power fraction of
/// the interfaces and polarizing elements for an incident state E (|E| = 1, transverse to the
/// initial direction k0), without the polarization-independent factors s (volume absorption,
/// absorber); weight = s ||P_T||_F^2 / 2 is the power for an unpolarized source. Directions are
/// global unit vectors; k is the ray direction after the trace.

#include <array>
#include <cstddef>
#include <span>
#include <vector>

#include "rtt/math/types.hpp"
#include "rtt/trace/ray_batch.hpp"

namespace rtt::py::polar_batch {

using CVec3 = Eigen::Vector3cd;

/// Tolerance of the input checks for unit length and transversality (|E| = 1, E . k0 = 0).
inline constexpr double kStateTolerance = 1e-12;

/// k . k0 at or above 1 - kSameDirection counts as "leaves in the incident direction"
/// (retardance() is defined only then).
inline constexpr double kSameDirection = 1e-12;

/// Initial direction of every ray: k0 = Re(P^T k). From P = P_T + k k0^T with k^T P_T = 0
/// (ADR 0021; Lam, Eqs. (3.2), (3.4)) follows P^T k = k0 for |k| = 1.
[[nodiscard]] std::vector<math::Vec3> initial_directions(const trace::RayBatch& rays);

/// The state `e` projected perpendicular to k0 of every ray and normalised:
/// (e - k (k . e)) / |...| with k = k0 / |k0|, a valid input of transmission() and stokes().
/// The projection is a convention (the far field of a dipole along e); it is singular for e
/// parallel to k0.
/// @throws std::invalid_argument if e is not finite or (nearly) parallel to a k0
///         (|projection| < kStateTolerance |e|)
[[nodiscard]] std::vector<CVec3> transverse_polarization(const trace::RayBatch& rays,
                                                         const CVec3& e);

/// Power of every ray for an unpolarized source: the weight column (copy). Dimensionless
/// (source normalised to 1).
[[nodiscard]] std::vector<double> transmission(const trace::RayBatch& rays);

/// Power of every ray for the incident state E (one for all rays or one per ray):
/// weight |P E|^2 / (||P_T||_F^2 / 2) with P_T = (I - k k^T) P, and 0 where P_T = 0.
/// Dimensionless (source normalised to 1).
/// @throws std::invalid_argument if `states` has neither 1 nor rays.size() entries (an empty
///         list is an error, not "unpolarized"), or a state is not finite, not of unit length or
///         not transverse to k0 (tolerance kStateTolerance)
[[nodiscard]] std::vector<double> transmission(const trace::RayBatch& rays,
                                               std::span<const CVec3> states);

/// Diattenuation of every ray (rtt::polar::diattenuation with k0 and k). With the
/// power-normalised P of the tracer the value is the power diattenuation, maximum^2 and
/// minimum^2 are power fractions without s (ADR 0021).
struct Diattenuations {
  std::vector<double> value;    ///< D in [0, 1]
  std::vector<double> maximum;  ///< largest transverse singular value of P
  std::vector<double> minimum;  ///< smallest transverse singular value of P
  std::vector<CVec3> axis;      ///< incident state with the largest transmission
};
[[nodiscard]] Diattenuations diattenuation(const trace::RayBatch& rays);

/// Retardance of every ray (rtt::polar::retardance of P about k) for rays that leave in their
/// incident direction (k . k0 >= 1 - kSameDirection); NaN value and fast axis otherwise. This is
/// the retardance of P including the geometric transformation: it equals the physical
/// retardance only if the transverse part of the path's Q is the identity (in-plane plates,
/// ideal thin elements); skew rays or out-of-plane folds that return to k0 add a geometric
/// rotation (e.g. delta = pi for a periscope with 90 deg image rotation). The physical
/// retardance of a whole path needs the accumulated Q, which the tracer does not keep (M6).
struct Retardances {
  std::vector<double> value;     ///< delta in [0, pi], rad
  std::vector<CVec3> fast_axis;  ///< eigenpolarization with the smaller phase
};
[[nodiscard]] Retardances retardance(const trace::RayBatch& rays);

/// Stokes parameters (s0, s1, s2, s3) of P E for every ray, in units of |E|^2 = 1 (power
/// fraction without s), in the basis e1 = `axis` projected perpendicular to k, e2 = k x e1
/// (rtt::polar::stokes; docs/architecture.md, Händigkeit und Stokes).
/// @throws std::invalid_argument as transmission() for `states`, or if `axis` is not finite or
///         (nearly) parallel to the direction k of a ray
[[nodiscard]] std::vector<std::array<double, 4>> stokes(const trace::RayBatch& rays,
                                                        std::span<const CVec3> states,
                                                        const math::Vec3& axis);

}  // namespace rtt::py::polar_batch
