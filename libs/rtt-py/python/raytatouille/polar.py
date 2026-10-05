"""Polarization quantities (rtt-polar) of traced rays and of single PRT matrices.

Semantics of a traced RayBatch (ADR 0021):

- ``RayBatch.prt`` holds the accumulated 3x3 polarization ray-tracing matrix P in global
  coordinates. P is POWER-NORMALISED: for an incident state E (|E| = 1, transverse to the
  initial direction k0) ``|P E|^2`` is the power fraction of the interfaces and polarizing
  elements, without the polarization-independent factors s (volume absorption, absorber). The
  phases are those of the field PRT matrix (Chipman, Lam); the field amplitude differs by the
  factors sqrt(c) of the interfaces passed. P k0 = k.
- ``RayBatch.weight`` is the power for an unpolarized source (source = 1):
  weight = s ||P_T||_F^2 / 2 with P_T = (I - k k^T) P. For a polarized state the power is
  ``transmission(rays, E)`` = weight |P E|^2 / (||P_T||_F^2 / 2).

The batch functions take a traced RayBatch and return NumPy arrays with one row per ray.
Polarization states are complex (3,) (one for all rays) or (N, 3) arrays, of unit length and
transverse to k0 (tolerance 1e-12); ``transverse_polarization`` projects a vector for every ray.
The physical retardance of a whole path (Q^-1 P with the accumulated Q) is not available yet
(M6); ``retardance`` is defined for rays that leave in their incident direction.

The single-matrix functions take NumPy arrays: unit vectors (3,) and matrices (3, 3).
Conventions as in C++ (docs/architecture.md, Polarisation und Fresnel): fields
~ exp(i(k.r - omega t)), Convention A (r_p = -r_s at normal incidence), basis s, p = k x s.
Invalid input raises ValueError.
"""

from __future__ import annotations

from typing import NamedTuple, overload

import numpy as np
import numpy.typing as npt

from . import _core
from ._core import RayBatch

__all__ = [
    "Diattenuation",
    "Diattenuations",
    "Retardance",
    "Retardances",
    "diattenuation",
    "geometric_transform",
    "initial_directions",
    "physical_retardance",
    "prt_matrix",
    "retardance",
    "stokes",
    "transmission",
    "transverse_polarization",
]

FloatArray = npt.NDArray[np.float64]
ComplexArray = npt.NDArray[np.complex128]


class Diattenuation(NamedTuple):
    """Diattenuation of one PRT matrix (Lam, Eq. (4.3)): D = (L1^2 - L2^2) / (L1^2 + L2^2) from
    the transverse singular values L1 >= L2. For the power-normalised P of the tracer it is the
    power diattenuation and L1^2, L2^2 are power fractions."""

    value: float
    maximum: float
    minimum: float
    axis: ComplexArray  #: incident state (3,) with the largest transmission


class Diattenuations(NamedTuple):
    """Diattenuation of every ray: arrays (N,) and axis (N, 3) complex."""

    value: FloatArray
    maximum: FloatArray
    minimum: FloatArray
    axis: ComplexArray


class Retardance(NamedTuple):
    """Retardance delta in [0, pi] rad and the fast axis (3,) complex (Lam, Eq. (4.4))."""

    value: float
    fast_axis: ComplexArray


class Retardances(NamedTuple):
    """Retardance of every ray: value (N,) in rad, fast_axis (N, 3); NaN where undefined."""

    value: FloatArray
    fast_axis: ComplexArray


def _real(a: npt.ArrayLike) -> FloatArray:
    return np.ascontiguousarray(a, dtype=np.float64)


def _complex(a: npt.ArrayLike) -> ComplexArray:
    return np.ascontiguousarray(a, dtype=np.complex128)


def initial_directions(rays: RayBatch) -> FloatArray:
    """Initial unit direction k0 of every ray, (N, 3): k0 = Re(P^T k) with the direction k after
    the trace (from P = P_T + k k0^T, k^T P_T = 0; ADR 0021)."""
    result: FloatArray = _core.polar_initial_directions(rays)
    return result


def transverse_polarization(rays: RayBatch, polarization: npt.ArrayLike) -> ComplexArray:
    """``polarization`` (3,) projected perpendicular to k0 of every ray and normalised, (N, 3)
    complex: a valid input of transmission() and stokes(). Raises ValueError if it is not finite
    or parallel to a k0."""
    result: ComplexArray = _core.polar_transverse_polarization(rays, _complex(polarization))
    return result


def transmission(rays: RayBatch, polarization: npt.ArrayLike | None = None) -> FloatArray:
    """Power of every ray, (N,), dimensionless (source = 1). Without ``polarization`` this is
    ``rays.weight`` (unpolarized source, copied); otherwise weight |P E|^2 / (||P_T||^2 / 2) for
    the state E, (3,) or (N, 3), and 0 where P_T = 0."""
    states = None if polarization is None else _complex(polarization)
    result: FloatArray = _core.polar_transmission(rays, states)
    return result


@overload
def diattenuation(rays: RayBatch) -> Diattenuations: ...


@overload
def diattenuation(rays: npt.ArrayLike, k_in: npt.ArrayLike,
                  k_out: npt.ArrayLike) -> Diattenuation: ...


def diattenuation(rays: RayBatch | npt.ArrayLike, k_in: npt.ArrayLike | None = None,
                  k_out: npt.ArrayLike | None = None) -> Diattenuations | Diattenuation:
    """Diattenuation of every ray of a traced RayBatch (with k0 and k), or of one PRT matrix
    ``p`` (3, 3) with P k_in = k_out: ``diattenuation(p, k_in, k_out)``."""
    if isinstance(rays, RayBatch):
        if k_in is not None or k_out is not None:
            raise TypeError("diattenuation(rays) takes no directions")
        value, maximum, minimum, axis = _core.polar_diattenuation(rays)
        return Diattenuations(value, maximum, minimum, axis)
    if k_in is None or k_out is None:
        raise TypeError("diattenuation(p, k_in, k_out) needs both directions")
    value1, maximum1, minimum1, axis1 = _core.polar_diattenuation_matrix(
        _complex(rays), _real(k_in), _real(k_out))
    return Diattenuation(value1, maximum1, minimum1, axis1)


@overload
def retardance(rays: RayBatch) -> Retardances: ...


@overload
def retardance(rays: npt.ArrayLike, k: npt.ArrayLike) -> Retardance: ...


def retardance(rays: RayBatch | npt.ArrayLike,
               k: npt.ArrayLike | None = None) -> Retardances | Retardance:
    """Retardance of every ray that leaves in its incident direction (k . k0 >= 1 - 1e-12),
    NaN otherwise; or of one matrix ``m`` (3, 3) with m k = k: ``retardance(m, k)``. Ideal
    polarizers make the decomposition ambiguous (rtt/polar/prt_analysis.hpp)."""
    if isinstance(rays, RayBatch):
        if k is not None:
            raise TypeError("retardance(rays) takes no direction")
        value, fast_axis = _core.polar_retardance(rays)
        return Retardances(value, fast_axis)
    if k is None:
        raise TypeError("retardance(m, k) needs the direction k")
    value1, fast_axis1 = _core.polar_retardance_matrix(_complex(rays), _real(k))
    return Retardance(value1, fast_axis1)


def physical_retardance(p: npt.ArrayLike, q: npt.ArrayLike, k_in: npt.ArrayLike) -> Retardance:
    """Physical retardance of one intercept or path: retardance of Q^-1 P (Lam, Sec. 4.5.1) for
    the PRT matrix p (3, 3) and the geometric transformation q (3, 3), both mapping k_in onto
    the same k_out."""
    value, fast_axis = _core.polar_physical_retardance(_complex(p), _real(q), _real(k_in))
    return Retardance(value, fast_axis)


def stokes(rays: RayBatch | npt.ArrayLike, polarization: npt.ArrayLike,
           axis: npt.ArrayLike) -> FloatArray:
    """Stokes parameters (s0, s1, s2, s3) in units of |E|^2 in the basis e1 = ``axis`` projected
    perpendicular to k, e2 = k x e1; s3 > 0 is right circular (docs/architecture.md, Händigkeit
    und Stokes).

    ``stokes(rays, polarization, axis)``: of P E for every ray, (N, 4), power fractions without
    s. ``stokes(e, axis, k)``: of one field e (3,) transverse to the unit vector k, (4,)."""
    if isinstance(rays, RayBatch):
        result: FloatArray = _core.polar_stokes(rays, _complex(polarization), _real(axis))
        return result
    single: FloatArray = _core.polar_stokes_vector(_complex(rays), _real(polarization),
                                                    _real(axis))
    return single


def prt_matrix(k_in: npt.ArrayLike, k_out: npt.ArrayLike, normal: npt.ArrayLike, a_s: complex,
               a_p: complex) -> ComplexArray:
    """PRT matrix (3, 3) of one interface with amplitudes a_s, a_p (Lam, Eqs. (3.1)-(3.9)):
    P = a_s s s^T + a_p p_out p_in^T + k_out k_in^T, s = k_in x N normalised, p = k x s. Unit
    vectors k_in, k_out and normal (either orientation)."""
    result: ComplexArray = _core.polar_prt_matrix(_real(k_in), _real(k_out), _real(normal),
                                                   complex(a_s), complex(a_p))
    return result


def geometric_transform(k_in: npt.ArrayLike, k_out: npt.ArrayLike, normal: npt.ArrayLike,
                        reflection: bool) -> FloatArray:
    """Geometric transformation Q (3, 3) of one interface: s s^T +/- p_out p_in^T +
    k_out k_in^T (+ for refraction, - for reflection; #57)."""
    result: FloatArray = _core.polar_geometric_transform(_real(k_in), _real(k_out),
                                                          _real(normal), reflection)
    return result
