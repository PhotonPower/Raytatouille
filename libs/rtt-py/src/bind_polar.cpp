// Polarization quantities for raytatouille.polar (#62). The functions here take exact dtypes
// and C-contiguous arrays; python/raytatouille/polar.py converts the user's input and wraps the
// results. Inputs are checked here (API boundary, ADR 0009): shapes, finite values, unit vectors.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/tuple.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "bindings.hpp"
#include "polar_batch.hpp"
#include "polar_matrix.hpp"
#include "rtt/polar/ideal.hpp"
#include "rtt/polar/prt_analysis.hpp"
#include "rtt/trace/ray_batch.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

using math::CMat3;
using math::Complex;
using math::Mat3;
using math::Vec3;
using polar_batch::CVec3;
using trace::RayBatch;

using RealIn = nb::ndarray<const double, nb::c_contig, nb::device::cpu>;
using ComplexIn = nb::ndarray<const Complex, nb::c_contig, nb::device::cpu>;
template <typename T>
using Out = nb::ndarray<nb::numpy, T>;

/// Tolerance of the unit-vector checks (as polar_batch.hpp).
constexpr double kUnit = polar_batch::kStateTolerance;

/// NumPy array of the given shape that owns `data` (a copy of the results).
template <typename T>
Out<T> owned(std::vector<T>&& values, std::initializer_list<std::size_t> shape) {
  auto data = std::make_unique<std::vector<T>>(std::move(values));
  T* pointer = data->data();
  // The capsule owns the copy from here on and frees it with the array.
  nb::capsule owner(data.get(), [](void* p) noexcept {
    const std::unique_ptr<std::vector<T>> owned_data(static_cast<std::vector<T>*>(p));
  });
  static_cast<void>(data.release());
  return Out<T>(pointer, shape.size(), shape.begin(), owner);
}

Out<double> vec3_rows(const std::vector<Vec3>& v) {
  std::vector<double> data;
  data.reserve(3 * v.size());
  for (const Vec3& x : v) data.insert(data.end(), {x.x(), x.y(), x.z()});
  return owned(std::move(data), {v.size(), 3});
}

Out<Complex> cvec3_rows(const std::vector<CVec3>& v) {
  std::vector<Complex> data;
  data.reserve(3 * v.size());
  for (const CVec3& x : v) data.insert(data.end(), {x.x(), x.y(), x.z()});
  return owned(std::move(data), {v.size(), 3});
}

Vec3 vec3(const RealIn& a, const char* name) {
  if (a.ndim() != 1 || a.shape(0) != 3) {
    throw std::invalid_argument(std::string(name) + " must have 3 elements");
  }
  const Vec3 v(a.data()[0], a.data()[1], a.data()[2]);
  if (!v.allFinite()) throw std::invalid_argument(std::string(name) + " is not finite");
  return v;
}

Vec3 unit(const RealIn& a, const char* name) {
  const Vec3 v = vec3(a, name);
  if (!(std::abs(v.norm() - 1.0) <= kUnit)) {
    throw std::invalid_argument(std::string(name) + " is not a unit vector (within 1e-12)");
  }
  return v;
}

CVec3 cvec3(const ComplexIn& a, const char* name) {
  if (a.ndim() != 1 || a.shape(0) != 3) {
    throw std::invalid_argument(std::string(name) + " must have 3 elements");
  }
  const CVec3 v(a.data()[0], a.data()[1], a.data()[2]);
  if (!v.allFinite()) throw std::invalid_argument(std::string(name) + " is not finite");
  return v;
}

/// (3,) or (N, 3) complex states.
std::vector<CVec3> states(const ComplexIn& a) {
  const bool one = a.ndim() == 1 && a.shape(0) == 3;
  const bool rows = a.ndim() == 2 && a.shape(1) == 3;
  if (!one && !rows) {
    throw std::invalid_argument("polarization must have shape (3,) or (N, 3) for N rays");
  }
  const std::size_t n = one ? 1 : a.shape(0);
  std::vector<CVec3> s(n);
  for (std::size_t i = 0; i < n; ++i) {
    s[i] = CVec3(a.data()[3 * i], a.data()[3 * i + 1], a.data()[3 * i + 2]);
  }
  return s;
}

template <typename Array, typename Matrix>
Matrix mat3(const Array& a, const char* name) {
  if (a.ndim() != 2 || a.shape(0) != 3 || a.shape(1) != 3) {
    throw std::invalid_argument(std::string(name) + " must be a 3x3 matrix");
  }
  Matrix m;
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) m(r, c) = a.data()[3 * r + c];
  }
  if (!m.allFinite()) throw std::invalid_argument(std::string(name) + " is not finite");
  return m;
}

Out<Complex> cmat3_out(const CMat3& m) {
  std::vector<Complex> data(9);
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) data[static_cast<std::size_t>(3 * r + c)] = m(r, c);
  }
  return owned(std::move(data), {3, 3});
}

Out<Complex> cvec3_out(const CVec3& v) {
  return owned(std::vector<Complex>{v.x(), v.y(), v.z()}, {3});
}

using DiattenuationTuple = std::tuple<double, double, double, Out<Complex>>;
using RetardanceTuple = std::tuple<double, Out<Complex>>;

RetardanceTuple retardance_tuple(const rtt::polar::Retardance& r) {
  return {r.value, cvec3_out(r.fast_axis)};
}

}  // namespace

void bind_polar(nb::module_& m) {
  // Batch functions on a traced RayBatch (polar_batch.hpp).
  m.def(
      "polar_initial_directions",
      [](const RayBatch& rays) { return vec3_rows(polar_batch::initial_directions(rays)); },
      "rays"_a);
  m.def(
      "polar_transverse_polarization",
      [](const RayBatch& rays, const ComplexIn& e) {
        return cvec3_rows(polar_batch::transverse_polarization(rays, cvec3(e, "polarization")));
      },
      "rays"_a, "polarization"_a);
  m.def(
      "polar_transmission",
      [](const RayBatch& rays, const std::optional<ComplexIn>& polarization) {
        if (!polarization) return owned(polar_batch::transmission(rays), {rays.size()});
        return owned(polar_batch::transmission(rays, states(*polarization)), {rays.size()});
      },
      "rays"_a, "polarization"_a.none());
  m.def(
      "polar_diattenuation",
      [](const RayBatch& rays) {
        polar_batch::Diattenuations d = polar_batch::diattenuation(rays);
        const std::size_t n = rays.size();
        return std::make_tuple(owned(std::move(d.value), {n}), owned(std::move(d.maximum), {n}),
                               owned(std::move(d.minimum), {n}), cvec3_rows(d.axis));
      },
      "rays"_a);
  m.def(
      "polar_retardance",
      [](const RayBatch& rays) {
        polar_batch::Retardances r = polar_batch::retardance(rays);
        return std::make_tuple(owned(std::move(r.value), {rays.size()}), cvec3_rows(r.fast_axis));
      },
      "rays"_a);
  m.def(
      "polar_stokes",
      [](const RayBatch& rays, const ComplexIn& polarization, const RealIn& axis) {
        const std::vector<std::array<double, 4>> s =
            polar_batch::stokes(rays, states(polarization), vec3(axis, "axis"));
        std::vector<double> data;
        data.reserve(4 * s.size());
        for (const auto& row : s) data.insert(data.end(), row.begin(), row.end());
        return owned(std::move(data), {s.size(), 4});
      },
      "rays"_a, "polarization"_a, "axis"_a);

  // Single matrices (rtt/polar/prt.hpp, prt_analysis.hpp, ideal.hpp).
  m.def(
      "polar_prt_matrix",
      [](const RealIn& k_in, const RealIn& k_out, const RealIn& normal, Complex a_s, Complex a_p) {
        if (!std::isfinite(a_s.real()) || !std::isfinite(a_s.imag()) ||
            !std::isfinite(a_p.real()) || !std::isfinite(a_p.imag())) {
          throw std::invalid_argument("a_s and a_p must be finite");
        }
        return cmat3_out(polar_matrix::prt_matrix(unit(k_in, "k_in"), unit(k_out, "k_out"),
                                                  unit(normal, "normal"), a_s, a_p));
      },
      "k_in"_a, "k_out"_a, "normal"_a, "a_s"_a, "a_p"_a);
  m.def(
      "polar_geometric_transform",
      [](const RealIn& k_in, const RealIn& k_out, const RealIn& normal, bool reflection) {
        const Mat3 q = polar_matrix::geometric_transform(unit(k_in, "k_in"), unit(k_out, "k_out"),
                                                         unit(normal, "normal"), reflection);
        std::vector<double> data(9);
        for (int r = 0; r < 3; ++r) {
          for (int c = 0; c < 3; ++c) data[static_cast<std::size_t>(3 * r + c)] = q(r, c);
        }
        return owned(std::move(data), {3, 3});
      },
      "k_in"_a, "k_out"_a, "normal"_a, "reflection"_a);
  m.def(
      "polar_diattenuation_matrix",
      [](const ComplexIn& p, const RealIn& k_in, const RealIn& k_out) -> DiattenuationTuple {
        const CMat3 pm = mat3<ComplexIn, CMat3>(p, "p");
        const Vec3 ki = unit(k_in, "k_in");
        const Vec3 ko = unit(k_out, "k_out");
        if (!((pm * ki.cast<Complex>() - ko.cast<Complex>()).norm() <= kUnit)) {
          throw std::invalid_argument("p must map k_in onto k_out (within 1e-12)");
        }
        const rtt::polar::Diattenuation d = rtt::polar::diattenuation(pm, ki, ko);
        return {d.value, d.maximum, d.minimum, cvec3_out(d.axis)};
      },
      "p"_a, "k_in"_a, "k_out"_a);
  m.def(
      "polar_retardance_matrix",
      [](const ComplexIn& matrix, const RealIn& k) -> RetardanceTuple {
        const CMat3 mm = mat3<ComplexIn, CMat3>(matrix, "m");
        const Vec3 kk = unit(k, "k");
        if (!((mm * kk.cast<Complex>() - kk.cast<Complex>()).norm() <= kUnit)) {
          throw std::invalid_argument("m must map k onto itself (within 1e-12)");
        }
        return retardance_tuple(rtt::polar::retardance(mm, kk));
      },
      "m"_a, "k"_a);
  m.def(
      "polar_physical_retardance",
      [](const ComplexIn& p, const RealIn& q, const RealIn& k_in) -> RetardanceTuple {
        const CMat3 pm = mat3<ComplexIn, CMat3>(p, "p");
        const Mat3 qm = mat3<RealIn, Mat3>(q, "q");
        const Vec3 ki = unit(k_in, "k_in");
        if (!((pm * ki.cast<Complex>() - (qm * ki).cast<Complex>()).norm() <= kUnit)) {
          throw std::invalid_argument("p and q must map k_in onto the same k_out (within 1e-12)");
        }
        return retardance_tuple(rtt::polar::physical_retardance(pm, qm, ki));
      },
      "p"_a, "q"_a, "k_in"_a);
  m.def(
      "polar_stokes_vector",
      [](const ComplexIn& e, const RealIn& axis, const RealIn& k) {
        const CVec3 ev = cvec3(e, "e");
        const Vec3 av = vec3(axis, "axis");
        const Vec3 kv = unit(k, "k");
        if (!(std::abs(kv.cast<Complex>().dot(ev)) <= kUnit * std::max(ev.norm(), 1.0))) {
          throw std::invalid_argument("e is not transverse to k (within 1e-12)");
        }
        const Vec3 projected = av - av.dot(kv) * kv;  // as rtt::polar::transverse_axis()
        if (!(projected.norm() >= rtt::polar::kAxisAlongK * av.norm()) || av.norm() == 0.0) {
          throw std::invalid_argument("axis is parallel to k");
        }
        const rtt::polar::Stokes<double> s = rtt::polar::stokes<double>(ev, av, kv);
        return owned(std::vector<double>{s.s0, s.s1, s.s2, s.s3}, {4});
      },
      "e"_a, "axis"_a, "k"_a);
}

}  // namespace rtt::py
