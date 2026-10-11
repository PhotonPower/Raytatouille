// Pupil function, its cells, refocus and the Zernike fit (ADR 0032; #202, #203) for Python
// (#209, P1a). Columns are read-only NumPy copies (analysis results, ADR 0002).

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include <complex>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

#include "bindings.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/analysis/zernike.hpp"
#include "rtt/math/types.hpp"
#include "rtt/trace/run_control.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

using analysis::PupilFunction;
using analysis::PupilMapping;

/// Read-only NumPy copy of a vector of doubles.
ReadOnlyArray<double> copy(const std::vector<double>& v) {
  return read_only_array<double>(v, [](double x) { return x; });
}

/// Read-only NumPy copy (x, y, z) of a vector.
ReadOnlyArray<double> copy(const math::Vec3& v) {
  const std::vector<double> xyz = {v.x(), v.y(), v.z()};
  return copy(xyz);
}

/// Read-only (rows, cols) array of doubles, row-major (a copy owned by the array).
nb::ndarray<nb::numpy, const double, nb::ndim<2>, nb::c_contig> matrix(std::vector<double> values,
                                                                       std::size_t rows,
                                                                       std::size_t cols) {
  auto data = std::make_unique<std::vector<double>>(std::move(values));
  const std::size_t shape[2] = {rows, cols};
  const nb::capsule owner(data.get(), [](void* q) noexcept {
    const std::unique_ptr<std::vector<double>> owned(static_cast<std::vector<double>*>(q));
  });
  const std::vector<double>* owned_by_capsule = data.release();  // freed by the capsule
  return {owned_by_capsule->data(), 2, shape, owner};
}

/// Read-only (n, 3, 3) array of the PRT matrices (9 complex per cell, row-major).
nb::ndarray<nb::numpy, const std::complex<double>, nb::ndim<3>, nb::c_contig> prt_array(
    const std::vector<std::complex<double>>& prt) {
  auto data = std::make_unique<std::vector<std::complex<double>>>(prt);
  const std::size_t shape[3] = {prt.size() / 9, 3, 3};
  const nb::capsule owner(data.get(), [](void* q) noexcept {
    const std::unique_ptr<std::vector<std::complex<double>>> owned(
        static_cast<std::vector<std::complex<double>>*>(q));
  });
  const std::vector<std::complex<double>>* owned_by_capsule = data.release();  // the capsule frees
  return {owned_by_capsule->data(), 3, shape, owner};
}

}  // namespace

void bind_pupil(nb::module_& m) {
  m.attr("MAX_PUPIL_GRID") = analysis::kMaxPupilGrid;
  m.attr("MAX_ZERNIKE_TERMS") = analysis::kMaxZernikeTerms;

  nb::class_<analysis::PupilCells>(
      m, "PupilCells",
      "N x N cells of the pupil function (ADR 0032, point 1), rows from py = -1, px fastest.")
      .def_ro("grid", &analysis::PupilCells::grid, "N.")
      .def_prop_ro(
          "fraction", [](const analysis::PupilCells& c) { return copy(c.fraction); },
          nb::rv_policy::reference, "Area fraction f in [0, 1] inside the unit circle (copy).")
      .def_prop_ro(
          "px", [](const analysis::PupilCells& c) { return copy(c.px); }, nb::rv_policy::reference,
          "Sample point x: centre, or centroid of the part inside (copy).")
      .def_prop_ro(
          "py", [](const analysis::PupilCells& c) { return copy(c.py); }, nb::rv_policy::reference,
          "Sample point y (copy).")
      .def("__len__", [](const analysis::PupilCells& c) { return c.fraction.size(); });
  m.def("pupil_cells", &analysis::pupil_cells, "grid"_a,
        "Cells of an N x N grid (ADR 0032, point 1); ValueError for N outside [2, "
        "MAX_PUPIL_GRID].");

  nb::class_<PupilMapping>(
      m, "PupilMapping",
      "Affine mapping s ~ s0 + M p of the pupil coordinates onto the direction sines (ADR 0032, "
      "point 4).")
      .def_prop_ro(
          "m",
          [](const PupilMapping& p) {
            return matrix(std::vector<double>(p.m.begin(), p.m.end()), 2, 2);
          },
          nb::rv_policy::reference, "M (2 x 2): rows s_x, s_y, columns px, py (copy).")
      .def_prop_ro(
          "s0",
          [](const PupilMapping& p) { return copy(std::vector<double>(p.s0.begin(), p.s0.end())); },
          nb::rv_policy::reference,
          "s of the chief ray of the evaluated wavelength; exactly 0 at the reference (copy).")
      .def_prop_ro(
          "centre", [](const PupilMapping& p) { return copy(p.centre); }, nb::rv_policy::reference,
          "Origin of the PSF plane C' (global, mm; copy).")
      .def_prop_ro(
          "e_x", [](const PupilMapping& p) { return copy(p.e_x); }, nb::rv_policy::reference,
          "Unit x axis of the PSF plane (global; copy).")
      .def_prop_ro(
          "e_y", [](const PupilMapping& p) { return copy(p.e_y); }, nb::rv_policy::reference,
          "Unit y axis, d_c x e_x (global; copy).")
      .def_prop_ro(
          "d_c", [](const PupilMapping& p) { return copy(p.d_c); }, nb::rv_policy::reference,
          "Chief ray of the reference wavelength on the image surface (global; copy).")
      .def_ro("delta_max", &PupilMapping::delta_max,
              "max |s - s0 - M p| / ||M||_2 over the arrived cells, pupil units.")
      .def_ro("delta_rms", &PupilMapping::delta_rms, "The same, RMS weighted with f.");

  nb::class_<PupilFunction>(
      m, "PupilFunction",
      "Pupil function of one field at one wavelength (ADR 0032, point 3). Check fraction > 0 "
      "first, then status: cells with f = 0 are not traced (status VIGNETTED, A = W = 0).")
      .def_prop_ro(
          "path", [](const PupilFunction& p) { return p.path.index; }, "Path index.")
      .def_ro("field", &PupilFunction::field)
      .def_ro("wavelength", &PupilFunction::wavelength)
      .def_ro("wavelength_um", &PupilFunction::wavelength_um, "Evaluated wavelength, um.")
      .def_ro("reference_wavelength_um", &PupilFunction::reference_wavelength_um,
              "Reference wavelength (unit of W), um.")
      .def_ro("grid", &PupilFunction::grid, "N.")
      .def_ro("focus_shift", &PupilFunction::focus_shift, "mm along the chief ray.")
      .def_ro("apodization", &PupilFunction::apodization, "Gaussian factor G.")
      .def_ro("lost_warning_fraction", &PupilFunction::lost_warning_fraction)
      .def_ro("mapping_warning", &PupilFunction::mapping_warning)
      .def_ro("sphere", &PupilFunction::sphere, "Reference sphere after focus_shift.")
      .def_prop_ro(
          "fraction", [](const PupilFunction& p) { return copy(p.fraction); },
          nb::rv_policy::reference, "Area fraction f per cell (copy).")
      .def_prop_ro(
          "px", [](const PupilFunction& p) { return copy(p.px); }, nb::rv_policy::reference,
          "Sample point x per cell (copy).")
      .def_prop_ro(
          "py", [](const PupilFunction& p) { return copy(p.py); }, nb::rv_policy::reference,
          "Sample point y per cell (copy).")
      .def_prop_ro(
          "amplitude", [](const PupilFunction& p) { return copy(p.amplitude); },
          nb::rv_policy::reference, "sqrt(weight) exp(-G rho^2); 0 if not arrived or f = 0 (copy).")
      .def_prop_ro(
          "w", [](const PupilFunction& p) { return copy(p.w); }, nb::rv_policy::reference,
          "W = OPL_ref - OPL, waves at the reference wavelength, > 0 leading (copy).")
      .def_prop_ro(
          "status",
          [](const PupilFunction& p) {
            return read_only_array<std::uint8_t>(
                p.status, [](trace::RayStatus s) { return status_value(s); });
          },
          nb::rv_policy::reference, "RayStatus values as uint8 (copy).")
      .def_prop_ro(
          "prt", [](const PupilFunction& p) { return prt_array(p.prt); }, nb::rv_policy::reference,
          "Accumulated PRT matrix per cell, shape (N^2, 3, 3), global (with_prt; else (0, 3, 3); "
          "copy).")
      .def_ro("mapping", &PupilFunction::mapping, "None for an afocal image space.")
      .def_ro("arrived_fraction", &PupilFunction::arrived_fraction,
              "sum f over the arrived cells / sum f over all cells.")
      .def_ro("losses", &PupilFunction::losses, "Traced cells by final status.")
      .def_ro("warnings", &PupilFunction::warnings,
              "rays.lost, stop.clips_beam, pupil.mapping_nonlinear (ADR 0023); "
              "raytatouille.analysis also issues them as RaytatouilleWarning.");

  m.def(
      "pupil_function",
      [](const compile::CompiledSystem& s, const PathArg& path, std::uint16_t field,
         std::optional<std::uint16_t> wavelength, int grid, double apodization, double focus_shift,
         trace::Aiming aiming, bool with_prt, double lost_warning_fraction, double mapping_warning,
         std::optional<int> threads, const std::optional<trace::CancelToken>& cancel,
         const std::optional<nb::callable>& progress) {
        const trace::RunControl control = run_control(cancel, progress);
        analysis::PupilFunctionOptions options;
        options.grid = grid;
        options.apodization = apodization;
        options.focus_shift = focus_shift;
        options.aiming = aiming;
        options.with_prt = with_prt;
        options.lost_warning_fraction = lost_warning_fraction;
        options.mapping_warning = mapping_warning;
        return released(threads, [&] {
          return analysis::pupil_function(s, path_id(s, path), field,
                                          wavelength_index(s, wavelength), options, control);
        });
      },
      "system"_a, "path"_a, "field"_a, "wavelength"_a.none(), "grid"_a, "apodization"_a,
      "focus_shift"_a, "aiming"_a, "with_prt"_a, "lost_warning_fraction"_a, "mapping_warning"_a,
      "threads"_a.none(), "cancel"_a.none(), "progress"_a.none(),
      "Pupil function of `field`; wavelength None means the reference.");
  m.def(
      "refocus",
      [](const compile::CompiledSystem& s, const PupilFunction& pupil, double focus_shift) {
        return analysis::refocus(s, pupil, focus_shift);
      },
      "system"_a, "pupil"_a, "focus_shift"_a, nb::call_guard<nb::gil_scoped_release>(),
      "The pupil function for another focus_shift from its end states, without tracing again; "
      "`system` must be the one that produced it.");

  nb::class_<analysis::NollIndex>(m, "NollIndex",
                                  "(n, m) of a Noll index; m >= 0, sine for odd j (Niu/Tian).")
      .def_ro("n", &analysis::NollIndex::n)
      .def_ro("m", &analysis::NollIndex::m);
  m.def("noll_to_nm", &analysis::noll_to_nm, "j"_a, "(n, m) of the Noll index j (Eq. (8)).");
  m.def("noll_from_nm", &analysis::noll_from_nm, "n"_a, "m"_a, "sine"_a,
        "Noll index of (n, m), sine or cosine (Eq. (9)).");
  m.def("zernike_radial", &analysis::zernike_radial, "n"_a, "m"_a, "rho"_a,
        "R_n^m(rho) (Eq. (4)), n <= 20.");
  m.def("zernike", &analysis::zernike, "j"_a, "rho"_a, "theta"_a,
        "Orthonormal Z_j(rho, theta) after Noll; theta from +x counter-clockwise.");

  nb::class_<analysis::ZernikeFit>(
      m, "ZernikeFit", "Zernike fit of a pupil function and the Strehl approximation (ADR 0032).")
      .def_prop_ro("path", [](const analysis::ZernikeFit& f) { return f.path.index; })
      .def_ro("field", &analysis::ZernikeFit::field)
      .def_ro("wavelength", &analysis::ZernikeFit::wavelength)
      .def_ro("wavelength_um", &analysis::ZernikeFit::wavelength_um)
      .def_ro("reference_wavelength_um", &analysis::ZernikeFit::reference_wavelength_um)
      .def_prop_ro(
          "coefficients", [](const analysis::ZernikeFit& f) { return copy(f.coefficients); },
          nb::rv_policy::reference, "c_j at index j - 1, waves at the reference (copy).")
      .def_ro("residual_rms", &analysis::ZernikeFit::residual_rms)
      .def_ro("sigma", &analysis::ZernikeFit::sigma, "Std of W (piston removed), waves.")
      .def_ro("sigma_bf", &analysis::ZernikeFit::sigma_bf,
              "Std about the diffraction focus (Z1 to Z4 removed), waves at the reference.")
      .def_ro("strehl", &analysis::ZernikeFit::strehl,
              "exp(-(2 pi sigma_bf lambda_ref / lambda)^2) (Marechal, good to about 0.1).")
      .def_ro("cells", &analysis::ZernikeFit::cells)
      .def_ro("arrived_fraction", &analysis::ZernikeFit::arrived_fraction)
      .def_ro("warnings", &analysis::ZernikeFit::warnings, "zernike.pupil_incomplete.");
  m.def(
      "zernike_fit",
      [](const PupilFunction& pupil, int terms, double incomplete_warning_fraction) {
        analysis::ZernikeOptions options;
        options.terms = terms;
        options.incomplete_warning_fraction = incomplete_warning_fraction;
        return analysis::zernike_fit(pupil, options);
      },
      "pupil"_a, "terms"_a, "incomplete_warning_fraction"_a,
      nb::call_guard<nb::gil_scoped_release>(),
      "Zernike fit with `terms` Noll terms, sigma, sigma_bf and the Strehl approximation.");
}

}  // namespace rtt::py
