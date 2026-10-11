// Jones pupil (ADR 0034, points 5 and 6; #208) for Python (#208 B2). Columns are read-only
// NumPy copies (analysis results, ADR 0002).

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
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
#include "rtt/analysis/jones.hpp"
#include "rtt/analysis/stokes.hpp"
#include "rtt/math/types.hpp"
#include "rtt/trace/run_control.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

using analysis::JonesBasis;
using analysis::JonesPupil;

/// Read-only NumPy copy of a vector of doubles.
ReadOnlyArray<double> copy(const std::vector<double>& v) {
  return read_only_array<double>(v, [](double x) { return x; });
}

/// Read-only NumPy copy (x, y, z) of a vector.
ReadOnlyArray<double> copy(const math::Vec3& v) {
  const std::vector<double> xyz = {v.x(), v.y(), v.z()};
  return copy(xyz);
}

/// Read-only (n, 2, 2) array of the Jones matrices (4 complex per point, row-major: rows the
/// output basis, columns the input basis).
nb::ndarray<nb::numpy, const std::complex<double>, nb::ndim<3>, nb::c_contig> jones_array(
    const std::vector<std::complex<double>>& j) {
  auto data = std::make_unique<std::vector<std::complex<double>>>(j);
  const std::size_t shape[3] = {j.size() / 4, 2, 2};
  const nb::capsule owner(data.get(), [](void* q) noexcept {
    const std::unique_ptr<std::vector<std::complex<double>>> owned(
        static_cast<std::vector<std::complex<double>>*>(q));
  });
  const std::vector<std::complex<double>>* owned_by_capsule = data.release();  // the capsule frees
  return {owned_by_capsule->data(), 3, shape, owner};
}

/// Options of both forms; sampling and aiming apply to the convenience form only.
analysis::JonesPupilOptions jones_options(JonesBasis basis,
                                          bool remove_global_phase,
                                          bool field_amplitude,
                                          double lost_warning_fraction,
                                          const trace::PupilSampling& sampling,
                                          trace::Aiming aiming) {
  analysis::JonesPupilOptions options;
  options.sampling = sampling;
  options.aiming = aiming;
  options.basis = basis;
  options.remove_global_phase = remove_global_phase;
  options.field_amplitude = field_amplitude;
  options.lost_warning_fraction = lost_warning_fraction;
  return options;
}

}  // namespace

void bind_jones(nb::module_& m) {
  nb::enum_<JonesBasis>(m, "JonesBasis", nb::is_arithmetic(),
                        "Output basis of the Jones pupil (ADR 0034, point 5).")
      .value("DOUBLE_POLE", JonesBasis::DoublePole,
             "R(k_c -> k) (x'_c, y'_c), singular only at -k_c (default).")
      .value("DIPOLE", JonesBasis::Dipole, "Lam's dipole basis with the poles at +-y'_c.");

  nb::class_<JonesPupil>(
      m, "JonesPupil",
      "Jones pupil of one field at one wavelength (ADR 0034, points 5 and 6): per pupil point "
      "the 2 x 2 Jones matrix J = exp(-2 pi i w lambda_ref / lambda) B_out^T P B_in and the maps. "
      "Check status first: a point that did not arrive has w = 0 and NaN in j and the maps; a "
      "point at a pole of a basis is ALIVE with NaN.")
      .def_prop_ro(
          "path", [](const JonesPupil& p) { return p.path.index; }, "Path index.")
      .def_ro("field", &JonesPupil::field)
      .def_ro("wavelength", &JonesPupil::wavelength)
      .def_ro("wavelength_um", &JonesPupil::wavelength_um, "Evaluated wavelength, um.")
      .def_ro("reference_wavelength_um", &JonesPupil::reference_wavelength_um,
              "Reference wavelength (unit of w), um.")
      .def_ro("basis", &JonesPupil::basis)
      .def_ro("remove_global_phase", &JonesPupil::remove_global_phase)
      .def_ro("field_amplitude", &JonesPupil::field_amplitude,
              "True: field amplitudes in the image-space medium; False: power-normalised.")
      .def_ro("global_phase", &JonesPupil::global_phase,
              "psi removed from every J, rad; 0 without remove_global_phase.")
      .def_prop_ro(
          "x_in", [](const JonesPupil& p) { return copy(p.x_in); }, nb::rv_policy::reference,
          "Chief frame at the input, x (global; copy).")
      .def_prop_ro(
          "y_in", [](const JonesPupil& p) { return copy(p.y_in); }, nb::rv_policy::reference,
          "Chief frame at the input, y (global; copy).")
      .def_prop_ro(
          "k_in", [](const JonesPupil& p) { return copy(p.k_in); }, nb::rv_policy::reference,
          "Start direction of the chief ray k0_c (global; copy).")
      .def_prop_ro(
          "x_out", [](const JonesPupil& p) { return copy(p.x_out); }, nb::rv_policy::reference,
          "Chief frame at the output, the x axis of the PSF plane (global; copy).")
      .def_prop_ro(
          "y_out", [](const JonesPupil& p) { return copy(p.y_out); }, nb::rv_policy::reference,
          "k_out x x_out (global; copy).")
      .def_prop_ro(
          "k_out", [](const JonesPupil& p) { return copy(p.k_out); }, nb::rv_policy::reference,
          "Direction of the chief ray on the image surface k_c (global; copy).")
      .def_ro("sphere", &JonesPupil::sphere,
              "Reference sphere of w; radius NaN in the form with start rays (no sphere).")
      .def_prop_ro(
          "px", [](const JonesPupil& p) { return copy(p.px); }, nb::rv_policy::reference,
          "Pupil coordinate x per point (copy).")
      .def_prop_ro(
          "py", [](const JonesPupil& p) { return copy(p.py); }, nb::rv_policy::reference,
          "Pupil coordinate y per point (copy).")
      .def_prop_ro(
          "status",
          [](const JonesPupil& p) {
            return read_only_array<std::uint8_t>(
                p.status, [](trace::RayStatus s) { return status_value(s); });
          },
          nb::rv_policy::reference, "RayStatus values as uint8 (copy).")
      .def_prop_ro(
          "w", [](const JonesPupil& p) { return copy(p.w); }, nb::rv_policy::reference,
          "Phase in waves at the reference wavelength, the sign of W; 0 if not arrived (copy).")
      .def_prop_ro(
          "j", [](const JonesPupil& p) { return jones_array(p.j); }, nb::rv_policy::reference,
          "Jones matrix per point, shape (n, 2, 2): rows the output, columns the input basis "
          "(copy).")
      .def_prop_ro(
          "transmission", [](const JonesPupil& p) { return copy(p.transmission); },
          nb::rv_policy::reference, "Power for an unpolarized source (weight; copy).")
      .def_prop_ro(
          "diattenuation", [](const JonesPupil& p) { return copy(p.diattenuation); },
          nb::rv_policy::reference, "D in [0, 1] (copy).")
      .def_prop_ro(
          "diattenuation_orientation",
          [](const JonesPupil& p) { return copy(p.diattenuation_orientation); },
          nb::rv_policy::reference,
          "Orientation psi in [0, pi) of the state of largest transmission, input basis, rad.")
      .def_prop_ro(
          "diattenuation_ellipticity",
          [](const JonesPupil& p) { return copy(p.diattenuation_ellipticity); },
          nb::rv_policy::reference,
          "Ellipticity chi in [-pi/4, pi/4] of that state, > 0 right-handed, rad.")
      .def_prop_ro(
          "retardance", [](const JonesPupil& p) { return copy(p.retardance); },
          nb::rv_policy::reference, "Physical retardance delta in [0, pi] from Q^T P, rad (copy).")
      .def_prop_ro(
          "retardance_orientation",
          [](const JonesPupil& p) { return copy(p.retardance_orientation); },
          nb::rv_policy::reference, "Orientation psi of the fast state, input basis, rad.")
      .def_prop_ro(
          "retardance_ellipticity",
          [](const JonesPupil& p) { return copy(p.retardance_ellipticity); },
          nb::rv_policy::reference, "Ellipticity chi of the fast state, rad.")
      .def_ro("losses", &JonesPupil::losses, "Launched rays by the status of their point.")
      .def_ro("warnings", &JonesPupil::warnings,
              "rays.lost, stop.clips_beam, pupil.basis_pole, pupil.field_amplitude_absorbing, "
              "pupil.chief_phase_unavailable (ADR 0023, 0034); raytatouille.analysis also issues "
              "them as RaytatouilleWarning.");

  m.def(
      "jones_pupil",
      [](const compile::CompiledSystem& s, const PathArg& path, std::uint16_t field,
         std::optional<std::uint16_t> wavelength, const trace::PupilSampling& sampling,
         trace::Aiming aiming, JonesBasis basis, bool remove_global_phase, bool field_amplitude,
         double lost_warning_fraction, std::optional<int> threads,
         const std::optional<trace::CancelToken>& cancel,
         const std::optional<nb::callable>& progress) {
        const trace::RunControl control = run_control(cancel, progress);
        const analysis::JonesPupilOptions options = jones_options(
            basis, remove_global_phase, field_amplitude, lost_warning_fraction, sampling, aiming);
        return released(threads, [&] {
          return analysis::jones_pupil(s, path_id(s, path), field, wavelength_index(s, wavelength),
                                       options, control);
        });
      },
      "system"_a, "path"_a, "field"_a, "wavelength"_a.none(), "sampling"_a, "aiming"_a, "basis"_a,
      "remove_global_phase"_a, "field_amplitude"_a, "lost_warning_fraction"_a, "threads"_a.none(),
      "cancel"_a.none(), "progress"_a.none(),
      "Jones pupil of `field` (convenience form); wavelength None means the reference.");
  m.def(
      "jones_pupil_rays",
      [](const compile::CompiledSystem& s, const PathArg& path, const trace::RayBatch& start,
         std::size_t chief, JonesBasis basis, bool remove_global_phase, bool field_amplitude,
         double lost_warning_fraction, std::optional<int> threads,
         const std::optional<trace::CancelToken>& cancel,
         const std::optional<nb::callable>& progress) {
        const trace::RunControl control = run_control(cancel, progress);
        const analysis::JonesPupilOptions options =
            jones_options(basis, remove_global_phase, field_amplitude, lost_warning_fraction,
                          trace::GridPupil{33}, trace::Aiming::Real);
        return released(threads, [&] {
          return analysis::jones_pupil(s, path_id(s, path), start, chief, options, control);
        });
      },
      "system"_a, "path"_a, "start"_a, "chief"_a, "basis"_a, "remove_global_phase"_a,
      "field_amplitude"_a, "lost_warning_fraction"_a, "threads"_a.none(), "cancel"_a.none(),
      "progress"_a.none(), "Jones pupil of the start rays `start` with chief ray `chief`.");

  nb::class_<analysis::StokesPupil>(
      m, "StokesPupil",
      "Stokes parameters of E_out = J E_in per point of a Jones pupil (ADR 0034, point 6, "
      "addendum #211): power of the polarized input, s3 > 0 right circular; NaN without J.")
      .def_prop_ro(
          "path", [](const analysis::StokesPupil& p) { return p.path.index; }, "Path index.")
      .def_ro("field", &analysis::StokesPupil::field)
      .def_ro("wavelength", &analysis::StokesPupil::wavelength)
      .def_ro("wavelength_um", &analysis::StokesPupil::wavelength_um)
      .def_ro("e1", &analysis::StokesPupil::e1, "Input component along x_in, normalised.")
      .def_ro("e2", &analysis::StokesPupil::e2, "Input component along y_in, normalised.")
      .def_prop_ro(
          "x_in", [](const analysis::StokesPupil& p) { return copy(p.x_in); },
          nb::rv_policy::reference, "Input chief frame, x (global; copy).")
      .def_prop_ro(
          "y_in", [](const analysis::StokesPupil& p) { return copy(p.y_in); },
          nb::rv_policy::reference, "Input chief frame, y (global; copy).")
      .def_prop_ro(
          "x_out", [](const analysis::StokesPupil& p) { return copy(p.x_out); },
          nb::rv_policy::reference, "Output chief frame, x (global; copy).")
      .def_prop_ro(
          "y_out", [](const analysis::StokesPupil& p) { return copy(p.y_out); },
          nb::rv_policy::reference, "Output chief frame, y (global; copy).")
      .def_prop_ro(
          "px", [](const analysis::StokesPupil& p) { return copy(p.px); }, nb::rv_policy::reference,
          "Pupil coordinate x per point (copy).")
      .def_prop_ro(
          "py", [](const analysis::StokesPupil& p) { return copy(p.py); }, nb::rv_policy::reference,
          "Pupil coordinate y per point (copy).")
      .def_prop_ro(
          "status",
          [](const analysis::StokesPupil& p) {
            return read_only_array<std::uint8_t>(
                p.status, [](trace::RayStatus s) { return status_value(s); });
          },
          nb::rv_policy::reference, "RayStatus values as uint8 (copy).")
      .def_prop_ro(
          "s0", [](const analysis::StokesPupil& p) { return copy(p.s0); }, nb::rv_policy::reference,
          "Power of the polarized input (copy).")
      .def_prop_ro(
          "s1", [](const analysis::StokesPupil& p) { return copy(p.s1); }, nb::rv_policy::reference,
          "|a|^2 - |b|^2, scaled as s0 (copy).")
      .def_prop_ro(
          "s2", [](const analysis::StokesPupil& p) { return copy(p.s2); }, nb::rv_policy::reference,
          "2 Re(a* b), scaled as s0 (copy).")
      .def_prop_ro(
          "s3", [](const analysis::StokesPupil& p) { return copy(p.s3); }, nb::rv_policy::reference,
          "2 Im(a b*), > 0 right circular, scaled as s0 (copy).")
      .def_prop_ro(
          "orientation", [](const analysis::StokesPupil& p) { return copy(p.orientation); },
          nb::rv_policy::reference, "Orientation psi in [0, pi) of the output ellipse, rad.")
      .def_prop_ro(
          "ellipticity", [](const analysis::StokesPupil& p) { return copy(p.ellipticity); },
          nb::rv_policy::reference, "Ellipticity chi in [-pi/4, pi/4], > 0 right-handed, rad.")
      .def_ro("warnings", &analysis::StokesPupil::warnings, "The warnings of the Jones pupil.");
  // STUB for the red run of #211 B: replaced by the call in the next commit.
  m.def(
      "stokes_pupil",
      [](const JonesPupil& /*jones*/, std::complex<double> /*e1*/, std::complex<double> /*e2*/) {
        return analysis::StokesPupil{};
      },
      "jones"_a, "e1"_a, "e2"_a, "Stokes map of `jones` for the input state (e1, e2).");
}

}  // namespace rtt::py
