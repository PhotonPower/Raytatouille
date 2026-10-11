// Diffraction PSF, through-focus Strehl, MTF and encircled energy (ADR 0033; #205, #206) for
// Python (#209, P1b). Arrays are read-only NumPy copies (analysis results, ADR 0002);
// raytatouille.analysis wraps the functions (system or pupil function, warnings).

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "bindings.hpp"
#include "rtt/analysis/psf.hpp"
#include "rtt/analysis/pupil.hpp"
#include "rtt/math/types.hpp"
#include "rtt/trace/run_control.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

using analysis::EeCentre;
using analysis::EncircledEnergy;
using analysis::Mtf;
using analysis::MtfCurve;
using analysis::Psf;
using analysis::PsfPeak;
using analysis::PsfWavelength;
using analysis::PupilFunction;
using analysis::ThroughFocus;
using Cx = std::complex<double>;

/// The wavelengths of a result as columns (PsfWavelengths).
struct PsfWavelengths {
  std::vector<PsfWavelength> points;
};

template <typename T, std::size_t D>
using Array = nb::ndarray<nb::numpy, const T, nb::ndim<D>, nb::c_contig>;

/// Read-only NumPy array of the given shape that owns a copy of `values` (row-major).
template <typename T, std::size_t D>
Array<T, D> owned(std::vector<T> values, const std::array<std::size_t, D>& shape) {
  auto data = std::make_unique<std::vector<T>>(std::move(values));
  const nb::capsule owner(data.get(), [](void* q) noexcept {
    const std::unique_ptr<std::vector<T>> freed(static_cast<std::vector<T>*>(q));
  });
  const std::vector<T>* by_capsule = data.release();  // freed by the capsule
  return {by_capsule->data(), D, shape.data(), owner};
}

Array<double, 1> copy(const std::vector<double>& v) {
  return owned<double, 1>(v, {v.size()});
}

Array<double, 1> copy(const math::Vec3& v) {
  return owned<double, 1>({v.x(), v.y(), v.z()}, {3});
}

Array<Cx, 1> copy(const std::vector<Cx>& v) {
  return owned<Cx, 1>(v, {v.size()});
}

/// A 2 x 2 matrix, row-major.
Array<double, 2> matrix2(const std::array<double, 4>& m) {
  return owned<double, 2>(std::vector<double>(m.begin(), m.end()), {2, 2});
}

/// A square side x side array of a row-major vector; (0, 0) for an empty vector.
template <typename T>
Array<T, 2> square(const std::vector<T>& v, std::size_t side) {
  if (v.empty()) return owned<T, 2>({}, {0, 0});
  return owned<T, 2>(v, {side, side});
}

/// Pupil-function options of the system forms (ADR 0032, point 3).
analysis::PupilFunctionOptions pupil_options(int grid,
                                             double apodization,
                                             double focus_shift,
                                             trace::Aiming aiming,
                                             double lost_warning_fraction,
                                             double mapping_warning) {
  analysis::PupilFunctionOptions o;
  o.grid = grid;
  o.apodization = apodization;
  o.focus_shift = focus_shift;
  o.aiming = aiming;
  o.lost_warning_fraction = lost_warning_fraction;
  o.mapping_warning = mapping_warning;
  return o;
}

using Wavelengths = std::optional<std::vector<std::uint16_t>>;

}  // namespace

void bind_psf(nb::module_& m) {
  auto wavelength_columns = columns<PsfWavelengths>(
      m, "PsfWavelengths",
      "The wavelengths of a PSF, MTF or encircled energy, one entry each; read-only NumPy "
      "copies.");
  column<std::uint16_t>(
      wavelength_columns, "wavelength", [](const PsfWavelength& w) { return w.wavelength; },
      "Index into the system wavelengths (copy).");
  column<double>(
      wavelength_columns, "weight", [](const PsfWavelength& w) { return w.weight; },
      "Model weight, normalised to sum 1 over the chosen wavelengths (copy).");
  column<double>(
      wavelength_columns, "energy", [](const PsfWavelength& w) { return w.energy; },
      "E = sum f A^2 of its pupil function (copy).");
  column<std::int32_t>(
      wavelength_columns, "length", [](const PsfWavelength& w) { return w.length; },
      "FFT length L_i (PSF) or 2 N (MTF, encircled energy) (copy).");
  column<double>(
      wavelength_columns, "outside_field", [](const PsfWavelength& w) { return w.outside_field; },
      "PSF: fraction of E outside the reference field (> 0 only for a longer wavelength); 0 "
      "otherwise (copy).");

  nb::class_<PsfPeak>(m, "PsfPeak", "Peak of the PSF lattice, not interpolated.")
      .def_ro("value", &PsfPeak::value, "Intensity of the cell (fraction of the energy).")
      .def_ro("k1", &PsfPeak::k1, "Lattice index along the first lattice vector.")
      .def_ro("k2", &PsfPeak::k2, "Lattice index along the second lattice vector.")
      .def_ro("x", &PsfPeak::x, "Position along e_x, mm.")
      .def_ro("y", &PsfPeak::y, "Position along e_y, mm.");

  nb::class_<Psf>(
      m, "Psf",
      "Diffraction PSF of one field (ADR 0033, points 2-4, 7): intensities on the lattice "
      "x(k) = lattice @ k (mm in e_x, e_y of the plane through C' normal to the chief ray), "
      "fraction of the pupil energy per cell.")
      .def_ro("field", &Psf::field)
      .def_ro("grid", &Psf::grid, "N.")
      .def_ro("length", &Psf::length, "L of the reference lattice.")
      .def_prop_ro(
          "centre", [](const Psf& p) { return copy(p.centre); }, nb::rv_policy::reference,
          "C', global, mm (copy).")
      .def_prop_ro(
          "e_x", [](const Psf& p) { return copy(p.e_x); }, nb::rv_policy::reference,
          "First axis of the PSF plane, global (copy).")
      .def_prop_ro(
          "e_y", [](const Psf& p) { return copy(p.e_y); }, nb::rv_policy::reference,
          "Second axis of the PSF plane, global (copy).")
      .def_prop_ro(
          "lattice", [](const Psf& p) { return matrix2(p.lattice); }, nb::rv_policy::reference,
          "B = lambda / (L a) M^-T (2 x 2), mm per lattice step: x = B[0, 0] k1 + B[0, 1] k2, "
          "y = B[1, 0] k1 + B[1, 1] k2 (copy).")
      .def_prop_ro(
          "intensity",
          [](const Psf& p) { return square(p.intensity, static_cast<std::size_t>(p.length)); },
          nb::rv_policy::reference,
          "(L, L) intensities, row k2 + L // 2, column k1 + L // 2 (k from -L // 2); fraction of "
          "the energy per cell (copy).")
      .def_ro("strehl", &Psf::strehl,
              "Strehl ratio at C' (polychromatic: the peak ratio of the summed images).")
      .def_ro("peak", &Psf::peak, "Largest lattice value.")
      .def_ro("edge_energy_deficit", &Psf::edge_energy_deficit,
              "1 - sum f^2 A^2 / sum f A^2, weighted.")
      .def_ro("polychromatic_scale_error", &Psf::polychromatic_scale_error,
              "max ||B_i - B_ref||_F / ||B_ref||_F over the wavelengths.")
      .def_prop_ro(
          "wavelengths", [](const Psf& p) { return PsfWavelengths{p.wavelengths}; },
          "The chosen wavelengths.")
      .def_ro("losses", &Psf::losses, "Traced cells by final status (ADR 0023).")
      .def_ro("warnings", &Psf::warnings,
              "rays.lost, stop.clips_beam, pupil.mapping_nonlinear, psf.aliasing, "
              "psf.pupil_undersampled; raytatouille.analysis also issues them as "
              "RaytatouilleWarning.");

  // Psf.positions(): pure Python in raytatouille.analysis.psf_positions (coordinator, P1b plan
  // M2); this method only forwards, so that it appears in the stub. Not part of the format.
  const nb::object psf_cls = m.attr("Psf");
  psf_cls.attr("positions") = nb::cpp_function(
      [](nb::handle self) {
        return nb::module_::import_("raytatouille.analysis").attr("psf_positions")(self);
      },
      nb::is_method(), nb::scope(psf_cls), nb::name("positions"),
      nb::sig("def positions(self) -> numpy.typing.NDArray[numpy.float64]"),
      "Positions of the lattice cells, shape (L, L, 2) in mm along e_x, e_y relative to C', in "
      "the order of intensity (raytatouille.analysis.psf_positions).");

  nb::class_<ThroughFocus>(m, "ThroughFocus",
                           "Strehl ratio through focus at one wavelength (ADR 0033, points 8, 11).")
      .def_ro("field", &ThroughFocus::field)
      .def_ro("wavelength", &ThroughFocus::wavelength, "Index into the system wavelengths.")
      .def_prop_ro(
          "shifts", [](const ThroughFocus& t) { return copy(t.shifts); }, nb::rv_policy::reference,
          "Focus shifts, mm along the chief ray (copy).")
      .def_prop_ro(
          "strehl", [](const ThroughFocus& t) { return copy(t.strehl); }, nb::rv_policy::reference,
          "|sum u|^2 / |sum f A|^2 per shift (copy).")
      .def_prop_ro(
          "arrived_fraction", [](const ThroughFocus& t) { return copy(t.arrived_fraction); },
          nb::rv_policy::reference,
          "Area of the cells that cross the shifted reference sphere, per shift (copy).")
      .def_ro("losses", &ThroughFocus::losses, "Traced cells by final status (focus shift 0).")
      .def_ro("warnings", &ThroughFocus::warnings,
              "rays.lost, stop.clips_beam, pupil.mapping_nonlinear; raytatouille.analysis also "
              "issues them as RaytatouilleWarning.");

  nb::class_<MtfCurve>(m, "MtfCurve", "One MTF curve along a pupil direction (ADR 0033, 5).")
      .def_prop_ro(
          "direction",
          [](const MtfCurve& c) {
            return owned<double, 1>(std::vector<double>(c.direction.begin(), c.direction.end()),
                                    {2});
          },
          nb::rv_policy::reference, "Unit direction of nu in (e_x, e_y) (copy).")
      .def_prop_ro(
          "frequency", [](const MtfCurve& c) { return copy(c.frequency); },
          nb::rv_policy::reference,
          "|nu|, cycles per mm, increasing; the last value lies at or above every cut-off (copy).")
      .def_prop_ro(
          "modulation", [](const MtfCurve& c) { return copy(c.modulation); },
          nb::rv_policy::reference, "|OTF|, 1 at nu = 0, 0 at the last frequency (copy).")
      .def_prop_ro(
          "otf", [](const MtfCurve& c) { return copy(c.otf); }, nb::rv_policy::reference,
          "Complex OTF = integral I exp(-2 pi i nu . x) / integral I with complex_otf; empty "
          "otherwise (copy).");

  nb::class_<Mtf>(m, "Mtf", "Diffraction MTF of one field (ADR 0033, point 5).")
      .def_ro("field", &Mtf::field)
      .def_ro("grid", &Mtf::grid, "N.")
      .def_prop_ro(
          "centre", [](const Mtf& t) { return copy(t.centre); }, nb::rv_policy::reference,
          "C', global, mm (copy).")
      .def_prop_ro(
          "e_x", [](const Mtf& t) { return copy(t.e_x); }, nb::rv_policy::reference,
          "First axis of the PSF plane, global (copy).")
      .def_prop_ro(
          "e_y", [](const Mtf& t) { return copy(t.e_y); }, nb::rv_policy::reference,
          "Second axis of the PSF plane, global (copy).")
      .def_ro("tangential", &Mtf::tangential, "Shifts in py (meridional for a field on +-y).")
      .def_ro("sagittal", &Mtf::sagittal, "Shifts in px.")
      .def_prop_ro(
          "frequency_matrix", [](const Mtf& t) { return matrix2(t.frequency_matrix); },
          nb::rv_policy::reference,
          "F = M a / lambda (2 x 2), cycles per mm per cell: nu = F @ l (copy).")
      .def_prop_ro(
          "map",
          [](const Mtf& t) { return square(t.map, static_cast<std::size_t>(2 * t.grid - 1)); },
          nb::rv_policy::reference,
          "With map: (2N - 1, 2N - 1) values of |OTF|, row l2 + N - 1, column l1 + N - 1; "
          "(0, 0) otherwise (copy).")
      .def_prop_ro(
          "otf_map",
          [](const Mtf& t) { return square(t.otf_map, static_cast<std::size_t>(2 * t.grid - 1)); },
          nb::rv_policy::reference,
          "With map and complex_otf: the OTF in the same order; (0, 0) otherwise (copy).")
      .def_ro("edge_energy_deficit", &Mtf::edge_energy_deficit,
              "1 - sum f^2 A^2 / sum f A^2, weighted.")
      .def_prop_ro(
          "wavelengths", [](const Mtf& t) { return PsfWavelengths{t.wavelengths}; },
          "The chosen wavelengths.")
      .def_ro("losses", &Mtf::losses, "Traced cells by final status (ADR 0023).")
      .def_ro("warnings", &Mtf::warnings,
              "rays.lost, stop.clips_beam, pupil.mapping_nonlinear, psf.pupil_undersampled; "
              "raytatouille.analysis also issues them as RaytatouilleWarning.");

  nb::enum_<EeCentre>(m, "EeCentre", "Centre of the encircled-energy circles.")
      .value("REFERENCE", EeCentre::Reference, "C', the centre of the reference sphere")
      .value("LINEAR_PHASE", EeCentre::LinearPhase,
             "the point of the phase of the first harmonic (exact for a tilt; not the centroid)");

  nb::class_<EncircledEnergy>(m, "EncircledEnergy",
                              "Encircled energy of one field (ADR 0033, point 6).")
      .def_ro("field", &EncircledEnergy::field)
      .def_ro("grid", &EncircledEnergy::grid, "N.")
      .def_prop_ro(
          "centre", [](const EncircledEnergy& e) { return copy(e.centre); },
          nb::rv_policy::reference, "Centre of the circles, global, mm (copy).")
      .def_prop_ro(
          "e_x", [](const EncircledEnergy& e) { return copy(e.e_x); }, nb::rv_policy::reference,
          "First axis of the PSF plane, global (copy).")
      .def_prop_ro(
          "e_y", [](const EncircledEnergy& e) { return copy(e.e_y); }, nb::rv_policy::reference,
          "Second axis of the PSF plane, global (copy).")
      .def_ro("offset", &EncircledEnergy::offset,
              "Centre relative to C' in the PSF plane, mm; 0 for REFERENCE.")
      .def_prop_ro(
          "radii", [](const EncircledEnergy& e) { return copy(e.radii); }, nb::rv_policy::reference,
          "Radii, mm, as given (copy).")
      .def_prop_ro(
          "energy", [](const EncircledEnergy& e) { return copy(e.energy); },
          nb::rv_policy::reference, "EE(r), fraction of the pupil energy (copy).")
      .def_ro("radius_max", &EncircledEnergy::radius_max,
              "Largest allowed radius, mm: a quarter of the shortest period.")
      .def_ro("edge_energy_deficit", &EncircledEnergy::edge_energy_deficit,
              "1 - sum f^2 A^2 / sum f A^2, weighted.")
      .def_prop_ro(
          "wavelengths", [](const EncircledEnergy& e) { return PsfWavelengths{e.wavelengths}; },
          "The chosen wavelengths.")
      .def_ro("losses", &EncircledEnergy::losses, "Traced cells by final status (ADR 0023).")
      .def_ro("warnings", &EncircledEnergy::warnings,
              "rays.lost, stop.clips_beam, pupil.mapping_nonlinear, psf.pupil_undersampled; "
              "raytatouille.analysis also issues them as RaytatouilleWarning.");

  // Functions; raytatouille.analysis wraps them. The control is built with the GIL held and
  // only the computation runs without it (released(), #83).
  m.def(
      "psf",
      [](const compile::CompiledSystem& s, const PathArg& path, std::uint16_t field,
         const Wavelengths& wavelengths, int grid, int padding, double apodization,
         double focus_shift, trace::Aiming aiming, double aliasing_factor,
         double undersampled_waves, double lost_warning_fraction, double mapping_warning,
         std::optional<int> threads, const std::optional<trace::CancelToken>& cancel,
         const std::optional<nb::callable>& progress) {
        const trace::RunControl control = run_control(cancel, progress);
        analysis::PsfOptions o;
        o.pupil = pupil_options(grid, apodization, focus_shift, aiming, lost_warning_fraction,
                                mapping_warning);
        o.padding = padding;
        o.wavelengths = wavelengths;
        o.aliasing_factor = aliasing_factor;
        o.undersampled_waves = undersampled_waves;
        const compile::PathId id = path_id(s, path);
        return released(threads, [&] { return analysis::psf(s, id, field, o, control); });
      },
      "system"_a, "path"_a, "field"_a, "wavelengths"_a.none(), "grid"_a, "padding"_a,
      "apodization"_a, "focus_shift"_a, "aiming"_a, "aliasing_factor"_a, "undersampled_waves"_a,
      "lost_warning_fraction"_a, "mapping_warning"_a, "threads"_a.none(), "cancel"_a.none(),
      "progress"_a.none(),
      "Diffraction PSF of `field`; wavelengths None means all with their model weights.");
  m.def(
      "psf_at",
      [](const PupilFunction& pupil,
         const nb::ndarray<const double, nb::shape<-1, 2>, nb::c_contig, nb::device::cpu>& points) {
        std::vector<analysis::Point2> x(points.shape(0));
        const double* v = points.data();
        for (std::size_t i = 0; i < x.size(); ++i) x[i] = {v[2 * i], v[2 * i + 1]};
        std::vector<double> values;
        {
          const nb::gil_scoped_release release;
          values = analysis::psf_at(pupil, x);
        }
        return copy(values);
      },
      "pupil"_a, "points"_a, nb::rv_policy::reference,
      "Intensity in Strehl units at points (n, 2) of the PSF plane, mm in e_x, e_y (copy).");
  m.def(
      "through_focus",
      [](const compile::CompiledSystem& s, const PathArg& path, std::uint16_t field,
         std::optional<std::uint16_t> wavelength, const std::vector<double>& shifts, int grid,
         double apodization, trace::Aiming aiming, double lost_warning_fraction,
         double mapping_warning, std::optional<int> threads,
         const std::optional<trace::CancelToken>& cancel,
         const std::optional<nb::callable>& progress) {
        const trace::RunControl control = run_control(cancel, progress);
        const analysis::PupilFunctionOptions o =
            pupil_options(grid, apodization, 0.0, aiming, lost_warning_fraction, mapping_warning);
        const compile::PathId id = path_id(s, path);
        const std::uint16_t wl = wavelength_index(s, wavelength);
        return released(
            threads, [&] { return analysis::through_focus(s, id, field, wl, shifts, o, control); });
      },
      "system"_a, "path"_a, "field"_a, "wavelength"_a.none(), "shifts"_a, "grid"_a, "apodization"_a,
      "aiming"_a, "lost_warning_fraction"_a, "mapping_warning"_a, "threads"_a.none(),
      "cancel"_a.none(), "progress"_a.none(),
      "Strehl ratio per focus shift; wavelength None means the reference.");
  m.def(
      "mtf",
      [](const compile::CompiledSystem& s, const PathArg& path, std::uint16_t field,
         const Wavelengths& wavelengths, int grid, double apodization, double focus_shift,
         trace::Aiming aiming, bool map, bool complex_otf, double undersampled_waves,
         double lost_warning_fraction, double mapping_warning, std::optional<int> threads,
         const std::optional<trace::CancelToken>& cancel,
         const std::optional<nb::callable>& progress) {
        const trace::RunControl control = run_control(cancel, progress);
        analysis::MtfOptions o;
        o.pupil = pupil_options(grid, apodization, focus_shift, aiming, lost_warning_fraction,
                                mapping_warning);
        o.wavelengths = wavelengths;
        o.outputs.map = map;
        o.outputs.complex_otf = complex_otf;
        o.undersampled_waves = undersampled_waves;
        const compile::PathId id = path_id(s, path);
        return released(threads, [&] { return analysis::mtf(s, id, field, o, control); });
      },
      "system"_a, "path"_a, "field"_a, "wavelengths"_a.none(), "grid"_a, "apodization"_a,
      "focus_shift"_a, "aiming"_a, "map"_a, "complex_otf"_a, "undersampled_waves"_a,
      "lost_warning_fraction"_a, "mapping_warning"_a, "threads"_a.none(), "cancel"_a.none(),
      "progress"_a.none(),
      "Diffraction MTF of `field`; wavelengths None means all with their model weights.");
  m.def(
      "mtf_of_pupil",
      [](const PupilFunction& pupil, bool map, bool complex_otf) {
        analysis::MtfOutputs outputs;
        outputs.map = map;
        outputs.complex_otf = complex_otf;
        return analysis::mtf(pupil, outputs);
      },
      "pupil"_a, "map"_a, "complex_otf"_a, nb::call_guard<nb::gil_scoped_release>(),
      "MTF of one pupil function (monochromatic).");
  m.def(
      "encircled_energy",
      [](const compile::CompiledSystem& s, const PathArg& path, std::uint16_t field,
         const std::vector<double>& radii, EeCentre centre, const Wavelengths& wavelengths,
         int grid, double apodization, double focus_shift, trace::Aiming aiming,
         double undersampled_waves, double lost_warning_fraction, double mapping_warning,
         std::optional<int> threads, const std::optional<trace::CancelToken>& cancel,
         const std::optional<nb::callable>& progress) {
        const trace::RunControl control = run_control(cancel, progress);
        analysis::EncircledEnergyOptions o;
        o.pupil = pupil_options(grid, apodization, focus_shift, aiming, lost_warning_fraction,
                                mapping_warning);
        o.wavelengths = wavelengths;
        o.centre = centre;
        o.undersampled_waves = undersampled_waves;
        const compile::PathId id = path_id(s, path);
        return released(
            threads, [&] { return analysis::encircled_energy(s, id, field, radii, o, control); });
      },
      "system"_a, "path"_a, "field"_a, "radii"_a, "centre"_a, "wavelengths"_a.none(), "grid"_a,
      "apodization"_a, "focus_shift"_a, "aiming"_a, "undersampled_waves"_a,
      "lost_warning_fraction"_a, "mapping_warning"_a, "threads"_a.none(), "cancel"_a.none(),
      "progress"_a.none(),
      "Encircled energy of `field` at `radii` (mm); wavelengths None means all.");
  m.def(
      "encircled_energy_of_pupil",
      [](const PupilFunction& pupil, const std::vector<double>& radii, EeCentre centre) {
        return analysis::encircled_energy(pupil, radii, centre);
      },
      "pupil"_a, "radii"_a, "centre"_a, nb::call_guard<nb::gil_scoped_release>(),
      "Encircled energy of one pupil function (monochromatic) at `radii` (mm).");
}

}  // namespace rtt::py
