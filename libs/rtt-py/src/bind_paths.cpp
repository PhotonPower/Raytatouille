// Path evaluation (#122) and ghost ranking (#124) from Python (#133): the result classes and the
// functions on a CompiledSystem or GhostSystem. raytatouille.analysis wraps the functions
// (System or CompiledSystem, sampling shorthand, warnings).

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "bindings.hpp"
#include "rtt/analysis/ghosts.hpp"
#include "rtt/analysis/paths.hpp"
#include "rtt/compile/ghosts.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sources.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

// Lists inside the results, bound with columns().
struct PathRays {
  std::vector<analysis::PathRay> points;
};
struct OplDifferencePoints {
  std::vector<analysis::OplDifferencePoint> points;
};
struct GhostEntries {
  std::vector<analysis::GhostEntry> points;
};

/// Path options of the convenience form (sampling, aiming) and of both forms (warning
/// threshold).
analysis::PathOptions path_options(const trace::PupilSampling& sampling,
                                   trace::Aiming aiming,
                                   double lost_warning_fraction) {
  analysis::PathOptions options;
  options.sampling = sampling;
  options.aiming = aiming;
  options.lost_warning_fraction = lost_warning_fraction;
  return options;
}

}  // namespace

void bind_paths(nb::module_& m) {
  using namespace analysis;

  auto path_rays = columns<PathRays>(
      m, "PathRays",
      "One entry per launched ray of a path transmission, in batch order; read-only NumPy "
      "copies.");
  column<double>(
      path_rays, "px", [](const PathRay& r) { return r.px; },
      "Normalised pupil x of the start ray (RayBatch.pupil_x) (copy).");
  column<double>(
      path_rays, "py", [](const PathRay& r) { return r.py; },
      "Normalised pupil y of the start ray (RayBatch.pupil_y) (copy).");
  column<double>(
      path_rays, "weight", [](const PathRay& r) { return r.weight; },
      "Final weight if the ray arrived, otherwise 0; dimensionless (copy).");
  column<std::uint8_t>(
      path_rays, "status", [](const PathRay& r) { return status_value(r.status); },
      "RayStatus values as uint8: ALIVE if the ray arrived, otherwise its final status "
      "(VIGNETTED if it stopped ALIVE on another surface) (copy).");

  nb::class_<PathTransmission>(
      m, "PathTransmission",
      "Transmission of one path for a set of start rays (#122, rtt/analysis/paths.hpp): weight "
      "is the power for an unpolarized source (ADR 0021).")
      .def_prop_ro(
          "path", [](const PathTransmission& t) { return t.path.index; }, "The evaluated path.")
      .def_ro("image_surface", &PathTransmission::image_surface,
              "Surface of the path's last event, index into surface_ids.")
      .def_prop_ro(
          "rays", [](const PathTransmission& t) { return PathRays{t.rays}; },
          "One entry per launched ray, in batch order.")
      .def_ro("rays_launched", &PathTransmission::rays_launched, "Rays of the start batch.")
      .def_ro("rays_arrived", &PathTransmission::rays_arrived,
              "Rays that reached the image surface ALIVE.")
      .def_ro("mean", &PathTransmission::mean,
              "Mean final weight over all launched rays, lost rays as 0: the transmitted power "
              "fraction of a uniformly illuminated pupil; with start weights other than 1 the "
              "apodized transmitted power. Dimensionless.")
      .def_ro("min", &PathTransmission::min,
              "Smallest final weight of the arrived rays (0 if none), dimensionless.")
      .def_ro("max", &PathTransmission::max,
              "Largest final weight of the arrived rays (0 if none), dimensionless.")
      .def_ro("losses", &PathTransmission::losses, "Launched rays by final status (ADR 0023).")
      .def_ro("warnings", &PathTransmission::warnings,
              "Warnings with codes: rays.lost above lost_warning_fraction, stop.clips_beam "
              "(ADR 0023); raytatouille.analysis also issues them as RaytatouilleWarning.");

  auto opl_points = columns<OplDifferencePoints>(
      m, "OplDifferencePoints",
      "OPL on path b minus OPL on path a per launched ray, in batch order; read-only NumPy "
      "copies. Check status before using delta (0 for rays that did not arrive on both "
      "paths).");
  column<double>(
      opl_points, "px", [](const OplDifferencePoint& p) { return p.px; },
      "Normalised pupil x of the start ray (copy).");
  column<double>(
      opl_points, "py", [](const OplDifferencePoint& p) { return p.py; },
      "Normalised pupil y of the start ray (copy).");
  column<double>(
      opl_points, "delta", [](const OplDifferencePoint& p) { return p.delta; },
      "OPL_b - OPL_a at the image surface, mm (copy).");
  column<std::uint8_t>(
      opl_points, "status", [](const OplDifferencePoint& p) { return status_value(p.status); },
      "RayStatus values as uint8: ALIVE if the ray arrived on both paths, otherwise its status "
      "on path a if lost there, else on path b (copy).");

  nb::class_<PathOplDifference>(
      m, "PathOplDifference",
      "Optical path difference OPL_b - OPL_a of two paths with the same image surface, per "
      "start ray, in mm (a path length, not a wavefront; #122).")
      .def_prop_ro(
          "path_a", [](const PathOplDifference& d) { return d.path_a.index; },
          "Reference path (subtracted).")
      .def_prop_ro(
          "path_b", [](const PathOplDifference& d) { return d.path_b.index; },
          "Path whose OPL is taken positive.")
      .def_ro("image_surface", &PathOplDifference::image_surface,
              "Common surface of the last events, index into surface_ids.")
      .def_prop_ro(
          "points", [](const PathOplDifference& d) { return OplDifferencePoints{d.points}; },
          "One entry per launched ray, in batch order.")
      .def_ro("chief", &PathOplDifference::chief,
              "delta of the first launched ray at pupil (0, 0) that arrived on both paths, mm; "
              "None if there is no such ray (it needs the pupil labels of the start rays).")
      .def_ro("losses_a", &PathOplDifference::losses_a, "Launched rays on path a by final status.")
      .def_ro("losses_b", &PathOplDifference::losses_b, "Launched rays on path b by final status.")
      .def_ro("warnings", &PathOplDifference::warnings,
              "Warnings of both paths, path a first (ADR 0023); raytatouille.analysis also "
              "issues them as RaytatouilleWarning.");

  auto ghost_entries = columns<GhostEntries>(
      m, "GhostEntries",
      "Ghosts of a ranking, by relative_irradiance descending (ties in ghost order); read-only "
      "NumPy copies, NaN where a diagnostic has no value.");
  column<std::uint32_t>(
      ghost_entries, "path", [](const GhostEntry& e) { return e.path.index; },
      "Index of the ghost path in GhostSystem.system (copy).");
  column<std::uint32_t>(
      ghost_entries, "surface_j", [](const GhostEntry& e) { return e.surface_j; },
      "First ghost reflection (back), index into surface_ids (copy).");
  column<std::uint32_t>(
      ghost_entries, "surface_i", [](const GhostEntry& e) { return e.surface_i; },
      "Second ghost reflection (forward), index into surface_ids (copy).");
  column<double>(
      ghost_entries, "power", [](const GhostEntry& e) { return e.power; },
      "P_g: mean final weight over the launched rays, dimensionless (copy).");
  column<double>(
      ghost_entries, "relative_power", [](const GhostEntry& e) { return e.relative_power; },
      "P_g / P_b, dimensionless (copy).");
  column<double>(
      ghost_entries, "rms_radius", [](const GhostEntry& e) { return e.rms_radius; },
      "r_g: weighted RMS radius on the image surface, mm; 0 without arrived rays (copy).");
  column<std::uint64_t>(
      ghost_entries, "rays_arrived", [](const GhostEntry& e) { return e.rays_arrived; },
      "Rays of the ghost that reached the image surface (copy).");
  column<double>(
      ghost_entries, "relative_irradiance",
      [](const GhostEntry& e) { return e.relative_irradiance; },
      "rho = (P_g / P_b) (r_b^2 + r0^2) / (r_g^2 + r0^2): the rank value, dimensionless "
      "(copy).");
  column<double>(
      ghost_entries, "focus_offset",
      [](const GhostEntry& e) { return e.focus_offset.value_or(kNan); },
      "Diagnostic: global z of the paraxial ghost focus minus the global z of the image surface "
      "vertex, mm, positive towards +z; NaN if the ghost leaves collimated or the base path has "
      "no finite entrance pupil (copy).");
  column<double>(
      ghost_entries, "paraxial_blur_radius",
      [](const GhostEntry& e) { return e.paraxial_blur_radius.value_or(kNan); },
      "Diagnostic: |height| of the same paraxial ray at the image surface vertex plane, mm; NaN "
      "without a finite entrance pupil (copy).");
  ghost_entries.def_prop_ro(
      "losses",
      [](const GhostEntries& g) {
        std::vector<RayLosses> losses;
        losses.reserve(g.points.size());
        for (const GhostEntry& e : g.points) losses.push_back(e.losses);
        return losses;
      },
      "Launched rays of each ghost by final status (ADR 0023), one RayLosses per entry (copy).");

  nb::class_<GhostRanking>(
      m, "GhostRanking",
      "Ghosts of a base path ranked by relative irradiance at the image (#124, "
      "rtt/analysis/ghosts.hpp): rho = (P_g / P_b) (r_b^2 + r0^2) / (r_g^2 + r0^2). rho depends "
      "on the resolution radius r0 (a detector model); the RMS radius is geometric.")
      .def_prop_ro(
          "base", [](const GhostRanking& r) { return r.base.index; },
          "The base path (useful image).")
      .def_ro("field", &GhostRanking::field, "Field index of the start rays.")
      .def_ro("wavelength", &GhostRanking::wavelength, "Wavelength index of the start rays.")
      .def_ro("base_power", &GhostRanking::base_power, "P_b, dimensionless.")
      .def_ro("base_rms_radius", &GhostRanking::base_rms_radius, "r_b, mm.")
      .def_ro("resolution_radius", &GhostRanking::resolution_radius, "r0 used, mm.")
      .def_prop_ro(
          "entries", [](const GhostRanking& r) { return GhostEntries{r.entries}; },
          "Ghosts by relative_irradiance descending, ties in ghost order.")
      .def_ro("warnings", &GhostRanking::warnings,
              "Warnings of the base path only (ADR 0023); raytatouille.analysis also issues "
              "them as RaytatouilleWarning.");

  // Functions; raytatouille.analysis wraps them. The control is built with the GIL held and
  // only the computation runs without it (released(), #83).
  m.def(
      "path_transmission_rays",
      [](const compile::CompiledSystem& s, const PathArg& path, const trace::RayBatch& rays,
         double lost_warning_fraction, std::optional<int> threads,
         const std::optional<trace::CancelToken>& cancel,
         const std::optional<nb::callable>& progress) {
        const trace::RunControl control = run_control(cancel, progress);
        PathOptions options;
        options.lost_warning_fraction = lost_warning_fraction;
        const compile::PathId id = path_id(s, path);
        return released(threads, [&] { return path_transmission(s, id, rays, options, control); });
      },
      "system"_a, "path"_a, "rays"_a, "lost_warning_fraction"_a, "threads"_a.none(),
      "cancel"_a.none(), "progress"_a.none(),
      "Transmission of `path` for the start rays `rays` (copied; main form).");
  m.def(
      "path_transmission",
      [](const compile::CompiledSystem& s, const PathArg& path, std::uint16_t field,
         std::optional<std::uint16_t> wavelength, const trace::PupilSampling& sampling,
         trace::Aiming aiming, double lost_warning_fraction, std::optional<int> threads,
         const std::optional<trace::CancelToken>& cancel,
         const std::optional<nb::callable>& progress) {
        const trace::RunControl control = run_control(cancel, progress);
        const PathOptions options = path_options(sampling, aiming, lost_warning_fraction);
        const compile::PathId id = path_id(s, path);
        const std::uint16_t wl = wavelength_index(s, wavelength);
        return released(threads,
                        [&] { return path_transmission(s, id, field, wl, options, control); });
      },
      "system"_a, "path"_a, "field"_a, "wavelength"_a.none(), "sampling"_a, "aiming"_a,
      "lost_warning_fraction"_a, "threads"_a.none(), "cancel"_a.none(), "progress"_a.none(),
      "Transmission of `path` for the rays of make_rays (convenience form); wavelength None "
      "means the reference.");
  m.def(
      "opl_difference_rays",
      [](const compile::CompiledSystem& s, const PathArg& path_a, const PathArg& path_b,
         const trace::RayBatch& rays, double lost_warning_fraction, std::optional<int> threads,
         const std::optional<trace::CancelToken>& cancel,
         const std::optional<nb::callable>& progress) {
        const trace::RunControl control = run_control(cancel, progress);
        PathOptions options;
        options.lost_warning_fraction = lost_warning_fraction;
        const compile::PathId a = path_id(s, path_a);
        const compile::PathId b = path_id(s, path_b);
        return released(threads, [&] { return opl_difference(s, a, b, rays, options, control); });
      },
      "system"_a, "path_a"_a, "path_b"_a, "rays"_a, "lost_warning_fraction"_a, "threads"_a.none(),
      "cancel"_a.none(), "progress"_a.none(),
      "OPL_b - OPL_a for the start rays `rays` (copied; main form).");
  m.def(
      "opl_difference",
      [](const compile::CompiledSystem& s, const PathArg& path_a, const PathArg& path_b,
         std::uint16_t field, std::optional<std::uint16_t> wavelength,
         const trace::PupilSampling& sampling, trace::Aiming aiming, double lost_warning_fraction,
         std::optional<int> threads, const std::optional<trace::CancelToken>& cancel,
         const std::optional<nb::callable>& progress) {
        const trace::RunControl control = run_control(cancel, progress);
        const PathOptions options = path_options(sampling, aiming, lost_warning_fraction);
        const compile::PathId a = path_id(s, path_a);
        const compile::PathId b = path_id(s, path_b);
        const std::uint16_t wl = wavelength_index(s, wavelength);
        return released(threads,
                        [&] { return opl_difference(s, a, b, field, wl, options, control); });
      },
      "system"_a, "path_a"_a, "path_b"_a, "field"_a, "wavelength"_a.none(), "sampling"_a,
      "aiming"_a, "lost_warning_fraction"_a, "threads"_a.none(), "cancel"_a.none(),
      "progress"_a.none(),
      "OPL_b - OPL_a for the rays of make_rays on path a (convenience form); wavelength None "
      "means the reference.");
  m.def(
      "ghost_ranking",
      [](const compile::GhostSystem& ghosts, std::uint16_t field,
         std::optional<std::uint16_t> wavelength, const trace::PupilSampling& sampling,
         trace::Aiming aiming, double resolution_radius, double lost_warning_fraction,
         std::optional<int> threads, const std::optional<trace::CancelToken>& cancel,
         const std::optional<nb::callable>& progress) {
        const trace::RunControl control = run_control(cancel, progress);
        GhostRankingOptions options;
        options.sampling = sampling;
        options.aiming = aiming;
        options.resolution_radius = resolution_radius;
        options.lost_warning_fraction = lost_warning_fraction;
        const std::uint16_t wl = wavelength_index(ghosts.system, wavelength);
        return released(threads,
                        [&] { return ghost_ranking(ghosts, field, wl, options, control); });
      },
      "ghosts"_a, "field"_a, "wavelength"_a.none(), "sampling"_a, "aiming"_a, "resolution_radius"_a,
      "lost_warning_fraction"_a, "threads"_a.none(), "cancel"_a.none(), "progress"_a.none(),
      "Ghosts of a GhostSystem ranked by relative irradiance; wavelength None means the "
      "reference.");
}

}  // namespace rtt::py
