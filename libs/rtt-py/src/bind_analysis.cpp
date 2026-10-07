#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "bindings.hpp"
#include "rtt/analysis/chromatic.hpp"
#include "rtt/analysis/field.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/model/system.hpp"
#include "rtt/paraxial/prescription.hpp"
#include "rtt/paraxial/seidel.hpp"
#include "rtt/trace/sources.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

// Lists inside analysis results, bound with columns() below.
struct Points2 {
  std::vector<analysis::Point2> points;
};
struct FanPoints {
  std::vector<analysis::FanPoint> points;
};
struct OpdPoints {
  std::vector<analysis::OpdPoint> points;
};
struct Foci {
  std::vector<analysis::FocusPosition> points;
};
struct DistortionSweep {
  std::vector<analysis::DistortionPoint> points;
};
struct FieldCurvatureSweep {
  std::vector<analysis::FieldCurvaturePoint> points;
};
struct SeidelSurfaces {
  std::vector<paraxial::SeidelSurface> points;
};
struct PrescriptionSurfaces {
  std::vector<paraxial::PrescriptionSurface> points;
};

constexpr double kNan = std::numeric_limits<double>::quiet_NaN();  // None in float arrays

/// Binds a list class; its attributes are read-only NumPy copies of the C++ result (ADR 0002:
/// analysis results are small), one array per member of the point type.
template <typename C>
nb::class_<C> columns(nb::module_& m, const char* name, const char* doc) {
  nb::class_<C> cls(m, name, doc);
  cls.def("__len__", [](const C& c) { return c.points.size(); });
  return cls;
}

/// Adds a read-only array attribute `name` with values get(point) of type T.
template <typename T, typename C, typename Get>
void column(nb::class_<C>& cls, const char* name, Get get, const char* doc) {
  // reference: the NumPy array views the copy and keeps its capsule alive (move would copy
  // again into a writable array; reference_internal needs an array without owner).
  cls.def_prop_ro(
      name, [get](const C& c) { return read_only_array<T>(c.points, get); },
      nb::rv_policy::reference, doc);
}

/// Element i with Python indexing (negative from the end).
/// @throws std::out_of_range (IndexError) outside -size <= i < size
template <typename T>
T at(const std::vector<T>& v, std::ptrdiff_t i) {
  const auto size = static_cast<std::ptrdiff_t>(v.size());
  if (i < -size || i >= size) throw std::out_of_range("index out of range");
  return v[static_cast<std::size_t>(i < 0 ? i + size : i)];
}

std::uint8_t status_value(trace::RayStatus s) {
  return static_cast<std::uint8_t>(s);
}

}  // namespace

void bind_analysis(nb::module_& m) {
  using namespace analysis;

  nb::class_<Point2>(m, "Point2", "Point in the local x, y plane of the image surface, mm.")
      .def_ro("x", &Point2::x)
      .def_ro("y", &Point2::y)
      .def("__repr__",
           [](const Point2& p) { return nb::str("Point2(x={!r}, y={!r})").format(p.x, p.y); });

  auto points2 = columns<Points2>(
      m, "Points2", "Points on the image surface (local x, y in mm), e.g. per wavelength.");
  column<double>(points2, "x", [](const Point2& p) { return p.x; }, "Local x, mm (copy).");
  column<double>(points2, "y", [](const Point2& p) { return p.y; }, "Local y, mm (copy).");

  nb::class_<SpotStatistics>(m, "SpotStatistics", "Weighted statistics of a spot diagram, mm.")
      .def_ro("centroid", &SpotStatistics::centroid, "c = sum w r / sum w.")
      .def_ro("rms_centroid", &SpotStatistics::rms_centroid, "RMS radius about the centroid.")
      .def_ro("rms_chief", &SpotStatistics::rms_chief, "RMS radius about the chief ray.")
      .def_ro("geo_centroid", &SpotStatistics::geo_centroid, "Max distance from the centroid.")
      .def_ro("geo_chief", &SpotStatistics::geo_chief, "Max distance from the chief ray.");

  nb::class_<SpotDiagram> spot_diagram(
      m, "SpotDiagram",
      "Spot diagram of one field. Coordinates in mm in the local x, y of the image surface (the "
      "surface of the last path event); the chief ray is the pupil-centre ray of the reference "
      "wavelength. x, y, wavelengths and weight are read-only NumPy copies.");
  spot_diagram.def_ro("field", &SpotDiagram::field, "Field index.")
      .def_ro("wavelength", &SpotDiagram::wavelength,
              "Wavelength index, or None for a polychromatic spot.")
      .def_ro("image_surface", &SpotDiagram::image_surface, "Index into surface_ids.")
      .def_ro("chief", &SpotDiagram::chief, "Chief ray on the image surface, mm.")
      .def_ro("stats", &SpotDiagram::stats)
      .def_ro("rays_launched", &SpotDiagram::rays_launched)
      .def_ro("rays_arrived", &SpotDiagram::rays_arrived)
      .def_ro("vignetted_fraction", &SpotDiagram::vignetted_fraction,
              "(launched - arrived) / launched, unweighted.")
      .def_prop_ro(
          "x",
          [](const SpotDiagram& s) {
            return read_only_array<double>(s.points, [](const SpotPoint& p) { return p.x; });
          },
          nb::rv_policy::reference, "Local x of the arrived rays, mm (copy).")
      .def_prop_ro(
          "y",
          [](const SpotDiagram& s) {
            return read_only_array<double>(s.points, [](const SpotPoint& p) { return p.y; });
          },
          nb::rv_policy::reference, "Local y of the arrived rays, mm (copy).")
      .def_prop_ro(
          "wavelengths",
          [](const SpotDiagram& s) {
            return read_only_array<std::uint16_t>(s.points,
                                                  [](const SpotPoint& p) { return p.wavelength; });
          },
          nb::rv_policy::reference, "Wavelength index of each arrived ray (copy).")
      .def_prop_ro(
          "weight",
          [](const SpotDiagram& s) {
            return read_only_array<double>(s.points, [](const SpotPoint& p) { return p.weight; });
          },
          nb::rv_policy::reference, "Statistical weight of each arrived ray (copy).");

  auto fan_points = columns<FanPoints>(
      m, "FanPoints",
      "One ray fan: transverse aberration relative to the chief ray. Check status before using "
      "ex, ey (0 for rays that did not arrive).");
  column<double>(
      fan_points, "p", [](const FanPoint& q) { return q.p; },
      "Normalised pupil coordinate (py tangential, px sagittal) (copy).");
  column<double>(
      fan_points, "ex", [](const FanPoint& q) { return q.ex; },
      "x - chief.x on the image surface, mm (copy).");
  column<double>(
      fan_points, "ey", [](const FanPoint& q) { return q.ey; },
      "y - chief.y on the image surface, mm (copy).");
  column<std::uint8_t>(
      fan_points, "status", [](const FanPoint& q) { return status_value(q.status); },
      "RayStatus values as uint8 (copy).");

  nb::class_<RayFan>(m, "RayFan", "Tangential and sagittal ray fans of one field.")
      .def_ro("field", &RayFan::field)
      .def_ro("wavelength", &RayFan::wavelength)
      .def_ro("image_surface", &RayFan::image_surface)
      .def_ro("chief", &RayFan::chief)
      .def_prop_ro(
          "tangential", [](const RayFan& f) { return FanPoints{f.tangential}; },
          "px = 0, py in [-1, 1]; ey(py) is the tangential aberration.")
      .def_prop_ro(
          "sagittal", [](const RayFan& f) { return FanPoints{f.sagittal}; },
          "py = 0, px in [-1, 1]; ex(px) is the sagittal aberration.");

  nb::class_<ReferenceSphere>(m, "ReferenceSphere", "Reference sphere in global coordinates.")
      .def_prop_ro(
          "centre",
          [](const ReferenceSphere& s) {
            const std::vector<double> xyz = {s.centre.x(), s.centre.y(), s.centre.z()};
            return read_only_array<double>(xyz, [](double v) { return v; });
          },
          nb::rv_policy::reference, "Chief-ray point on the image surface (x, y, z), mm (copy).")
      .def_ro("radius", &ReferenceSphere::radius, "Distance to the exit-pupil centre, mm.");

  auto opd_points = columns<OpdPoints>(
      m, "OpdPoints",
      "OPD samples in waves at the reference wavelength, W > 0 leading (Wyant & Creath). Check "
      "status before using w (0 for rays that did not arrive).");
  column<double>(
      opd_points, "px", [](const OpdPoint& q) { return q.px; }, "Normalised pupil x (copy).");
  column<double>(
      opd_points, "py", [](const OpdPoint& q) { return q.py; }, "Normalised pupil y (copy).");
  column<double>(opd_points, "w", [](const OpdPoint& q) { return q.w; }, "OPD, waves (copy).");
  column<std::uint8_t>(
      opd_points, "status", [](const OpdPoint& q) { return status_value(q.status); },
      "RayStatus values as uint8 (copy).");

  nb::class_<OpdMap>(m, "OpdMap", "OPD map over the pupil (GridPupil order).")
      .def_ro("field", &OpdMap::field)
      .def_ro("wavelength", &OpdMap::wavelength)
      .def_ro("sphere", &OpdMap::sphere)
      .def_prop_ro(
          "points", [](const OpdMap& o) { return OpdPoints{o.points}; },
          "Grid points (rows from py = -1, px fastest).")
      .def_ro("rms", &OpdMap::rms, "Standard deviation of W (piston removed), waves.")
      .def_ro("pv", &OpdMap::pv, "Max W - min W, waves.")
      .def_ro("arrived", &OpdMap::arrived)
      .def_ro("vignetted", &OpdMap::vignetted);

  nb::class_<OpdFan>(m, "OpdFan", "Tangential (px = 0) and sagittal (py = 0) OPD fans.")
      .def_ro("field", &OpdFan::field)
      .def_ro("wavelength", &OpdFan::wavelength)
      .def_ro("sphere", &OpdFan::sphere)
      .def_prop_ro("tangential", [](const OpdFan& o) { return OpdPoints{o.tangential}; })
      .def_prop_ro("sagittal", [](const OpdFan& o) { return OpdPoints{o.sagittal}; });

  auto foci = columns<Foci>(m, "Foci", "Focus positions per system wavelength.");
  column<std::uint16_t>(
      foci, "wavelength", [](const FocusPosition& f) { return f.wavelength; },
      "Wavelength index (copy).");
  column<double>(
      foci, "paraxial_z", [](const FocusPosition& f) { return f.paraxial_z; },
      "Global z of the paraxial focus, mm (copy).");
  column<double>(
      foci, "real_z", [](const FocusPosition& f) { return f.real_z; },
      "Global z of the axis crossing of the real zone ray, mm (copy).");

  nb::class_<LongitudinalColour>(
      m, "LongitudinalColour",
      "Longitudinal colour: focus(first) - focus(second) along the image-space propagation, mm.")
      .def_ro("pair", &LongitudinalColour::pair, nb::rv_policy::copy, "Pair used (a copy).")
      .def_prop_ro("foci", [](const LongitudinalColour& l) { return Foci{l.foci}; })
      .def_ro("paraxial", &LongitudinalColour::paraxial, "Paraxial difference, mm.")
      .def_ro("real", &LongitudinalColour::real, "Difference of the real zone rays, mm.");

  nb::class_<LateralColour>(m, "LateralColour",
                            "Chief ray per wavelength on the image surface and its offset from "
                            "the reference wavelength, mm.")
      .def_ro("field", &LateralColour::field)
      .def_prop_ro("chief", [](const LateralColour& l) { return Points2{l.chief}; })
      .def_prop_ro("offset", [](const LateralColour& l) { return Points2{l.offset}; });

  nb::class_<DistortionPoint>(m, "DistortionPoint", "Distortion at one field value.")
      .def_ro("fraction", &DistortionPoint::fraction, "Relative field.")
      .def_ro("field", &DistortionPoint::field, nb::rv_policy::copy,
              "Field value in the system's field units (a copy).")
      .def_ro("real_height", &DistortionPoint::real_height, "Signed real chief-ray height, mm.")
      .def_ro("paraxial_height", &DistortionPoint::paraxial_height,
              "Signed paraxial chief-ray height at the image-surface vertex plane, mm.")
      .def_ro("percent", &DistortionPoint::percent, "D = (real - paraxial) / paraxial, %.");
  auto distortion_sweep = columns<DistortionSweep>(
      m, "DistortionSweep",
      "Distortion over the field sweep (relative field 0 ... 1 along +y of the largest field).");
  distortion_sweep.def(
      "__getitem__", [](const DistortionSweep& d, std::ptrdiff_t i) { return at(d.points, i); },
      "i"_a);
  column<double>(
      distortion_sweep, "fraction", [](const DistortionPoint& p) { return p.fraction; },
      "Relative field (copy).");
  column<double>(
      distortion_sweep, "field_x", [](const DistortionPoint& p) { return p.field.x; },
      "Field x in the system's field units (copy).");
  column<double>(
      distortion_sweep, "field_y", [](const DistortionPoint& p) { return p.field.y; },
      "Field y in the system's field units (copy).");
  column<double>(
      distortion_sweep, "real_height", [](const DistortionPoint& p) { return p.real_height; },
      "mm (copy).");
  column<double>(
      distortion_sweep, "paraxial_height",
      [](const DistortionPoint& p) { return p.paraxial_height; }, "mm (copy).");
  column<double>(
      distortion_sweep, "percent", [](const DistortionPoint& p) { return p.percent; }, "% (copy).");

  nb::class_<FieldCurvaturePoint>(
      m, "FieldCurvaturePoint",
      "Tangential and sagittal focus at one field value, from the image-surface vertex along the "
      "image-space propagation, mm.")
      .def_ro("fraction", &FieldCurvaturePoint::fraction)
      .def_ro("field", &FieldCurvaturePoint::field, nb::rv_policy::copy,
              "Field value in the system's field units (a copy).")
      .def_ro("tangential", &FieldCurvaturePoint::tangential)
      .def_ro("sagittal", &FieldCurvaturePoint::sagittal)
      .def_ro("astigmatism", &FieldCurvaturePoint::astigmatism, "tangential - sagittal, mm.");
  auto field_curvature_sweep =
      columns<FieldCurvatureSweep>(m, "FieldCurvatureSweep", "Field curvature over the sweep.");
  field_curvature_sweep.def(
      "__getitem__", [](const FieldCurvatureSweep& d, std::ptrdiff_t i) { return at(d.points, i); },
      "i"_a);
  column<double>(
      field_curvature_sweep, "fraction", [](const FieldCurvaturePoint& p) { return p.fraction; },
      "(copy)");
  column<double>(
      field_curvature_sweep, "field_x", [](const FieldCurvaturePoint& p) { return p.field.x; },
      "(copy)");
  column<double>(
      field_curvature_sweep, "field_y", [](const FieldCurvaturePoint& p) { return p.field.y; },
      "(copy)");
  column<double>(
      field_curvature_sweep, "tangential",
      [](const FieldCurvaturePoint& p) { return p.tangential; }, "mm (copy).");
  column<double>(
      field_curvature_sweep, "sagittal", [](const FieldCurvaturePoint& p) { return p.sagittal; },
      "mm (copy).");
  column<double>(
      field_curvature_sweep, "astigmatism",
      [](const FieldCurvaturePoint& p) { return p.astigmatism; }, "mm (copy).");

  // Seidel sums (rtt-paraxial, seidel.hpp).
  using paraxial::ChromaticPair;
  using paraxial::RayStart;
  using paraxial::SeidelSurface;
  using paraxial::SeidelTerms;
  nb::class_<ChromaticPair>(m, "ChromaticPair",
                            "Wavelength indices for colour terms: value(first) - value(second).")
      .def(
          "__init__",
          [](ChromaticPair* p, std::uint16_t first, std::uint16_t second) {
            new (p) ChromaticPair{first, second};
          },
          "first"_a, "second"_a)
      .def_rw("first", &ChromaticPair::first)
      .def_rw("second", &ChromaticPair::second);
  nb::class_<SeidelTerms>(m, "SeidelTerms",
                          "Seidel sums and chromatic terms in mm (W040 = S_I / 8, W131 = S_II / "
                          "2, W222 = S_III / 2, W220 = (S_IV + S_III) / 4, W311 = S_V / 2).")
      .def_ro("s1", &SeidelTerms::s1, "S_I, spherical aberration.")
      .def_ro("s2", &SeidelTerms::s2, "S_II, coma.")
      .def_ro("s3", &SeidelTerms::s3, "S_III, astigmatism.")
      .def_ro("s4", &SeidelTerms::s4, "S_IV, Petzval field curvature.")
      .def_ro("s5", &SeidelTerms::s5, "S_V, distortion.")
      .def_ro("c_l", &SeidelTerms::c_l, "C_L, longitudinal colour; 0 without a pair.")
      .def_ro("c_t", &SeidelTerms::c_t, "C_T, transverse colour; 0 without a pair.");
  nb::class_<RayStart>(m, "RayStart", "Paraxial start of a ray in object space.")
      .def_ro("z", &RayStart::z, "Global z of the start plane, mm.")
      .def_ro("y", &RayStart::y, "Height, mm.")
      .def_ro("u", &RayStart::u, "Slope dy/dz.");
  auto seidel_surfaces = columns<SeidelSurfaces>(
      m, "SeidelSurfaces", "Contributions per event of the path, in path order, mm.");
  column<std::uint32_t>(
      seidel_surfaces, "surface", [](const SeidelSurface& q) { return q.surface; },
      "Index into surface_ids (copy).");
  column<double>(
      seidel_surfaces, "y", [](const SeidelSurface& q) { return q.y; },
      "Marginal ray height at the vertex plane, mm (copy).");
  column<double>(
      seidel_surfaces, "y_bar", [](const SeidelSurface& q) { return q.y_bar; },
      "Chief ray height at the vertex plane, mm (copy).");
  column<double>(
      seidel_surfaces, "a", [](const SeidelSurface& q) { return q.a; },
      "Refraction invariant A of the marginal ray (copy).");
  column<double>(
      seidel_surfaces, "a_bar", [](const SeidelSurface& q) { return q.a_bar; },
      "Refraction invariant of the chief ray (copy).");
  column<double>(
      seidel_surfaces, "lagrange", [](const SeidelSurface& q) { return q.lagrange; },
      "Lagrange invariant after the event, mm (copy).");
  column<double>(
      seidel_surfaces, "s1", [](const SeidelSurface& q) { return q.terms.s1; }, "(copy)");
  column<double>(
      seidel_surfaces, "s2", [](const SeidelSurface& q) { return q.terms.s2; }, "(copy)");
  column<double>(
      seidel_surfaces, "s3", [](const SeidelSurface& q) { return q.terms.s3; }, "(copy)");
  column<double>(
      seidel_surfaces, "s4", [](const SeidelSurface& q) { return q.terms.s4; }, "(copy)");
  column<double>(
      seidel_surfaces, "s5", [](const SeidelSurface& q) { return q.terms.s5; }, "(copy)");
  column<double>(
      seidel_surfaces, "c_l", [](const SeidelSurface& q) { return q.terms.c_l; }, "(copy)");
  column<double>(
      seidel_surfaces, "c_t", [](const SeidelSurface& q) { return q.terms.c_t; }, "(copy)");
  nb::class_<paraxial::Seidel>(
      m, "Seidel",
      "Seidel sums of a path from the paraxial marginal ray and the chief ray of the maximum "
      "field (conventions in rtt/paraxial/seidel.hpp).")
      .def_prop_ro("surfaces", [](const paraxial::Seidel& s) { return SeidelSurfaces{s.surfaces}; })
      .def_ro("sum", &paraxial::Seidel::sum, "Sum over all events, mm.")
      .def_ro("lagrange", &paraxial::Seidel::lagrange, "H = n (u_bar y - u y_bar), mm.")
      .def_ro("marginal", &paraxial::Seidel::marginal)
      .def_ro("chief", &paraxial::Seidel::chief)
      .def_ro("chromatic", &paraxial::Seidel::chromatic, nb::rv_policy::copy,
              "Pair used for C_L, C_T (a copy), or None.");
  m.def(
      "seidel",
      [](const compile::CompiledSystem& s, const PathArg& path,
         std::optional<std::uint16_t> wavelength, std::optional<ChromaticPair> pair) {
        return paraxial::seidel(s, path_id(s, path), wavelength_index(s, wavelength), pair);
      },
      "system"_a, "path"_a, "wavelength"_a.none(), "pair"_a.none(),
      "Seidel sums; wavelength None means the reference, pair None gives C_L = C_T = 0.");

  using paraxial::Prescription;
  using paraxial::PrescriptionSurface;
  auto prescription_surfaces = columns<PrescriptionSurfaces>(
      m, "PrescriptionSurfaces",
      "Paraxial marginal and chief ray per event of the path, in path order; NaN where a ray "
      "is not defined (conventions in rtt/paraxial/prescription.hpp).");
  column<std::uint32_t>(
      prescription_surfaces, "surface", [](const PrescriptionSurface& q) { return q.surface; },
      "Index into surface_ids (copy).");
  column<double>(
      prescription_surfaces, "z", [](const PrescriptionSurface& q) { return q.z; },
      "Global z of the surface vertex, mm (copy).");
  column<double>(
      prescription_surfaces, "n", [](const PrescriptionSurface& q) { return q.n; },
      "Signed index after the event, negative while light travels towards -z (copy).");
  column<double>(
      prescription_surfaces, "y",
      [](const PrescriptionSurface& q) { return q.marginal ? q.marginal->y : kNan; },
      "Marginal ray height at the vertex plane, mm (copy).");
  column<double>(
      prescription_surfaces, "u",
      [](const PrescriptionSurface& q) { return q.marginal ? q.marginal->u : kNan; },
      "Marginal ray slope dy/dz after the event (copy).");
  column<double>(
      prescription_surfaces, "i",
      [](const PrescriptionSurface& q) { return q.marginal ? q.marginal->i : kNan; },
      "Paraxial angle of incidence of the marginal ray, i = u + y c with u before the "
      "event, rad (copy).");
  column<double>(
      prescription_surfaces, "y_bar",
      [](const PrescriptionSurface& q) { return q.chief ? q.chief->y : kNan; },
      "Chief ray height at the vertex plane, mm (copy).");
  column<double>(
      prescription_surfaces, "u_bar",
      [](const PrescriptionSurface& q) { return q.chief ? q.chief->u : kNan; },
      "Chief ray slope dy/dz after the event (copy).");
  column<double>(
      prescription_surfaces, "i_bar",
      [](const PrescriptionSurface& q) { return q.chief ? q.chief->i : kNan; },
      "Paraxial angle of incidence of the chief ray, rad (copy).");
  column<double>(
      prescription_surfaces, "lagrange",
      [](const PrescriptionSurface& q) { return q.lagrange.value_or(kNan); },
      "Lagrange invariant n' (u_bar' y - u' y_bar) after the event, mm (copy).");
  nb::class_<Prescription>(
      m, "Prescription",
      "Paraxial prescription data of a path: marginal and chief ray per event and system data "
      "(conventions in rtt/paraxial/prescription.hpp). Lateral and angular magnification are "
      "in first_order.")
      .def_prop_ro("surfaces",
                   [](const Prescription& p) { return PrescriptionSurfaces{p.surfaces}; })
      .def_ro("total_track", &Prescription::total_track,
              "Sum of |dz| between consecutive event vertices, mm (unfolded with mirrors).")
      .def_ro("object_distance", &Prescription::object_distance,
              "First event vertex minus object z, mm; None for an object at infinity.")
      .def_ro("paraxial_working_f_number", &Prescription::paraxial_working_f_number,
              "1 / (2 |n' u'|) of the marginal ray in image space (any cone of light, also "
              "afocal with a finite object); None if u' = 0 (afocal, object at infinity) or "
              "without a marginal ray.")
      .def_ro("paraxial_image_na", &Prescription::paraxial_image_na,
              "|n' u'| of the marginal ray in image space; None without a marginal ray.")
      .def_ro("lagrange_invariant", &Prescription::lagrange_invariant,
              "H = n (u_bar y - u y_bar) in object space, mm; None without both rays.")
      .def_ro("marginal_start", &Prescription::marginal_start, nb::rv_policy::copy,
              "Marginal ray in object space (a copy), or None.")
      .def_ro("chief_start", &Prescription::chief_start, nb::rv_policy::copy,
              "Chief ray in object space (a copy), or None.")
      .def_ro("first_order", &Prescription::first_order, "First-order data of the path.");
  m.def(
      "prescription",
      [](const compile::CompiledSystem& s, const PathArg& path,
         std::optional<std::uint16_t> wavelength) {
        return paraxial::prescription(s, path_id(s, path), wavelength_index(s, wavelength));
      },
      "system"_a, "path"_a, "wavelength"_a.none(),
      "Paraxial prescription data; wavelength None means the reference.");

  // Analysis functions on a CompiledSystem; raytatouille.analysis wraps them (System or
  // CompiledSystem, sampling shorthand). The GIL is released, `threads` limits the workers.
  const auto release = nb::call_guard<nb::gil_scoped_release>();
  m.def(
      "spot",
      [](const compile::CompiledSystem& s, const PathArg& path, std::uint16_t field,
         std::optional<std::uint16_t> wavelength, const trace::PupilSampling& sampling,
         trace::Aiming aiming, std::optional<int> threads) {
        const SpotOptions options{sampling, aiming};
        return with_threads(threads, [&] {
          return analysis::spot(s, path_id(s, path), field, wavelength, options);
        });
      },
      "system"_a, "path"_a, "field"_a, "wavelength"_a.none(), "sampling"_a, "aiming"_a,
      "threads"_a.none(), release,
      "Spot diagram of `field`; wavelength None means polychromatic (model weights).");
  m.def(
      "ray_fan",
      [](const compile::CompiledSystem& s, const PathArg& path, std::uint16_t field,
         std::optional<std::uint16_t> wavelength, int points, trace::Aiming aiming,
         std::optional<int> threads) {
        const FanOptions options{points, aiming};
        return with_threads(threads, [&] {
          return analysis::ray_fan(s, path_id(s, path), field, wavelength_index(s, wavelength),
                                   options);
        });
      },
      "system"_a, "path"_a, "field"_a, "wavelength"_a.none(), "points"_a, "aiming"_a,
      "threads"_a.none(), release, "Ray fans of `field`; wavelength None means the reference.");
  m.def(
      "opd_map",
      [](const compile::CompiledSystem& s, const PathArg& path, std::uint16_t field,
         std::optional<std::uint16_t> wavelength, int grid, trace::Aiming aiming,
         std::optional<int> threads) {
        OpdOptions options;
        options.grid = grid;
        options.aiming = aiming;
        return with_threads(threads, [&] {
          return analysis::opd_map(s, path_id(s, path), field, wavelength_index(s, wavelength),
                                   options);
        });
      },
      "system"_a, "path"_a, "field"_a, "wavelength"_a.none(), "grid"_a, "aiming"_a,
      "threads"_a.none(), release, "OPD map of `field`; wavelength None means the reference.");
  m.def(
      "opd_fan",
      [](const compile::CompiledSystem& s, const PathArg& path, std::uint16_t field,
         std::optional<std::uint16_t> wavelength, int points, trace::Aiming aiming,
         std::optional<int> threads) {
        OpdOptions options;
        options.fan_points = points;
        options.aiming = aiming;
        return with_threads(threads, [&] {
          return analysis::opd_fan(s, path_id(s, path), field, wavelength_index(s, wavelength),
                                   options);
        });
      },
      "system"_a, "path"_a, "field"_a, "wavelength"_a.none(), "points"_a, "aiming"_a,
      "threads"_a.none(), release, "OPD fans of `field`; wavelength None means the reference.");
  m.def(
      "longitudinal_colour",
      [](const compile::CompiledSystem& s, const PathArg& path,
         std::optional<paraxial::ChromaticPair> pair, double zone, trace::Aiming aiming,
         std::optional<int> threads) {
        const ChromaticOptions options{pair, zone, aiming};
        return with_threads(
            threads, [&] { return analysis::longitudinal_colour(s, path_id(s, path), options); });
      },
      "system"_a, "path"_a, "pair"_a.none(), "zone"_a, "aiming"_a, "threads"_a.none(), release,
      "Longitudinal colour; pair None means first and last system wavelength.");
  m.def(
      "lateral_colour",
      [](const compile::CompiledSystem& s, const PathArg& path, std::uint16_t field,
         trace::Aiming aiming, std::optional<int> threads) {
        return with_threads(
            threads, [&] { return analysis::lateral_colour(s, path_id(s, path), field, aiming); });
      },
      "system"_a, "path"_a, "field"_a, "aiming"_a, "threads"_a.none(), release,
      "Lateral colour of `field`.");
  m.def(
      "distortion",
      [](const compile::CompiledSystem& s, const PathArg& path,
         std::optional<std::uint16_t> wavelength, int samples, trace::Aiming aiming,
         std::optional<int> threads) {
        const FieldSweepOptions options{samples, aiming};
        return with_threads(threads, [&] {
          return DistortionSweep{
              analysis::distortion(s, path_id(s, path), wavelength_index(s, wavelength), options)};
        });
      },
      "system"_a, "path"_a, "wavelength"_a.none(), "samples"_a, "aiming"_a, "threads"_a.none(),
      release, "Distortion over the field sweep; wavelength None means the reference.");
  m.def(
      "distortion_at",
      [](const compile::CompiledSystem& s, const PathArg& path, const model::Field& field,
         std::optional<std::uint16_t> wavelength, trace::Aiming aiming,
         std::optional<int> threads) {
        return with_threads(threads, [&] {
          return analysis::distortion_at(s, path_id(s, path), field,
                                         wavelength_index(s, wavelength), aiming);
        });
      },
      "system"_a, "path"_a, "field"_a, "wavelength"_a.none(), "aiming"_a, "threads"_a.none(),
      release, "Distortion at one field value; wavelength None means the reference.");
  m.def(
      "field_curvature",
      [](const compile::CompiledSystem& s, const PathArg& path,
         std::optional<std::uint16_t> wavelength, int samples, double delta, trace::Aiming aiming,
         std::optional<int> threads) {
        const FieldCurvatureOptions options{samples, delta, aiming};
        return with_threads(threads, [&] {
          return FieldCurvatureSweep{analysis::field_curvature(
              s, path_id(s, path), wavelength_index(s, wavelength), options)};
        });
      },
      "system"_a, "path"_a, "wavelength"_a.none(), "samples"_a, "delta"_a, "aiming"_a,
      "threads"_a.none(), release,
      "Field curvature over the field sweep; wavelength None means the reference.");
  m.def(
      "field_curvature_at",
      [](const compile::CompiledSystem& s, const PathArg& path, const model::Field& field,
         std::optional<std::uint16_t> wavelength, double delta, trace::Aiming aiming,
         std::optional<int> threads) {
        FieldCurvatureOptions options;
        options.delta = delta;
        options.aiming = aiming;
        return with_threads(threads, [&] {
          return analysis::field_curvature_at(s, path_id(s, path), field,
                                              wavelength_index(s, wavelength), options);
        });
      },
      "system"_a, "path"_a, "field"_a, "wavelength"_a.none(), "delta"_a, "aiming"_a,
      "threads"_a.none(), release,
      "Field curvature at one field value; wavelength None means the reference.");
}

}  // namespace rtt::py
