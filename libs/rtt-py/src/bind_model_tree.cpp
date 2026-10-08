// Read access to the model tree (#82, ADR 0024 point 1): immutable typed copies of every model
// type below System. The classes have no constructor and only read-only attributes; System
// hands out copies (bind_model.cpp), so an object stays valid when the System changes or is
// freed.

#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <string>
#include <variant>
#include <vector>

#include "bindings.hpp"
#include "rtt/model/element.hpp"
#include "rtt/model/param.hpp"
#include "rtt/model/path.hpp"
#include "rtt/model/pose.hpp"
#include "rtt/model/surface.hpp"
#include "rtt/model/system.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

/// A read-only class: value equality, not hashable (it holds lists), no constructor.
template <typename T>
nb::class_<T> read_only_class(nb::module_& m, const char* name, const char* doc) {
  nb::class_<T> cls(m, name, doc);
  cls.def("__eq__", &equal<T>, "other"_a);
  cls.attr("__hash__") = nb::none();
  return cls;
}

/// Children of an assembly as the alternatives themselves (Assembly or Element).
std::vector<std::variant<model::Assembly, model::Element>> children(const model::Assembly& a) {
  std::vector<std::variant<model::Assembly, model::Element>> out;
  out.reserve(a.children.size());
  for (const model::Node& n : a.children) {
    out.push_back(std::visit(
        [](const auto& v) -> std::variant<model::Assembly, model::Element> { return v; }, n.value));
  }
  return out;
}

}  // namespace

void bind_model_tree(nb::module_& m) {
  nb::enum_<model::ElementKind>(m, "ElementKind", "Kind of an element (file key \"type\").")
      .value("LENS", model::ElementKind::Lens,
             ">= 2 surfaces, one material per segment (cemented groups allowed)")
      .value("MIRROR", model::ElementKind::Mirror,
             ">= 1 surface, substrate material optional (shorthand `material` only)")
      .value("PLATE", model::ElementKind::Plate,
             ">= 2 plane surfaces, one material per segment (plates, cubes, prisms)")
      .value("THIN_ELEMENT", model::ElementKind::ThinElement,
             "exactly 1 surface without thickness (ideal polarizer, retarder, splitter)")
      .value("STOP", model::ElementKind::Stop, "exactly 1 surface, aperture stop")
      .value("DETECTOR", model::ElementKind::Detector, "exactly 1 surface");

  nb::enum_<model::EventKind>(m, "EventKind", "What a ray does at the surface of an event.")
      .value("REFRACT", model::EventKind::Refract, "refraction into the next medium")
      .value("REFLECT", model::EventKind::Reflect, "reflection")
      .value("TRANSMIT", model::EventKind::Transmit,
             "pass without refraction (thin element, detector)")
      .value("ORDINARY", model::EventKind::Ordinary, "ordinary ray in a uniaxial crystal")
      .value("EXTRAORDINARY", model::EventKind::Extraordinary,
             "extraordinary ray in a uniaxial crystal")
      .value("DIFFRACT", model::EventKind::Diffract, "diffraction into the event's order");

  nb::enum_<model::FieldType>(m, "FieldType", "Meaning of the field values x and y.")
      .value("ANGLE_DEG", model::FieldType::AngleDeg, "field angle in degree")
      .value("OBJECT_HEIGHT", model::FieldType::ObjectHeight, "object height in mm")
      .value("PARAXIAL_IMAGE_HEIGHT", model::FieldType::ParaxialImageHeight,
             "paraxial image height in mm");

  nb::enum_<model::SystemApertureType>(m, "SystemApertureType",
                                       "What the value of the system aperture means.")
      .value("ENTRANCE_PUPIL_DIAMETER", model::SystemApertureType::EntrancePupilDiameter,
             "entrance pupil diameter in mm (file: \"epd\")")
      .value("IMAGE_SPACE_F_NUMBER", model::SystemApertureType::ImageSpaceFNumber,
             "image-space F-number, dimensionless (file: \"image_fnumber\")")
      .value("OBJECT_SPACE_NA", model::SystemApertureType::ObjectSpaceNA,
             "object-space numerical aperture, dimensionless (file: \"object_na\")")
      .value("STOP_SIZE", model::SystemApertureType::StopSize,
             "the stop aperture as defined; the value is not used (file: \"stop_size\")");

  read_only_class<model::Param>(
      m, "Param",
      "Numeric design parameter (read-only copy). Its unit is that of the attribute that holds "
      "it.")
      .def_ro("value", &model::Param::value, "Current value.")
      .def_ro("variable", &model::Param::variable, "True if the optimizer may change it.")
      .def_ro("pickup", &model::Param::pickup,
              "Pickup expression (evaluated from M5 on); None for an independent value.")
      .def("__repr__", [](const model::Param& p) {
        nb::str text = nb::str("Param(value={!r}").format(p.value);
        if (p.variable) text = nb::str("{}, variable=True").format(text);
        if (p.pickup) text = nb::str("{}, pickup={!r}").format(text, *p.pickup);
        return nb::str("{})").format(text);
      });

  read_only_class<model::Pose>(
      m, "Pose",
      "Placement of a node in its parent coordinate system (read-only copy): translation, then "
      "intrinsic X -> Y -> Z rotation about the pivot, p_parent = position + pivot + "
      "R (p_child - pivot) with R = Rx Ry Rz. Right-handed coordinates, +z is the optical "
      "axis, y the meridional direction.")
      .def_ro("position", &model::Pose::position,
              "Translation x, y, z in mm, in parent coordinates (Param each).")
      .def_ro("rotation_deg", &model::Pose::rotation_deg,
              "Rotation angles about x, y, z in degree, intrinsic X -> Y -> Z (Param each).")
      .def_ro("pivot", &model::Pose::pivot,
              "Pivot point x, y, z in mm, in local coordinates; tolerances tilt about it.");

  // ---------------------------------------------------------------- shape -----
  read_only_class<model::Plane>(m, "Plane", "Flat surface z = 0 (read-only copy).");
  read_only_class<model::Conic>(
      m, "Conic",
      "Conic of revolution, a sphere for conic = 0 (read-only copy). The radius is positive if "
      "the centre of curvature lies on the +z side of the vertex.")
      .def_ro("radius", &model::Conic::radius, "Vertex radius of curvature in mm (Param).")
      .def_ro("conic", &model::Conic::conic, "Conic constant, dimensionless (Param).");
  read_only_class<model::EvenAsphere>(
      m, "EvenAsphere",
      "Conic plus even polynomial z = conic(r) + A4 r^4 + A6 r^6 + ... (read-only copy).")
      .def_ro("radius", &model::EvenAsphere::radius,
              "Vertex radius of curvature in mm, sign as for Conic (Param).")
      .def_ro("conic", &model::EvenAsphere::conic, "Conic constant, dimensionless (Param).")
      .def_ro("coefficients", &model::EvenAsphere::coefficients,
              "A4, A6, A8, ... in mm^(1-2k) for the term r^(2k) (Param each).");
  read_only_class<model::ZernikeSag>(
      m, "ZernikeSag",
      "Additive Zernike sag term in Noll ordering, coefficients[0] is j = 1 (read-only copy).")
      .def_ro("normalization_radius", &model::ZernikeSag::normalization_radius,
              "Normalization radius in mm (Param).")
      .def_ro("coefficients", &model::ZernikeSag::coefficients,
              "Zernike coefficients in mm (Param each).");
  read_only_class<model::ShapeStack>(
      m, "ShapeStack",
      "Total sag = base + sum of the terms, in the local xy plane (read-only copy).")
      .def_ro("base", &model::ShapeStack::base, "Base shape: Plane, Conic or EvenAsphere.")
      .def_ro("terms", &model::ShapeStack::terms, "Additive terms (ZernikeSag).");

  // ------------------------------------------------------------- aperture -----
  read_only_class<model::CircularAperture>(m, "CircularAperture",
                                           "Circular aperture, an annulus if inner_radius > 0 "
                                           "(read-only copy).")
      .def_ro("radius", &model::CircularAperture::radius, "Outer radius in mm.")
      .def_ro("inner_radius", &model::CircularAperture::inner_radius,
              "Inner radius in mm; 0 for a full disc.");
  read_only_class<model::RectangularAperture>(m, "RectangularAperture",
                                              "Rectangular aperture (read-only copy).")
      .def_ro("half_width_x", &model::RectangularAperture::half_width_x, "Half width in x, mm.")
      .def_ro("half_width_y", &model::RectangularAperture::half_width_y, "Half width in y, mm.");
  read_only_class<model::EllipticalAperture>(m, "EllipticalAperture",
                                             "Elliptical aperture (read-only copy).")
      .def_ro("semi_axis_x", &model::EllipticalAperture::semi_axis_x, "Semi-axis in x, mm.")
      .def_ro("semi_axis_y", &model::EllipticalAperture::semi_axis_y, "Semi-axis in y, mm.");

  // --------------------------------------------------------------- phases -----
  read_only_class<model::LinearGrating>(
      m, "LinearGrating",
      "Straight-line grating; orientation_deg = 0 means grooves parallel to the local y axis "
      "(read-only copy).")
      .def_ro("lines_per_mm", &model::LinearGrating::lines_per_mm,
              "Groove density in lines per mm (Param).")
      .def_ro("orientation_deg", &model::LinearGrating::orientation_deg,
              "Groove orientation in degree.");
  read_only_class<model::RadialPhase>(
      m, "RadialPhase",
      "Rotationally symmetric phase phi = sum c_k rho^(2k), rho = r / normalization_radius "
      "(read-only copy).")
      .def_ro("normalization_radius", &model::RadialPhase::normalization_radius,
              "Normalization radius in mm (Param).")
      .def_ro("coefficients", &model::RadialPhase::coefficients,
              "c_1, c_2, ... in radian; coefficients[0] belongs to rho^2 (Param each).");

  // ---------------------------------------------------------- interaction -----
  read_only_class<model::Fresnel>(
      m, "Fresnel",
      "Uncoated interface, Fresnel coefficients from the two media; the default (read-only "
      "copy).");
  read_only_class<model::IdealMirror>(
      m, "IdealMirror",
      "Lossless mirror, R = 1 for s and p, no phase difference (read-only copy).");
  read_only_class<model::IdealAntiReflection>(
      m, "IdealAntiReflection", "Lossless anti-reflection interface, T = 1 (read-only copy).");
  read_only_class<model::Absorber>(m, "Absorber", "Fully absorbing surface (read-only copy).");
  read_only_class<model::IdealBeamSplitter>(
      m, "IdealBeamSplitter",
      "Lossless beam splitter, T = 1 - R per polarization (read-only copy).")
      .def_ro("reflectance_s", &model::IdealBeamSplitter::reflectance_s,
              "Intensity reflectance for s polarization, dimensionless.")
      .def_ro("reflectance_p", &model::IdealBeamSplitter::reflectance_p,
              "Intensity reflectance for p polarization, dimensionless.");
  read_only_class<model::CoatingRef>(
      m, "CoatingRef", "Thin-film coating defined in a coating catalog (read-only copy).")
      .def_ro("name", &model::CoatingRef::name,
              "Catalog reference of the coating, e.g. \"DEMO:AR_MGF2\".");
  read_only_class<model::IdealPolarizer>(
      m, "IdealPolarizer",
      "Ideal linear polarizer; the axis is projected perpendicular to the ray direction "
      "(read-only copy).")
      .def_ro("transmission_axis", &model::IdealPolarizer::transmission_axis,
              "Transmission axis x, y, z in element coordinates, dimensionless.")
      .def_ro("extinction_ratio", &model::IdealPolarizer::extinction_ratio,
              "T_min / T_max, dimensionless; 0 is a perfect polarizer.");
  read_only_class<model::IdealRetarder>(
      m, "IdealRetarder",
      "Ideal linear retarder; a positive retardance delays the slow axis (read-only copy).")
      .def_ro("fast_axis", &model::IdealRetarder::fast_axis,
              "Fast axis x, y, z in element coordinates, dimensionless.")
      .def_ro("retardance_waves", &model::IdealRetarder::retardance_waves,
              "Retardance in waves (0.25 for a quarter-wave plate).");

  // -------------------------------------------------------------- surface -----
  read_only_class<model::Surface>(m, "Surface",
                                  "Surface with its layer stack (read-only copy): shape, "
                                  "aperture, phases, interaction.")
      .def_prop_ro(
          "id", [](const model::Surface& s) { return s.id.str(); },
          "Stable surface id, unique in the system; paths refer to it.")
      .def_ro("pose", &model::Surface::pose, "Placement in the element.")
      .def_ro("shape", &model::Surface::shape, "Sag of the surface.")
      .def_ro("aperture", &model::Surface::aperture,
              "Clear aperture (CircularAperture, RectangularAperture, EllipticalAperture); None "
              "for an unlimited surface.")
      .def_ro("phases", &model::Surface::phases, "Phase layers (LinearGrating, RadialPhase).")
      .def_ro("interaction", &model::Surface::interaction,
              "Interaction (Fresnel, IdealMirror, IdealAntiReflection, IdealBeamSplitter, "
              "CoatingRef, IdealPolarizer, IdealRetarder, Absorber).")
      .def("__repr__",
           [](const model::Surface& s) { return nb::str("Surface(id={!r})").format(s.id.str()); });

  // -------------------------------------------------------------- element -----
  read_only_class<model::Element>(
      m, "Element",
      "A physical body (read-only copy). An element with N surfaces has N - 1 segments; segment i "
      "is the glass between surfaces i and i + 1. The material is given either as `material` "
      "(one catalog reference for all segments) or as `segment_materials` (one per segment).")
      .def_ro("name", &model::Element::name, "Name, unique among assemblies and elements.")
      .def_ro("kind", &model::Element::kind, "Kind of the element.")
      .def_ro("pose", &model::Element::pose, "Placement in the parent assembly.")
      .def_ro("material", &model::Element::material,
              "Catalog reference for all segments, e.g. \"SCHOTT:N-BK7\"; None if the element "
              "has no material or uses segment_materials.")
      .def_ro("segment_materials", &model::Element::segment_materials,
              "Catalog reference per segment; empty if `material` is used or the element has no "
              "material.")
      .def_ro("surfaces", &model::Element::surfaces, "Surfaces in order.")
      .def("segment_material", &model::Element::segment_material, "segment"_a,
           "Material of segment `segment` (between surfaces `segment` and `segment + 1`): the "
           "list entry or the shorthand `material`. None if the element has no material or the "
           "segment does not exist.")
      .def("__repr__", [](const model::Element& e) {
        return nb::str("Element(name={!r}, kind={})").format(e.name, nb::cast(e.kind));
      });

  read_only_class<model::Assembly>(
      m, "Assembly", "Group of nodes moved and toleranced as one rigid body (read-only copy).")
      .def_ro("name", &model::Assembly::name, "Name, unique among assemblies and elements.")
      .def_ro("pose", &model::Assembly::pose, "Placement in the parent assembly.")
      .def_prop_ro("children", &children,
                   "Child nodes in order: Assembly or Element. Each access copies the subtrees; "
                   "keep the list in a variable instead of reading it again in a loop.")
      .def("__repr__", [](const model::Assembly& a) {
        return nb::str("Assembly(name={!r}, children={})").format(a.name, a.children.size());
      });

  // ---------------------------------------------------------------- paths -----
  read_only_class<model::Event>(m, "Event", "One event of a ray path (read-only copy).")
      .def_prop_ro(
          "surface", [](const model::Event& e) { return e.surface.str(); }, "Id of the surface.")
      .def_ro("kind", &model::Event::kind, "What the ray does at the surface.")
      .def_ro("order", &model::Event::order,
              "Diffraction order; only meaningful for EventKind.DIFFRACT.");

  read_only_class<model::Path>(
      m, "Path",
      "Ray path (read-only copy): either automatic (all surfaces in tree order; refracting at "
      "lenses and plates, reflecting at mirrors, transmitting at stops, detectors and thin "
      "elements) or a list of events.")
      .def_ro("name", &model::Path::name, "Name, unique among the paths.")
      .def_ro("automatic", &model::Path::automatic,
              "True for the automatic path (file: \"events\": \"auto\"); `events` is then empty.")
      .def_ro("events", &model::Path::events, "Events in order; empty for an automatic path.")
      .def("__repr__", [](const model::Path& p) {
        return nb::str("Path(name={!r}, automatic={!r}, events={})")
            .format(p.name, p.automatic, p.events.size());
      });

  // ------------------------------------------------------- system parts -----
  read_only_class<model::FieldSet>(m, "FieldSet", "Field points of the system (read-only copy).")
      .def_ro("type", &model::FieldSet::type,
              "Meaning of x and y: angle in degree or height in mm.")
      .def_prop_ro(
          "points", [](const model::FieldSet& f) { return f.points; },
          "Field points (copies; Field is mutable, but changing one does not change this "
          "FieldSet or the System).");

  read_only_class<model::SystemAperture>(m, "SystemAperture",
                                         "Aperture of the system (read-only copy).")
      .def_ro("type", &model::SystemAperture::type, "What `value` means.")
      .def_ro("value", &model::SystemAperture::value,
              "mm for ENTRANCE_PUPIL_DIAMETER, dimensionless for IMAGE_SPACE_F_NUMBER and "
              "OBJECT_SPACE_NA, not used for STOP_SIZE (Param).");

  read_only_class<model::ObjectSpace>(m, "ObjectSpace", "Position of the object (read-only copy).")
      .def_ro("at_infinity", &model::ObjectSpace::at_infinity, "True for an object at infinity.")
      .def_ro("distance", &model::ObjectSpace::distance,
              "Distance from the object to the global origin along -z in mm; used only if not "
              "at_infinity (Param).");
}

}  // namespace rtt::py
