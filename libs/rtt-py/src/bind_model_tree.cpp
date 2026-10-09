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
#include "rtt/model/optimization.hpp"
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

/// target, weight and configuration of an operand as attributes of its own (file keys of the
/// same names, ADR 0030 point 2).
template <typename T>
void common_attributes(nb::class_<T>& cls) {
  cls.def_prop_ro(
         "target", [](const T& o) { return o.common.target; },
         "Target in the unit of the operand's value.")
      .def_prop_ro(
          "weight", [](const T& o) { return o.common.weight; },
          "Weight w >= 0 of the residual sqrt(w) (value - target).")
      .def_prop_ro(
          "configuration", [](const T& o) { return o.common.configuration; },
          "Name of a configuration; None: configuration 0.");
}

/// The merit function (ADR 0030, points 2-4; #162 part B). The defaults are those of
/// rtt/model/optimization.hpp.
void bind_optimization(nb::module_& m) {
  nb::enum_<model::FirstOrderQuantity>(m, "FirstOrderQuantity",
                                       "Paraxial quantity of a first-order operand.")
      .value("EFL", model::FirstOrderQuantity::Efl, "effective focal length, mm (file: \"efl\")")
      .value("BFL", model::FirstOrderQuantity::Bfl, "back focal length, mm (file: \"bfl\")")
      .value("IMAGE_F_NUMBER", model::FirstOrderQuantity::ImageFNumber,
             "image-space F-number, dimensionless (file: \"image_fnumber\")")
      .value("MAGNIFICATION", model::FirstOrderQuantity::Magnification,
             "paraxial magnification, dimensionless (file: \"magnification\")");
  nb::enum_<model::RayCoordinate>(m, "RayCoordinate", "Coordinate of a ray operand.")
      .value("X", model::RayCoordinate::X, "x in mm (file: \"ray_x\")")
      .value("Y", model::RayCoordinate::Y, "y in mm (file: \"ray_y\")");
  nb::enum_<model::SpotReference>(m, "SpotReference", "Reference point of a spot RMS.")
      .value("CENTROID", model::SpotReference::Centroid, "the centroid (file: \"centroid\")")
      .value("CHIEF", model::SpotReference::Chief, "the chief ray (file: \"chief\")");

  auto first_order = read_only_class<model::FirstOrderOperand>(
      m, "FirstOrderOperand",
      "A paraxial quantity of a path (file types efl, bfl, image_fnumber, magnification; "
      "read-only copy).");
  common_attributes(first_order);
  first_order.def_ro("quantity", &model::FirstOrderOperand::quantity, "The quantity.")
      .def_ro("path", &model::FirstOrderOperand::path, "Name of a path.")
      .def_ro("wavelength", &model::FirstOrderOperand::wavelength,
              "Wavelength index; None: the reference wavelength.");

  auto ray = read_only_class<model::RayOperand>(
      m, "RayOperand",
      "x or y of a real ray in mm in the local coordinates of a surface (file types ray_x, "
      "ray_y; read-only copy).");
  common_attributes(ray);
  ray.def_ro("coordinate", &model::RayOperand::coordinate, "X or Y.")
      .def_ro("path", &model::RayOperand::path, "Name of a path.")
      .def_prop_ro(
          "surface", [](const model::RayOperand& o) { return o.surface.str(); },
          "Id of a surface on the path.")
      .def_ro("occurrence", &model::RayOperand::occurrence,
              "Which event at the surface (0-based) if the path meets it more than once; None: "
              "it meets it once.")
      .def_ro("field", &model::RayOperand::field, "Field index.")
      .def_ro("px", &model::RayOperand::px,
              "Normalised pupil x (unit circle: rim of the paraxial entrance pupil).")
      .def_ro("py", &model::RayOperand::py, "Normalised pupil y, +y meridional.")
      .def_ro("wavelength", &model::RayOperand::wavelength,
              "Wavelength index; None: the reference wavelength.");

  auto spot = read_only_class<model::SpotRmsOperand>(
      m, "SpotRmsOperand",
      "RMS spot radius of one field in mm (file type spot_rms; read-only copy).");
  common_attributes(spot);
  spot.def_ro("path", &model::SpotRmsOperand::path, "Name of a path.")
      .def_ro("field", &model::SpotRmsOperand::field, "Field index.")
      .def_ro("wavelength", &model::SpotRmsOperand::wavelength,
              "Wavelength index; None: the reference wavelength, or all if polychromatic.")
      .def_ro("polychromatic", &model::SpotRmsOperand::polychromatic,
              "True: all wavelengths with their weights.")
      .def_ro("reference", &model::SpotRmsOperand::reference, "Centroid or chief ray.")
      .def_ro("rings", &model::SpotRmsOperand::rings, "Rings of the hexapolar pupil.");

  auto opd = read_only_class<model::OpdRmsOperand>(
      m, "OpdRmsOperand",
      "RMS wavefront error of one field in waves (file type opd_rms; read-only copy).");
  common_attributes(opd);
  opd.def_ro("path", &model::OpdRmsOperand::path, "Name of a path.")
      .def_ro("field", &model::OpdRmsOperand::field, "Field index.")
      .def_ro("wavelength", &model::OpdRmsOperand::wavelength,
              "Wavelength index; None: the reference wavelength.")
      .def_ro("grid", &model::OpdRmsOperand::grid, "Points per side of the pupil grid.");

  auto value = read_only_class<model::ParamValueOperand>(
      m, "ParamValueOperand",
      "Value of a parameter row in a configuration (file type param_value; read-only copy).");
  common_attributes(value);
  value.def_ro("parameter", &model::ParamValueOperand::parameter, "Name of a parameter row.");

  read_only_class<model::SpotGenerator>(
      m, "SpotGenerator",
      "Spot generator: two residuals per ray, target 0 (file type rms_spot; read-only copy).")
      .def_ro("path", &model::SpotGenerator::path, "Name of a path.")
      .def_ro("configuration", &model::SpotGenerator::configuration,
              "Name of a configuration; None: configuration 0.")
      .def_ro("fields", &model::SpotGenerator::fields, "Field indices; None: all fields.")
      .def_ro("wavelengths", &model::SpotGenerator::wavelengths,
              "Wavelength indices; None: all wavelengths.")
      .def_ro("reference", &model::SpotGenerator::reference, "Centroid or chief ray.")
      .def_ro("rings", &model::SpotGenerator::rings, "Gauss-Legendre rings in rho^2.")
      .def_ro("arms", &model::SpotGenerator::arms, "Arms per ring.")
      .def_ro("weight", &model::SpotGenerator::weight, "Weight >= 0.");

  read_only_class<model::WavefrontGenerator>(
      m, "WavefrontGenerator",
      "Wavefront generator: one residual per ray in waves, target 0 (file type rms_wavefront; "
      "read-only copy).")
      .def_ro("path", &model::WavefrontGenerator::path, "Name of a path.")
      .def_ro("configuration", &model::WavefrontGenerator::configuration,
              "Name of a configuration; None: configuration 0.")
      .def_ro("fields", &model::WavefrontGenerator::fields, "Field indices; None: all fields.")
      .def_ro("wavelengths", &model::WavefrontGenerator::wavelengths,
              "Wavelength indices; None: all wavelengths.")
      .def_ro("rings", &model::WavefrontGenerator::rings, "Gauss-Legendre rings in rho^2.")
      .def_ro("arms", &model::WavefrontGenerator::arms, "Arms per ring.")
      .def_ro("weight", &model::WavefrontGenerator::weight, "Weight >= 0.");

  read_only_class<model::Optimization>(m, "Optimization",
                                       "The merit function (ADR 0030; read-only copy).")
      .def_ro("operands", &model::Optimization::operands, "Operands in file order.")
      .def_ro("generators", &model::Optimization::generators, "Generators in file order.");
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
             "extraordinary ray in a uniaxial crystal");

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
             "image-space F-number at infinite conjugates, EPD = |EFL| / F#, also for a finite "
             "object; dimensionless (file: \"image_fnumber\")")
      .value("OBJECT_SPACE_NA", model::SystemApertureType::ObjectSpaceNA,
             "object-space numerical aperture NA, dimensionless, finite objects only, read "
             "paraxially: the paraxial marginal slope from the axial object point is u = NA / n, "
             "n the index of object space (Greivenkamp, OPTI-502, Sec. 9, p. 9-34: "
             "NA = n sin U ~ n u). With paraxial aiming the marginal ray has tan U = NA / n; with "
             "real aiming it hits the paraxial stop edge R_s, and its tan U differs from NA / n "
             "by the pupil aberration; "
             "file: \"object_na\"")
      .value("STOP_SIZE", model::SystemApertureType::StopSize,
             "the stop aperture as defined; the value is not used (file: \"stop_size\")");

  read_only_class<model::Param>(
      m, "Param",
      "Numeric design parameter (read-only copy). Its unit is that of the attribute that holds "
      "it. Either a value (unbound) or bound to a row of System.parameters (ADR 0029).")
      .def_prop_ro(
          "value",
          [](const model::Param& p) -> std::optional<double> {
            if (p.is_bound()) return std::nullopt;
            return p.value;
          },
          "Current value; None for a Param bound to a parameter row (its value comes from the "
          "table, ADR 0029).")
      .def_ro("variable", &model::Param::variable,
              "True if the optimizer may change it (unbound Params only).")
      .def_ro("min", &model::Param::min,
              "Lower bound for the optimizer, unit of the value; None if not set.")
      .def_ro("max", &model::Param::max,
              "Upper bound for the optimizer, unit of the value; None if not set.")
      .def_ro("param", &model::Param::param,
              "Name of the parameter row it is bound to (file: {\"param\": ...}); None for a "
              "value.")
      .def("__repr__", [](const model::Param& p) {
        if (p.is_bound()) return nb::str("Param(param={!r})").format(*p.param);
        nb::str text = nb::str("Param(value={!r}").format(p.value);
        if (p.variable) text = nb::str("{}, variable=True").format(text);
        if (p.min) text = nb::str("{}, min={!r}").format(text, *p.min);
        if (p.max) text = nb::str("{}, max={!r}").format(text, *p.max);
        return nb::str("{})").format(text);
      });

  nb::enum_<model::PoseReference>(
      m, "PoseReference",
      "What a Pose is relative to (ADR 0028): global(X) = global(reference) * pose(X).")
      .value("ABSOLUTE", model::PoseReference::Absolute, "the parent node (file: \"absolute\")")
      .value("RELATIVE_TO_PRECEDING", model::PoseReference::RelativeToPreceding,
             "the last surface before the node in tree order (file: "
             "\"relative_to_preceding\")")
      .value("RELATIVE_TO_SIBLING", model::PoseReference::RelativeToSibling,
             "the preceding sibling in the same assembly or element (file: "
             "\"relative_to_sibling\")");

  nb::enum_<model::PoseOrder>(m, "PoseOrder",
                              "Order of translation and rotation of a Pose (ADR 0028).")
      .value("TRANSLATE_FIRST", model::PoseOrder::TranslateFirst,
             "p_ref = t + P + R (p - P), pivot P in the node's coordinates (file: "
             "\"translate_first\")")
      .value("ROTATE_FIRST", model::PoseOrder::RotateFirst,
             "p_ref = Q + R (t + p - Q), pivot Q in the reference frame: a coordinate break "
             "(file: \"rotate_first\")");

  read_only_class<model::Pose>(
      m, "Pose",
      "Placement of a node relative to its reference (read-only copy, ADR 0028). With order "
      "TRANSLATE_FIRST: translation, then intrinsic X -> Y -> Z rotation about the pivot, "
      "p_ref = position + pivot + R (p_child - pivot) with R = Rx Ry Rz; with ROTATE_FIRST the "
      "rotation about the pivot (in the reference frame) comes first and the translation runs "
      "along the rotated axes. Right-handed coordinates, +z is the optical axis, y the "
      "meridional direction.")
      .def_ro("reference", &model::Pose::reference,
              "What the pose is relative to (PoseReference); ABSOLUTE is the parent node.")
      .def_ro("order", &model::Pose::order, "Order of translation and rotation (PoseOrder).")
      .def_ro("position", &model::Pose::position,
              "Translation x, y, z in mm, in the reference frame (TRANSLATE_FIRST) or along the "
              "rotated axes (ROTATE_FIRST) (Param each).")
      .def_ro("rotation_deg", &model::Pose::rotation_deg,
              "Rotation angles about x, y, z in degree, intrinsic X -> Y -> Z (Param each).")
      .def_ro("pivot", &model::Pose::pivot,
              "Pivot point x, y, z in mm, in the node's coordinates (TRANSLATE_FIRST) or in the "
              "reference frame (ROTATE_FIRST); tolerances tilt about it.");

  read_only_class<model::Configuration>(
      m, "Configuration", "A configuration: one column of the parameter table (ADR 0029).")
      .def_ro("name", &model::Configuration::name, "Name, unique.");

  read_only_class<model::ParameterRow>(
      m, "ParameterRow",
      "Row of the parameter table (read-only copy, ADR 0029): exactly one of value (all "
      "configurations), values (one per configuration) and expression (earlier rows); the "
      "other two are None. Its unit is that of the Params bound to it.")
      .def_ro("name", &model::ParameterRow::name, "Name, [A-Za-z_][A-Za-z0-9_]*, unique.")
      .def_prop_ro(
          "value",
          [](const model::ParameterRow& r) -> std::optional<double> {
            if (const auto* v = std::get_if<double>(&r.form)) return *v;
            return std::nullopt;
          },
          "The value in all configurations, or None.")
      .def_prop_ro(
          "values",
          [](const model::ParameterRow& r) -> std::optional<std::vector<double>> {
            if (const auto* v = std::get_if<std::vector<double>>(&r.form)) return *v;
            return std::nullopt;
          },
          "One value per configuration (order of System.configurations), or None.")
      .def_prop_ro(
          "expression",
          [](const model::ParameterRow& r) -> std::optional<std::string> {
            if (const auto* e = std::get_if<model::ParameterExpression>(&r.form)) return e->text;
            return std::nullopt;
          },
          "Expression of earlier rows, or None.")
      .def_ro("variable", &model::ParameterRow::variable,
              "True if the optimizer may change it (never for an expression).")
      .def_ro("min", &model::ParameterRow::min, "Lower bound, or None.")
      .def_ro("max", &model::ParameterRow::max, "Upper bound, or None.")
      .def("__repr__", [](const model::ParameterRow& r) {
        return nb::str("ParameterRow({!r})").format(r.name);
      });

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
              "A4, A6, A8, ...: A_n in mm^(1 - n) for the term A_n r^n (Param each).");
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
  read_only_class<model::IdealLens>(
      m, "IdealLens",
      "Ideal lens without thickness (ADR 0031): images the plane at s onto the plane at s' with "
      "1/s' = 1/s + 1/f for every ray, not only paraxially. The OPD is exact only for objects in "
      "the plane of object_distance (read-only copy).")
      .def_ro("focal_length", &model::IdealLens::focal_length,
              "f in mm along the propagation in the surrounding medium, f > 0 converges (Param).")
      .def_ro("object_distance", &model::IdealLens::object_distance,
              "Design conjugate for the optical path in mm, positive for a real object in front "
              "(Param); None: object at infinity.");
  read_only_class<model::IdealCylinderLens>(
      m, "IdealCylinderLens",
      "Ideal cylinder lens (ADR 0031): the ideal lens in the direction of power "
      "(-sin psi, cos psi, 0) only; the slope along the axis (cos psi, sin psi, 0) stays "
      "(read-only copy).")
      .def_ro("focal_length", &model::IdealCylinderLens::focal_length,
              "f in mm along the propagation, f > 0 converges (Param).")
      .def_ro("axis_deg", &model::IdealCylinderLens::axis_deg,
              "Cylinder axis psi in degree, local surface frame, from the x to the y axis.")
      .def_ro("object_distance", &model::IdealCylinderLens::object_distance,
              "Design conjugate for the optical path in mm (Param); None: object at infinity.");

  // -------------------------------------------------------------- surface -----
  read_only_class<model::DiffractionEfficiency>(
      m, "DiffractionEfficiency",
      "Power fraction of one diffraction order at a surface with phase layers (ADR 0025; "
      "read-only copy).")
      .def_ro("order", &model::DiffractionEfficiency::order, "Diffraction order.")
      .def_ro("efficiency", &model::DiffractionEfficiency::efficiency,
              "Power fraction of the order, 0 ... 1, dimensionless.");
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
      .def_ro("diffraction_efficiency", &model::Surface::diffraction_efficiency,
              "Efficiency per diffraction order (DiffractionEfficiency); None: every order has "
              "efficiency 1. If given, orders not listed have efficiency 0.")
      .def_ro("interaction", &model::Surface::interaction,
              "Interaction (Fresnel, IdealMirror, IdealAntiReflection, IdealBeamSplitter, "
              "CoatingRef, IdealPolarizer, IdealRetarder, Absorber, IdealLens, "
              "IdealCylinderLens).")
      .def("__repr__",
           [](const model::Surface& s) { return nb::str("Surface(id={!r})").format(s.id.str()); });

  // -------------------------------------------------------------- element -----
  read_only_class<model::CrystalMaterial>(
      m, "CrystalMaterial",
      "Uniaxial crystal (ADR 0026): catalog references for the principal indices n_O and n_E "
      "(read-only copy).")
      .def_ro("ordinary", &model::CrystalMaterial::ordinary,
              "Catalog reference for n_O, e.g. \"BIREFRINGENT:CALCITE\".")
      .def_ro("extraordinary", &model::CrystalMaterial::extraordinary,
              "Catalog reference for n_E, e.g. \"BIREFRINGENT:CALCITE-E\".");
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
      .def_ro("crystal", &model::Element::crystal,
              "Uniaxial crystal for all segments (CrystalMaterial); None for an isotropic "
              "material or none. Then `material` and `segment_materials` are empty.")
      .def_ro("optic_axis", &model::Element::optic_axis,
              "Optic axis x, y, z of the crystal in element coordinates, dimensionless, not "
              "normalised; None without crystal.")
      .def_ro("surfaces", &model::Element::surfaces, "Surfaces in order.")
      .def("segment_material", &model::Element::segment_material, "segment"_a,
           "Material of segment `segment` (between surfaces `segment` and `segment + 1`): the "
           "list entry or the shorthand `material`. None if the element has no isotropic "
           "material (none, or a crystal) or the segment does not exist.")
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
              "Diffraction order (ADR 0025); non-zero only at a surface with a phase layer.");

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

  bind_optimization(m);
}

}  // namespace rtt::py
