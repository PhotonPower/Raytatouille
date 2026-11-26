// Paraxial model of paths with plane fold mirrors (#357, ADR 0054). A plane mirror maps a point
// to its mirror image at the mirror plane; unfolded at each reflection the ray runs straight on
// (J. E. Greivenkamp, OPTI-502 lecture notes (2019), Sec. 2, pp. 2-2 and 2-11). So a folded
// path has the first-order data of its unfolded counterpart: here m1/singlet_const (stop, a
// singlet, the image 106.363 mm behind the origin) with plane mirrors between the lens and the
// image. Acceptance (decision 3 on #357): EFL, BFL, positions and diameters of the pupils to
// 1e-12 relative against the unfolded system; unfolded paths stay bitwise.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <iomanip>
#include <numbers>
#include <optional>
#include <sstream>
#include <string>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/math/types.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/prescription.hpp"
#include "rtt/paraxial/seidel.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::math::Vec3;
using rtt::paraxial::FirstOrder;
using rtt::paraxial::ParaxialError;

namespace {

constexpr double kImageZ = 106.363;  ///< image of m1/singlet_const, unfolded
constexpr double kMirrorZ = 49.0;    ///< first fold mirror, 40 mm behind the lens

std::string num(double v) {
  std::ostringstream s;
  s << std::setprecision(17) << v;
  return s.str();
}

std::string vec(double x, double y, double z) {
  return "[" + num(x) + ", " + num(y) + ", " + num(z) + "]";
}

/// m1/singlet_const with `middle` (mirror elements as JSON) between the lens and the image, the
/// image at `image_pose`, the events `events` after the lens and the aperture `aperture`.
std::string singlet(const std::string& middle,
                    const std::string& image_pose,
                    const std::string& events,
                    const std::string& aperture = R"({"type": "epd", "value": 20.0})",
                    bool stop_element = true) {
  const std::string stop = stop_element
                               ? R"({"type": "stop", "name": "stop", "surfaces": [{"id": "STO",
          "aperture": {"type": "circular", "radius": 10.0}}]},)"
                               : "";
  const std::string stop_event = stop_element ? R"({"surface": "STO", "kind": "transmit"},)" : "";
  return R"({"schema_version": "0.11.0", "name": "fold", "units": {"length": "mm",
    "wavelength": "um"}, "wavelengths": [{"um": 0.5876, "reference": true}],
    "aperture": )" +
         aperture + R"(, "fields": {"points": [{}, {"y": 5.0}]},
    "root": {"type": "assembly", "name": "root", "children": [)" +
         stop + R"(
      {"type": "lens", "name": "L1", "pose": {"position": [0.0, 0.0, 5.0]},
       "material": "CONST:1.5168", "surfaces": [
        {"id": "L1.S1", "shape": {"base": {"type": "conic", "radius": 51.68}},
         "aperture": {"type": "circular", "radius": 12.7}},
        {"id": "L1.S2", "pose": {"position": [0.0, 0.0, 4.0]},
         "aperture": {"type": "circular", "radius": 12.7}}]},)" +
         middle + R"(
      {"type": "detector", "name": "image", "pose": )" +
         image_pose + R"(, "surfaces": [{"id": "IMG"}]}]},
    "paths": [{"name": "main", "events": [)" +
         stop_event + R"(
      {"surface": "L1.S1", "kind": "refract"}, {"surface": "L1.S2", "kind": "refract"},)" +
         events + R"(
      {"surface": "IMG", "kind": "transmit"}]}]})";
}

std::string mirror(const std::string& name,
                   const Vec3& position,
                   double tilt_deg,
                   const std::string& shape = {},
                   const std::string& aperture = {}) {
  return R"({"type": "mirror", "name": ")" + name + R"(", "pose": {"position": )" +
         vec(position.x(), position.y(), position.z()) + R"(, "rotation_deg": [)" + num(tilt_deg) +
         R"(, 0.0, 0.0]}, "surfaces": [{"id": ")" + name + "\"" + shape + aperture + "}]},";
}

std::string pose(const Vec3& position, double tilt_deg) {
  return R"({"position": )" + vec(position.x(), position.y(), position.z()) +
         R"(, "rotation_deg": [)" + num(tilt_deg) + ", 0.0, 0.0]}";
}

std::string reflect(const std::string& name) {
  return R"({"surface": ")" + name + R"(", "kind": "reflect"},)";
}

CompiledSystem compiled(const std::string& json) {
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(rtt::io::parse_system(json), lib);
}

const double kPi = std::numbers::pi;

/// One fold mirror at z = 49 tilted by `tilt_deg` about x: its normal is Rx(t) e_z =
/// (0, -sin t, cos t), so +z reflects into d = (0, sin 2t, -cos 2t); the image lies on d at the
/// unfolded distance, facing along d (rotation 180 deg + 2t about x).
CompiledSystem one_fold(double tilt_deg) {
  const double t = tilt_deg * kPi / 180.0;
  const Vec3 d(0.0, std::sin(2.0 * t), -std::cos(2.0 * t));
  const Vec3 image = Vec3(0.0, 0.0, kMirrorZ) + (kImageZ - kMirrorZ) * d;
  return compiled(singlet(mirror("M", Vec3(0.0, 0.0, kMirrorZ), tilt_deg),
                          pose(image, 180.0 + 2.0 * tilt_deg), reflect("M")));
}

/// A periscope: the mirror at z = 49 tilted by 45 deg sends the light to +y, a second one 10 mm
/// further, tilted by 45 deg, back to +z; two folds, even parity.
CompiledSystem periscope() {
  const Vec3 m2(0.0, 10.0, kMirrorZ);
  return compiled(singlet(mirror("M1", Vec3(0.0, 0.0, kMirrorZ), 45.0) + mirror("M2", m2, 45.0),
                          pose(Vec3(0.0, 10.0, kImageZ - 10.0), 0.0),
                          reflect("M1") + reflect("M2")));
}

CompiledSystem unfolded() {
  return compiled(singlet("", pose(Vec3(0.0, 0.0, kImageZ), 0.0), ""));
}

bool close(double a, double b) {
  return std::abs(a - b) <= 1e-12 * std::max(1.0, std::abs(b));
}

bool close(const std::optional<double>& a, const std::optional<double>& b) {
  return a.has_value() == b.has_value() && (!a || close(*a, *b));
}

/// The first-order data of a folded path against those of the unfolded singlet.
void same_first_order(const FirstOrder& f, const FirstOrder& u) {
  CHECK(f.image_direction == u.image_direction);
  CHECK(close(f.power, u.power));
  CHECK(close(f.efl, u.efl));
  CHECK(close(f.bfl, u.bfl));
  CHECK(close(f.ffl, u.ffl));
  CHECK(close(f.front_focal_z, u.front_focal_z));
  CHECK(close(f.rear_focal_z, u.rear_focal_z));
  CHECK(close(f.front_principal_z, u.front_principal_z));
  CHECK(close(f.rear_principal_z, u.rear_principal_z));
  CHECK(close(f.image_z, u.image_z));
  CHECK(close(f.angular_magnification, u.angular_magnification));
  REQUIRE(f.entrance_pupil);
  REQUIRE(f.exit_pupil);
  CHECK(close(f.entrance_pupil->z, u.entrance_pupil->z));
  CHECK(close(f.entrance_pupil->diameter, u.entrance_pupil->diameter));
  CHECK(close(f.exit_pupil->z, u.exit_pupil->z));
  CHECK(close(f.exit_pupil->diameter, u.exit_pupil->diameter));
}

/// Determinant of a 3x3 matrix (cofactor expansion along the first row).
double det(const rtt::math::Mat3& m) {
  return m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)) -
         m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
         m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
}

/// Global point of the unfolded image-space point (0, 0, z).
Vec3 global_on_axis(const FirstOrder& f, double z) {
  return f.image_space.linear * Vec3(0.0, 0.0, z) + f.image_space.offset;
}

template <typename F>
std::string paraxial_error(F&& f) {
  try {
    f();
  } catch (const ParaxialError& e) {
    return e.what();
  }
  FAIL("no ParaxialError");
  return {};
}

}  // namespace

TEST_CASE("a plane fold mirror unfolds to the singlet", "[paraxial][fold][357]") {
  const CompiledSystem flat = unfolded();
  const FirstOrder u = rtt::paraxial::first_order(flat, PathId{0}, 0);
  REQUIRE(u.folds == 0);
  REQUIRE(u.image_space.linear == rtt::math::Mat3::Identity());
  REQUIRE(u.image_space.offset == Vec3::Zero());
  // content guard: a real lens with finite focus and pupils
  REQUIRE(u.efl);
  REQUIRE(*u.efl > 90.0);
  for (const double tilt : {45.0, 30.0, 60.0}) {
    INFO("tilt " << tilt << " deg");
    const CompiledSystem cs = one_fold(tilt);
    const FirstOrder f = rtt::paraxial::first_order(cs, PathId{0}, 0);
    REQUIRE(f.folds == 1);
    same_first_order(f, u);
    // the image space frame: a reflection (parity change, Greivenkamp p. 2-4) that maps the
    // unfolded axis behind the mirror onto the reflected axis
    CHECK(std::abs(det(f.image_space.linear) + 1.0) <= 1e-14);
    const double t = tilt * kPi / 180.0;
    const Vec3 d(0.0, std::sin(2.0 * t), -std::cos(2.0 * t));
    const Vec3 expected = Vec3(0.0, 0.0, kMirrorZ) + (*f.image_z - kMirrorZ) * d;
    CHECK((global_on_axis(f, *f.image_z) - expected).norm() <= 1e-12);
    // the detector of the file lies on the unfolded axis at its unfolded distance
    CHECK((global_on_axis(f, kImageZ) - cs.surfaces().back().to_global.translation()).norm() <=
          1e-12);
    // Seidel sums: the plane mirror adds nothing
    const auto sf = rtt::paraxial::seidel(cs, PathId{0}, 0);
    const auto su = rtt::paraxial::seidel(flat, PathId{0}, 0);
    CHECK(close(sf.sum.s1, su.sum.s1));
    CHECK(close(sf.sum.s2, su.sum.s2));
    CHECK(close(sf.sum.s3, su.sum.s3));
    CHECK(close(sf.sum.s4, su.sum.s4));
    CHECK(close(sf.sum.s5, su.sum.s5));
    REQUIRE(sf.surfaces.size() == su.surfaces.size() + 1);
    CHECK(sf.surfaces[3].terms.s1 == 0.0);  // the mirror (events: STO, L1.S1, L1.S2, M, IMG)
    CHECK(std::abs(su.sum.s1) > 1e-6);      // content guard: a real spherical aberration
  }
}

TEST_CASE("two fold mirrors keep the parity", "[paraxial][fold][357]") {
  const FirstOrder u = rtt::paraxial::first_order(unfolded(), PathId{0}, 0);
  const CompiledSystem cs = periscope();
  const FirstOrder f = rtt::paraxial::first_order(cs, PathId{0}, 0);
  REQUIRE(f.folds == 2);
  same_first_order(f, u);
  CHECK(std::abs(det(f.image_space.linear) - 1.0) <= 1e-14);
  // behind the periscope the axis runs along +z again, 10 mm higher and 10 mm shorter
  CHECK((global_on_axis(f, 100.0) - Vec3(0.0, 10.0, 90.0)).norm() <= 1e-12);
  // the traced rays: z along the unfolded axis, the same heights as without the mirrors
  const auto rf = rtt::paraxial::trace_ray(cs, PathId{0}, 0, -10.0, 1.0, 0.02);
  const auto ru = rtt::paraxial::trace_ray(unfolded(), PathId{0}, 0, -10.0, 1.0, 0.02);
  REQUIRE(rf.size() == ru.size() + 2);
  CHECK(close(rf.back().z, ru.back().z));
  CHECK(close(rf.back().y, ru.back().y));
  CHECK(close(rf.back().u, ru.back().u));
  CHECK(rf.back().n > 0.0);
}

TEST_CASE("what a fold does not cover stays a ParaxialError", "[paraxial][fold][357]") {
  const std::string image = pose(Vec3(0.0, kImageZ - kMirrorZ, kMirrorZ), 270.0);
  // a curved mirror tilted against the axis: off-axis, parabasal in v2
  CHECK(
      paraxial_error([&] {
        (void)rtt::paraxial::first_order(
            compiled(singlet(mirror("M", Vec3(0.0, 0.0, kMirrorZ), 45.0,
                                    R"(, "shape": {"base": {"type": "conic", "radius": -500.0}})"),
                             image, reflect("M"))),
            PathId{0}, 0);
      }).find("tilted against the z axis") != std::string::npos);
  // a plane mirror with a fourth-order term (A4 of an odd asphere with c = 0)
  CHECK(paraxial_error([&] {
          (void)rtt::paraxial::first_order(
              compiled(singlet(mirror("M", Vec3(0.0, 0.0, kMirrorZ), 45.0,
                                      R"(, "shape": {"base": {"type": "odd_asphere",
                                       "curvature": 0.0, "conic": 0.0,
                                       "coefficients": [0.0, 0.0, 0.0, 1e-6]}})"),
                               image, reflect("M"))),
              PathId{0}, 0);
        }).find("tilted against the z axis") != std::string::npos);
  // the stop on the fold mirror (aperture.stop, ADR 0039)
  CHECK(
      paraxial_error([&] {
        (void)rtt::paraxial::first_order(
            compiled(singlet(mirror("M", Vec3(0.0, 0.0, kMirrorZ), 45.0, {},
                                    R"(, "aperture": {"type": "circular", "radius": 20.0})"),
                             image, reflect("M"),
                             R"({"type": "epd", "value": 20.0, "stop": {"surface": "M"}})", false)),
            PathId{0}, 0);
      }).find("stop on a fold mirror") != std::string::npos);
  // per section (ADR 0054, point 3)
  const CompiledSystem cs = one_fold(45.0);
  CHECK(paraxial_error([&] {
          (void)rtt::paraxial::first_order_section(cs, PathId{0}, 0, rtt::paraxial::Section::y);
        }).find("per section") != std::string::npos);
}

TEST_CASE("the prescription of a folded path is that of the unfolded one",
          "[paraxial][fold][357]") {
  // PR B1 (ADR 0054, point 5): the rows are the unfolded ones plus the fold mirror, which leaves
  // the rays unchanged; z along the unfolded axis.
  const auto pf = rtt::paraxial::prescription(one_fold(45.0), PathId{0}, 0);
  const auto pu = rtt::paraxial::prescription(unfolded(), PathId{0}, 0);
  REQUIRE(pf.surfaces.size() == pu.surfaces.size() + 1);
  CHECK(close(pf.total_track, pu.total_track));
  CHECK(close(pf.paraxial_working_f_number, pu.paraxial_working_f_number));
  const auto same_row = [](const rtt::paraxial::PrescriptionSurface& a,
                           const rtt::paraxial::PrescriptionSurface& b) {
    CHECK(close(a.z, b.z));
    CHECK(a.n == b.n);
    REQUIRE(a.marginal);
    REQUIRE(b.marginal);
    CHECK(close(a.marginal->y, b.marginal->y));
    CHECK(close(a.marginal->u, b.marginal->u));
    REQUIRE(a.chief);
    REQUIRE(b.chief);
    CHECK(close(a.chief->y, b.chief->y));
    CHECK(close(a.chief->u, b.chief->u));
  };
  for (std::size_t k = 0; k < 3; ++k) same_row(pf.surfaces[k], pu.surfaces[k]);  // STO, L1
  same_row(pf.surfaces[4], pu.surfaces[3]);                                      // IMG
  // the fold row: at z = 49 on the unfolded axis, the marginal ray passes unchanged
  CHECK(close(pf.surfaces[3].z, kMirrorZ));
  REQUIRE(pf.surfaces[3].marginal);
  CHECK(pf.surfaces[3].marginal->u == pf.surfaces[2].marginal->u);
  REQUIRE(pf.surfaces[2].marginal->u != 0.0);  // content guard: a converging beam
}

TEST_CASE("a concave mirror at normal incidence behind a fold", "[paraxial][fold][357]") {
  // A plane fold at z = 20 to +y, then a concave mirror R = -100 mm 30 mm further at normal
  // incidence, which sends the light back (n changes sign, ADR 0054, point 1); unfolded it is the
  // same mirror 50 mm behind the stop with the image at z = 2.
  const std::string folded =
      R"({"schema_version": "0.11.0", "name": "x", "units": {"length": "mm", "wavelength": "um"}, "wavelengths": [{"um": 0.5876, "reference": true}], "aperture": {"type": "epd", "value": 10.0}, "fields": {"points": [{}, {"y": 1.0}]}, "root": {"type": "assembly", "name": "system", "children": [{"type": "stop", "name": "stop", "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 5.0}}]}, {"type": "mirror", "name": "M1", "pose": {"position": [0.0, 0.0, 20.0], "rotation_deg": [45.0, 0.0, 0.0]}, "surfaces": [{"id": "M1", "aperture": {"type": "circular", "radius": 10.0}}]}, {"type": "mirror", "name": "M2", "pose": {"position": [0.0, 30.0, 20.0], "rotation_deg": [-90.0, 0.0, 0.0]}, "surfaces": [{"id": "M2", "shape": {"base": {"type": "conic", "radius": -100.0}}, "aperture": {"type": "circular", "radius": 8.0}}]}, {"type": "detector", "name": "image", "pose": {"position": [0.0, -18.0, 20.0], "rotation_deg": [90.0, 0.0, 0.0]}, "surfaces": [{"id": "IMG"}]}]}, "paths": [{"name": "main", "events": [{"surface": "STO", "kind": "transmit"}, {"surface": "M1", "kind": "reflect"}, {"surface": "M2", "kind": "reflect"}, {"surface": "IMG", "kind": "transmit"}]}]})";
  const std::string flat =
      R"({"schema_version": "0.11.0", "name": "x", "units": {"length": "mm", "wavelength": "um"}, "wavelengths": [{"um": 0.5876, "reference": true}], "aperture": {"type": "epd", "value": 10.0}, "fields": {"points": [{}, {"y": 1.0}]}, "root": {"type": "assembly", "name": "system", "children": [{"type": "stop", "name": "stop", "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 5.0}}]}, {"type": "mirror", "name": "M2", "pose": {"position": [0.0, 0.0, 50.0]}, "surfaces": [{"id": "M2", "shape": {"base": {"type": "conic", "radius": -100.0}}, "aperture": {"type": "circular", "radius": 8.0}}]}, {"type": "detector", "name": "image", "pose": {"position": [0.0, 0.0, 2.0], "rotation_deg": [180.0, 0.0, 0.0]}, "surfaces": [{"id": "IMG"}]}]}, "paths": [{"name": "main", "events": [{"surface": "STO", "kind": "transmit"}, {"surface": "M2", "kind": "reflect"}, {"surface": "IMG", "kind": "transmit"}]}]})";
  const CompiledSystem cf = compiled(folded);
  const CompiledSystem cu = compiled(flat);
  const FirstOrder f = rtt::paraxial::first_order(cf, PathId{0}, 0);
  const FirstOrder u = rtt::paraxial::first_order(cu, PathId{0}, 0);
  REQUIRE(f.folds == 1);
  REQUIRE(u.image_direction == -1);  // content guard: light leaves towards -z
  same_first_order(f, u);
  REQUIRE(f.efl);
  CHECK(std::abs(*f.efl - 50.0) < 0.1);  // |R| / 2, times n of AIR (absolute indices)
  const auto sf = rtt::paraxial::seidel(cf, PathId{0}, 0);
  const auto su = rtt::paraxial::seidel(cu, PathId{0}, 0);
  CHECK(close(sf.sum.s1, su.sum.s1));
  CHECK(close(sf.sum.s3, su.sum.s3));
  CHECK(std::abs(su.sum.s3) > 1e-9);  // content guard: off-axis field, astigmatism
}
