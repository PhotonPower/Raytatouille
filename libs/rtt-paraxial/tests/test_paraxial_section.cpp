// First-order data per section x-z and y-z (#315, ADR 0042, point 2; acceptance 1 (paraxial
// part), 6 and 7 (paraxial part)): first_order_section, trace_ray_section and is_anamorphic.
//
// Oracle: an independent 2x2 matrix computation in this file (Expected, expected()), built from
// the curvatures per section that the hand computation gives for each test system, not from the
// code under test. Elementary matrices of the y-nu trace in reduced form (y, nu = n u):
// transfer [[1, t/n], [0, 1]], refraction [[1, 0], [-phi, 1]] with phi = (n' - n) c
// (Greivenkamp, OPTI-201/202 lecture notes, Sec. 9, p. 9-2, docs/quellen.md, row "Paraxialer
// y-nu-Trace"); an ideal lens adds phi = |n| / f (ADR 0031, point 6; docs/quellen.md, row
// "Ideale Linse paraxial"). Cardinal points from the matrix [[A, B], [C, D]] between the first
// and last vertex: Phi = -C, EFL = 1/Phi, BFL = A n'/Phi (Greivenkamp, Sec. 9, p. 9-12:
// phi = -n'u'_k/y_1, BFD = -y_k/u'_k). Entrance pupil = paraxial image of the stop through the
// part in front of it: with [[A, B], [C, D]] from an object-space plane z_0 to the stop plane,
// z_ep = z_0 + B n/A and r_ep = r_stop/|A| (docs/quellen.md, row "Pupillen, Winkel- und
// Lateralvergrößerung"; the same derivation as in paraxial.cpp, here written as matrix
// products).
//
// Curvature per section (ADR 0042, point 2): the surface lies in the global frame with
// p_g = o + M p_l, M = diag(Q, sigma); the global Hessian is H_g = sigma Q H_l Q^T, so a
// rotation by 90 degree about z swaps the sections and a local axis -z (sigma = -1, tilt by
// 180 degree about x, Q = diag(1, -1)) reverses the sign of both curvatures. The ideal
// cylinder lens has the power tensor phi b_g b_g^T with b_g = Q b_l, b_l = (-sin psi, cos psi).
//
// Tolerances: relative 1e-13 for a handful of y-nu steps that round by a few u = 1.1e-16 each.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <optional>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/paraxial/paraxial.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::paraxial::FirstOrder;
using rtt::paraxial::ParaxialError;
using rtt::paraxial::RayAtEvent;
using rtt::paraxial::Section;

namespace {

constexpr double kN = 1.5168;  // CONST:1.5168
constexpr double kRel = 1e-13;

// ---------------------------------------------------------------- test systems (JSON) -----

/// "%.17g" reads back as the same double.
std::string json_number(double v) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.17g", v);
  return buf;
}

std::string pose(double z, double rx = 0.0, double ry = 0.0, double rz = 0.0) {
  return R"({"position": [0.0, 0.0, )" + json_number(z) + R"(], "rotation_deg": [)" +
         json_number(rx) + ", " + json_number(ry) + ", " + json_number(rz) + "]}";
}

std::string stop(double z, double radius) {
  return R"({"type": "stop", "name": "stop", "pose": )" + pose(z) +
         R"(, "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": )" +
         json_number(radius) + "}}]}";
}

/// Lens `name` at global z with thickness t: front surface `<name>.S1` with `shape` (JSON of
/// shape.base) and the pose `front_pose` relative to the element, rear surface plane.
std::string lens(const std::string& name,
                 const std::string& element_pose,
                 const std::string& shape,
                 const std::string& front_pose,
                 double t) {
  return R"({"type": "lens", "name": ")" + name + R"(", "pose": )" + element_pose +
         R"(, "material": "CONST:1.5168", "surfaces": [{"id": ")" + name + R"(.S1", "pose": )" +
         front_pose + R"(, "shape": {"base": )" + shape +
         R"(}, "aperture": {"type": "rectangular", "half_width_x": 12.0, "half_width_y": 12.0}},)" +
         R"({"id": ")" + name + R"(.S2", "pose": )" + pose(t) +
         R"(, "aperture": {"type": "rectangular", "half_width_x": 12.0, "half_width_y": 12.0}}]})";
}

std::string thin(const std::string& element_pose, const std::string& interaction) {
  return R"({"type": "thin_element", "name": "lens", "pose": )" + element_pose +
         R"(, "surfaces": [{"id": "IL", "aperture": {"type": "circular", "radius": 12.0}, )" +
         R"("interaction": )" + interaction + "}]}";
}

std::string detector(double z) {
  return R"({"type": "detector", "name": "image", "pose": )" + pose(z) +
         R"(, "surfaces": [{"id": "IMG"}]})";
}

/// System in VACUUM at 0.5876 um, object at infinity, path "main" = "auto".
CompiledSystem build(const std::vector<std::string>& children,
                     const std::string& aperture = R"({"type": "epd", "value": 10.0})") {
  std::string list;
  for (std::size_t i = 0; i < children.size(); ++i) {
    if (i > 0) list += ", ";
    list += children[i];
  }
  const std::string json =
      R"({"schema_version": "0.9.0", "name": "section test",)"
      R"( "units": {"length": "mm", "wavelength": "um"}, "environment": {"medium": "VACUUM"},)"
      R"( "wavelengths": [{"um": 0.5876, "reference": true}], "aperture": )" +
      aperture + R"(, "fields": {"points": [{}]}, "root": {"type": "assembly", "name": "system",)" +
      R"( "children": [)" + list + R"(]}, "paths": [{"name": "main", "events": "auto"}]})";
  return rtt::compile::compile(rtt::io::parse_system(json), rtt::material::MaterialLibrary{});
}

CompiledSystem load_reference(const std::string& file) {
  return rtt::compile::compile(rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + file),
                               rtt::material::MaterialLibrary{});
}

// ------------------------------------------------- independent matrix computation -----

struct Abcd {
  double a = 1.0, b = 0.0, c = 0.0, d = 1.0;
};

/// l * r (r acts first).
Abcd operator*(const Abcd& l, const Abcd& r) {
  return {l.a * r.a + l.b * r.c, l.a * r.b + l.b * r.d, l.c * r.a + l.d * r.c,
          l.c * r.b + l.d * r.d};
}

Abcd transfer(double t, double n) {
  return {1.0, t / n, 0.0, 1.0};
}
Abcd refraction(double phi) {
  return {1.0, 0.0, -phi, 1.0};
}

/// One event of a path as the hand computation sees it in one section.
struct Ev {
  double z = 0.0;         ///< global z of the vertex, mm
  double c = 0.0;         ///< curvature in this section, global orientation, 1/mm
  double n_after = 1.0;   ///< index after the event (all positive: no mirrors here)
  double thin_phi = 0.0;  ///< power of an ideal (cylinder) lens in this section, 1/mm
};

struct Expected {
  std::optional<double> efl;
  std::optional<double> bfl;
  std::optional<double> ffl;                ///< first vertex -> F, positive upstream
  std::optional<double> front_principal_z;  ///< global z of H
  std::optional<double> rear_principal_z;   ///< global z of H'
  std::optional<double> z_ep;
  std::optional<double> r_ep;
};

/// Matrix from the plane z_in (index n_in) through events [first, last) to the plane z_out.
Abcd matrix(const std::vector<Ev>& ev,
            std::size_t first,
            std::size_t last,
            double z_in,
            double n_in,
            double z_out) {
  Abcd m;
  double z = z_in;
  double n = n_in;
  for (std::size_t i = first; i < last; ++i) {
    m = transfer(ev[i].z - z, n) * m;
    m = refraction(ev[i].c * (ev[i].n_after - n) + n * ev[i].thin_phi) * m;
    z = ev[i].z;
    n = ev[i].n_after;
  }
  return transfer(z_out - z, n) * m;
}

/// efl, bfl and entrance pupil of `ev` (object space n = 1, object at infinity). `stop` is the
/// index of the stop event (a Stop element: it only transmits); `r_stop` its radius for the
/// system aperture stop_size, none for an EPD (then r_ep = epd / 2).
Expected expected(const std::vector<Ev>& ev,
                  std::size_t stop,
                  std::optional<double> r_stop,
                  double epd = 10.0) {
  Expected e;
  // Vertices: the first and last events that change the ray (index step or ideal lens).
  std::optional<std::size_t> v1;
  std::size_t vk = 0;
  double n = 1.0;
  for (std::size_t i = 0; i < ev.size(); ++i) {
    if (ev[i].n_after != n || ev[i].thin_phi != 0.0) {
      if (!v1) v1 = i;
      vk = i;
    }
    n = ev[i].n_after;
  }
  // Without a vertex nothing changes the ray: afocal (the afocal section of an ideal cylinder
  // lens in air).
  if (v1) {
    const Abcd m = matrix(ev, *v1, vk + 1, ev[*v1].z, 1.0, ev[vk].z);
    const double phi = -m.c;
    if (std::abs(phi) > 1e-14) {
      // Front side (Greivenkamp, Sec. 9, p. 9-14: f_F = -n/phi, FFD = -y_1/u_1): a ray from F
      // with slope u reaches V1 with (y, nu) = ((z_V1 - z_F) u, u) and leaves parallel:
      // C (z_V1 - z_F) + D = 0, so ffl = z_V1 - z_F = D/phi (n = 1 in object space).
      // Principal planes at the focal lengths n/phi and n'/phi from F and F' (p. 9-12, 9-14).
      const double n_img = ev[vk].n_after;
      e.efl = 1.0 / phi;
      e.bfl = m.a * n_img / phi;
      e.ffl = m.d / phi;
      e.front_principal_z = ev[*v1].z - *e.ffl + 1.0 / phi;
      e.rear_principal_z = ev[vk].z + *e.bfl - n_img / phi;
    }
  }
  const Abcd front = matrix(ev, 0, stop, ev.front().z, 1.0, ev[stop].z);
  e.z_ep = ev.front().z + front.b * 1.0 / front.a;
  e.r_ep = r_stop ? *r_stop / std::abs(front.a) : epd / 2.0;
  return e;
}

bool close(double a, double b, double rel = kRel) {
  return std::abs(a - b) <= rel * std::max(std::abs(b), 1.0);
}

void check(const FirstOrder& fo, const Expected& e) {
  REQUIRE(fo.efl.has_value() == e.efl.has_value());
  if (e.efl) {
    INFO("efl " << *fo.efl << " expected " << *e.efl);
    CHECK(close(*fo.efl, *e.efl));
    REQUIRE(fo.bfl.has_value());
    INFO("bfl " << *fo.bfl << " expected " << *e.bfl);
    CHECK(close(*fo.bfl, *e.bfl));
    REQUIRE(fo.ffl.has_value());
    INFO("ffl " << *fo.ffl << " expected " << *e.ffl);
    CHECK(close(*fo.ffl, *e.ffl));
    CHECK(close(fo.front_principal_z.value(), *e.front_principal_z));
    CHECK(close(fo.rear_principal_z.value(), *e.rear_principal_z));
  } else {
    CHECK_FALSE(fo.bfl.has_value());
    CHECK_FALSE(fo.ffl.has_value());
  }
  REQUIRE(fo.entrance_pupil.has_value());
  REQUIRE(fo.entrance_pupil->z.has_value());
  INFO("z_ep " << *fo.entrance_pupil->z << " expected " << *e.z_ep);
  CHECK(close(*fo.entrance_pupil->z, *e.z_ep));
  REQUIRE(fo.entrance_pupil->diameter.has_value());
  CHECK(close(*fo.entrance_pupil->diameter, 2.0 * *e.r_ep));
}

FirstOrder section(const CompiledSystem& cs, Section s, PathId path = PathId{0}) {
  return rtt::paraxial::first_order_section(cs, path, cs.reference_wavelength(), s);
}

/// Runs f and returns the ParaxialError it throws.
template <typename F>
ParaxialError paraxial_error(F&& f) {
  try {
    f();
  } catch (const ParaxialError& e) {
    return e;
  }
  FAIL("no ParaxialError");
  return ParaxialError("unreachable");
}

void same(const std::optional<rtt::paraxial::Pupil>& a,
          const std::optional<rtt::paraxial::Pupil>& b) {
  REQUIRE(a.has_value() == b.has_value());
  if (!a || !b) return;
  CHECK(a->z == b->z);
  CHECK(a->diameter == b->diameter);
}

/// Bitwise equality of all fields (std::optional<double> compares the values with ==; a
/// difference in the sign of a zero would pass, the CHECKs on signbit below catch it).
void same(const FirstOrder& a, const FirstOrder& b) {
  CHECK(a.object_index == b.object_index);
  CHECK(a.image_index == b.image_index);
  CHECK(a.image_direction == b.image_direction);
  CHECK(a.power == b.power);
  CHECK(std::signbit(a.power) == std::signbit(b.power));
  CHECK(a.efl == b.efl);
  CHECK(a.front_focal_length == b.front_focal_length);
  CHECK(a.rear_focal_length == b.rear_focal_length);
  CHECK(a.ffl == b.ffl);
  CHECK(a.bfl == b.bfl);
  CHECK(a.front_focal_z == b.front_focal_z);
  CHECK(a.rear_focal_z == b.rear_focal_z);
  CHECK(a.front_principal_z == b.front_principal_z);
  CHECK(a.rear_principal_z == b.rear_principal_z);
  CHECK(a.image_z == b.image_z);
  CHECK(a.lateral_magnification == b.lateral_magnification);
  CHECK(a.angular_magnification == b.angular_magnification);
  same(a.entrance_pupil, b.entrance_pupil);
  same(a.exit_pupil, b.exit_pupil);
}

// Curvatures of the test surfaces.
constexpr double kRy = 51.68;  // f = R / (n - 1) = 100 mm
constexpr double kRx = 77.52;  // f = 150 mm

const std::string kBiconic =
    R"({"type": "biconic", "radius_x": 77.52, "radius_y": 51.68, "conic_x": -0.5})";
const std::string kCylinder = R"({"type": "cylinder", "radius": 51.68})";

/// Events of the singlet "stop at 0, lens at 5 (front curvature c, t = 4), image at z_img".
std::vector<Ev> singlet_events(double c, double z_img = 110.0) {
  return {{0.0, 0.0, 1.0}, {5.0, c, kN}, {9.0, 0.0, 1.0}, {z_img, 0.0, 1.0}};
}

}  // namespace

TEST_CASE("plano-convex cylinder lens in VACUUM: y section R/(n-1), x section afocal",
          "[paraxial][section]") {
  // ADR 0042, acceptance 1 (paraxial part): curved side R = 51.68 in the y section in front,
  // plane rear side, t = 4, n = 1.5168. Two components (Greivenkamp, OPTI-201/202, Sec. 7,
  // p. 7-4, docs/quellen.md row "Zwei Komponenten"): P = P1 + P2 - (t/n) P1 P2 with P2 = 0, so
  // efl = R/(n - 1) = 100 mm independent of t; H' lies t/n before the rear vertex
  // (Greivenkamp, Sec. 9, p. 9-12), so bfl = efl - t/n.
  const CompiledSystem cs =
      build({stop(0.0, 10.0), lens("L1", pose(5.0), kCylinder, pose(0.0), 4.0), detector(110.0)});
  const FirstOrder y = section(cs, Section::y);
  REQUIRE(y.efl.has_value());
  CHECK(close(*y.efl, kRy / (kN - 1.0)));
  CHECK(close(*y.efl, 100.0, 1e-12));
  CHECK(close(y.bfl.value(), 100.0 - 4.0 / kN, 1e-12));
  const FirstOrder x = section(cs, Section::x);
  CHECK_FALSE(x.efl.has_value());
  CHECK_FALSE(x.bfl.has_value());
  CHECK_FALSE(x.rear_focal_z.has_value());
  CHECK(x.power == 0.0);
  // Both sections against the matrix computation, entrance pupil = the stop at z = 0.
  check(y, expected(singlet_events(1.0 / kRy), 0, std::nullopt));
  check(x, expected(singlet_events(0.0), 0, std::nullopt));
}

TEST_CASE("cylinder lens rotated by 90 degree about z: power in x, y section afocal",
          "[paraxial][section]") {
  // H_g = Q H_l Q^T with Q the rotation by 90 degree: (0, c) -> (c, 0). cos(90 degree) rounds to
  // 6.1e-17, so the y section gets the curvature 3.7e-33 c from the rounding of Q alone; it must
  // still count as afocal (ADR 0042, point 2: a section without power is afocal).
  const CompiledSystem cs =
      build({stop(0.0, 10.0), lens("L1", pose(5.0, 0.0, 0.0, 90.0), kCylinder, pose(0.0), 4.0),
             detector(110.0)});
  check(section(cs, Section::x), expected(singlet_events(1.0 / kRy), 0, std::nullopt));
  const FirstOrder y = section(cs, Section::y);
  CHECK_FALSE(y.efl.has_value());
  CHECK_FALSE(y.bfl.has_value());
  check(y, expected(singlet_events(0.0), 0, std::nullopt));
}

TEST_CASE("crossed cylinders (90 degree) against the matrix computation", "[paraxial][section]") {
  // L1: cylinder R = 51.68 at z = 5 (power in y); L2: cylinder R = 38.76 at z = 30, rotated by
  // 90 degree about z (power in x). Stop at z = 0.
  const CompiledSystem cs =
      build({stop(0.0, 10.0), lens("L1", pose(5.0), kCylinder, pose(0.0), 4.0),
             lens("L2", pose(30.0, 0.0, 0.0, 90.0), R"({"type": "cylinder", "radius": 38.76})",
                  pose(0.0), 4.0),
             detector(110.0)});
  const auto events = [](double c1, double c2) {
    return std::vector<Ev>{{0.0, 0.0, 1.0}, {5.0, c1, kN},    {9.0, 0.0, 1.0},
                           {30.0, c2, kN},  {34.0, 0.0, 1.0}, {110.0, 0.0, 1.0}};
  };
  const Expected ex = expected(events(0.0, 1.0 / 38.76), 0, std::nullopt);
  const Expected ey = expected(events(1.0 / kRy, 0.0), 0, std::nullopt);
  // Content guard: the sections differ (f_x = 75 at z = 30, f_y = 100 at z = 5).
  REQUIRE(std::abs(*ex.efl - *ey.efl) > 1.0);
  REQUIRE(std::abs(*ex.bfl - *ey.bfl) > 1.0);
  check(section(cs, Section::x), ex);
  check(section(cs, Section::y), ey);
}

TEST_CASE("biconic singlet against the matrix computation", "[paraxial][section]") {
  // Front biconic R_x = 77.52, R_y = 51.68, k_x = -0.5 (the conic does not enter first order).
  const CompiledSystem cs =
      build({stop(0.0, 10.0), lens("L1", pose(5.0), kBiconic, pose(0.0), 4.0), detector(110.0)});
  const FirstOrder x = section(cs, Section::x);
  const FirstOrder y = section(cs, Section::y);
  check(x, expected(singlet_events(1.0 / kRx), 0, std::nullopt));
  check(y, expected(singlet_events(1.0 / kRy), 0, std::nullopt));
  CHECK(close(x.efl.value(), 150.0, 1e-12));
  CHECK(close(y.efl.value(), 100.0, 1e-12));
}

TEST_CASE("toroid in front of the stop: entrance pupil per section", "[paraxial][section]") {
  // Toroid: y profile R = 51.68, radius of rotation 77.52, so c_y = 1/51.68 and c_x = 1/77.52
  // (ADR 0042, point 1, table). Lens at z = 5 (t = 4), stop at z = 20 with radius 4 and the
  // system aperture stop_size: the entrance pupil is the image of the stop through the lens and
  // differs per section in position and size.
  const CompiledSystem cs = build(
      {lens("L1", pose(5.0), R"({"type": "toroid", "radius": 51.68, "radius_of_rotation": 77.52})",
            pose(0.0), 4.0),
       stop(20.0, 4.0), detector(110.0)},
      R"({"type": "stop_size"})");
  const auto events = [](double c) {
    return std::vector<Ev>{{5.0, c, kN}, {9.0, 0.0, 1.0}, {20.0, 0.0, 1.0}, {110.0, 0.0, 1.0}};
  };
  const Expected ex = expected(events(1.0 / kRx), 2, 4.0);
  const Expected ey = expected(events(1.0 / kRy), 2, 4.0);
  // Content guard: z_ep,x != z_ep,y (ADR 0042, acceptance 6) and different pupil sizes.
  REQUIRE(std::abs(*ex.z_ep - *ey.z_ep) > 0.1);
  REQUIRE(std::abs(*ex.r_ep - *ey.r_ep) > 0.01);
  const FirstOrder x = section(cs, Section::x);
  const FirstOrder y = section(cs, Section::y);
  check(x, ex);
  check(y, ey);
  CHECK(*x.entrance_pupil->z != *y.entrance_pupil->z);
}

TEST_CASE("ideal cylinder lens: psi = 0, 90 degree and psi = 0 on rotated surfaces",
          "[paraxial][section]") {
  // f = 50 at z = 0 in VACUUM, stop at z = -10 (radius 5, EPD 10). phi = |n|/f along
  // b_g = Q b_l, b_l = (-sin psi, cos psi) (ADR 0042, point 2; ADR 0031, point 6): section s
  // has phi b_g,s^2. The focal section has efl = bfl = f; the other is afocal.
  struct Case {
    const char* name;
    std::string element_pose;
    double psi_deg;
    Section focal;
  };
  const std::vector<Case> cases = {
      {"psi = 0", pose(0.0), 0.0, Section::y},
      {"psi = 90", pose(0.0), 90.0, Section::x},
      {"psi = 0, surface rotated 90 about z", pose(0.0, 0.0, 0.0, 90.0), 0.0, Section::x},
      // sigma = -1: no sign enters (b b^T is even in b, phi counts along the propagation).
      {"psi = 0, surface tilted 180 about x", pose(0.0, 180.0, 0.0, 0.0), 0.0, Section::y},
  };
  for (const Case& c : cases) {
    INFO(c.name);
    const CompiledSystem cs =
        build({stop(-10.0, 5.0),
               thin(c.element_pose,
                    R"({"type": "ideal_cylinder_lens", "focal_length": 50.0, "axis_deg": )" +
                        json_number(c.psi_deg) + "}"),
               detector(60.0)});
    const auto events = [](double phi) {
      return std::vector<Ev>{{-10.0, 0.0, 1.0}, {0.0, 0.0, 1.0, phi}, {60.0, 0.0, 1.0}};
    };
    const Section other = c.focal == Section::x ? Section::y : Section::x;
    const FirstOrder focal = section(cs, c.focal);
    const FirstOrder afocal = section(cs, other);
    check(focal, expected(events(1.0 / 50.0), 0, std::nullopt));
    CHECK(close(focal.efl.value(), 50.0));
    CHECK(close(focal.bfl.value(), 50.0));
    check(afocal, expected(events(0.0), 0, std::nullopt));
    CHECK_FALSE(afocal.efl.has_value());
  }
}

TEST_CASE("sign of the transform: biconic rotated by 90 degree about z and with axis -z",
          "[paraxial][section]") {
  SECTION("rotated by 90 degree about z: the sections swap") {
    // Q = [[0, -1], [1, 0]]: H_g = Q diag(c_x, c_y) Q^T = diag(c_y, c_x).
    const CompiledSystem cs =
        build({stop(0.0, 10.0), lens("L1", pose(5.0), kBiconic, pose(0.0, 0.0, 0.0, 90.0), 4.0),
               detector(110.0)});
    check(section(cs, Section::x), expected(singlet_events(1.0 / kRy), 0, std::nullopt));
    check(section(cs, Section::y), expected(singlet_events(1.0 / kRx), 0, std::nullopt));
  }
  SECTION("tilted by 180 degree about x (sigma = -1): both curvatures change sign") {
    // Q = diag(1, -1), sigma = -1: H_g = -diag(c_x, c_y). The biconic with R_x = -77.52 and
    // R_y = -51.68 then acts as the convex front of the biconic singlet above.
    const CompiledSystem cs = build(
        {stop(0.0, 10.0),
         lens("L1", pose(5.0),
              R"({"type": "biconic", "radius_x": -77.52, "radius_y": -51.68, "conic_x": -0.5})",
              pose(0.0, 180.0, 0.0, 0.0), 4.0),
         detector(110.0)});
    check(section(cs, Section::x), expected(singlet_events(1.0 / kRx), 0, std::nullopt));
    check(section(cs, Section::y), expected(singlet_events(1.0 / kRy), 0, std::nullopt));
    CHECK(close(section(cs, Section::x).efl.value(), 150.0, 1e-12));
  }
}

TEST_CASE("trace_ray_section against the y-nu trace by hand", "[paraxial][section]") {
  // Plano-convex cylinder lens (power in y). Ray from z = -1 with height 2 and slope 0.01 in
  // each section; by hand: transfer y' = y + (z' - z) nu/n, refraction nu' = nu - y c (n' - n)
  // (Greivenkamp, Sec. 9, p. 9-2).
  const CompiledSystem cs =
      build({stop(0.0, 10.0), lens("L1", pose(5.0), kCylinder, pose(0.0), 4.0), detector(110.0)});
  for (const Section s : {Section::x, Section::y}) {
    INFO((s == Section::x ? "x" : "y"));
    const std::vector<Ev> ev = singlet_events(s == Section::y ? 1.0 / kRy : 0.0);
    const std::vector<RayAtEvent> rays =
        rtt::paraxial::trace_ray_section(cs, PathId{0}, 0, -1.0, 2.0, 0.01, s);
    REQUIRE(rays.size() == ev.size());
    double z = -1.0;
    double y = 2.0;
    double n = 1.0;
    double nu = 0.01;
    for (std::size_t i = 0; i < ev.size(); ++i) {
      INFO("event " << i);
      y += (ev[i].z - z) * nu / n;
      nu -= y * ev[i].c * (ev[i].n_after - n);
      z = ev[i].z;
      n = ev[i].n_after;
      CHECK(rays[i].z == ev[i].z);
      CHECK(close(rays[i].y, y));
      CHECK(close(rays[i].u, nu / n));
      CHECK(rays[i].n == n);
    }
    // Content guard: the y section bends the ray, the x section does not.
    if (s == Section::y) {
      CHECK(rays.back().u < 0.0);
    } else {
      CHECK(close(rays.back().u, 0.01));
    }
  }
}

TEST_CASE("rotationally symmetric paths: both sections bit for bit as first_order",
          "[paraxial][section]") {
  // ADR 0042, point 2: rotationally symmetric surfaces take c exactly as first_order, without
  // the rounding of Q Q^T, so every field of FirstOrder is bitwise equal.
  struct Ref {
    const char* file;
    const char* path;
  };
  for (const Ref& r :
       {Ref{"m1/singlet_const.rtt.json", "main"},
        Ref{"m1/two_lenses_stop_between.rtt.json", nullptr},
        Ref{"m1/paraboloid_mirror.rtt.json", nullptr},
        Ref{"m1/telecentric_singlet.rtt.json", nullptr},
        Ref{"m2/paraboloid_stop.rtt.json", nullptr}, Ref{"m2/telecentric_4f.rtt.json", nullptr},
        Ref{"r2/ideal_lens.rtt.json", "ideal lens"}}) {
    INFO(r.file);
    const CompiledSystem cs = load_reference(r.file);
    const PathId path = r.path ? *cs.find_path(r.path) : PathId{0};
    // Content guard: the path is focal except for the afocal 4f system.
    const FirstOrder ref = rtt::paraxial::first_order(cs, path, cs.reference_wavelength());
    REQUIRE(ref.efl.has_value() == (std::string(r.file) != "m2/telecentric_4f.rtt.json"));
    CHECK_FALSE(rtt::paraxial::is_anamorphic(cs, path));
    for (std::uint16_t w = 0; w < cs.wavelengths_um().size(); ++w) {
      INFO("wavelength " << w);
      const FirstOrder fo = rtt::paraxial::first_order(cs, path, w);
      for (const Section s : {Section::x, Section::y}) {
        same(rtt::paraxial::first_order_section(cs, path, w, s), fo);
        const std::vector<RayAtEvent> a = rtt::paraxial::trace_ray(cs, path, w, -3.0, 1.5, 0.02);
        const std::vector<RayAtEvent> b =
            rtt::paraxial::trace_ray_section(cs, path, w, -3.0, 1.5, 0.02, s);
        REQUIRE(a.size() == b.size());
        for (std::size_t i = 0; i < a.size(); ++i) {
          CHECK(a[i].z == b[i].z);
          CHECK(a[i].y == b[i].y);
          CHECK(a[i].u == b[i].u);
          CHECK(a[i].n == b[i].n);
        }
      }
    }
  }
}

TEST_CASE("principal axes not aligned: ParaxialError at the surface", "[paraxial][section]") {
  // ADR 0042, acceptance 7: a biconic and an ideal cylinder lens rotated by 45 degree have
  // |H_g,xy| (resp. |phi b_x b_y|) far above kSymmetryTolerance.
  const auto check_error = [](const CompiledSystem& cs, const std::string& id,
                              const std::string& pointer) {
    for (const Section s : {Section::x, Section::y}) {
      for (const ParaxialError& e :
           {paraxial_error([&] { (void)rtt::paraxial::first_order_section(cs, PathId{0}, 0, s); }),
            paraxial_error([&] {
              (void)rtt::paraxial::trace_ray_section(cs, PathId{0}, 0, -20.0, 1.0, 0.0, s);
            })}) {
        INFO(e.what());
        CHECK(std::string(e.what()).find("principal axes not aligned with x and y") !=
              std::string::npos);
        REQUIRE(e.surface() == rtt::model::SurfaceId(id));
        CHECK(e.location() == pointer);
      }
    }
  };
  check_error(
      build({stop(0.0, 10.0), lens("L1", pose(5.0), kBiconic, pose(0.0, 0.0, 0.0, 45.0), 4.0),
             detector(110.0)}),
      "L1.S1", "/root/children/1/surfaces/0");
  check_error(build({stop(-10.0, 5.0),
                     thin(pose(0.0), R"({"type": "ideal_cylinder_lens", "focal_length": 50.0,)"
                                     R"( "axis_deg": 45.0})"),
                     detector(60.0)}),
              "IL", "/root/children/1/surfaces/0");
  check_error(build({stop(-10.0, 5.0),
                     thin(pose(0.0, 0.0, 0.0, 45.0),
                          R"({"type": "ideal_cylinder_lens", "focal_length": 50.0})"),
                     detector(60.0)}),
              "IL", "/root/children/1/surfaces/0");
}

TEST_CASE("first_order and trace_ray still reject anamorphic paths", "[paraxial][section]") {
  // ADR 0042, point 2 ("Unverändert"): same messages as before.
  const CompiledSystem biconic =
      build({stop(0.0, 10.0), lens("L1", pose(5.0), kBiconic, pose(0.0), 4.0), detector(110.0)});
  for (const ParaxialError& e :
       {paraxial_error([&] { (void)rtt::paraxial::first_order(biconic, PathId{0}, 0); }),
        paraxial_error([&] { (void)rtt::paraxial::trace_ray(biconic, PathId{0}, 0, -1, 1, 0); })}) {
    INFO(e.what());
    CHECK(std::string(e.what()).find("not rotationally symmetric in second order") !=
          std::string::npos);
    CHECK(e.location() == "/root/children/1/surfaces/0");
  }
  const CompiledSystem ideal = load_reference("r2/ideal_lens.rtt.json");
  const PathId cyl = *ideal.find_path("cylinder lens");
  const ParaxialError e = paraxial_error([&] { (void)rtt::paraxial::first_order(ideal, cyl, 0); });
  CHECK(std::string(e.what()).find("the ideal cylinder lens is not rotationally symmetric") !=
        std::string::npos);
  // The section API takes it (psi = 0: power in y, x afocal).
  CHECK(rtt::paraxial::first_order_section(ideal, cyl, 0, Section::y).efl.has_value());
  CHECK_FALSE(rtt::paraxial::first_order_section(ideal, cyl, 0, Section::x).efl.has_value());
}

TEST_CASE("is_anamorphic", "[paraxial][section]") {
  const CompiledSystem ideal = load_reference("r2/ideal_lens.rtt.json");
  CHECK_FALSE(rtt::paraxial::is_anamorphic(ideal, *ideal.find_path("ideal lens")));
  CHECK(rtt::paraxial::is_anamorphic(ideal, *ideal.find_path("cylinder lens")));
  CHECK(rtt::paraxial::is_anamorphic(
      build({stop(0.0, 10.0), lens("L1", pose(5.0), kBiconic, pose(0.0), 4.0), detector(110.0)}),
      PathId{0}));
  CHECK(rtt::paraxial::is_anamorphic(
      build({stop(0.0, 10.0), lens("L1", pose(5.0), kCylinder, pose(0.0), 4.0), detector(110.0)}),
      PathId{0}));
  // Equal sections: the Hessian is (c, 0, c), not anamorphic.
  CHECK_FALSE(rtt::paraxial::is_anamorphic(
      build({stop(0.0, 10.0),
             lens("L1", pose(5.0), R"({"type": "biconic", "radius_x": 51.68, "radius_y": 51.68})",
                  pose(0.0), 4.0),
             detector(110.0)}),
      PathId{0}));
  CHECK_THROWS_AS(rtt::paraxial::is_anamorphic(ideal, PathId{7}), ParaxialError);
}

TEST_CASE("ideal cylinder lens in its section without power is no vertex", "[paraxial][section]") {
  // Review of C1 (P2-2): for psi = 90 degree, b_l = (-1, cos 90 degree) = (-1, 6.1e-17), and
  // psi = 0 on a surface rotated by 90 degree gives b_g,y of the same size, so b_y^2 = 3.7e-33.
  // That factor is rounding, not power: in the y section the lens must not count as a vertex
  // that changes the ray, or bfl, ffl and the principal planes would be measured from it. The
  // oracle (expected()) has no power there (thin_phi = 0). Singlet: sphere R = 51.68 at z = 5,
  // plane rear side at z = 9, stop at z = 0; the cylinder lens (f = 50) behind it at z = 20 or
  // in front of the first vertex at z = 2.
  const std::string sphere = R"({"type": "conic", "radius": 51.68})";
  struct Case {
    const char* name;
    std::string element_pose;
    double psi_deg;
  };
  for (const double z_lens : {20.0, 2.0}) {
    for (const Case& c : {Case{"psi = 90", pose(z_lens), 90.0},
                          Case{"psi = 0, rotated 90 about z", pose(z_lens, 0.0, 0.0, 90.0), 0.0}}) {
      INFO(c.name << ", lens at z = " << z_lens);
      const std::string cylinder =
          thin(c.element_pose, R"({"type": "ideal_cylinder_lens", "focal_length": 50.0,)"
                               R"( "axis_deg": )" +
                                   json_number(c.psi_deg) + "}");
      const std::string singlet = lens("L1", pose(5.0), sphere, pose(0.0), 4.0);
      const CompiledSystem cs = z_lens > 9.0
                                    ? build({stop(0.0, 10.0), singlet, cylinder, detector(110.0)})
                                    : build({stop(0.0, 10.0), cylinder, singlet, detector(110.0)});
      const auto events = [z_lens](double phi) {
        std::vector<Ev> ev = {{0.0, 0.0, 1.0}, {5.0, 1.0 / kRy, kN}, {9.0, 0.0, 1.0}};
        const Ev il{z_lens, 0.0, 1.0, phi};
        ev.insert(z_lens > 9.0 ? ev.end() : ev.begin() + 1, il);
        ev.push_back({110.0, 0.0, 1.0});
        return ev;
      };
      // Content guard: in y the singlet alone (bfl from z = 9, ffl from z = 5); in x both.
      const Expected ey = expected(events(0.0), 0, std::nullopt);
      const Expected ex = expected(events(1.0 / 50.0), 0, std::nullopt);
      REQUIRE(close(*ey.bfl, 100.0 - 4.0 / kN, 1e-12));
      REQUIRE(std::abs(*ex.efl - *ey.efl) > 1.0);
      check(section(cs, Section::y), ey);
      check(section(cs, Section::x), ex);
    }
  }
}
