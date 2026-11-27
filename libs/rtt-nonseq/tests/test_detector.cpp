// Sources, detectors and the ray path report of the non-sequential tracer (ADR 0057, points 6-8,
// Nachtrag PR 4): acceptance cases 3, 4 and 6 of point 11 and the rules of the fields, the report
// and the hit records.

#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/nonseq/scene.hpp"
#include "rtt/nonseq/source.hpp"
#include "rtt/nonseq/trace.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "test_support.hpp"

using rtt::compile::CompiledSystem;
using rtt::material::MaterialLibrary;
using rtt::math::Isometry3;
using rtt::math::Mat3;
using rtt::math::Vec3;
using rtt::model::CircularAperture;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Node;
using rtt::model::Pose;
using rtt::model::RectangularAperture;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using namespace rtt::nonseq;

namespace {

constexpr double kPi = std::numbers::pi;

CompiledSystem compiled(const System& s) {
  const MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

Surface detector_surface(double hx, double hy) {
  Surface image;
  image.id = SurfaceId("IMG");
  image.aperture = RectangularAperture{hx, hy};
  return image;
}

/// Only a detector: plane, rectangular aperture of half widths hx, hy at z = d.
System detector_only(double d, double hx, double hy) {
  System s = rtt::model::test::make_singlet();
  s.root.children.clear();
  s.root.children.push_back({Element{
      "image", ElementKind::Detector, Pose::along_z(d), std::nullopt, {detector_surface(hx, hy)}}});
  return s;
}

/// A plane-parallel plate of CONST:1.5 at z = 5 ... 9 before a detector of half width 20 at
/// z = 106.363.
System plate_system() {
  System s = detector_only(106.363, 20.0, 20.0);
  Surface a;
  a.id = SurfaceId("P.S1");
  a.aperture = CircularAperture{12.7, 0.0};
  Surface b = a;
  b.id = SurfaceId("P.S2");
  b.pose = Pose::along_z(4.0);
  s.root.children.insert(
      s.root.children.begin(),
      Node{Element{"P", ElementKind::Plate, Pose::along_z(5.0), std::string("CONST:1.5"), {a, b}}});
  return s;
}

std::uint32_t face_of(const Scene& scene, FaceKind kind) {
  for (std::uint32_t f = 0; f < scene.faces().size(); ++f) {
    if (scene.faces()[f].kind == kind) return f;
  }
  return kNoFace;
}

std::uint32_t body_face(const Scene& scene, const CompiledSystem& system, const char* id) {
  for (std::uint32_t f = 0; f < scene.faces().size(); ++f) {
    const Face& face = scene.faces()[f];
    if (face.kind == FaceKind::Body && system.surfaces()[face.system_surface].id.str() == id)
      return f;
  }
  return kNoFace;
}

void check_closed(const NonseqResult& r) {
  for (const Balance& b : r.balance) {
    CHECK(std::abs(b.sinks() - b.source) <= 1e-12 * (1.0 + b.source));
  }
}

Source collimated_rect(double hx, double hy, std::uint64_t n, double z) {
  Source s;
  s.emitter = CollimatedEmitter{Rectangle{hx, hy}};
  s.pose = Isometry3(Mat3::Identity(), Vec3(0.0, 0.0, z));
  s.power = 2.0;
  s.rays = n;
  s.sampling = Sampling::Grid;
  return s;
}

}  // namespace

TEST_CASE("Irradiance: a collimated grid on matching pixels gives Phi / (M A_pixel) each",
          "[nonseq][detector]") {
  // A rectangle 4 mm x 4 mm on a 4 x 4 grid (n = 4, M = 16) shines on a detector of the same size
  // with 4 x 4 pixels of 1 mm^2: one ray per pixel, E = (2 W / 16) / 1 mm^2 = 0.125 W/mm^2.
  const CompiledSystem system = compiled(detector_only(50.0, 2.0, 2.0));
  const Scene scene = Scene::build(system);
  const std::uint32_t det = face_of(scene, FaceKind::Detector);
  const std::vector<Source> src{collimated_rect(2.0, 2.0, 16, -1.0)};
  const std::vector<Detector> dets{IrradianceDetector{det, 4, 4}};
  const NonseqResult r = trace(system, scene, src, dets);
  REQUIRE(r.detectors.size() == 1);
  const DetectorField& f = r.detectors[0];
  CHECK(f.rows == 4);
  CHECK(f.cols == 4);
  REQUIRE(f.value.size() == 16);
  REQUIRE(f.cell_size.size() == 1);
  CHECK(std::abs(f.cell_size[0] - 1.0) < 1e-15);
  for (std::size_t k = 0; k < 16; ++k) {
    CHECK(std::abs(f.value[k] - 0.125) < 1e-15);
    CHECK(f.hits[k] == 1);
    CHECK(f.weights[k] == 1.0);
  }
  CHECK(std::abs(r.balance[0].detected[det] - 2.0) < 1e-15);
  CHECK(r.balance[0].source == 2.0);
  CHECK(r.rays == 16);
}

TEST_CASE("Irradiance: pixel rows run along local y, columns along local x", "[nonseq][detector]") {
  // A thin collimated strip at x in [1, 2], y in [-2, -1.5] lands only in column 3 (x in [1, 2])
  // and row 0 (y in [-2, -1]).
  const CompiledSystem system = compiled(detector_only(50.0, 2.0, 2.0));
  const Scene scene = Scene::build(system);
  const std::uint32_t det = face_of(scene, FaceKind::Detector);
  Source s = collimated_rect(0.5, 0.25, 4, -1.0);
  s.pose = Isometry3(Mat3::Identity(), Vec3(1.5, -1.75, -1.0));
  const NonseqResult r = trace(system, scene, std::vector<Source>{s},
                               std::vector<Detector>{IrradianceDetector{det, 4, 4}});
  const DetectorField& f = r.detectors[0];
  for (std::uint32_t row = 0; row < 4; ++row) {
    for (std::uint32_t col = 0; col < 4; ++col) {
      const std::uint64_t h = f.hits[row * 4 + col];
      CHECK(h == ((row == 0 && col == 3) ? 4U : 0U));
    }
  }
}

TEST_CASE("Irradiance detectors need a plane with a rectangular aperture", "[nonseq][detector]") {
  const auto codes = [](const System& sys) {
    const CompiledSystem system = compiled(sys);
    const Scene scene = Scene::build(system);
    const std::vector<Detector> dets{IrradianceDetector{face_of(scene, FaceKind::Detector), 2, 2}};
    try {
      (void)trace(system, scene, std::vector<Source>{collimated_rect(1.0, 1.0, 4, -1.0)}, dets);
    } catch (const SceneError& e) {
      std::vector<std::string> c;
      for (const auto& d : e.diagnostics()) c.push_back(d.code + " " + d.location);
      return c;
    }
    return std::vector<std::string>{};
  };
  System circ = detector_only(50.0, 2.0, 2.0);
  std::get<Element>(circ.root.children[0].value).surfaces[0].aperture = CircularAperture{2.0, 0.0};
  const auto c1 = codes(circ);
  REQUIRE(c1.size() == 1);
  CHECK(c1[0] == "nonseq.detector_aperture /root/children/0/surfaces/0/aperture");
  System curved = detector_only(50.0, 2.0, 2.0);
  std::get<Element>(curved.root.children[0].value).surfaces[0].shape.base =
      rtt::model::Conic{rtt::model::Param(-100.0), rtt::model::Param(0.0)};
  const auto c2 = codes(curved);
  REQUIRE(c2.size() == 1);
  CHECK(c2[0] == "nonseq.detector_shape /root/children/0/surfaces/0/shape");
  CHECK(codes(detector_only(50.0, 2.0, 2.0)).empty());
}

TEST_CASE("Invalid detectors and options are rejected", "[nonseq][detector]") {
  const CompiledSystem system = compiled(plate_system());
  const Scene scene = Scene::build(system);
  const std::uint32_t det = face_of(scene, FaceKind::Detector);
  const std::vector<Source> src{collimated_rect(1.0, 1.0, 4, -1.0)};
  const auto bad = [&](const Detector& d) {
    CHECK_THROWS_AS(trace(system, scene, src, std::vector<Detector>{d}), std::invalid_argument);
  };
  bad(IrradianceDetector{0, 2, 2});  // a plate face
  bad(IrradianceDetector{kNoFace, 2, 2});
  bad(IrradianceDetector{det, 0, 2});
  bad(IntensityDetector{Isometry3(), 0.0, 2, 2});
  bad(IntensityDetector{Isometry3(), kPi + 1e-9, 2, 2});
  bad(IntensityDetector{Isometry3(), 1.0, 2, 0});
  TraceOptions o;
  o.report_min_share = 1.5;
  CHECK_THROWS_AS(trace(system, scene, src, std::vector<Detector>{}, o), std::invalid_argument);
  Source s = src[0];
  s.power = -1.0;
  CHECK_THROWS_AS(trace(system, scene, std::vector<Source>{s}, std::vector<Detector>{}),
                  std::invalid_argument);
}

TEST_CASE("Source power splits by the wavelength weights", "[nonseq][detector]") {
  System sys = detector_only(50.0, 5.0, 5.0);
  sys.wavelengths = {{0.5, 1.0, true}, {0.6, 3.0, false}, {0.7, 4.0, false}};
  const CompiledSystem system = compiled(sys);
  const Scene scene = Scene::build(system);
  const std::uint32_t det = face_of(scene, FaceKind::Detector);
  Source s = collimated_rect(1.0, 1.0, 10, -1.0);  // grid: M = 16
  s.wavelengths = {2, 1};
  const NonseqResult r = trace(system, scene, std::vector<Source>{s},
                               std::vector<Detector>{IrradianceDetector{det, 1, 1}});
  REQUIRE(r.balance.size() == 3);
  // Phi_lambda = 2 W * g / (3 + 4).
  CHECK(r.balance[0].source == 0.0);
  CHECK(std::abs(r.balance[1].source - 2.0 * 3.0 / 7.0) < 1e-14);
  CHECK(std::abs(r.balance[2].source - 2.0 * 4.0 / 7.0) < 1e-14);
  CHECK(r.rays == 32);
  const DetectorField& f = r.detectors[0];
  REQUIRE(f.value.size() == 3);
  CHECK(f.hits[0] == 0);
  CHECK(f.hits[1] == 16);
  CHECK(f.hits[2] == 16);
  CHECK(std::abs(f.value[1] * 100.0 - 2.0 * 3.0 / 7.0) < 1e-13);  // A_pixel = 100 mm^2
  check_closed(r);
}

TEST_CASE("case 3: Lambert source before a parallel detector, E = E0 cos^4 theta",
          "[nonseq][detector][acceptance]") {
  // NBS TN 910-1, Eq. (2.15): dPhi = L cos theta_1 dA_1 cos theta_2 dA_2 / D^2. For a small
  // source of area A parallel to the detector at distance d, theta_1 = theta_2 = theta and
  // D = d / cos theta give E(theta) = L A cos^4 theta / d^2 (ADR 0057, Quellen). With
  // L = Phi / (pi A) (Eq. (A2-14), max_angle 90 deg), E0 = Phi / (pi d^2).
  //
  // Random sampling (#450, decision 9; Nachtrag PR 4, a), N rays: the power in a pixel is
  // binomial with p = P_pix / Phi, so the tolerance is 4 sigma = 4 sqrt(N p (1 - p)) Phi / N,
  // plus the finite source: a source point at offset s shifts every hit by s, so the pixel mean
  // of E moves by at most a max|grad E| = a E0 (4/d) u (1 + u^2)^-3 at u^2 = 1/5, i.e.
  // 1.0352 a E0 / d. The pixel mean of E0 cos^4 theta is a midpoint sum over 100 x 100 points;
  // its error (h/d)^2 ~ 1e-6 relative is far below 4 sigma.
  const double d = 50.0;
  const double half = 50.0;
  const double a = 0.05;
  const double phi = 1.0;
  const std::uint64_t n = 400000;
  const CompiledSystem system = compiled(detector_only(d, half, half));
  const Scene scene = Scene::build(system);
  const std::uint32_t det = face_of(scene, FaceKind::Detector);
  Source s;
  s.emitter = LambertEmitter{Disc{a}, kPi / 2};
  s.power = phi;
  s.rays = n;
  TraceOptions o;
  o.seed = 7;
  const std::uint32_t pixels = 10;
  const NonseqResult r = trace(system, scene, std::vector<Source>{s},
                               std::vector<Detector>{IrradianceDetector{det, pixels, pixels}}, o);
  const DetectorField& f = r.detectors[0];
  const double e0 = phi / (kPi * d * d);
  const double w = 2.0 * half / pixels;
  const double area = w * w;
  REQUIRE(std::abs(f.cell_size[0] - area) < 1e-12);
  const int m = 100;
  double worst = 0.0;
  for (std::uint32_t row = 0; row < pixels; ++row) {
    for (std::uint32_t col = 0; col < pixels; ++col) {
      double mean = 0.0;
      for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
          const double x = -half + (col + (i + 0.5) / m) * w;
          const double y = -half + (row + (j + 0.5) / m) * w;
          const double c2 = d * d / (d * d + x * x + y * y);  // cos^2 theta
          mean += e0 * c2 * c2;
        }
      }
      mean /= m * m;
      const double p = mean * area / phi;
      const double sigma_e =
          std::sqrt(static_cast<double>(n) * p * (1.0 - p)) * phi / static_cast<double>(n) / area;
      const double tol = 4.0 * sigma_e + 1.0352 * a * e0 / d;
      const double got = f.value[row * pixels + col];
      INFO("pixel " << row << ", " << col << ": E = " << got << ", expected " << mean);
      CHECK(std::abs(got - mean) <= tol);
      worst = std::max(worst, std::abs(got - mean) / tol);
    }
  }
  // Content guard: the field falls off from the centre to the corner as cos^4 theta.
  CHECK(f.value[0] < 0.5 * f.value[4 * pixels + 4]);
  // Detected power on the detector plus escaped power closes the balance.
  check_closed(r);
  CHECK(std::abs(r.balance[0].detected[det] + r.balance[0].escaped - phi) < 1e-12);
}

TEST_CASE("Lambert grid: n^4 rays, Phi / M each, all power on a large detector (NBS (A2-14))",
          "[nonseq][detector]") {
  // max_angle 30 deg: every ray reaches a detector of half width 50 mm at d = 50 mm
  // (tan 30 deg * 50 + 1 < 50), so the detected power is Phi = pi L A sin^2 theta_max.
  const CompiledSystem system = compiled(detector_only(50.0, 50.0, 50.0));
  const Scene scene = Scene::build(system);
  const std::uint32_t det = face_of(scene, FaceKind::Detector);
  Source s;
  s.emitter = LambertEmitter{Rectangle{1.0, 1.0}, kPi / 6};
  s.power = 0.3;
  s.rays = 17;  // n = 3, M = 81
  s.sampling = Sampling::Grid;
  const NonseqResult r = trace(system, scene, std::vector<Source>{s},
                               std::vector<Detector>{IrradianceDetector{det, 1, 1}});
  CHECK(r.rays == 81);
  CHECK(r.detectors[0].hits[0] == 81);
  CHECK(std::abs(r.balance[0].detected[det] - 0.3) < 1e-15);
  CHECK(r.balance[0].escaped == 0.0);
}

TEST_CASE("case 6: polar intensity of a point source, I = Phi / (2 pi (1 - cos theta_max))",
          "[nonseq][detector][acceptance]") {
  // PBRT 4ed, §A.5.4: uniform directions in the cone; NBS TN 910-1, Eq. (A2-10): the cone has
  // 2 pi (1 - cos theta_max) sr, so I is the same in every cell. Grid with n = 1000 (M = 10^6):
  // phi = 2 pi xi_2 at the cell centres never lies on a boundary of the n_phi = 8 cells
  // (n / n_phi = 125 columns each), and 1 - cos theta is linear in xi_1, so a theta cell holds
  // the centres of an interval of xi_1: at most one row more or less than n times its length.
  // One row in one phi cell carries Phi / (n n_phi), so |I - I_exact| <= Phi / (n n_phi DeltaOmega)
  // (plus rounding).
  const double tmax = kPi / 3;
  const double phi = 1.5;
  const CompiledSystem system = compiled(detector_only(-10.0, 5.0, 5.0));  // behind the source
  const Scene scene = Scene::build(system);
  Source s;
  s.emitter = PointEmitter{tmax};
  s.power = phi;
  s.rays = 1000000;
  s.sampling = Sampling::Grid;
  const std::uint32_t nt = 10;
  const std::uint32_t np = 8;
  TraceOptions o;
  o.path_report = false;
  const std::vector<Detector> dets{IntensityDetector{Isometry3(), tmax, nt, np},
                                   IntensityDetector{Isometry3(), 2.0 * tmax, nt, np},
                                   IntensityDetector{Isometry3(), tmax / 2.0, 5, np}};
  const NonseqResult r = trace(system, scene, std::vector<Source>{s}, dets, o);
  CHECK(std::abs(r.balance[0].escaped - phi) < 1e-12);
  const double expected = phi / (2.0 * kPi * (1.0 - std::cos(tmax)));
  const DetectorField& f = r.detectors[0];
  REQUIRE(f.rows == nt);
  REQUIRE(f.cols == np);
  REQUIRE(f.cell_size.size() == nt);
  std::uint64_t hits = 0;
  for (std::uint32_t row = 0; row < nt; ++row) {
    const double t1 = row * tmax / nt;
    const double t2 = (row + 1) * tmax / nt;
    const double domega = (2.0 * kPi / np) * (std::cos(t1) - std::cos(t2));
    CHECK(std::abs(f.cell_size[row] - domega) < 1e-15);
    const double tol = phi / (1000.0 * np * domega) + 1e-12 * expected;
    for (std::uint32_t col = 0; col < np; ++col) {
      INFO("cell " << row << ", " << col);
      CHECK(std::abs(f.value[row * np + col] - expected) <= tol);
      hits += f.hits[row * np + col];
    }
  }
  CHECK(hits == 1000000);
  // The wider detector sees nothing beyond theta_max of the source.
  const DetectorField& g = r.detectors[1];
  for (std::uint32_t row = nt / 2; row < nt; ++row) {
    for (std::uint32_t col = 0; col < np; ++col) CHECK(g.hits[row * np + col] == 0);
  }
  CHECK(std::abs(g.value[0] - expected) <= phi / (1000.0 * np * g.cell_size[0]) + 1e-12);
  // The narrower detector counts only theta <= 30 deg: the rows with (k + 1/2)/1000 <=
  // (1 - cos 30 deg)/(1 - cos 60 deg) = sin^2(15 deg)/sin^2(30 deg) = 0.26795, k <= 267, i.e.
  // 268 rows of 1000 rays.
  std::uint64_t narrow = 0;
  for (const std::uint64_t h : r.detectors[2].hits) narrow += h;
  CHECK(narrow == 268000);
}

TEST_CASE("Intensity cells are measured in the detector pose", "[nonseq][detector]") {
  // Detector axis tilted by 40 deg about x: local -> global R = rot_x(40 deg). A ray along
  // global +z has local direction R^T z = (0, sin 40, cos 40): theta = 40 deg (row 1 of
  // 3 x 30 deg), phi = 90 deg (column 0 of 3 x 120 deg).
  const CompiledSystem system = compiled(detector_only(-10.0, 5.0, 5.0));
  const Scene scene = Scene::build(system);
  const double al = 40.0 * kPi / 180.0;
  const Mat3 rot =
      (Mat3() << 1, 0, 0, 0, std::cos(al), -std::sin(al), 0, std::sin(al), std::cos(al)).finished();
  const std::vector<Detector> dets{
      IntensityDetector{Isometry3(rot, Vec3(5.0, 0.0, 0.0)), kPi / 2, 3, 3}};
  const NonseqResult r =
      trace(system, scene, std::vector<Source>{collimated_rect(1.0, 1.0, 4, 0.0)}, dets);
  const DetectorField& f = r.detectors[0];
  for (std::uint32_t k = 0; k < 9; ++k) CHECK(f.hits[k] == (k == 3 ? 4U : 0U));
}

TEST_CASE("A transparent irradiance detector counts both directions", "[nonseq][detector]") {
  System sys = detector_only(50.0, 2.0, 2.0);
  Surface m;
  m.id = SurfaceId("M");
  m.aperture = CircularAperture{10.0, 0.0};
  m.interaction = rtt::model::IdealMirror{};
  sys.root.children.push_back(
      {Element{"M", ElementKind::Mirror, Pose::along_z(60.0), std::nullopt, {m}}});
  const CompiledSystem system = compiled(sys);
  SceneOptions so;
  so.transparent_detectors = {0};
  const Scene scene = Scene::build(system, so);
  const std::uint32_t det = face_of(scene, FaceKind::Detector);
  const NonseqResult r =
      trace(system, scene, std::vector<Source>{collimated_rect(2.0, 2.0, 4, -1.0)},
            std::vector<Detector>{IrradianceDetector{det, 2, 2}});
  const DetectorField& f = r.detectors[0];
  for (std::size_t k = 0; k < 4; ++k) {
    CHECK(f.hits[k] == 2);
    CHECK(std::abs(f.value[k] - 2.0 * (2.0 / 4.0) / 4.0) < 1e-15);  // 2 passes, A_pixel 4 mm^2
  }
  CHECK(std::abs(r.balance[0].passed[det] - 4.0) < 1e-15);
  CHECK(r.balance[0].detected[det] == 0.0);
  CHECK(std::abs(r.balance[0].escaped - 2.0) < 1e-15);
}

TEST_CASE("case 4: point source through a lens hits the detector where the sequential trace does",
          "[nonseq][detector][acceptance]") {
  // Ray by ray (Nachtrag Review PR 1, f): the same start rays traced sequentially along the
  // automatic path and non-sequentially with full splitting; ideal_anti_reflection on the lens
  // faces leaves only the transmitted branch (Reflect has weight 0), so each start ray gives
  // exactly one hit. Both run the same intersect_surface and apply_event; the limit 1e-9 mm is
  // geom::kDefaultTMin, the self-hit bound of the scene.
  System sys = rtt::model::test::make_singlet();
  sys.root.children.erase(sys.root.children.begin());  // no stop
  auto& lens = std::get<Element>(sys.root.children[0].value);
  lens.material = "CONST:1.5168";  // N-BK7 at the d line, without the catalog
  for (Surface& surf : lens.surfaces) surf.interaction = rtt::model::IdealAntiReflection{};
  std::get<Element>(sys.root.children[1].value).surfaces[0].aperture =
      RectangularAperture{20.0, 20.0};
  const CompiledSystem system = compiled(sys);
  const Scene scene = Scene::build(system);
  const std::uint32_t det = face_of(scene, FaceKind::Detector);

  Source s;
  s.emitter = PointEmitter{8.0 * kPi / 180.0};
  s.pose = Isometry3(Mat3::Identity(), Vec3(0.0, 0.0, -50.0));
  s.power = 1.0;
  s.rays = 400;
  s.sampling = Sampling::Grid;
  const std::uint64_t m = source_ray_count(s);
  REQUIRE(m == 400);
  TraceOptions o;
  o.splitting = Splitting::Full;
  o.record_hits = true;
  const NonseqResult r = trace(system, scene, std::vector<Source>{s},
                               std::vector<Detector>{IrradianceDetector{det, 4, 4}}, o);
  REQUIRE(r.hits.size() == m);

  rtt::trace::RayBatch rays(m);
  for (std::uint64_t i = 0; i < m; ++i) {
    const SourceRay sr = source_ray(s, 0, i, o.seed);
    rays.pos_x()[i] = sr.position.x();
    rays.pos_y()[i] = sr.position.y();
    rays.pos_z()[i] = sr.position.z();
    rays.dir_x()[i] = sr.direction.x();
    rays.dir_y()[i] = sr.direction.y();
    rays.dir_z()[i] = sr.direction.z();
  }
  const rtt::trace::SequentialTracer tracer;
  (void)tracer.trace(system, rtt::compile::PathId{0}, rays);
  double spread = 0.0;
  for (const DetectorHit& h : r.hits) {
    REQUIRE(h.index < m);
    const std::size_t i = h.index;
    REQUIRE(rays.status()[i] == rtt::trace::RayStatus::Alive);
    const Vec3 seq(rays.pos_x()[i], rays.pos_y()[i], rays.pos_z()[i]);
    INFO("ray " << i);
    CHECK((h.point - seq).norm() <= 1e-9);
    spread = std::max(spread, std::hypot(h.point.x(), h.point.y()));
  }
  // Content guard: the hits are spread over the detector, not all at one point.
  CHECK(spread > 0.1);
  CHECK(std::abs(r.balance[0].detected[det] - 1.0) < 1e-12);
}

TEST_CASE("Ray path report of a Fresnel plate", "[nonseq][report]") {
  const CompiledSystem system = compiled(plate_system());
  const Scene scene = Scene::build(system);
  const std::uint32_t det = face_of(scene, FaceKind::Detector);
  const std::uint32_t s1 = body_face(scene, system, "P.S1");
  const std::uint32_t s2 = body_face(scene, system, "P.S2");
  REQUIRE(s1 != kNoFace);
  REQUIRE(s2 != kNoFace);
  const Source src = collimated_rect(1.0, 1.0, 1, -1.0);  // one ray along the axis
  TraceOptions o;
  o.splitting = Splitting::Full;
  o.min_relative_power = 0.0;
  o.report_min_share = 0.01;
  const NonseqResult r = trace(system, scene, std::vector<Source>{src}, std::vector<Detector>{}, o);
  REQUIRE(r.paths.size() == 1);
  const PathReport& rep = r.paths[0];
  CHECK(rep.source == 2.0);
  // Detector, escaped, and the branches cut at max_depth (power ~ R^64, last).
  REQUIRE(rep.endpoints.size() == 3);
  CHECK(rep.endpoints[2].end == PathEnd::Cut);
  CHECK(rep.endpoints[2].power < 1e-80);
  CHECK(rep.endpoints[2].paths.empty());
  const double n0 = system.media()[system.environment_medium()].index[0].real();
  const double rr = std::pow((1.5 - n0) / (1.5 + n0), 2);
  using rtt::model::EventKind;
  // Detector first (T = (1 - R)/(1 + R) ~ 0.92), then escaped.
  const PathEndpoint& e0 = rep.endpoints[0];
  CHECK(e0.end == PathEnd::Detected);
  CHECK(e0.face == det);
  CHECK(std::abs(e0.power - r.balance[0].detected[det]) < 1e-15);
  REQUIRE(e0.paths.size() == 1);  // R^2 T^2 ~ 0.0015 < 0.01
  CHECK(e0.paths[0].steps ==
        std::vector<PathStep>{{s1, EventKind::Refract, 0}, {s2, EventKind::Refract, 0}});
  CHECK(std::abs(e0.paths[0].power - 2.0 * (1 - rr) * (1 - rr)) < 1e-14);
  CHECK(std::abs(e0.paths[0].share - (1 - rr) * (1 - rr)) < 1e-14);
  CHECK(e0.paths[0].rays == 1);
  const PathEndpoint& e1 = rep.endpoints[1];
  CHECK(e1.end == PathEnd::Escaped);
  CHECK(e1.face == kNoFace);
  CHECK(std::abs(e1.power - r.balance[0].escaped) < 1e-15);
  REQUIRE(e1.paths.size() == 2);
  CHECK(e1.paths[0].steps == std::vector<PathStep>{{s1, EventKind::Reflect, 0}});
  CHECK(std::abs(e1.paths[0].share - rr) < 1e-15);
  CHECK(e1.paths[1].steps == std::vector<PathStep>{{s1, EventKind::Refract, 0},
                                                   {s2, EventKind::Reflect, 0},
                                                   {s1, EventKind::Refract, 0}});
  CHECK(std::abs(e1.paths[1].share - rr * (1 - rr) * (1 - rr)) < 1e-15);
  CHECK(e1.rays > 2);  // the unlisted paths still count
}

TEST_CASE("Monte Carlo path report: every start ray ends on one path", "[nonseq][report]") {
  const CompiledSystem system = compiled(plate_system());
  const Scene scene = Scene::build(system);
  Source src = collimated_rect(1.0, 1.0, 3000, -1.0);
  src.sampling = Sampling::Random;
  TraceOptions o;
  o.report_min_share = 0.0;
  const NonseqResult r = trace(system, scene, std::vector<Source>{src}, std::vector<Detector>{}, o);
  std::uint64_t rays = 0;
  double power = 0.0;
  for (const PathEndpoint& e : r.paths[0].endpoints) {
    rays += e.rays;
    power += e.power;
    std::uint64_t listed = 0;
    for (const RayPath& p : e.paths) listed += p.rays;
    CHECK(listed == e.rays);  // report_min_share 0 lists all
    for (std::size_t k = 1; k < e.paths.size(); ++k)
      CHECK(e.paths[k - 1].power >= e.paths[k].power);
  }
  CHECK(rays == 3000);
  // No absorption: the end points carry the whole source.
  CHECK(std::abs(power - 2.0) < 1e-12);
  TraceOptions off = o;
  off.path_report = false;
  CHECK(trace(system, scene, std::vector<Source>{src}, std::vector<Detector>{}, off).paths.empty());
}

TEST_CASE("Sources, fields and report are bitwise the same for 1, 2 and 8 threads",
          "[nonseq][detector]") {
  // ADR 0057, point 5 and case 7: blocks of 4096 start rays, summed in block order.
  const CompiledSystem system = compiled(plate_system());
  const Scene scene = Scene::build(system);
  const std::uint32_t det = face_of(scene, FaceKind::Detector);
  Source lam;
  lam.emitter = LambertEmitter{Disc{3.0}, 0.3};
  lam.pose = Isometry3(Mat3::Identity(), Vec3(0.0, 0.0, -1.0));
  lam.power = 1.0;
  lam.rays = 10000;
  Source pt;
  pt.emitter = PointEmitter{0.2};
  pt.power = 0.5;
  pt.rays = 5000;
  pt.sampling = Sampling::Grid;
  const std::vector<Source> sources{lam, pt};
  const std::vector<Detector> dets{IrradianceDetector{det, 8, 8},
                                   IntensityDetector{Isometry3(), kPi, 6, 4}};
  TraceOptions o;
  o.seed = 99;
  o.record_hits = true;
  const auto run = [&](int threads) {
    tbb::task_arena arena(threads);
    return arena.execute([&] { return trace(system, scene, sources, dets, o); });
  };
  const NonseqResult a = run(1);
  REQUIRE(a.detectors[0].hits[27] + a.detectors[0].hits[36] > 0);
  for (const int t : {2, 8}) {
    const NonseqResult b = run(t);
    for (std::size_t k = 0; k < a.detectors.size(); ++k) {
      CHECK(a.detectors[k].value == b.detectors[k].value);
      CHECK(a.detectors[k].hits == b.detectors[k].hits);
      CHECK(a.detectors[k].weights == b.detectors[k].weights);
    }
    REQUIRE(a.hits.size() == b.hits.size());
    for (std::size_t k = 0; k < a.hits.size(); ++k) {
      CHECK(a.hits[k].index == b.hits[k].index);
      CHECK(a.hits[k].point == b.hits[k].point);
    }
    REQUIRE(a.paths.size() == b.paths.size());
    for (std::size_t s = 0; s < a.paths.size(); ++s) {
      REQUIRE(a.paths[s].endpoints.size() == b.paths[s].endpoints.size());
      for (std::size_t e = 0; e < a.paths[s].endpoints.size(); ++e) {
        const PathEndpoint& x = a.paths[s].endpoints[e];
        const PathEndpoint& y = b.paths[s].endpoints[e];
        CHECK(x.power == y.power);
        CHECK(x.rays == y.rays);
        REQUIRE(x.paths.size() == y.paths.size());
        for (std::size_t p = 0; p < x.paths.size(); ++p) {
          CHECK(x.paths[p].steps == y.paths[p].steps);
          CHECK(x.paths[p].power == y.paths[p].power);
        }
      }
    }
    CHECK(a.balance[0].detected == b.balance[0].detected);
    CHECK(a.balance[1].escaped == b.balance[1].escaped);
  }
  check_closed(a);
}

TEST_CASE("case 5 with a source outside: every ray through the entry face reaches the end",
          "[nonseq][detector][acceptance]") {
  // A rod of CONST:1.5 (radius 2, z = 5 ... 55) with ideal_anti_reflection on both end faces and
  // a point source 1 mm before the entry face, cone 60 deg: every ray meets the entry face
  // (tan 60 deg * 1 mm < 2 mm). Inside, sin theta = sin theta_air / 1.5 < 1/1.5, so the angle at
  // the mantle, 90 deg - theta, exceeds the critical angle asin(1/1.5) (de Greve, Eq. (25)) and
  // every reflection there is total; the exit face sees theta < 41.8 deg and transmits fully
  // (Nachtrag Review PR 1, f). The detector 1e-6 mm behind the exit face gets the whole power.
  System sys = detector_only(55.0 + 1e-6, 2.0 + 1e-5, 2.0 + 1e-5);
  Surface a;
  a.id = SurfaceId("R.S1");
  a.aperture = CircularAperture{2.0, 0.0};
  a.interaction = rtt::model::IdealAntiReflection{};
  Surface b = a;
  b.id = SurfaceId("R.S2");
  b.pose = Pose::along_z(50.0);
  sys.root.children.insert(
      sys.root.children.begin(),
      Node{Element{"R", ElementKind::Plate, Pose::along_z(5.0), std::string("CONST:1.5"), {a, b}}});
  const CompiledSystem system = compiled(sys);
  const Scene scene = Scene::build(system);
  const std::uint32_t det = face_of(scene, FaceKind::Detector);
  Source s;
  s.emitter = PointEmitter{kPi / 3};
  s.pose = Isometry3(Mat3::Identity(), Vec3(0.0, 0.0, 4.0));
  s.power = 1.0;
  s.rays = 400;
  s.sampling = Sampling::Grid;
  TraceOptions o;
  o.splitting = Splitting::Full;
  const NonseqResult r = trace(system, scene, std::vector<Source>{s},
                               std::vector<Detector>{IrradianceDetector{det, 1, 1}}, o);
  CHECK(std::abs(r.balance[0].detected[det] - 1.0) < 1e-12);
  CHECK(r.detectors[0].hits[0] == 400);
  // Content guard: the steeper rays are reflected at the mantle on the way.
  CHECK(r.events > 3 * 400);
  check_closed(r);
}
