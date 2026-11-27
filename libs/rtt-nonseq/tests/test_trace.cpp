// Trace core of the non-sequential tracer (ADR 0057, PR 3): acceptance cases 1, 2, 5 and 7 of
// point 11 (Fresnel plate with full splitting and with Monte Carlo, light guide with TIR,
// balance and determinism), and the branch rules of point 4 at mirrors, beam splitters,
// transparent detectors and a crystal.

#include <oneapi/tbb/task_arena.h>

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/nonseq/scene.hpp"
#include "rtt/nonseq/trace.hpp"
#include "test_support.hpp"

using rtt::compile::CompiledSystem;
using rtt::material::MaterialLibrary;
using rtt::math::Vec3;
using rtt::model::CircularAperture;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::Node;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::nonseq::Balance;
using rtt::nonseq::FaceKind;
using rtt::nonseq::NonseqResult;
using rtt::nonseq::RaySource;
using rtt::nonseq::Scene;
using rtt::nonseq::Splitting;
using rtt::nonseq::StartRay;
using rtt::nonseq::TraceOptions;

namespace {

CompiledSystem compiled(const System& s) {
  const MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

/// make_singlet without the stop: `body` replaces the lens (index 0 afterwards), the detector
/// without aperture at z = 106.363 stays (index 1).
System with_body(const Element& body) {
  System s = rtt::model::test::make_singlet();
  s.root.children.erase(s.root.children.begin());
  s.root.children[0] = Node{body};
  return s;
}

/// A plane-parallel plate (CONST:1.5 by default) with Fresnel surfaces at z = 5 and
/// z = 5 + thickness.
Element plate(double thickness, double radius, const char* material = "CONST:1.5") {
  Surface a;
  a.id = SurfaceId("P.S1");
  a.aperture = CircularAperture{radius, 0.0};
  Surface b = a;
  b.id = SurfaceId("P.S2");
  b.pose = Pose::along_z(thickness);
  return Element{"P", ElementKind::Plate, Pose::along_z(5.0), material, {a, b}};
}

std::uint32_t face_of(const Scene& scene, FaceKind kind) {
  for (std::uint32_t f = 0; f < scene.faces().size(); ++f) {
    if (scene.faces()[f].kind == kind) return f;
  }
  return rtt::nonseq::kNoFace;
}

double sum(const std::vector<double>& v) {
  double s = 0.0;
  for (const double x : v) s += x;
  return s;
}

/// Power reflectances |r_s|^2 and |r_p|^2 of an interface n1 -> n2 at incidence angle i
/// (S. J. Byrnes, arXiv:1603.02720v5, Fresnel coefficients of Eqs. (1)-(6), docs/quellen.md;
/// the sign convention of r_p does not enter |r_p|^2).
std::pair<double, double> fresnel_r(double n1, double n2, double i) {
  const double ci = std::cos(i);
  const double st = n1 / n2 * std::sin(i);
  const double ct = std::sqrt(1.0 - st * st);
  const double rs = (n1 * ci - n2 * ct) / (n1 * ci + n2 * ct);
  const double rp = (n2 * ci - n1 * ct) / (n2 * ci + n1 * ct);
  return {rs * rs, rp * rp};
}

/// Real indices of the environment (AIR, ~1.00027) and of the first medium of element 0 at
/// wavelength 0.
std::pair<double, double> indices(const CompiledSystem& s) {
  return {s.media()[s.environment_medium()].index[0].real(),
          s.media()[s.elements()[0].media[0]].index[0].real()};
}

/// Reflectance of one face at normal incidence, ((n - n0) / (n + n0))^2 (Byrnes, as above).
double normal_r(const CompiledSystem& s) {
  const auto [n0, n] = indices(s);
  return (n - n0) * (n - n0) / ((n + n0) * (n + n0));
}

/// Every sink of every balance adds up to its source (ADR 0057, point 8).
void check_closed(const NonseqResult& r) {
  for (const Balance& b : r.balance) {
    CHECK(std::abs(b.sinks() - b.source) <= 1e-12 * (1.0 + b.source));
  }
}

}  // namespace

TEST_CASE("case 1: Fresnel plate, full splitting, T = (1 - R) / (1 + R)", "[nonseq][trace]") {
  // Whittaker and Gehring, arXiv:1012.3677, Eqs. (1)-(3) (docs/quellen.md): a non-absorbing
  // plate with equal faces transmits T = (1 - R)/(1 + R) and reflects 2R/(1 + R), R the
  // reflectance of one face, the incoherent sum of all passes. With min_relative_power = 0 and
  // max_depth 64 the series is cut after R^62 < 1e-80 (Nachtrag Review PR 1, f): the limit is
  // the rounding of the sum, a few 1e-16.
  const CompiledSystem system = compiled(with_body(plate(4.0, 12.7)));
  const Scene scene = Scene::build(system);
  TraceOptions o;
  o.splitting = Splitting::Full;
  o.min_relative_power = 0.0;
  const std::uint32_t det = face_of(scene, FaceKind::Detector);

  SECTION("normal incidence: R = ((n - n0)/(n + n0))^2, about 0.04") {
    const std::vector<RaySource> src{{{{Vec3(0.0, 0.0, -1.0), Vec3::UnitZ(), 1.0, 0}}}};
    const NonseqResult r = rtt::nonseq::trace(system, scene, src, o);
    const Balance& b = r.balance[0];
    const double R = normal_r(system);
    CHECK(std::abs(b.detected[det] - (1.0 - R) / (1.0 + R)) < 1e-14);
    CHECK(std::abs(b.escaped - 2.0 * R / (1.0 + R)) < 1e-14);
    CHECK(b.cut < 1e-60);
    CHECK(r.events > 10);  // guard: the reflections were followed
    check_closed(r);
  }
  SECTION("45 degree incidence: s and p separately, the unpolarized source the mean") {
    // A wide thin plate (radius 100, 1 mm), so that the passes do not reach the mantle; the
    // planes of incidence of both faces coincide, s and p stay eigenstates, and the weight of
    // an unpolarized ray is the mean of the s and p powers (ADR 0021).
    const CompiledSystem wide = compiled(with_body(plate(1.0, 100.0)));
    const Scene ws = Scene::build(wide);
    const double i = std::numbers::pi / 4.0;
    const std::vector<RaySource> src{
        {{{Vec3(-4.0, 0.0, 1.0), Vec3(std::sin(i), 0.0, std::cos(i)), 1.0, 0}}}};
    const NonseqResult r = rtt::nonseq::trace(wide, ws, src, o);
    const auto [n0, n] = indices(wide);
    const auto [rs, rp] = fresnel_r(n0, n, i);
    const double t = 0.5 * ((1.0 - rs) / (1.0 + rs) + (1.0 - rp) / (1.0 + rp));
    CHECK(std::abs(r.balance[0].detected[face_of(ws, FaceKind::Detector)] - t) < 1e-13);
    check_closed(r);
  }
}

TEST_CASE("case 2: the same plate with Monte Carlo converges within 4 sigma", "[nonseq][trace]") {
  // Every Fresnel event loses nothing, so every ray keeps weight 1 and ends either at the
  // detector or escaped: the detected power is a binomial count, sigma = sqrt(T (1 - T) / N).
  const CompiledSystem system = compiled(with_body(plate(4.0, 12.7)));
  const Scene scene = Scene::build(system);
  constexpr std::size_t n = 20000;
  RaySource src;
  src.rays.assign(n, StartRay{Vec3(0.0, 0.0, -1.0), Vec3::UnitZ(), 1.0 / n, 0});
  TraceOptions o;
  o.seed = 20261127;
  const NonseqResult r = rtt::nonseq::trace(system, scene, std::vector<RaySource>{src}, o);
  const double R = normal_r(system);
  const double t = (1.0 - R) / (1.0 + R);
  const double sigma = std::sqrt(t * (1.0 - t) / n);
  const double got = r.balance[0].detected[face_of(scene, FaceKind::Detector)];
  CHECK(std::abs(got - t) < 4.0 * sigma);
  CHECK(r.balance[0].escaped > 0.0);  // guard: both branches were chosen
  CHECK(r.rays == n);
  check_closed(r);
}

TEST_CASE("case 5: light guide, TIR at the mantle carries every steep enough ray",
          "[nonseq][trace]") {
  // A rod of CONST:1.5, radius 2, length 50, ideal_anti_reflection at both ends; the detector
  // sits at the exit face (1e-6 mm behind it, radius 2 + 1e-5 mm) (Nachtrag Review PR 1, f). A
  // ray inside at angle theta to the axis meets the mantle at 90 deg - theta: beyond the
  // critical angle asin(1/1.5) = 41.8 deg (sin_t^2 > 1, de Greve, Eq. (24)) it is totally
  // reflected without loss; the end face at theta < 41.8 deg transmits fully. For theta = 70 deg
  // the mantle angle 20 deg is below it, and most of the power leaves through the mantle.
  Element rod = plate(50.0, 2.0);
  for (auto& s : rod.surfaces) s.interaction = rtt::model::IdealAntiReflection{};
  System sys = with_body(rod);
  Surface det;
  det.id = SurfaceId("D");
  det.aperture = CircularAperture{2.0 + 1e-5, 0.0};
  rtt::model::test::element(sys, 1).pose = Pose::along_z(55.0 + 1e-6);
  rtt::model::test::element(sys, 1).surfaces = {det};
  const CompiledSystem system = compiled(sys);
  const Scene scene = Scene::build(system);
  const std::uint32_t d = face_of(scene, FaceKind::Detector);
  TraceOptions o;
  o.splitting = Splitting::Full;
  for (const double deg : {10.0, 20.0, 30.0, 40.0}) {
    const double th = deg * std::numbers::pi / 180.0;
    const std::vector<RaySource> src{
        {{{Vec3(0.0, 0.0, 6.0), Vec3(std::sin(th), 0.0, std::cos(th)), 1.0, 0}}}};
    const NonseqResult r = rtt::nonseq::trace(system, scene, src, o);
    // Every reflection and the AR exit keep the weight up to rounding (a few 1e-16 per event).
    CHECK(std::abs(r.balance[0].detected[d] - 1.0) < 1e-12);
    CHECK(r.events >= 2);
    check_closed(r);
  }
  // 45 deg: totally reflected at the mantle and at the AR end faces too, so the ray stays
  // trapped until max_depth cuts it. At TIR the ideal anti-reflection surface reflects like a
  // Fresnel surface (ADR 0057, Nachtrag PR 3, b): nothing is absorbed, everything is cut.
  {
    const double t45 = std::numbers::pi / 4.0;
    const std::vector<RaySource> src{
        {{{Vec3(0.0, 0.0, 6.0), Vec3(std::sin(t45), 0.0, std::cos(t45)), 1.0, 0}}}};
    const NonseqResult r = rtt::nonseq::trace(system, scene, src, o);
    CHECK(sum(r.balance[0].absorbed) < 1e-12);
    CHECK(std::abs(r.balance[0].cut - 1.0) < 1e-12);
    check_closed(r);
  }
  const double th = 70.0 * std::numbers::pi / 180.0;
  const std::vector<RaySource> steep{
      {{{Vec3(0.0, 0.0, 6.0), Vec3(std::sin(th), 0.0, std::cos(th)), 1.0, 0}}}};
  const NonseqResult r = rtt::nonseq::trace(system, scene, steep, o);
  CHECK(r.balance[0].detected[d] < 1e-3);
  CHECK(r.balance[0].escaped > 0.9);
  check_closed(r);
}

TEST_CASE("case 7: balance closes and the result is bitwise equal for 1, 2 and 8 threads",
          "[nonseq][trace]") {
  // Three blocks of 4096 start rays with Monte Carlo on the plate, two wavelengths of the same
  // index (CONST) sharing the balance by index; the sums run in block order (point 5).
  System sys = with_body(plate(4.0, 12.7));
  sys.wavelengths.push_back({0.6563, 1.0, false});
  const CompiledSystem system = compiled(sys);
  const Scene scene = Scene::build(system);
  RaySource src;
  for (std::size_t i = 0; i < 10000; ++i) {
    const double x = -6.0 + 12.0 * static_cast<double>(i % 100) / 99.0;
    const double y = -6.0 + 12.0 * std::floor(static_cast<double>(i) / 100.0) / 99.0;
    src.rays.push_back(
        {Vec3(x, y, -1.0), Vec3(0.02 * x, 0.0, 1.0), 1e-4, static_cast<std::uint16_t>(i % 2)});
  }
  TraceOptions o;
  o.seed = 7;
  const std::vector<RaySource> sources{src};
  const auto run = [&](int threads) {
    oneapi::tbb::task_arena arena(threads);
    return arena.execute([&] { return rtt::nonseq::trace(system, scene, sources, o); });
  };
  const NonseqResult a = run(1);
  const NonseqResult b = run(2);
  const NonseqResult c = run(8);
  REQUIRE(a.balance.size() == 2);
  for (const NonseqResult* other : {&b, &c}) {
    CHECK(other->events == a.events);
    for (std::size_t k = 0; k < a.balance.size(); ++k) {
      CHECK(other->balance[k].source == a.balance[k].source);
      CHECK(other->balance[k].volume == a.balance[k].volume);
      CHECK(other->balance[k].not_modelled == a.balance[k].not_modelled);
      CHECK(other->balance[k].passed == a.balance[k].passed);
      CHECK(other->balance[k].escaped == a.balance[k].escaped);
      CHECK(other->balance[k].detected == a.balance[k].detected);
      CHECK(other->balance[k].absorbed == a.balance[k].absorbed);
      CHECK(other->balance[k].cut == a.balance[k].cut);
    }
  }
  CHECK(a.rays == 10000);
  CHECK(sum(a.balance[0].detected) > 0.3);  // guard: the plate transmits most of the light
  check_closed(a);
}

TEST_CASE("mirror: the front reflects, the back absorbs", "[nonseq][trace]") {
  Surface m;
  m.id = SurfaceId("M");
  m.aperture = CircularAperture{5.0, 0.0};
  const System sys = with_body(Element{"M", ElementKind::Mirror, Pose::along_z(5.0), {}, {m}});
  const CompiledSystem system = compiled(sys);
  const Scene scene = Scene::build(system);
  const std::uint32_t mirror = face_of(scene, FaceKind::Mirror);
  // Front: coming along +z from the -z side; an ideal conductor reflects everything back.
  const std::vector<RaySource> front{{{{Vec3(0.0, 0.0, 0.0), Vec3::UnitZ(), 1.0, 0}}}};
  const NonseqResult f = rtt::nonseq::trace(system, scene, front, {});
  CHECK(std::abs(f.balance[0].escaped - 1.0) < 1e-15);
  // Back: coming along -z from the +z side.
  const std::vector<RaySource> back{{{{Vec3(0.0, 0.0, 10.0), -Vec3::UnitZ(), 1.0, 0}}}};
  const NonseqResult k = rtt::nonseq::trace(system, scene, back, {});
  CHECK(k.balance[0].absorbed[mirror] == 1.0);
  check_closed(f);
  check_closed(k);
}

TEST_CASE("thin ideal beam splitter splits in the environment (ADR 0057, Nachtrag PR 3, a)",
          "[nonseq][trace]") {
  Surface s;
  s.id = SurfaceId("BS");
  s.aperture = CircularAperture{5.0, 0.0};
  s.interaction = rtt::model::IdealBeamSplitter{0.3, 0.3};
  const System sys =
      with_body(Element{"BS", ElementKind::ThinElement, Pose::along_z(5.0), {}, {s}});
  const CompiledSystem system = compiled(sys);
  const Scene scene = Scene::build(system);
  TraceOptions o;
  o.splitting = Splitting::Full;
  const std::vector<RaySource> src{{{{Vec3(0.0, 0.0, 0.0), Vec3::UnitZ(), 1.0, 0}}}};
  const NonseqResult r = rtt::nonseq::trace(system, scene, src, o);
  CHECK(std::abs(r.balance[0].detected[face_of(scene, FaceKind::Detector)] - 0.7) < 1e-15);
  CHECK(std::abs(r.balance[0].escaped - 0.3) < 1e-15);
  check_closed(r);
}

TEST_CASE("transparent detector records both directions and lets the ray pass", "[nonseq][trace]") {
  // The plate of case 1 with a transparent detector in front of it at z = 2: it sees the
  // incoming power 1 and the reflected power 2R/(1 + R) on its way back.
  System sys = with_body(plate(4.0, 12.7));
  Surface d;
  d.id = SurfaceId("T");
  d.aperture = CircularAperture{20.0, 0.0};
  sys.root.children.insert(
      sys.root.children.begin(),
      Node{Element{"front", ElementKind::Detector, Pose::along_z(2.0), {}, {d}}});
  const CompiledSystem system = compiled(sys);
  rtt::nonseq::SceneOptions so;
  so.transparent_detectors = {0};
  const Scene scene = Scene::build(system, so);
  TraceOptions o;
  o.splitting = Splitting::Full;
  o.min_relative_power = 0.0;
  const std::vector<RaySource> src{{{{Vec3(0.0, 0.0, -1.0), Vec3::UnitZ(), 1.0, 0}}}};
  const NonseqResult r = rtt::nonseq::trace(system, scene, src, o);
  // The plate is element 1 here (the front detector is element 0).
  const double n0 = system.media()[system.environment_medium()].index[0].real();
  const double n = system.media()[system.elements()[1].media[0]].index[0].real();
  const double R = (n - n0) * (n - n0) / ((n + n0) * (n + n0));
  CHECK(std::abs(r.balance[0].passed[0] - (1.0 + 2.0 * R / (1.0 + R))) < 1e-14);
  CHECK(std::abs(r.balance[0].escaped - 2.0 * R / (1.0 + R)) < 1e-14);
  check_closed(r);
}

TEST_CASE("crystal plate: entry into o and e, exit by Refract only", "[nonseq][trace]") {
  // ADR 0057, point 4, and Nachtrag Review PR 1, b: no reflection at entry (projection model,
  // ADR 0026) and none inside; at normal incidence both modes leave through the back face.
  Element p = plate(4.0, 12.7);
  p.material.reset();
  p.crystal = rtt::model::CrystalMaterial{"CONST:1.6584", "CONST:1.4864"};
  p.optic_axis = std::array<double, 3>{1.0, 0.0, 1.0};
  System sys = with_body(p);
  // compile rejects an automatic path through a crystal; the scene does not use the paths.
  sys.paths = {{"det", false, {{SurfaceId("IMG"), rtt::model::EventKind::Transmit}}}};
  const CompiledSystem system = compiled(sys);
  const Scene scene = Scene::build(system);
  TraceOptions o;
  o.splitting = Splitting::Full;
  const std::vector<RaySource> src{{{{Vec3(0.0, 0.0, -1.0), Vec3::UnitZ(), 1.0, 0}}}};
  const NonseqResult r = rtt::nonseq::trace(system, scene, src, o);
  const Balance& b = r.balance[0];
  CHECK(b.not_modelled == 0.0);
  CHECK(r.events == 3);  // entry with two branches (one event), one exit per branch
  // The projection model has no loss at entry or exit: everything arrives (rounding).
  CHECK(std::abs(b.detected[face_of(scene, FaceKind::Detector)] - 1.0) < 1e-14);
  check_closed(r);
}

TEST_CASE("trace rejects invalid input", "[nonseq][trace]") {
  const CompiledSystem system = compiled(with_body(plate(4.0, 12.7)));
  const Scene scene = Scene::build(system);
  const auto bad = [&](StartRay ray, TraceOptions o = {}) {
    const std::vector<RaySource> src{{{ray}}};
    CHECK_THROWS_AS(rtt::nonseq::trace(system, scene, src, o), std::invalid_argument);
  };
  bad({Vec3::Zero(), Vec3::UnitZ(), 1.0, 1});   // wavelength out of range
  bad({Vec3::Zero(), Vec3::UnitZ(), -1.0, 0});  // negative power
  bad({Vec3::Zero(), Vec3::Zero(), 1.0, 0});    // zero direction
  TraceOptions zero_depth;
  zero_depth.max_depth = 0;
  bad({Vec3::Zero(), Vec3::UnitZ(), 1.0, 0}, zero_depth);
}

TEST_CASE("a cancelled run returns the finished blocks without throwing", "[nonseq][trace]") {
  const CompiledSystem system = compiled(with_body(plate(4.0, 12.7)));
  const Scene scene = Scene::build(system);
  RaySource src;
  src.rays.assign(5000, StartRay{Vec3(0.0, 0.0, -1.0), Vec3::UnitZ(), 1.0, 0});
  rtt::trace::RunControl control;
  control.cancel = rtt::trace::CancelToken{};
  control.cancel->request_cancel();  // before the first block: nothing is traced
  const NonseqResult r =
      rtt::nonseq::trace(system, scene, std::vector<RaySource>{src}, {}, control);
  CHECK(r.cancelled);
  CHECK(r.rays == 0);
  CHECK(r.balance[0].source == 0.0);
}

TEST_CASE("absorbing plate: volume absorption in the balance", "[nonseq][trace]") {
  // CONST:1.5 with kappa = 2e-5, 4 mm, normal incidence, full splitting. Each pass transmits
  // tau = exp(-4 pi kappa d / lambda) (Byrnes, Eqs. (1), (2), as move_to_hit), so the incoherent
  // sum gives T = (1 - R)^2 tau / (1 - R^2 tau^2) and R_tot = R + (1 - R)^2 R tau^2 /
  // (1 - R^2 tau^2) (the series of Whittaker and Gehring with tau per pass), R = |r|^2 with the
  // complex index. The formula takes 1 - R as the power into the absorbing medium; the
  // tracer's interface powers (ADR 0021) differ from it at order kappa^2 = 4e-10, hence 1e-9.
  constexpr double kappa = 2e-5;
  const CompiledSystem system = compiled(with_body(plate(4.0, 12.7, "CONST:1.5,2e-5")));
  const Scene scene = Scene::build(system);
  REQUIRE(system.media()[system.elements()[0].media[0]].index[0].imag() == kappa);
  TraceOptions o;
  o.splitting = Splitting::Full;
  o.min_relative_power = 0.0;
  const std::vector<RaySource> src{{{{Vec3(0.0, 0.0, -1.0), Vec3::UnitZ(), 1.0, 0}}}};
  const NonseqResult r = rtt::nonseq::trace(system, scene, src, o);
  const double lambda_mm = system.wavelengths_um()[0] * 1e-3;
  const double n0 = system.media()[system.environment_medium()].index[0].real();
  const auto nn = system.media()[system.elements()[0].media[0]].index[0];
  const double R = std::norm((nn - n0) / (nn + n0));
  const double tau = std::exp(-4.0 * std::numbers::pi * kappa * 4.0 / lambda_mm);
  const double t = (1.0 - R) * (1.0 - R) * tau / (1.0 - R * R * tau * tau);
  const double r_tot = R + (1.0 - R) * (1.0 - R) * R * tau * tau / (1.0 - R * R * tau * tau);
  const Balance& b = r.balance[0];
  CHECK(std::abs(b.detected[face_of(scene, FaceKind::Detector)] - t) < 1e-9);
  CHECK(std::abs(b.escaped - r_tot) < 1e-9);
  CHECK(std::abs(b.volume - (1.0 - t - r_tot)) < 1e-9);
  CHECK(b.volume > 0.5);  // guard: most of the light is absorbed in the glass
  check_closed(r);
}

TEST_CASE("mirror with substrate: Fresnel against the substrate, the rest absorbed",
          "[nonseq][trace]") {
  // ADR 0057, Nachtrag PR 3, e: the front of a mirror with material reflects with Fresnel
  // against its substrate, R = ((n - n0)/(n + n0))^2 at normal incidence (Byrnes, as above);
  // the substrate is no body, the transmitted rest is absorbed at the face.
  Surface m;
  m.id = SurfaceId("M");
  m.aperture = CircularAperture{5.0, 0.0};
  const CompiledSystem system =
      compiled(with_body(Element{"M", ElementKind::Mirror, Pose::along_z(5.0), "CONST:1.5", {m}}));
  const Scene scene = Scene::build(system);
  const double n0 = system.media()[system.environment_medium()].index[0].real();
  const double n = system.media()[system.elements()[0].media[0]].index[0].real();
  const double R = (n - n0) * (n - n0) / ((n + n0) * (n + n0));
  const std::vector<RaySource> src{{{{Vec3(0.0, 0.0, 0.0), Vec3::UnitZ(), 1.0, 0}}}};
  TraceOptions o;
  o.splitting = Splitting::Full;
  const NonseqResult r = rtt::nonseq::trace(system, scene, src, o);
  CHECK(std::abs(r.balance[0].escaped - R) < 1e-15);
  CHECK(std::abs(r.balance[0].absorbed[face_of(scene, FaceKind::Mirror)] - (1.0 - R)) < 1e-15);
  check_closed(r);
}

TEST_CASE("thin grating with only order 0 listed: its efficiency applies", "[nonseq][trace]") {
  // ADR 0057, point 4: every listed order is a branch; a pass between equal media without
  // event would skip the efficiency 0.5 of the single listed order 0.
  Surface g;
  g.id = SurfaceId("G");
  g.aperture = CircularAperture{5.0, 0.0};
  g.phases = {rtt::model::LinearGrating{rtt::model::Param(300.0), 90.0}};
  g.diffraction_efficiency = std::vector<rtt::model::DiffractionEfficiency>{{0, 0.5}};
  const CompiledSystem system =
      compiled(with_body(Element{"G", ElementKind::ThinElement, Pose::along_z(5.0), {}, {g}}));
  const Scene scene = Scene::build(system);
  TraceOptions o;
  o.splitting = Splitting::Full;
  const std::vector<RaySource> src{{{{Vec3(0.0, 0.0, 0.0), Vec3::UnitZ(), 1.0, 0}}}};
  const NonseqResult r = rtt::nonseq::trace(system, scene, src, o);
  CHECK(std::abs(r.balance[0].detected[face_of(scene, FaceKind::Detector)] - 0.5) < 1e-15);
  CHECK(r.events == 1);
  check_closed(r);
}

TEST_CASE("crystal rod with AR ends: TIR inside the crystal is not modelled", "[nonseq][trace]") {
  // ADR 0057, Nachtrag Review PR 1, b: inside a crystal only the exit; at TIR the ray ends as
  // not modelled. The ray enters through the mantle and meets the AR end face beyond the
  // critical angle; the TIR of the AR surface (Nachtrag PR 3, b) must not apply inside a
  // crystal (before the fix this read past the end of the branch list).
  Element p = plate(50.0, 2.0);
  for (auto& s : p.surfaces) s.interaction = rtt::model::IdealAntiReflection{};
  p.material.reset();
  p.crystal = rtt::model::CrystalMaterial{"CONST:1.6584", "CONST:1.4864"};
  p.optic_axis = std::array<double, 3>{1.0, 0.0, 1.0};
  System sys = with_body(p);
  sys.paths = {{"det", false, {{SurfaceId("IMG"), rtt::model::EventKind::Transmit}}}};
  const CompiledSystem system = compiled(sys);
  const Scene scene = Scene::build(system);
  TraceOptions o;
  o.splitting = Splitting::Full;
  const std::vector<RaySource> src{
      {{{Vec3(-3.0, 0.0, 53.0), Vec3(1.0, 0.0, 1.0).normalized(), 1.0, 0}}}};
  const NonseqResult r = rtt::nonseq::trace(system, scene, src, o);
  // Both modes end at the end face by TIR: entry (one event) and one TIR per mode.
  CHECK(std::abs(r.balance[0].not_modelled - 1.0) < 1e-14);
  CHECK(r.events == 3);
  check_closed(r);
}
