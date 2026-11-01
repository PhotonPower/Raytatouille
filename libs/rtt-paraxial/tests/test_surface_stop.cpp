// Paraxial pupils with the stop on a surface (#53, ADR 0039, point 4): the stop plane is the vertex
// plane of the stop surface; the entrance pupil images it through the steps before the surface,
// the exit pupil through the surface itself and the steps after it (the stop lies in front of
// the refraction at its own surface).
//
// Systems: m2/cooke_triplet_surface_stop (stop on L2.S1) and the workaround of #34,
// m2/cooke_triplet (a plane stop element delta = 1e-8 mm behind the vertex of L2.S1, in the
// glass, on an explicit path).

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <map>
#include <string>
#include <variant>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::SurfaceId;
using rtt::model::System;

namespace {

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

const MaterialLibrary& schott_m2() {
  struct Lib {  // MaterialLibrary is not movable
    MaterialLibrary lib;
    Lib() { lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf"); }
  };
  static const Lib l;
  return l.lib;
}

const CompiledSystem& surface_stop() {
  static const CompiledSystem cs =
      rtt::compile::compile(load("m2/cooke_triplet_surface_stop.rtt.json"), schott_m2());
  return cs;
}

const CompiledSystem& workaround() {
  static const CompiledSystem cs =
      rtt::compile::compile(load("m2/cooke_triplet.rtt.json"), schott_m2());
  return cs;
}

/// Image of the stop through L2.S1, L2.S2, L3.S1, L3.S2 by an independent y-nu trace
/// (Greivenkamp, OPTI-201/202, Sec. 9, p. 9-2, docs/quellen.md: refraction n'u' = nu - y phi with
/// phi = c (n' - n), transfer y' = y + t u'), starting at
/// the vertex plane of L2.S1 in the medium in front of it. Radii and vertices from the file
/// (L2 at z = 10.691, L2.S2 0.975 behind; L3 at z = 16.488, L3.S2 3.127 behind), indices from
/// the compiled media of the path at the reference wavelength. `refract_at_stop` = false is the
/// content guard, the trace of the code before ADR 0039, point 4: the backward trace started
/// behind the stop event, as for a stop element, but in the medium in front of the stop (air), so
/// it left out the refraction at L2.S1 and took air for the glass of L2.
struct ExitPupil {
  double z = 0.0;
  double diameter = 0.0;
};

ExitPupil independent_exit_pupil(const CompiledSystem& cs, bool refract_at_stop) {
  const std::uint16_t wl = cs.reference_wavelength();
  const auto& events = cs.path(PathId{0}).events;  // L1.S1, L1.S2, L2.S1, L2.S2, L3.S1, L3.S2, IMG
  auto n = [&](std::uint32_t medium) { return cs.media()[medium].index[wl].real(); };
  const std::array<double, 4> radius{-24.456, 21.896, 86.759, -20.4942};
  const std::array<double, 4> z{10.691, 10.691 + 0.975, 16.488, 16.488 + 3.127};
  const std::array<double, 4> n_before{
      n(events[2].medium_before),
      refract_at_stop ? n(events[3].medium_before) : n(events[2].medium_before),
      n(events[4].medium_before), n(events[5].medium_before)};
  const std::array<double, 4> n_after{n(events[2].medium_after), n(events[3].medium_after),
                                      n(events[4].medium_after), n(events[5].medium_after)};
  // Two rays from the stop plane: the chief ray from the centre (y = 0, u = 0.1) and the rim ray
  // (y = r = 10, u = 0), r the radius of the aperture of L2.S1.
  struct Ray {
    double y;
    double u;
  };
  Ray chief{0.0, 0.1};
  Ray rim{10.0, 0.0};
  for (std::size_t k = 0; k < 4; ++k) {
    if (k > 0) {
      const double t = z[k] - z[k - 1];
      chief.y += t * chief.u;
      rim.y += t * rim.u;
    }
    if (k == 0 && !refract_at_stop) continue;
    const double phi = (n_after[k] - n_before[k]) / radius[k];
    chief.u = (n_before[k] * chief.u - chief.y * phi) / n_after[k];
    rim.u = (n_before[k] * rim.u - rim.y * phi) / n_after[k];
  }
  // The exit pupil lies where the chief ray crosses the axis; its radius is the rim ray's height
  // there (the rim ray leaves the stop rim parallel to the axis, so it marks the image of the rim).
  const double z_xp = z[3] - chief.y / chief.u;
  const double h = rim.y + (z_xp - z[3]) * rim.u;
  return {z_xp, 2.0 * std::abs(h)};
}

/// Paraxial stop radius R_s of the triplet with EPD 12.5 mm and the object at infinity: the
/// marginal ray enters parallel at y = 6.25 mm, passes the rim of the entrance pupil and therefore
/// the rim of the stop (conjugate planes). Independent y-nu trace through L1 (L1.S1 at z = 0,
/// R = 23.713; L1.S2 4.831 behind, R = 7331.288) to the vertex plane of L2.S1 (z = 10.691), same
/// equations as above.
double independent_stop_radius(const CompiledSystem& cs) {
  const std::uint16_t wl = cs.reference_wavelength();
  const auto& events = cs.path(PathId{0}).events;
  auto n = [&](std::uint32_t medium) { return cs.media()[medium].index[wl].real(); };
  const std::array<double, 2> radius{23.713, 7331.288};
  const std::array<double, 2> z{0.0, 4.831};
  double y = 6.25;
  double u = 0.0;
  for (std::size_t k = 0; k < 2; ++k) {
    if (k > 0) y += (z[k] - z[k - 1]) * u;
    const double nb = n(events[k].medium_before);
    const double na = n(events[k].medium_after);
    u = (nb * u - y * (na - nb) / radius[k]) / na;
  }
  return std::abs(y + (10.691 - z[1]) * u);
}

}  // namespace

TEST_CASE("surface stop: first order as with the stop element", "[paraxial][aperture-stop]") {
  // f, BFL and FFL do not depend on the stop. The workaround splits the gap L2.S1 -> L2.S2 at
  // STO into two transfers, which may round differently: tolerance 1e-12 relative (a few ulp of
  // a ~50 mm value per step, seven steps).
  const std::uint16_t wl = surface_stop().reference_wavelength();
  const auto a = rtt::paraxial::first_order(surface_stop(), PathId{0}, wl);
  const auto b = rtt::paraxial::first_order(workaround(), PathId{0}, wl);
  REQUIRE(a.rear_focal_length);
  REQUIRE(a.bfl);
  REQUIRE(a.ffl);
  REQUIRE(b.rear_focal_length);
  REQUIRE(b.bfl);
  REQUIRE(b.ffl);
  CHECK(std::abs(*a.rear_focal_length - *b.rear_focal_length) <= 1e-12 * *b.rear_focal_length);
  CHECK(std::abs(*a.bfl - *b.bfl) <= 1e-12 * std::abs(*b.bfl));
  CHECK(std::abs(*a.ffl - *b.ffl) <= 1e-12 * std::abs(*b.ffl));
}

TEST_CASE("surface stop: entrance pupil of the triplet against #34", "[paraxial][aperture-stop]") {
  // Independent y-nu trace of the stop centre through L2.S1 and L1 (Greivenkamp, OPTI-201/202
  // Sec. 9, p. 9-2; Python script in the PR of #34, without rtt): z_EP = 11.679209257540240 mm
  // with the stop in the vertex plane of L2.S1. Tolerance 1e-9 mm (double rounding ~1e-14 mm).
  // Content guard: the workaround (stop 1e-8 mm behind the vertex) gives 11.679209268268448 mm,
  // 1.07e-8 mm away, so a stop in the wrong plane fails.
  const std::uint16_t wl = surface_stop().reference_wavelength();
  const auto fo = rtt::paraxial::first_order(surface_stop(), PathId{0}, wl);
  REQUIRE(fo.entrance_pupil);
  REQUIRE(fo.entrance_pupil->z);
  REQUIRE(fo.entrance_pupil->diameter);
  const double z_ep = *fo.entrance_pupil->z;
  INFO(std::setprecision(17) << "z_EP " << z_ep);
  CHECK(std::abs(z_ep - 11.679209257540240) <= 1e-9);
  CHECK(std::abs(z_ep - 11.679209268268448) > 5e-9);
  CHECK(std::abs(*fo.entrance_pupil->diameter - 12.5) <= 1e-12);
}

TEST_CASE("surface stop: exit pupil through the stop surface itself", "[paraxial][aperture-stop]") {
  // Against the independent traces above, tolerance 1e-9 mm (double y-nu traces over at most
  // four surfaces, rounding ~1e-13 mm).
  const CompiledSystem& cs = surface_stop();
  const std::uint16_t wl = cs.reference_wavelength();
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, wl);
  REQUIRE(fo.exit_pupil);
  REQUIRE(fo.exit_pupil->z);
  REQUIRE(fo.exit_pupil->diameter);
  REQUIRE(fo.entrance_pupil);
  REQUIRE(fo.entrance_pupil->diameter);
  const ExitPupil xp = independent_exit_pupil(cs, true);
  INFO(std::setprecision(17) << "z_XP " << *fo.exit_pupil->z << " (independent " << xp.z
                             << "), D_XP " << *fo.exit_pupil->diameter);
  CHECK(std::abs(*fo.exit_pupil->z - xp.z) <= 1e-9);
  // The pupils are images of the paraxial stop radius R_s (EPD 12.5 mm), not of the aperture
  // radius 10: D_XP = xp.diameter * R_s / 10, both traces linear in the rim height.
  const double r_s = independent_stop_radius(cs);
  INFO(std::setprecision(17) << "R_s " << r_s);
  CHECK(std::abs(*fo.exit_pupil->diameter - xp.diameter * r_s / 10.0) <= 1e-9);
  // Content guard: the trace of the old code (see independent_exit_pupil) puts the exit pupil
  // clearly elsewhere. A refraction in the stop plane alone would change neither position nor size
  // (it only turns rays at the points of the stop, like a field lens in an image plane); the
  // old code went wrong by continuing with the index in front of the stop.
  const ExitPupil wrong = independent_exit_pupil(cs, false);
  INFO(std::setprecision(17) << "old code: z_XP " << wrong.z);
  CHECK(std::abs(wrong.z - xp.z) > 1e-3);
}

TEST_CASE("surface stop: Seidel sums against the workaround", "[paraxial][aperture-stop]") {
  // The workaround's stop sits u_bar delta higher on the chief ray; the stop shifting parameter
  // gives |dS| < 0.65e-9 mm for every value of every surface (derivation in
  // rtt-analysis/tests/test_m2_acceptance.cpp, comment "Stop:"). Plus 1e-12 mm for rounding.
  const std::uint16_t wl = surface_stop().reference_wavelength();
  const auto a = rtt::paraxial::seidel(surface_stop(), PathId{0}, wl);
  const auto b = rtt::paraxial::seidel(workaround(), PathId{0}, wl);
  std::map<std::string, rtt::paraxial::SeidelTerms> by_id;
  for (const auto& s : b.surfaces) by_id[workaround().surfaces()[s.surface].id.str()] = s.terms;
  std::size_t compared = 0;
  for (const auto& s : a.surfaces) {
    const std::string id = surface_stop().surfaces()[s.surface].id.str();
    INFO(id);
    REQUIRE(by_id.contains(id));
    const auto& t = by_id[id];
    constexpr double kBound = 0.65e-9 + 1e-12;
    CHECK(std::abs(s.terms.s1 - t.s1) <= kBound);
    CHECK(std::abs(s.terms.s2 - t.s2) <= kBound);
    CHECK(std::abs(s.terms.s3 - t.s3) <= kBound);
    CHECK(std::abs(s.terms.s4 - t.s4) <= kBound);
    CHECK(std::abs(s.terms.s5 - t.s5) <= kBound);
    ++compared;
  }
  CHECK(compared == 7);  // six refracting surfaces and IMG; STO only in the workaround
}

TEST_CASE("surface stop: a mirror as the stop", "[paraxial][aperture-stop]") {
  // m2/paraboloid_stop without its stop element, the stop on the mirror M (vertex at z = 0):
  // nothing lies in front of the stop, so the entrance pupil is the stop plane itself, and the
  // exit pupil is its image through the mirror at zero distance, the same plane.
  System s = load("m2/paraboloid_stop.rtt.json");
  std::erase_if(s.root.children, [](const rtt::model::Node& node) {
    const auto* e = std::get_if<rtt::model::Element>(&node.value);
    return e != nullptr && e->kind == rtt::model::ElementKind::Stop;
  });
  s.aperture.stop = SurfaceId("M");
  const CompiledSystem cs = rtt::compile::compile(s, MaterialLibrary{});
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, cs.reference_wavelength());
  REQUIRE(fo.entrance_pupil);
  REQUIRE(fo.entrance_pupil->z);
  REQUIRE(fo.entrance_pupil->diameter);
  CHECK(*fo.entrance_pupil->z == 0.0);
  CHECK(std::abs(*fo.entrance_pupil->diameter - 60.0) <= 1e-12);
  REQUIRE(fo.exit_pupil);
  REQUIRE(fo.exit_pupil->z);
  CHECK(std::abs(*fo.exit_pupil->z) <= 1e-12);
}

TEST_CASE("surface stop: stop_size takes the aperture of the stop surface",
          "[paraxial][aperture-stop]") {
  // ADR 0039, points 1 and 4: with stop_size the entrance pupil is the image of the aperture of
  // L2.S1 (radius 10 mm). Pupil and stop are conjugate and the map is linear in the height, so
  // D_EP = 12.5 mm * 10 / R_s with R_s the independent stop radius of the 12.5 mm pupil (above).
  System s = load("m2/cooke_triplet_surface_stop.rtt.json");
  s.aperture.type = rtt::model::SystemApertureType::StopSize;
  s.aperture.value = rtt::model::Param();
  const CompiledSystem cs = rtt::compile::compile(s, schott_m2());
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, cs.reference_wavelength());
  REQUIRE(fo.entrance_pupil);
  REQUIRE(fo.entrance_pupil->diameter);
  const double expected = 12.5 * 10.0 / independent_stop_radius(surface_stop());
  INFO(std::setprecision(17) << "D_EP " << *fo.entrance_pupil->diameter << ", expected "
                             << expected);
  CHECK(std::abs(*fo.entrance_pupil->diameter - expected) <= 1e-9);
  // Content guard: the aperture radius matters (10 mm, not the file's EPD 12.5 mm).
  CHECK(std::abs(*fo.entrance_pupil->diameter - 12.5) > 1.0);
}
