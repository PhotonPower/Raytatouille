// OPD at the points of an arbitrary pupil sampling (#168, ADR 0030 addendum): the conventions of
// opd_map (reference sphere, chief ray of the same wavelength, waves at the reference
// wavelength), without statistics, so that every sampling is allowed, also the Gaussian one.

#include <algorithm>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/trace/run_control.hpp"
#include "rtt/trace/sources.hpp"

using rtt::analysis::OpdOptions;
using rtt::analysis::OpdPoint;
using rtt::analysis::OpdPupilPoints;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::trace::RayStatus;

namespace {

/// tests/reference/m1/singlet_const.rtt.json: a plano-convex CONST singlet, fields 0, 1, 2;
/// `lens_radius` > 0 replaces the aperture of the first lens surface (12.7 mm).
CompiledSystem singlet(double lens_radius = 0.0) {
  rtt::model::System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  if (lens_radius > 0.0) {
    std::get<rtt::model::Element>(s.root.children[1].value).surfaces[0].aperture =
        rtt::model::CircularAperture{lens_radius, 0.0};
  }
  return rtt::compile::compile(s, rtt::material::MaterialLibrary{});
}

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

bool same(const OpdPoint& a, const OpdPoint& b) {
  return same(a.px, b.px) && same(a.py, b.py) && same(a.w, b.w) && a.status == b.status;
}

}  // namespace

TEST_CASE("opd_points: a grid gives the points of opd_map bit for bit", "[opd][points]") {
  // opd_map runs on the same function: same sphere, same points, same losses.
  const CompiledSystem cs = singlet();
  OpdOptions options;
  options.grid = 9;
  for (const std::uint16_t field : {std::uint16_t{0}, std::uint16_t{2}}) {
    INFO("field " << field);
    const rtt::analysis::OpdMap map = rtt::analysis::opd_map(cs, PathId{0}, field, 0, options);
    const OpdPupilPoints p =
        rtt::analysis::opd_points(cs, PathId{0}, field, 0, rtt::trace::GridPupil{9}, options);
    CHECK(p.field == field);
    CHECK(p.wavelength == 0);
    REQUIRE(p.points.size() == map.points.size());
    for (std::size_t i = 0; i < p.points.size(); ++i) CHECK(same(p.points[i], map.points[i]));
    CHECK(same(p.sphere.radius, map.sphere.radius));
    CHECK(p.sphere.centre == map.sphere.centre);
    CHECK(p.losses.launched == map.losses.launched);
    if (field == 2) {
      // Content: off axis the wavefront is not flat.
      double largest = 0.0;
      for (const OpdPoint& q : p.points) largest = std::max(largest, std::abs(q.w));
      CHECK(largest > 1e-3);
    }
  }
}

TEST_CASE("opd_points: the gauss sampling, point by point as single rays", "[opd][points]") {
  // Every point of GaussPupil gives the W of the same point aimed alone (SinglePupilPoint): the
  // chief ray and the reference sphere do not depend on the sampling.
  const CompiledSystem cs = singlet();
  const rtt::trace::GaussPupil gauss{3, 6};
  const OpdPupilPoints p = rtt::analysis::opd_points(cs, PathId{0}, 2, 0, gauss);
  const std::vector<rtt::trace::PupilPoint> pupil = rtt::trace::pupil_points(gauss);
  REQUIRE(p.points.size() == pupil.size());
  bool nonzero = false;
  for (std::size_t i = 0; i < pupil.size(); ++i) {
    INFO("point " << i);
    const OpdPupilPoints single = rtt::analysis::opd_points(
        cs, PathId{0}, 2, 0, rtt::trace::SinglePupilPoint{pupil[i].px, pupil[i].py});
    REQUIRE(single.points.size() == 1);
    CHECK(same(p.points[i], single.points[0]));
    CHECK(p.points[i].status == RayStatus::Alive);
    nonzero = nonzero || p.points[i].w != 0.0;
  }
  CHECK(nonzero);
}

TEST_CASE("opd_points: lost rays keep W = 0 and their status, without an exception",
          "[opd][points]") {
  // An aperture of 5 mm on L1.S1 (beam radius 10 mm, collimated, 5 mm behind the stop) lets the
  // chief ray (the reference of W, not a gauss point) and the inner ring (rho = 0.34) arrive and
  // stops the outer rings (rho = 0.71, 0.94); opd_points has no statistics and reports them as
  // points.
  const CompiledSystem cs = singlet(5.0);
  const OpdPupilPoints p =
      rtt::analysis::opd_points(cs, PathId{0}, 0, 0, rtt::trace::GaussPupil{3, 6});
  std::size_t lost = 0;
  for (const OpdPoint& q : p.points) {
    if (q.status == RayStatus::Alive) continue;
    ++lost;
    CHECK(q.w == 0.0);
  }
  CHECK(lost == 12);
  CHECK(p.points.size() == 18);
  CHECK(p.losses.launched == p.points.size());
}

TEST_CASE("opd_points: run control gives the same points; input errors", "[opd][points]") {
  const CompiledSystem cs = singlet();
  const rtt::trace::GaussPupil gauss{2, 6};
  const OpdPupilPoints plain = rtt::analysis::opd_points(cs, PathId{0}, 1, 0, gauss);
  rtt::trace::RunControl control;
  control.cancel = rtt::trace::CancelToken{};
  const OpdPupilPoints controlled =
      rtt::analysis::opd_points(cs, PathId{0}, 1, 0, gauss, OpdOptions{}, control);
  REQUIRE(plain.points.size() == controlled.points.size());
  for (std::size_t i = 0; i < plain.points.size(); ++i) {
    CHECK(same(plain.points[i], controlled.points[i]));
  }
  CHECK_THROWS_AS(rtt::analysis::opd_points(cs, PathId{5}, 0, 0, gauss), std::invalid_argument);
  CHECK_THROWS_AS(rtt::analysis::opd_points(cs, PathId{0}, 0, 9, gauss), std::invalid_argument);
  CHECK_THROWS_AS(rtt::analysis::opd_points(cs, PathId{0}, 0, 0, rtt::trace::GaussPupil{0, 6}),
                  std::invalid_argument);
}
