// Explicit pupil point lists in the analyses (ADR 0032, point 2; Refs #202): opd_points and spot
// give bitwise the same result for the list of a sampling's points as for the sampling itself.
// This is what the pupil function (ADR 0032, point 2) relies on.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <string>
#include <vector>

#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/sources.hpp"

using rtt::trace::PointListPupil;

namespace {

bool same_bits(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
}

/// Singlet with a constant index (no catalogue); field 1 is off axis, so the wavefront is not flat.
rtt::compile::CompiledSystem singlet() {
  return rtt::compile::compile(
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json"),
      rtt::material::MaterialLibrary{});
}

PointListPupil list_of(const rtt::trace::PupilSampling& sampling) {
  return PointListPupil{rtt::trace::pupil_points(sampling)};
}

}  // namespace

TEST_CASE("point list pupil: opd_points of the list equals opd_points of the grid bitwise",
          "[opd][point_list]") {
  const rtt::compile::CompiledSystem cs = singlet();
  const rtt::trace::GridPupil grid{9};
  const auto a = rtt::analysis::opd_points(cs, rtt::compile::PathId{0}, 1, 0, grid);
  const auto b = rtt::analysis::opd_points(cs, rtt::compile::PathId{0}, 1, 0, list_of(grid));
  REQUIRE(a.points.size() == b.points.size());
  // Content guard: the points arrive and the off-axis wavefront is not flat.
  std::size_t arrived = 0;
  double max_w = 0.0;
  for (std::size_t i = 0; i < a.points.size(); ++i) {
    INFO("point " << i);
    CHECK(same_bits(a.points[i].px, b.points[i].px));
    CHECK(same_bits(a.points[i].py, b.points[i].py));
    CHECK(same_bits(a.points[i].w, b.points[i].w));
    CHECK(a.points[i].status == b.points[i].status);
    if (a.points[i].status == rtt::trace::RayStatus::Alive) ++arrived;
    max_w = std::max(max_w, std::abs(a.points[i].w));
  }
  CHECK(arrived > a.points.size() / 2);
  CHECK(max_w > 1e-3);
  CHECK(same_bits(a.sphere.radius, b.sphere.radius));
}

TEST_CASE("point list pupil: spot accepts the list and equals the hexapolar spot bitwise",
          "[spot][point_list]") {
  // spot rejects only GaussPupil (its points need quadrature weights); a list is accepted.
  const rtt::compile::CompiledSystem cs = singlet();
  rtt::analysis::SpotOptions hexa;
  hexa.sampling = rtt::trace::HexapolarPupil{4};
  rtt::analysis::SpotOptions list;
  list.sampling = list_of(hexa.sampling);
  const auto a = rtt::analysis::spot(cs, rtt::compile::PathId{0}, 1, 0, hexa);
  const auto b = rtt::analysis::spot(cs, rtt::compile::PathId{0}, 1, 0, list);
  REQUIRE(a.points.size() == b.points.size());
  for (std::size_t i = 0; i < a.points.size(); ++i) {
    INFO("point " << i);
    CHECK(same_bits(a.points[i].x, b.points[i].x));
    CHECK(same_bits(a.points[i].y, b.points[i].y));
    CHECK(same_bits(a.points[i].weight, b.points[i].weight));
  }
  CHECK(same_bits(a.stats.rms_centroid, b.stats.rms_centroid));
  // Content guard: a real spot with several points and a finite size.
  CHECK(a.points.size() > 10);
  CHECK(a.stats.rms_centroid > 0.0);
}
