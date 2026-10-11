// Explicit pupil point lists (ADR 0032, point 2; Refs #202): pupil_points gives the list as is,
// make_rays aims every point exactly like a SinglePupilPoint, and an empty list or a non-finite
// point is an input error.

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sources.hpp"

using rtt::trace::PointListPupil;
using rtt::trace::pupil_points;
using rtt::trace::PupilPoint;
using rtt::trace::RayBatch;

namespace {

bool same_bits(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
}

}  // namespace

TEST_CASE("point list pupil: pupil_points returns the list unchanged, in its order",
          "[sources][point_list]") {
  // Includes a point outside the unit circle, a duplicate and -0.0, which must keep its sign bit.
  const PointListPupil list{{{0.25, -0.5}, {1.5, 0.0}, {0.25, -0.5}, {-0.0, 0.75}}};
  const std::vector<PupilPoint> p = pupil_points(list);
  REQUIRE(p.size() == list.points.size());
  for (std::size_t i = 0; i < p.size(); ++i) {
    INFO("point " << i);
    CHECK(same_bits(p[i].px, list.points[i].px));
    CHECK(same_bits(p[i].py, list.points[i].py));
  }
}

TEST_CASE("point list pupil: an empty list or a non-finite point is invalid_argument",
          "[sources][point_list]") {
  using Catch::Matchers::ContainsSubstring;
  CHECK_THROWS_AS(pupil_points(PointListPupil{}), std::invalid_argument);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  CHECK_THROWS_WITH(pupil_points(PointListPupil{{{0.0, 0.0}, {nan, 0.0}}}),
                    ContainsSubstring("point 1"));
  CHECK_THROWS_WITH(pupil_points(PointListPupil{{{0.0, 0.0}, {0.1, 0.2}, {0.0, -inf}}}),
                    ContainsSubstring("point 2"));
}

TEST_CASE("point list pupil: make_rays aims each point exactly like a SinglePupilPoint",
          "[sources][point_list]") {
  // Two fields: the order is field-major, then the list order (as for every sampling).
  const rtt::model::System s =
      rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/m1/singlet_const.rtt.json");
  const rtt::material::MaterialLibrary lib;
  const rtt::compile::CompiledSystem cs = rtt::compile::compile(s, lib);
  REQUIRE(cs.fields().points.size() >= 2);
  const std::array<std::uint16_t, 2> fields{0, 1};
  const PointListPupil list{{{0.0, 0.0}, {0.3, -0.7}, {-0.6, 0.2}, {0.0, 1.0}}};
  const RayBatch rays = rtt::trace::make_rays(cs, rtt::compile::PathId{0}, fields, 0, list);
  REQUIRE(rays.size() == fields.size() * list.points.size());
  // Content guard: the rays are aimed (Alive) and start at different points across the pupil (the
  // object is at infinity, so the directions of one field are equal).
  bool starts_differ = false;
  for (std::size_t f = 0; f < fields.size(); ++f) {
    for (std::size_t k = 0; k < list.points.size(); ++k) {
      const std::size_t i = f * list.points.size() + k;
      INFO("field " << fields[f] << ", point " << k);
      const std::array<std::uint16_t, 1> one{fields[f]};
      const RayBatch single =
          rtt::trace::make_rays(cs, rtt::compile::PathId{0}, one, 0,
                                rtt::trace::SinglePupilPoint{list.points[k].px, list.points[k].py});
      CHECK(rays.status()[i] == rtt::trace::RayStatus::Alive);
      CHECK(same_bits(rays.pos_x()[i], single.pos_x()[0]));
      CHECK(same_bits(rays.pos_y()[i], single.pos_y()[0]));
      CHECK(same_bits(rays.pos_z()[i], single.pos_z()[0]));
      CHECK(same_bits(rays.dir_x()[i], single.dir_x()[0]));
      CHECK(same_bits(rays.dir_y()[i], single.dir_y()[0]));
      CHECK(same_bits(rays.dir_z()[i], single.dir_z()[0]));
      CHECK(same_bits(rays.pupil_x()[i], list.points[k].px));
      CHECK(same_bits(rays.pupil_y()[i], list.points[k].py));
      CHECK(rays.field()[i] == fields[f]);
      CHECK(rays.status()[i] == single.status()[0]);
      if (k > 0 && rays.pos_y()[i] != rays.pos_y()[i - k]) starts_differ = true;
    }
  }
  CHECK(starts_differ);
}
