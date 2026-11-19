// Grid sag in the analyses (#316, ADR 0051, point 9, PR 4): a grid has no paraxial model (order
// 1), so aiming and OPD run with the reference of a nominal system (ADR 0041, point 3), and
// without it they throw ParaxialError.
//
// Independent check: with exact derivatives the bicubic Hermite interpolant reproduces every
// polynomial of degree <= 3 in x and in y exactly (the cubic Hermite polynomial is exact for
// cubics, Burden/Faires section 3.4, and the tensor product keeps that per axis). A grid of the
// astigmatic sag a (x^2 - y^2) / R^2 is therefore the same surface as the xy_polynomial with
// c_x2 = a, c_y2 = -a; the OPD of both mirrors must agree up to rounding and the Newton stopping
// criterion: |F| <= 1e-12 mm per intersection on each side moves the OPL by at most 2e-12 mm per
// surface, so |dW| <= 2 * 2e-12 * 1e3 / lambda = 7e-9 waves; the test uses 1e-8.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "rtt/analysis/nominal_reference.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;

namespace {

constexpr double kLambda = 0.5876;  // um
constexpr double kR = 10.0;         // normalisation radius of the XY polynomial, mm
constexpr double kA = 1e-4;         // mm

CompiledSystem compiled(const std::string& json) {
  return rtt::compile::compile(rtt::io::parse_system(json), rtt::material::MaterialLibrary{});
}

/// Collimated light through a stop (radius 10) onto a plane mirror with `terms`, back through an
/// ideal lens f = 5000 at z = -20 to the focus at z = -5020.
std::string mirror_system(const std::string& terms) {
  return R"({
  "schema_version": "0.11.0",
  "units": {"length": "mm", "wavelength": "um"},
  "environment": {"medium": "VACUUM"},
  "wavelengths": [{"um": 0.5876, "reference": true}],
  "aperture": {"type": "epd", "value": 20.0},
  "fields": {"points": [{}]},
  "root": {"type": "assembly", "name": "system", "children": [
    {"type": "stop", "name": "stop", "pose": {"position": [0.0, 0.0, -1.0]},
     "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 10.0}}]},
    {"type": "mirror", "name": "M", "surfaces": [
      {"id": "M", "shape": {"base": {"type": "plane"}, "terms": [)" +
         terms + R"(]},
       "aperture": {"type": "circular", "radius": 12.0}}]},
    {"type": "thin_element", "name": "IL", "pose": {"position": [0.0, 0.0, -20.0]},
     "surfaces": [{"id": "IL", "aperture": {"type": "circular", "radius": 15.0},
                   "interaction": {"type": "ideal_lens", "focal_length": 5000.0}}]},
    {"type": "detector", "name": "F", "pose": {"position": [0.0, 0.0, -5020.0]},
     "surfaces": [{"id": "F"}]}
  ]},
  "paths": [{"name": "main", "events": "auto"}]
})";
}

/// 51 x 51 points with spacing 0.5 on [-12.5, 12.5]^2 of a (x^2 - y^2) / R^2 with its exact
/// derivatives.
std::string grid_term() {
  std::ostringstream z;
  std::ostringstream zx;
  std::ostringstream zy;
  std::ostringstream zxy;
  for (std::ostringstream* s : {&z, &zx, &zy, &zxy}) s->precision(17);
  constexpr int kN = 51;
  for (int j = 0; j < kN; ++j) {
    for (int i = 0; i < kN; ++i) {
      const double x = (i - 25) * 0.5;
      const double y = (j - 25) * 0.5;
      const char* sep = i + j > 0 ? ", " : "";
      z << sep << kA * (x * x - y * y) / (kR * kR);
      zx << sep << 2.0 * kA * x / (kR * kR);
      zy << sep << -2.0 * kA * y / (kR * kR);
      zxy << sep << 0.0;
    }
  }
  return R"({"type": "grid_sag", "nx": 51, "ny": 51, "spacing_x": 0.5, "spacing_y": 0.5, "sag": [)" +
         z.str() + R"(], "dz_dx": [)" + zx.str() + R"(], "dz_dy": [)" + zy.str() +
         R"(], "d2z_dxdy": [)" + zxy.str() + "]}";
}

std::string xy_term() {
  std::ostringstream s;
  s.precision(17);
  s << R"({"type": "xy_polynomial", "normalization_radius": )" << kR
    << R"(, "coefficients": [0.0, 0.0, )" << kA << ", 0.0, " << -kA << "]}";
  return s.str();
}

}  // namespace

TEST_CASE("a grid sag traces and gives the OPD with the nominal reference (ADR 0051, point 9)",
          "[grid][opd]") {
  const CompiledSystem grid = compiled(mirror_system(grid_term()));
  const CompiledSystem xy = compiled(mirror_system(xy_term()));
  const CompiledSystem nominal = compiled(mirror_system(R"({"type": "xy_polynomial",
      "normalization_radius": 10.0, "coefficients": [0.0]})"));
  // No paraxial model (order 1): first order and an analysis without reference throw.
  CHECK_THROWS_AS((void)rtt::paraxial::first_order(grid, PathId{0}, 0),
                  rtt::paraxial::ParaxialError);
  CHECK_THROWS_AS(rtt::analysis::opd_points(grid, PathId{0}, 0, 0, rtt::trace::GridPupil{16}),
                  rtt::paraxial::ParaxialError);
  rtt::analysis::OpdOptions o;
  o.reference = rtt::analysis::nominal_reference(nominal, PathId{0});
  const auto wg = rtt::analysis::opd_points(grid, PathId{0}, 0, 0, rtt::trace::GridPupil{24}, o);
  const auto wx = rtt::analysis::opd_points(xy, PathId{0}, 0, 0, rtt::trace::GridPupil{24}, o);
  REQUIRE(wg.points.size() == wx.points.size());
  std::size_t alive = 0;
  double spread = 0.0;  // content guard: the astigmatism is in the OPD
  for (std::size_t k = 0; k < wg.points.size(); ++k) {
    const auto& g = wg.points[k];
    const auto& x = wx.points[k];
    REQUIRE(g.status == x.status);
    if (g.status != rtt::trace::RayStatus::Alive) continue;
    ++alive;
    CHECK(std::abs(g.w - x.w) <= 1e-8);
    spread = std::max(spread, std::abs(g.w));
  }
  CHECK(alive > 300);
  // dW = -2 a (x^2 - y^2) / R^2: up to 2e3 kA / kLambda = 0.34 waves at the rim (minus piston).
  CHECK(spread > 0.25 * 2e3 * kA / kLambda);
}
