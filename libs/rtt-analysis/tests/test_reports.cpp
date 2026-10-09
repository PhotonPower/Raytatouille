// Reports as data (#177): raytrace report against RayPaths (#80), an axial ray through the
// vertices, the system data report against the paraxial prescription, and the dimension
// report against closed forms (centre thickness, edge thickness of a plano-convex lens).

#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "rtt/analysis/reports.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/paraxial/prescription.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::analysis::ApertureKind;
using rtt::analysis::DimensionReport;
using rtt::analysis::RaytraceReport;
using rtt::analysis::RaytraceRow;
using rtt::analysis::SegmentDimensions;
using rtt::analysis::SystemReport;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::material::MaterialLibrary;
using rtt::model::Element;
using rtt::model::System;
using rtt::trace::RayBatch;
using rtt::trace::RayStatus;

namespace {

System load(const std::string& relative) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + relative);
}

/// tests/reference/m1/singlet_const.rtt.json: stop at z = 0, plano-convex lens CONST:1.5168 at
/// z = 5 (L1.S1 R = 51.68 mm, L1.S2 plane 4 mm behind), both apertures 12.7 mm, image at
/// z = 106.363.
CompiledSystem singlet() {
  return rtt::compile::compile(load("m1/singlet_const.rtt.json"), MaterialLibrary{});
}

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

std::uint32_t surface(const CompiledSystem& cs, const std::string& id) {
  for (std::uint32_t i = 0; i < cs.surfaces().size(); ++i) {
    if (cs.surfaces()[i].id.str() == id) return i;
  }
  throw std::logic_error("no surface " + id);
}

}  // namespace

TEST_CASE("raytrace report: global columns are bitwise those of RayPaths", "[reports]") {
  const CompiledSystem cs = singlet();
  const std::uint16_t fields[] = {0, 1, 2};
  const RayBatch start =
      rtt::trace::make_rays(cs, PathId{0}, fields, cs.reference_wavelength(),
                            rtt::trace::HexapolarPupil{3}, rtt::trace::Aiming::Real);
  // One ray steeply off axis, so that it is lost and its rows end early.
  RayBatch with_lost = start;
  with_lost.dir_x()[1] = 0.0;
  with_lost.dir_y()[1] = 0.8;
  with_lost.dir_z()[1] = 0.6;

  const RaytraceReport report = rtt::analysis::raytrace_report(cs, PathId{0}, with_lost);
  RayBatch traced = with_lost;
  rtt::trace::RayPaths paths;
  static_cast<void>(rtt::trace::SequentialTracer{}.trace(cs, PathId{0}, traced, paths));

  REQUIRE(report.rays == with_lost.size());
  REQUIRE(report.slots == paths.slots);
  std::size_t expected_rows = 0;
  for (std::size_t r = 0; r < paths.ray_count(); ++r) expected_rows += paths.count[r];
  REQUIRE(report.rows.size() == expected_rows);
  REQUIRE(paths.count[1] < paths.slots);  // the off-axis ray stopped

  std::size_t i = 0;
  for (std::size_t r = 0; r < paths.ray_count(); ++r) {
    for (std::size_t s = 0; s < paths.count[r]; ++s, ++i) {
      INFO("ray " << r << ", slot " << s);
      const RaytraceRow& row = report.rows[i];
      CHECK(row.ray == paths.ray_indices[r]);
      CHECK(row.slot == s);
      const auto p = paths.position_at(r, s);
      const auto d = paths.direction_at(r, s);
      CHECK((same(row.x, p.x()) && same(row.y, p.y()) && same(row.z, p.z())));
      CHECK((same(row.dx, d.x()) && same(row.dy, d.y()) && same(row.dz, d.z())));
      CHECK(same(row.opl, paths.opl[r * paths.slots + s]));
      CHECK(same(row.weight, paths.weight[r * paths.slots + s]));
      CHECK(row.status == paths.status[r * paths.slots + s]);
      if (s == 0) {
        CHECK(row.surface == rtt::trace::kNoSurface);
        CHECK(std::isnan(row.local_x));
      } else {
        CHECK(row.surface == paths.event_surfaces[s - 1]);
      }
    }
  }
  // The start batch is not changed.
  CHECK(same(with_lost.dir_y()[1], 0.8));
}

TEST_CASE("raytrace report: an axial ray meets every surface at its vertex", "[reports]") {
  const CompiledSystem cs = singlet();
  RayBatch start(1);  // on the axis along +z from the origin
  const RaytraceReport report = rtt::analysis::raytrace_report(cs, PathId{0}, start);
  REQUIRE(report.rows.size() == report.slots);  // the ray arrives
  for (std::size_t s = 1; s < report.slots; ++s) {
    const RaytraceRow& row = report.rows[s];
    INFO("surface " << cs.surfaces()[row.surface].id.str());
    const auto vertex = cs.surfaces()[row.surface].to_global.apply_point(rtt::math::Vec3::Zero());
    CHECK(std::abs(row.x - vertex.x()) <= 1e-12);
    CHECK(std::abs(row.y - vertex.y()) <= 1e-12);
    CHECK(std::abs(row.z - vertex.z()) <= 1e-12);
    CHECK(std::abs(row.local_x) <= 1e-12);
    CHECK(std::abs(row.local_y) <= 1e-12);
    CHECK(std::abs(row.local_z) <= 1e-12);
    CHECK(std::abs(row.local_dz - 1.0) <= 1e-12);
    CHECK(row.status == RayStatus::Alive);
  }
}

TEST_CASE("raytrace report: invalid input", "[reports]") {
  const CompiledSystem cs = singlet();
  CHECK_THROWS_AS(rtt::analysis::raytrace_report(cs, PathId{7}, RayBatch(1)),
                  std::invalid_argument);
  CHECK_THROWS_AS(rtt::analysis::raytrace_report(cs, PathId{0}, RayBatch(0)),
                  std::invalid_argument);
}

TEST_CASE("system report: settings and the paraxial prescription", "[reports]") {
  const CompiledSystem cs = singlet();
  const SystemReport r = rtt::analysis::system_report(cs, PathId{0}, 1);
  CHECK(r.wavelength == 1);
  CHECK(r.wavelengths_um == cs.wavelengths_um());
  CHECK(r.reference_wavelength == cs.reference_wavelength());
  CHECK(r.field_count == 3);
  CHECK(r.surface_count == cs.surfaces().size());
  CHECK(r.event_count == cs.paths()[0].events.size());
  REQUIRE(r.stop.has_value());
  CHECK(cs.surfaces()[*r.stop].id.str() == "STO");
  REQUIRE(r.prescription.has_value());
  CHECK(r.warnings.empty());
  // The same numbers as paraxial::prescription.
  const rtt::paraxial::Prescription p = rtt::paraxial::prescription(cs, PathId{0}, 1);
  CHECK(same(r.prescription->first_order.efl, p.first_order.efl));
  CHECK(same(r.prescription->total_track, p.total_track));
  REQUIRE(r.prescription->surfaces.size() == p.surfaces.size());
  CHECK_THROWS_AS(rtt::analysis::system_report(cs, PathId{0}, 9), std::invalid_argument);
  CHECK_THROWS_AS(rtt::analysis::system_report(cs, PathId{3}, 0), std::invalid_argument);
}

TEST_CASE("system report: a path without paraxial data keeps the settings", "[reports]") {
  // Diffraction order +1 at the grating: not rotationally symmetric for the paraxial trace.
  const CompiledSystem cs =
      rtt::compile::compile(load("m4/grating_transmission.rtt.json"), MaterialLibrary{});
  const PathId order = *cs.find_path("order +1");
  const SystemReport r = rtt::analysis::system_report(cs, order, 0);
  CHECK_FALSE(r.prescription.has_value());
  REQUIRE(r.warnings.size() == 1);
  CHECK(r.warnings[0].code == "report.paraxial_unavailable");
  CHECK(r.warnings[0].severity == rtt::model::Severity::Warning);
  CHECK(r.event_count == 3);
}

TEST_CASE("dimension report: plano-convex singlet against closed forms", "[reports]") {
  const CompiledSystem cs = singlet();
  const DimensionReport d = rtt::analysis::dimension_report(cs);
  REQUIRE(d.segments.size() == 1);  // only the lens has two surfaces
  const SegmentDimensions& s = d.segments[0];
  CHECK(s.first_surface == surface(cs, "L1.S1"));
  CHECK(cs.elements()[s.element].name == "L1");
  CHECK(s.coaxial);
  // Centre thickness: the vertex distance of the model, L1.S2 at z = 4 in the lens.
  CHECK(std::abs(s.centre_thickness - 4.0) <= 1e-12);
  CHECK(s.semi_diameter_first == 12.7);
  CHECK(s.semi_diameter_second == 12.7);
  CHECK(s.aperture_first == ApertureKind::Circular);
  CHECK(s.diameter == 25.4);
  // Edge thickness at h = 12.7: t + sag2(h) - sag1(h), sag of a sphere R - sqrt(R^2 - h^2)
  // (k = 0), sag2 = 0 for the plane.
  const double r = 51.68;
  const double h = 12.7;
  const double expected = 4.0 - (r - std::sqrt(r * r - h * h));
  CHECK(std::abs(s.edge_thickness - expected) <= 1e-12);
  CHECK(expected > 2.0);  // a real edge, about 2.4 mm
}

TEST_CASE("dimension report: larger semi-diameter, circumscribed apertures, tilt", "[reports]") {
  System sys = load("m1/singlet_const.rtt.json");
  Element& lens = std::get<Element>(sys.root.children[1].value);
  // L1.S1 circular 12.7, L1.S2 a rectangle 6 x 8 (half diagonal 5): the edge sits at 12.7.
  lens.surfaces[1].aperture = rtt::model::RectangularAperture{3.0, 4.0};
  CompiledSystem cs = rtt::compile::compile(sys, MaterialLibrary{});
  SegmentDimensions s = rtt::analysis::dimension_report(cs).segments.at(0);
  CHECK(s.semi_diameter_second == 5.0);
  CHECK(s.aperture_second == ApertureKind::Rectangular);
  CHECK(s.diameter == 25.4);
  // Ellipse: the major semi-axis; no aperture: NaN and no edge.
  lens.surfaces[1].aperture = rtt::model::EllipticalAperture{2.0, 14.0};
  cs = rtt::compile::compile(sys, MaterialLibrary{});
  s = rtt::analysis::dimension_report(cs).segments.at(0);
  CHECK(s.semi_diameter_second == 14.0);
  CHECK(s.diameter == 28.0);
  lens.surfaces[1].aperture.reset();
  cs = rtt::compile::compile(sys, MaterialLibrary{});
  s = rtt::analysis::dimension_report(cs).segments.at(0);
  CHECK(s.aperture_second == ApertureKind::None);
  CHECK(std::isnan(s.semi_diameter_second));
  CHECK(s.semi_diameter_first == 12.7);
  CHECK(std::isnan(s.edge_thickness));
  // A tilted second surface: not coaxial, the thicknesses are NaN.
  lens.surfaces[1].aperture = rtt::model::CircularAperture{12.7, 0.0};
  lens.surfaces[1].pose.rotation_deg[0] = rtt::model::Param(5.0);
  cs = rtt::compile::compile(sys, MaterialLibrary{});
  s = rtt::analysis::dimension_report(cs).segments.at(0);
  CHECK_FALSE(s.coaxial);
  CHECK(std::isnan(s.centre_thickness));
  CHECK(std::isnan(s.edge_thickness));
}
