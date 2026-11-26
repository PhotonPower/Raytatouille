// Paraxial aiming and spot diagrams of paths with plane fold mirrors (#357, ADR 0054, PR A).
// The aiming needs only object-space data (the entrance pupil), which lie before the first fold;
// the analyses that need image-space positions in global coordinates (the reference sphere, the
// paraxial image height) give an error for folded paths until PR B.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>

#include "rtt/analysis/chromatic.hpp"
#include "rtt/analysis/field.hpp"
#include "rtt/analysis/nominal_reference.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/system.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/sequential.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;

namespace {

rtt::model::System load(const std::string& name) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + name);
}

/// `text` with the first occurrence of `what` replaced by `with`.
std::string replaced(std::string text, const std::string& what, const std::string& with) {
  const std::size_t at = text.find(what);
  REQUIRE(at != std::string::npos);
  return text.replace(at, what.size(), with);
}

/// m6/fold45_focus with a stop of radius 5 mm 1 mm before the ideal lens (not in its plane: a
/// ray would meet the lens at distance 0): the file has none, and aiming needs one (EPD 10 mm as
/// in the file).
rtt::model::System fold45_with_stop() {
  const std::ifstream in(std::string(RTT_REFERENCE_DIR) + "/m6/fold45_focus.rtt.json");
  REQUIRE(in);
  std::ostringstream buffer;
  buffer << in.rdbuf();
  std::string text = replaced(buffer.str(), R"("children": [)", R"("children": [
      {"type": "stop", "name": "stop", "pose": {"position": [0.0, 0.0, -1.0]},
       "surfaces": [{"id": "STO",
       "aperture": {"type": "circular", "radius": 5.0}}]},)");
  text = replaced(text, R"("events": [)", R"("events": [{"surface": "STO", "kind": "transmit"},)");
  return rtt::io::parse_system(text);
}

CompiledSystem compiled(const rtt::model::System& s) {
  const rtt::material::MaterialLibrary lib;
  return rtt::compile::compile(s, lib);
}

rtt::analysis::SpotOptions options(rtt::trace::Aiming aiming) {
  rtt::analysis::SpotOptions o;
  o.aiming = aiming;
  return o;
}

}  // namespace

TEST_CASE("paraxial aiming through a 45 deg fold, against rays set by hand",
          "[analysis][fold][357]") {
  // m6/fold45_focus with a stop (fold45_with_stop()): an ideal lens f = 50 mm, EPD 10 mm, the
  // plane fold at z = 20, the detector 2 mm behind the focus. On axis the paraxially aimed rays
  // are the collimated rays one sets by hand (as the M6 tests do): parallel to z through the
  // pupil, so both land on the same points of the detector.
  const CompiledSystem cs = compiled(fold45_with_stop());
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, 0);
  REQUIRE(fo.folds == 1);
  REQUIRE(fo.entrance_pupil);
  REQUIRE(fo.entrance_pupil->diameter == 10.0);
  const std::uint16_t field = 0;
  rtt::trace::RayBatch aimed =
      rtt::trace::make_rays(cs, PathId{0}, std::span<const std::uint16_t>(&field, 1), 0,
                            rtt::trace::HexapolarPupil{6}, rtt::trace::Aiming::Paraxial);
  REQUIRE(aimed.size() == 127);
  rtt::trace::RayBatch hand(aimed.size());
  double r_max = 0.0;
  for (std::size_t i = 0; i < aimed.size(); ++i) {
    REQUIRE(aimed.dir_x()[i] == 0.0);
    REQUIRE(aimed.dir_y()[i] == 0.0);
    hand.pos_x()[i] = aimed.pos_x()[i];
    hand.pos_y()[i] = aimed.pos_y()[i];
    hand.pos_z()[i] = -2.0;
    r_max = std::max(r_max, std::hypot(aimed.pos_x()[i], aimed.pos_y()[i]));
  }
  CHECK(std::abs(r_max - 5.0) <= 1e-12);  // the rim of the entrance pupil
  const rtt::trace::SequentialTracer tracer;
  (void)tracer.trace(cs, PathId{0}, aimed);
  (void)tracer.trace(cs, PathId{0}, hand);
  double spread = 0.0;
  for (std::size_t i = 0; i < aimed.size(); ++i) {
    INFO("ray " << i);
    REQUIRE(aimed.status()[i] == rtt::trace::RayStatus::Alive);
    REQUIRE(hand.status()[i] == rtt::trace::RayStatus::Alive);
    CHECK(std::abs(aimed.pos_x()[i] - hand.pos_x()[i]) <= 1e-12);
    CHECK(std::abs(aimed.pos_y()[i] - hand.pos_y()[i]) <= 1e-12);
    CHECK(std::abs(aimed.pos_z()[i] - hand.pos_z()[i]) <= 1e-12);
    spread = std::max(spread, std::abs(aimed.pos_z()[i] - 20.0));
  }
  // content guard: the detector lies along y behind the fold; the defocus spreads the rays over
  // about 2 / 50 of the pupil in z there
  CHECK(spread > 0.1);
  const auto s = rtt::analysis::spot(cs, PathId{0}, 0, 0, options(rtt::trace::Aiming::Paraxial));
  REQUIRE(s.rays_arrived == 127);
}

TEST_CASE("the spot behind a fold equals that of the unfolded singlet", "[analysis][fold][357]") {
  // m6/singlet_fold45 is m1/singlet_const with a plane 45 deg mirror 40 mm behind the lens
  // (decision 3 on #357: spot RMS to 1e-12 mm).
  const CompiledSystem folded = compiled(load("m6/singlet_fold45.rtt.json"));
  const CompiledSystem flat = compiled(load("m1/singlet_const.rtt.json"));
  for (const auto aiming : {rtt::trace::Aiming::Paraxial, rtt::trace::Aiming::Real}) {
    for (std::uint16_t field = 0; field < 3; ++field) {
      INFO("field " << field << ", aiming " << static_cast<int>(aiming));
      const auto a = rtt::analysis::spot(folded, PathId{0}, field, 1, options(aiming));
      const auto b = rtt::analysis::spot(flat, PathId{0}, field, 1, options(aiming));
      REQUIRE(a.rays_arrived == b.rays_arrived);
      CHECK(std::abs(a.stats.rms_centroid - b.stats.rms_centroid) <= 1e-12);
      CHECK(std::abs(a.stats.geo_centroid - b.stats.geo_centroid) <= 1e-12);
      CHECK(b.stats.rms_centroid > 1e-3);  // content guard: the singlet has aberrations
    }
  }
}

TEST_CASE("image-space analyses behind a fold equal those of the unfolded singlet",
          "[analysis][fold][357]") {
  // PR B1 (ADR 0054, point 5): the reference sphere about the exit pupil mapped by image_space,
  // fields as paraxial image height in the unfolded image plane, and longitudinal colour,
  // distortion and field curvature measured along the unfolded axis. m6/singlet_fold45 against
  // m1/singlet_const (decision 3 on #357: 1e-12).
  const rtt::model::System s = load("m6/singlet_fold45.rtt.json");
  const CompiledSystem cs = compiled(s);
  const CompiledSystem flat = compiled(load("m1/singlet_const.rtt.json"));
  const auto near = [](double a, double b) {
    return std::abs(a - b) <= 1e-12 * std::max(1.0, std::abs(b));
  };
  // OPD on the reference sphere, field 2
  const auto of = rtt::analysis::opd_map(cs, PathId{0}, 2, 1);
  const auto ou = rtt::analysis::opd_map(flat, PathId{0}, 2, 1);
  REQUIRE(of.arrived == ou.arrived);
  CHECK(std::abs(of.rms - ou.rms) <= 1e-9);  // waves
  CHECK(std::abs(of.pv - ou.pv) <= 1e-9);
  CHECK(near(of.sphere.radius, ou.sphere.radius));
  CHECK(ou.rms > 0.01);  // content guard: the singlet has aberrations
  // the exit pupil of a nominal reference: the same point as make_reference takes
  const auto e = rtt::analysis::nominal_reference(cs, PathId{0}).exit_pupil();
  REQUIRE(e);
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, cs.reference_wavelength());
  // geometrically: the exit pupil of the unfolded singlet at z_XP, mirrored at the 45 deg
  // mirror at z = 49 into -y (the mirror normal is (0, -sin 45, cos 45))
  const double z_xp =
      *rtt::paraxial::first_order(flat, PathId{0}, flat.reference_wavelength()).exit_pupil->z;
  CHECK((*e - rtt::math::Vec3(0.0, -(49.0 - z_xp), 49.0)).norm() <= 1e-12);
  CHECK(fo.folds == 1);
  // longitudinal colour
  const auto lf = rtt::analysis::longitudinal_colour(cs, PathId{0});
  const auto lu = rtt::analysis::longitudinal_colour(flat, PathId{0});
  REQUIRE(lf.foci.size() == lu.foci.size());
  for (std::size_t k = 0; k < lf.foci.size(); ++k) {
    CHECK(near(lf.foci[k].paraxial_z, lu.foci[k].paraxial_z));
    CHECK(near(lf.foci[k].real_z, lu.foci[k].real_z));
  }
  CHECK(std::abs(lu.foci[1].real_z - lu.foci[1].paraxial_z) > 1e-3);  // spherical aberration
  // distortion and field curvature
  const auto df = rtt::analysis::distortion(cs, PathId{0}, 1);
  const auto du = rtt::analysis::distortion(flat, PathId{0}, 1);
  REQUIRE(df.size() == du.size());
  for (std::size_t k = 0; k < df.size(); ++k) {
    INFO("distortion point " << k);
    CHECK(near(df[k].real_height, du[k].real_height));
    CHECK(near(df[k].paraxial_height, du[k].paraxial_height));
    CHECK(std::abs(df[k].percent - du[k].percent) <= 1e-10);
  }
  CHECK(std::abs(du.back().percent) > 1e-3);  // content guard
  const auto cf = rtt::analysis::field_curvature(cs, PathId{0}, 1);
  const auto cu = rtt::analysis::field_curvature(flat, PathId{0}, 1);
  REQUIRE(cf.size() == cu.size());
  for (std::size_t k = 0; k < cf.size(); ++k) {
    INFO("field curvature point " << k);
    // field curvature intersects two neighbouring rays (ADR 0032): the rounding of their
    // positions (about 1e-14 mm, a few ulp more behind the fold) grows by the inverse of their
    // small angle; measured up to 1.2e-8 mm at the field edge
    CHECK(std::abs(cf[k].tangential - cu[k].tangential) <= 1e-7);
    CHECK(std::abs(cf[k].sagittal - cu[k].sagittal) <= 1e-7);
  }
  CHECK(std::abs(cu.back().astigmatism) > 1e-3);  // content guard
  // fields as paraxial image height (in the unfolded image plane)
  rtt::model::System h = s;
  h.fields.type = rtt::model::FieldType::ParaxialImageHeight;
  rtt::model::System hu = load("m1/singlet_const.rtt.json");
  hu.fields.type = rtt::model::FieldType::ParaxialImageHeight;
  const auto sh =
      rtt::analysis::spot(compiled(h), PathId{0}, 1, 1, options(rtt::trace::Aiming::Paraxial));
  const auto su =
      rtt::analysis::spot(compiled(hu), PathId{0}, 1, 1, options(rtt::trace::Aiming::Paraxial));
  REQUIRE(sh.rays_arrived == su.rays_arrived);
  CHECK(std::abs(sh.stats.rms_centroid - su.stats.rms_centroid) <= 1e-12);
  CHECK(su.stats.rms_centroid > 1e-3);
}

TEST_CASE("the fold mirror as the first surface bounds the start plane", "[analysis][fold][357]") {
  // The mirror at z = 5, tilted by 45 deg, is the first surface; the stop lies 15 mm behind it
  // along +y (not 10 mm: the top ray would meet the mirror in the stop plane, at distance 0), so
  // the entrance pupil is at the unfolded z = 20. The lower half of the mirror
  // reaches z = -5 for the bundle of radius 10 mm; the start plane must lie before it, else the
  // rays start behind the mirror and miss it.
  const CompiledSystem cs = compiled(rtt::io::parse_system(
      R"({"schema_version": "0.11.0", "name": "fold before the stop", "units": {"length": "mm", "wavelength": "um"}, "wavelengths": [{"um": 0.5876, "reference": true}], "aperture": {"type": "epd", "value": 20.0}, "fields": {"points": [{}]}, "root": {"type": "assembly", "name": "system", "children": [{"type": "mirror", "name": "M", "pose": {"position": [0.0, 0.0, 5.0], "rotation_deg": [45.0, 0.0, 0.0]}, "surfaces": [{"id": "M", "aperture": {"type": "circular", "radius": 20.0}}]}, {"type": "stop", "name": "stop", "pose": {"position": [0.0, 15.0, 5.0], "rotation_deg": [-90.0, 0.0, 0.0]}, "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 10.0}}]}, {"type": "detector", "name": "image", "pose": {"position": [0.0, 50.0, 5.0], "rotation_deg": [-90.0, 0.0, 0.0]}, "surfaces": [{"id": "IMG"}]}]}, "paths": [{"name": "main", "events": [{"surface": "M", "kind": "reflect"}, {"surface": "STO", "kind": "transmit"}, {"surface": "IMG", "kind": "transmit"}]}]})"));
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, 0);
  REQUIRE(fo.folds == 1);
  REQUIRE(fo.entrance_pupil);
  REQUIRE(fo.entrance_pupil->z);
  CHECK(std::abs(*fo.entrance_pupil->z - 20.0) <= 1e-12);
  for (const auto aiming : {rtt::trace::Aiming::Paraxial, rtt::trace::Aiming::Real}) {
    INFO("aiming " << static_cast<int>(aiming));
    const auto s = rtt::analysis::spot(cs, PathId{0}, 0, 0, options(aiming));
    CHECK(s.rays_launched == 127);
    CHECK(s.rays_arrived == 127);
  }
}

TEST_CASE("a periscope over x and y with oblique fields", "[analysis][fold][357]") {
  // m1/singlet_const with a mirror at z = 49 into +y and a second one 20 mm further into +x
  // (two folds, even parity, a linear map that is not symmetric, so linear and its transpose
  // differ), fields (0, 0), (3, 4), (5, 0), (0, 5) deg: distortion, field curvature and OPD equal
  // those of the unfolded singlet with the same fields.
  const CompiledSystem cs =
      compiled(rtt::io::parse_system(
          R"({"schema_version": "0.11.0", "name": "periscope x-y", "units": {"length": "mm", "wavelength": "um"}, "wavelengths": [{"um": 0.4861}, {"um": 0.5876, "reference": true}, {"um": 0.6563}], "aperture": {"type": "epd", "value": 20.0}, "fields": {"points": [{}, {"x": 3.0, "y": 4.0}, {"x": 5.0}, {"y": 5.0}]}, "root": {"type": "assembly", "name": "system", "children": [{"type": "stop", "name": "stop", "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 10.0}}]}, {"type": "lens", "name": "L1", "pose": {"position": [0.0, 0.0, 5.0]}, "material": "CONST:1.5168", "surfaces": [{"id": "L1.S1", "shape": {"base": {"type": "conic", "radius": 51.68}}, "aperture": {"type": "circular", "radius": 12.7}}, {"id": "L1.S2", "pose": {"position": [0.0, 0.0, 4.0]}, "aperture": {"type": "circular", "radius": 12.7}}]}, {"type": "mirror", "name": "M1", "pose": {"position": [0.0, 0.0, 49.0], "rotation_deg": [45.0, 0.0, 0.0]}, "surfaces": [{"id": "M1"}]}, {"type": "mirror", "name": "M2", "pose": {"position": [0.0, 20.0, 49.0], "rotation_deg": [90.0, 45.0, 0.0]}, "surfaces": [{"id": "M2"}]}, {"type": "detector", "name": "image", "pose": {"position": [37.363, 20.0, 49.0], "rotation_deg": [0.0, 90.0, 0.0]}, "surfaces": [{"id": "IMG"}]}]}, "paths": [{"name": "main", "events": [{"surface": "STO", "kind": "transmit"}, {"surface": "L1.S1", "kind": "refract"}, {"surface": "L1.S2", "kind": "refract"}, {"surface": "M1", "kind": "reflect"}, {"surface": "M2", "kind": "reflect"}, {"surface": "IMG", "kind": "transmit"}]}]})"));
  const CompiledSystem flat = compiled(rtt::io::parse_system(
      R"({"schema_version": "0.11.0", "name": "Plano-convex singlet f = 100 mm, constant index 1.5168", "units": {"length": "mm", "wavelength": "um"}, "wavelengths": [{"um": 0.4861}, {"um": 0.5876, "reference": true}, {"um": 0.6563}], "aperture": {"type": "epd", "value": 20.0}, "fields": {"points": [{}, {"x": 3.0, "y": 4.0}, {"x": 5.0}, {"y": 5.0}]}, "root": {"type": "assembly", "name": "system", "children": [{"type": "stop", "name": "stop", "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 10.0}}]}, {"type": "lens", "name": "L1", "pose": {"position": [0.0, 0.0, 5.0]}, "material": "CONST:1.5168", "surfaces": [{"id": "L1.S1", "shape": {"base": {"type": "conic", "radius": 51.68}}, "aperture": {"type": "circular", "radius": 12.7}}, {"id": "L1.S2", "pose": {"position": [0.0, 0.0, 4.0]}, "aperture": {"type": "circular", "radius": 12.7}}]}, {"type": "detector", "name": "image", "pose": {"position": [0.0, 0.0, 106.363]}, "surfaces": [{"id": "IMG"}]}]}, "paths": [{"name": "main", "events": "auto"}]})"));
  const auto fo = rtt::paraxial::first_order(cs, PathId{0}, 1);
  REQUIRE(fo.folds == 2);
  CHECK((fo.image_space.linear - fo.image_space.linear.transpose()).cwiseAbs().maxCoeff() > 0.5);
  // UnfoldedFrame: there and back, and the image vertex at the unfolded (0, 0, 106.363)
  const rtt::math::Vec3 p(1.0, -2.0, 70.0);
  CHECK((fo.image_space.to_unfolded(fo.image_space.to_global(p)) - p).norm() <= 1e-12);
  CHECK((fo.image_space.to_unfolded(cs.surfaces().back().to_global.translation()) -
         rtt::math::Vec3(0.0, 0.0, 106.363))
            .norm() <= 1e-12);
  for (std::uint16_t f = 1; f < 4; ++f) {
    INFO("field " << f);
    const auto& field = cs.fields().points[f];
    const auto df = rtt::analysis::distortion_at(cs, PathId{0}, field, 1);
    const auto du = rtt::analysis::distortion_at(flat, PathId{0}, field, 1);
    CHECK(std::abs(df.real_height - du.real_height) <= 1e-12);
    CHECK(std::abs(df.paraxial_height - du.paraxial_height) <= 1e-12);
    const auto cf = rtt::analysis::field_curvature_at(cs, PathId{0}, field, 1);
    const auto cu = rtt::analysis::field_curvature_at(flat, PathId{0}, field, 1);
    CHECK(std::abs(cf.tangential - cu.tangential) <= 1e-7);
    CHECK(std::abs(cf.sagittal - cu.sagittal) <= 1e-7);
    const auto of = rtt::analysis::opd_map(cs, PathId{0}, f, 1);
    const auto ou = rtt::analysis::opd_map(flat, PathId{0}, f, 1);
    CHECK(std::abs(of.rms - ou.rms) <= 1e-9);
    CHECK(std::abs(du.real_height) > 1.0);  // content guard: off axis
  }
}
