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

TEST_CASE("image-space analyses of folded paths wait for PR B", "[analysis][fold][357]") {
  const rtt::model::System s = load("m6/singlet_fold45.rtt.json");
  const CompiledSystem cs = compiled(s);
  // the reference sphere about the exit pupil (OPD, PSF, MTF)
  try {
    (void)rtt::analysis::opd_map(cs, PathId{0}, 0, 1);
    FAIL("no exception");
  } catch (const std::invalid_argument& e) {
    CHECK(std::string(e.what()).find("folded path") != std::string::npos);
  }
  // fields as paraxial image height
  rtt::model::System h = s;
  h.fields.type = rtt::model::FieldType::ParaxialImageHeight;
  const CompiledSystem ch = compiled(h);
  try {
    (void)rtt::analysis::spot(ch, PathId{0}, 1, 1, options(rtt::trace::Aiming::Paraxial));
    FAIL("no exception");
  } catch (const std::invalid_argument& e) {
    CHECK(std::string(e.what()).find("folded path") != std::string::npos);
  }
}
