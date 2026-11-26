// The optimization on folded paths (#357, ADR 0054, point 5, PR B1): first_order has data for a
// plane fold mirror, so the merit function runs on such a path as on its unfolded counterpart,
// without and with the overlay of a tolerance trial (the nominal basis skips the fold as a
// vertex for the BFL of a trial).

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>

#include "rtt/io/json_io.hpp"
#include "rtt/model/system.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/overlay.hpp"

namespace {

/// The reference file `name` with an EFL and a BFL operand and an RMS spot generator.
rtt::model::System with_merit(const std::string& name) {
  const std::ifstream in(std::string(RTT_REFERENCE_DIR) + "/" + name);
  REQUIRE(in);
  std::ostringstream buffer;
  buffer << in.rdbuf();
  std::string text = buffer.str();
  const std::size_t end = text.rfind('}');
  REQUIRE(end != std::string::npos);
  text.insert(end, R"(, "optimization": {
    "operands": [{"type": "efl", "path": "main", "target": 100.0, "weight": 1.0},
                 {"type": "bfl", "path": "main", "target": 90.0, "weight": 1.0},
                 {"type": "image_fnumber", "path": "main", "target": 5.0, "weight": 1.0},
                 {"type": "opd_rms", "path": "main", "field": 2, "target": 0.0, "weight": 1.0}],
    "generators": [{"type": "rms_spot", "path": "main"}]})");
  return rtt::io::parse_system(text);
}

rtt::optim::MeritEvaluation evaluate(const rtt::model::System& s, bool overlay) {
  const rtt::material::MaterialLibrary lib;
  std::optional<rtt::optim::TrialOverlay> o;
  if (overlay) o.emplace(rtt::optim::TrialOverlay{{}, rtt::optim::nominal_basis(s, lib, nullptr)});
  const rtt::optim::MeritFunction merit(s, lib, nullptr, o ? &*o : nullptr);
  return merit.evaluate(merit.start());
}

}  // namespace

TEST_CASE("operands on a folded path equal those of the unfolded one", "[optim][fold][357]") {
  const rtt::model::System folded = with_merit("m6/singlet_fold45.rtt.json");
  const rtt::model::System flat = with_merit("m1/singlet_const.rtt.json");
  for (const bool overlay : {false, true}) {
    INFO("overlay " << overlay);
    const auto f = evaluate(folded, overlay);
    const auto u = evaluate(flat, overlay);
    REQUIRE(f.values.size() == u.values.size());
    REQUIRE(f.residuals.size() == u.residuals.size());
    for (std::size_t i = 0; i < f.values.size(); ++i) {
      INFO("value " << i);
      CHECK(std::abs(f.values[i] - u.values[i]) <= 1e-9 * std::max(1.0, std::abs(u.values[i])));
    }
    // The spot residuals are the x and y of each ray in the local frame of the detector, which
    // the fold turns (the detector faces +y): the same magnitudes, signs may differ.
    for (std::size_t i = 0; i < f.residuals.size(); ++i) {
      INFO("residual " << i);
      CHECK(std::abs(std::abs(f.residuals[i]) - std::abs(u.residuals[i])) <= 1e-9);
    }
    // content guard: EFL near 100 mm, BFL 97.4 mm (the singlet of m1), F/5, aberrations
    CHECK(std::abs(u.values[0] - 100.0) < 2.0);
    CHECK(std::abs(u.values[1] - 97.4) < 0.5);
    CHECK(std::abs(u.values[2] - 5.0) < 0.2);
    CHECK(u.values[3] > 0.01);
  }
}

TEST_CASE("the BFL of a trial behind a fold before a concave mirror", "[optim][fold][357]") {
  // The fold at z = 20 into +y, then a concave mirror R = -100 mm at normal incidence, which is
  // the last vertex (it changes the sign of the paraxial index, the fold does not): with the
  // overlay of a trial the BFL is measured from the concave mirror, as in the unfolded system.
  const auto with_bfl = [](std::string text) {
    const std::size_t end = text.rfind('}');
    REQUIRE(end != std::string::npos);
    text.insert(end, R"(, "optimization": {"operands": [
      {"type": "efl", "path": "main", "target": 50.0, "weight": 1.0},
      {"type": "bfl", "path": "main", "target": 50.0, "weight": 1.0}]})");
    return rtt::io::parse_system(text);
  };
  const rtt::model::System folded = with_bfl(
      R"({"schema_version": "0.11.0", "name": "x", "units": {"length": "mm", "wavelength": "um"}, "wavelengths": [{"um": 0.5876, "reference": true}], "aperture": {"type": "epd", "value": 10.0}, "fields": {"points": [{}, {"y": 1.0}]}, "root": {"type": "assembly", "name": "system", "children": [{"type": "stop", "name": "stop", "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 5.0}}]}, {"type": "mirror", "name": "M1", "pose": {"position": [0.0, 0.0, 20.0], "rotation_deg": [45.0, 0.0, 0.0]}, "surfaces": [{"id": "M1", "aperture": {"type": "circular", "radius": 10.0}}]}, {"type": "mirror", "name": "M2", "pose": {"position": [0.0, 30.0, 20.0], "rotation_deg": [-90.0, 0.0, 0.0]}, "surfaces": [{"id": "M2", "shape": {"base": {"type": "conic", "radius": -100.0}}, "aperture": {"type": "circular", "radius": 8.0}}]}, {"type": "detector", "name": "image", "pose": {"position": [0.0, -18.0, 20.0], "rotation_deg": [90.0, 0.0, 0.0]}, "surfaces": [{"id": "IMG"}]}]}, "paths": [{"name": "main", "events": [{"surface": "STO", "kind": "transmit"}, {"surface": "M1", "kind": "reflect"}, {"surface": "M2", "kind": "reflect"}, {"surface": "IMG", "kind": "transmit"}]}]})");
  const rtt::model::System flat = with_bfl(
      R"({"schema_version": "0.11.0", "name": "x", "units": {"length": "mm", "wavelength": "um"}, "wavelengths": [{"um": 0.5876, "reference": true}], "aperture": {"type": "epd", "value": 10.0}, "fields": {"points": [{}, {"y": 1.0}]}, "root": {"type": "assembly", "name": "system", "children": [{"type": "stop", "name": "stop", "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 5.0}}]}, {"type": "mirror", "name": "M2", "pose": {"position": [0.0, 0.0, 50.0]}, "surfaces": [{"id": "M2", "shape": {"base": {"type": "conic", "radius": -100.0}}, "aperture": {"type": "circular", "radius": 8.0}}]}, {"type": "detector", "name": "image", "pose": {"position": [0.0, 0.0, 2.0], "rotation_deg": [180.0, 0.0, 0.0]}, "surfaces": [{"id": "IMG"}]}]}, "paths": [{"name": "main", "events": [{"surface": "STO", "kind": "transmit"}, {"surface": "M2", "kind": "reflect"}, {"surface": "IMG", "kind": "transmit"}]}]})");
  for (const bool overlay : {false, true}) {
    INFO("overlay " << overlay);
    const auto f = evaluate(folded, overlay);
    const auto u = evaluate(flat, overlay);
    REQUIRE(f.values.size() == 2);
    CHECK(std::abs(f.values[0] - u.values[0]) <= 1e-9);
    CHECK(std::abs(f.values[1] - u.values[1]) <= 1e-9);
    CHECK(std::abs(u.values[1] - 50.0) < 0.1);  // content guard: |R| / 2 from the mirror
  }
}
