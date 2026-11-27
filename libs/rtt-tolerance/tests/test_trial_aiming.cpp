// One aiming for the criteria and the compensators of a trial (#291, decision 3; ADR 0030,
// addendum #291): TrialOptions::aiming, real by default. m1/two_lenses_stop_between (the stop
// between the lenses, so the aimings differ) with a spot criterion at the field edge.

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <span>

#include "rtt/analysis/spot.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/tolerance/trial.hpp"
#include "rtt/trace/sources.hpp"

TEST_CASE("the criteria of a trial aim as TrialOptions::aiming", "[tolerance][291]") {
  const rtt::model::System s = rtt::io::parse_system(
      R"({"schema_version": "0.12.0", "name": "Two singlets, stop between, stop-size aperture", "units": {"length": "mm", "wavelength": "um"}, "wavelengths": [{"um": 0.5876, "reference": true}], "aperture": {"type": "stop_size"}, "fields": {"points": [{}, {"y": 5.0}, {"x": -3.0, "y": 4.0}]}, "root": {"type": "assembly", "name": "system", "children": [{"type": "lens", "name": "L1", "pose": {"position": [0.0, 0.0, 10.0]}, "material": "CONST:1.5168", "surfaces": [{"id": "L1.S1", "shape": {"base": {"type": "conic", "radius": 60.0}}}, {"id": "L1.S2", "pose": {"position": [0.0, 0.0, 5.0]}, "shape": {"base": {"type": "conic", "radius": -60.0}}}]}, {"type": "stop", "name": "stop", "pose": {"position": [0.0, 0.0, 20.0]}, "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 4.0}}]}, {"type": "lens", "name": "L2", "pose": {"position": [0.0, 0.0, 25.0]}, "material": "CONST:1.5168", "surfaces": [{"id": "L2.S1", "shape": {"base": {"type": "conic", "radius": 45.0}}}, {"id": "L2.S2", "pose": {"position": [0.0, 0.0, 4.0]}}]}, {"type": "detector", "name": "image", "pose": {"position": [0.0, 0.0, 80.0]}, "surfaces": [{"id": "IMG"}]}]}, "paths": [{"name": "main", "events": "auto"}], "tolerances": {"items": [{"id": "l2-dz", "node": "L2", "kind": "shift_z", "limit": 0.05, "distribution": "uniform"}], "criteria": [{"id": "spot-edge", "operand": {"type": "spot_rms", "path": "main", "field": 2}, "max": 1.0}]}})");
  const rtt::material::MaterialLibrary lib;
  const rtt::compile::CompiledSystem cs = rtt::compile::compile(s, lib);
  double values[2] = {0.0, 0.0};
  int k = 0;
  for (const auto aiming : {rtt::trace::Aiming::Real, rtt::trace::Aiming::Paraxial}) {
    INFO("aiming " << static_cast<int>(aiming));
    rtt::tolerance::TrialOptions options;
    options.aiming = aiming;
    const rtt::tolerance::TrialContext context =
        rtt::tolerance::prepare(s, lib, nullptr, s.tolerances.items, options);
    const auto& nominal = context.nominal();
    REQUIRE(nominal.criteria.size() == 1);
    REQUIRE(nominal.criteria[0].value);
    // Phi_0 is the spot of the nominal system with that aiming (hexapolar, 6 rings)
    rtt::analysis::SpotOptions o;
    o.aiming = aiming;
    const auto spot =
        rtt::analysis::spot(cs, rtt::compile::PathId{0}, 2, cs.reference_wavelength(), o);
    CHECK(*nominal.criteria[0].value == spot.stats.rms_centroid);
    values[k++] = *nominal.criteria[0].value;
  }
  CHECK(std::abs(values[0] - values[1]) > 1e-9);  // content guard: the aimings differ
}

TEST_CASE("the compensator run of a trial aims as TrialOptions::aiming", "[tolerance][291]") {
  // the image distance compensates a shift of L2 on the spot at the field edge; the run with
  // paraxial aiming ends elsewhere than the one with real aiming, and TrialOptions::aiming
  // replaces compensators.aiming (decision 3 on #291)
  const rtt::model::System s = rtt::io::parse_system(
      R"({"schema_version": "0.12.0", "name": "Two singlets, stop between, stop-size aperture", "units": {"length": "mm", "wavelength": "um"}, "wavelengths": [{"um": 0.5876, "reference": true}], "aperture": {"type": "stop_size"}, "fields": {"points": [{}, {"y": 5.0}, {"x": -3.0, "y": 4.0}]}, "root": {"type": "assembly", "name": "system", "children": [{"type": "lens", "name": "L1", "pose": {"position": [0.0, 0.0, 10.0]}, "material": "CONST:1.5168", "surfaces": [{"id": "L1.S1", "shape": {"base": {"type": "conic", "radius": 60.0}}}, {"id": "L1.S2", "pose": {"position": [0.0, 0.0, 5.0]}, "shape": {"base": {"type": "conic", "radius": -60.0}}}]}, {"type": "stop", "name": "stop", "pose": {"position": [0.0, 0.0, 20.0]}, "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 4.0}}]}, {"type": "lens", "name": "L2", "pose": {"position": [0.0, 0.0, 25.0]}, "material": "CONST:1.5168", "surfaces": [{"id": "L2.S1", "shape": {"base": {"type": "conic", "radius": 45.0}}}, {"id": "L2.S2", "pose": {"position": [0.0, 0.0, 4.0]}}]}, {"type": "detector", "name": "image", "pose": {"position": [0.0, 0.0, 80.0]}, "surfaces": [{"id": "IMG"}]}]}, "paths": [{"name": "main", "events": "auto"}], "tolerances": {"items": [{"id": "l2-dz", "node": "L2", "kind": "shift_z", "limit": 0.5, "distribution": "uniform"}], "compensators": [{"node": "image", "param": "/pose/position/2", "travel_min": -5.0, "travel_max": 5.0}], "compensator_merit": {"operands": [{"type": "spot_rms", "path": "main", "field": 2, "target": 0.0}]}, "criteria": [{"id": "spot-edge", "operand": {"type": "spot_rms", "path": "main", "field": 2}, "max": 1.0}]}})");
  const rtt::material::MaterialLibrary lib;
  const auto end = [&](rtt::trace::Aiming aiming, rtt::trace::Aiming compensator_aiming) {
    rtt::tolerance::TrialOptions options;
    options.aiming = aiming;
    options.compensators.aiming = compensator_aiming;
    const rtt::tolerance::TrialContext context =
        rtt::tolerance::prepare(s, lib, nullptr, s.tolerances.items, options);
    const double shift[] = {0.5};
    const rtt::tolerance::TrialResult r = rtt::tolerance::evaluate(context, shift);
    REQUIRE(r.status == rtt::tolerance::TrialStatus::Valid);
    REQUIRE(r.compensators.size() == 1);
    return r.compensators[0].end;
  };
  using rtt::trace::Aiming;
  const double real = end(Aiming::Real, Aiming::Real);
  const double paraxial = end(Aiming::Paraxial, Aiming::Paraxial);
  CHECK(real != paraxial);  // content guard: another bundle, another best image distance
  CHECK(std::abs(real - paraxial) < 1.0);
  CHECK(end(Aiming::Real, Aiming::Paraxial) == real);  // compensators.aiming is replaced
}
