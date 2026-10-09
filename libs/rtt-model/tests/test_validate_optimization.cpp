// validate() for the merit function of schema 0.4.0 (ADR 0030, points 2-4; #162 part B): names,
// indices, surfaces on the path, weights, sampling and finite numbers. The evaluation (and
// merit.operand_unsupported) is #167.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "rtt/model/validate.hpp"
#include "test_support.hpp"

using namespace rtt::model;
using rtt::model::test::make_singlet;

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

/// True if `d` has exactly one diagnostic with this code, and it is at `location`.
bool only_at(const std::vector<Diagnostic>& d,
             const std::string& code,
             const std::string& location) {
  const auto n =
      std::count_if(d.begin(), d.end(), [&](const Diagnostic& x) { return x.code == code; });
  INFO(code << " at " << location << ": " << n << " diagnostics with this code");
  for (const Diagnostic& x : d) UNSCOPED_INFO(to_string(x));
  return n == 1 && std::any_of(d.begin(), d.end(), [&](const Diagnostic& x) {
           return x.code == code && x.location == location;
         });
}

OperandCommon target(double t) {
  OperandCommon c;
  c.target = t;
  return c;
}

/// The singlet of test_support.hpp (paths "main", automatic; 2 fields, 1 wavelength) with a
/// parameter row "D", configurations "near" and "far", and one operand of every type plus both
/// generators, all valid.
System with_merit() {
  System s = make_singlet();
  s.configurations = {{"near"}, {"far"}};
  ParameterRow d;
  d.name = "D";
  d.form = 5.0;
  s.parameters = {d};
  FirstOrderOperand efl{target(100.0), FirstOrderQuantity::Efl, "main", std::nullopt};
  efl.common.configuration = "far";
  RayOperand ray{
      target(0.0), RayCoordinate::Y, "main", SurfaceId("IMG"), std::nullopt, 1, 0.0, 1.0, 0};
  SpotRmsOperand spot{target(0.0), "main", 1, std::nullopt, true, SpotReference::Chief, 4};
  OpdRmsOperand opd{target(0.0), "main", 0, 0, 17};
  ParamValueOperand value{target(5.0), "D"};
  value.common.configuration = "near";
  value.common.weight = 0.0;
  s.optimization.operands = {efl, ray, spot, opd, value};
  SpotGenerator spots;
  spots.path = "main";
  spots.fields = std::vector<std::uint16_t>{0, 1};
  spots.wavelengths = std::vector<std::uint16_t>{0};
  WavefrontGenerator waves;
  waves.path = "main";
  waves.configuration = "near";
  s.optimization.generators = {spots, waves};
  return s;
}

template <class T>
T& operand(System& s, std::size_t i) {
  return std::get<T>(s.optimization.operands[i]);
}

template <class T>
T& generator(System& s, std::size_t i) {
  return std::get<T>(s.optimization.generators[i]);
}

}  // namespace

TEST_CASE("merit function: one operand of every type and both generators are valid",
          "[validate][optimization]") {
  const System s = with_merit();
  const std::vector<Diagnostic> d = validate(s);
  for (const Diagnostic& x : d) UNSCOPED_INFO(to_string(x));
  REQUIRE(d.empty());
  REQUIRE_FALSE(s.optimization.empty());
  REQUIRE(Optimization{}.empty());
}

TEST_CASE("merit function: unknown names", "[validate][optimization]") {
  System s = with_merit();
  operand<FirstOrderOperand>(s, 0).path = "other";
  REQUIRE(only_at(validate(s), "merit.unknown_path", "/optimization/operands/0/path"));

  s = with_merit();
  generator<WavefrontGenerator>(s, 1).path = "";
  REQUIRE(only_at(validate(s), "merit.unknown_path", "/optimization/generators/1/path"));

  s = with_merit();
  operand<ParamValueOperand>(s, 4).common.configuration = "middle";
  REQUIRE(only_at(validate(s), "merit.unknown_configuration",
                  "/optimization/operands/4/configuration"));

  // Without a configurations section every name is unknown (only the nominal configuration).
  s = with_merit();
  s.configurations.clear();
  const std::vector<Diagnostic> d = validate(s);
  CHECK(std::count_if(d.begin(), d.end(), [](const Diagnostic& x) {
          return x.code == "merit.unknown_configuration";
        }) == 3);

  s = with_merit();
  operand<ParamValueOperand>(s, 4).parameter = "E";
  REQUIRE(only_at(validate(s), "merit.unknown_parameter", "/optimization/operands/4/parameter"));
}

TEST_CASE("merit function: the surface of a ray operand must be on its path once",
          "[validate][optimization]") {
  System s = with_merit();
  operand<RayOperand>(s, 1).surface = SurfaceId("NONE");
  REQUIRE(only_at(validate(s), "merit.surface_not_on_path", "/optimization/operands/1/surface"));

  // An automatic path meets every surface once: occurrence 0 is fine, 1 is not.
  s = with_merit();
  operand<RayOperand>(s, 1).occurrence = 0;
  REQUIRE(validate(s).empty());
  operand<RayOperand>(s, 1).occurrence = 1;
  REQUIRE(only_at(validate(s), "merit.surface_not_on_path", "/optimization/operands/1/occurrence"));

  // An explicit double pass meets L1.S1 twice: without occurrence ambiguous, with 0 or 1 fine.
  s = with_merit();
  s.paths.push_back({"double",
                     false,
                     {{SurfaceId("STO"), EventKind::Transmit, 0},
                      {SurfaceId("L1.S1"), EventKind::Refract, 0},
                      {SurfaceId("L1.S2"), EventKind::Reflect, 0},
                      {SurfaceId("L1.S1"), EventKind::Refract, 0}}});
  RayOperand& ray = operand<RayOperand>(s, 1);
  ray.path = "double";
  ray.surface = SurfaceId("L1.S1");
  REQUIRE(only_at(validate(s), "merit.surface_ambiguous", "/optimization/operands/1/surface"));
  ray.occurrence = 1;
  REQUIRE(validate(s).empty());
  ray.occurrence = 2;
  REQUIRE(only_at(validate(s), "merit.surface_not_on_path", "/optimization/operands/1/occurrence"));
  // A surface of the system that this path does not meet.
  ray.occurrence.reset();
  ray.surface = SurfaceId("IMG");
  REQUIRE(only_at(validate(s), "merit.surface_not_on_path", "/optimization/operands/1/surface"));
}

TEST_CASE("merit function: field and wavelength indices", "[validate][optimization]") {
  // The singlet has fields 0 and 1 and wavelength 0.
  System s = with_merit();
  operand<RayOperand>(s, 1).field = 2;
  REQUIRE(only_at(validate(s), "merit.index_out_of_range", "/optimization/operands/1/field"));

  s = with_merit();
  operand<FirstOrderOperand>(s, 0).wavelength = 1;
  REQUIRE(only_at(validate(s), "merit.index_out_of_range", "/optimization/operands/0/wavelength"));

  s = with_merit();
  operand<OpdRmsOperand>(s, 3).field = 7;
  REQUIRE(only_at(validate(s), "merit.index_out_of_range", "/optimization/operands/3/field"));

  s = with_merit();
  generator<SpotGenerator>(s, 0).fields = std::vector<std::uint16_t>{0, 5};
  REQUIRE(only_at(validate(s), "merit.index_out_of_range", "/optimization/generators/0/fields/1"));

  s = with_merit();
  generator<WavefrontGenerator>(s, 1).wavelengths = std::vector<std::uint16_t>{1};
  REQUIRE(
      only_at(validate(s), "merit.index_out_of_range", "/optimization/generators/1/wavelengths/0"));
}

TEST_CASE("merit function: selections and wavelength of a polychromatic spot",
          "[validate][optimization]") {
  // The file cannot hold these (read errors), a model built in code can.
  System s = with_merit();
  generator<SpotGenerator>(s, 0).fields = std::vector<std::uint16_t>{};
  REQUIRE(only_at(validate(s), "merit.selection_empty", "/optimization/generators/0/fields"));

  s = with_merit();
  generator<WavefrontGenerator>(s, 1).wavelengths = std::vector<std::uint16_t>{};
  REQUIRE(only_at(validate(s), "merit.selection_empty", "/optimization/generators/1/wavelengths"));

  s = with_merit();
  operand<SpotRmsOperand>(s, 2).wavelength = 0;  // polychromatic is set
  REQUIRE(only_at(validate(s), "merit.polychromatic_wavelength",
                  "/optimization/operands/2/wavelength"));
}

TEST_CASE("merit function: weights, sampling and finite numbers", "[validate][optimization]") {
  for (const double w : {-1.0, kNaN, kInf}) {
    INFO("weight " << w);
    System s = with_merit();
    operand<OpdRmsOperand>(s, 3).common.weight = w;
    REQUIRE(only_at(validate(s), "merit.weight_invalid", "/optimization/operands/3/weight"));
    s = with_merit();
    generator<SpotGenerator>(s, 0).weight = w;
    REQUIRE(only_at(validate(s), "merit.weight_invalid", "/optimization/generators/0/weight"));
  }

  System s = with_merit();
  operand<SpotRmsOperand>(s, 2).rings = 0;
  REQUIRE(only_at(validate(s), "merit.sampling_invalid", "/optimization/operands/2/rings"));
  s = with_merit();
  operand<OpdRmsOperand>(s, 3).grid = -1;
  REQUIRE(only_at(validate(s), "merit.sampling_invalid", "/optimization/operands/3/grid"));
  s = with_merit();
  generator<SpotGenerator>(s, 0).rings = 0;
  REQUIRE(only_at(validate(s), "merit.sampling_invalid", "/optimization/generators/0/rings"));
  s = with_merit();
  generator<WavefrontGenerator>(s, 1).arms = 0;
  REQUIRE(only_at(validate(s), "merit.sampling_invalid", "/optimization/generators/1/arms"));

  s = with_merit();
  operand<FirstOrderOperand>(s, 0).common.target = kNaN;
  REQUIRE(only_at(validate(s), "value.not_finite", "/optimization/operands/0/target"));
  s = with_merit();
  operand<RayOperand>(s, 1).px = kInf;
  REQUIRE(only_at(validate(s), "value.not_finite", "/optimization/operands/1/px"));
  s = with_merit();
  operand<RayOperand>(s, 1).py = kNaN;
  REQUIRE(only_at(validate(s), "value.not_finite", "/optimization/operands/1/py"));
}
