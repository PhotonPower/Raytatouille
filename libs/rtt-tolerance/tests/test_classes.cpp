// Tolerance classes (ADR 0035, point 5; ADR 0036, point 9; #260): the class table and expand().
//
// The numbers of the class table are the decision of the maintainer (plan #260, section 3;
// ADR 0036, addendum #260); the row abbe is relative: a share of nu_d of the glass (SCHOTT
// TIE-29, table 5-1, steps 4/3/1), so expand needs the materials.

#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <array>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/perturbation.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/model/tolerance.hpp"
#include "rtt/model/tolerance_classes.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/tolerance/classes.hpp"
#include "rtt/tolerance/expand.hpp"
#include "rtt/tolerance/random.hpp"

namespace fs = std::filesystem;

using rtt::model::ClassSkip;
using rtt::model::ToleranceAnchor;
using rtt::model::ToleranceClass;
using rtt::model::ToleranceDistribution;
using rtt::model::ToleranceItem;
using rtt::model::ToleranceKind;
using rtt::model::ToleranceLevel;
using rtt::tolerance::class_limit;
using rtt::tolerance::expand;
using rtt::tolerance::ItemOrigin;
using rtt::tolerance::ToleranceItems;

namespace {

using K = ToleranceKind;
using L = ToleranceLevel;

struct Combo {
  ToleranceKind kind;
  ToleranceLevel level;
  auto operator<=>(const Combo&) const = default;
};

/// The 15 combinations a class generates (ADR 0035, point 5): surface 3, element 5 plus the two
/// material kinds, assembly 5.
constexpr std::array<Combo, 15> kCombos{{
    {K::CurvatureRelative, L::Surface},
    {K::Conic, L::Surface},
    {K::Thickness, L::Surface},
    {K::Thickness, L::Element},
    {K::DecenterX, L::Element},
    {K::DecenterY, L::Element},
    {K::TiltXDeg, L::Element},
    {K::TiltYDeg, L::Element},
    {K::Index, L::Element},
    {K::Abbe, L::Element},
    {K::DecenterX, L::Assembly},
    {K::DecenterY, L::Assembly},
    {K::ShiftZ, L::Assembly},
    {K::TiltXDeg, L::Assembly},
    {K::TiltYDeg, L::Assembly},
}};

constexpr std::array<ToleranceClass, 3> kClasses{ToleranceClass::Coarse, ToleranceClass::Normal,
                                                 ToleranceClass::Precise};

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

rtt::model::System reference() {
  return rtt::io::load_system(fs::path(RTT_REFERENCE_DIR) / "m7" / "achromat_tolerances.rtt.json");
}

std::vector<fs::path> reference_files() {
  std::vector<fs::path> files;
  for (const auto& entry : fs::recursive_directory_iterator(RTT_REFERENCE_DIR)) {
    const std::string name = entry.path().filename().string();
    if (entry.is_regular_file() && name.ends_with(".rtt.json")) files.push_back(entry.path());
  }
  std::sort(files.begin(), files.end());
  return files;
}

bool has_errors(const rtt::model::System& system) {
  const std::vector<rtt::model::Diagnostic> d = rtt::model::validate(system);
  return std::any_of(d.begin(), d.end(), [](const rtt::model::Diagnostic& x) {
    return x.severity == rtt::model::Severity::Error;
  });
}

std::map<std::string, ToleranceItem> by_id(const ToleranceItems& r) {
  std::map<std::string, ToleranceItem> out;
  for (const ToleranceItem& i : r.items) out.emplace(i.id, i);
  return out;
}

/// By value: bound to a reference, GCC takes the result for one into the temporary string of a
/// literal argument (-Wdangling-reference, AGENTS.md).
rtt::model::ClassTolerance find_class(const ToleranceItems& r, const std::string& id) {
  const auto it = std::find_if(r.classes.items.begin(), r.classes.items.end(),
                               [&](const rtt::model::ClassTolerance& c) { return c.id == id; });
  REQUIRE(it != r.classes.items.end());
  return *it;
}

/// The glass catalogue of the reference achromat (N-BK7, F2).
const rtt::material::MaterialLibrary& schott() {
  struct Schott {  // MaterialLibrary is not movable
    rtt::material::MaterialLibrary lib;
    Schott() { lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf"); }
  };
  static const Schott instance;
  return instance.lib;
}

/// nu_d of `reference` in the environment of `system`, at run time from the catalogue: the
/// material at the vacuum wavelengths of the lines d, F, C (TIE-29, Eq. (2.1-1)), computed here
/// without the helper of rtt-compile.
double nu_d(const rtt::model::System& system, const std::string& reference) {
  const auto m = schott().resolve(reference);
  const auto at = [&](double air_um) {
    return m
        ->index(rtt::compile::line_vacuum_wavelength_um(air_um), system.environment.temperature_c,
                system.environment.pressure_atm)
        .real();
  };
  const double nd = at(rtt::compile::kLineDAirUm);
  const double nf = at(rtt::compile::kLineFAirUm);
  const double nc = at(rtt::compile::kLineCAirUm);
  return (nd - 1.0) / (nf - nc);
}

/// The glass of segment `segment` of L1 in the reference file.
std::string l1_glass(std::uint32_t segment) {
  return segment == 0 ? "SCHOTT:N-BK7" : "SCHOTT:F2";
}

rtt::model::ToleranceTarget node(const std::string& name) {
  rtt::model::ToleranceTarget t;
  t.anchor = ToleranceAnchor::Node;
  t.name = name;
  return t;
}

}  // namespace

TEST_CASE("classes: every generated combination has a limit, coarse > normal > precise",
          "[tolerance][classes]") {
  for (const Combo& c : kCombos) {
    INFO("kind " << rtt::model::tolerance_kind_name(c.kind) << ", level "
                 << static_cast<int>(c.level));
    const double coarse = class_limit(ToleranceClass::Coarse, c.kind, c.level);
    const double normal = class_limit(ToleranceClass::Normal, c.kind, c.level);
    const double precise = class_limit(ToleranceClass::Precise, c.kind, c.level);
    CHECK(std::isfinite(coarse));
    CHECK(precise > 0.0);
    CHECK(coarse > normal);
    CHECK(normal > precise);
    CHECK_THROWS_AS((void)class_limit(ToleranceClass::None, c.kind, c.level),
                    std::invalid_argument);
  }
  // The list above is the one of rtt-model: the reference file generates exactly these
  // combinations (skipped ones included), so a kind added there fails here.
  const ToleranceItems r = expand(reference(), schott());
  std::set<Combo> generated;
  for (const rtt::model::ClassTolerance& c : r.classes.items) generated.insert({c.kind, c.level});
  CHECK(generated == std::set<Combo>(kCombos.begin(), kCombos.end()));
}

TEST_CASE("classes: the values decided by the maintainer", "[tolerance][classes]") {
  // ADR 0036, addendum #260 (decision of the maintainer on plan #260, section 3). Units: mm,
  // degree; curvature_relative, conic and index dimensionless; abbe a share of nu_d.
  struct Value {
    Combo combo;
    std::array<double, 3> limits;  // coarse, normal, precise
  };
  const std::array<Value, 15> table{{
      {{K::CurvatureRelative, L::Surface}, {5e-3, 2e-3, 5e-4}},
      {{K::Conic, L::Surface}, {0.05, 0.02, 0.005}},
      {{K::Thickness, L::Surface}, {0.25, 0.1, 0.025}},
      {{K::Thickness, L::Element}, {0.25, 0.1, 0.05}},
      {{K::DecenterX, L::Element}, {0.25, 0.1, 0.05}},
      {{K::DecenterY, L::Element}, {0.25, 0.1, 0.05}},
      {{K::TiltXDeg, L::Element}, {0.2, 0.1, 0.05}},
      {{K::TiltYDeg, L::Element}, {0.2, 0.1, 0.05}},
      {{K::Index, L::Element}, {1e-3, 5e-4, 2e-4}},
      {{K::Abbe, L::Element}, {0.008, 0.005, 0.002}},  // TIE-29 steps 4, 3, 1
      {{K::DecenterX, L::Assembly}, {0.25, 0.1, 0.05}},
      {{K::DecenterY, L::Assembly}, {0.25, 0.1, 0.05}},
      {{K::ShiftZ, L::Assembly}, {0.25, 0.1, 0.05}},
      {{K::TiltXDeg, L::Assembly}, {0.2, 0.1, 0.05}},
      {{K::TiltYDeg, L::Assembly}, {0.2, 0.1, 0.05}},
  }};
  for (const Value& v : table) {
    INFO("kind " << rtt::model::tolerance_kind_name(v.combo.kind) << ", level "
                 << static_cast<int>(v.combo.level));
    for (std::size_t c = 0; c < kClasses.size(); ++c) {
      CHECK(same(class_limit(kClasses[c], v.combo.kind, v.combo.level), v.limits[c]));
    }
  }
}

TEST_CASE("classes: abbe_share is the abbe row of the table", "[tolerance][classes]") {
  // #261 (P3-2 of the review of #275, decided by the maintainer): the share of nu_d under its
  // own name; class_limit keeps giving the same number.
  const std::array<double, 3> shares{0.008, 0.005, 0.002};  // TIE-29 steps 4, 3, 1
  for (std::size_t c = 0; c < kClasses.size(); ++c) {
    CHECK(same(rtt::tolerance::abbe_share(kClasses[c]), shares[c]));
    CHECK(same(rtt::tolerance::abbe_share(kClasses[c]),
               class_limit(kClasses[c], K::Abbe, L::Element)));
  }
  CHECK_THROWS_AS((void)rtt::tolerance::abbe_share(ToleranceClass::None), std::invalid_argument);
  CHECK_THROWS_AS((void)rtt::tolerance::abbe_share(static_cast<ToleranceClass>(9)),
                  std::invalid_argument);
}

TEST_CASE("classes: kinds no class generates at a level have no limit", "[tolerance][classes]") {
  const std::array<Combo, 9> invalid{{
      {K::Radius, L::Surface},
      {K::Param, L::Surface},
      {K::TiltXDeg, L::Surface},
      {K::TiltZDeg, L::Element},
      {K::ShiftZ, L::Element},
      {K::Parameter, L::Element},
      {K::Thickness, L::Assembly},
      {K::Index, L::Assembly},
      {K::TiltZDeg, L::Assembly},
  }};
  for (const ToleranceClass c : kClasses) {
    for (const Combo& x : invalid) {
      INFO("kind " << rtt::model::tolerance_kind_name(x.kind) << ", level "
                   << static_cast<int>(x.level));
      CHECK_THROWS_AS((void)class_limit(c, x.kind, x.level), std::invalid_argument);
    }
  }
  CHECK_THROWS_AS((void)class_limit(static_cast<ToleranceClass>(9), K::Thickness, L::Element),
                  std::invalid_argument);
  CHECK_THROWS_AS(
      (void)class_limit(ToleranceClass::Normal, static_cast<ToleranceKind>(99), L::Element),
      std::invalid_argument);
  CHECK_THROWS_AS(
      (void)class_limit(ToleranceClass::Normal, K::Thickness, static_cast<ToleranceLevel>(7)),
      std::invalid_argument);
}

TEST_CASE("expand: class, class override, single values and explicit items",
          "[tolerance][classes]") {
  // tests/reference/m7/achromat_tolerances.rtt.json: class normal; objective precise (and so L1
  // and its surfaces); image none; single values curvature_relative 0.005 at L1.S2 and index
  // 0.001 uniform at L1 segment 1; explicit S3-wedge, S1-conic (param on /shape/base/conic of
  // L1.S1) and R1. L1.S1 has its radius bound to the row R1 and conic 0; STO is plane; stop
  // has no material in the beam.
  const rtt::model::System system = reference();
  REQUIRE(system.tolerances.tolerance_class == ToleranceClass::Normal);
  REQUIRE(system.tolerances.items.size() == 3);
  const ToleranceItems r = expand(system, schott());
  // Content guard: a crown and a flint (V(d) of the NM lines in tests/catalogs/schott.agf).
  REQUIRE(std::abs(nu_d(system, l1_glass(0)) - 64.17) < 0.1);
  REQUIRE(std::abs(nu_d(system, l1_glass(1)) - 36.37) < 0.1);

  struct Expected {
    std::string id;
    ToleranceClass tolerance_class;
    ToleranceLevel level;
  };
  constexpr ToleranceClass n = ToleranceClass::Normal;
  constexpr ToleranceClass p = ToleranceClass::Precise;
  const std::vector<Expected> generated{
      {"class:node:stop:thickness", n, L::Element},
      {"class:node:stop:decenter_x", n, L::Element},
      {"class:node:stop:decenter_y", n, L::Element},
      {"class:node:stop:tilt_x_deg", n, L::Element},
      {"class:node:stop:tilt_y_deg", n, L::Element},
      {"class:node:objective:decenter_x", p, L::Assembly},
      {"class:node:objective:decenter_y", p, L::Assembly},
      {"class:node:objective:shift_z", p, L::Assembly},
      {"class:node:objective:tilt_x_deg", p, L::Assembly},
      {"class:node:objective:tilt_y_deg", p, L::Assembly},
      {"class:node:L1:thickness", p, L::Element},
      {"class:node:L1:decenter_x", p, L::Element},
      {"class:node:L1:decenter_y", p, L::Element},
      {"class:node:L1:tilt_x_deg", p, L::Element},
      {"class:node:L1:tilt_y_deg", p, L::Element},
      {"class:node:L1:index:0", p, L::Element},
      {"class:node:L1:index:1", p, L::Element},
      {"class:node:L1:abbe:0", p, L::Element},
      {"class:node:L1:abbe:1", p, L::Element},
      {"class:surface:L1.S2:curvature_relative", p, L::Surface},
      {"class:surface:L1.S2:thickness", p, L::Surface},
      {"class:surface:L1.S3:curvature_relative", p, L::Surface},
      {"class:surface:L1.S3:thickness", p, L::Surface},
  };
  REQUIRE(r.items.size() == generated.size() + 3);
  REQUIRE(r.origins.size() == r.items.size());

  for (std::size_t i = 0; i < generated.size(); ++i) {
    const Expected& e = generated[i];
    const ToleranceItem& item = r.items[i];
    INFO("item " << i << ": " << item.id << ", expected " << e.id);
    REQUIRE(item.id == e.id);
    REQUIRE(r.origins[i].source == ItemOrigin::Source::Class);
    REQUIRE(r.origins[i].index < r.classes.items.size());
    const rtt::model::ClassTolerance& c = r.classes.items[r.origins[i].index];
    CHECK(c.id == e.id);
    CHECK(c.generated());
    CHECK(c.tolerance_class == e.tolerance_class);
    CHECK(c.level == e.level);
    CHECK(item.target == c.target);
    CHECK(item.kind == c.kind);
    // Limit and distribution: the single value, else the class table and normal.
    if (e.id == "class:surface:L1.S2:curvature_relative") {
      CHECK(same(item.limit, 0.005));
      CHECK(item.distribution == ToleranceDistribution::Normal);
    } else if (e.id == "class:node:L1:index:1") {
      CHECK(same(item.limit, 0.001));
      CHECK(item.distribution == ToleranceDistribution::Uniform);
    } else if (c.kind == K::Abbe) {
      // A share of nu_d of the glass of the segment (TIE-29, table 5-1).
      REQUIRE(c.target.segment.has_value());
      const double vd = nu_d(system, l1_glass(c.target.segment.value_or(0)));
      CHECK(same(item.limit, class_limit(e.tolerance_class, c.kind, e.level) * vd));
      CHECK(item.distribution == ToleranceDistribution::Normal);
    } else {
      CHECK(same(item.limit, class_limit(e.tolerance_class, c.kind, e.level)));
      CHECK(item.distribution == ToleranceDistribution::Normal);
    }
  }
  // The explicit items follow unchanged, in file order.
  for (std::size_t i = 0; i < 3; ++i) {
    const std::size_t k = generated.size() + i;
    CHECK(r.items[k] == system.tolerances.items[i]);
    CHECK(r.origins[k] == ItemOrigin{ItemOrigin::Source::Explicit, i});
  }

  // Skipped and replaced kinds stay in `classes`, with the reason.
  CHECK(find_class(r, "class:node:stop:index").skip == ClassSkip::NoMaterialInBeam);
  CHECK(find_class(r, "class:surface:STO:curvature_relative").skip == ClassSkip::Plane);
  CHECK(find_class(r, "class:surface:L1.S1:curvature_relative").skip == ClassSkip::Bound);
  const rtt::model::ClassTolerance conic = find_class(r, "class:surface:L1.S1:conic");
  CHECK(conic.skip == ClassSkip::Sphere);
  CHECK(conic.replaced_by == std::optional<std::size_t>(1));  // S1-conic
  // image is under the class none: nothing listed.
  CHECK(std::none_of(r.classes.items.begin(), r.classes.items.end(),
                     [](const rtt::model::ClassTolerance& c) { return c.target.name == "image"; }));
  REQUIRE(r.classes.value_overrides.size() == 2);
  CHECK(r.classes.value_overrides[0].effect == rtt::model::OverrideEffect::Applied);
  CHECK(r.classes.value_overrides[1].effect == rtt::model::OverrideEffect::Applied);
}

TEST_CASE("expand: without a class only the explicit items", "[tolerance][classes]") {
  rtt::model::System system = reference();
  system.tolerances.tolerance_class = ToleranceClass::None;
  std::erase_if(system.tolerances.overrides, [](const rtt::model::ToleranceOverride& o) {
    return std::holds_alternative<rtt::model::ToleranceClassOverride>(o);
  });
  REQUIRE(!has_errors(system));
  const ToleranceItems r = expand(system, rtt::material::MaterialLibrary{});  // no glass needed
  CHECK(r.classes.items.empty());
  REQUIRE(r.items.size() == 3);
  for (std::size_t i = 0; i < 3; ++i) {
    CHECK(r.items[i] == system.tolerances.items[i]);
    CHECK(r.origins[i] == ItemOrigin{ItemOrigin::Source::Explicit, i});
  }
}

TEST_CASE("expand: an added override leaves the other items and their draws unchanged",
          "[tolerance][classes]") {
  const rtt::model::System before = reference();
  rtt::model::System after = before;
  // L1 (and its surfaces) coarse instead of precise; a single value at the objective.
  after.tolerances.overrides.emplace_back(
      rtt::model::ToleranceClassOverride{"L1", ToleranceClass::Coarse});
  after.tolerances.overrides.emplace_back(rtt::model::ToleranceValueOverride{
      node("objective"), K::DecenterX, 0.07, ToleranceDistribution::Uniform});
  REQUIRE(!has_errors(after));
  const ToleranceItems a = expand(before, schott());
  const ToleranceItems b = expand(after, schott());
  const std::map<std::string, ToleranceItem> ia = by_id(a);
  const std::map<std::string, ToleranceItem> ib = by_id(b);
  REQUIRE(ia.size() == ib.size());

  const auto affected = [](const std::string& id) {
    return id.starts_with("class:node:L1:") || id.starts_with("class:surface:L1.") ||
           id == "class:node:objective:decenter_x";
  };
  const std::uint64_t seed = 20261012;  // monte_carlo.seed of the reference file
  std::size_t unchanged = 0;
  std::size_t changed = 0;
  for (const auto& [id, x] : ia) {
    INFO(id);
    const auto it = ib.find(id);
    REQUIRE(it != ib.end());
    const ToleranceItem& y = it->second;
    if (!affected(id)) {
      ++unchanged;
      CHECK(x == y);
      bool draws_same = true;
      for (std::uint32_t t = 0; t < 200; ++t) {
        draws_same = draws_same && same(rtt::tolerance::perturbation(seed, t, x),
                                        rtt::tolerance::perturbation(seed, t, y));
      }
      CHECK(draws_same);
      continue;
    }
    ++changed;
    CHECK(x.target == y.target);
    CHECK(x.kind == y.kind);
    if (id == "class:node:objective:decenter_x") {
      CHECK(same(y.limit, 0.07));
      CHECK(y.distribution == ToleranceDistribution::Uniform);
    } else if (id == "class:surface:L1.S2:curvature_relative") {
      CHECK(same(y.limit, 0.005));  // the single value goes before the class override
    } else if (id == "class:node:L1:index:1") {
      CHECK(same(y.limit, 0.001));
    } else {
      const rtt::model::ClassTolerance c = find_class(b, id);
      CHECK(c.tolerance_class == ToleranceClass::Coarse);
      const double share = class_limit(ToleranceClass::Coarse, y.kind, c.level);
      CHECK(same(y.limit, y.kind == K::Abbe
                              ? share * nu_d(after, l1_glass(c.target.segment.value_or(0)))
                              : share));
      CHECK(y.limit > x.limit);  // coarse instead of precise
    }
  }
  // stop 5, objective 4, explicit 3; L1 9, L1.S2 2, L1.S3 2, objective:decenter_x 1.
  CHECK(unchanged == 12U);
  CHECK(changed == 14U);
}

TEST_CASE("expand: the draws of a class set are reproducible for every number of threads",
          "[tolerance][classes]") {
  const rtt::model::System system = reference();
  const ToleranceItems first = expand(system, schott());
  const ToleranceItems second = expand(system, schott());
  REQUIRE(first.items == second.items);
  REQUIRE(first.items.size() > system.tolerances.items.size());  // the class generates items
  const std::uint64_t seed = 20261012;
  constexpr std::uint32_t kN = 500;
  const std::size_t m = first.items.size();
  std::vector<double> serial;
  for (std::uint32_t t = 0; t < kN; ++t) {
    const std::vector<double> v = rtt::tolerance::perturbations(seed, t, first.items);
    serial.insert(serial.end(), v.begin(), v.end());
  }
  // Content guard: the values are not all equal (e.g. all 0).
  CHECK(std::adjacent_find(serial.begin(), serial.end(), std::not_equal_to<>()) != serial.end());
  for (const int threads : {1, 4, oneapi::tbb::task_arena::automatic}) {
    INFO("threads " << threads);
    std::vector<double> out(serial.size());
    oneapi::tbb::task_arena arena(threads);
    arena.execute([&] {
      oneapi::tbb::parallel_for(std::uint32_t{0}, kN, [&](std::uint32_t t) {
        const std::vector<double> v = rtt::tolerance::perturbations(seed, t, second.items);
        std::copy(v.begin(), v.end(), out.begin() + static_cast<std::ptrdiff_t>(t * m));
      });
    });
    bool all_same = true;
    for (std::size_t i = 0; i < serial.size(); ++i) all_same = all_same && same(out[i], serial[i]);
    CHECK(all_same);
  }
}

TEST_CASE("expand: a class never makes a valid reference file invalid", "[tolerance][classes]") {
  // ADR 0035, point 5: class alone never invalidates a valid file (up to a hash collision of two
  // generated ids). Every reference file that compiles, with each class; the glasses of a file
  // come from the catalogue it compiles with (as in the compile tests of the overlay).
  rtt::material::MaterialLibrary m1;
  m1.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  rtt::material::MaterialLibrary m2;
  m2.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf");
  rtt::coating::CoatingLibrary coatings;
  coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/demo.json");
  coatings.add_catalog(std::string(RTT_CATALOG_DIR) + "/coatings/m3.json");
  std::size_t files = 0;
  std::size_t with_generated = 0;
  for (const fs::path& file : reference_files()) {
    rtt::model::System system = rtt::io::load_system(file);
    const rtt::material::MaterialLibrary* lib = nullptr;
    for (const rtt::material::MaterialLibrary* candidate : {&m1, &m2}) {
      try {
        (void)rtt::compile::compile(system, *candidate, coatings);
      } catch (const rtt::compile::CompileError&) {
        continue;
      }
      lib = candidate;
      break;
    }
    if (lib == nullptr) continue;
    ++files;
    for (const ToleranceClass c : kClasses) {
      INFO(file.string() << ", class " << static_cast<int>(c));
      system.tolerances.tolerance_class = c;
      ToleranceItems r;
      REQUIRE_NOTHROW(r = expand(system, *lib));
      REQUIRE(r.origins.size() == r.items.size());
      const auto generated = static_cast<std::size_t>(
          std::count_if(r.classes.items.begin(), r.classes.items.end(),
                        [](const rtt::model::ClassTolerance& x) { return x.generated(); }));
      CHECK(r.items.size() == generated + system.tolerances.items.size());
      for (const ToleranceItem& i : r.items) {
        CHECK(std::isfinite(i.limit));
        CHECK(i.limit > 0.0);
      }
      if (c == ToleranceClass::Normal && generated > 0) ++with_generated;
    }
  }
  CHECK(files >= 20U);
  CHECK(with_generated == files);  // every system has at least one element
}

TEST_CASE("expand: the glass of a generated abbe item must resolve", "[tolerance][classes]") {
  // Decision of the maintainer on #260: an unknown material is a CompileError (material.unknown)
  // at the material of the element, as in compile; only the abbe items need the glass.
  const rtt::model::System system = reference();
  try {
    (void)expand(system, rtt::material::MaterialLibrary{});
    FAIL("no CompileError");
  } catch (const rtt::compile::CompileError& e) {
    std::vector<std::pair<std::string, std::string>> got;
    for (const rtt::model::Diagnostic& d : e.diagnostics()) got.emplace_back(d.code, d.location);
    // L1 is child 0 of the objective (child 1 of the root); abbe:0 and abbe:1.
    const std::vector<std::pair<std::string, std::string>> expected{
        {"material.unknown", "/root/children/1/children/0/material/0"},
        {"material.unknown", "/root/children/1/children/0/material/1"}};
    CHECK(got == expected);
  }
}

TEST_CASE("expand: abbe of a glass without dispersion is skipped", "[tolerance][classes]") {
  // Decision of the maintainer on #260: without dispersion (n_F = n_C) there is no nu_d; the
  // class skips abbe there with the reason ConstAbbe, as for a CONST material.
  class Flat final : public rtt::material::Material {
   public:
    [[nodiscard]] rtt::math::Complex index(double, double, double) const override {
      return {1.6, 0.0};
    }
  };
  rtt::material::MaterialLibrary lib;
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  lib.add("TEST:FLAT", std::make_shared<const Flat>());
  rtt::model::System system = reference();
  auto& objective = std::get<rtt::model::Assembly>(system.root.children[1].value);
  auto& l1 = std::get<rtt::model::Element>(objective.children[0].value);
  REQUIRE(l1.name == "L1");
  REQUIRE(l1.segment_materials.size() == 2);
  l1.segment_materials[0] = "TEST:FLAT";
  const ToleranceItems r = expand(system, lib);
  const rtt::model::ClassTolerance flat = find_class(r, "class:node:L1:abbe:0");
  CHECK(flat.skip == ClassSkip::ConstAbbe);
  CHECK_FALSE(flat.generated());
  const std::map<std::string, ToleranceItem> items = by_id(r);
  CHECK(items.count("class:node:L1:abbe:0") == 0U);
  // Content guard: the flint keeps its abbe item, the flat glass its index item.
  REQUIRE(items.count("class:node:L1:abbe:1") == 1U);
  CHECK(
      same(items.at("class:node:L1:abbe:1").limit,
           class_limit(ToleranceClass::Precise, K::Abbe, L::Element) * nu_d(system, l1_glass(1))));
  CHECK(items.count("class:node:L1:index:0") == 1U);
}

TEST_CASE("expand: an invalid tolerance set throws CompileError with all errors",
          "[tolerance][classes]") {
  rtt::model::System system = reference();
  // thickness at surface 0 (no class generates it) and an unknown target.
  rtt::model::ToleranceTarget s0;
  s0.anchor = ToleranceAnchor::Surface;
  s0.name = "L1.S1";
  system.tolerances.overrides.emplace_back(
      rtt::model::ToleranceValueOverride{s0, K::Thickness, 0.01, ToleranceDistribution::Normal});
  system.tolerances.overrides.emplace_back(rtt::model::ToleranceValueOverride{
      node("nowhere"), K::DecenterX, 0.01, ToleranceDistribution::Normal});
  try {
    (void)expand(system, schott());
    FAIL("no CompileError");
  } catch (const rtt::compile::CompileError& e) {
    std::vector<std::pair<std::string, std::string>> got;
    for (const rtt::model::Diagnostic& d : e.diagnostics()) {
      CHECK(d.severity == rtt::model::Severity::Error);
      got.emplace_back(d.code, d.location);
    }
    const std::pair<std::string, std::string> invalid{"tolerances.target_invalid",
                                                      "/tolerances/overrides/5/kind"};
    const std::pair<std::string, std::string> unknown{"tolerances.unknown_target",
                                                      "/tolerances/overrides/6/node"};
    CHECK(std::find(got.begin(), got.end(), invalid) != got.end());
    CHECK(std::find(got.begin(), got.end(), unknown) != got.end());
  }
}
