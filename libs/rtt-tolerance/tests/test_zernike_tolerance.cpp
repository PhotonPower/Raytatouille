// Zernike surface errors as a tolerance kind (#316, ADR 0051, points 7 and 9, PR 5): the trial
// appends a Noll term through compile::Perturbation (decision of the maintainer in #316); its
// OPD is bit for bit that of the same term appended by hand, with the reference of the nominal
// system (ADR 0041, point 3), and coma has no paraxial model without that reference.
//
// The system: a biconvex N-BK7 singlet (L.S1 R = 50 with aperture radius 6, L.S2 R = -50) behind
// a stop, detector at the paraxial focus region; EPD 10; RMS wavefront of field 0 as criterion.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "rtt/analysis/nominal_reference.hpp"
#include "rtt/analysis/opd.hpp"
#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/model/tolerance.hpp"
#include "rtt/model/validate.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/tolerance/monte_carlo.hpp"
#include "rtt/tolerance/sensitivity.hpp"
#include "rtt/tolerance/trial.hpp"
#include "rtt/trace/sources.hpp"

using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::model::Criterion;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::FieldType;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::model::ToleranceAnchor;
using rtt::model::ToleranceDistribution;
using rtt::model::ToleranceItem;
using rtt::model::ToleranceKind;
using rtt::tolerance::TrialContext;
using rtt::tolerance::TrialStatus;

namespace {

constexpr double kRadius = 6.0;  // aperture radius of L.S1 = normalization radius, mm

const rtt::material::MaterialLibrary& library() {
  struct Schott {  // MaterialLibrary is not movable
    rtt::material::MaterialLibrary lib;
    Schott() { lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf"); }
  };
  static const Schott s;
  return s.lib;
}

Surface surface(const std::string& id, double z_mm, rtt::model::BaseShape base) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  s.shape.base = std::move(base);
  return s;
}

System singlet() {
  System s;
  s.name = "singlet";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  s.root.name = "root";
  Surface stop = surface("STO", 0.0, rtt::model::Plane{});
  stop.aperture = rtt::model::CircularAperture{8.0, 0.0};
  Surface front = surface("L.S1", 0.0, rtt::model::Conic{Param(50.0), Param(0.0)});
  front.aperture = rtt::model::CircularAperture{kRadius, 0.0};
  s.root.children = {
      {Element{"stop", ElementKind::Stop, Pose::along_z(-2.0), std::nullopt, {stop}}},
      {Element{"L",
               ElementKind::Lens,
               Pose::along_z(0.0),
               "SCHOTT:N-BK7",
               {front, surface("L.S2", 5.0, rtt::model::Conic{Param(-50.0), Param(0.0)})}}},
      {Element{"D",
               ElementKind::Detector,
               Pose::along_z(48.0),
               std::nullopt,
               {surface("IMG", 0.0, rtt::model::Plane{})}}}};
  s.paths = {{"main", true, {}}};
  rtt::model::OpdRmsOperand opd;
  opd.path = "main";
  Criterion c;
  c.id = "opd";
  c.operand = opd;
  c.max = 10.0;
  s.tolerances.criteria = {c};
  return s;
}

ToleranceItem zernike_item(const std::string& id, int term, double limit) {
  ToleranceItem i;
  i.id = id;
  i.target.anchor = ToleranceAnchor::Surface;
  i.target.name = "L.S1";
  i.target.term = static_cast<std::uint16_t>(term);
  i.kind = ToleranceKind::Zernike;
  i.limit = limit;
  i.distribution = ToleranceDistribution::Uniform;
  return i;
}

/// The nominal singlet with a_j appended by hand to L.S1 (normalization radius kRadius).
System by_hand(int term, double a) {
  System s = singlet();
  std::vector<Param> c(static_cast<std::size_t>(term - 1), Param(0.0));
  c.push_back(Param(a));
  // In place: a temporary of the term variant trips a false -Wmaybe-uninitialized in GCC 13 at
  // -O3 (the destructor of the GridSag alternative).
  std::get<Element>(s.root.children[1].value)
      .surfaces[0]
      .shape.terms.emplace_back(std::in_place_type<rtt::model::ZernikeSag>, Param(kRadius), c);
  return s;
}

CompiledSystem compiled(const System& s, const rtt::compile::Perturbation& p = {}) {
  return rtt::compile::compile(s, library(), rtt::coating::CoatingLibrary{}, 0, p);
}

std::vector<double> opd(const CompiledSystem& cs, const rtt::analysis::NominalReference& ref) {
  rtt::analysis::OpdOptions o;
  o.reference = ref;
  const auto points = rtt::analysis::opd_points(cs, PathId{0}, 0, 0, rtt::trace::GridPupil{16}, o);
  std::vector<double> w;
  for (const auto& p : points.points) {
    REQUIRE(p.status == rtt::trace::RayStatus::Alive);
    w.push_back(p.w);
  }
  return w;
}

}  // namespace

TEST_CASE("a Zernike tolerance is the term appended by hand, bit for bit (ADR 0051, point 9)",
          "[tolerance][zernike]") {
  const System nominal = singlet();
  REQUIRE_FALSE(rtt::model::has_errors(rtt::model::validate(nominal)));
  const CompiledSystem cn = compiled(nominal);
  const rtt::analysis::NominalReference ref = rtt::analysis::nominal_reference(cn, PathId{0});
  const std::vector<double> w0 = opd(cn, ref);
  for (const int j : {4, 7}) {
    const double limit = 2e-4;  // mm
    const TrialContext ctx = rtt::tolerance::prepare(
        nominal, library(), nullptr, std::vector<ToleranceItem>{zernike_item("z", j, limit)});
    REQUIRE(ctx.nominal().status == TrialStatus::Valid);
    for (const double delta : {limit, -limit}) {
      INFO("j " << j << ", delta " << delta);
      const rtt::tolerance::Applied a = rtt::tolerance::apply(ctx, std::vector<double>{delta});
      CHECK(a.system == nominal);  // the model keeps the tree of the nominal
      const std::vector<double> wt = opd(compiled(a.system, a.perturbation), ref);
      const std::vector<double> wh = opd(compiled(by_hand(j, delta)), ref);
      REQUIRE(wt.size() == wh.size());
      CHECK(wt == wh);  // bit for bit
      CHECK(wt != w0);  // content guard
      // The trial itself is valid with the nominal reference (ADR 0041, point 3).
      CHECK(rtt::tolerance::evaluate(ctx, std::vector<double>{delta}).status == TrialStatus::Valid);
    }
  }
}

TEST_CASE("the sign of a Zernike defocus error follows the rule of ADR 0051, point 9",
          "[tolerance][zernike]") {
  // a_4 > 0 on L.S1 (air n = 1 to glass n' = 1.52): dW = -(n - n')(dz - dz(0)) > 0 and grows with
  // rho^2 (Z_4 = sqrt(3)(2 rho^2 - 1)); the change at the rim exceeds that at the centre.
  const System nominal = singlet();
  const CompiledSystem cn = compiled(nominal);
  const rtt::analysis::NominalReference ref = rtt::analysis::nominal_reference(cn, PathId{0});
  rtt::analysis::OpdOptions o;
  o.reference = ref;
  const auto at = [&](const CompiledSystem& cs) {
    return rtt::analysis::opd_points(cs, PathId{0}, 0, 0, rtt::trace::GridPupil{16}, o);
  };
  const auto base = at(cn);
  const auto with = at(compiled(by_hand(4, 1e-4)));
  double centre = 0.0;
  double rim = 0.0;
  double centre_rho = 2.0;
  double rim_rho = -1.0;
  for (std::size_t k = 0; k < base.points.size(); ++k) {
    const double rho = std::hypot(base.points[k].px, base.points[k].py);
    const double dw = with.points[k].w - base.points[k].w;
    if (rho < centre_rho) {
      centre_rho = rho;
      centre = dw;
    }
    if (rho <= 1.0 && rho > rim_rho) {
      rim_rho = rho;
      rim = dw;
    }
  }
  CHECK(rim - centre > 0.0);
  // The same through the tolerance: the sensitivity of the RMS wavefront is not zero.
  const TrialContext ctx = rtt::tolerance::prepare(
      nominal, library(), nullptr, std::vector<ToleranceItem>{zernike_item("z", 4, 1e-4)});
  const auto s = rtt::tolerance::sensitivity(ctx);
  REQUIRE(s.criteria.size() == 1);
  REQUIRE(s.criteria[0].items.size() == 1);
  const auto& c = s.criteria[0].items[0];
  REQUIRE(c.plus.has_value());
  REQUIRE(c.minus.has_value());
  CHECK(c.plus.value_or(0.0) != c.minus.value_or(0.0));
}

TEST_CASE("coma as a tolerance runs with the nominal reference only (ADR 0041, point 3)",
          "[tolerance][zernike]") {
  const TrialContext ctx = rtt::tolerance::prepare(
      singlet(), library(), nullptr, std::vector<ToleranceItem>{zernike_item("z", 7, 2e-4)});
  const rtt::tolerance::Applied a = rtt::tolerance::apply(ctx, std::vector<double>{2e-4});
  const CompiledSystem trial = compiled(a.system, a.perturbation);
  CHECK_THROWS_AS((void)rtt::paraxial::first_order(trial, PathId{0}, 0),
                  rtt::paraxial::ParaxialError);
  const rtt::tolerance::TrialResult r = rtt::tolerance::evaluate(ctx, std::vector<double>{2e-4});
  INFO(r.reason);
  CHECK(r.status == TrialStatus::Valid);
  REQUIRE(r.criteria.size() == 1);
  CHECK(r.criteria[0].value.has_value());
}

TEST_CASE("Zernike items of one surface share one term; Monte Carlo is deterministic",
          "[tolerance][zernike]") {
  const std::vector<ToleranceItem> items{zernike_item("a", 4, 1e-4), zernike_item("b", 9, 5e-5),
                                         zernike_item("c", 6, 1e-4)};
  const TrialContext ctx = rtt::tolerance::prepare(singlet(), library(), nullptr, items);
  const rtt::tolerance::Applied a =
      rtt::tolerance::apply(ctx, std::vector<double>{1e-4, -5e-5, 2e-5});
  const auto it = a.perturbation.zernike.find(SurfaceId("L.S1"));
  REQUIRE(it != a.perturbation.zernike.end());
  REQUIRE(it->second.size() == 9);
  CHECK(it->second[3] == 1e-4);
  CHECK(it->second[5] == 2e-5);
  CHECK(it->second[8] == -5e-5);
  CHECK(it->second[0] == 0.0);
  rtt::tolerance::MonteCarloOptions o;
  o.seed = 316;
  o.trials = 6;
  const auto m1 = rtt::tolerance::monte_carlo(ctx, o);
  const auto m2 = rtt::tolerance::monte_carlo(ctx, o);
  REQUIRE(m1.criteria.size() == 1);
  CHECK(m1.criteria[0].mean == m2.criteria[0].mean);
  CHECK(m1.yield.good == m2.yield.good);
}

TEST_CASE("Zernike tolerance errors (ADR 0051, point 7)", "[tolerance][zernike]") {
  const auto codes = [](const System& s) {
    std::vector<std::pair<std::string, std::string>> out;
    for (const auto& d : rtt::model::validate(s)) out.emplace_back(d.code, d.location);
    return out;
  };
  const auto has = [&](const System& s, const std::string& code, const std::string& location) {
    const auto list = codes(s);
    return std::find(list.begin(), list.end(), std::pair{code, location}) != list.end();
  };
  const std::string item = "/tolerances/items/0";
  System s = singlet();
  for (const int j : {1, 232}) {
    INFO("term " << j);
    s.tolerances.items = {zernike_item("z", j, 1e-4)};
    CHECK(has(s, "tolerances.target_invalid", item + "/term"));
  }
  s.tolerances.items = {zernike_item("z", 4, 1e-4)};
  s.tolerances.items[0].target.term.reset();
  CHECK(has(s, "tolerances.target_invalid", item));
  // Without an aperture the semi-diameter is not defined.
  s.tolerances.items = {zernike_item("z", 4, 1e-4)};
  s.tolerances.items[0].target.name = "L.S2";
  CHECK(has(s, "tolerances.target_invalid", item + "/surface"));
  // term only with the kind zernike.
  s.tolerances.items = {zernike_item("z", 4, 1e-4)};
  s.tolerances.items[0].kind = ToleranceKind::Conic;
  CHECK(has(s, "tolerances.target_invalid", item + "/term"));
  // Zernike needs a surface.
  s.tolerances.items = {zernike_item("z", 4, 1e-4)};
  s.tolerances.items[0].target.anchor = ToleranceAnchor::Node;
  s.tolerances.items[0].target.name = "L";
  CHECK(has(s, "tolerances.target_invalid", item + "/kind"));
  // No class generates it: a single value is an error (use an explicit item).
  s.tolerances.items.clear();
  rtt::model::ToleranceValueOverride v;
  v.target.anchor = ToleranceAnchor::Surface;
  v.target.name = "L.S1";
  v.kind = ToleranceKind::Zernike;
  v.limit = 1e-4;
  s.tolerances.overrides = {v};
  CHECK(has(s, "tolerances.target_invalid", "/tolerances/overrides/0/kind"));
  // The same coefficient twice is one quantity twice.
  s.tolerances.overrides.clear();
  s.tolerances.items = {zernike_item("a", 4, 1e-4), zernike_item("b", 4, 1e-4)};
  CHECK(has(s, "tolerances.quantity_duplicate", "/tolerances/items/1"));
  s.tolerances.items = {zernike_item("a", 4, 1e-4), zernike_item("b", 5, 1e-4)};
  CHECK_FALSE(has(s, "tolerances.quantity_duplicate", "/tolerances/items/1"));
}

TEST_CASE("a Zernike item leaves Phi_0 and the other items bit for bit unchanged",
          "[tolerance][zernike]") {
  // At delta = 0 the overlay appends nothing (like an identity motion): the nominal keeps the
  // closed-form intersection of its conic.
  ToleranceItem conic;
  conic.id = "k";
  conic.target.anchor = ToleranceAnchor::Surface;
  conic.target.name = "L.S2";
  conic.kind = ToleranceKind::Conic;
  conic.limit = 0.01;
  const TrialContext without =
      rtt::tolerance::prepare(singlet(), library(), nullptr, std::vector<ToleranceItem>{conic});
  const TrialContext with = rtt::tolerance::prepare(
      singlet(), library(), nullptr, std::vector<ToleranceItem>{conic, zernike_item("z", 4, 1e-4)});
  REQUIRE(without.nominal().criteria[0].value.has_value());
  CHECK(with.nominal().criteria[0].value == without.nominal().criteria[0].value);
  CHECK(rtt::tolerance::evaluate(with, std::vector<double>{0.01, 0.0}).criteria[0].value ==
        rtt::tolerance::evaluate(without, std::vector<double>{0.01}).criteria[0].value);
  const rtt::tolerance::Applied zero = rtt::tolerance::apply(with, std::vector<double>{0.0, 0.0});
  const CompiledSystem cs = compiled(zero.system, zero.perturbation);
  const std::optional<std::uint32_t> s1 = cs.find_surface(SurfaceId("L.S1"));
  REQUIRE(s1.has_value());
  CHECK(std::holds_alternative<rtt::geom::Conic<double>>(cs.surfaces()[s1.value_or(0)].shape));
}
