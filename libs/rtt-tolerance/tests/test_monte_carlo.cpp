// Monte Carlo (ADR 0036, point 7; #259). The systems are those of test_sensitivity.cpp: ideal
// lenses (ADR 0031) in VACUUM, where the image position is linear in a decenter, so the
// distribution of a criterion follows from that of the perturbation. Every band is derived in the
// comment of its test case before the run; the seeds are fixed and were not chosen.
//
// Bands of the statistics (docs/quellen.md):
// - mean: |mean - mu| <= 4 sigma / sqrt(N) (NIST e-Handbook, Sec. 1.3.5.2, with the known sigma);
// - variance: |s^2 - sigma^2| <= 4 sqrt(Var(s^2)), Var(s^2) = (mu_4 - (N - 3) / (N - 1) sigma^4) /
// N
//   (MathWorld, "Sample Variance Distribution", Eq. (22), for the estimator with N - 1). This
//   assumes a nearly normal s^2; it is not a Chebyshev band;
// - yield: the true fraction lies in the Wilson interval with z = z_{1 - 1e-4 / 2} = 3.8906,
//   computed here from erfc (NIST DLMF, Eq. (7.2.2)). With 95 % a correct code would fail for one
//   seed in twenty; with 1e-4 for one in 10^4.
// Bound for traced coordinates: kEps = 1e-12 mm, as in test_sensitivity.cpp.

#include <oneapi/tbb/task_arena.h>

#include <algorithm>
#include <atomic>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "monte_carlo_detail.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/model.hpp"
#include "rtt/tolerance/monte_carlo.hpp"
#include "rtt/tolerance/random.hpp"
#include "rtt/tolerance/sensitivity.hpp"
#include "rtt/tolerance/trial.hpp"
#include "rtt/trace/run_control.hpp"

using rtt::model::Criterion;
using rtt::model::Element;
using rtt::model::ElementKind;
using rtt::model::FieldType;
using rtt::model::Param;
using rtt::model::Pose;
using rtt::model::RayCoordinate;
using rtt::model::Surface;
using rtt::model::SurfaceId;
using rtt::model::System;
using rtt::model::ToleranceDistribution;
using rtt::model::ToleranceItem;
using rtt::model::ToleranceKind;
using rtt::tolerance::CriterionStatistics;
using rtt::tolerance::MonteCarloOptions;
using rtt::tolerance::MonteCarloResult;
using rtt::tolerance::MonteCarloStatus;
using rtt::tolerance::TrialContext;
using rtt::tolerance::TrialRecord;
using rtt::tolerance::TrialStatus;

namespace {

constexpr double kEps = 1e-12;  // mm, see the head of the file

bool same(double a, double b) {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

Surface plane(const std::string& id, double z_mm, double aperture_radius) {
  Surface s;
  s.id = SurfaceId(id);
  s.pose = Pose::along_z(z_mm);
  s.aperture = rtt::model::CircularAperture{aperture_radius, 0.0};
  return s;
}

/// As in test_sensitivity.cpp: in VACUUM at 0.5876 um, EPD 10 mm, the stop "stop" (STO) at
/// z = -1, the ideal lens "lens" (IL, f) at z = 0; object at infinity (fields 0 and 0.5 degree)
/// or at `object_distance` (object heights 0 and 1 mm, also the lens's design conjugate).
System lens_system(double f, std::optional<double> object_distance = std::nullopt) {
  System s;
  s.name = "monte_carlo";
  s.environment.medium = "VACUUM";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {rtt::model::SystemApertureType::EntrancePupilDiameter, Param(10.0)};
  if (object_distance) {
    s.object.at_infinity = false;
    s.object.distance = Param(*object_distance);
    s.fields = {FieldType::ObjectHeight, {{0.0, 0.0, 1.0}, {0.0, 1.0, 1.0}}};
  } else {
    s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 0.5, 1.0}}};
  }
  s.root.name = "root";
  Surface lens = plane("IL", 0.0, 30.0);
  std::optional<Param> conjugate;
  if (object_distance) conjugate = Param(*object_distance);
  lens.interaction = rtt::model::IdealLens{Param(f), conjugate};
  s.root.children = {
      {Element{
          "stop", ElementKind::Stop, Pose::along_z(-1.0), std::nullopt, {plane("STO", 0.0, 10.0)}}},
      {Element{"lens", ElementKind::ThinElement, Pose::along_z(0.0), std::nullopt, {lens}}}};
  return s;
}

/// Adds the detector "detector" (IMG) at z and the automatic path "main".
void add_detector(System& s, double z, double aperture_radius = 50.0) {
  s.root.children.push_back({Element{"detector",
                                     ElementKind::Detector,
                                     Pose::along_z(z),
                                     std::nullopt,
                                     {plane("IMG", 0.0, aperture_radius)}}});
  s.paths = {{"main", true, {}}};
}

ToleranceItem item(const std::string& id,
                   const std::string& node,
                   ToleranceKind kind,
                   double limit,
                   ToleranceDistribution distribution = ToleranceDistribution::Uniform) {
  ToleranceItem i;
  i.id = id;
  i.target.anchor = rtt::model::ToleranceAnchor::Node;
  i.target.name = node;
  i.kind = kind;
  i.limit = limit;
  i.distribution = distribution;
  return i;
}

/// A ray coordinate on `surface` of `path` (ADR 0030, point 3) with |Phi - Phi_0| <= max_change.
Criterion ray(const std::string& id,
              RayCoordinate coordinate,
              std::uint16_t field,
              double max_change,
              const std::string& surface = "IMG") {
  rtt::model::RayOperand op;
  op.coordinate = coordinate;
  op.path = "main";
  op.surface = SurfaceId(surface);
  op.field = field;
  op.py = 0.0;
  Criterion c;
  c.id = id;
  c.operand = op;
  c.max_change = max_change;
  return c;
}

/// The rms spot of field 0 about its centroid, 6 hexapolar rings, with `max`.
Criterion spot(double max) {
  rtt::model::SpotRmsOperand op;
  op.path = "main";
  op.field = 0;
  op.rings = 6;
  op.reference = rtt::model::SpotReference::Centroid;
  Criterion c;
  c.id = "spot";
  c.operand = op;
  c.max = max;
  return c;
}

const rtt::material::MaterialLibrary& materials() {
  static const rtt::material::MaterialLibrary library;
  return library;
}

TrialContext context(const System& s) {
  return rtt::tolerance::prepare(s, materials(), nullptr, s.tolerances.items);
}

/// The value of `x`; fails the test case if there is none (clang-tidy, #277).
template <typename T>
const T& required(const std::optional<T>& x) {
  const T* p = x ? &*x : nullptr;
  REQUIRE(p != nullptr);
  return *p;
}

MonteCarloOptions with(std::uint64_t seed, std::uint64_t trials) {
  MonteCarloOptions o;
  o.seed = seed;
  o.trials = trials;
  return o;
}

/// z with erfc(z / sqrt(2)) = alpha, i.e. z_{1 - alpha/2} of the standard normal distribution:
/// P(|X| > z) = 2 (1 - Phi(z)) = erfc(z / sqrt(2)) with Phi(x) = erfc(-x / sqrt(2)) / 2 (NIST
/// DLMF, Eqs. (7.2.1), (7.2.2), substitution t = x / sqrt(2)). Bisection; erfc falls strictly.
double z_of(double alpha) {
  double lo = 0.0;
  double hi = 10.0;
  for (int i = 0; i < 200; ++i) {
    const double mid = (lo + hi) / 2.0;
    (std::erfc(mid / std::numbers::sqrt2) > alpha ? lo : hi) = mid;
  }
  return lo;
}

/// z_{1 - 1e-4 / 2}, see the head of the file; 3.890591886 only as a check.
double z_tests() {
  const volatile double alpha = 1e-4;
  const double z = z_of(alpha);
  REQUIRE(std::abs(z - 3.890591886) <= 1e-8);
  return z;
}

/// The criterion values of criterion k in the Valid trials, in trial order.
std::vector<double> valid_values(const MonteCarloResult& r, std::size_t k) {
  std::vector<double> v;
  for (const TrialRecord& t : r.trials) {
    if (t.result.status != TrialStatus::Valid) continue;
    REQUIRE(k < t.result.criteria.size());
    v.push_back(required(t.result.criteria[k].value));
  }
  return v;
}

/// Checks the statistics of a uniform or truncated normal sample against the bands of the head
/// of the file: true mean 0, true variance `variance`, fourth central moment `mu4`.
void check_bands(const CriterionStatistics& c, double variance, double mu4) {
  const auto n = static_cast<double>(c.count);
  REQUIRE(c.count >= 2);
  const double mean = required(c.mean);
  const double s = required(c.standard_deviation);
  const double var_s2 = (mu4 - (n - 3.0) / (n - 1.0) * variance * variance) / n;
  INFO("mean " << mean << ", s^2 " << s * s << " against " << variance);
  CHECK(std::abs(mean) <= 4.0 * std::sqrt(variance / n));
  CHECK(std::abs(s * s - variance) <= 4.0 * std::sqrt(var_s2));
}

}  // namespace

TEST_CASE("monte carlo: nearest rank with the integer catch", "[tolerance][monte_carlo]") {
  // ADR 0036, point 7: r = ceil(p n), but 0.07 * 100 = 7.000000000000001 in double gives rank 7,
  // not 8.
  using rtt::tolerance::detail::nearest_rank;
  REQUIRE(0.07 * 100.0 > 7.0);  // the trap exists in this build
  CHECK(nearest_rank(0.07, 100) == 7);
  REQUIRE(0.07 * 10000.0 > 700.0);
  CHECK(nearest_rank(0.07, 10000) == 700);
  CHECK(nearest_rank(0.3, 10) == 3);   // 0.3 * 10 is exactly 3 in double
  CHECK(nearest_rank(0.35, 10) == 4);  // 3.5: ceil
  CHECK(nearest_rank(0.5, 200) == 100);
  CHECK(nearest_rank(0.99, 200) == 198);
  CHECK(nearest_rank(1.0, 7) == 7);
  CHECK(nearest_rank(1e-9, 10) == 1);  // ceil(1e-8) = 1
  CHECK(nearest_rank(1e-300, 10) == 1);
  CHECK(nearest_rank(0.5, 1) == 1);
}

TEST_CASE("monte carlo: the Wilson interval", "[tolerance][monte_carlo]") {
  // NIST e-Handbook, Sec. 7.2.4.1 (docs/quellen.md, with z = z_{1 - alpha/2} > 0):
  // [p + z^2/(2n) -+ z sqrt(p (1 - p)/n + z^2/(4n^2))] / (1 + z^2/n).
  using rtt::tolerance::detail::wilson;
  const double z = 1.959963984540054;
  const auto bound = [z](double p, double n, double sign) {
    return (p + z * z / (2.0 * n) +
            sign * z * std::sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n))) /
           (1.0 + z * z / n);
  };
  const rtt::tolerance::Yield half = wilson(5, 10, z);
  CHECK(half.good == 5);
  CHECK(half.total == 10);
  CHECK(required(half.fraction) == 0.5);
  CHECK(std::abs(required(half.lower) - bound(0.5, 10.0, -1.0)) <= 1e-15);
  CHECK(std::abs(required(half.upper) - bound(0.5, 10.0, 1.0)) <= 1e-15);
  CHECK(std::abs(required(half.lower) + required(half.upper) - 1.0) <= 1e-15);  // symmetric
  // By hand: p = 1/2 gives 1/2 -+ z sqrt(1/40 + z^2/400) / (1 + z^2/10); p = 0 gives the upper
  // bound z^2 / (n + z^2), p = 1 the lower bound n / (n + z^2) (40 digits, Python decimal).
  CHECK(std::abs(required(half.lower) - 0.2365930905125640) <= 1e-15);
  CHECK(std::abs(required(half.upper) - 0.7634069094874360) <= 1e-15);
  CHECK(std::abs(required(wilson(0, 10, z).upper) - 0.2775327998628892) <= 1e-15);
  CHECK(std::abs(required(wilson(10, 10, z).lower) - 0.7224672001371108) <= 1e-15);
  // Edges: exactly 0 and 1 (the formula gives about 1e-18 off in double).
  const rtt::tolerance::Yield none = wilson(0, 10, z);
  CHECK(required(none.fraction) == 0.0);
  CHECK(required(none.lower) == 0.0);
  CHECK(std::abs(required(none.upper) - bound(0.0, 10.0, 1.0)) <= 1e-15);
  const rtt::tolerance::Yield all = wilson(10, 10, z);
  CHECK(required(all.fraction) == 1.0);
  CHECK(required(all.upper) == 1.0);
  CHECK(std::abs(required(all.lower) - bound(1.0, 10.0, -1.0)) <= 1e-15);
  // No trials: no values.
  const rtt::tolerance::Yield empty = wilson(0, 0, z);
  CHECK(empty.total == 0);
  CHECK(!empty.fraction.has_value());
  CHECK(!empty.lower.has_value());
  CHECK(!empty.upper.has_value());
}

TEST_CASE("monte carlo: mean and standard deviation in two passes", "[tolerance][monte_carlo]") {
  // ADR 0036, point 7: first the mean, then the squared deviations. For 1e8 + {1, 2, 3} the
  // deviations -1, 0, 1 are exact, so s = 1 exactly; sum x^2 - n mean^2 would lose all digits
  // (x^2 = 1e16: ulp 2).
  using rtt::tolerance::detail::moments;
  const std::vector<double> x = {1e8 + 1.0, 1e8 + 2.0, 1e8 + 3.0};
  const rtt::tolerance::detail::Moments m = moments(x);
  CHECK(required(m.mean) == 1e8 + 2.0);
  CHECK(required(m.standard_deviation) == 1.0);
  const std::vector<double> one = {4.0};
  CHECK(required(moments(one).mean) == 4.0);
  CHECK(!moments(one).standard_deviation.has_value());
  CHECK(!moments(std::vector<double>{}).mean.has_value());
}

TEST_CASE("monte carlo: options, seed and number of trials", "[tolerance][monte_carlo]") {
  // ADR 0035, point 8 and addendum #259 (D6): option before file; no trials -> invalid_argument;
  // the seed falls back to 0 and is in the result.
  System s = lens_system(50.0, 100.0);
  add_detector(s, 100.0);
  s.tolerances.items = {item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1)};
  s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 0, 0.1)};
  const TrialContext no_defaults = context(s);
  CHECK_THROWS_AS(rtt::tolerance::monte_carlo(no_defaults), std::invalid_argument);
  MonteCarloOptions o;
  o.trials = 3;
  const MonteCarloResult zero = rtt::tolerance::monte_carlo(no_defaults, o);
  CHECK(zero.seed == 0);
  CHECK(zero.requested == 3);
  REQUIRE(zero.trials.size() == 3);
  // Trial n carries perturbations(seed, n, items) (ADR 0036, point 5; D2).
  for (std::uint32_t n = 0; n < 3; ++n) {
    INFO("trial " << n);
    CHECK(zero.trials[n].trial == n);
    const std::vector<double> expected = rtt::tolerance::perturbations(0, n, s.tolerances.items);
    REQUIRE(zero.trials[n].values.size() == 1);
    CHECK(same(zero.trials[n].values[0], expected[0]));
    CHECK(expected[0] != 0.0);
  }

  s.tolerances.monte_carlo = rtt::model::MonteCarloDefaults{7, 4};
  const TrialContext defaults = context(s);
  const MonteCarloResult file = rtt::tolerance::monte_carlo(defaults);
  CHECK(file.seed == 7);
  CHECK(file.requested == 4);
  CHECK(file.trials.size() == 4);
  const MonteCarloResult option = rtt::tolerance::monte_carlo(defaults, with(9, 2));
  CHECK(option.seed == 9);
  CHECK(option.trials.size() == 2);
  REQUIRE(!option.trials.empty());
  CHECK(
      same(option.trials[0].values[0], rtt::tolerance::perturbations(9, 0, s.tolerances.items)[0]));
  // Content guard: the file's seed 7 and the option's seed 9 draw differently.
  REQUIRE(!file.trials.empty());
  CHECK(same(file.trials[0].values[0], rtt::tolerance::perturbations(7, 0, s.tolerances.items)[0]));
  CHECK(!same(option.trials[0].values[0], file.trials[0].values[0]));
  // The context gives the file's defaults (TrialContext::monte_carlo_defaults()).
  const auto& from_file = defaults.monte_carlo_defaults();
  CHECK(from_file == s.tolerances.monte_carlo);
  CHECK(!no_defaults.monte_carlo_defaults().has_value());

  const auto rejects = [&](const MonteCarloOptions& options) {
    CHECK_THROWS_AS(rtt::tolerance::monte_carlo(defaults, options), std::invalid_argument);
  };
  MonteCarloOptions bad = with(1, 0);
  rejects(bad);
  bad.trials = (std::uint64_t{1} << 32) + 1;
  rejects(bad);
  for (const double p : {0.0, -0.5, 1.5, std::numeric_limits<double>::quiet_NaN()}) {
    INFO("percentile " << p);
    bad = with(1, 2);
    bad.percentiles = {0.5, p};
    rejects(bad);
  }
  for (const double z : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
    INFO("z " << z);
    bad = with(1, 2);
    bad.z = z;
    rejects(bad);
  }

  // Items that prepare() takes but perturbation() rejects: thrown before any trial (D6).
  std::vector<ToleranceItem> zero_limit = s.tolerances.items;
  zero_limit[0].limit = 0.0;
  std::vector<ToleranceItem> unknown = s.tolerances.items;
  unknown[0].distribution = static_cast<ToleranceDistribution>(99);
  for (const std::vector<ToleranceItem>& items : {zero_limit, unknown}) {
    const TrialContext c = rtt::tolerance::prepare(s, materials(), nullptr, items);
    int reports = 0;
    rtt::trace::RunControl control;
    control.min_interval = std::chrono::milliseconds{0};
    control.progress = [&](const rtt::trace::Progress&) { ++reports; };
    CHECK_THROWS_AS(rtt::tolerance::monte_carlo(c, with(1, 4), &control), std::invalid_argument);
    CHECK(reports == 0);
  }
}

TEST_CASE("monte carlo: statistics of a uniform decenter against the analytic values",
          "[tolerance][monte_carlo]") {
  // Ideal lens f = 50, object at 2 f: m = -1, so ray_y of the chief ray on the detector is
  // Y = (1 - m) X = 2 X for the decenter X (test_sensitivity.cpp, "1 - m"), exact within kEps.
  // X uniform on [-T, T] with T = 0.1 (ADR 0036, point 5): Y uniform on [-0.2, 0.2], mean 0,
  // sigma^2 = 0.2^2 / 3 (NIST 1.3.6.6.2), mu_4 = 9/5 sigma^4 (kurtosis 9/5, ibid.).
  // Criterion |Y - Phi_0| <= L = 0.1 with Phi_0 = 0: true yield L / 0.2 = 0.5.
  System s = lens_system(50.0, 100.0);
  add_detector(s, 100.0);
  s.tolerances.items = {item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1)};
  s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 0, 0.1)};
  const std::uint64_t trials = 200;
  const MonteCarloResult r = rtt::tolerance::monte_carlo(context(s), with(259, trials));
  REQUIRE(r.status == MonteCarloStatus::Completed);
  REQUIRE(r.trials.size() == trials);
  REQUIRE(r.criteria.size() == 1);
  const CriterionStatistics& c = r.criteria[0];
  CHECK(c.id == "chief");
  CHECK(c.count == trials);
  CHECK(c.with_value == trials);
  CHECK(c.without_value == 0);
  CHECK(c.note.empty());
  // Every trial: Y = 2 X.
  for (const TrialRecord& t : r.trials) {
    INFO("trial " << t.trial);
    REQUIRE(t.result.status == TrialStatus::Valid);
    CHECK(std::abs(required(t.result.criteria[0].value) - 2.0 * t.values[0]) <= kEps);
  }
  const double sigma2 = 0.2 * 0.2 / 3.0;
  check_bands(c, sigma2, 9.0 / 5.0 * sigma2 * sigma2);
  CHECK(required(c.min) >= -0.2 - kEps);
  CHECK(required(c.max) <= 0.2 + kEps);
  CHECK(required(c.min) < required(c.max));

  // Percentiles: x_(r) of the sorted values with r = 100, 180, 190, 198 (0.9 * 200 and
  // 0.95 * 200 are exact; 0.5 and 0.99 too).
  std::vector<double> sorted = valid_values(r, 0);
  std::sort(sorted.begin(), sorted.end());
  REQUIRE(c.percentiles.size() == 4);
  const std::size_t ranks[] = {100, 180, 190, 198};
  const double ps[] = {0.5, 0.9, 0.95, 0.99};
  for (std::size_t i = 0; i < 4; ++i) {
    INFO("p " << ps[i]);
    CHECK(c.percentiles[i].p == ps[i]);
    CHECK(c.percentiles[i].rank == ranks[i]);
    CHECK(same(c.percentiles[i].value, sorted[ranks[i] - 1]));
  }
  CHECK(same(required(c.min), sorted.front()));
  CHECK(same(required(c.max), sorted.back()));

  // Yield: one criterion, so the total equals that of the criterion; the true 0.5 lies in the
  // Wilson interval with z_tests() (a separate run with the same seed: only z differs).
  CHECK(c.yield.total == trials);
  CHECK(r.yield.good == c.yield.good);
  std::size_t good = 0;
  for (const TrialRecord& t : r.trials) good += t.good ? 1 : 0;
  CHECK(r.yield.good == good);
  CHECK(good > 0);
  CHECK(good < trials);
  MonteCarloOptions strict = with(259, trials);
  strict.z = z_tests();
  const MonteCarloResult wide = rtt::tolerance::monte_carlo(context(s), strict);
  CHECK(wide.yield.good == good);
  CHECK(required(wide.yield.lower) <= 0.5);
  CHECK(required(wide.yield.upper) >= 0.5);
  // ... and the default 95 % interval is narrower.
  CHECK(required(r.yield.upper) - required(r.yield.lower) <
        required(wide.yield.upper) - required(wide.yield.lower));
}

TEST_CASE("monte carlo: RSS against Monte Carlo, linear with a truncated normal",
          "[tolerance][monte_carlo]") {
  // The decenter of the previous case, now truncated normal: sigma_p = T / 2 = 0.05 cut at
  // +-2 sigma_p (ADR 0036, point 5). Y = 2 X has the parent sigma 0.1 and (docs/quellen.md,
  // "Varianz der bei +-a sigma abgeschnittenen Normalverteilung", a = 2, Z = erf(a / sqrt(2)),
  // phi(a) = exp(-a^2/2) / sqrt(2 pi)):
  //   variance sigma_t^2 = (1 - 2 a phi / Z) 0.1^2 = 0.773741 * 0.01,
  //   mu_4 = (3 Z - 2 a^3 phi - 6 a phi) / Z * 0.1^4 = 1.416189e-4 (both in parent units).
  // sensitivity() predicts sigma_RSS = c |delta| = c * 0.2 (#258), the same sigma_t; the sample
  // variance of Monte Carlo must lie in the band around it.
  System s = lens_system(50.0, 100.0);
  add_detector(s, 100.0);
  s.tolerances.items = {
      item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1, ToleranceDistribution::Normal)};
  s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 0, 0.1)};
  const TrialContext ctx = context(s);
  const volatile double a_in = 2.0;
  const double a = a_in;
  const double phi = std::exp(-a * a / 2.0) / std::sqrt(2.0 * std::numbers::pi);
  const double big_z = std::erf(a / std::numbers::sqrt2);
  const double variance = (1.0 - 2.0 * a * phi / big_z) * 0.01;
  const double mu4 = (3.0 * big_z - 2.0 * a * a * a * phi - 6.0 * a * phi) / big_z * 1e-4;
  REQUIRE(std::abs(variance - 0.773741e-2) <= 1e-8);
  REQUIRE(std::abs(mu4 - 1.416189e-4) <= 1e-10);

  const rtt::tolerance::ToleranceSensitivity rss = rtt::tolerance::sensitivity(ctx);
  REQUIRE(rss.criteria.size() == 1);
  const double predicted = required(rss.criteria[0].predicted_sigma);
  CHECK(std::abs(predicted * predicted - variance) <= 1e-15);

  const MonteCarloResult r = rtt::tolerance::monte_carlo(ctx, with(25901, 200));
  REQUIRE(r.criteria.size() == 1);
  REQUIRE(r.criteria[0].count == 200);
  check_bands(r.criteria[0], predicted * predicted, mu4);
  // Truncated at +-2 sigma_p: |Y| <= 0.2.
  CHECK(required(r.criteria[0].max) <= 0.2 + kEps);
  CHECK(required(r.criteria[0].min) >= -0.2 - kEps);
}

TEST_CASE("monte carlo: RSS against Monte Carlo, quadratic", "[tolerance][monte_carlo]") {
  // test_sensitivity.cpp, "q = k^2 T^2": a detector shift x gives the spot Phi = k |x| with
  // Phi_0 = 0, so Phi^2 = k^2 x^2. For x uniform on [-T, T]: E[Phi^2] = k^2 T^2 / 3, which is the
  // prediction E[Phi^2]_RSS = Phi_0^2 + c^2 q with c^2 = 1/3 and q = k^2 T^2 (#258). The mean of
  // Phi^2 over N trials lies within 4 s(Phi^2) / sqrt(N) of it (NIST 1.3.5.2, with the sample
  // standard deviation of Phi^2 in place of the unknown sigma).
  System s = lens_system(100.0);
  add_detector(s, 100.0);
  s.tolerances.items = {item("det.dz", "detector", ToleranceKind::ShiftZ, 0.5)};
  s.tolerances.criteria = {spot(0.01)};
  const TrialContext ctx = context(s);
  const rtt::tolerance::ToleranceSensitivity rss = rtt::tolerance::sensitivity(ctx);
  REQUIRE(rss.criteria.size() == 1);
  const double predicted = required(rss.criteria[0].predicted_mean_square);
  REQUIRE(predicted > 0.0);

  const MonteCarloResult r = rtt::tolerance::monte_carlo(ctx, with(25902, 200));
  REQUIRE(r.criteria.size() == 1);
  REQUIRE(r.criteria[0].count == 200);
  std::vector<double> squares = valid_values(r, 0);
  for (double& x : squares) x *= x;
  const rtt::tolerance::detail::Moments m = rtt::tolerance::detail::moments(squares);
  const double mean = required(m.mean);
  const double s2 = required(m.standard_deviation);
  INFO("mean of Phi^2 " << mean << " against " << predicted);
  CHECK(std::abs(mean - predicted) <= 4.0 * s2 / std::sqrt(200.0));
  // Content guard: the spot varies.
  CHECK(required(r.criteria[0].max) > 2.0 * required(r.criteria[0].min));
}

namespace {

/// Addendum #259, D7: the warnings of the result against those of the trials.
void check_warnings(const MonteCarloResult& r, std::uint64_t trials) {
  // D7: one warning per code and location, in the order of their first appearance (trial order,
  // then the order within that trial), with the number of trials and the first trial.
  struct Seen {
    std::size_t count = 0;
    std::uint32_t first = 0;
    std::size_t order = 0;
    std::string message;  ///< of the first trial
  };
  using Key = std::pair<std::string, std::string>;
  std::map<Key, Seen> seen;
  for (const TrialRecord& t : r.trials) {
    std::vector<Key> in_trial;
    std::vector<std::string> messages;
    for (const rtt::model::Diagnostic& w : t.result.warnings) {
      const Key key{w.code, w.location};
      if (std::find(in_trial.begin(), in_trial.end(), key) == in_trial.end()) {
        in_trial.push_back(key);
        messages.push_back(w.message);
      }
    }
    for (std::size_t j = 0; j < in_trial.size(); ++j) {
      const std::size_t order = seen.size();
      ++seen.try_emplace(in_trial[j], Seen{0, t.trial, order, messages[j]}).first->second.count;
    }
  }
  REQUIRE(!seen.empty());  // content guard: the trials report warnings
  REQUIRE(r.warnings.size() == seen.size());
  for (std::size_t i = 0; i < r.warnings.size(); ++i) {
    const rtt::model::Diagnostic& w = r.warnings[i];
    INFO(w.code << " at " << w.location << ": " << w.message);
    const auto it = seen.find({w.code, w.location});
    REQUIRE(it != seen.end());
    CHECK(it->second.order == i);
    const std::string counted = "in " + std::to_string(it->second.count) + " of " +
                                std::to_string(trials) + " trials, first trial " +
                                std::to_string(it->second.first);
    CHECK(w.message == it->second.message + " (" + counted + ")");
  }
}

}  // namespace

TEST_CASE("monte carlo: invalid trials and keep_trials", "[tolerance][monte_carlo]") {
  // The 2 f system with a detector of radius 0.15 mm: the chief ray lands at Y = 2 X with X
  // uniform on [-0.1, 0.1], so it misses the detector for |X| > 0.075, in about a quarter of the
  // trials. ray_y on IMG is then not defined and the trial is Invalid (#257); ray_y on the stop
  // (in front of the lens) keeps its value. Addendum #259, D5: the statistics of both criteria
  // count only the Valid trials, their yields all N trials; keep_trials off keeps the trials that
  // are not Valid.
  System s = lens_system(50.0, 100.0);
  add_detector(s, 100.0, 0.15);
  s.tolerances.items = {item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1)};
  s.tolerances.criteria = {ray("image", RayCoordinate::Y, 0, 0.1),
                           ray("stop", RayCoordinate::Y, 0, 1.0, "STO")};
  const TrialContext ctx = context(s);
  const std::uint64_t trials = 100;
  const MonteCarloResult r = rtt::tolerance::monte_carlo(ctx, with(25903, trials));
  REQUIRE(r.status == MonteCarloStatus::Completed);
  REQUIRE(r.trials.size() == trials);
  std::size_t valid = 0;
  std::size_t stop_good = 0;
  std::size_t good = 0;
  for (const TrialRecord& t : r.trials) {
    INFO("trial " << t.trial);
    REQUIRE(t.result.criteria.size() == 2);
    const bool lost = std::abs(t.values[0]) > 0.075 + 1e-9;
    const bool inside = std::abs(t.values[0]) < 0.075 - 1e-9;
    if (lost) CHECK(t.result.status == TrialStatus::Invalid);
    if (inside) CHECK(t.result.status == TrialStatus::Valid);
    if (t.result.status == TrialStatus::Valid) {
      ++valid;
    } else {
      CHECK(!t.result.reason.empty());
      CHECK(!t.result.criteria[0].value.has_value());
      CHECK(!t.good);
    }
    CHECK(t.result.criteria[1].value.has_value());  // the stop is met in every trial
    stop_good += t.result.criteria[1].satisfied.value_or(false) ? 1 : 0;
    good += t.good ? 1 : 0;
    CHECK(t.good == (t.result.status == TrialStatus::Valid &&
                     t.result.criteria[0].satisfied.value_or(false) &&
                     t.result.criteria[1].satisfied.value_or(false)));
  }
  REQUIRE(valid > 0);
  REQUIRE(valid < trials);  // content guard: there are Invalid trials
  REQUIRE(r.criteria.size() == 2);
  const CriterionStatistics& image = r.criteria[0];
  const CriterionStatistics& stop = r.criteria[1];
  CHECK(image.with_value == valid);
  CHECK(image.without_value == trials - valid);
  CHECK(stop.with_value == trials);
  CHECK(stop.without_value == 0);
  CHECK(image.count == valid);
  CHECK(stop.count == valid);  // the stop values of Invalid trials stay out
  const rtt::tolerance::detail::Moments m = rtt::tolerance::detail::moments(valid_values(r, 1));
  CHECK(same(required(stop.mean), required(m.mean)));
  CHECK(stop.yield.total == trials);
  CHECK(stop.yield.good == stop_good);  // counted over all N trials
  CHECK(stop_good == trials);
  CHECK(image.yield.total == trials);
  CHECK(image.yield.good <= valid);
  CHECK(r.yield.good == good);
  CHECK(r.yield.total == trials);

  // D7 (see the next case): a chief ray that misses the detector makes the trial Invalid with
  // a reason (#257); it is not a warning.
  CHECK(r.warnings.empty());

  // keep_trials off: only the trials that are not Valid, the same values; statistics unchanged.
  MonteCarloOptions lean = with(25903, trials);
  lean.keep_trials = false;
  const MonteCarloResult l = rtt::tolerance::monte_carlo(ctx, lean);
  REQUIRE(l.trials.size() == trials - valid);
  std::size_t j = 0;
  for (const TrialRecord& t : r.trials) {
    if (t.result.status == TrialStatus::Valid) continue;
    REQUIRE(j < l.trials.size());
    CHECK(l.trials[j].trial == t.trial);
    CHECK(same(l.trials[j].values[0], t.values[0]));
    CHECK(l.trials[j].result.reason == t.result.reason);
    ++j;
  }
  REQUIRE(l.criteria.size() == 2);
  CHECK(same(required(l.criteria[1].mean), required(stop.mean)));
  CHECK(l.yield.good == r.yield.good);
  CHECK(l.warnings.size() == r.warnings.size());
}

TEST_CASE("monte carlo: the warnings of the trials", "[tolerance][monte_carlo]") {
  // Addendum #259, D7. The detector (radius 1.05 mm) stands 10 mm behind the focus of f = 50 mm
  // with EPD 10 mm: the axial beam covers a disc of radius 5 * 10 / 50 = 1 mm. A decentre of the
  // detector uniform on [-1, 1] keeps the chief ray on it (it would miss for |dy| > 1.05, and
  // the spot is then not defined), so every trial is Valid. The spot warns "rays.lost" above half
  // of the rays lost (SpotOptions::lost_warning_fraction). The part of the disc outside the
  // detector (area of two overlapping circles) is 0.40 at |dy| = 0.7 and 0.55 at |dy| = 0.95, one
  // half at about 0.87. Of the 127 rays of 6 hexapolar rings, 0.41 are lost at |dy| = 0.7 and
  // 0.55 at |dy| = 0.95; their fraction crosses one half at about 0.84. Some trials warn, some
  // not.
  System s = lens_system(50.0);
  add_detector(s, 60.0, 1.05);
  s.tolerances.items = {item("det.dy", "detector", ToleranceKind::DecenterY, 1.0)};
  s.tolerances.criteria = {spot(1.0)};
  const TrialContext ctx = context(s);
  REQUIRE(ctx.nominal().status == TrialStatus::Valid);
  REQUIRE(ctx.nominal().warnings.empty());  // the nominal beam lies inside the detector
  const std::uint64_t trials = 50;
  const MonteCarloResult r = rtt::tolerance::monte_carlo(ctx, with(25904, trials));
  REQUIRE(r.status == MonteCarloStatus::Completed);
  REQUIRE(r.trials.size() == trials);
  std::size_t clipped = 0;
  for (const TrialRecord& t : r.trials) {
    INFO("trial " << t.trial << ", dy " << t.values[0]);
    CHECK(t.result.status == TrialStatus::Valid);
    CHECK(t.result.criteria[0].value.has_value());
    if (std::abs(t.values[0]) < 0.7) CHECK(t.result.warnings.empty());
    if (std::abs(t.values[0]) > 0.95) CHECK(!t.result.warnings.empty());
    clipped += t.result.warnings.empty() ? 0 : 1;
  }
  REQUIRE(clipped > 0);       // content guard: some trials warn
  REQUIRE(clipped < trials);  // and some do not, so the count in the message is not N
  check_warnings(r, trials);
}

TEST_CASE("monte carlo: a criterion without valid values", "[tolerance][monte_carlo]") {
  // A detector of radius 1e-3 mm: the chief ray of field 1 (0.5 degree on a 100 mm image
  // distance, Y about 0.87 mm) misses it in every trial; n = 0 gives no values and the note.
  System s = lens_system(100.0);
  add_detector(s, 100.0, 1e-3);
  s.tolerances.items = {item("lens.dy", "lens", ToleranceKind::DecenterY, 0.01)};
  s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 1, 0.1)};
  const MonteCarloResult r = rtt::tolerance::monte_carlo(context(s), with(1, 5));
  REQUIRE(r.criteria.size() == 1);
  const CriterionStatistics& c = r.criteria[0];
  CHECK(c.count == 0);
  CHECK(c.without_value == 5);
  CHECK(c.note == "no valid values");
  CHECK(!c.mean.has_value());
  CHECK(!c.standard_deviation.has_value());
  CHECK(!c.min.has_value());
  CHECK(!c.max.has_value());
  CHECK(c.percentiles.empty());
  CHECK(c.yield.good == 0);
  CHECK(required(c.yield.lower) == 0.0);
  CHECK(r.trials.size() == 5);
}

TEST_CASE("monte carlo: cancelled runs keep the finished trials", "[tolerance][monte_carlo]") {
  // ADR 0036, points 4 and 7: the finished trials in ascending order, N their number, status
  // Cancelled.
  System s = lens_system(50.0, 100.0);
  add_detector(s, 100.0);
  s.tolerances.items = {item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1)};
  s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 0, 0.1)};
  const TrialContext ctx = context(s);
  const auto check_nothing = [](const MonteCarloResult& r) {
    CHECK(r.status == MonteCarloStatus::Cancelled);
    CHECK(r.trials.empty());
    CHECK(r.yield.total == 0);
    CHECK(!r.yield.fraction.has_value());
    REQUIRE(r.criteria.size() == 1);
    CHECK(r.criteria[0].count == 0);
    CHECK(r.criteria[0].note == "no valid values");
  };
  SECTION("before the first trial") {
    rtt::trace::RunControl control;
    control.cancel = rtt::trace::CancelToken{};
    control.cancel->request_cancel();
    check_nothing(rtt::tolerance::monte_carlo(ctx, with(1, 50), &control));
  }
  SECTION("already in prepare (nominal() cancelled, #257)") {
    rtt::trace::RunControl control;
    control.cancel = rtt::trace::CancelToken{};
    control.cancel->request_cancel();
    const TrialContext cancelled =
        rtt::tolerance::prepare(s, materials(), nullptr, s.tolerances.items, {}, &control);
    REQUIRE(cancelled.nominal().status == TrialStatus::Cancelled);
    check_nothing(rtt::tolerance::monte_carlo(cancelled, with(1, 50)));
  }
  SECTION("after the first progress report") {
    rtt::trace::RunControl control;
    control.cancel = rtt::trace::CancelToken{};
    control.min_interval = std::chrono::milliseconds{0};
    std::atomic<int> reports{0};
    control.progress = [&](const rtt::trace::Progress& p) {
      CHECK(p.stage == "monte_carlo");
      if (++reports == 1) control.cancel->request_cancel();
    };
    const std::uint64_t trials = 200;
    const MonteCarloResult r = rtt::tolerance::monte_carlo(ctx, with(1, trials), &control);
    CHECK(reports >= 1);
    CHECK(r.status == MonteCarloStatus::Cancelled);
    REQUIRE(!r.trials.empty());
    CHECK(r.trials.size() < trials);
    CHECK(r.yield.total == r.trials.size());
    REQUIRE(r.criteria.size() == 1);
    CHECK(r.criteria[0].count == r.trials.size());
    for (std::size_t i = 1; i < r.trials.size(); ++i) {
      CHECK(r.trials[i - 1].trial < r.trials[i].trial);
    }
    // Each finished trial is the one of an uncancelled run.
    const MonteCarloResult full = rtt::tolerance::monte_carlo(ctx, with(1, trials));
    REQUIRE(full.trials.size() == trials);
    for (const TrialRecord& t : r.trials) {
      CHECK(same(required(t.result.criteria[0].value),
                 required(full.trials[t.trial].result.criteria[0].value)));
    }
  }
}

TEST_CASE("monte carlo: the overload with a System expands and keeps the expansion",
          "[tolerance][monte_carlo]") {
  System s = lens_system(50.0, 100.0);
  add_detector(s, 100.0);
  s.tolerances.items = {item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1)};
  s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 0, 0.1)};
  s.tolerances.monte_carlo = rtt::model::MonteCarloDefaults{3, 4};
  const MonteCarloResult r = rtt::tolerance::monte_carlo(s, materials(), nullptr);
  const rtt::tolerance::ToleranceItems& expansion = required(r.expansion);
  CHECK(expansion.items == r.items);
  REQUIRE(r.items.size() == 1);
  CHECK(r.items[0].id == "lens.dy");
  CHECK(r.trials.size() == 4);
  // The same as the overload with a context.
  const MonteCarloResult direct = rtt::tolerance::monte_carlo(context(s));
  CHECK(!direct.expansion.has_value());
  REQUIRE(direct.criteria.size() == 1);
  CHECK(same(required(direct.criteria[0].mean), required(r.criteria[0].mean)));
}

TEST_CASE("monte carlo: bitwise the same for 1, 4 and all threads", "[tolerance][monte_carlo]") {
  // ADR 0036, point 8: every number of threads and every repetition; another seed changes the
  // values.
  System s = lens_system(50.0, 100.0);
  add_detector(s, 100.0);
  s.tolerances.items = {
      item("lens.dy", "lens", ToleranceKind::DecenterY, 0.1),
      item("lens.dz", "lens", ToleranceKind::ShiftZ, 0.2, ToleranceDistribution::Normal),
      item("det.dz", "detector", ToleranceKind::ShiftZ, 0.1, ToleranceDistribution::Endpoint)};
  rtt::model::SpotRmsOperand op;
  op.path = "main";
  op.field = 1;
  Criterion quadratic;
  quadratic.id = "spot";
  quadratic.operand = op;
  quadratic.max = 0.05;
  s.tolerances.criteria = {ray("chief", RayCoordinate::Y, 1, 0.03), quadratic};
  const TrialContext ctx = context(s);
  const auto run = [&](int threads, std::uint64_t seed) {
    oneapi::tbb::task_arena arena(threads);
    return arena.execute([&] { return rtt::tolerance::monte_carlo(ctx, with(seed, 40)); });
  };
  const auto values = [](const MonteCarloResult& r) {
    std::vector<double> v;
    const auto add = [&](const std::optional<double>& x) { v.push_back(x.value_or(-1e300)); };
    for (const TrialRecord& t : r.trials) {
      for (const double x : t.values) v.push_back(x);
      for (const auto& c : t.result.criteria) add(c.value);
      v.push_back(t.good ? 1.0 : 0.0);
    }
    for (const CriterionStatistics& c : r.criteria) {
      for (const auto* x : {&c.mean, &c.standard_deviation, &c.min, &c.max}) add(*x);
      for (const auto& p : c.percentiles) v.push_back(p.value);
      add(c.yield.fraction);
      add(c.yield.lower);
      add(c.yield.upper);
    }
    add(r.yield.fraction);
    add(r.yield.lower);
    add(r.yield.upper);
    return v;
  };
  const MonteCarloResult first = run(1, 11);
  // Content guard: a complete run with varying values.
  REQUIRE(first.status == MonteCarloStatus::Completed);
  REQUIRE(first.criteria.size() == 2);
  CHECK(first.criteria[0].count > 2);
  CHECK(required(first.criteria[0].standard_deviation) > 0.0);
  CHECK(required(first.criteria[1].standard_deviation) > 0.0);
  const std::vector<double> one = values(first);
  REQUIRE(one.size() > 200);
  for (const int threads : {1, 4, oneapi::tbb::task_arena::automatic}) {
    INFO(threads << " threads");
    const std::vector<double> other = values(run(threads, 11));
    REQUIRE(other.size() == one.size());
    bool all_same = true;
    for (std::size_t i = 0; i < one.size(); ++i) all_same = all_same && same(one[i], other[i]);
    CHECK(all_same);
  }
  const std::vector<double> reseeded = values(run(4, 12));
  REQUIRE(reseeded.size() == one.size());
  bool any_differs = false;
  for (std::size_t i = 0; i < one.size(); ++i)
    any_differs = any_differs || !same(one[i], reseeded[i]);
  CHECK(any_differs);
}
