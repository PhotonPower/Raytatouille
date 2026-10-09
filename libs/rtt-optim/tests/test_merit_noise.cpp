// Noise of a merit evaluation (ADR 0030, point 9; #167): the measurement behind kMeritPrecision.
// Hidden ([.noise]), not run in CI; run it by hand with
//   rtt_optim_tests "[.noise]" -s
// For every reference system and path that the merit function accepts, all unbound radius and
// z-position Params are variables (plus the variable rows and Params already in the file), and
// the residuals are the quantities the operands use: EFL, and per field the real ray x and y on
// the image surface for a hexapolar pupil with 3 rings, the RMS spot (centroid) and the RMS OPD.
// Per variable j the residuals are evaluated at theta + i delta e_j, i = -4..4, with
// delta = 1e-6 max(|theta_j|, 1); the noise of residual k is sigma_k = |Delta^8 f_k| / sqrt(12870)
// (the eighth difference of independent noise of variance sigma^2 has variance
// sum_i C(8, i)^2 sigma^2 = C(16, 8) sigma^2 = 12870 sigma^2; the smooth part, f^(8) delta^8, is
// far below). eps_f of a case = max_{k, j} sigma_k / ||f(theta_0)||_inf. A residual that is not
// defined at one of the 9 points drops out for that variable and is counted.

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/model/optimization.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/optim/merit.hpp"
#include "rtt/optim/optimize.hpp"
#include "rtt/trace/sources.hpp"

namespace fs = std::filesystem;
using rtt::model::System;

namespace {

/// Merit function of the measurement for `path` of `s` (see the file comment).
void add_operands(System& s, const rtt::compile::CompiledSystem& cs, const std::string& path) {
  const rtt::compile::CompiledPath& p = cs.path(*cs.find_path(path));
  const std::uint32_t image = p.events.back().surface;
  std::uint32_t occurrence = 0;
  for (std::size_t e = 0; e + 1 < p.events.size(); ++e) {
    occurrence += p.events[e].surface == image ? 1U : 0U;
  }
  rtt::model::Optimization& merit = s.optimization;
  merit = {};
  rtt::model::FirstOrderOperand efl;
  efl.path = path;
  merit.operands.emplace_back(efl);
  const std::vector<rtt::trace::PupilPoint> pupil =
      rtt::trace::pupil_points(rtt::trace::HexapolarPupil{3});
  for (std::size_t f = 0; f < s.fields.points.size(); ++f) {
    const auto field = static_cast<std::uint16_t>(f);
    for (const rtt::trace::PupilPoint& q : pupil) {
      for (const auto c : {rtt::model::RayCoordinate::X, rtt::model::RayCoordinate::Y}) {
        rtt::model::RayOperand ray;
        ray.coordinate = c;
        ray.path = path;
        ray.surface = cs.surfaces()[image].id;
        if (occurrence > 0) ray.occurrence = occurrence;
        ray.field = field;
        ray.px = q.px;
        ray.py = q.py;
        merit.operands.emplace_back(ray);
      }
    }
    rtt::model::SpotRmsOperand spot;
    spot.path = path;
    spot.field = field;
    spot.rings = 3;
    merit.operands.emplace_back(spot);
    rtt::model::OpdRmsOperand opd;
    opd.path = path;
    opd.field = field;
    opd.grid = 17;
    merit.operands.emplace_back(opd);
  }
}

/// Makes every unbound radius and z position a variable (bounds are kept, the measurement
/// evaluates the external values directly).
void add_variables(System& s) {
  rtt::model::for_each_param(s, [](std::string_view pointer, rtt::model::Param& p) {
    if (!p.is_bound() && (pointer.ends_with("/radius") || pointer.ends_with("/position/2"))) {
      p.variable = true;
    }
  });
}

struct Case {
  std::string name;
  std::size_t variables = 0;
  std::size_t residuals = 0;
  double eps = 0.0;  ///< max sigma / ||f0||_inf
  std::size_t worst_residual = 0;
  std::size_t worst_variable = 0;
  std::size_t dropped = 0;  ///< (residual, variable) pairs without 9 valid values
};

constexpr std::array<double, 9> kBinomial8 = {1, -8, 28, -56, 70, -56, 28, -8, 1};

Case measure(const std::string& name, const System& s, const rtt::material::MaterialLibrary& lib) {
  const rtt::optim::MeritFunction merit(s, lib, nullptr);
  const std::vector<double> p0 = merit.start();
  const rtt::optim::MeritEvaluation start = merit.evaluate(p0);
  Case c{name, p0.size(), start.residuals.size()};
  double norm = 0.0;
  for (const double r : start.residuals) {
    if (std::isfinite(r)) norm = std::max(norm, std::abs(r));
  }
  REQUIRE(norm > 0.0);
  double worst = 0.0;
  for (std::size_t j = 0; j < p0.size(); ++j) {
    const double delta = 1e-6 * std::max(std::abs(p0[j]), 1.0);
    std::vector<std::vector<double>> f;
    for (int i = -4; i <= 4; ++i) {
      std::vector<double> p = p0;
      p[j] += i * delta;
      try {
        f.push_back(merit.evaluate(p).residuals);
      } catch (const std::exception&) {
        f.emplace_back(start.residuals.size(), std::nan(""));
      }
    }
    for (std::size_t k = 0; k < start.residuals.size(); ++k) {
      double d8 = 0.0;
      bool valid = std::isfinite(start.residuals[k]);
      for (std::size_t i = 0; i < 9; ++i) {
        valid = valid && std::isfinite(f[i][k]);
        d8 += kBinomial8[i] * f[i][k];
      }
      if (!valid) {
        ++c.dropped;
        continue;
      }
      const double sigma = std::abs(d8) / std::sqrt(12870.0);
      if (sigma > worst) {
        worst = sigma;
        c.worst_residual = k;
        c.worst_variable = j;
      }
    }
  }
  c.eps = worst / norm;
  return c;
}

}  // namespace

TEST_CASE("merit noise: relative noise of a merit evaluation (manual, #167)", "[.noise]") {
  const fs::path root(RTT_REFERENCE_DIR);
  rtt::material::MaterialLibrary schott;
  schott.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  rtt::material::MaterialLibrary schott_m2;
  schott_m2.add_catalog(std::string(RTT_CATALOG_DIR) + "/m2/schott.agf");
  std::vector<fs::path> files;
  for (const auto& dir : fs::directory_iterator(root)) {
    if (!dir.is_directory() || !dir.path().filename().string().starts_with("m")) continue;
    for (const auto& e : fs::directory_iterator(dir.path())) {
      if (e.path().string().ends_with(".rtt.json")) files.push_back(e.path());
    }
  }
  std::sort(files.begin(), files.end());
  double overall = 0.0;
  std::cout << std::setprecision(3);
  for (const fs::path& file : files) {
    const std::string rel = fs::relative(file, root).generic_string();
    System base;
    try {
      base = rtt::io::load_system(file);
    } catch (const std::exception& e) {
      std::cout << "SKIP " << rel << ": load: " << e.what() << "\n";
      continue;
    }
    const rtt::material::MaterialLibrary& lib = rel.starts_with("m2/") ? schott_m2 : schott;
    std::optional<rtt::compile::CompiledSystem> cs;
    try {
      cs.emplace(rtt::compile::compile(base, lib));
    } catch (const std::exception& e) {
      std::cout << "SKIP " << rel << ": compile: " << e.what() << "\n";
      continue;
    }
    for (const auto& path : cs->paths()) {
      const std::string name = rel + " [" + path.name + "]";
      System s = base;
      add_variables(s);
      add_operands(s, *cs, path.name);
      try {
        const Case c = measure(name, s, lib);
        std::cout << "CASE " << c.name << ": n = " << c.variables << ", m = " << c.residuals
                  << ", eps_f = " << c.eps << " (residual " << c.worst_residual << ", variable "
                  << c.worst_variable << "), dropped " << c.dropped << "\n";
        overall = std::max(overall, c.eps);
      } catch (const std::exception& e) {
        std::cout << "SKIP " << name << ": " << e.what() << "\n";
      }
    }
  }
  const double proposal = std::pow(10.0, std::ceil(std::log10(10.0 * overall)));
  std::cout << "MAX eps_f = " << overall << ", kMeritPrecision proposal = " << proposal
            << " (in use: " << rtt::optim::kMeritPrecision << ")\n";
  CHECK(overall > 0.0);
}
