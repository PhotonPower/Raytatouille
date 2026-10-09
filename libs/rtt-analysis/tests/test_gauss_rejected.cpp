// Analyses that average over their rays without quadrature weights reject the Gaussian pupil
// sampling (#168, coordinator's decision (a)): an unweighted RMS or mean over Gauss points is not
// the mean over the pupil. They take the other samplings unchanged.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <stdexcept>
#include <string>

#include "rtt/analysis/ghosts.hpp"
#include "rtt/analysis/paths.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/coating/catalog.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/compile/ghosts.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/trace/sources.hpp"

using Catch::Matchers::ContainsSubstring;
using rtt::compile::CompiledSystem;
using rtt::compile::PathId;
using rtt::trace::GaussPupil;

namespace {

rtt::model::System load(const char* file) {
  return rtt::io::load_system(std::string(RTT_REFERENCE_DIR) + "/" + file);
}

/// The message of the std::invalid_argument of `f`, or "<no exception>".
template <class F>
std::string rejection(F f) {
  try {
    f();
  } catch (const std::invalid_argument& e) {
    return e.what();
  }
  return "<no exception>";
}

}  // namespace

TEST_CASE("spot rejects the gauss sampling and keeps the others", "[gauss][spot]") {
  const rtt::material::MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(load("m1/singlet_const.rtt.json"), lib);
  rtt::analysis::SpotOptions options;
  options.sampling = GaussPupil{};
  CHECK_THAT(rejection([&] { (void)rtt::analysis::spot(cs, PathId{0}, 0, 0, options); }),
             ContainsSubstring("gauss"));
  options.sampling = rtt::trace::HexapolarPupil{2};
  CHECK(rtt::analysis::spot(cs, PathId{0}, 0, 0, options).rays_launched == 19);
}

TEST_CASE("the path analyses reject the gauss sampling", "[gauss][paths]") {
  const rtt::material::MaterialLibrary lib;
  const CompiledSystem cs = rtt::compile::compile(load("m1/singlet_const.rtt.json"), lib);
  rtt::analysis::PathOptions options;
  options.sampling = GaussPupil{2, 6};
  CHECK_THAT(
      rejection([&] { (void)rtt::analysis::path_transmission(cs, PathId{0}, 0, 0, options); }),
      ContainsSubstring("gauss"));
  CHECK_THAT(rejection([&] {
               (void)rtt::analysis::opl_difference(cs, PathId{0}, PathId{0}, 0, 0, options);
             }),
             ContainsSubstring("gauss"));
  options.sampling = rtt::trace::HexapolarPupil{1};
  CHECK(rtt::analysis::path_transmission(cs, PathId{0}, 0, 0, options).rays.size() == 7);
}

TEST_CASE("the ghost ranking rejects the gauss sampling", "[gauss][ghosts]") {
  const rtt::material::MaterialLibrary lib;
  const rtt::coating::CoatingLibrary coatings;
  const rtt::compile::GhostSystem g =
      rtt::compile::compile_with_ghosts(load("m4/ghost_plates.rtt.json"), "main", lib, coatings);
  rtt::analysis::GhostRankingOptions options;
  options.sampling = GaussPupil{};
  CHECK_THAT(rejection([&] { (void)rtt::analysis::ghost_ranking(g, 0, 0, options); }),
             ContainsSubstring("gauss"));
}
