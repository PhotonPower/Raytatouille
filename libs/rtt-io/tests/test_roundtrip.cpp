#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>

#include "rtt/io/json_io.hpp"
#include "rtt/model/validate.hpp"
#include "test_support.hpp"

namespace fs = std::filesystem;

namespace {

std::vector<fs::path> reference_files() {
  std::vector<fs::path> files;
  for (const auto& entry : fs::recursive_directory_iterator(RTT_REFERENCE_DIR)) {
    const std::string name = entry.path().filename().string();
    if (entry.is_regular_file() && name.size() > 9 && name.ends_with(".rtt.json")) {
      files.push_back(entry.path());
    }
  }
  std::sort(files.begin(), files.end());
  return files;
}

std::string read_text(const fs::path& p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

}  // namespace

TEST_CASE("reference files exist", "[io]") {
  REQUIRE(reference_files().size() >= 3U);
}

TEST_CASE("file -> model -> file is byte-identical for canonical files", "[io][roundtrip]") {
  for (const fs::path& file : reference_files()) {
    INFO(file.string());
    const std::string original = read_text(file);
    const rtt::model::System system = rtt::io::parse_system(original);
    REQUIRE(rtt::io::to_json(system) == original);
  }
}

TEST_CASE("reference files are semantically valid", "[io]") {
  for (const fs::path& file : reference_files()) {
    INFO(file.string());
    const auto diagnostics = rtt::model::validate(rtt::io::load_system(file));
    for (const auto& d : diagnostics) UNSCOPED_INFO(rtt::model::to_string(d));
    REQUIRE_FALSE(rtt::model::has_errors(diagnostics));
  }
}

TEST_CASE("model -> file -> model preserves the model", "[io][roundtrip]") {
  rtt::model::System s = rtt::model::test::make_singlet();
  // Exercise non-default parameter features.
  auto& lens = rtt::model::test::element(s, 1);
  lens.pose.position[2].variable = true;
  lens.surfaces[1].pose.position[2].pickup = "L1.S1.thickness";
  lens.surfaces[0].shape.terms.push_back(
      rtt::model::ZernikeSag{rtt::model::Param(12.7), {0.0, 1e-4, rtt::model::Param(-2e-5)}});
  s.object = {false, rtt::model::Param(250.0)};
  s.environment.temperature_c = 22.5;

  const std::string text = rtt::io::to_json(s);
  const rtt::model::System back = rtt::io::parse_system(text);
  REQUIRE(back == s);
  REQUIRE(rtt::io::to_json(back) == text);
}

TEST_CASE("save and load through a file", "[io]") {
  const rtt::model::System s = rtt::model::test::make_singlet();
  const fs::path tmp = fs::temp_directory_path() / "rtt_io_roundtrip_test.rtt.json";
  rtt::io::save_system(s, tmp);
  REQUIRE(rtt::io::load_system(tmp) == s);
  fs::remove(tmp);
}

TEST_CASE("non-finite numbers are rejected on write", "[io]") {
  rtt::model::System s = rtt::model::test::make_singlet();
  s.wavelengths[0].um = std::numeric_limits<double>::quiet_NaN();
  REQUIRE_THROWS_AS(rtt::io::to_json(s), std::invalid_argument);
}
