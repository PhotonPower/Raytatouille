// The C++ result writer (ADR 0052, point 1): floats as Python's repr and strings as json.dumps
// with ensure_ascii=False. The expected texts were printed by Python 3 (repr(float.fromhex(h))
// and json.dumps(..., ensure_ascii=False)); test_result_json.py compares whole results.

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "rtt/io/result_json.hpp"
#include "rtt/io/result_names.hpp"

TEST_CASE("floats are written as Python's repr (ADR 0052, point 1)", "[io][result_json]") {
  const std::pair<double, const char*> cases[] = {
      {0x0.0p+0, "0.0"},
      {-0x0.0p+0, "-0.0"},
      {0x1.0p+0, "1.0"},
      {0x1.9p+6, "100.0"},
      {0x1.999999999999ap-4, "0.1"},
      {0x1.4f8b588e368f1p-17, "1e-05"},
      {0x1.a36e2eb1c432dp-14, "0.0001"},
      {0x1.a36e2eb1c432cp-14, "9.999999999999999e-05"},
      {0x1.1c37937e08p+53, "1e+16"},
      {0x1.1c37937e07fffp+53, "9999999999999998.0"},
      {0x1.b69b4ba630f35p+56, "1.2345678901234568e+17"},
      {0x0.0000000000001p-1022, "5e-324"},
      {0x1.fffffffffffffp+1023, "1.7976931348623157e+308"},
      {0x1.0p-1022, "2.2250738585072014e-308"},
      {0x1.edd2f1a9fbe77p+6, "123.456"},
      {-0x1.421f5f40d8376p-23, "-1.5e-07"},
      {0x1.0f0cf064dd592p+73, "1e+22"},
      {0x1.56e1fc2f8f359p-997, "1e-300"},
      {0x1.3333333333334p-2, "0.30000000000000004"},
      {0x1.550f7dca7p+51, "3000000000000000.0"},
      {0x1.5ee2a2eb5a5c4p+53, "1.2345678901234568e+16"},
  };
  for (const auto& [value, expected] : cases) {
    INFO(expected);
    CHECK(rtt::io::python_float_repr(value) == expected);
  }
}

TEST_CASE("strings, NaN and the envelope as json.dumps (ADR 0052, point 1)", "[io][result_json]") {
  rtt::optim::OptimResult r;
  r.status = rtt::optim::LmStatus::Cancelled;
  // Quote, backslash, newline, U+0001, DEL, e acute and U+2028 in UTF-8.
  r.patch = "a\"b\\c\n\x01\x7f\xc3\xa9\xe2\x80\xa8";
  rtt::optim::OptimIteration it;
  it.k = 1;
  it.rho = std::numeric_limits<double>::quiet_NaN();
  it.step_norm = std::numeric_limits<double>::infinity();
  r.history = {it};
  const std::string text = rtt::io::to_result_json(r);
  CHECK(text.starts_with(
      R"({"format": "raytatouille-result", "schema_version": "0.1.23", "type": "OptimResult", )"
      // An ordinary literal: MSVC reads \u in a raw string literal as a universal character name.
      "\"data\": {\"status\": \"CANCELLED\", \"patch\": \"a\\\"b\\\\c\\n\\u0001"
      "\x7f\xc3\xa9\xe2\x80\xa8\""));
  CHECK(text.find(R"("rho": {"dtype": "float64", "shape": [1], "values": ["NaN"]})") !=
        std::string::npos);
  CHECK(text.find(R"("step_norm": {"dtype": "float64", "shape": [1], "values": ["Infinity"]})") !=
        std::string::npos);
  CHECK(text.ends_with(R"("failed_evaluations": 0}})"));  // no final newline
}

TEST_CASE("every enum member of the format has one distinct name (ADR 0052, point 1)",
          "[io][result_json]") {
  const auto check = [](const std::vector<std::string_view>& names) {
    for (const std::string_view n : names) CHECK(n != "?");
    const std::set<std::string_view> distinct(names.begin(), names.end());
    CHECK(distinct.size() == names.size());
  };
  using rtt::optim::LmStatus;
  check({rtt::io::result_name(LmStatus::ConvergedGradient),
         rtt::io::result_name(LmStatus::ConvergedStep),
         rtt::io::result_name(LmStatus::ConvergedMerit),
         rtt::io::result_name(LmStatus::MaxIterations), rtt::io::result_name(LmStatus::Cancelled),
         rtt::io::result_name(LmStatus::Failed)});
  using rtt::diagnostics::Severity;
  check({rtt::io::result_name(Severity::Error), rtt::io::result_name(Severity::Warning)});
  std::vector<std::string_view> status;
  for (std::size_t s = 0; s < rtt::trace::kRayStatusCount; ++s) {
    status.push_back(rtt::io::result_name(static_cast<rtt::trace::RayStatus>(s)));
  }
  check(status);
}

TEST_CASE("NaN as a scalar and -infinity in an array as the format writes them",
          "[io][result_json]") {
  // json.dumps of to_dict: a scalar NaN as {"float": "NaN"}, an array element -inf as
  // "-Infinity" (ADR 0023, point 3).
  rtt::analysis::ZernikeFit f;
  f.strehl = std::numeric_limits<double>::quiet_NaN();
  f.coefficients = {1.5, -std::numeric_limits<double>::infinity()};
  const std::string text = rtt::io::to_result_json(f);
  CHECK(text.find(R"("strehl": {"float": "NaN"})") != std::string::npos);
  CHECK(
      text.find(
          R"("coefficients": {"dtype": "float64", "shape": [2], "values": [1.5, "-Infinity"]})") !=
      std::string::npos);
  CHECK_THROWS_AS(rtt::io::python_float_repr(f.strehl), std::invalid_argument);
}
