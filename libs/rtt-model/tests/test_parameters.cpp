// Evaluation of the parameter table (ADR 0029): row values per configuration, evaluating
// diagnostics with their pointers, resolve_parameters and for_each_param. Expected values are
// closed forms; the expression grammar itself is tested in test_expression.cpp.

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "rtt/model/parameters.hpp"
#include "test_support.hpp"

using namespace rtt::model;
using Catch::Matchers::ContainsSubstring;
using rtt::model::test::element;
using rtt::model::test::make_singlet;

namespace {

ParameterRow row(std::string name, ParameterForm form) {
  ParameterRow r;
  r.name = std::move(name);
  r.form = std::move(form);
  return r;
}

ParameterRow expr(std::string name, std::string text) {
  return row(std::move(name), ParameterExpression{std::move(text)});
}

/// The zoom of ADR 0029, point 8: TOTAL 40, G [20, 5] (variable), B = TOTAL - G.
System zoom() {
  System s = make_singlet();
  s.configurations = {{"wide"}, {"tele"}};
  ParameterRow g = row("G", std::vector<double>{20.0, 5.0});
  g.variable = true;
  g.min = 2.0;
  g.max = 30.0;
  s.parameters = {row("TOTAL", 40.0), g, expr("B", "TOTAL - G")};
  return s;
}

/// Diagnostics of evaluate_parameters with this code.
std::vector<Diagnostic> with_code(const std::vector<Diagnostic>& d, std::string_view code) {
  std::vector<Diagnostic> out;
  std::copy_if(d.begin(), d.end(), std::back_inserter(out),
               [&](const Diagnostic& x) { return x.code == code; });
  return out;
}

std::vector<Diagnostic> evaluating(const System& s) {
  std::vector<Diagnostic> d;
  static_cast<void>(evaluate_parameters(s, &d));
  return d;
}

}  // namespace

TEST_CASE("configurations: count and lookup", "[parameters]") {
  System s = make_singlet();
  CHECK(configuration_count(s) == 1);  // only the nominal configuration
  CHECK_FALSE(find_configuration(s, "").has_value());
  s.configurations = {{"wide"}, {"tele"}};
  CHECK(configuration_count(s) == 2);
  CHECK(find_configuration(s, "tele") == 1U);
  CHECK_FALSE(find_configuration(s, "Tele").has_value());  // case-sensitive
  s.parameters = {row("A", 1.0), row("B", 2.0), row("A", 3.0)};
  CHECK(find_parameter(s, "B") == 1U);
  CHECK(find_parameter(s, "A") == 0U);  // the first of a duplicate name
  CHECK_FALSE(find_parameter(s, "C").has_value());
}

TEST_CASE("table: the zoom of ADR 0029 per configuration", "[parameters]") {
  std::vector<Diagnostic> d;
  const ParameterValues v = evaluate_parameters(zoom(), &d);
  CHECK(d.empty());
  REQUIRE(v.rows == 3);
  REQUIRE(v.configurations == 2);
  REQUIRE(v.values.size() == 6);
  CHECK(v.at(0, 0) == 40.0);  // TOTAL, the same in both
  CHECK(v.at(0, 1) == 40.0);
  CHECK(v.at(1, 0) == 20.0);  // G
  CHECK(v.at(1, 1) == 5.0);
  CHECK(v.at(2, 0) == 20.0);  // B = TOTAL - G
  CHECK(v.at(2, 1) == 35.0);
}

TEST_CASE("table: chains over earlier rows, nominal configuration", "[parameters]") {
  System s = make_singlet();
  s.parameters = {row("R1", 51.68), expr("R2", "-R1"), expr("H", "(R1 - R2) / 4"),
                  expr("K", "H * H - 1e2"), row("V", std::vector<double>{3.0})};
  std::vector<Diagnostic> d;
  const ParameterValues v = evaluate_parameters(s, &d);
  CHECK(d.empty());
  REQUIRE(v.configurations == 1);
  CHECK(v.at(1, 0) == -51.68);
  CHECK(v.at(2, 0) == (51.68 - -51.68) / 4.0);
  CHECK(v.at(3, 0) == v.at(2, 0) * v.at(2, 0) - 100.0);
  CHECK(v.at(4, 0) == 3.0);
  CHECK(evaluate_parameters(System{}).rows == 0);  // no table
}

TEST_CASE("table: errors of expressions with pointer and position", "[parameters]") {
  System s = make_singlet();
  s.parameters = {row("A", 1.0),      expr("B", "A + later"),  row("later", 2.0),
                  expr("C", "C + 1"), expr("D", "A + nobody"), expr("E", "2 *"),
                  expr("F", "A + 1")};
  const std::vector<Diagnostic> d = evaluating(s);
  REQUIRE(d.size() == 4);  // B, C, D, E; F is fine
  CHECK(d[0].code == "parameters.forward_reference");
  CHECK(d[0].location == "/parameters/1/expression");
  CHECK_THAT(d[0].message, ContainsSubstring("later"));
  CHECK_THAT(d[0].message, ContainsSubstring("character 5"));
  CHECK(d[1].code == "parameters.forward_reference");  // the row itself
  CHECK(d[1].location == "/parameters/3/expression");
  CHECK(d[2].code == "parameters.unknown_name");
  CHECK(d[2].location == "/parameters/4/expression");
  CHECK_THAT(d[2].message, ContainsSubstring("nobody"));
  CHECK(d[3].code == "parameters.expression_syntax");
  CHECK(d[3].location == "/parameters/5/expression");
  for (const Diagnostic& x : d) CHECK(x.severity == Severity::Error);
  const ParameterValues v = evaluate_parameters(s);
  CHECK(std::isnan(v.at(1, 0)));
  CHECK(v.at(6, 0) == 2.0);
}

TEST_CASE("table: not finite in exactly one configuration", "[parameters]") {
  System s = make_singlet();
  s.configurations = {{"wide"}, {"tele"}};
  s.parameters = {row("X", std::vector<double>{1.0, 0.0}), expr("Y", "1 / X"), expr("Z", "Y + 1")};
  const std::vector<Diagnostic> d = evaluating(s);
  // Only Y in "tele"; Z uses Y and gets no diagnostic of its own.
  REQUIRE(d.size() == 1);
  CHECK(d[0].code == "parameters.not_finite");
  CHECK(d[0].location == "/parameters/1/expression");
  CHECK_THAT(d[0].message, ContainsSubstring("'tele'"));
  CHECK_THAT(d[0].message, !ContainsSubstring("'wide'"));
  const ParameterValues v = evaluate_parameters(s);
  CHECK(v.at(1, 0) == 1.0);
  CHECK(std::isinf(v.at(1, 1)));
  CHECK(std::isnan(v.at(2, 0)));  // Z: its input row is undefined
  CHECK(std::isnan(v.at(2, 1)));

  // In the nominal configuration the message says so.
  System n = make_singlet();
  n.parameters = {row("X", 0.0), expr("Y", "X / X")};
  const std::vector<Diagnostic> dn = evaluating(n);
  REQUIRE(dn.size() == 1);
  CHECK_THAT(dn[0].message, ContainsSubstring("nominal"));
}

TEST_CASE("table: undefined inputs give no follow-up diagnostics", "[parameters]") {
  System s = make_singlet();
  s.configurations = {{"a"}, {"b"}};
  s.parameters = {row("W", std::vector<double>{1.0}),  // values_count: one value for two
                  row("N", std::numeric_limits<double>::quiet_NaN()),  // value.not_finite
                  expr("P", "W + 1"),
                  expr("Q", "N * 0"),
                  expr("S", "1 +"),
                  expr("T", "S + P")};
  const std::vector<Diagnostic> d = evaluating(s);
  // Only the syntax error of S; the structural errors of W and N belong to validate().
  REQUIRE(d.size() == 1);
  CHECK(d[0].code == "parameters.expression_syntax");
  CHECK(d[0].location == "/parameters/4/expression");
  const ParameterValues v = evaluate_parameters(s);
  for (std::size_t r = 0; r < v.rows; ++r) {
    for (std::size_t k = 0; k < v.configurations; ++k) CHECK(std::isnan(v.at(r, k)));
  }
}

TEST_CASE("validate reports the evaluating codes once", "[parameters][validate]") {
  System s = zoom();
  s.parameters.push_back(expr("Q", "B / (G - 5)"));  // division by zero in "tele"
  s.parameters.push_back(expr("R", "TOTAL + unknown"));
  const std::vector<Diagnostic> d = validate(s);
  CHECK(with_code(d, "parameters.not_finite").size() == 1);
  CHECK(with_code(d, "parameters.unknown_name").size() == 1);
  CHECK(test::has_error_at(d, "/parameters/3/expression"));
  CHECK(test::has_error_at(d, "/parameters/4/expression"));
  CHECK(validate(zoom()).empty());
}

TEST_CASE("resolve: bound Params get the value of their configuration", "[parameters]") {
  System s = zoom();
  element(s, 1).pose.position[2] = Param::bound("G");
  element(s, 2).pose.position[2] = Param::bound("B");
  element(s, 0).pose.position[0] = Param(0.5);  // unbound: unchanged

  for (std::size_t k = 0; k < 2; ++k) {
    INFO("configuration " << k);
    System r = resolve_parameters(s, k);
    const Param& g = element(r, 1).pose.position[2];
    const Param& b = element(r, 2).pose.position[2];
    CHECK(g.value == (k == 0 ? 20.0 : 5.0));
    CHECK(b.value == (k == 0 ? 20.0 : 35.0));
    CHECK_FALSE(g.is_bound());
    CHECK_FALSE(g.variable);
    CHECK_FALSE(g.min.has_value());
    CHECK(element(r, 0).pose.position[0] == Param(0.5));
    // The table and the configurations stay.
    CHECK(r.parameters == s.parameters);
    CHECK(r.configurations == s.configurations);
    // The same system with the values entered by hand.
    System hand = s;
    element(hand, 1).pose.position[2] = Param(k == 0 ? 20.0 : 5.0);
    element(hand, 2).pose.position[2] = Param(k == 0 ? 20.0 : 35.0);
    CHECK(r == hand);
    CHECK(validate(r).empty());
  }
}

TEST_CASE("resolve: without bound Params the copy equals the input", "[parameters]") {
  CHECK(resolve_parameters(make_singlet(), 0) == make_singlet());
  CHECK(resolve_parameters(zoom(), 1) == zoom());
}

TEST_CASE("resolve: invalid configuration index or table", "[parameters]") {
  CHECK_THROWS_AS(resolve_parameters(make_singlet(), 1), std::invalid_argument);
  CHECK_THROWS_AS(resolve_parameters(zoom(), 2), std::invalid_argument);

  System s = zoom();
  s.parameters.push_back(expr("Q", "1 / (G - 5)"));                  // not finite in "tele"
  CHECK_THROWS_AS(resolve_parameters(s, 0), std::invalid_argument);  // any table error

  s = zoom();
  element(s, 1).pose.position[2] = Param::bound("nobody");  // param.unknown_parameter
  CHECK_THROWS_AS(resolve_parameters(s, 0), std::invalid_argument);

  s = zoom();
  s.configurations = {{"wide"}, {"wide"}};  // configurations.name_duplicate
  CHECK_THROWS_AS(resolve_parameters(s, 0), std::invalid_argument);

  // An error elsewhere in the model does not concern the table.
  s = zoom();
  s.wavelengths.clear();
  CHECK_NOTHROW(resolve_parameters(s, 0));
}

TEST_CASE("for_each_param: every Param of the singlet in edit-form order", "[parameters]") {
  System s = make_singlet();
  std::vector<std::string> pointers;
  for_each_param(static_cast<const System&>(s),
                 [&](std::string_view p, const Param&) { pointers.emplace_back(p); });
  // object distance, aperture, then 8 poses (root, 3 elements, 4 surfaces) with 6 Params each,
  // and radius and conic of L1.S1.
  REQUIRE(pointers.size() == 2 + 8 * 6 + 2);
  CHECK(pointers[0] == "/object/distance");
  CHECK(pointers[1] == "/aperture/value");
  CHECK(pointers[2] == "/root/pose/position/0");
  CHECK(pointers[5] == "/root/pose/rotation_deg/0");
  CHECK(pointers[8] == "/root/children/0/pose/position/0");
  CHECK(pointers[14] == "/root/children/0/surfaces/0/pose/position/0");
  CHECK(std::find(pointers.begin(), pointers.end(),
                  "/root/children/1/surfaces/0/shape/base/radius") != pointers.end());
  CHECK(pointers.back() == "/root/children/2/surfaces/0/pose/rotation_deg/2");

  // The non-const visitor changes the model.
  for_each_param(s, [](std::string_view p, Param& param) {
    if (p == "/root/children/1/surfaces/0/shape/base/radius") param.variable = true;
  });
  CHECK(std::get<Conic>(element(s, 1).surfaces[0].shape.base).radius.variable);
}
