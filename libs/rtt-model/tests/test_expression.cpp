// Expressions of the parameter table (ADR 0029, point 2): grammar, correctly rounded literals,
// evaluation order without FMA contraction, errors with code and position, and the limits.

#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "expression.hpp"

using namespace rtt::model::detail;

namespace {

/// A table with the rows a, b, c (indices 0, 1, 2) before the current row "self" (3) and a
/// later row "later" (4).
NameLookup lookup(std::string_view name) {
  static const std::map<std::string, std::size_t, std::less<>> earlier = {
      {"a", 0}, {"b", 1}, {"c", 2}};
  if (const auto it = earlier.find(name); it != earlier.end()) {
    return {NameStatus::Found, it->second};
  }
  if (name == "self" || name == "later") return {NameStatus::Later, 0};
  return {NameStatus::Unknown, 0};
}

const std::vector<double> kRows = {2.0, 3.0, 4.0};

/// Parses and evaluates `text` with a = 2, b = 3, c = 4; fails the test on a parse error.
double value_of(std::string_view text, const std::vector<double>& rows = kRows) {
  const auto parsed = parse_expression(text, lookup);
  if (const auto* error = std::get_if<ExpressionError>(&parsed)) {
    FAIL("unexpected error for '" << text << "': " << error->code << ": " << error->message);
  }
  return evaluate(std::get<Expression>(parsed), rows);
}

/// The parse error of `text`; fails the test if `text` parses.
ExpressionError error_of(std::string_view text) {
  const auto parsed = parse_expression(text, lookup);
  if (!std::holds_alternative<ExpressionError>(parsed)) FAIL("'" << text << "' parsed");
  return std::get<ExpressionError>(parsed);
}

std::uint64_t bits(double x) {
  return std::bit_cast<std::uint64_t>(x);
}

}  // namespace

TEST_CASE("expression: precedence, parentheses and left associativity", "[expression]") {
  CHECK(value_of("2 + 3 * 4") == 14.0);
  CHECK(value_of("(2 + 3) * 4") == 20.0);
  CHECK(value_of("8 / 4 / 2") == 1.0);   // (8 / 4) / 2, not 8 / (4 / 2) = 4
  CHECK(value_of("10 - 4 - 3") == 3.0);  // (10 - 4) - 3, not 10 - (4 - 3) = 9
  CHECK(value_of("a - b - c") == -5.0);  // (2 - 3) - 4
  CHECK(value_of("a * b + c") == 10.0);
  CHECK(value_of("a + b * c") == 14.0);
  CHECK(value_of("c / a * b") == 6.0);  // (4 / 2) * 3
  CHECK(value_of("((a))") == 2.0);
  CHECK(value_of("\ta\t+  b ") == 5.0);  // spaces and tabs between tokens
  CHECK(value_of("a") == 2.0);
  CHECK(value_of("0") == 0.0);
}

TEST_CASE("expression: unary minus is exact negation and nests", "[expression]") {
  CHECK(value_of("-2 * -3") == 6.0);
  CHECK(value_of("--a") == 2.0);
  CHECK(value_of("---a") == -2.0);
  CHECK(value_of("-(a + b)") == -5.0);
  CHECK(value_of("a - -b") == 5.0);
  CHECK(value_of("-a * b") == -6.0);
  // Negation of zero gives -0.0, bit for bit.
  CHECK(bits(value_of("-0")) == bits(-0.0));
  CHECK(bits(value_of("-(a - a)")) == bits(-0.0));
}

TEST_CASE("expression: literals follow the JSON number rule", "[expression]") {
  CHECK(value_of("1e-3") == 1e-3);
  CHECK(value_of("2.5E+2") == 250.0);
  CHECK(value_of("1E2") == 100.0);
  CHECK(value_of("0.5") == 0.5);
  CHECK(value_of("0e5") == 0.0);
  CHECK(value_of("123456789") == 123456789.0);
}

TEST_CASE("expression: literals are correctly rounded (ties to even)", "[expression]") {
  // 0.1 is the double 0x3FB999999999999A.
  CHECK(bits(value_of("0.1")) == 0x3FB999999999999AULL);
  // 2^53 + 1 lies exactly halfway between 2^53 and 2^53 + 2; ties to even gives 2^53.
  CHECK(value_of("9007199254740993") == 9007199254740992.0);
  // 1 + 2^-53 lies exactly halfway between 1 and 1 + 2^-52; ties to even gives 1.
  CHECK(value_of("1.00000000000000011102230246251565404236316680908203125") == 1.0);
  // Just above that halfway point the literal rounds up to 1 + 2^-52.
  CHECK(value_of("1.00000000000000011102230246251565404236316680908203126") ==
        1.0 + std::ldexp(1.0, -52));
  // 2^53 + 3 lies halfway between 2^53 + 2 and 2^53 + 4; ties to even gives 2^53 + 4.
  CHECK(value_of("9007199254740995") == 9007199254740996.0);
}

TEST_CASE("expression: no FMA contraction in the evaluation", "[expression]") {
  // a = 1 + 2^-30 and b = 1 - 2^-30 are exact; a * b = 1 - 2^-60 rounds to 1, so a * b - 1 is
  // exactly 0. A fused multiply-add would give -2^-60 instead.
  const std::vector<double> rows = {1.000000000931322574615478515625,
                                    0.999999999068677425384521484375, 1.0};
  REQUIRE(rows[0] == 1.0 + std::ldexp(1.0, -30));
  REQUIRE(rows[1] == 1.0 - std::ldexp(1.0, -30));
  CHECK(bits(value_of("a * b - c", rows)) == bits(0.0));
  CHECK(bits(value_of("a * b - 1", rows)) == bits(0.0));
  // The same with the literals in the text.
  CHECK(bits(value_of("1.000000000931322574615478515625 * 0.999999999068677425384521484375 - 1")) ==
        bits(0.0));
}

TEST_CASE("expression: non-finite results come back to the caller", "[expression]") {
  const std::vector<double> rows = {1.0, 0.0, 1e308};
  CHECK(std::isinf(value_of("a / b", rows)));
  CHECK(std::isnan(value_of("b / b", rows)));
  CHECK(std::isinf(value_of("c * 10", rows)));
}

TEST_CASE("expression: names resolve to earlier rows only", "[expression]") {
  const ExpressionError unknown = error_of("a + d");
  CHECK(unknown.code == "parameters.unknown_name");
  CHECK(unknown.position == 4);
  const ExpressionError later = error_of("2 * later");
  CHECK(later.code == "parameters.forward_reference");
  CHECK(later.position == 4);
  const ExpressionError self = error_of("self + 1");
  CHECK(self.code == "parameters.forward_reference");
  CHECK(self.position == 0);
  // Names are case sensitive.
  CHECK(error_of("A").code == "parameters.unknown_name");
  // Syntax errors come before name errors: "sqrt" is unknown, but the call is the error.
  CHECK(error_of("sqrt(2)").code == "parameters.expression_syntax");
}

TEST_CASE("expression: names are [A-Za-z_][A-Za-z0-9_]*", "[expression]") {
  const ResolveName names = [](std::string_view name) -> NameLookup {
    if (name == "_x1") return {NameStatus::Found, 0};
    if (name == "Zoom_2b") return {NameStatus::Found, 1};
    return {NameStatus::Unknown, 0};
  };
  const auto parsed = parse_expression("_x1 * Zoom_2b", names);
  REQUIRE(std::holds_alternative<Expression>(parsed));
  CHECK(evaluate(std::get<Expression>(parsed), std::vector<double>{7.0, 3.0}) == 21.0);
  // The lookup sees each name exactly as written.
  std::vector<std::string> seen;
  const ResolveName record = [&seen](std::string_view name) -> NameLookup {
    seen.emplace_back(name);
    return {NameStatus::Found, 0};
  };
  REQUIRE(std::holds_alternative<Expression>(parse_expression("a_1+B*(c)", record)));
  CHECK(seen == std::vector<std::string>{"a_1", "B", "c"});
}

TEST_CASE("expression: syntax errors with position", "[expression]") {
  struct Case {
    std::string_view text;
    std::size_t position;
  };
  const Case cases[] = {
      {"", 0},                       // empty
      {"   ", 3},                    // only blanks: the end
      {"2 +", 3},                    // missing operand at the end
      {"(2", 2},                     // missing ")"
      {"2)", 1},                     // stray ")"
      {"2 3", 2},                    // two operands
      {"a b", 2},     {"2a", 1},     // number directly followed by a name
      {".5", 0},                     // no leading digit
      {"1.", 0},                     // no digit after the point
      {"01", 1},                     // leading zero followed by a digit
      {"+2", 0},                     // no unary plus
      {"2 * +3", 4},  {"a(", 1},     // no function calls
      {"sqrt(2)", 4}, {"2 ^ 3", 2},  // unknown operator
      {"1e", 0},                     // exponent without digits
      {"1e+", 0},     {"2 * * 3", 4},
      {"()", 1},      {"2\n+ 3", 1},  // only spaces and tabs separate tokens
      {"1e400", 0},                   // out of the range of double
      {"Ä", 0},                       // ASCII only
  };
  for (const Case& c : cases) {
    INFO("text: '" << c.text << "'");
    const ExpressionError e = error_of(c.text);
    CHECK(e.code == "parameters.expression_syntax");
    CHECK(e.position == c.position);
    CHECK(!e.message.empty());
  }
}

TEST_CASE("expression: depth limit counts parentheses and unary minus", "[expression]") {
  const auto nested = [](std::size_t n) { return std::string(n, '(') + "a" + std::string(n, ')'); };
  CHECK(value_of(nested(kMaxExpressionDepth)) == 2.0);
  CHECK(error_of(nested(kMaxExpressionDepth + 1)).code == "parameters.expression_syntax");

  const auto minus = [](std::size_t n) { return std::string(n, '-') + "a"; };
  CHECK(value_of(minus(kMaxExpressionDepth)) == 2.0);  // 64 is even
  CHECK(error_of(minus(kMaxExpressionDepth + 1)).code == "parameters.expression_syntax");

  // Both count together: 32 parentheses around 32 minus signs is 64, one more is too deep.
  const std::string mixed =
      std::string(32, '(') + std::string(32, '-') + "a" + std::string(32, ')');
  CHECK(value_of(mixed) == 2.0);
  CHECK(error_of("-" + mixed).code == "parameters.expression_syntax");

  // Depth is nesting, not count: many factors side by side are fine.
  std::string sum = "a";
  for (int i = 0; i < 100; ++i) sum += " + (a)";
  CHECK(value_of(sum) == 202.0);
}

TEST_CASE("expression: length limit", "[expression]") {
  // "a" followed by " + a" until the length is reached exactly.
  std::string text = "a";
  while (text.size() + 4 <= kMaxExpressionLength) text += " + a";
  text += std::string(kMaxExpressionLength - text.size(), ' ');
  REQUIRE(text.size() == kMaxExpressionLength);
  CHECK(std::isfinite(value_of(text)));
  const ExpressionError e = error_of(text + " ");
  CHECK(e.code == "parameters.expression_syntax");
  CHECK(e.position == kMaxExpressionLength);
}

TEST_CASE("expression: the postfix program and its stack size", "[expression]") {
  const auto parsed = parse_expression("a - (b - c) * 2", lookup);
  REQUIRE(std::holds_alternative<Expression>(parsed));
  const Expression& e = std::get<Expression>(parsed);
  // a b c - 2 * -
  REQUIRE(e.code.size() == 7);
  CHECK(e.code[0].op == OpCode::Row);
  CHECK(e.code[0].row == 0);
  CHECK(e.code[3].op == OpCode::Subtract);
  CHECK(e.code[5].op == OpCode::Multiply);
  CHECK(e.code[6].op == OpCode::Subtract);
  CHECK(e.stack_size == 3);
  CHECK(evaluate(e, kRows) == 4.0);  // 2 - (3 - 4) * 2
}
