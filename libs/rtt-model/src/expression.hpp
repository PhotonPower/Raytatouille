#pragma once

/// @file expression.hpp
/// Expressions of the parameter table (ADR 0029, point 2), internal to rtt-model: a parser that
/// turns the text of a derived row into a small postfix program, and its evaluation in double.
///
/// Grammar (ADR 0029, point 2):
///
///     expression = term { ("+" | "-") term } ;
///     term       = factor { ("*" | "/") factor } ;
///     factor     = "-" factor | primary ;
///     primary    = number | name | "(" expression ")" ;
///
/// Numbers follow the JSON number rule without a sign; names are [A-Za-z_][A-Za-z0-9_]*.
/// Spaces and tabs separate tokens. The parser never throws and never recurses deeper than
/// kMaxExpressionDepth levels of factor.

#include <cstddef>
#include <cstdint>
#include <functional>
#include <span>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace rtt::model::detail {

/// At most this many characters in an expression (ADR 0029, point 2).
inline constexpr std::size_t kMaxExpressionLength = 1000;
/// At most this many nested factors: every "(" and every unary "-" counts (ADR 0029, point 2).
inline constexpr std::size_t kMaxExpressionDepth = 64;

/// Result of looking up a name used in an expression.
enum class NameStatus : std::uint8_t {
  Found,    ///< an earlier row; `row` is its index
  Unknown,  ///< no row of this name
  Later,    ///< the row itself or a later one (forward reference)
};

struct NameLookup {
  NameStatus status = NameStatus::Unknown;
  std::size_t row = 0;  ///< index of the row if Found
};

/// Looks up a name in the table, as seen from the row whose expression is parsed.
using ResolveName = std::function<NameLookup(std::string_view name)>;

enum class OpCode : std::uint8_t { Constant, Row, Add, Subtract, Multiply, Divide, Negate };

struct Instruction {
  OpCode op = OpCode::Constant;
  double constant = 0.0;  ///< for Constant
  std::size_t row = 0;    ///< for Row: index into the row values of evaluate()
};

/// A parsed expression: a postfix program, evaluated left to right with a stack.
struct Expression {
  std::vector<Instruction> code;
  std::size_t stack_size = 0;  ///< largest stack depth that evaluate() needs
};

/// A parse error. `code` is a diagnostic code of the registry (parameters.expression_syntax,
/// parameters.unknown_name or parameters.forward_reference); `position` is the 0-based offset
/// of the offending token in the text; `message` names the token and the 1-based position.
struct ExpressionError {
  std::string code;
  std::size_t position = 0;
  std::string message;
};

/// Parses `text` and resolves its names with `resolve`. Never throws (except std::bad_alloc).
[[nodiscard]] std::variant<Expression, ExpressionError> parse_expression(
    std::string_view text, const ResolveName& resolve);

/// Evaluates `expression` in double with the values of the rows (index = Instruction::row):
/// every operation in the order of the grammar, left-associative, unary minus as exact
/// negation, without reassociation and without FMA contraction (expression_eval.cpp is compiled
/// with -ffp-contract=off; ADR 0029, point 2). The result may be infinite or NaN; the caller
/// reports that as parameters.not_finite.
/// @pre every Row instruction indexes into `rows`
[[nodiscard]] double evaluate(const Expression& expression, std::span<const double> rows);

}  // namespace rtt::model::detail
