// Parser for the expressions of the parameter table (ADR 0029, point 2).

#include "expression.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <variant>
#include <vector>

namespace rtt::model::detail {
namespace {

constexpr std::string_view kSyntax = "parameters.expression_syntax";

bool is_digit(char c) noexcept {
  return c >= '0' && c <= '9';
}
bool is_letter(char c) noexcept {
  return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}
bool is_name_start(char c) noexcept {
  return is_letter(c) || c == '_';
}
bool is_name_char(char c) noexcept {
  return is_name_start(c) || is_digit(c);
}
bool is_blank(char c) noexcept {
  return c == ' ' || c == '\t';
}

enum class TokenKind : std::uint8_t {
  Number,
  Name,
  Plus,
  Minus,
  Star,
  Slash,
  Open,
  Close,
  End,
};

struct Token {
  TokenKind kind = TokenKind::End;
  std::size_t position = 0;
  std::string_view text;  ///< the characters of the token (Number and Name)
};

/// A name used in the expression, resolved after the whole text has parsed, so that syntax
/// errors come first.
struct NameUse {
  std::string_view name;
  std::size_t position = 0;
  std::size_t instruction = 0;  ///< index of its Row instruction
};

std::string describe(const Token& t) {
  switch (t.kind) {
    case TokenKind::End:
      return "the end of the expression";
    case TokenKind::Number:
    case TokenKind::Name:
      return "'" + std::string(t.text) + "'";
    case TokenKind::Plus:
      return "'+'";
    case TokenKind::Minus:
      return "'-'";
    case TokenKind::Star:
      return "'*'";
    case TokenKind::Slash:
      return "'/'";
    case TokenKind::Open:
      return "'('";
    case TokenKind::Close:
      return "')'";
  }
  return "?";
}

/// Recursive-descent parser after the grammar of expression.hpp. Errors are returned, never
/// thrown; the first error stops parsing.
class Parser {
 public:
  explicit Parser(std::string_view text) : text_(text) {}

  std::variant<Expression, ExpressionError> run(const ResolveName& resolve) {
    if (text_.size() > kMaxExpressionLength) {
      return fail(kMaxExpressionLength,
                  "longer than " + std::to_string(kMaxExpressionLength) + " characters");
    }
    next();
    if (!error_) parse_expression(0);
    if (!error_ && current_.kind != TokenKind::End) {
      fail(current_.position, "unexpected " + describe(current_));
    }
    if (error_) return std::move(*error_);
    for (const NameUse& use : names_) {
      const NameLookup found = resolve(use.name);
      switch (found.status) {
        case NameStatus::Found:
          expression_.code[use.instruction].row = found.row;
          break;
        case NameStatus::Unknown:
          return ExpressionError{"parameters.unknown_name", use.position,
                                 "unknown name '" + std::string(use.name) + "' at character " +
                                     std::to_string(use.position + 1)};
        case NameStatus::Later:
          return ExpressionError{"parameters.forward_reference", use.position,
                                 "'" + std::string(use.name) + "' at character " +
                                     std::to_string(use.position + 1) +
                                     " is not an earlier row of the table"};
      }
    }
    return std::move(expression_);
  }

 private:
  ExpressionError fail(std::size_t position, const std::string& what) {
    if (!error_) {
      error_ = ExpressionError{std::string(kSyntax), position,
                               what + " at character " + std::to_string(position + 1)};
    }
    return *error_;
  }

  /// Reads the next token into current_; a malformed token sets error_.
  void next() {
    while (pos_ < text_.size() && is_blank(text_[pos_])) ++pos_;
    const std::size_t start = pos_;
    if (pos_ == text_.size()) {
      current_ = {TokenKind::End, start, {}};
      return;
    }
    const char c = text_[pos_];
    const auto single = [&](TokenKind kind) {
      ++pos_;
      current_ = {kind, start, text_.substr(start, 1)};
    };
    switch (c) {
      case '+':
        return single(TokenKind::Plus);
      case '-':
        return single(TokenKind::Minus);
      case '*':
        return single(TokenKind::Star);
      case '/':
        return single(TokenKind::Slash);
      case '(':
        return single(TokenKind::Open);
      case ')':
        return single(TokenKind::Close);
      default:
        break;
    }
    if (is_digit(c)) {
      read_number(start);
      return;
    }
    if (is_name_start(c)) {
      while (pos_ < text_.size() && is_name_char(text_[pos_])) ++pos_;
      current_ = {TokenKind::Name, start, text_.substr(start, pos_ - start)};
      return;
    }
    fail(start, "unexpected character");
    current_ = {TokenKind::End, start, {}};
  }

  /// number = int [ "." digit { digit } ] [ ("e" | "E") [ "+" | "-" ] digit { digit } ],
  /// int = "0" | nonzero { digit } (JSON, RFC 8259, section 6, without the sign).
  void read_number(std::size_t start) {
    const auto digits = [&] {
      const std::size_t first = pos_;
      while (pos_ < text_.size() && is_digit(text_[pos_])) ++pos_;
      return pos_ > first;
    };
    if (text_[pos_] == '0') {
      ++pos_;  // a leading zero stands alone; a following digit is the next token
    } else {
      digits();
    }
    if (pos_ < text_.size() && text_[pos_] == '.') {
      ++pos_;
      if (!digits()) {
        fail(start, "malformed number");
        current_ = {TokenKind::End, start, {}};
        return;
      }
    }
    if (pos_ < text_.size() && (text_[pos_] == 'e' || text_[pos_] == 'E')) {
      ++pos_;
      if (pos_ < text_.size() && (text_[pos_] == '+' || text_[pos_] == '-')) ++pos_;
      if (!digits()) {
        fail(start, "malformed number");
        current_ = {TokenKind::End, start, {}};
        return;
      }
    }
    current_ = {TokenKind::Number, start, text_.substr(start, pos_ - start)};
  }

  void emit(OpCode op, double constant = 0.0) {
    expression_.code.push_back({op, constant, 0});
    if (op == OpCode::Constant || op == OpCode::Row) {
      ++stack_;
      expression_.stack_size = std::max(expression_.stack_size, stack_);
    } else if (op != OpCode::Negate) {
      --stack_;  // binary: two operands become one
    }
  }

  /// expression = term { ("+" | "-") term }
  void parse_expression(std::size_t depth) {
    parse_term(depth);
    while (!error_ && (current_.kind == TokenKind::Plus || current_.kind == TokenKind::Minus)) {
      const OpCode op = current_.kind == TokenKind::Plus ? OpCode::Add : OpCode::Subtract;
      next();
      if (error_) return;
      parse_term(depth);
      if (!error_) emit(op);
    }
  }

  /// term = factor { ("*" | "/") factor }
  void parse_term(std::size_t depth) {
    parse_factor(depth);
    while (!error_ && (current_.kind == TokenKind::Star || current_.kind == TokenKind::Slash)) {
      const OpCode op = current_.kind == TokenKind::Star ? OpCode::Multiply : OpCode::Divide;
      next();
      if (error_) return;
      parse_factor(depth);
      if (!error_) emit(op);
    }
  }

  /// factor = "-" factor | primary; primary = number | name | "(" expression ")".
  /// `depth` counts the enclosing "(" and unary "-" (ADR 0029, point 2).
  void parse_factor(std::size_t depth) {
    if (error_) return;
    const Token t = current_;
    switch (t.kind) {
      case TokenKind::Minus: {
        if (depth + 1 > kMaxExpressionDepth) {
          fail(t.position, "nested deeper than " + std::to_string(kMaxExpressionDepth));
          return;
        }
        next();
        if (error_) return;
        parse_factor(depth + 1);
        if (!error_) emit(OpCode::Negate);
        return;
      }
      case TokenKind::Open: {
        if (depth + 1 > kMaxExpressionDepth) {
          fail(t.position, "nested deeper than " + std::to_string(kMaxExpressionDepth));
          return;
        }
        next();
        if (error_) return;
        parse_expression(depth + 1);
        if (error_) return;
        if (current_.kind != TokenKind::Close) {
          fail(current_.position, "expected ')' instead of " + describe(current_));
          return;
        }
        next();
        return;
      }
      case TokenKind::Number: {
        // std::from_chars: correctly rounded (to nearest, ties to even) and independent of the
        // locale (ADR 0029, point 2).
        double value = 0.0;
        const auto [end, ec] = std::from_chars(t.text.data(), t.text.data() + t.text.size(), value);
        if (ec != std::errc() || end != t.text.data() + t.text.size() || !std::isfinite(value)) {
          fail(t.position, "number " + describe(t) + " out of the range of double");
          return;
        }
        next();
        if (!error_) emit(OpCode::Constant, value);
        return;
      }
      case TokenKind::Name: {
        names_.push_back({t.text, t.position, expression_.code.size()});
        emit(OpCode::Row);
        next();
        return;
      }
      default:
        fail(t.position, "unexpected " + describe(t));
        return;
    }
  }

  std::string_view text_;
  std::size_t pos_ = 0;
  Token current_;
  Expression expression_;
  std::size_t stack_ = 0;
  std::vector<NameUse> names_;
  std::optional<ExpressionError> error_;
};

}  // namespace

std::variant<Expression, ExpressionError> parse_expression(std::string_view text,
                                                           const ResolveName& resolve) {
  return Parser(text).run(resolve);
}

}  // namespace rtt::model::detail
