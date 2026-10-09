// Evaluation of the expressions of the parameter table (ADR 0029, point 2). This translation unit
// is compiled without FMA contraction (-ffp-contract=off, see CMakeLists.txt), and it holds
// nothing but the evaluation, so that no caller with other options inlines the arithmetic:
// every operation is rounded on its own, in the order of the postfix program, which gives the
// same bits on every platform.

#include <cstddef>
#include <span>
#include <vector>

#include "expression.hpp"

namespace rtt::model::detail {

double evaluate(const Expression& expression, std::span<const double> rows) {
  std::vector<double> stack;
  stack.reserve(expression.stack_size);
  for (const Instruction& in : expression.code) {
    switch (in.op) {
      case OpCode::Constant:
        stack.push_back(in.constant);
        break;
      case OpCode::Row:
        stack.push_back(rows[in.row]);
        break;
      case OpCode::Negate:
        stack.back() = -stack.back();
        break;
      default: {
        const double right = stack.back();
        stack.pop_back();
        double& left = stack.back();
        switch (in.op) {
          case OpCode::Add:
            left = left + right;
            break;
          case OpCode::Subtract:
            left = left - right;
            break;
          case OpCode::Multiply:
            left = left * right;
            break;
          case OpCode::Divide:
            left = left / right;
            break;
          default:
            break;
        }
        break;
      }
    }
  }
  return stack.empty() ? 0.0 : stack.back();
}

}  // namespace rtt::model::detail
