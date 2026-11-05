#include "vertex_rule.hpp"

#include <variant>

namespace rtt::paraxial::detail {

VertexRule vertex_rule(const geom::VertexExpansion<double>& e) {
  VertexRule r;
  const auto no_model = [&r](const std::string& why) {
    r.first_order_error = why;
    r.fourth_order_error = why;
    return r;
  };
  if (e.order < 2) {
    return no_model("no paraxial model at the vertex (the sag is not twice differentiable there)");
  }
  if (e.slope.first != 0.0 || e.slope.second != 0.0) {
    return no_model(
        "slope at the vertex (tilt terms of the sag); no rotationally symmetric "
        "paraxial model");
  }
  if (e.hessian[0] != e.hessian[2] || e.hessian[1] != 0.0) {
    return no_model("not rotationally symmetric in second order; no paraxial model");
  }
  r.value = e.value;
  r.c = e.hessian[0];
  if (e.order < 4) {
    r.fourth_order_error = "the sag has no expansion to fourth order at the vertex";
  } else if (!e.rotational4) {
    r.fourth_order_error = "not rotationally symmetric to fourth order";
  } else {
    r.conic = e.conic;
    r.a4 = e.a4;
  }
  return r;
}

VertexRule vertex_rule(const compile::CompiledSurface& s) {
  return vertex_rule(
      std::visit([](const auto& shape) { return shape.vertex_expansion(); }, s.shape));
}

double vertex_z(double origin_z, double axis_z, double value) {
  if (value == 0.0) return origin_z;
  return axis_z > 0.0 ? origin_z + value : origin_z - value;
}

}  // namespace rtt::paraxial::detail
