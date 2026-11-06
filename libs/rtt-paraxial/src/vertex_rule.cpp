#include "vertex_rule.hpp"

#include <algorithm>
#include <cmath>
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

geom::VertexExpansion<double> vertex_expansion(const compile::CompiledSurface& s) {
  return std::visit([](const auto& shape) { return shape.vertex_expansion(); }, s.shape);
}

VertexRule vertex_rule(const compile::CompiledSurface& s) {
  return vertex_rule(vertex_expansion(s));
}

bool is_anamorphic(const geom::VertexExpansion<double>& e) {
  return e.order >= 2 && (e.hessian[0] != e.hessian[2] || e.hessian[1] != 0.0);
}

SectionFrame section_frame(const math::Mat3& rotation) {
  return {rotation(0, 0), rotation(0, 1), rotation(1, 0), rotation(1, 1),
          rotation(2, 2) > 0.0 ? 1.0 : -1.0};
}

SectionRule section_rule(const geom::VertexExpansion<double>& e,
                         const SectionFrame& frame,
                         Section section) {
  SectionRule r;
  if (!is_anamorphic(e) || e.slope.first != 0.0 || e.slope.second != 0.0) {
    // Rotationally symmetric in second order, or rejected by the vertex rule (order < 2, slope
    // at the vertex) with its message: c exactly as first_order() (ADR 0042, point 2).
    const VertexRule v = vertex_rule(e);
    if (v.first_order_error) {
      r.error = *v.first_order_error + " (ADR 0041)";
      return r;
    }
    r.value = v.value;
    r.c = frame.sigma > 0.0 ? v.c : -v.c;
    r.c_scale = r.c;
    return r;
  }
  // Global Hessian H_g = sigma Q H_l Q^T (ADR 0042, point 2, derived there from
  // z_g - o_z = sigma z_l(Q^T (x_g - o_x, y_g - o_y))): with rows q_x = (q00, q01) and
  // q_y = (q10, q11) of Q, H_g,ab = sigma q_a H_l q_b^T.
  const double hxx = e.hessian[0];
  const double hxy = e.hessian[1];
  const double hyy = e.hessian[2];
  const double gxx =
      frame.sigma * (frame.q00 * frame.q00 * hxx + 2.0 * frame.q00 * frame.q01 * hxy +
                     frame.q01 * frame.q01 * hyy);
  const double gyy =
      frame.sigma * (frame.q10 * frame.q10 * hxx + 2.0 * frame.q10 * frame.q11 * hxy +
                     frame.q11 * frame.q11 * hyy);
  const double gxy = frame.sigma * (frame.q00 * frame.q10 * hxx +
                                    (frame.q00 * frame.q11 + frame.q01 * frame.q10) * hxy +
                                    frame.q01 * frame.q11 * hyy);
  if (!(std::abs(gxy) <= kSymmetryTolerance * std::max(std::abs(gxx), std::abs(gyy)))) {
    r.error = "principal axes not aligned with x and y; no paraxial model per section (ADR 0042)";
    return r;
  }
  r.value = e.value;
  r.c = section == Section::x ? gxx : gyy;
  r.c_scale = std::abs(hxx) + 2.0 * std::abs(hxy) + std::abs(hyy);
  return r;
}

CylinderLensSection cylinder_lens_section(const math::Vec3& power_axis,
                                          const SectionFrame& frame,
                                          Section section) {
  CylinderLensSection r;
  // b_g = Q b_l (ADR 0042, point 2); b_l lies in the local x-y plane (b_l,z = 0).
  const double bx = frame.q00 * power_axis.x() + frame.q01 * power_axis.y();
  const double by = frame.q10 * power_axis.x() + frame.q11 * power_axis.y();
  if (!(std::abs(bx * by) <= kSymmetryTolerance)) {
    r.error =
        "principal axes not aligned with x and y: the direction of power of the ideal cylinder "
        "lens lies along neither; no paraxial model per section (ADR 0042)";
    return r;
  }
  r.factor = section == Section::x ? bx * bx : by * by;
  return r;
}

double vertex_z(double origin_z, double axis_z, double value) {
  if (value == 0.0) return origin_z;
  return axis_z > 0.0 ? origin_z + value : origin_z - value;
}

}  // namespace rtt::paraxial::detail
