#include "rtt/analysis/field.hpp"

namespace rtt::analysis {

std::vector<DistortionPoint> distortion(const compile::CompiledSystem& /*system*/,
                                        compile::PathId /*path*/,
                                        std::uint16_t /*wavelength*/,
                                        const FieldSweepOptions& /*options*/) {
  return {};  // TODO(issue-31): implement
}

DistortionPoint distortion_at(const compile::CompiledSystem& /*system*/,
                              compile::PathId /*path*/,
                              const model::Field& /*field*/,
                              std::uint16_t /*wavelength*/,
                              trace::Aiming /*aiming*/) {
  return {};  // TODO(issue-31): implement
}

std::vector<FieldCurvaturePoint> field_curvature(const compile::CompiledSystem& /*system*/,
                                                 compile::PathId /*path*/,
                                                 std::uint16_t /*wavelength*/,
                                                 const FieldCurvatureOptions& /*options*/) {
  return {};  // TODO(issue-31): implement
}

FieldCurvaturePoint field_curvature_at(const compile::CompiledSystem& /*system*/,
                                       compile::PathId /*path*/,
                                       const model::Field& /*field*/,
                                       std::uint16_t /*wavelength*/,
                                       const FieldCurvatureOptions& /*options*/) {
  return {};  // TODO(issue-31): implement
}

}  // namespace rtt::analysis
