#include "rtt/trace/sources.hpp"

namespace rtt::trace {

std::vector<PupilPoint> pupil_points(const PupilSampling& /*sampling*/) {
  return {};  // TODO(issue-8): implement
}

AimedRay aim_ray(const compile::CompiledSystem& /*system*/,
                 compile::PathId /*path*/,
                 std::uint16_t /*field*/,
                 std::uint16_t /*wavelength*/,
                 double /*px*/,
                 double /*py*/,
                 Aiming /*aiming*/) {
  return {};  // TODO(issue-8): implement
}

RayBatch make_rays(const compile::CompiledSystem& /*system*/,
                   compile::PathId /*path*/,
                   std::span<const std::uint16_t> /*fields*/,
                   std::uint16_t /*wavelength*/,
                   const PupilSampling& /*sampling*/,
                   Aiming /*aiming*/) {
  return {};  // TODO(issue-8): implement
}

}  // namespace rtt::trace
