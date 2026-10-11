// STUB for the red run of #205 (tests first): links, but computes nothing.
#include "rtt/analysis/psf.hpp"

namespace rtt::analysis {

Psf psf(const compile::CompiledSystem& /*system*/,
        compile::PathId /*path*/,
        std::uint16_t /*field*/,
        const PsfOptions& /*options*/) {
  return {};
}

Psf psf(const compile::CompiledSystem& /*system*/,
        compile::PathId /*path*/,
        std::uint16_t /*field*/,
        const PsfOptions& /*options*/,
        const trace::RunControl& /*control*/) {
  return {};
}

std::vector<double> psf_at(const PupilFunction& /*pupil*/, std::span<const Point2> x) {
  return std::vector<double>(x.size(), 0.0);
}

std::vector<double> through_focus(const compile::CompiledSystem& /*system*/,
                                  compile::PathId /*path*/,
                                  std::uint16_t /*field*/,
                                  std::uint16_t /*wavelength*/,
                                  std::span<const double> shifts,
                                  const PupilFunctionOptions& /*options*/) {
  return std::vector<double>(shifts.size(), 0.0);
}

}  // namespace rtt::analysis
