#include "rtt/analysis/paths.hpp"

// STUB for the red run of the #122 tests: every function returns an empty result.

namespace rtt::analysis {

PathTransmission path_transmission(const compile::CompiledSystem& /*system*/,
                                   compile::PathId path,
                                   const trace::RayBatch& /*start*/,
                                   const PathOptions& /*options*/) {
  PathTransmission out;
  out.path = path;
  return out;
}

PathTransmission path_transmission(const compile::CompiledSystem& /*system*/,
                                   compile::PathId path,
                                   std::uint16_t /*field*/,
                                   std::uint16_t /*wavelength*/,
                                   const PathOptions& /*options*/) {
  PathTransmission out;
  out.path = path;
  return out;
}

PathOplDifference opl_difference(const compile::CompiledSystem& /*system*/,
                                 compile::PathId path_a,
                                 compile::PathId path_b,
                                 const trace::RayBatch& /*start*/,
                                 const PathOptions& /*options*/) {
  PathOplDifference out;
  out.path_a = path_a;
  out.path_b = path_b;
  return out;
}

PathOplDifference opl_difference(const compile::CompiledSystem& /*system*/,
                                 compile::PathId path_a,
                                 compile::PathId path_b,
                                 std::uint16_t /*field*/,
                                 std::uint16_t /*wavelength*/,
                                 const PathOptions& /*options*/) {
  PathOplDifference out;
  out.path_a = path_a;
  out.path_b = path_b;
  return out;
}

PathTransmission path_transmission(const compile::CompiledSystem& system,
                                   compile::PathId path,
                                   const trace::RayBatch& start,
                                   const PathOptions& options,
                                   const trace::RunControl& /*control*/) {
  return path_transmission(system, path, start, options);
}

PathTransmission path_transmission(const compile::CompiledSystem& system,
                                   compile::PathId path,
                                   std::uint16_t field,
                                   std::uint16_t wavelength,
                                   const PathOptions& options,
                                   const trace::RunControl& /*control*/) {
  return path_transmission(system, path, field, wavelength, options);
}

PathOplDifference opl_difference(const compile::CompiledSystem& system,
                                 compile::PathId path_a,
                                 compile::PathId path_b,
                                 const trace::RayBatch& start,
                                 const PathOptions& options,
                                 const trace::RunControl& /*control*/) {
  return opl_difference(system, path_a, path_b, start, options);
}

PathOplDifference opl_difference(const compile::CompiledSystem& system,
                                 compile::PathId path_a,
                                 compile::PathId path_b,
                                 std::uint16_t field,
                                 std::uint16_t wavelength,
                                 const PathOptions& options,
                                 const trace::RunControl& /*control*/) {
  return opl_difference(system, path_a, path_b, field, wavelength, options);
}

}  // namespace rtt::analysis
