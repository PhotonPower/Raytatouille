#include "rtt/model/tolerance_classes.hpp"

// Stub for the red run of #253 part B.

namespace rtt::model {

std::string_view tolerance_kind_name(ToleranceKind /*kind*/) noexcept {
  return {};
}

std::uint64_t tolerance_id_hash(std::string_view /*id*/) noexcept {
  return 0;
}

ClassTolerances class_tolerances(const System& /*system*/) {
  return {};
}

}  // namespace rtt::model
