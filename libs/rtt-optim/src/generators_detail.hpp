#pragma once

/// @file generators_detail.hpp
/// generator_residuals() with the nominal reference of a trial (ADR 0030, addendum #256),
/// private to rtt-optim.

#include <span>

#include "rtt/analysis/nominal_reference.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/model/optimization.hpp"
#include "rtt/optim/generators.hpp"

namespace rtt::optim::detail {

/// generator_residuals() (generators.hpp); with `nominal` the rays and the reference rays are
/// aimed with its aiming data and the wavefront has its reference sphere (#255, #265). Null:
/// exactly generator_residuals().
void generator_residuals(const model::Generator& generator,
                         const compile::CompiledSystem& system,
                         const analysis::NominalReference* nominal,
                         std::span<double> out,
                         GeneratorStats& stats);

}  // namespace rtt::optim::detail
