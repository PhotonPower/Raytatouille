#include "rtt/analysis/stokes.hpp"

#include <complex>

namespace rtt::analysis {

// STUB for the red run of #211: replaced by the implementation in the next commit.
StokesPupil stokes_pupil(const JonesPupil& /*jones*/,
                         std::complex<double> /*e1*/,
                         std::complex<double> /*e2*/) {
  return {};
}

}  // namespace rtt::analysis
