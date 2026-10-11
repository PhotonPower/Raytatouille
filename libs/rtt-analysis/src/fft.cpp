// STUB for the red run of #204 (tests first): checks the sizes but does not transform.
#include "fft.hpp"

#include <stdexcept>

namespace rtt::analysis::detail {

void fft2(std::span<std::complex<double>> data,
          std::size_t rows,
          std::size_t cols,
          FftDirection /*direction*/) {
  if (rows == 0 || cols == 0 || data.size() != rows * cols) {
    throw std::invalid_argument("fft2: rows and cols must be positive and match data.size()");
  }
}

void fft1(std::span<std::complex<double>> data, FftDirection direction) {
  fft2(data, 1, data.size(), direction);
}

}  // namespace rtt::analysis::detail
