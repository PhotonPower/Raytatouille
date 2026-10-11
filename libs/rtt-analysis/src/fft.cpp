// The only translation unit that includes pocketfft (ADR 0033, point 1; ODR): without its thread
// pool, without SIMD paths that depend on -m flags, and without the plan cache, so the dependency
// holds no mutable global state (rule 4) and the result depends only on the input.
//
// GCC 16 at -O2 reports -Wnull-dereference inside the header (pocketfft_hdronly.h:1495 and :266,
// after inlining into fft2): arr<T> allocates with malloc, which may return null for a size of 0,
// a path that cannot occur there (the radix ip >= 2). The SYSTEM include does not cover this
// middle-end warning, so the suppression covers this translation unit, which holds only the
// wrapper; limited to GCC, as in rtt-polar/src/prt_analysis.cpp.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnull-dereference"
#endif

#define POCKETFFT_NO_MULTITHREADING
#define POCKETFFT_NO_VECTORS
#define POCKETFFT_CACHE_SIZE 0

#include "fft.hpp"

#include <pocketfft_hdronly.h>

#include <cstddef>
#include <stdexcept>

namespace rtt::analysis::detail {

void fft2(std::span<std::complex<double>> data,
          std::size_t rows,
          std::size_t cols,
          FftDirection direction) {
  if (rows == 0 || cols == 0 || data.size() != rows * cols) {
    throw std::invalid_argument("fft2: rows and cols must be positive and match data.size()");
  }
  const pocketfft::shape_t shape{rows, cols};
  // Strides in bytes, row-major.
  const auto row_stride = static_cast<std::ptrdiff_t>(cols * sizeof(std::complex<double>));
  const pocketfft::stride_t stride{row_stride,
                                   static_cast<std::ptrdiff_t>(sizeof(std::complex<double>))};
  const pocketfft::shape_t axes{0, 1};
  // pocketfft FORWARD is sum_j x_j exp(-2 pi i j k / n) per axis (checked in test_fft.cpp); in
  // place, factor 1, one thread.
  pocketfft::c2c(shape, stride, stride, axes, direction == FftDirection::Forward, data.data(),
                 data.data(), 1.0, 1);
}

void fft1(std::span<std::complex<double>> data, FftDirection direction) {
  fft2(data, 1, data.size(), direction);
}

}  // namespace rtt::analysis::detail

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
