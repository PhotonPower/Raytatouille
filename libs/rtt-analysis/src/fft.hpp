#pragma once

/// @file fft.hpp
/// Internal discrete Fourier transform of rtt-analysis (ADR 0033, point 1): a thin wrapper
/// around pocketfft, which only fft.cpp includes. Not part of the public API.

#include <complex>
#include <cstddef>
#include <span>

namespace rtt::analysis::detail {

/// Direction of the transform.
enum class FftDirection {
  Forward,   ///< y_k = sum_j x_j exp(-2 pi i j k / n) along each axis
  Backward,  ///< y_k = sum_j x_j exp(+2 pi i j k / n) along each axis
};

/// Two-dimensional complex DFT in place, unnormalised: Backward(Forward(x)) = rows * cols * x.
/// `data` is row-major, rows x cols (index r * cols + c); both axes are transformed. Any length
/// (pocketfft: mixed radix, Bluestein for large prime factors). Single-threaded, without global
/// state; the result depends only on the input (bitwise the same on every call and thread).
/// @throws std::invalid_argument if rows or cols is 0 or data.size() != rows * cols
void fft2(std::span<std::complex<double>> data,
          std::size_t rows,
          std::size_t cols,
          FftDirection direction);

/// One-dimensional complex DFT in place: fft2 with one row.
/// @throws std::invalid_argument if data is empty
void fft1(std::span<std::complex<double>> data, FftDirection direction);

}  // namespace rtt::analysis::detail
