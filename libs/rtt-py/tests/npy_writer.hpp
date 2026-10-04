#pragma once

// .npy output of rtt_py_reference (test program of rtt-py, not part of the library).

#include <bit>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace rtt::py::reference {

static_assert(std::endian::native == std::endian::little, ".npy descriptors below are '<'");

/// Writes a C-contiguous array as NumPy .npy file, format version 1.0.
inline void write_npy(const std::filesystem::path& file,
                      std::string_view descr,
                      const std::vector<std::size_t>& shape,
                      std::span<const char> bytes) {
  std::string dims;
  for (const std::size_t d : shape) dims += std::to_string(d) + ", ";
  if (shape.size() > 1) dims.resize(dims.size() - 2);  // "(N,)" for 1-D, "(N, 3, 3)" else
  std::string header = "{'descr': '" + std::string(descr) + "', 'fortran_order': False, " +
                       "'shape': (" + dims + "), }";
  // Magic (6) + version (2) + header length (2) + header + '\n' is a multiple of 64.
  const std::size_t total = 10 + header.size() + 1;
  header.append((64 - total % 64) % 64, ' ');
  header += '\n';
  std::ofstream out(file, std::ios::binary);
  out.write("\x93NUMPY\x01\x00", 8);
  const auto length = static_cast<std::uint16_t>(header.size());
  const char length_bytes[2] = {static_cast<char>(length & 0xffU), static_cast<char>(length >> 8U)};
  out.write(length_bytes, 2);
  out.write(header.data(), static_cast<std::streamsize>(header.size()));
  out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
  if (!out) throw std::runtime_error("cannot write " + file.string());
}

template <typename T>
void write_column(const std::filesystem::path& file,
                  std::string_view descr,
                  std::span<const T> values) {
  const std::span<const std::byte> bytes = std::as_bytes(values);
  write_npy(file, descr, {values.size()},
            {reinterpret_cast<const char*>(bytes.data()), bytes.size()});
}

/// Analysis cases of test_bitwise_analysis.py; writes <out>/<case>.<name>.npy.
void run_analysis_cases(const std::filesystem::path& reference_dir,
                        const std::filesystem::path& catalog_dir,
                        const std::filesystem::path& out,
                        int threads);

}  // namespace rtt::py::reference
