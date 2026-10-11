#include "rtt/analysis/psf.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "common.hpp"
#include "fft.hpp"
#include "rtt/diagnostics/codes.hpp"

namespace rtt::analysis {

namespace {

using compile::CompiledSystem;
using compile::PathId;
using Complex = std::complex<double>;

/// Largest FFT length (ADR 0033, point 11: memory guard). L^2 complex values are 1 GiB at 8192;
/// psf() keeps one transform at a time (the pupil's is freed before the ideal one of
/// psf.aliasing) plus the L^2 intensities of the reference lattice (0.5 GiB), so about 1.5 GiB at
/// the limit (review of #233).
constexpr int kMaxFftLength = 8192;

std::string path_pointer(PathId path) {
  return "/paths/" + std::to_string(path.index);
}

void check_options(const PsfOptions& o) {
  if (o.padding < 2) throw std::invalid_argument("analysis: PSF padding must be >= 2");
  if (!(o.aliasing_factor >= 1.0) || !std::isfinite(o.aliasing_factor)) {
    throw std::invalid_argument("analysis: aliasing_factor must be finite and >= 1");
  }
  if (!(o.undersampled_waves > 0.0) || !std::isfinite(o.undersampled_waves)) {
    throw std::invalid_argument("analysis: undersampled_waves must be finite and > 0");
  }
}

/// The chosen wavelengths with their model weights normalised to sum 1 over the
/// choice.
std::vector<std::pair<std::uint16_t, double>> chosen_wavelengths(
    const CompiledSystem& system, const std::optional<std::vector<std::uint16_t>>& list) {
  std::vector<std::uint16_t> indices;
  if (list) {
    if (list->empty()) throw std::invalid_argument("analysis: the wavelength list is empty");
    for (const std::uint16_t w : *list) {
      detail::check_wavelength(system, w);
      if (std::find(indices.begin(), indices.end(), w) != indices.end()) {
        throw std::invalid_argument("analysis: wavelength " + std::to_string(w) +
                                    " is chosen more than once");
      }
      indices.push_back(w);
    }
  } else {
    for (std::size_t i = 0; i < system.wavelengths_um().size(); ++i) {
      indices.push_back(static_cast<std::uint16_t>(i));
    }
  }
  double sum = 0.0;
  for (const std::uint16_t w : indices) sum += system.wavelength_weights()[w];
  if (!(sum > 0.0)) throw std::invalid_argument("analysis: the chosen wavelengths have weight 0");
  std::vector<std::pair<std::uint16_t, double>> out;
  out.reserve(indices.size());
  for (const std::uint16_t w : indices) out.emplace_back(w, system.wavelength_weights()[w] / sum);
  return out;
}

/// Sums over the cells of a pupil function (ADR 0033, point 4).
struct PupilSums {
  double fa = 0.0;    ///< sum f A
  double fa2 = 0.0;   ///< E = sum f A^2, the exact pupil energy
  double f2a2 = 0.0;  ///< sum f^2 A^2 = sum |u|^2
};

PupilSums pupil_sums(const PupilFunction& p) {
  PupilSums s;
  for (std::size_t j = 0; j < p.fraction.size(); ++j) {
    const double fa = p.fraction[j] * p.amplitude[j];
    s.fa += fa;
    s.fa2 += fa * p.amplitude[j];
    s.f2a2 += fa * fa;
  }
  return s;
}

/// Cell value u_j = f A exp(-2 pi i W lambda_ref / lambda) (ADR 0032, point 3;
/// ADR 0033, 2).
Complex cell_value(const PupilFunction& p, std::size_t j) {
  const double fa = p.fraction[j] * p.amplitude[j];
  if (fa == 0.0) return {0.0, 0.0};
  const double phase =
      -2.0 * std::numbers::pi * p.w[j] * p.reference_wavelength_um / p.wavelength_um;
  return std::polar(fa, phase);
}

/// The checks every evaluation of a pupil function needs: a mapping and an
/// arrived cell.
PupilSums require_evaluable(const PupilFunction& p) {
  if (!p.mapping) {
    throw AnalysisError(
        "analysis: the image space is afocal (no mapping of "
        "the pupil onto the image plane); "
        "no PSF",
        path_pointer(p.path));
  }
  const PupilSums s = pupil_sums(p);
  if (!(s.fa > 0.0)) {
    throw AnalysisError("analysis: no cell of the pupil function arrived; no PSF",
                        path_pointer(p.path));
  }
  return s;
}

/// The mapping of a pupil function that require_evaluable() accepted (checked again, so that
/// every access is visibly checked).
const PupilMapping& mapping_of(const PupilFunction& p) {
  if (!p.mapping) throw std::logic_error("analysis: PSF without the checked mapping");
  return *p.mapping;
}

/// sum_j u_j in cell order: U at x = 0 up to the factor |det M| a^2 (ADR 0033, point 2).
Complex sum_of_cells(const PupilFunction& p) {
  Complex sum{0.0, 0.0};
  for (std::size_t j = 0; j < p.fraction.size(); ++j) sum += cell_value(p, j);
  return sum;
}

/// |sum u|^2 / |sum f A|^2: the Strehl ratio of the discrete pupil (ADR 0033,
/// point 4; Wyant & Creath Eq. (64)), the intensity at C' relative to the pupil
/// without aberration.
double strehl_of(const PupilFunction& p, const PupilSums& s) {
  return std::norm(sum_of_cells(p)) / (s.fa * s.fa);
}

/// Centre coordinate (2 i + 1 - N) / N of cell column or row i (pupil.hpp): the
/// lattice of the FFT (ADR 0033, point 3). W and A were sampled at the centroid
/// of a boundary cell (ADR 0032, point 1); the transform places every cell at
/// its centre.
double cell_centre(int i, int n) {
  return static_cast<double>(2 * i + 1 - n) / static_cast<double>(n);
}

/// Largest phase step |Delta(W lambda_ref / lambda)| between horizontally or
/// vertically neighbouring arrived cells with f > 0, waves at lambda (ADR 0033,
/// point 3).
double largest_phase_step(const PupilFunction& p) {
  const auto n = static_cast<std::size_t>(p.grid);
  const double scale = p.reference_wavelength_um / p.wavelength_um;
  const auto usable = [&](std::size_t j) {
    return p.fraction[j] > 0.0 && p.status[j] == trace::RayStatus::Alive;
  };
  double largest = 0.0;
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t j = k * n + i;
      if (!usable(j)) continue;
      if (i + 1 < n && usable(j + 1)) {
        largest = std::max(largest, std::abs(p.w[j + 1] - p.w[j]) * scale);
      }
      if (k + 1 < n && usable(j + n)) {
        largest = std::max(largest, std::abs(p.w[j + n] - p.w[j]) * scale);
      }
    }
  }
  return largest;
}

/// The pupil in the top left N x N corner of an L x L array (rows py, columns
/// px), zero padded, transformed forward: U at the lattice points up to a phase
/// (ADR 0033, point 3), index (k2 mod L) * L + (k1 mod L).
std::vector<Complex> transformed(const PupilFunction& p, std::size_t length, bool ideal) {
  const auto n = static_cast<std::size_t>(p.grid);
  std::vector<Complex> data(length * length, Complex(0.0, 0.0));
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t j = k * n + i;
      data[k * length + i] =
          ideal ? Complex(p.fraction[j] * p.amplitude[j], 0.0) : cell_value(p, j);
    }
  }
  detail::fft2(data, length, length, detail::FftDirection::Forward);
  return data;
}

/// Fraction of sum |U|^2 in the outer max(1, L / 16) rows and columns of the
/// centred lattice (ADR 0033, point 3, psf.aliasing).
double border_fraction(const std::vector<Complex>& u, std::size_t length) {
  const std::size_t border = std::max<std::size_t>(1, length / 16);
  const std::size_t half = length / 2;
  const auto outer = [&](std::size_t idx) {
    const std::size_t s = (idx + half) % length;  // position in the centred lattice
    return s < border || s >= length - border;
  };
  double all = 0.0;
  double edge = 0.0;
  for (std::size_t r = 0; r < length; ++r) {
    for (std::size_t c = 0; c < length; ++c) {
      const double v = std::norm(u[r * length + c]);
      all += v;
      if (outer(r) || outer(c)) edge += v;
    }
  }
  return all > 0.0 ? edge / all : 0.0;
}

/// Lattice matrix B = lambda / (L a) M^-T, row-major, mm (ADR 0033, point 3).
std::array<double, 4> lattice_matrix(const PupilFunction& p, std::size_t length) {
  const auto& m = mapping_of(p).m;
  const double det = m[0] * m[3] - m[1] * m[2];
  const double a = 2.0 / static_cast<double>(p.grid);
  const double scale = p.wavelength_um * 1e-3 / (static_cast<double>(length) * a * det);
  // M^-T = (1 / det) [m3, -m2; -m1, m0]
  return {scale * m[3], -scale * m[2], -scale * m[1], scale * m[0]};
}

double frobenius(const std::array<double, 4>& b) {
  return std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2] + b[3] * b[3]);
}

/// One wavelength of a PSF in the making.
struct Part {
  PupilFunction pupil;
  PupilSums sums;
  double weight = 0.0;
  std::size_t length = 0;
};

void add_warning(std::vector<model::Diagnostic>& out,
                 diagnostics::DiagnosticCode code,
                 PathId path,
                 std::string message) {
  out.push_back({code.severity(), path_pointer(path), std::move(message), std::string(code.str())});
}

}  // namespace

Psf psf(const compile::CompiledSystem& system,
        compile::PathId path,
        std::uint16_t field,
        const PsfOptions& options) {
  return psf(system, path, field, options, trace::RunControl{});
}

Psf psf(const compile::CompiledSystem& system,
        compile::PathId path,
        std::uint16_t field,
        const PsfOptions& options,
        const trace::RunControl& control) {
  detail::check_path(system, path);
  check_options(options);
  const auto chosen = chosen_wavelengths(system, options.wavelengths);
  const int n = options.pupil.grid;
  // The memory guard before any ray: L = Q N of the reference (the larger L_i
  // of longer wavelengths are checked below once M is known).
  if (static_cast<long long>(options.padding) * n > kMaxFftLength) {
    throw std::invalid_argument("analysis: PSF length L = padding * grid exceeds " +
                                std::to_string(kMaxFftLength));
  }

  // Pupil functions, in the order of the choice (ADR 0033, point 7).
  std::vector<Part> parts;
  parts.reserve(chosen.size());
  for (const auto& [w, weight] : chosen) {
    Part part;
    part.pupil = pupil_function(system, path, field, w, options.pupil, control);
    part.sums = require_evaluable(part.pupil);
    part.weight = weight;
    parts.push_back(std::move(part));
  }
  // The lattice of the reference: the system reference wavelength if chosen,
  // else the first.
  std::size_t ref = 0;
  for (std::size_t i = 0; i < parts.size(); ++i) {
    if (parts[i].pupil.wavelength == system.reference_wavelength()) ref = i;
  }
  const PupilFunction& p_ref = parts[ref].pupil;
  const auto det_of = [](const PupilFunction& p) {
    const auto& m = mapping_of(p).m;
    return std::abs(m[0] * m[3] - m[1] * m[2]);
  };
  // L_i = round(L_ref (lambda_i / lambda_ref) sqrt(|det M_ref| / |det M_i|)),
  // to the nearest integer (ADR 0033, point 7 with the addendum of #205: the
  // bound 1 / (2 L_i) of the scale error needs the nearest integer); L_ref
  // grows in steps of 2 until every L_i >= 2 N - 1.
  std::size_t l_ref = static_cast<std::size_t>(options.padding) * static_cast<std::size_t>(n);
  const auto length_of = [&](const Part& part, std::size_t l) {
    const double ratio = part.pupil.wavelength_um / p_ref.wavelength_um *
                         std::sqrt(det_of(p_ref) / det_of(part.pupil));
    return static_cast<std::size_t>(std::llround(static_cast<double>(l) * ratio));
  };
  const auto min_length = static_cast<std::size_t>(2 * n - 1);
  for (;;) {
    bool enough = true;
    for (const Part& part : parts) enough = enough && length_of(part, l_ref) >= min_length;
    if (enough) break;
    l_ref += 2;
  }
  for (Part& part : parts) {
    part.length = length_of(part, l_ref);
    if (part.length > static_cast<std::size_t>(kMaxFftLength)) {
      throw std::invalid_argument("analysis: PSF length L_i = " + std::to_string(part.length) +
                                  " of wavelength " + std::to_string(part.pupil.wavelength) +
                                  " exceeds " + std::to_string(kMaxFftLength));
    }
  }
  parts[ref].length = l_ref;

  Psf result;
  result.field = field;
  result.grid = n;
  result.length = static_cast<int>(l_ref);
  result.centre = mapping_of(p_ref).centre;
  result.e_x = mapping_of(p_ref).e_x;
  result.e_y = mapping_of(p_ref).e_y;
  result.lattice = lattice_matrix(p_ref, l_ref);
  result.intensity.assign(l_ref * l_ref, 0.0);

  double sum_we = 0.0;
  for (const Part& part : parts) sum_we += part.weight * part.sums.fa2;
  const auto half_ref = static_cast<long long>(l_ref / 2);
  double worst_aliasing = 0.0;
  std::optional<std::uint16_t> aliasing_wavelength;
  double worst_step = 0.0;
  std::optional<std::uint16_t> step_wavelength;
  // Polychromatic Strehl ratio as the peak ratio of the summed images (review of #233): the
  // intensity at k = 0, sum omega_i |sum u|_i^2 / (L_i^2 E_i), over that of the pupils without
  // aberration, sum omega_i (sum f A)_i^2 / (L_i^2 E_i), omega_i = w_i E_i / sum w E.
  double peak_sum = 0.0;
  double ideal_sum = 0.0;
  for (const Part& part : parts) {
    const PupilFunction& p = part.pupil;
    const std::size_t l = part.length;
    const double l2 = static_cast<double>(l) * static_cast<double>(l);
    double outside = 0.0;
    double edge = 0.0;
    {
      const std::vector<Complex> u = transformed(p, l, false);
      // ADR 0033, point 4: I_k = |U_k|^2 / (L^2 E), weighted with w_i E_i / sum w E.
      const double scale = part.weight * part.sums.fa2 / sum_we / (l2 * part.sums.fa2);
      // Every cell of the own lattice, k in [-(L_i / 2), L_i - L_i / 2): inside the reference
      // field it is added, outside (only for L_i > L) its energy is reported as outside_field.
      const auto li = static_cast<long long>(l);
      const auto centred = [li](std::size_t idx) {
        const auto k = static_cast<long long>(idx);
        return k < li - li / 2 ? k : k - li;
      };
      const auto inside = [half_ref, l_ref](long long k) {
        return k >= -half_ref && k < static_cast<long long>(l_ref) - half_ref;
      };
      const double unit = 1.0 / (l2 * part.sums.fa2);
      for (std::size_t row = 0; row < l; ++row) {
        const long long k2 = centred(row);
        for (std::size_t col = 0; col < l; ++col) {
          const long long k1 = centred(col);
          const double v = std::norm(u[row * l + col]);
          if (inside(k1) && inside(k2)) {
            const auto target = static_cast<std::size_t>(k2 + half_ref) * l_ref +
                                static_cast<std::size_t>(k1 + half_ref);
            result.intensity[target] += scale * v;
          } else {
            outside += unit * v;
          }
        }
      }
      edge = border_fraction(u, l);
    }  // the pupil's transform is freed before the ideal one: one L^2 array at a time
    // Warnings (ADR 0033, point 3).
    const double edge_ideal = border_fraction(transformed(p, l, true), l);
    if (edge > options.aliasing_factor * edge_ideal) {
      const double excess = edge_ideal > 0.0 ? edge / edge_ideal : edge;
      if (!aliasing_wavelength || excess > worst_aliasing) {
        worst_aliasing = excess;
        aliasing_wavelength = p.wavelength;
      }
    }
    const double step = largest_phase_step(p);
    if (step > options.undersampled_waves && (!step_wavelength || step > worst_step)) {
      worst_step = step;
      step_wavelength = p.wavelength;
    }
    // Weighted scalars.
    const double we = part.weight * part.sums.fa2 / sum_we;
    peak_sum += we * std::norm(sum_of_cells(p)) / (l2 * part.sums.fa2);
    ideal_sum += we * part.sums.fa * part.sums.fa / (l2 * part.sums.fa2);
    result.edge_energy_deficit += we * (1.0 - part.sums.f2a2 / part.sums.fa2);
    if (&part != &parts[ref]) {
      const std::array<double, 4> b = lattice_matrix(p, l);
      std::array<double, 4> d{};
      for (std::size_t c = 0; c < 4; ++c) d[c] = b[c] - result.lattice[c];
      result.polychromatic_scale_error =
          std::max(result.polychromatic_scale_error, frobenius(d) / frobenius(result.lattice));
    }
    result.wavelengths.push_back(
        {p.wavelength, part.weight, part.sums.fa2, static_cast<int>(l), outside});
  }

  // Monochromatic exactly |sum u|^2 / |sum f A|^2, bit for bit as through_focus().
  result.strehl =
      parts.size() == 1 ? strehl_of(parts[0].pupil, parts[0].sums) : peak_sum / ideal_sum;

  // Peak of the lattice, the first maximum in row-major order.
  std::size_t best = 0;
  for (std::size_t c = 1; c < result.intensity.size(); ++c) {
    if (result.intensity[c] > result.intensity[best]) best = c;
  }
  result.peak.value = result.intensity[best];
  result.peak.k1 = static_cast<int>(static_cast<long long>(best % l_ref) - half_ref);
  result.peak.k2 = static_cast<int>(static_cast<long long>(best / l_ref) - half_ref);
  result.peak.x = result.lattice[0] * result.peak.k1 + result.lattice[1] * result.peak.k2;
  result.peak.y = result.lattice[2] * result.peak.k1 + result.lattice[3] * result.peak.k2;

  // Losses and warnings over every traced cell of every wavelength, counted as
  // the pupil functions count them (f > 0, final status), so one wavelength
  // gives its own losses.
  detail::LossCounter losses(system, path, options.pupil.lost_warning_fraction);
  const model::Diagnostic* mapping_warning = nullptr;
  double mapping_delta = -1.0;
  for (const Part& part : parts) {
    const PupilFunction& p = part.pupil;
    for (std::size_t j = 0; j < p.fraction.size(); ++j) {
      if (p.fraction[j] > 0.0) losses.add_ray(p.status[j], p.end.last_surface[j]);
    }
    for (const model::Diagnostic& d : p.warnings) {
      if (d.code == "pupil.mapping_nonlinear" && mapping_of(p).delta_max > mapping_delta) {
        mapping_warning = &d;
        mapping_delta = mapping_of(p).delta_max;
      }
    }
  }
  result.losses = losses.result();
  result.warnings = losses.warnings();
  if (mapping_warning) result.warnings.push_back(*mapping_warning);
  if (aliasing_wavelength) {
    constexpr diagnostics::DiagnosticCode kCode = "psf.aliasing";
    add_warning(result.warnings, kCode, path,
                "the PSF is wider than its field and folds in (wavelength " +
                    std::to_string(*aliasing_wavelength) +
                    "): the energy in the outer L/16 rows "
                    "and columns is " +
                    std::to_string(worst_aliasing) +
                    " times that of the pupil without aberration; increase the grid N "
                    "(ADR 0033, point 3)");
  }
  if (step_wavelength) {
    constexpr diagnostics::DiagnosticCode kCode = "psf.pupil_undersampled";
    add_warning(result.warnings, kCode, path,
                "the pupil grid is too coarse for the aberration (wavelength " +
                    std::to_string(*step_wavelength) + "): neighbouring cells differ by " +
                    std::to_string(worst_step) + " waves, more than " +
                    std::to_string(options.undersampled_waves) + " (ADR 0033, point 3)");
  }
  return result;
}

std::vector<double> psf_at(const PupilFunction& pupil, std::span<const Point2> x) {
  const PupilSums s = require_evaluable(pupil);
  const auto& m = mapping_of(pupil).m;
  const double lambda_mm = pupil.wavelength_um * 1e-3;
  const int n = pupil.grid;
  // Cell values and their directions s_j = M c_j at the cell centres c_j (the
  // FFT lattice).
  std::vector<Complex> u;
  std::vector<double> sx;
  std::vector<double> sy;
  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < n; ++i) {
      const auto j =
          static_cast<std::size_t>(k) * static_cast<std::size_t>(n) + static_cast<std::size_t>(i);
      const Complex v = cell_value(pupil, j);
      if (v == Complex(0.0, 0.0)) continue;
      const double cx = cell_centre(i, n);
      const double cy = cell_centre(k, n);
      u.push_back(v);
      sx.push_back(m[0] * cx + m[1] * cy);
      sy.push_back(m[2] * cx + m[3] * cy);
    }
  }
  std::vector<double> out;
  out.reserve(x.size());
  for (const Point2& point : x) {
    Complex sum{0.0, 0.0};
    for (std::size_t j = 0; j < u.size(); ++j) {
      const double phase =
          -2.0 * std::numbers::pi * (point.x * sx[j] + point.y * sy[j]) / lambda_mm;
      sum += u[j] * std::polar(1.0, phase);
    }
    out.push_back(std::norm(sum) / (s.fa * s.fa));
  }
  return out;
}

ThroughFocus through_focus(const compile::CompiledSystem& system,
                           compile::PathId path,
                           std::uint16_t field,
                           std::uint16_t wavelength,
                           std::span<const double> shifts,
                           const PupilFunctionOptions& options) {
  return through_focus(system, path, field, wavelength, shifts, options, trace::RunControl{});
}

ThroughFocus through_focus(const compile::CompiledSystem& system,
                           compile::PathId path,
                           std::uint16_t field,
                           std::uint16_t wavelength,
                           std::span<const double> shifts,
                           const PupilFunctionOptions& options,
                           const trace::RunControl& control) {
  for (const double d : shifts) {
    if (!std::isfinite(d)) throw std::invalid_argument("analysis: focus shifts must be finite");
  }
  PupilFunctionOptions base = options;
  base.focus_shift = 0.0;
  const PupilFunction p0 = pupil_function(system, path, field, wavelength, base, control);
  // Afocal as in psf() (review of #233): the reference of an afocal image space is not settled.
  static_cast<void>(require_evaluable(p0));
  ThroughFocus result;
  result.field = field;
  result.wavelength = wavelength;
  result.shifts.assign(shifts.begin(), shifts.end());
  result.strehl.reserve(shifts.size());
  result.arrived_fraction.reserve(shifts.size());
  result.losses = p0.losses;
  result.warnings = p0.warnings;
  for (const double d : shifts) {
    const PupilFunction p = refocus(system, p0, d);
    const PupilSums s = pupil_sums(p);
    if (!(s.fa > 0.0)) {
      throw AnalysisError(
          "analysis: no cell of the pupil function arrived at focus shift " + std::to_string(d),
          path_pointer(path));
    }
    result.strehl.push_back(strehl_of(p, s));
    result.arrived_fraction.push_back(p.arrived_fraction);
  }
  return result;
}

// STUB for the red run of #206 (tests first): MTF and encircled energy link, compute nothing.

Mtf mtf(const compile::CompiledSystem& /*system*/,
        compile::PathId /*path*/,
        std::uint16_t /*field*/,
        const MtfOptions& /*options*/) {
  return {};
}

Mtf mtf(const compile::CompiledSystem& /*system*/,
        compile::PathId /*path*/,
        std::uint16_t /*field*/,
        const MtfOptions& /*options*/,
        const trace::RunControl& /*control*/) {
  return {};
}

Mtf mtf(const PupilFunction& /*pupil*/, const MtfOutputs& /*outputs*/) {
  return {};
}

EncircledEnergy encircled_energy(const compile::CompiledSystem& /*system*/,
                                 compile::PathId /*path*/,
                                 std::uint16_t /*field*/,
                                 std::span<const double> /*radii*/,
                                 const EncircledEnergyOptions& /*options*/) {
  return {};
}

EncircledEnergy encircled_energy(const compile::CompiledSystem& /*system*/,
                                 compile::PathId /*path*/,
                                 std::uint16_t /*field*/,
                                 std::span<const double> /*radii*/,
                                 const EncircledEnergyOptions& /*options*/,
                                 const trace::RunControl& /*control*/) {
  return {};
}

EncircledEnergy encircled_energy(const PupilFunction& /*pupil*/,
                                 std::span<const double> /*radii*/,
                                 EeCentre /*centre*/) {
  return {};
}

}  // namespace rtt::analysis
