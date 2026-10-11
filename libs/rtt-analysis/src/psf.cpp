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

// ---------------------------------------------------------------------------------------------
// MTF and encircled energy (ADR 0033, points 5-7, with the erratum on the signs; #206).
// ---------------------------------------------------------------------------------------------

namespace {

/// Autocorrelation C(l) = sum_j u_(j+l) conj(u_j) of the cell values for every shift l (l1 along
/// px, l2 along py), by FFT with L = 2 N >= 2 N - 1, so that no shift wraps around (ADR 0033,
/// point 5): Backward(|Forward(u)|^2) = L^2 C (derived from the two DFT sums).
class Autocorrelation {
 public:
  Autocorrelation() = default;

  explicit Autocorrelation(const PupilFunction& p)
      : n_(p.grid), length_(2 * static_cast<std::size_t>(p.grid)) {
    const auto n = static_cast<std::size_t>(p.grid);
    c_.assign(length_ * length_, Complex(0.0, 0.0));
    for (std::size_t k = 0; k < n; ++k) {
      for (std::size_t i = 0; i < n; ++i) c_[k * length_ + i] = cell_value(p, k * n + i);
    }
    detail::fft2(c_, length_, length_, detail::FftDirection::Forward);
    for (Complex& v : c_) v = Complex(std::norm(v), 0.0);
    detail::fft2(c_, length_, length_, detail::FftDirection::Backward);
    const double scale = 1.0 / (static_cast<double>(length_) * static_cast<double>(length_));
    for (Complex& v : c_) v *= scale;
  }

  /// C(l1, l2); exactly 0 for |l1| >= N or |l2| >= N, where the shifted pupil has no overlap.
  [[nodiscard]] Complex at(int l1, int l2) const {
    if (std::abs(l1) >= n_ || std::abs(l2) >= n_) return {0.0, 0.0};
    const auto l = static_cast<long long>(length_);
    const auto row = static_cast<std::size_t>(((l2 % l) + l) % l);
    const auto col = static_cast<std::size_t>(((l1 % l) + l) % l);
    return c_[row * length_ + col];
  }

 private:
  int n_ = 0;
  std::size_t length_ = 0;
  std::vector<Complex> c_;
};

/// One wavelength of an MTF or EE in the making.
struct Spectral {
  PupilFunction pupil;
  PupilSums sums;
  double weight = 0.0;  ///< normalised model weight w_i
  Autocorrelation c;    ///< of the pupil function, set with it
};

/// F = M a / lambda, row-major, cycles per mm per cell (ADR 0033, point 5).
std::array<double, 4> frequency_matrix_of(const PupilFunction& p) {
  const auto& m = mapping_of(p).m;
  const double scale = 2.0 / static_cast<double>(p.grid) / (p.wavelength_um * 1e-3);
  return {scale * m[0], scale * m[1], scale * m[2], scale * m[3]};
}

/// OTF(nu(l)) = conj(C(l)) / E, with C(0) := E (ADR 0033, point 5 and the erratum): 1 at l = 0.
Complex otf_at(const Autocorrelation& c, const PupilSums& s, int l1, int l2) {
  if (l1 == 0 && l2 == 0) return {1.0, 0.0};
  return std::conj(c.at(l1, l2)) / s.fa2;
}

/// The curve along column `col` of F (0: sagittal, l = (l, 0); 1: tangential, l = (0, l)),
/// samples l = 0 ... last (last >= N; from N on the OTF is 0), with the complex OTF.
MtfCurve curve_of(const Autocorrelation& c,
                  const PupilSums& s,
                  const std::array<double, 4>& f,
                  int last,
                  int col) {
  MtfCurve curve;
  const double fx = col == 0 ? f[0] : f[1];
  const double fy = col == 0 ? f[2] : f[3];
  const double step = std::hypot(fx, fy);
  curve.direction = {fx / step, fy / step};
  for (int l = 0; l <= last; ++l) {
    const Complex otf = col == 0 ? otf_at(c, s, l, 0) : otf_at(c, s, 0, l);
    curve.frequency.push_back(static_cast<double>(l) * step);
    curve.otf.push_back(otf);
    curve.modulation.push_back(std::abs(otf));
  }
  return curve;
}

/// The complex OTF of `curve` at |nu|, linearly interpolated (real and imaginary part); 0 at and
/// above the cut-off, the last sample (ADR 0033, point 7).
Complex interpolated(const MtfCurve& curve, double nu) {
  if (nu >= curve.frequency.back()) return {0.0, 0.0};
  const auto hi = std::upper_bound(curve.frequency.begin(), curve.frequency.end(), nu);
  const auto k = static_cast<std::size_t>(hi - curve.frequency.begin());
  const double g = (nu - curve.frequency[k - 1]) / (curve.frequency[k] - curve.frequency[k - 1]);
  return (1.0 - g) * curve.otf[k - 1] + g * curve.otf[k];
}

/// The part of the choice that serves as reference: the system reference wavelength if chosen,
/// otherwise the first (as psf()).
std::size_t reference_part(const CompiledSystem& system, const std::vector<Spectral>& parts) {
  std::size_t ref = 0;
  for (std::size_t i = 0; i < parts.size(); ++i) {
    if (parts[i].pupil.wavelength == system.reference_wavelength()) ref = i;
  }
  return ref;
}

/// The pupil functions of the chosen wavelengths, each checked (mapping, an arrived cell).
std::vector<Spectral> spectral_parts(const CompiledSystem& system,
                                     PathId path,
                                     std::uint16_t field,
                                     const std::optional<std::vector<std::uint16_t>>& wavelengths,
                                     const PupilFunctionOptions& pupil_options,
                                     const trace::RunControl& control) {
  const auto chosen = chosen_wavelengths(system, wavelengths);
  std::vector<Spectral> parts;
  parts.reserve(chosen.size());
  for (const auto& [w, weight] : chosen) {
    Spectral part;
    part.pupil = pupil_function(system, path, field, w, pupil_options, control);
    part.sums = require_evaluable(part.pupil);
    part.weight = weight;
    part.c = Autocorrelation(part.pupil);
    parts.push_back(std::move(part));
  }
  return parts;
}

/// Losses and warnings over the pupil functions, counted as they count them, plus
/// pupil.mapping_nonlinear (the largest delta) and psf.pupil_undersampled (as psf()).
void losses_and_warnings(const CompiledSystem& system,
                         PathId path,
                         const std::vector<Spectral>& parts,
                         double lost_warning_fraction,
                         double undersampled_waves,
                         RayLosses& losses_out,
                         std::vector<model::Diagnostic>& warnings_out) {
  detail::LossCounter losses(system, path, lost_warning_fraction);
  const model::Diagnostic* mapping_warning = nullptr;
  double mapping_delta = -1.0;
  double worst_step = 0.0;
  std::optional<std::uint16_t> step_wavelength;
  for (const Spectral& part : parts) {
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
    const double step = largest_phase_step(p);
    if (step > undersampled_waves && (!step_wavelength || step > worst_step)) {
      worst_step = step;
      step_wavelength = p.wavelength;
    }
  }
  losses_out = losses.result();
  warnings_out = losses.warnings();
  if (mapping_warning) warnings_out.push_back(*mapping_warning);
  if (step_wavelength) {
    constexpr diagnostics::DiagnosticCode kCode = "psf.pupil_undersampled";
    add_warning(warnings_out, kCode, path,
                "the pupil grid is too coarse for the aberration (wavelength " +
                    std::to_string(*step_wavelength) + "): neighbouring cells differ by " +
                    std::to_string(worst_step) + " waves, more than " +
                    std::to_string(undersampled_waves) + " (ADR 0033, point 3)");
  }
}

void check_undersampled(double waves) {
  if (!(waves > 0.0) || !std::isfinite(waves)) {
    throw std::invalid_argument("analysis: undersampled_waves must be finite and > 0");
  }
}

/// The MTF of the parts (one or several wavelengths), without losses and warnings.
Mtf mtf_of(const std::vector<Spectral>& parts, std::size_t ref, const MtfOutputs& outputs) {
  const Spectral& r = parts[ref];
  const PupilFunction& pr = r.pupil;
  const int n = pr.grid;
  Mtf result;
  result.field = pr.field;
  result.grid = n;
  result.centre = mapping_of(pr).centre;
  result.e_x = mapping_of(pr).e_x;
  result.e_y = mapping_of(pr).e_y;
  result.frequency_matrix = frequency_matrix_of(pr);
  double sum_we = 0.0;
  for (const Spectral& part : parts) sum_we += part.weight * part.sums.fa2;
  if (parts.size() == 1) {
    result.sagittal = curve_of(r.c, r.sums, result.frequency_matrix, n, 0);
    result.tangential = curve_of(r.c, r.sums, result.frequency_matrix, n, 1);
  } else {
    // ADR 0033, point 7: the complex OTF of each wavelength on its own frequencies, interpolated
    // linearly in |nu| onto those of the reference, summed with w_i E_i. The axis runs in the
    // steps of the reference up to the first sample at or above the largest cut-off of all
    // wavelengths (a shorter wavelength reaches beyond that of the reference; review of #239),
    // so the last value is 0 for every wavelength. The direction of the reference is taken for
    // all: exact where M_i e is parallel to M_ref e (centred systems with the field on +-y).
    for (int col = 0; col < 2; ++col) {
      const auto step_of = [col](const std::array<double, 4>& f) {
        return col == 0 ? std::hypot(f[0], f[2]) : std::hypot(f[1], f[3]);
      };
      const double step = step_of(result.frequency_matrix);
      double cutoff = static_cast<double>(n) * step;
      for (const Spectral& part : parts) {
        cutoff =
            std::max(cutoff, static_cast<double>(n) * step_of(frequency_matrix_of(part.pupil)));
      }
      int last = n;
      while (static_cast<double>(last) * step < cutoff) ++last;
      MtfCurve sum = curve_of(r.c, r.sums, result.frequency_matrix, last, col);
      std::vector<Complex> acc(sum.otf.size(), Complex(0.0, 0.0));
      for (const Spectral& part : parts) {
        const double we = part.weight * part.sums.fa2;
        if (&part == &r) {
          for (std::size_t l = 0; l < acc.size(); ++l) acc[l] += we * sum.otf[l];
          continue;
        }
        const MtfCurve own = curve_of(part.c, part.sums, frequency_matrix_of(part.pupil), n, col);
        for (std::size_t l = 0; l < acc.size(); ++l) {
          acc[l] += we * interpolated(own, sum.frequency[l]);
        }
      }
      for (std::size_t l = 0; l < acc.size(); ++l) {
        sum.otf[l] = acc[l] / sum_we;
        sum.modulation[l] = std::abs(sum.otf[l]);
      }
      (col == 0 ? result.sagittal : result.tangential) = std::move(sum);
    }
  }
  if (outputs.map) {  // monochromatic only (checked by the callers)
    const std::size_t side = 2 * static_cast<std::size_t>(n) - 1;
    result.map.reserve(side * side);
    if (outputs.complex_otf) result.otf_map.reserve(side * side);
    for (int l2 = -(n - 1); l2 <= n - 1; ++l2) {
      for (int l1 = -(n - 1); l1 <= n - 1; ++l1) {
        const Complex otf = otf_at(r.c, r.sums, l1, l2);
        result.map.push_back(std::abs(otf));
        if (outputs.complex_otf) result.otf_map.push_back(otf);
      }
    }
  }
  if (!outputs.complex_otf) {
    result.sagittal.otf.clear();
    result.tangential.otf.clear();
  }
  for (const Spectral& part : parts) {
    result.edge_energy_deficit +=
        part.weight * part.sums.fa2 / sum_we * (1.0 - part.sums.f2a2 / part.sums.fa2);
    result.wavelengths.push_back({part.pupil.wavelength, part.weight, part.sums.fa2, 2 * n, 0.0});
  }
  return result;
}

/// Largest EE radius of one wavelength: a quarter of the shortest period lambda / a |M^-T e_i|
/// of the periodic PSF (ADR 0033, point 6).
double radius_max_of(const PupilFunction& p) {
  const std::array<double, 4> b = lattice_matrix(p, 1);  // lambda / a M^-T
  return 0.25 * std::min(std::hypot(b[0], b[2]), std::hypot(b[1], b[3]));
}

/// EeCentre::LinearPhase: x0 with (a / lambda) M^T x0 = (arg C(e1), arg C(e2)) / (2 pi), i.e.
/// x0 = (lambda / a) M^-T (arg C(e1), arg C(e2)) / (2 pi); exact for a moved PSF, since
/// I(x - t) has C(l) exp(+2 pi i nu(l) . t) (ADR 0033, erratum).
Point2 linear_phase_centre(const Spectral& part) {
  const std::array<double, 4> b = lattice_matrix(part.pupil, 1);
  const double phi1 = std::arg(part.c.at(1, 0)) / (2.0 * std::numbers::pi);
  const double phi2 = std::arg(part.c.at(0, 1)) / (2.0 * std::numbers::pi);
  return {b[0] * phi1 + b[1] * phi2, b[2] * phi1 + b[3] * phi2};
}

/// EE(r) of one wavelength about x0 (ADR 0033, point 6 with the erratum):
/// (|det M| a^2 / lambda^2) sum_l Re[C(l) exp(-2 pi i nu(l) . x0)] K_r(nu(l)) / E with
/// K_r(nu) = r J1(2 pi r |nu|) / |nu|, K_r(0) = pi r^2, C(0) = sum |u|^2; serial in grid order.
std::vector<double> energy_of(const Spectral& part,
                              const Point2& x0,
                              std::span<const double> radii) {
  const PupilFunction& p = part.pupil;
  const int n = p.grid;
  const auto& m = mapping_of(p).m;
  const double a = 2.0 / static_cast<double>(n);
  const double lambda = p.wavelength_um * 1e-3;
  const double scale = std::abs(m[0] * m[3] - m[1] * m[2]) * a * a / (lambda * lambda);
  const std::array<double, 4> f = frequency_matrix_of(p);
  std::vector<double> sums(radii.size(), 0.0);
  for (int l2 = -(n - 1); l2 <= n - 1; ++l2) {
    for (int l1 = -(n - 1); l1 <= n - 1; ++l1) {
      const double nx = f[0] * l1 + f[1] * l2;
      const double ny = f[2] * l1 + f[3] * l2;
      const double nu = std::hypot(nx, ny);
      double re = 0.0;
      if (l1 == 0 && l2 == 0) {
        re = part.sums.f2a2;  // C(0) = sum |u|^2, the discrete sum (ADR 0033, point 6)
      } else {
        const Complex shift = std::polar(1.0, -2.0 * std::numbers::pi * (nx * x0.x + ny * x0.y));
        re = (part.c.at(l1, l2) * shift).real();
      }
      for (std::size_t i = 0; i < radii.size(); ++i) {
        const double r = radii[i];
        const double kernel =
            nu == 0.0 ? std::numbers::pi * r * r
                      : r * std::cyl_bessel_j(1.0, 2.0 * std::numbers::pi * r * nu) / nu;
        sums[i] += re * kernel;
      }
    }
  }
  for (double& s : sums) s *= scale / part.sums.fa2;
  return sums;
}

void check_radii(std::span<const double> radii) {
  for (const double r : radii) {
    if (!std::isfinite(r) || r < 0.0) {
      throw std::invalid_argument("analysis: encircled-energy radii must be finite and >= 0");
    }
  }
}

/// EE of the parts about EeCentre `centre` (one common centre for all wavelengths).
EncircledEnergy energy_of_parts(const std::vector<Spectral>& parts,
                                std::size_t ref,
                                std::span<const double> radii,
                                EeCentre centre) {
  const PupilFunction& pr = parts[ref].pupil;
  EncircledEnergy result;
  result.field = pr.field;
  result.grid = pr.grid;
  result.e_x = mapping_of(pr).e_x;
  result.e_y = mapping_of(pr).e_y;
  result.radii.assign(radii.begin(), radii.end());
  double sum_we = 0.0;
  for (const Spectral& part : parts) sum_we += part.weight * part.sums.fa2;
  result.radius_max = std::numeric_limits<double>::infinity();
  for (const Spectral& part : parts) {
    result.radius_max = std::min(result.radius_max, radius_max_of(part.pupil));
  }
  for (const double r : radii) {
    if (r > result.radius_max) {
      throw std::invalid_argument("analysis: encircled-energy radius " + std::to_string(r) +
                                  " mm above radius_max " + std::to_string(result.radius_max) +
                                  " mm (a quarter of the period of the PSF; ADR 0033, point 6)");
    }
  }
  if (centre == EeCentre::LinearPhase) {
    for (const Spectral& part : parts) {
      const double we = part.weight * part.sums.fa2;
      const Point2 x0 = linear_phase_centre(part);
      result.offset.x += we * x0.x;
      result.offset.y += we * x0.y;
    }
    result.offset.x /= sum_we;
    result.offset.y /= sum_we;
  }
  result.centre =
      mapping_of(pr).centre + result.offset.x * result.e_x + result.offset.y * result.e_y;
  result.energy.assign(radii.size(), 0.0);
  for (const Spectral& part : parts) {
    const double we = part.weight * part.sums.fa2;
    const std::vector<double> e = energy_of(part, result.offset, radii);
    for (std::size_t i = 0; i < radii.size(); ++i) result.energy[i] += we * e[i];
    result.edge_energy_deficit += we / sum_we * (1.0 - part.sums.f2a2 / part.sums.fa2);
    result.wavelengths.push_back(
        {part.pupil.wavelength, part.weight, part.sums.fa2, 2 * part.pupil.grid, 0.0});
  }
  for (double& e : result.energy) e /= sum_we;
  return result;
}

/// The single part of a pupil function given directly (weight 1).
std::vector<Spectral> single_part(const PupilFunction& pupil) {
  std::vector<Spectral> parts(1);
  parts[0].pupil = pupil;
  parts[0].sums = require_evaluable(parts[0].pupil);
  parts[0].weight = 1.0;
  parts[0].c = Autocorrelation(parts[0].pupil);
  return parts;
}

}  // namespace

Mtf mtf(const compile::CompiledSystem& system,
        compile::PathId path,
        std::uint16_t field,
        const MtfOptions& options) {
  return mtf(system, path, field, options, trace::RunControl{});
}

Mtf mtf(const compile::CompiledSystem& system,
        compile::PathId path,
        std::uint16_t field,
        const MtfOptions& options,
        const trace::RunControl& control) {
  detail::check_path(system, path);
  check_undersampled(options.undersampled_waves);
  if (options.outputs.map) {
    const std::size_t count =
        options.wavelengths ? options.wavelengths->size() : system.wavelengths_um().size();
    if (count != 1) {
      throw std::invalid_argument("analysis: the 2D MTF map needs exactly one wavelength");
    }
  }
  const std::vector<Spectral> parts =
      spectral_parts(system, path, field, options.wavelengths, options.pupil, control);
  Mtf result = mtf_of(parts, reference_part(system, parts), options.outputs);
  losses_and_warnings(system, path, parts, options.pupil.lost_warning_fraction,
                      options.undersampled_waves, result.losses, result.warnings);
  return result;
}

Mtf mtf(const PupilFunction& pupil, const MtfOutputs& outputs) {
  const std::vector<Spectral> parts = single_part(pupil);
  Mtf result = mtf_of(parts, 0, outputs);
  result.losses = pupil.losses;
  result.warnings = pupil.warnings;
  return result;
}

EncircledEnergy encircled_energy(const compile::CompiledSystem& system,
                                 compile::PathId path,
                                 std::uint16_t field,
                                 std::span<const double> radii,
                                 const EncircledEnergyOptions& options) {
  return encircled_energy(system, path, field, radii, options, trace::RunControl{});
}

EncircledEnergy encircled_energy(const compile::CompiledSystem& system,
                                 compile::PathId path,
                                 std::uint16_t field,
                                 std::span<const double> radii,
                                 const EncircledEnergyOptions& options,
                                 const trace::RunControl& control) {
  detail::check_path(system, path);
  check_undersampled(options.undersampled_waves);
  check_radii(radii);
  const std::vector<Spectral> parts =
      spectral_parts(system, path, field, options.wavelengths, options.pupil, control);
  EncircledEnergy result =
      energy_of_parts(parts, reference_part(system, parts), radii, options.centre);
  losses_and_warnings(system, path, parts, options.pupil.lost_warning_fraction,
                      options.undersampled_waves, result.losses, result.warnings);
  return result;
}

EncircledEnergy encircled_energy(const PupilFunction& pupil,
                                 std::span<const double> radii,
                                 EeCentre centre) {
  check_radii(radii);
  const std::vector<Spectral> parts = single_part(pupil);
  EncircledEnergy result = energy_of_parts(parts, 0, radii, centre);
  result.losses = pupil.losses;
  result.warnings = pupil.warnings;
  return result;
}

}  // namespace rtt::analysis
