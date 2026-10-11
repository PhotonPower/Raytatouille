#include "rtt/analysis/pupil.hpp"

#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "common.hpp"
#include "reference.hpp"
#include "rtt/diagnostics/codes.hpp"

namespace rtt::analysis {

namespace {

/// Antiderivative of h(t) = sqrt(1 - t^2) on [-1, 1]: P(t) = (t sqrt(1 - t^2) + asin t) / 2.
double area_antiderivative(double t) {
  return 0.5 * (t * std::sqrt(1.0 - t * t) + std::asin(t));
}

/// Antiderivative of t h(t): -(1 - t^2)^(3/2) / 3.
double moment_antiderivative(double t) {
  const double r = 1.0 - t * t;
  return -r * std::sqrt(r) / 3.0;
}

/// Antiderivative of h(t)^2 = 1 - t^2: t - t^3 / 3.
double square_antiderivative(double t) {
  return t - t * t * t / 3.0;
}

/// Area and first moments of the part of [x0, x1] x [y0, y1] inside the unit circle.
struct Part {
  double area = 0.0;
  double mx = 0.0;  ///< integral of x over the part
  double my = 0.0;  ///< integral of y over the part
};

/// Closed form of ADR 0032, point 1 (derived there): for t in [x0, x1] with |t| < 1 the cut in y
/// has length l(t) = max(0, min(y1, h) - max(y0, -h)), h = sqrt(1 - t^2). The kinks
/// t = +-sqrt(1 - y0^2), +-sqrt(1 - y1^2) split [x0, x1] into pieces on which the upper end is
/// either y1 or h and the lower end either y0 or -h; the form is decided at the midpoint m of the
/// piece: upper h if h(m) <= y1, lower -h if -h(m) >= y0. The "<=" matters where an edge y = +-1
/// touches the circle in t = 0 (odd N): there h <= y1 on the whole piece with equality in one
/// point only. A piece with upper <= lower at m lies outside (its sign cannot change inside the
/// piece, the crossings are kinks). Area = int (upper - lower) dt, x moment = int t (upper -
/// lower) dt, y moment = int (upper^2 - lower^2) / 2 dt, each with the antiderivatives above.
Part part_inside(double x0, double x1, double y0, double y1) {
  Part part;
  const double lo = std::max(x0, -1.0);
  const double hi = std::min(x1, 1.0);
  if (!(lo < hi)) return part;
  // Ends and at most four kinks (a vector: std::sort on a partly filled std::array trips
  // -Warray-bounds of GCC 16 at -O2).
  std::vector<double> t;
  t.reserve(6);
  t.push_back(lo);
  for (const double y : {y0, y1}) {
    if (std::abs(y) >= 1.0) continue;
    const double kink = std::sqrt(1.0 - y * y);
    for (const double k : {-kink, kink}) {
      if (k > lo && k < hi) t.push_back(k);
    }
  }
  t.push_back(hi);
  std::sort(t.begin(), t.end());
  for (std::size_t j = 0; j + 1 < t.size(); ++j) {
    const double ta = t[j];
    const double tb = t[j + 1];
    if (!(ta < tb)) continue;
    const double m = 0.5 * (ta + tb);
    const double hm = std::sqrt(1.0 - m * m);
    const bool upper_h = hm <= y1;
    const bool lower_h = -hm >= y0;
    const double upper_m = upper_h ? hm : y1;
    const double lower_m = lower_h ? -hm : y0;
    if (!(upper_m > lower_m)) continue;
    const double dt = tb - ta;
    const double dh = area_antiderivative(tb) - area_antiderivative(ta);
    const double dth = moment_antiderivative(tb) - moment_antiderivative(ta);
    const double dt2 = 0.5 * (tb * tb - ta * ta);
    const double dh2 = square_antiderivative(tb) - square_antiderivative(ta);
    const double upper = upper_h ? dh : y1 * dt;
    const double lower = lower_h ? -dh : y0 * dt;
    const double upper_t = upper_h ? dth : y1 * dt2;
    const double lower_t = lower_h ? -dth : y0 * dt2;
    const double upper_sq = upper_h ? dh2 : y1 * y1 * dt;
    const double lower_sq = lower_h ? dh2 : y0 * y0 * dt;
    part.area += upper - lower;
    part.mx += upper_t - lower_t;
    part.my += 0.5 * (upper_sq - lower_sq);
  }
  return part;
}

}  // namespace

PupilCells pupil_cells(int grid) {
  if (grid < 2 || grid > kMaxPupilGrid) {
    throw std::invalid_argument("analysis: pupil grid must lie in [2, " +
                                std::to_string(kMaxPupilGrid) + "], got " + std::to_string(grid));
  }
  const int n = grid;
  const auto size = static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
  PupilCells cells;
  cells.grid = n;
  cells.fraction.resize(size);
  cells.px.resize(size);
  cells.py.resize(size);
  const double nd = static_cast<double>(n);
  const double a = 2.0 / nd;
  // Edges and centres as (2 i - N) / N and (2 i + 1 - N) / N: the mirror cell has the negated
  // numerator, so the values are mirror images bit for bit.
  const auto edge = [nd, n](int i) { return static_cast<double>(2 * i - n) / nd; };
  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < n; ++i) {
      const std::size_t j =
          static_cast<std::size_t>(k) * static_cast<std::size_t>(n) + static_cast<std::size_t>(i);
      const double cx = static_cast<double>(2 * i + 1 - n) / nd;
      const double cy = static_cast<double>(2 * k + 1 - n) / nd;
      cells.px[j] = cx;
      cells.py[j] = cy;
      double x0 = edge(i);
      double x1 = edge(i + 1);
      double y0 = edge(k);
      double y1 = edge(k + 1);
      // The disk is convex: full if the farthest corner lies inside, empty if the nearest point
      // of the cell lies outside.
      const double far_x = std::max(std::abs(x0), std::abs(x1));
      const double far_y = std::max(std::abs(y0), std::abs(y1));
      const double near_x = x0 <= 0.0 && x1 >= 0.0 ? 0.0 : std::min(std::abs(x0), std::abs(x1));
      const double near_y = y0 <= 0.0 && y1 >= 0.0 ? 0.0 : std::min(std::abs(y0), std::abs(y1));
      if (far_x * far_x + far_y * far_y <= 1.0) {
        cells.fraction[j] = 1.0;
        continue;
      }
      if (near_x * near_x + near_y * near_y >= 1.0) {
        cells.fraction[j] = 0.0;
        continue;
      }
      // Boundary cell. Cells on the negative side are computed mirrored, so that mirror cells
      // give the same numbers (x -> -x, y -> -y).
      const double sx = x1 <= 0.0 ? -1.0 : 1.0;
      const double sy = y1 <= 0.0 ? -1.0 : 1.0;
      if (sx < 0.0) {
        x0 = -edge(i + 1);
        x1 = -edge(i);
      }
      if (sy < 0.0) {
        y0 = -edge(k + 1);
        y1 = -edge(k);
      }
      // Integrated along the coordinate nearer to the axis: the antiderivatives are then evaluated
      // away from t = +-1, where asin amplifies the rounding of t by 1 / sqrt(1 - t^2) (up to
      // 5e-13 in f at N = 64 the other way round). The swapped cell (x <-> y) takes the same
      // path with the roles exchanged, so f and the centroid are symmetric bit for bit.
      // A cell on a diagonal (the same x and y range after mirroring) is its own swap image, so
      // its two moments are equal; their mean puts the centroid on the diagonal bit for bit (the
      // two formulas alone differ by up to 3e-14 in the centroid at N = 64).
      Part part;
      if (x0 == y0 && x1 == y1) {
        part = part_inside(x0, x1, y0, y1);
        part.mx = part.my = 0.5 * (part.mx + part.my);
      } else if (std::max(std::abs(x0), std::abs(x1)) > std::max(std::abs(y0), std::abs(y1))) {
        const Part swapped = part_inside(y0, y1, x0, x1);
        part = {swapped.area, swapped.my, swapped.mx};
      } else {
        part = part_inside(x0, x1, y0, y1);
      }
      if (!(part.area > 0.0)) {
        cells.fraction[j] = 0.0;  // touches the circle in a point only
        continue;
      }
      cells.fraction[j] = std::min(1.0, part.area / (a * a));
      cells.px[j] = sx * (part.mx / part.area);
      cells.py[j] = sy * (part.my / part.area);
    }
  }
  return cells;
}

namespace {

using compile::CompiledSystem;
using compile::PathId;

/// Step of the parabasal rays of the mapping in pupil units (ADR 0032, point 4): the same step as
/// the field-curvature analysis (#31); the central difference has an error O(step^2).
constexpr double kMappingStep = 1e-3;
/// |x_B x d_c| at or below this: the chief ray grazes the image surface, no frame (point 4).
constexpr double kGrazingLimit = 1e-6;
/// |det M| at or below this times n'^2: afocal image space, no mapping (point 4).
constexpr double kAfocalDeterminant = 1e-12;

std::string path_pointer(PathId path) {
  return "/paths/" + std::to_string(path.index);
}

void check_options(const PupilFunctionOptions& o) {
  // grid is checked by pupil_cells; lost_warning_fraction by the LossCounter.
  if (!std::isfinite(o.apodization) || !(o.apodization >= 0.0)) {
    throw std::invalid_argument("analysis: apodization must be finite and >= 0");
  }
  if (!std::isfinite(o.focus_shift)) {
    throw std::invalid_argument("analysis: focus_shift must be finite");
  }
  if (!(o.mapping_warning >= 0.0 && o.mapping_warning <= 1.0)) {
    throw std::invalid_argument("analysis: mapping_warning must lie in [0, 1]");
  }
  if (!(o.lost_warning_fraction >= 0.0 && o.lost_warning_fraction <= 1.0)) {
    throw std::invalid_argument("analysis: lost_warning_fraction must lie in [0, 1]");
  }
}

/// ADR 0032, point 9 (H6): the OPD of a crystal mode needs the index of the mode (#35).
void reject_crystal_modes(const CompiledSystem& system, PathId path) {
  for (const auto& event : system.path(path).events) {
    if (event.kind == model::EventKind::Ordinary || event.kind == model::EventKind::Extraordinary) {
      throw AnalysisError("analysis: paths with crystal modes are not supported (#35)",
                          path_pointer(path));
    }
  }
}

math::Vec3 end_pos(const PupilEndStates& e, std::size_t j) {
  return {e.pos_x[j], e.pos_y[j], e.pos_z[j]};
}

math::Vec3 end_dir(const PupilEndStates& e, std::size_t j) {
  return {e.dir_x[j], e.dir_y[j], e.dir_z[j]};
}

/// Direction sines s = n' (d . e_x, d . e_y) (ADR 0032, point 4).
std::array<double, 2> direction_sines(const math::Vec3& d, const PupilMapping& m, double n_image) {
  return {n_image * d.dot(m.e_x), n_image * d.dot(m.e_y)};
}

/// W, status and amplitude of every cell, the sphere, the losses and arrived_fraction from the
/// end states against `ref`; the same code for pupil_function() and refocus(), so both agree bit
/// for bit. The warnings of the losses come first, then the mapping (finish_mapping()).
void evaluate(const CompiledSystem& system, PupilFunction& pf, const detail::Reference& ref) {
  const PupilEndStates& e = pf.end;
  pf.sphere = ref.sphere;
  const double opl_chief =
      detail::chief_opl(e.chief_pos, e.chief_dir, e.chief_opl, ref, e.image_index);
  detail::LossCounter losses(system, pf.path, pf.lost_warning_fraction);
  double sum_all = 0.0;
  double sum_arrived = 0.0;
  for (std::size_t j = 0; j < pf.fraction.size(); ++j) {
    sum_all += pf.fraction[j];
    pf.w[j] = 0.0;
    pf.amplitude[j] = 0.0;
    if (pf.fraction[j] == 0.0) {
      pf.status[j] = trace::RayStatus::Vignetted;  // not traced, not counted
      continue;
    }
    const bool on_image =
        e.ray_status[j] == trace::RayStatus::Alive && e.last_surface[j] == e.image_surface;
    const OpdPoint q = detail::point_against(end_pos(e, j), end_dir(e, j), e.opl[j], on_image,
                                             e.ray_status[j], opl_chief, ref, e.image_index);
    pf.status[j] = q.status;
    losses.add_ray(q.status, e.last_surface[j]);
    if (q.status != trace::RayStatus::Alive) continue;
    pf.w[j] = q.w;
    // ADR 0032, point 3 and 8: A = sqrt(weight) exp(-G rho^2), rho of the sample point (Ansys
    // Zemax OpticStudio User Guide, Release 2025 R1, "Apodization Type").
    const double rho2 = pf.px[j] * pf.px[j] + pf.py[j] * pf.py[j];
    pf.amplitude[j] = std::sqrt(e.weight[j]) * std::exp(-pf.apodization * rho2);
    sum_arrived += pf.fraction[j];
  }
  pf.arrived_fraction = sum_all > 0.0 ? sum_arrived / sum_all : 0.0;
  pf.losses = losses.result();
  pf.warnings = losses.warnings();
}

/// delta_max and delta_rms of the mapping over the arrived cells (ADR 0032, point 4; delta_rms
/// weighted with f like the other pupil statistics), centre = C', and the warning.
void finish_mapping(PupilFunction& pf) {
  if (!pf.mapping) return;
  PupilMapping& m = *pf.mapping;
  m.centre = pf.sphere.centre;
  // ||M||_2: the largest singular value of the 2 x 2 matrix, sqrt((T + sqrt(T^2 - 4 D^2)) / 2)
  // with T the squared Frobenius norm and D the determinant (eigenvalues of M^T M).
  const double tr = m.m[0] * m.m[0] + m.m[1] * m.m[1] + m.m[2] * m.m[2] + m.m[3] * m.m[3];
  const double det = m.m[0] * m.m[3] - m.m[1] * m.m[2];
  const double norm = std::sqrt(0.5 * (tr + std::sqrt(std::max(0.0, tr * tr - 4.0 * det * det))));
  double max = 0.0;
  double sum_f = 0.0;
  double sum_fd2 = 0.0;
  for (std::size_t j = 0; j < pf.fraction.size(); ++j) {
    if (pf.fraction[j] == 0.0 || pf.status[j] != trace::RayStatus::Alive) continue;
    const std::array<double, 2> s = direction_sines(end_dir(pf.end, j), m, pf.end.image_index);
    const double rx = s[0] - m.s0[0] - (m.m[0] * pf.px[j] + m.m[1] * pf.py[j]);
    const double ry = s[1] - m.s0[1] - (m.m[2] * pf.px[j] + m.m[3] * pf.py[j]);
    const double d = std::hypot(rx, ry) / norm;
    max = std::max(max, d);
    sum_f += pf.fraction[j];
    sum_fd2 += pf.fraction[j] * d * d;
  }
  m.delta_max = max;
  m.delta_rms = sum_f > 0.0 ? std::sqrt(sum_fd2 / sum_f) : 0.0;
  if (m.delta_max > pf.mapping_warning) {
    constexpr diagnostics::DiagnosticCode kCode = "pupil.mapping_nonlinear";
    pf.warnings.push_back(
        {kCode.severity(), path_pointer(pf.path),
         "pupil coordinates map nonlinearly onto the direction sines: delta_max = " +
             std::to_string(m.delta_max) + " pupil units above " +
             std::to_string(pf.mapping_warning) + " (ADR 0032, point 4)",
         std::string(kCode.str())});
  }
}

/// The affine mapping s ~ s0 + M p (ADR 0032, point 4), without the deltas; none for an afocal
/// image space.
std::optional<PupilMapping> make_mapping(const CompiledSystem& system,
                                         PathId path,
                                         std::uint16_t field,
                                         std::uint16_t wavelength,
                                         trace::Aiming aiming,
                                         const detail::Reference& ref,
                                         const math::Vec3& chief_dir) {
  const std::uint32_t image = detail::image_surface(system, path);
  const double n_image = detail::image_index(system, path, wavelength);
  PupilMapping m;
  m.d_c = ref.chief_dir;
  const math::Vec3 x_b = system.surfaces()[image].to_global.apply_vector(math::Vec3::UnitX());
  m.e_x = (x_b - x_b.dot(m.d_c) * m.d_c).normalized();
  m.e_y = m.d_c.cross(m.e_x);
  const trace::PointListPupil parabasal{
      {{kMappingStep, 0.0}, {-kMappingStep, 0.0}, {0.0, kMappingStep}, {0.0, -kMappingStep}}};
  const trace::RayBatch rays =
      detail::trace_rays(system, path, field, wavelength, parabasal, aiming);
  std::array<std::array<double, 2>, 4> s{};
  for (std::size_t i = 0; i < 4; ++i) {
    if (!detail::arrived(rays, i, image)) {
      detail::throw_lost(system, rays, i, field,
                         "analysis: a parabasal ray of the pupil mapping does not reach the image "
                         "surface");
    }
    s[i] = direction_sines({rays.dir_x()[i], rays.dir_y()[i], rays.dir_z()[i]}, m, n_image);
  }
  // Central differences, columns px and py: error O(kMappingStep^2).
  const double h2 = 2.0 * kMappingStep;
  m.m = {(s[0][0] - s[1][0]) / h2, (s[2][0] - s[3][0]) / h2, (s[0][1] - s[1][1]) / h2,
         (s[2][1] - s[3][1]) / h2};
  const double det = m.m[0] * m.m[3] - m.m[1] * m.m[2];
  if (!(std::abs(det) > kAfocalDeterminant * n_image * n_image)) return std::nullopt;
  // s0: 0 at the reference wavelength by construction (the chief ray of d_c itself).
  if (wavelength != system.reference_wavelength()) m.s0 = direction_sines(chief_dir, m, n_image);
  return m;
}

}  // namespace

PupilFunction pupil_function(const compile::CompiledSystem& system,
                             compile::PathId path,
                             std::uint16_t field,
                             std::uint16_t wavelength,
                             const PupilFunctionOptions& options) {
  return pupil_function(system, path, field, wavelength, options, trace::RunControl{});
}

PupilFunction pupil_function(const compile::CompiledSystem& system,
                             compile::PathId path,
                             std::uint16_t field,
                             std::uint16_t wavelength,
                             const PupilFunctionOptions& options,
                             const trace::RunControl& control) {
  detail::check_path(system, path);
  detail::check_wavelength(system, wavelength);
  check_options(options);
  PupilCells cells = pupil_cells(options.grid);  // checks the grid
  reject_crystal_modes(system, path);            // before any ray is traced

  const detail::Reference ref =
      detail::make_reference(system, path, field, options.aiming, options.focus_shift);
  const std::uint32_t image = detail::image_surface(system, path);
  // ADR 0032, point 4: the frame of the PSF plane needs x_B not parallel to d_c.
  const math::Vec3 x_b = system.surfaces()[image].to_global.apply_vector(math::Vec3::UnitX());
  if (x_b.cross(ref.chief_dir).norm() <= kGrazingLimit) {
    throw AnalysisError(
        "analysis: the chief ray grazes the image surface (its x axis is parallel to the chief "
        "ray); the plane of the pupil mapping has no frame",
        system.surfaces()[image].location);
  }

  PupilFunction pf;
  pf.path = path;
  pf.field = field;
  pf.wavelength = wavelength;
  pf.wavelength_um = system.wavelengths_um()[wavelength];
  pf.reference_wavelength_um = system.wavelengths_um()[system.reference_wavelength()];
  pf.grid = options.grid;
  pf.focus_shift = options.focus_shift;
  pf.apodization = options.apodization;
  pf.lost_warning_fraction = options.lost_warning_fraction;
  pf.mapping_warning = options.mapping_warning;
  const std::size_t size = cells.fraction.size();
  pf.fraction = std::move(cells.fraction);
  pf.px = std::move(cells.px);
  pf.py = std::move(cells.py);
  pf.amplitude.assign(size, 0.0);
  pf.w.assign(size, 0.0);
  pf.status.assign(size, trace::RayStatus::Vignetted);

  // End states: the chief ray of this wavelength and the cells with f > 0 as a point list
  // (ADR 0032, point 2), traced like opd_points().
  PupilEndStates& e = pf.end;
  e.image_surface = image;
  e.image_index = detail::image_index(system, path, wavelength);
  e.reference_pos = ref.chief_pos;
  e.reference_dir = ref.chief_dir;
  e.exit_pupil = ref.exit_pupil;
  const trace::RayBatch chief = detail::chief_ray(system, path, field, wavelength, options.aiming);
  e.chief_pos = {chief.pos_x()[0], chief.pos_y()[0], chief.pos_z()[0]};
  e.chief_dir = {chief.dir_x()[0], chief.dir_y()[0], chief.dir_z()[0]};
  e.chief_opl = chief.opl()[0];
  trace::PointListPupil list;
  for (std::size_t j = 0; j < size; ++j) {
    if (pf.fraction[j] > 0.0) list.points.push_back({pf.px[j], pf.py[j]});
  }
  const trace::RayBatch rays =
      detail::trace_rays(system, path, field, wavelength, list, options.aiming, control);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (auto* column :
       {&e.pos_x, &e.pos_y, &e.pos_z, &e.dir_x, &e.dir_y, &e.dir_z, &e.opl, &e.weight}) {
    column->assign(size, nan);
  }
  e.ray_status.assign(size, trace::RayStatus::Vignetted);
  e.last_surface.assign(size, trace::kNoSurface);
  if (options.with_prt) pf.prt.assign(9 * size, std::complex<double>(0.0, 0.0));
  std::size_t i = 0;
  for (std::size_t j = 0; j < size; ++j) {
    if (pf.fraction[j] == 0.0) continue;
    e.pos_x[j] = rays.pos_x()[i];
    e.pos_y[j] = rays.pos_y()[i];
    e.pos_z[j] = rays.pos_z()[i];
    e.dir_x[j] = rays.dir_x()[i];
    e.dir_y[j] = rays.dir_y()[i];
    e.dir_z[j] = rays.dir_z()[i];
    e.opl[j] = rays.opl()[i];
    e.weight[j] = rays.weight()[i];
    e.ray_status[j] = rays.status()[i];
    e.last_surface[j] = rays.last_surface()[i];
    if (options.with_prt) {
      for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
          pf.prt[9 * j + static_cast<std::size_t>(3 * row + col)] = rays.prt(row, col)[i];
        }
      }
    }
    ++i;
  }

  evaluate(system, pf, ref);
  pf.mapping = make_mapping(system, path, field, wavelength, options.aiming, ref, e.chief_dir);
  finish_mapping(pf);
  return pf;
}

PupilFunction refocus(const compile::CompiledSystem& system,
                      const PupilFunction& pupil,
                      double focus_shift) {
  if (!std::isfinite(focus_shift)) {
    throw std::invalid_argument("analysis: focus_shift must be finite");
  }
  detail::check_path(system, pupil.path);
  // Review of #227 (f2): a system other than the one of `pupil` would make the loss counter read
  // surfaces out of range. The cheap guards catch the usual mix-up: the same image surface, the
  // field and wavelength in range, every last surface a surface of `system`.
  detail::check_wavelength(system, pupil.wavelength);
  const auto surfaces = system.surfaces().size();
  if (pupil.field >= system.fields().points.size() ||
      pupil.end.image_surface != detail::image_surface(system, pupil.path) ||
      std::any_of(
          pupil.end.last_surface.begin(), pupil.end.last_surface.end(),
          [surfaces](std::uint32_t s) { return s != trace::kNoSurface && s >= surfaces; })) {
    throw std::invalid_argument(
        "analysis: refocus needs the system that produced the pupil function");
  }
  PupilFunction pf = pupil;
  pf.focus_shift = focus_shift;
  const detail::Reference ref =
      detail::reference_at(pf.end.reference_pos, pf.end.reference_dir, pf.end.exit_pupil,
                           focus_shift, pf.reference_wavelength_um * 1e-3);
  evaluate(system, pf, ref);
  finish_mapping(pf);
  return pf;
}

}  // namespace rtt::analysis
