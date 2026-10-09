#include "rtt/optim/levenberg_marquardt.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/partitioner.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Sources (docs/quellen.md): K. Madsen, H. B. Nielsen, O. Tingleff, "Methods for Non-Linear Least
// Squares Problems", 2nd ed., IMM DTU 2004 (MNT); H. B. Nielsen, "Damping Parameter in Marquardt's
// Method", IMM-REP-1999-05 (1999); J. J. Moré, B. S. Garbow, K. E. Hillstrom, "User Guide for
// MINPACK-1", ANL-80-74 (1980). Decisions: ADR 0030.

namespace rtt::optim {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Nielsen 1999, eq. (2.5) with beta = nu_0 = 2, gamma = 3, p = 3 (the same as MNT eq. (2.21)).
constexpr double kBeta = 2.0;
constexpr double kGamma = 3.0;

bool same_bits(double a, double b) noexcept {
  return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

bool finite(double x) noexcept {
  return std::isfinite(x);
}

// 1/2 f^T f (MNT eq. (3.1b)), summed in index order.
double half_square(std::span<const double> f) noexcept {
  double sum = 0.0;
  for (const double v : f) sum += v * v;
  return 0.5 * sum;
}

// Column-major m x n matrix.
struct Columns {
  std::size_t m = 0;
  std::size_t n = 0;
  std::vector<double> values;
  [[nodiscard]] double operator()(std::size_t i, std::size_t j) const { return values[j * m + i]; }
  double& operator()(std::size_t i, std::size_t j) { return values[j * m + i]; }
};

void check_input(const Residuals& residuals,
                 std::size_t m,
                 const std::vector<double>& p0,
                 const std::vector<Bounds>& bounds,
                 const LmOptions& o) {
  const auto fail = [](const char* what) {
    throw std::invalid_argument(std::string("levenberg_marquardt: ") + what);
  };
  if (!residuals) fail("no residual function");
  if (m == 0) fail("m must be at least 1");
  if (p0.empty()) fail("at least one variable is needed");
  if (bounds.size() != p0.size()) fail("one Bounds entry per variable is needed");
  for (std::size_t j = 0; j < p0.size(); ++j) {
    const Bounds& b = bounds[j];
    if (!finite(p0[j])) fail("start values must be finite");
    if ((b.min && !finite(*b.min)) || (b.max && !finite(*b.max))) fail("bounds must be finite");
    if (b.min && b.max && !(*b.min < *b.max)) fail("min must be less than max");
    if ((b.min && p0[j] < *b.min) || (b.max && p0[j] > *b.max)) {
      fail("start values must lie within the bounds");
    }
  }
  if (o.max_iterations < 0) fail("max_iterations must not be negative");
  if (!finite(o.ftol) || o.ftol < 0.0) fail("ftol must be finite and >= 0");
  if (!finite(o.xtol) || o.xtol < 0.0) fail("xtol must be finite and >= 0");
  if (!finite(o.gtol) || o.gtol < 0.0) fail("gtol must be finite and >= 0");
  if (!finite(o.tau) || !(o.tau > 0.0)) fail("tau must be finite and > 0");
  if (!finite(o.function_precision) || !(o.function_precision > 0.0)) {
    fail("function_precision must be finite and > 0");
  }
}

// True if the stage was interrupted: a cancellation request (returned as true) or an exception of
// the progress callback (rethrown by finish()).
bool interrupted(trace::RunMonitor& monitor) {
  if (!monitor.stop()) return false;
  try {
    monitor.finish();
  } catch (const trace::Cancelled&) {
    return true;
  }
  return true;
}

// Solves A x = b for a symmetric n x n matrix A (row-major) by the Cholesky factorization
// A = C^T C with C upper triangular, MNT Appendix A, Algorithm A.4 (with its test for positive
// definiteness: d = a_kk - sum c_ik^2 > 0), then forward substitution C^T z = b and back
// substitution C x = z (MNT p. 52). Returns false if A is not positive definite or x is not
// finite. Sums run in index order (ADR 0030, point 11).
bool cholesky_solve(std::vector<double> a,
                    std::size_t n,
                    std::span<const double> b,
                    std::span<double> x) {
  // a is overwritten with C in its upper triangle.
  for (std::size_t k = 0; k < n; ++k) {
    double d = a[k * n + k];
    for (std::size_t i = 0; i < k; ++i) d -= a[i * n + k] * a[i * n + k];
    if (!(d > 0.0)) return false;
    const double ckk = std::sqrt(d);
    a[k * n + k] = ckk;
    for (std::size_t j = k + 1; j < n; ++j) {
      double s = a[k * n + j];
      for (std::size_t i = 0; i < k; ++i) s -= a[i * n + j] * a[i * n + k];
      a[k * n + j] = s / ckk;
    }
  }
  for (std::size_t k = 0; k < n; ++k) {  // C^T z = b, z in x
    double s = b[k];
    for (std::size_t i = 0; i < k; ++i) s -= a[i * n + k] * x[i];
    x[k] = s / a[k * n + k];
  }
  for (std::size_t k = n; k-- > 0;) {  // C x = z
    double s = x[k];
    for (std::size_t j = k + 1; j < n; ++j) s -= a[k * n + j] * x[j];
    x[k] = s / a[k * n + k];
  }
  return std::all_of(x.begin(), x.end(), finite);
}

class Run {
 public:
  Run(const Residuals& residuals,
      std::size_t m,
      std::vector<double> p0,
      std::vector<Bounds> bounds,
      const LmOptions& options,
      const trace::RunControl& control)
      : residuals_(residuals),
        m_(m),
        n_(p0.size()),
        p0_(std::move(p0)),
        bounds_(std::move(bounds)),
        options_(options),
        control_(control),
        optimize_(control, static_cast<std::size_t>(options.max_iterations), "optimize") {
    theta0_.resize(n_);
    for (std::size_t j = 0; j < n_; ++j) theta0_[j] = to_internal(p0_[j], bounds_[j]);
    theta_ = theta0_;
  }

  LmResult run();

 private:
  enum class JacobianOutcome : std::uint8_t { Ok, Cancelled, Failed };

  // External values of theta. A variable whose internal value equals the start keeps its input
  // value bitwise, so the round trip p -> theta -> p never writes 1-ulp changes (ADR 0030, point
  // 12).
  void external(std::span<const double> theta, std::span<double> p) const {
    for (std::size_t j = 0; j < n_; ++j) {
      p[j] = same_bits(theta[j], theta0_[j]) ? p0_[j] : to_external(theta[j], bounds_[j]);
    }
  }

  // One evaluation at external values p; false for an invalid result, an exception or a non-finite
  // residual (ADR 0030, point 10).
  bool evaluate(std::span<const double> p, std::span<double> f) const noexcept {
    try {
      if (!residuals_(p, f)) return false;
    } catch (...) {
      return false;
    }
    return std::all_of(f.begin(), f.end(), finite);
  }

  JacobianOutcome jacobian();
  void normal_equations();
  LmResult finish(LmStatus status);

  const Residuals& residuals_;
  std::size_t m_;
  std::size_t n_;
  std::vector<double> p0_;
  std::vector<Bounds> bounds_;
  const LmOptions& options_;
  const trace::RunControl& control_;
  trace::RunMonitor optimize_;

  std::vector<double> theta0_;
  std::vector<double> theta_;  // last accepted state
  std::vector<double> f_;      // residuals there
  double merit_ = kNaN;
  Columns j_;              // Jacobian in theta at theta_
  std::vector<double> d_;  // scaling D (MINPACK-1 section 2.5), fixed after the start
  std::vector<double> a_;  // J^T J, n x n row-major
  std::vector<double> g_;  // J^T f
  std::size_t evaluations_ = 0;
  std::size_t failed_evaluations_ = 0;
  int k_ = 0;
  std::optional<std::size_t> failed_variable_;
  std::vector<LmIteration> history_;
};

// Central difference in theta (ADR 0007, ADR 0030 point 9), all 2n evaluations in parallel.
//
// Step h_j = eps_f^(1/3) max(|theta_j|, 1). Derivation (no literature source, docs/quellen.md
// "hergeleitet"): the central difference (f(x + h) - f(x - h)) / (2h) has the truncation error h^2
// |f'''| / 6 and, with residuals known to the absolute noise eps_f |f| (eps_f relative to the fixed
// scale ||f(theta_0)||_inf), the rounding error eps_f |f| / h. Their sum is smallest for h = (3
// eps_f |f| / |f'''|)^(1/3), so h is proportional to eps_f^(1/3) when |f| / |f'''| is of the order
// max(|theta|, 1)^3. Model: the forward difference of MINPACK-1, section 2.4, with h = sqrt(EPSFCN)
// |x| (sqrt for the forward difference). The denominator is the difference of the two points as
// represented, (theta + h) - (theta - h), so the rounding of theta +- h does not enter the slope.
Run::JacobianOutcome Run::jacobian() {
  const std::size_t count = 2 * n_;
  std::vector<double> plus(n_);
  std::vector<double> minus(n_);
  const double factor = std::cbrt(options_.function_precision);
  for (std::size_t j = 0; j < n_; ++j) {
    const double h = factor * std::max(std::abs(theta_[j]), 1.0);
    plus[j] = theta_[j] + h;
    minus[j] = theta_[j] - h;
  }
  // Every evaluation writes only its own slot (ADR 0004): values, state.
  std::vector<double> values(count * m_);
  enum : std::uint8_t { kSkipped = 0, kValid = 1, kInvalid = 2 };
  std::vector<std::uint8_t> state(count, kSkipped);
  trace::RunMonitor monitor(control_, count, "jacobian");
  const auto body = [&](const oneapi::tbb::blocked_range<std::size_t>& range) {
    std::vector<double> theta = theta_;
    std::vector<double> p(n_);
    for (std::size_t idx = range.begin(); idx != range.end(); ++idx) {
      if (monitor.stop() || optimize_.stop()) continue;  // checked before every evaluation
      const std::size_t j = idx / 2;
      theta[j] = idx % 2 == 0 ? plus[j] : minus[j];
      external(theta, p);
      const std::span<double> slot(values.data() + idx * m_, m_);
      state[idx] = evaluate(p, slot) ? kValid : kInvalid;
      theta[j] = theta_[j];
      monitor.add(1);
    }
  };
  const oneapi::tbb::blocked_range<std::size_t> all(0, count, 1);
  if (control_.active()) {
    // Grain 1: the cancellation is checked before every evaluation (ADR 0004, addendum #83).
    oneapi::tbb::parallel_for(all, body, oneapi::tbb::simple_partitioner());
  } else {
    oneapi::tbb::parallel_for(all, body, oneapi::tbb::static_partitioner());
  }
  // Counters in index order (integers, ADR 0030 point 11).
  bool skipped = false;
  for (const std::uint8_t s : state) {
    if (s == kSkipped) {
      skipped = true;
    } else {
      ++evaluations_;
      if (s == kInvalid) ++failed_evaluations_;
    }
  }
  if (skipped) {
    if (interrupted(monitor) || interrupted(optimize_)) return JacobianOutcome::Cancelled;
  } else {
    try {
      monitor.finish();  // the final report of the stage
    } catch (const trace::Cancelled&) {
      // A request after the last evaluation of the stage: the run ends here with the last
      // accepted state, like a request before an evaluation.
      return JacobianOutcome::Cancelled;
    }
  }
  for (std::size_t j = 0; j < n_; ++j) {
    if (state[2 * j] != kValid || state[2 * j + 1] != kValid) {
      failed_variable_ = j;  // smallest index (ADR 0030, point 10)
      return JacobianOutcome::Failed;
    }
  }
  j_ = Columns{m_, n_, std::vector<double>(m_ * n_)};
  for (std::size_t j = 0; j < n_; ++j) {
    const double denominator = plus[j] - minus[j];
    const double* fp = values.data() + 2 * j * m_;
    const double* fm = values.data() + (2 * j + 1) * m_;
    for (std::size_t i = 0; i < m_; ++i) j_(i, j) = (fp[i] - fm[i]) / denominator;
  }
  return JacobianOutcome::Ok;
}

// A = J^T J and g = J^T f, summed in index order.
void Run::normal_equations() {
  a_.assign(n_ * n_, 0.0);
  g_.assign(n_, 0.0);
  for (std::size_t r = 0; r < n_; ++r) {
    double gr = 0.0;
    for (std::size_t i = 0; i < m_; ++i) gr += j_(i, r) * f_[i];
    g_[r] = gr;
    for (std::size_t c = r; c < n_; ++c) {
      double sum = 0.0;
      for (std::size_t i = 0; i < m_; ++i) sum += j_(i, r) * j_(i, c);
      a_[r * n_ + c] = sum;
      a_[c * n_ + r] = sum;
    }
  }
}

LmResult Run::finish(LmStatus status) {
  LmResult r;
  r.status = status;
  r.p.resize(n_);
  external(theta_, r.p);
  r.changed.resize(n_);
  r.at_bound.resize(n_);
  for (std::size_t j = 0; j < n_; ++j) {
    r.changed[j] = same_bits(theta_[j], theta0_[j]) ? 0 : 1;
    r.at_bound[j] = at_bound(r.p[j], bounds_[j]) ? 1 : 0;
  }
  r.f = f_;
  r.F = f_.empty() ? kNaN : merit_;
  r.iterations = k_;
  r.evaluations = evaluations_;
  r.failed_evaluations = failed_evaluations_;
  r.failed_variable = failed_variable_;
  r.history = std::move(history_);
  if (status != LmStatus::Cancelled) {
    // The final report; an exception of the callback propagates.
    try {
      optimize_.finish();
    } catch (const trace::Cancelled&) {
      // A request after the last evaluation: the state is complete, the run counts as
      // cancelled (the call honours every request, RunControl).
      r.status = LmStatus::Cancelled;
    }
  }
  return r;
}

// MNT Algorithm 3.16 in the scaled variables z = D theta, with the order of the stopping tests of
// ADR 0030, point 8 (steps 1-7).
LmResult Run::run() {
  // Step 1: start.
  if (interrupted(optimize_)) return finish(LmStatus::Cancelled);
  f_.resize(m_);
  ++evaluations_;
  if (!evaluate(p0_, f_)) {
    throw std::invalid_argument("levenberg_marquardt: residuals invalid at the start point");
  }
  merit_ = half_square(f_);
  switch (jacobian()) {
    case JacobianOutcome::Cancelled:
      return finish(LmStatus::Cancelled);
    case JacobianOutcome::Failed:
      return finish(LmStatus::Failed);
    case JacobianOutcome::Ok:
      break;
  }
  // Scaling D: column norms of the start Jacobian, 1 for a zero column (MINPACK-1, section 2.5,
  // MODE = 1); fixed for the run (ADR 0030, point 6).
  d_.resize(n_);
  for (std::size_t j = 0; j < n_; ++j) {
    double sum = 0.0;
    for (std::size_t i = 0; i < m_; ++i) sum += j_(i, j) * j_(i, j);
    d_[j] = sum > 0.0 ? std::sqrt(sum) : 1.0;
  }
  normal_equations();
  // ||D^-1 g||_inf <= gtol (MNT eq. (3.15a) in z).
  const auto gradient_converged = [this] {
    double norm = 0.0;
    for (std::size_t j = 0; j < n_; ++j) {
      norm = std::max(norm, std::abs(g_[j]) / d_[j]);
    }
    return norm <= options_.gtol;
  };
  if (gradient_converged()) return finish(LmStatus::ConvergedGradient);
  // mu_0 = tau max_i (J_z^T J_z)_ii (MNT eq. (3.14)); nu_0 = beta (Nielsen eq. (2.5)).
  double max_diagonal = 0.0;
  for (std::size_t j = 0; j < n_; ++j) {
    max_diagonal = std::max(max_diagonal, a_[j * n_ + j] / (d_[j] * d_[j]));
  }
  double mu = options_.tau * max_diagonal;
  double nu = kBeta;

  std::vector<double> h(n_);
  std::vector<double> minus_g(n_);
  std::vector<double> theta_new(n_);
  std::vector<double> p_new(n_);
  std::vector<double> f_new(m_);
  while (true) {
    // Step 7: k >= k_max.
    if (k_ >= options_.max_iterations) return finish(LmStatus::MaxIterations);
    // Step 2: solve; k counts every solve.
    ++k_;
    LmIteration it{k_, merit_, mu, kNaN, kNaN, false, evaluations_};
    const auto rejected = [&] {
      mu *= nu;  // Nielsen eq. (2.5), else branch
      nu *= 2.0;
    };
    const auto record = [&] {
      it.F = merit_;
      it.evaluations = evaluations_;
      history_.push_back(it);
      optimize_.add(1);
    };
    // (J^T J + mu D^T D) h = -g (MNT eq. (3.13) in z) via the normal equations (MNT example 3.6)
    // and the Cholesky factorization of MNT Algorithm A.4; a failure is a rejected trial.
    std::vector<double> system = a_;
    for (std::size_t j = 0; j < n_; ++j) system[j * n_ + j] += mu * d_[j] * d_[j];
    for (std::size_t j = 0; j < n_; ++j) minus_g[j] = -g_[j];
    const bool solved = cholesky_solve(std::move(system), n_, minus_g, h);
    if (!solved) {  // a rejected trial without an evaluation
      rejected();
      record();
      continue;
    }
    // Step 3: the step test before f(theta + h) (MNT eq. (3.15b) in z); h is not applied.
    double step_norm_sq = 0.0;
    double theta_norm_sq = 0.0;
    for (std::size_t j = 0; j < n_; ++j) {
      const double dh = d_[j] * h[j];
      const double dt = d_[j] * theta_[j];
      step_norm_sq += dh * dh;
      theta_norm_sq += dt * dt;
    }
    it.step_norm = std::sqrt(step_norm_sq);
    if (it.step_norm <= options_.xtol * (std::sqrt(theta_norm_sq) + options_.xtol)) {
      record();
      return finish(LmStatus::ConvergedStep);
    }
    // Step 4: the trial, after the cancellation check (ADR 0030, point 11).
    if (interrupted(optimize_)) {
      history_.push_back(it);
      return finish(LmStatus::Cancelled);
    }
    for (std::size_t j = 0; j < n_; ++j) theta_new[j] = theta_[j] + h[j];
    external(theta_new, p_new);
    ++evaluations_;
    if (!evaluate(p_new, f_new)) {
      ++failed_evaluations_;
      rejected();
      record();
      continue;
    }
    const double merit_new = half_square(f_new);
    // L(0) - L(h) = 1/2 h^T (mu D^T D h - g) (MNT p. 25 in z); positive in exact arithmetic, a
    // non-positive value from rounding counts as no predicted gain.
    double predicted = 0.0;
    for (std::size_t j = 0; j < n_; ++j) {
      predicted += h[j] * (mu * d_[j] * d_[j] * h[j] - g_[j]);
    }
    predicted *= 0.5;
    const double merit_old = merit_;
    // Gain ratio (MNT eq. (2.18)).
    const double rho = predicted > 0.0 ? (merit_old - merit_new) / predicted : 0.0;
    it.rho = rho;
    it.accepted = rho > 0.0;
    if (it.accepted) {
      theta_ = theta_new;
      f_ = f_new;
      merit_ = merit_new;
      // Nielsen eq. (2.5): mu := mu max{1/gamma, 1 - (beta - 1)(2 rho - 1)^p}, nu := beta.
      const double t = 2.0 * rho - 1.0;
      mu *= std::max(1.0 / kGamma, 1.0 - (kBeta - 1.0) * t * t * t);
      nu = kBeta;
    } else {
      rejected();
    }
    // Step 5: F-convergence after every valid trial (MINPACK-1 section 2.3, LMDIF with squared
    // norms; without the ACTRED = -1 rule of LMDIF, ADR 0030 point 8).
    const bool merit_converged = merit_old > 0.0 &&
                                 std::abs(merit_old - merit_new) / merit_old <= options_.ftol &&
                                 predicted / merit_old <= options_.ftol && rho <= 2.0;
    // Step 6: a new Jacobian after an accepted step, then the gradient test.
    if (it.accepted && !merit_converged) {
      switch (jacobian()) {
        case JacobianOutcome::Cancelled:
          record();
          return finish(LmStatus::Cancelled);
        case JacobianOutcome::Failed:
          record();
          return finish(LmStatus::Failed);
        case JacobianOutcome::Ok:
          normal_equations();
          break;
      }
    }
    record();
    if (merit_converged) return finish(LmStatus::ConvergedMerit);
    if (it.accepted && gradient_converged()) return finish(LmStatus::ConvergedGradient);
  }
}

}  // namespace

LmResult levenberg_marquardt(const Residuals& residuals,
                             std::size_t m,
                             std::vector<double> p0,
                             std::vector<Bounds> bounds,
                             const LmOptions& options,
                             const trace::RunControl& control) {
  check_input(residuals, m, p0, bounds, options);
  Run run(residuals, m, std::move(p0), std::move(bounds), options, control);
  return run.run();
}

}  // namespace rtt::optim
