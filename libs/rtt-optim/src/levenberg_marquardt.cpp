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
// MINPACK-1", ANL-80-74 (1980); for the equality constraints S. Boyd, L. Vandenberghe,
// "Introduction to Applied Linear Algebra" (VMLS), Cambridge 2018, sections 16.1-16.3, eqs.
// (16.1), (16.4), (16.5), (16.12). Decisions: ADR 0030 and ADR 0056. In comments, g alone means the
// constraints; the member g_ is J^T f (as without constraints).

namespace rtt::optim {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Nielsen 1999, eq. (2.5) with beta = nu_0 = 2, gamma = 3, p = 3 (the same as MNT eq. (2.21)).
constexpr double kBeta = 2.0;
constexpr double kGamma = 3.0;

// ADR 0056, point 3: at most this many Gauss-Newton corrections back onto g = 0 per trial.
constexpr int kMaxProjection = 10;

// Smallest pivot of the Cholesky factorization of the p x p constraint matrices, relative to the
// diagonal element: below it the rows of C count as linearly dependent (ADR 0056, point 2). Two
// equal rows leave a pivot of the order of the rounding, eps_M a_kk; 64 eps_M separates that from
// constraints whose gradients still differ by an angle of about 1e-7 rad. It covers rounding only:
// two constraints that are dependent analytically but computed differently differ in their
// central differences by the truncation error O(h^2) (h = 1e-3 max(|theta|, 1) with the
// function_precision of optimize), which may exceed it; the dependence then shows only when the
// steps have made the differences agree (review of #196).
constexpr double kRankTolerance = 64.0 * 2.220446049250313e-16;

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
  const std::size_t p = o.constraint_tolerances.size();
  if (p > 0) {
    if (p >= m) fail("at least one residual besides the constraints is needed");
    if (p > p0.size()) fail("at most as many constraints as variables are allowed");
    for (const double t : o.constraint_tolerances) {
      if (!finite(t) || !(t > 0.0)) fail("constraint tolerances must be finite and > 0");
    }
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
//
// With min_pivot > 0 a pivot d <= min_pivot a_kk (a_kk the diagonal element before the
// factorization) also counts as a failure: the rank test of the constraint matrices (ADR 0056).
bool cholesky_solve(std::vector<double> a,
                    std::size_t n,
                    std::span<const double> b,
                    std::span<double> x,
                    double min_pivot = 0.0) {
  // a is overwritten with C in its upper triangle.
  for (std::size_t k = 0; k < n; ++k) {
    const double akk = a[k * n + k];
    double d = akk;
    for (std::size_t i = 0; i < k; ++i) d -= a[i * n + k] * a[i * n + k];
    if (min_pivot > 0.0 ? !(d > min_pivot * akk) : !(d > 0.0)) return false;
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
      const trace::RunControl& control,
      const LmAcceptedCallback& on_accepted)
      : residuals_(residuals),
        m_(m),
        n_(p0.size()),
        p0_(std::move(p0)),
        bounds_(std::move(bounds)),
        options_(options),
        control_(control),
        on_accepted_(on_accepted),
        optimize_(control, static_cast<std::size_t>(options.max_iterations), "optimize"),
        p_(options.constraint_tolerances.size()),
        mf_(m - options.constraint_tolerances.size()) {
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
  void scaling();

  // --- Equality constraints (ADR 0056) ---------------------------------------------------------
  enum class Kkt : std::uint8_t { Ok, NotPositive, Dependent };
  enum class Projection : std::uint8_t { Feasible, Infeasible, Invalid, Dependent, Cancelled };
  LmResult run_constrained();
  [[nodiscard]] double objective(std::span<const double> f) const noexcept {
    return half_square(f.first(mf_));
  }
  [[nodiscard]] bool feasible(std::span<const double> f) const noexcept;
  [[nodiscard]] double violation(std::span<const double> f) const noexcept;
  [[nodiscard]] std::size_t worst_constraint(std::span<const double> f) const noexcept;
  [[nodiscard]] std::vector<double> constraint_metric() const;
  Kkt kkt_step(double mu, std::span<double> h, std::span<double> w) const;
  [[nodiscard]] bool correct(std::vector<double>& theta, std::span<const double> f) const;
  Projection project(std::vector<double>& theta, std::vector<double>& f);
  Projection project_start();
  [[nodiscard]] bool lagrangian_gradient_converged() const;
  LmResult finish(LmStatus status);
  bool report_accepted(const LmIteration& it);

  const Residuals& residuals_;
  std::size_t m_;
  std::size_t n_;
  std::vector<double> p0_;
  std::vector<Bounds> bounds_;
  const LmOptions& options_;
  const trace::RunControl& control_;
  const LmAcceptedCallback& on_accepted_;
  trace::RunMonitor optimize_;
  std::size_t p_;   // number of constraints (the last p_ residuals)
  std::size_t mf_;  // number of residuals of F (the first m_ - p_)

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
  std::optional<std::size_t> infeasible_constraint_;
  bool dependent_constraints_ = false;
  std::vector<LmIteration> history_;
};

// Central difference in theta (ADR 0007, ADR 0030 point 9), all 2n evaluations in parallel; one
// invalid side gives the one-sided difference with f(theta) (addendum #289, below).
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
    // Both monitors are finished, so an exception of the progress callback in either stage
    // propagates even when a cancellation was requested at the same time.
    const bool stage_cancelled = interrupted(monitor);
    const bool run_cancelled = interrupted(optimize_);
    if (stage_cancelled || run_cancelled) return JacobianOutcome::Cancelled;
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
    if (state[2 * j] != kValid && state[2 * j + 1] != kValid) {
      failed_variable_ = j;  // smallest index (ADR 0030, point 10; addendum #289)
      return JacobianOutcome::Failed;
    }
  }
  // Exactly one invalid side (ADR 0030, addendum #289): the one-sided difference with the stored
  // f(theta) of the last accepted state and the same h (MINPACK-1, section 2.4: the forward
  // difference; error O(h) instead of O(h^2)). The denominator is again the difference of the
  // two points as represented. The choice of the side depends only on the validity of the 2n
  // evaluations, each deterministic in its own slot.
  j_ = Columns{m_, n_, std::vector<double>(m_ * n_)};
  for (std::size_t j = 0; j < n_; ++j) {
    const bool plus_valid = state[2 * j] == kValid;
    const bool minus_valid = state[2 * j + 1] == kValid;
    const double* fp = plus_valid ? values.data() + 2 * j * m_ : f_.data();
    const double* fm = minus_valid ? values.data() + (2 * j + 1) * m_ : f_.data();
    const double denominator =
        (plus_valid ? plus[j] : theta_[j]) - (minus_valid ? minus[j] : theta_[j]);
    for (std::size_t i = 0; i < m_; ++i) j_(i, j) = (fp[i] - fm[i]) / denominator;
  }
  return JacobianOutcome::Ok;
}

// Scaling D: column norms of the start Jacobian (all m rows, the constraints included), 1 for a
// zero column (MINPACK-1, section 2.5, MODE = 1); fixed for the run (ADR 0030, point 6).
void Run::scaling() {
  d_.resize(n_);
  for (std::size_t j = 0; j < n_; ++j) {
    double sum = 0.0;
    for (std::size_t i = 0; i < m_; ++i) sum += j_(i, j) * j_(i, j);
    d_[j] = sum > 0.0 ? std::sqrt(sum) : 1.0;
  }
}

// A = J^T J and g = J^T f, summed in index order, over the residuals of F (all m without
// constraints; the rows of C are not part of them, ADR 0056).
void Run::normal_equations() {
  a_.assign(n_ * n_, 0.0);
  g_.assign(n_, 0.0);
  for (std::size_t r = 0; r < n_; ++r) {
    double gr = 0.0;
    for (std::size_t i = 0; i < mf_; ++i) gr += j_(i, r) * f_[i];
    g_[r] = gr;
    for (std::size_t c = r; c < n_; ++c) {
      double sum = 0.0;
      for (std::size_t i = 0; i < mf_; ++i) sum += j_(i, r) * j_(i, c);
      a_[r * n_ + c] = sum;
      a_[c * n_ + r] = sum;
    }
  }
}

// Reports the accepted state to on_accepted_ (ADR 0030, addendum #378); false if the callback
// threw, which ends the run like a cancellation. Sequential: called from the solver loop only.
bool Run::report_accepted(const LmIteration& it) {
  if (!on_accepted_) return true;
  std::vector<double> p(n_);
  external(theta_, p);
  std::vector<std::uint8_t> changed(n_);
  for (std::size_t j = 0; j < n_; ++j) changed[j] = same_bits(theta_[j], theta0_[j]) ? 0 : 1;
  try {
    on_accepted_(it, p, changed);
  } catch (...) {
    return false;  // decision of #378: the run ends as Cancelled, the exception is not rethrown
  }
  return true;
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
  r.infeasible_constraint = infeasible_constraint_;
  r.dependent_constraints = dependent_constraints_;
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
  if (p_ > 0) return run_constrained();
  // Step 1: start.
  if (interrupted(optimize_)) return finish(LmStatus::Cancelled);
  f_.resize(m_);
  ++evaluations_;
  if (!evaluate(p0_, f_)) {
    throw std::invalid_argument("levenberg_marquardt: residuals invalid at the start point");
  }
  merit_ = half_square(f_);
  // f = 0 already: g = J^T f = 0 for every J, so the gradient test holds without differences
  // (ADR 0030, point 9).
  if (std::all_of(f_.begin(), f_.end(), [](double v) { return v == 0.0; })) {
    return finish(LmStatus::ConvergedGradient);
  }
  switch (jacobian()) {
    case JacobianOutcome::Cancelled:
      return finish(LmStatus::Cancelled);
    case JacobianOutcome::Failed:
      return finish(LmStatus::Failed);
    case JacobianOutcome::Ok:
      break;
  }
  scaling();
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
      // The live state (#378), before the Jacobian of the new state. `it` carries F and the
      // evaluations up to here; record() below sets the same F (merit_ does not change).
      it.F = merit_;
      it.evaluations = evaluations_;
      if (!report_accepted(it)) {
        record();
        return finish(LmStatus::Cancelled);
      }
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

// --- Equality constraints (ADR 0056) -----------------------------------------------------------
// Notation: f the first mf_ residuals (rows 0..mf_-1 of the Jacobian, J), g the last p_ residuals
// (rows mf_..m_-1, C), in the internal variables theta; z = D theta are the scaled variables.

bool Run::feasible(std::span<const double> f) const noexcept {
  for (std::size_t i = 0; i < p_; ++i) {
    if (!(std::abs(f[mf_ + i]) <= options_.constraint_tolerances[i])) return false;
  }
  return true;
}

// max_i |g_i| / tolerance_i: <= 1 exactly when feasible.
double Run::violation(std::span<const double> f) const noexcept {
  double v = 0.0;
  for (std::size_t i = 0; i < p_; ++i) {
    v = std::max(v, std::abs(f[mf_ + i]) / options_.constraint_tolerances[i]);
  }
  return v;
}

// The first index of the largest |g_i| / tolerance_i.
std::size_t Run::worst_constraint(std::span<const double> f) const noexcept {
  std::size_t worst = 0;
  double v = -1.0;
  for (std::size_t i = 0; i < p_; ++i) {
    const double vi = std::abs(f[mf_ + i]) / options_.constraint_tolerances[i];
    if (vi > v) {
      v = vi;
      worst = i;
    }
  }
  return worst;
}

// M = C D^-2 C^T (p x p, row-major), summed in index order: the Gram matrix of the rows of C in the
// scaled variables, C_z = C D^-1.
std::vector<double> Run::constraint_metric() const {
  std::vector<double> m(p_ * p_);
  for (std::size_t r = 0; r < p_; ++r) {
    for (std::size_t c = r; c < p_; ++c) {
      double sum = 0.0;
      for (std::size_t j = 0; j < n_; ++j) {
        sum += j_(mf_ + r, j) * j_(mf_ + c, j) / (d_[j] * d_[j]);
      }
      m[r * p_ + c] = sum;
      m[c * p_ + r] = sum;
    }
  }
  return m;
}

// The damped step (ADR 0056, point 2): h minimizes ||J h + f||^2 + mu ||D h||^2 subject to
// C h = -g. That is VMLS eq. (16.1) with A = [J; sqrt(mu) D], b = [-f; 0] and d = -g; its KKT
// equations (16.4), with the first block row divided by 2 and w = z/2, are
//   (J^T J + mu D^T D) h + C^T w = -J^T f,   C h = -g.
// Solved by block elimination (derivation, docs/quellen.md "hergeleitet"): with H = J^T J + mu D^T
// D (positive definite for mu > 0), x0 = H^-1 (-J^T f) and Y = H^-1 C^T the first row gives
// h = x0 - Y w, and the second S w = C x0 + g with S = C Y = C H^-1 C^T, which is positive
// definite exactly when C has independent rows (VMLS (16.5); the second condition holds because
// of mu D^T D). Cholesky (MNT Algorithm A.4) for H and S; a failure of H is a rejected trial as
// without constraints, a failure of S (rank test kRankTolerance) means dependent constraints.
// w are the multipliers of the Lagrangian F + w^T g: at a solution (h -> 0) J^T f + C^T w = 0.
Run::Kkt Run::kkt_step(double mu, std::span<double> h, std::span<double> w) const {
  std::vector<double> hm = a_;
  for (std::size_t j = 0; j < n_; ++j) hm[j * n_ + j] += mu * d_[j] * d_[j];
  std::vector<double> rhs(n_);
  std::vector<double> x0(n_);
  for (std::size_t j = 0; j < n_; ++j) rhs[j] = -g_[j];
  if (!cholesky_solve(hm, n_, rhs, x0)) return Kkt::NotPositive;
  std::vector<double> y(n_ * p_);  // column c of Y at y[c * n_ ...]
  for (std::size_t c = 0; c < p_; ++c) {
    for (std::size_t j = 0; j < n_; ++j) rhs[j] = j_(mf_ + c, j);
    if (!cholesky_solve(hm, n_, rhs, std::span<double>(y.data() + c * n_, n_))) {
      return Kkt::NotPositive;
    }
  }
  std::vector<double> s(p_ * p_);
  for (std::size_t r = 0; r < p_; ++r) {
    for (std::size_t c = r; c < p_; ++c) {
      // S is symmetric in exact arithmetic; the mean of both products keeps it so after rounding.
      double rc = 0.0;
      double cr = 0.0;
      for (std::size_t j = 0; j < n_; ++j) {
        rc += j_(mf_ + r, j) * y[c * n_ + j];
        cr += j_(mf_ + c, j) * y[r * n_ + j];
      }
      s[r * p_ + c] = 0.5 * (rc + cr);
      s[c * p_ + r] = s[r * p_ + c];
    }
  }
  std::vector<double> b(p_);
  for (std::size_t r = 0; r < p_; ++r) {
    double sum = f_[mf_ + r];
    for (std::size_t j = 0; j < n_; ++j) sum += j_(mf_ + r, j) * x0[j];
    b[r] = sum;
  }
  if (!cholesky_solve(std::move(s), p_, b, w, kRankTolerance)) return Kkt::Dependent;
  for (std::size_t j = 0; j < n_; ++j) {
    double sum = x0[j];
    for (std::size_t c = 0; c < p_; ++c) sum -= y[c * n_ + j] * w[c];
    h[j] = sum;
  }
  return std::all_of(h.begin(), h.end(), finite) ? Kkt::Ok : Kkt::NotPositive;
}

// One correction back onto g = 0 (ADR 0056, point 3) with the matrix C of the current Jacobian:
// the least-norm solution of the linearized constraints C Delta = -g in the scaled variables,
// VMLS eq. (16.12) x = C^T (C C^T)^-1 d with C_z = C D^-1 and d = -g:
//   Delta theta = D^-1 C_z^T (C_z C_z^T)^-1 (-g) = D^-2 C^T M^-1 (-g),  M = C D^-2 C^T.
// False if M is singular (dependent constraints, rank test kRankTolerance).
bool Run::correct(std::vector<double>& theta, std::span<const double> f) const {
  std::vector<double> minus_g(p_);
  std::vector<double> u(p_);
  for (std::size_t i = 0; i < p_; ++i) minus_g[i] = -f[mf_ + i];
  if (!cholesky_solve(constraint_metric(), p_, minus_g, u, kRankTolerance)) return false;
  for (std::size_t j = 0; j < n_; ++j) {
    double sum = 0.0;
    for (std::size_t i = 0; i < p_; ++i) sum += j_(mf_ + i, j) * u[i];
    theta[j] += sum / (d_[j] * d_[j]);
  }
  return true;
}

// Projection of a trial (ADR 0056, point 3): from theta (with its residuals f) the simplified
// Gauss-Newton corrections with the matrix C of the last accepted state, until every
// |g_i| <= tolerance_i, at most kMaxProjection corrections. After the KKT step the violation is
// of second order in h, so the fixed C suffices. An increase of the violation (no convergence)
// or an invalid evaluation ends it without success.
Run::Projection Run::project(std::vector<double>& theta, std::vector<double>& f) {
  std::vector<double> p(n_);
  double last = violation(f);
  for (int step = 0; step < kMaxProjection; ++step) {
    if (last <= 1.0) return Projection::Feasible;
    if (!correct(theta, f)) return Projection::Dependent;
    if (interrupted(optimize_)) return Projection::Cancelled;
    external(theta, p);
    ++evaluations_;
    if (!evaluate(p, f)) {
      ++failed_evaluations_;
      return Projection::Invalid;
    }
    const double v = violation(f);
    if (!(v < last)) return Projection::Infeasible;  // v >= last > 1: no convergence
    last = v;
  }
  return last <= 1.0 ? Projection::Feasible : Projection::Infeasible;
}

// Projection of an infeasible start (ADR 0056, point 5, as decided on #196): Gauss-Newton with a
// new Jacobian before every correction after the first (the start Jacobian is current), so the
// convergence is quadratic also far from g = 0, at most kMaxProjection corrections. Works on
// theta_ and f_; without success they are reset to the start.
Run::Projection Run::project_start() {
  const std::vector<double> f_start = f_;
  const auto reset = [&](Projection outcome) {
    theta_ = theta0_;
    f_ = f_start;
    return outcome;
  };
  std::vector<double> p(n_);
  double last = violation(f_);
  for (int step = 0; step < kMaxProjection; ++step) {
    if (last <= 1.0) return Projection::Feasible;
    if (step > 0) {
      switch (jacobian()) {
        case JacobianOutcome::Cancelled:
          return reset(Projection::Cancelled);
        case JacobianOutcome::Failed:
          failed_variable_.reset();  // reported as an infeasible start, not a failed column
          return reset(Projection::Invalid);
        case JacobianOutcome::Ok:
          break;
      }
    }
    if (!correct(theta_, f_)) return reset(Projection::Dependent);
    if (interrupted(optimize_)) return reset(Projection::Cancelled);
    external(theta_, p);
    ++evaluations_;
    if (!evaluate(p, f_)) {
      ++failed_evaluations_;
      return reset(Projection::Invalid);
    }
    const double v = violation(f_);
    if (!(v < last)) return reset(Projection::Infeasible);  // v >= last > 1: no convergence
    last = v;
  }
  return last <= 1.0 ? Projection::Feasible : reset(Projection::Infeasible);
}

// The gradient test with the gradient of the Lagrangian (ADR 0056, point 4):
// ||D^-1 (J^T f + C^T w)||_inf <= gtol with the multipliers w that minimize that norm in the
// 2-norm, w = -M^-1 C D^-2 J^T f (the normal equations of min ||D^-1 (J^T f + C^T w)||). False if
// M is singular (the KKT step reports that).
bool Run::lagrangian_gradient_converged() const {
  std::vector<double> rhs(p_);
  for (std::size_t i = 0; i < p_; ++i) {
    double sum = 0.0;
    for (std::size_t j = 0; j < n_; ++j) sum += j_(mf_ + i, j) * g_[j] / (d_[j] * d_[j]);
    rhs[i] = -sum;
  }
  std::vector<double> w(p_);
  if (!cholesky_solve(constraint_metric(), p_, rhs, w, kRankTolerance)) return false;
  double norm = 0.0;
  for (std::size_t j = 0; j < n_; ++j) {
    double sum = g_[j];
    for (std::size_t i = 0; i < p_; ++i) sum += j_(mf_ + i, j) * w[i];
    norm = std::max(norm, std::abs(sum) / d_[j]);
  }
  return norm <= options_.gtol;
}

// The run with constraints (ADR 0056): MNT Algorithm 3.16 as run(), with the KKT step instead of
// the damped normal equations, the projection of every trial onto g = 0 before the gain ratio and
// F = 1/2 ||f||^2 over the residuals of F only.
LmResult Run::run_constrained() {
  // Step 1: start.
  if (interrupted(optimize_)) return finish(LmStatus::Cancelled);
  f_.resize(m_);
  ++evaluations_;
  if (!evaluate(p0_, f_)) {
    throw std::invalid_argument("levenberg_marquardt: residuals invalid at the start point");
  }
  merit_ = objective(f_);
  if (feasible(f_) && std::all_of(f_.begin(), f_.begin() + static_cast<std::ptrdiff_t>(mf_),
                                  [](double v) { return v == 0.0; })) {
    return finish(LmStatus::ConvergedGradient);
  }
  const auto new_jacobian = [this]() -> std::optional<LmStatus> {
    switch (jacobian()) {
      case JacobianOutcome::Cancelled:
        return LmStatus::Cancelled;
      case JacobianOutcome::Failed:
        return LmStatus::Failed;
      case JacobianOutcome::Ok:
        break;
    }
    return std::nullopt;
  };
  if (const auto end = new_jacobian()) return finish(*end);
  scaling();
  // An infeasible start is projected first (ADR 0056, point 5).
  if (!feasible(f_)) {
    switch (project_start()) {
      case Projection::Feasible:
        break;
      case Projection::Cancelled:
        return finish(LmStatus::Cancelled);
      case Projection::Dependent:
        dependent_constraints_ = true;
        return finish(LmStatus::Failed);
      case Projection::Infeasible:
      case Projection::Invalid:
        infeasible_constraint_ = worst_constraint(f_);
        return finish(LmStatus::Failed);
    }
    merit_ = objective(f_);
    if (const auto end = new_jacobian()) return finish(*end);
  }
  normal_equations();
  if (lagrangian_gradient_converged()) return finish(LmStatus::ConvergedGradient);
  // mu_0 = tau max_i (J_z^T J_z)_ii (MNT eq. (3.14)) over the rows of F; tau if they are all zero
  // in the scaled variables (only the constraints move then).
  double max_diagonal = 0.0;
  for (std::size_t j = 0; j < n_; ++j) {
    max_diagonal = std::max(max_diagonal, a_[j * n_ + j] / (d_[j] * d_[j]));
  }
  double mu = options_.tau * (max_diagonal > 0.0 ? max_diagonal : 1.0);
  double nu = kBeta;

  std::vector<double> h(n_);
  std::vector<double> w(p_);
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
    switch (kkt_step(mu, h, w)) {
      case Kkt::Ok:
        break;
      case Kkt::NotPositive:  // a rejected trial without an evaluation
        rejected();
        record();
        continue;
      case Kkt::Dependent:
        dependent_constraints_ = true;
        record();
        return finish(LmStatus::Failed);
    }
    // Step 3: the step test (MNT eq. (3.15b) in z); h is not applied.
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
    // Step 4: the trial and its projection onto g = 0, after the cancellation check.
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
    bool valid = true;
    switch (project(theta_new, f_new)) {
      case Projection::Feasible:
        break;
      case Projection::Cancelled:
        it.F = merit_;
        it.evaluations = evaluations_;  // the trial and corrections so far (cumulative)
        history_.push_back(it);
        return finish(LmStatus::Cancelled);
      case Projection::Dependent:
        dependent_constraints_ = true;
        record();
        return finish(LmStatus::Failed);
      case Projection::Infeasible:
      case Projection::Invalid:
        valid = false;
        break;
    }
    if (!valid) {
      rejected();
      record();
      continue;
    }
    const double merit_new = objective(f_new);
    // Gain ratio on the Lagrangian l = F + w^T g with the multipliers w of the KKT step (ADR 0056,
    // point 4). Both states are feasible only up to the tolerances, and inside that band F
    // changes along the normal of the constraints by about |f| tolerance / |grad g|, to first
    // order; near a constrained minimum with a large residual that is as large as the true gain,
    // and F alone would accept or reject by chance. In l that first-order change cancels
    // (grad F + C^T w = 0 at the solution). The model is L(h) = 1/2 ||f + J h||^2 with the
    // linearized constraint g + C h = 0 (MNT eq. (2.18)). With the first KKT row,
    // -h^T J^T f = ||J h||^2 + mu ||D h||^2 + (C h)^T w = ||J h||^2 + mu ||D h||^2 - g^T w, so
    //   predicted = L(0) - L(h) + w^T g = -h^T J^T f - 1/2 ||J h||^2 + w^T g
    //             = 1/2 ||J h||^2 + mu ||D h||^2,
    // positive and without cancellation (as 1/2 h^T (mu D^T D h - J^T f) without constraints,
    // MNT p. 25), and
    //   actual    = F - F_new + w^T (g - g_new)
    // without the difference of the large values F and F_new (F changes near the solution only to
    // second order, below the rounding of F): by a^2 - b^2 = (a - b)(a + b),
    // F - F_new = 1/2 sum_i (f_i - f_new_i)(f_i + f_new_i). That removes the rounding of the
    // subtraction F - F_new (eps_M F), not the error of the evaluations themselves: residuals
    // known to delta_i leave an error of about sum_i |f_i| |delta_i - delta_new_i| (with the
    // function_precision of optimize about 1e-9 F); the essential gain is the cancellation of
    // the normal part in l. A non-positive predicted value (rounding) counts as no predicted
    // gain.
    double predicted = 0.0;
    for (std::size_t i = 0; i < mf_; ++i) {
      double r = 0.0;
      for (std::size_t j = 0; j < n_; ++j) r += j_(i, j) * h[j];
      predicted += r * r;
    }
    predicted *= 0.5;
    for (std::size_t j = 0; j < n_; ++j) predicted += mu * d_[j] * d_[j] * h[j] * h[j];
    double actual = 0.0;
    for (std::size_t i = 0; i < mf_; ++i) actual += (f_[i] - f_new[i]) * (f_[i] + f_new[i]);
    actual *= 0.5;
    for (std::size_t i = 0; i < p_; ++i) actual += w[i] * (f_[mf_ + i] - f_new[mf_ + i]);
    const double merit_old = merit_;
    const double rho = predicted > 0.0 ? actual / predicted : 0.0;
    it.rho = rho;
    it.accepted = rho > 0.0;
    if (it.accepted) {
      theta_ = theta_new;
      f_ = f_new;
      merit_ = merit_new;
      const double t = 2.0 * rho - 1.0;
      mu *= std::max(1.0 / kGamma, 1.0 - (kBeta - 1.0) * t * t * t);
      nu = kBeta;
      it.F = merit_;
      it.evaluations = evaluations_;
      if (!report_accepted(it)) {
        record();
        return finish(LmStatus::Cancelled);
      }
    } else {
      rejected();
    }
    // Step 5: F-convergence after every valid trial, as run(), with the actual decrease above.
    const bool merit_converged = merit_old > 0.0 && std::abs(actual) / merit_old <= options_.ftol &&
                                 predicted / merit_old <= options_.ftol && rho <= 2.0;
    // Step 6: a new Jacobian after an accepted step, then the gradient test of the Lagrangian.
    if (it.accepted && !merit_converged) {
      if (const auto end = new_jacobian()) {
        record();
        return finish(*end);
      }
      normal_equations();
    }
    record();
    if (merit_converged) return finish(LmStatus::ConvergedMerit);
    if (it.accepted && lagrangian_gradient_converged()) {
      return finish(LmStatus::ConvergedGradient);
    }
  }
}

}  // namespace

LmResult levenberg_marquardt(const Residuals& residuals,
                             std::size_t m,
                             std::vector<double> p0,
                             std::vector<Bounds> bounds,
                             const LmOptions& options,
                             const trace::RunControl& control,
                             const LmAcceptedCallback& on_accepted) {
  check_input(residuals, m, p0, bounds, options);
  Run run(residuals, m, std::move(p0), std::move(bounds), options, control, on_accepted);
  return run.run();
}

}  // namespace rtt::optim
