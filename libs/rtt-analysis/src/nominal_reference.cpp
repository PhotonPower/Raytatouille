#include "rtt/analysis/nominal_reference.hpp"

#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

#include "rtt/compile/errors.hpp"
#include "rtt/paraxial/paraxial.hpp"

namespace rtt::analysis {

namespace detail {

/// Data of a NominalReference: per wavelength the aiming data or the exception of their
/// computation, and the exit pupil or its exception.
struct NominalData {
  compile::PathId path;
  std::vector<std::optional<trace::AimReference>> aim;
  std::vector<std::exception_ptr> aim_error;
  std::optional<math::Vec3> exit_pupil;
  std::exception_ptr exit_pupil_error;
};

}  // namespace detail

compile::PathId NominalReference::path() const noexcept {
  return data_->path;
}

const trace::AimReference& NominalReference::aim(std::uint16_t wavelength) const {
  if (wavelength >= data_->aim.size()) {
    throw std::invalid_argument("analysis: wavelength index " + std::to_string(wavelength) +
                                " does not exist in the nominal reference");
  }
  const std::optional<trace::AimReference>& aim = data_->aim[wavelength];
  if (!aim.has_value()) std::rethrow_exception(data_->aim_error[wavelength]);
  return aim.value();
}

std::optional<math::Vec3> NominalReference::exit_pupil() const {
  if (data_->exit_pupil_error) std::rethrow_exception(data_->exit_pupil_error);
  return data_->exit_pupil;
}

NominalReference nominal_reference(const compile::CompiledSystem& nominal, compile::PathId path) {
  if (path.index >= nominal.paths().size()) {
    throw std::invalid_argument("analysis: path index " + std::to_string(path.index) +
                                " does not exist");
  }
  auto data = std::make_shared<detail::NominalData>();
  data->path = path;
  const std::size_t count = nominal.wavelengths_um().size();
  data->aim.resize(count);
  data->aim_error.resize(count);
  for (std::size_t w = 0; w < count; ++w) {
    try {
      data->aim[w] = trace::aim_reference(nominal, path, static_cast<std::uint16_t>(w));
    } catch (...) {
      data->aim_error[w] = std::current_exception();
    }
  }
  // The exit pupil as the reference sphere without reference takes it (make_reference in
  // reference.cpp): first_order at the reference wavelength, E = (0, 0, z_XP).
  try {
    compile::require_stop(nominal, path);
    const paraxial::FirstOrder fo =
        paraxial::first_order(nominal, path, nominal.reference_wavelength());
    // With a stop on the path first_order() always gives the exit pupil.
    if (!fo.exit_pupil) throw std::logic_error("analysis: no exit pupil despite a stop");
    const std::optional<double>& z = fo.exit_pupil->z;
    if (z.has_value()) data->exit_pupil = math::Vec3(0.0, 0.0, z.value());
  } catch (...) {
    data->exit_pupil_error = std::current_exception();
  }
  return NominalReference(std::move(data));
}

}  // namespace rtt::analysis
