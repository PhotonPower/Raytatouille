#include "overlay_basis.hpp"

#include <cstddef>
#include <cstdint>
#include <exception>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/prescription.hpp"

namespace rtt::optim {

namespace {

compile::CompiledSystem compile_one(const model::System& system,
                                    const material::MaterialLibrary& materials,
                                    const coating::CoatingLibrary* coatings,
                                    std::size_t configuration) {
  return coatings != nullptr ? compile::compile(system, materials, *coatings, configuration)
                             : compile::compile(system, materials, configuration);
}

/// The last event of `path` that changes the ray at `wavelength`, by the rule of
/// rtt::paraxial::first_order (paraxial.hpp, "Vertices"): a reflection, a change of index or an
/// ideal lens; transmit events (stop, detector, thin elements without power) do not count.
std::optional<std::size_t> last_vertex(const compile::CompiledSystem& cs,
                                       compile::PathId path,
                                       std::uint16_t wavelength) {
  std::optional<std::size_t> vertex;
  const std::vector<compile::CompiledEvent>& events = cs.path(path).events;
  for (std::size_t e = 0; e < events.size(); ++e) {
    const compile::CompiledEvent& event = events[e];
    const double before = cs.media()[event.medium_before].index[wavelength].real();
    const double after = cs.media()[event.medium_after].index[wavelength].real();
    if (event.kind == model::EventKind::Reflect || after != before ||
        cs.surfaces()[event.surface].ideal_lens.has_value()) {
      vertex = e;
    }
  }
  return vertex;
}

std::string message_of(const std::exception_ptr& error) {
  try {
    std::rethrow_exception(error);
  } catch (const std::exception& e) {
    return e.what();
  } catch (...) {
    return "unknown error";
  }
}

}  // namespace

namespace detail {

model::System variant_of(const model::System& system) {
  model::System v = system;
  v.object.at_infinity = true;
  v.fields = {model::FieldType::AngleDeg, {{0.0, 0.0, 1.0}}};
  v.aperture = {model::SystemApertureType::StopSize, model::Param(0.0)};
  return v;
}

}  // namespace detail

std::size_t NominalBasis::configurations() const noexcept {
  return data_->configurations.size();
}

NominalBasis nominal_basis(const model::System& nominal,
                           const material::MaterialLibrary& materials,
                           const coating::CoatingLibrary* coatings) {
  auto data = std::make_shared<detail::BasisData>();
  data->wavelengths = nominal.wavelengths.size();
  const std::size_t count = model::configuration_count(nominal);
  const model::System variant = detail::variant_of(nominal);
  for (std::size_t k = 0; k < count; ++k) {
    const compile::CompiledSystem cs = compile_one(nominal, materials, coatings, k);
    detail::ConfigurationData config;
    config.object_at_infinity = nominal.object.at_infinity;
    // The variant of a finite object, compiled once per configuration; its failure only
    // concerns efl and bfl.
    std::optional<compile::CompiledSystem> inf;
    std::string inf_error;
    if (!config.object_at_infinity) {
      try {
        inf.emplace(compile_one(variant, materials, coatings, k));
      } catch (...) {
        inf_error = "the variant with the object at infinity does not compile: " +
                    message_of(std::current_exception());
      }
    }
    for (std::uint32_t p = 0; p < cs.paths().size(); ++p) {
      const compile::PathId path{p};
      detail::PathData pd{cs.path(path).name, analysis::nominal_reference(cs, path), {}};
      for (std::size_t w = 0; w < cs.wavelengths_um().size(); ++w) {
        const auto wl = static_cast<std::uint16_t>(w);
        detail::WavelengthData wd;
        try {
          const paraxial::FirstOrder fo = paraxial::first_order(cs, path, wl);
          wd.efl = fo.efl;
          wd.magnification = fo.lateral_magnification;
          wd.last_vertex = last_vertex(cs, path, wl);
          wd.f_number = paraxial::prescription(cs, path, wl).paraxial_working_f_number.has_value();
        } catch (...) {
          wd.unsupported = message_of(std::current_exception());
        }
        if (config.object_at_infinity) {
          try {
            static_cast<void>(pd.reference.aim(wl));
          } catch (...) {
            wd.infinite_unsupported = message_of(std::current_exception());
          }
        } else if (!inf) {
          wd.infinite_unsupported = inf_error;
        } else {
          try {
            wd.infinite_aim = trace::aim_reference(*inf, path, wl);
          } catch (...) {
            wd.infinite_unsupported = message_of(std::current_exception());
          }
        }
        pd.wavelengths.push_back(std::move(wd));
      }
      config.paths.push_back(std::move(pd));
    }
    data->configurations.push_back(std::move(config));
  }
  return NominalBasis(std::move(data));
}

}  // namespace rtt::optim
