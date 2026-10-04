#include "rtt/material/material.hpp"

#include <charconv>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <system_error>
#include <utility>

namespace rtt::material {
namespace {

constexpr std::string_view kConstPrefix = "CONST:";

/// Medium with a fixed complex index n + i*kappa, independent of wavelength and temperature.
/// Sign convention n + i*kappa, kappa >= 0 for absorption: docs/architecture.md, "Konventionen".
class ConstantMaterial final : public Material {
 public:
  explicit ConstantMaterial(math::Complex index) : index_(index) {}

  [[nodiscard]] math::Complex index(double /*wavelength_um*/,
                                    double /*temperature_c*/,
                                    double /*pressure_atm*/) const override {
    return index_;
  }

 private:
  math::Complex index_;
};

/// Parses the whole of `text` as a finite double in C-locale syntax (no whitespace, no '+').
std::optional<double> parse_finite(std::string_view text) {
  double value = 0.0;
  const char* const first = text.data();
  const char* const last = first + text.size();
  const auto [ptr, ec] = std::from_chars(first, last, value);
  if (ec != std::errc{} || ptr != last || !std::isfinite(value)) return std::nullopt;
  return value;
}

[[noreturn]] void reject(std::string_view reference, std::string_view reason) {
  throw UnknownMaterial("invalid material reference '" + std::string(reference) +
                        "': " + std::string(reason));
}

/// Parses "<n>" or "<n>,<kappa>" (the part after "CONST:").
math::Complex parse_constant(std::string_view reference, std::string_view args) {
  const auto comma = args.find(',');
  const std::string_view n_text = args.substr(0, comma);
  const std::optional<double> n = parse_finite(n_text);
  if (!n) reject(reference, "expected CONST:<n> or CONST:<n>,<kappa> with finite numbers");
  if (*n <= 0.0) reject(reference, "refractive index n must be positive");
  if (comma == std::string_view::npos) return {*n, 0.0};

  const std::optional<double> kappa = parse_finite(args.substr(comma + 1));
  if (!kappa) reject(reference, "expected CONST:<n> or CONST:<n>,<kappa> with finite numbers");
  if (*kappa < 0.0) reject(reference, "extinction coefficient kappa must be >= 0");
  // Adding +0.0 maps -0.0 to +0.0 (IEEE 754-2019, 6.3): the sign of a zero imaginary part
  // selects the branch of complex square roots downstream, so kappa >= 0 includes the sign.
  return {*n, *kappa + 0.0};
}

std::shared_ptr<const Material> make_material(std::string_view reference) {
  if (reference == "VACUUM") return std::make_shared<const ConstantMaterial>(1.0);
  // M1 placeholder: air as n = 1. Ciddor air and glass data relative to air follow in M2.
  if (reference == "AIR") return std::make_shared<const ConstantMaterial>(1.0);
  if (reference.starts_with(kConstPrefix)) {
    return std::make_shared<const ConstantMaterial>(
        parse_constant(reference, reference.substr(kConstPrefix.size())));
  }
  throw UnknownMaterial(
      "unknown material reference '" + std::string(reference) +
      "': supported are VACUUM, AIR, CONST:<n>, CONST:<n>,<kappa> and "
      "materials registered with MaterialLibrary::add; glass catalogs follow in M2");
}

}  // namespace

std::shared_ptr<const Material> MaterialLibrary::resolve(std::string_view reference) const {
  const std::scoped_lock lock(mutex_);
  if (const auto it = cache_.find(reference); it != cache_.end()) return it->second;
  auto material = make_material(reference);
  cache_.emplace(std::string(reference), material);
  return material;
}

void MaterialLibrary::add(std::string name, std::shared_ptr<const Material> material) {
  if (name.empty()) throw std::invalid_argument("MaterialLibrary::add: empty material name");
  if (name == "VACUUM" || name == "AIR" || name.starts_with(kConstPrefix)) {
    throw std::invalid_argument("MaterialLibrary::add: '" + name + "' is a reserved name");
  }
  if (material == nullptr) {
    throw std::invalid_argument("MaterialLibrary::add: material '" + name + "' is null");
  }
  const std::scoped_lock lock(mutex_);
  if (cache_.contains(name)) {
    throw std::invalid_argument("MaterialLibrary::add: '" + name + "' is already in use");
  }
  cache_.emplace(std::move(name), std::move(material));
}

}  // namespace rtt::material
