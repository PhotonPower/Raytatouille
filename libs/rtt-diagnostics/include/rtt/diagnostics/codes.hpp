#pragma once

/// @file codes.hpp
/// Registry of the stable diagnostic codes (ADR 0022), layer Basis: every library can name a
/// code without depending on rtt-model.
///
/// A code is a lower-case dotted name "<group>.<what>", e.g. "material.unknown". Codes are
/// stable across versions: a code is never renamed, reused or given another meaning; new codes
/// are only added. An interim code added between two releases (its summary says "interim") may
/// be removed again before the next release, because it was never published as stable (ADR
/// 0022). Every code is listed in docs/diagnostics.md.

#include <array>
#include <cstdint>
#include <stdexcept>
#include <string_view>

namespace rtt::diagnostics {

/// Severity of a diagnostic. Errors stop compile(); warnings are reported as data.
enum class Severity : std::uint8_t { Error, Warning };

/// One entry of the registry.
struct CodeInfo {
  std::string_view code;      ///< dotted name, e.g. "material.unknown"
  Severity severity;          ///< the severity every diagnostic with this code has
  std::string_view producer;  ///< "validate", "compile", "analysis", "agf", "edit" or "io"
                              ///< (docs/diagnostics.md)
  std::string_view summary;   ///< one-line meaning, as in docs/diagnostics.md
};

/// All registered codes, sorted by code. Iterable at run time (documentation check, Python).
inline constexpr std::array kCodes = std::to_array<CodeInfo>({
    {"agf.duplicate_glass", Severity::Warning, "agf",
     "glass repeated with identical data; the repetition is ignored"},
    {"agf.duplicate_glass_conflict", Severity::Warning, "agf",
     "glass defined twice with different data; it is ambiguous and cannot be resolved"},
    {"agf.preamble_skipped", Severity::Warning, "agf", "text line before the first record skipped"},
    {"agf.stray_line", Severity::Warning, "agf", "single word before the next NM record skipped"},
    {"aperture.na_not_physical", Severity::Warning, "validate", "object-space NA >= 1 in air"},
    {"aperture.stop_missing", Severity::Error, "validate",
     "aperture type stop_size without a stop element"},
    {"aperture.value_invalid", Severity::Error, "validate",
     "system aperture value not finite or <= 0"},
    {"bounds.invalid", Severity::Error, "validate",
     "min not below max, or bounds at a derived parameter row"},
    {"bounds.value_outside", Severity::Warning, "validate",
     "value of a Param or row outside its bounds"},
    {"coating.design_wavelength_out_of_range", Severity::Error, "compile",
     "QWOT design wavelength outside the valid range of the layer material"},
    {"coating.layer_material_unknown", Severity::Error, "compile",
     "material of a coating layer cannot be resolved"},
    {"coating.no_substrate", Severity::Error, "compile",
     "coating on a lens, plate or mirror without material"},
    {"coating.not_allowed", Severity::Error, "compile",
     "coating on a thin element, stop or detector (no substrate)"},
    {"coating.substrate_ambiguous", Severity::Error, "compile",
     "coating on an inner surface between two segments"},
    {"coating.thickness_invalid", Severity::Error, "compile", "layer thickness cannot be computed"},
    {"coating.unknown", Severity::Error, "compile", "coating reference cannot be resolved"},
    {"coating.wavelength_out_of_range", Severity::Error, "compile",
     "system wavelength outside the valid range of a layer material"},
    {"configurations.name_duplicate", Severity::Error, "validate", "configuration name used twice"},
    {"configurations.name_invalid", Severity::Error, "validate",
     "empty or blank configuration name"},
    {"crystal.absorbing", Severity::Error, "compile",
     "crystal part with kappa != 0 at a system wavelength (not supported in M4)"},
    {"crystal.interaction_unsupported", Severity::Error, "compile",
     "interaction other than fresnel or ideal_anti_reflection at a crystal surface"},
    {"crystal.kind_not_allowed", Severity::Error, "validate",
     "crystal material at an element that is not a lens or a plate"},
    {"crystal.material_conflict", Severity::Error, "validate",
     "crystal and isotropic material both set (API only)"},
    {"crystal.optic_axis_invalid", Severity::Error, "validate", "optic axis zero or not finite"},
    {"crystal.optic_axis_missing", Severity::Error, "validate", "crystal without optic axis"},
    {"crystal.optic_axis_not_allowed", Severity::Error, "validate",
     "optic axis at an element without crystal"},
    {"crystal.unsupported", Severity::Error, "compile",
     "crystal case that cannot be traced yet (ADR 0026, point 4)"},
    {"edit.base_not_representable", Severity::Error, "edit",
     "the system to patch cannot be written in the edit form (non-finite number, material and "
     "segment list or crystal both set, invalid UTF-8)"},
    {"edit.invalid_value", Severity::Error, "edit",
     "the patched edit form is not a valid system file (type, key, enum value)"},
    {"edit.patch_invalid", Severity::Error, "edit",
     "not a valid JSON Patch document or operation (RFC 6902)"},
    {"edit.path_not_found", Severity::Error, "edit",
     "patch target, source or parent missing, or the parent is no object or array"},
    {"edit.read_only", Severity::Error, "edit",
     "patch writes to /schema_version, /units or the whole document"},
    {"edit.test_failed", Severity::Error, "edit", "a test operation of the patch failed"},
    {"element.material_both", Severity::Error, "validate",
     "one material for all segments and a segment list at the same time"},
    {"element.material_empty", Severity::Error, "validate", "empty material reference"},
    {"element.material_list_not_allowed", Severity::Error, "validate",
     "mirror with a list of segment materials"},
    {"element.material_missing", Severity::Error, "validate", "lens or plate without material"},
    {"element.material_not_allowed", Severity::Error, "validate",
     "thin element, stop or detector with a material"},
    {"element.plate_surface_not_plane", Severity::Error, "validate",
     "plate surface is not a plane"},
    {"element.segment_count", Severity::Error, "validate",
     "number of segment materials does not match the surfaces"},
    {"element.surface_count", Severity::Error, "validate",
     "wrong number of surfaces for the element kind"},
    {"environment.medium_empty", Severity::Error, "validate", "empty environment medium"},
    {"environment.pressure_invalid", Severity::Error, "validate", "pressure not finite or < 0 atm"},
    {"environment.temperature_invalid", Severity::Error, "validate",
     "temperature not finite or not above absolute zero"},
    {"fields.coordinate_invalid", Severity::Error, "validate", "field coordinate not finite"},
    {"fields.empty", Severity::Error, "validate", "no field point"},
    {"fields.weight_invalid", Severity::Error, "validate", "field weight not finite or < 0"},
    {"interaction.axis_invalid", Severity::Error, "validate",
     "polarizer or retarder axis zero or not finite"},
    {"interaction.coating_name_empty", Severity::Error, "validate", "empty coating reference"},
    {"interaction.extinction_ratio_invalid", Severity::Error, "validate",
     "polarizer extinction ratio outside [0, 1]"},
    {"interaction.reflectance_invalid", Severity::Error, "validate",
     "beam splitter reflectance outside [0, 1]"},
    {"interaction.retardance_invalid", Severity::Error, "validate", "retardance not finite"},
    {"io.pickup_dropped", Severity::Warning, "io",
     "pickup of a file before schema 0.4 dropped on reading"},
    {"material.unknown", Severity::Error, "compile", "material reference cannot be resolved"},
    {"material.wavelength_out_of_range", Severity::Error, "compile",
     "system wavelength outside the valid range of a material on a path"},
    {"node.name_duplicate", Severity::Error, "validate", "assembly or element name used twice"},
    {"node.name_empty", Severity::Error, "validate", "empty assembly or element name"},
    {"object.distance_invalid", Severity::Error, "validate",
     "finite object distance not finite or <= 0 mm"},
    {"param.bound_conflict", Severity::Error, "validate", "bound Param with variable or bounds"},
    {"param.unknown_parameter", Severity::Error, "validate",
     "Param bound to an unknown parameter row"},
    {"param.unresolved", Severity::Error, "compile",
     "bound Param not evaluated yet (interim, removed with #165)"},
    {"parameters.name_duplicate", Severity::Error, "validate", "parameter row name used twice"},
    {"parameters.name_invalid", Severity::Error, "validate",
     "parameter row name not of the form [A-Za-z_][A-Za-z0-9_]*"},
    {"parameters.values_count", Severity::Error, "validate",
     "values of a row not one per configuration"},
    {"parameters.variable_expression", Severity::Error, "validate",
     "derived parameter row marked variable"},
    {"paths.crystal_mode_required", Severity::Error, "compile",
     "a crystal entered with Refract instead of Ordinary or Extraordinary"},
    {"paths.empty", Severity::Error, "validate", "no path"},
    {"paths.events_empty", Severity::Error, "validate", "explicit path without events"},
    {"paths.inner_surface_ambiguous", Severity::Error, "compile",
     "inner surface between different materials entered from outside the element"},
    {"paths.mangin_mirror_automatic", Severity::Error, "compile",
     "mirror with substrate and several surfaces on an automatic path"},
    {"paths.mode_without_crystal", Severity::Error, "compile",
     "Ordinary or Extraordinary at an event whose medium after is not a crystal"},
    {"paths.name_duplicate", Severity::Error, "validate", "path name used twice"},
    {"paths.name_empty", Severity::Error, "validate", "empty path name"},
    {"paths.order_not_allowed", Severity::Error, "validate",
     "diffraction order at a surface without phase layer"},
    {"paths.uniform_plate_automatic", Severity::Error, "compile",
     "plate of one material with more than 2 surfaces on an automatic path"},
    {"paths.unknown_surface", Severity::Error, "validate", "event at an unknown surface id"},
    {"phase.lines_per_mm_invalid", Severity::Error, "validate",
     "grating line density not finite or <= 0"},
    {"phase.radius_invalid", Severity::Error, "validate",
     "phase normalization radius not finite or <= 0 mm"},
    {"pose.no_preceding", Severity::Error, "validate",
     "relative_to_preceding without a preceding surface"},
    {"pose.no_sibling", Severity::Error, "validate",
     "relative_to_sibling without a preceding sibling"},
    {"pose.reference_unsupported", Severity::Error, "compile",
     "relative pose not evaluated yet (interim, removed with #163)"},
    {"pose.relative_first_surface", Severity::Error, "validate",
     "first surface of an element placed relatively"},
    {"rays.lost", Severity::Warning, "analysis",
     "more rays lost than the threshold of the analysis (default 50 %)"},
    {"shape.asphere_without_coefficients", Severity::Warning, "validate",
     "even asphere without coefficients"},
    {"shape.radius_invalid", Severity::Error, "validate", "radius zero or not finite"},
    {"shape.zernike_radius_invalid", Severity::Error, "validate",
     "Zernike normalization radius not finite or <= 0 mm"},
    {"shape.zernike_unsupported", Severity::Error, "compile",
     "Zernike sag terms are not supported yet"},
    {"shape.zernike_without_coefficients", Severity::Warning, "validate",
     "Zernike term without coefficients"},
    {"stop.aperture_missing", Severity::Error, "validate", "stop without aperture"},
    {"stop.clips_beam", Severity::Warning, "analysis",
     "rays of the sampling vignetted at the stop: the stop clips the beam"},
    {"stop.multiple", Severity::Error, "validate", "more than one stop element"},
    {"stop.not_on_path", Severity::Warning, "compile", "a path does not visit the stop element"},
    {"surface.efficiency_invalid", Severity::Error, "validate",
     "diffraction efficiency list empty, value outside [0, 1], order twice, or no phase layer"},
    {"surface.id_duplicate", Severity::Error, "validate", "surface id used twice"},
    {"surface.id_empty", Severity::Error, "validate", "empty surface id"},
    {"surface_aperture.half_width_invalid", Severity::Error, "validate",
     "rectangular aperture half width not finite or <= 0 mm"},
    {"surface_aperture.inner_radius_invalid", Severity::Error, "validate",
     "inner radius not finite, < 0 or >= radius"},
    {"surface_aperture.radius_invalid", Severity::Error, "validate",
     "circular aperture radius not finite or <= 0 mm"},
    {"surface_aperture.semi_axis_invalid", Severity::Error, "validate",
     "elliptical aperture semi axis not finite or <= 0 mm"},
    {"value.not_finite", Severity::Error, "validate", "number of the model not finite"},
    {"wavelengths.empty", Severity::Error, "validate", "no wavelength"},
    {"wavelengths.reference_count", Severity::Error, "validate",
     "not exactly one reference wavelength"},
    {"wavelengths.too_many", Severity::Error, "compile", "more than 65535 wavelengths"},
    {"wavelengths.value_invalid", Severity::Error, "validate", "wavelength not finite or <= 0 um"},
    {"wavelengths.weight_invalid", Severity::Error, "validate",
     "wavelength weight not finite or < 0"},
});

/// Registry entry of `code`, or nullptr if it is not registered.
[[nodiscard]] constexpr const CodeInfo* find_code(std::string_view code) noexcept {
  for (const CodeInfo& info : kCodes) {
    if (info.code == code) return &info;
  }
  return nullptr;
}

/// A registered code, checked at compile time: DiagnosticCode c = "material.unknown"; does not
/// compile for a code that is not in kCodes. The producers of diagnostics take this type, so
/// every diagnostic has a registered code and the severity of its registry entry.
class DiagnosticCode {
 public:
  // NOLINTNEXTLINE(google-explicit-constructor): string literals convert on purpose
  consteval DiagnosticCode(const char* code) : code_(code), severity_(Severity::Error) {
    const CodeInfo* info = find_code(code_);
    if (info == nullptr) throw std::invalid_argument("unregistered diagnostic code");
    code_ = info->code;
    severity_ = info->severity;
  }

  [[nodiscard]] constexpr std::string_view str() const noexcept { return code_; }
  [[nodiscard]] constexpr Severity severity() const noexcept { return severity_; }

 private:
  std::string_view code_;
  Severity severity_;
};

}  // namespace rtt::diagnostics
