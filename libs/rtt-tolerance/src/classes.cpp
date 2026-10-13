#include "rtt/tolerance/classes.hpp"

#include <array>
#include <stdexcept>
#include <string>

namespace rtt::tolerance {

namespace {

using K = model::ToleranceKind;
using L = model::ToleranceLevel;

/// One row of the class table: the limits of coarse, normal and precise.
struct Row {
  K kind;
  L level;
  double coarse;
  double normal;
  double precise;
};

/// The class table (ADR 0036, addendum #260, decided by the maintainer). Anchors: J. H. Burge,
/// "Tolerancing Optical Systems", OPTI 521, slide 18 (typical machining +-0.25 mm, precision
/// machining +-0.025 mm), slides 23 and 24 (null corrector: lens spacing, decenter and tilt in
/// mounting 50 um, lens thickness 25 um, radius 25 um at R = 41.6 mm and 50 um at R = 129.7 mm,
/// index 2e-4). normal is the geometric mean of coarse and precise, rounded to the nearest value
/// of the 1-2-5 series on a log scale. Tilts: an edge offset over a diameter of 50 mm. conic and
/// coarse index are own values. abbe is relative, a share of nu_d: SCHOTT TIE-29, section 5,
/// table 5-1, steps 4, 3 and 1 (+-0.8, 0.5, 0.2 %).
constexpr std::array<Row, 15> kTable{{
    {K::CurvatureRelative, L::Surface, 5e-3, 2e-3, 5e-4},
    {K::Conic, L::Surface, 0.05, 0.02, 0.005},
    {K::Thickness, L::Surface, 0.25, 0.1, 0.025},  // mm, centre thickness
    {K::Thickness, L::Element, 0.25, 0.1, 0.05},   // mm, position or air space
    {K::DecenterX, L::Element, 0.25, 0.1, 0.05},   // mm
    {K::DecenterY, L::Element, 0.25, 0.1, 0.05},   // mm
    {K::TiltXDeg, L::Element, 0.2, 0.1, 0.05},     // degree
    {K::TiltYDeg, L::Element, 0.2, 0.1, 0.05},     // degree
    {K::Index, L::Element, 1e-3, 5e-4, 2e-4},
    {K::Abbe, L::Element, 0.008, 0.005, 0.002},  // share of nu_d
    {K::DecenterX, L::Assembly, 0.25, 0.1, 0.05},  // mm
    {K::DecenterY, L::Assembly, 0.25, 0.1, 0.05},  // mm
    {K::ShiftZ, L::Assembly, 0.25, 0.1, 0.05},     // mm
    {K::TiltXDeg, L::Assembly, 0.2, 0.1, 0.05},    // degree
    {K::TiltYDeg, L::Assembly, 0.2, 0.1, 0.05},    // degree
}};

}  // namespace

double class_limit(model::ToleranceClass tolerance_class,
                   model::ToleranceKind kind,
                   model::ToleranceLevel level) {
  for (const Row& r : kTable) {
    if (r.kind != kind || r.level != level) continue;
    switch (tolerance_class) {
      case model::ToleranceClass::Coarse:
        return r.coarse;
      case model::ToleranceClass::Normal:
        return r.normal;
      case model::ToleranceClass::Precise:
        return r.precise;
      case model::ToleranceClass::None:
        break;
    }
    throw std::invalid_argument("tolerance: class_limit needs the class coarse, normal or precise");
  }
  throw std::invalid_argument("tolerance: no class generates the kind '" +
                              std::string(model::tolerance_kind_name(kind)) + "' at this level");
}

}  // namespace rtt::tolerance
