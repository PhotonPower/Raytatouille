#pragma once

#include <algorithm>
#include <string>
#include <vector>

#include "rtt/model/model.hpp"

namespace rtt::model::test {

/// Minimal consistent system: stop, plano-convex singlet, detector, automatic path.
inline System make_singlet() {
  System s;
  s.name = "singlet";
  s.wavelengths = {{0.5876, 1.0, true}};
  s.aperture = {SystemApertureType::EntrancePupilDiameter, Param(20.0)};
  s.fields = {FieldType::AngleDeg, {{0.0, 0.0, 1.0}, {0.0, 5.0, 1.0}}};
  s.root.name = "system";

  Surface stop_surface;
  stop_surface.id = SurfaceId("STO");
  stop_surface.aperture = CircularAperture{10.0, 0.0};
  s.root.children.push_back({Element{"stop", ElementKind::Stop, {}, std::nullopt, {stop_surface}}});

  Surface s1;
  s1.id = SurfaceId("L1.S1");
  s1.shape.base = Conic{Param(51.68), Param(0.0)};
  s1.aperture = CircularAperture{12.7, 0.0};
  Surface s2;
  s2.id = SurfaceId("L1.S2");
  s2.pose = Pose::along_z(4.0);
  s2.aperture = CircularAperture{12.7, 0.0};
  s.root.children.push_back(
      {Element{"L1", ElementKind::Lens, Pose::along_z(5.0), "SCHOTT:N-BK7", {s1, s2}}});

  Surface image;
  image.id = SurfaceId("IMG");
  s.root.children.push_back(
      {Element{"image", ElementKind::Detector, Pose::along_z(106.363), std::nullopt, {image}}});

  s.paths = {{"main", true, {}}};
  return s;
}

inline Element& element(System& s, std::size_t i) {
  return std::get<Element>(s.root.children[i].value);
}

/// True if any error diagnostic has exactly this location.
inline bool has_error_at(const std::vector<Diagnostic>& d, const std::string& location) {
  return std::any_of(d.begin(), d.end(), [&](const Diagnostic& x) {
    return x.severity == Severity::Error && x.location == location;
  });
}

}  // namespace rtt::model::test
