// Vertices of the polygon apertures (ADR 0050, point 2) and the centred-circle test of the
// obscurations (ADR 0050, point 3).

#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <variant>
#include <vector>

#include "rtt/model/surface.hpp"

namespace rtt::model {

std::vector<std::array<double, 2>> aperture_vertices(const Aperture& aperture) {
  if (const auto* p = std::get_if<PolygonAperture>(&aperture)) return p->vertices;
  const auto* r = std::get_if<RegularPolygonAperture>(&aperture);
  if (r == nullptr || r->sides < 3) return {};
  // phi_k = theta_0 + 2 pi k / n, first vertex on +x for theta_0 = 0; positive angles turn
  // from +x towards +y (ADR 0050, point 2; MathWorld "Regular Polygon").
  const double theta0 = r->rotation_deg * std::numbers::pi / 180.0;
  const double radius = r->circumradius.value;
  std::vector<std::array<double, 2>> out;
  out.reserve(r->sides);
  for (std::uint16_t k = 0; k < r->sides; ++k) {
    const double phi =
        theta0 + 2.0 * std::numbers::pi * static_cast<double>(k) / static_cast<double>(r->sides);
    out.push_back({radius * std::cos(phi), radius * std::sin(phi)});
  }
  return out;
}

bool is_centred_circle(const Obscuration& obscuration) noexcept {
  const auto* c = std::get_if<CircularObscuration>(&obscuration);
  return c != nullptr && c->decenter_x.value == 0.0 && c->decenter_y.value == 0.0;
}

}  // namespace rtt::model
