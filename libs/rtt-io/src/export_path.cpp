#include "export_path.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "frames.hpp"
#include "rtt/compile/compiled_system.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/isometry.hpp"
#include "rtt/math/units.hpp"
#include "rtt/model/parameters.hpp"
#include "rtt/model/path.hpp"
#include "rtt/model/pose.hpp"

namespace rtt::io {

ExportError::ExportError(std::string code, std::string surface, const std::string& message)
    : std::runtime_error((surface.empty() ? std::string() : surface + ": ") + message + " [" +
                         code + "]"),
      code_(std::move(code)),
      surface_(std::move(surface)) {}

namespace detail {

void export_fail(diagnostics::DiagnosticCode code,
                 const std::string& surface,
                 const std::string& message) {
  throw ExportError(std::string(code.str()), surface, message);
}

model::Diagnostic export_warning(diagnostics::DiagnosticCode code,
                                 std::string location,
                                 std::string message) {
  return {model::Severity::Warning, std::move(location), std::move(message),
          std::string(code.str())};
}

namespace {

using model::Element;
using model::ElementKind;
using model::EventKind;

/// Threshold of the coaxial subset (ADR 0045, point 2): mm for the vertex, dimensionless on the
/// unit vectors of the frame. Chosen like the symmetry check of paraxial::first_order, whose
/// constant lies in a private header of rtt-paraxial.
constexpr double kCoaxial = 1e-12;

/// A surface of the model tree in tree order (pre-order), as CompiledSystem::surfaces().
struct TreeSurface {
  const model::Surface* surface = nullptr;
  const Element* element = nullptr;
  /// t of the pose that places this surface relative to the surface before it, for a canonical
  /// system: the element pose for the first surface of an element, else the surface pose.
  double t = 0.0;
};

void collect(const model::Assembly& a, std::vector<TreeSurface>& out) {
  for (const model::Node& n : a.children) {
    if (const auto* e = std::get_if<Element>(&n.value)) {
      for (std::size_t i = 0; i < e->surfaces.size(); ++i) {
        const model::Pose& p = i == 0 ? e->pose : e->surfaces[i].pose;
        out.push_back({&e->surfaces[i], e, p.position[2].value});
      }
    } else {
      collect(std::get<model::Assembly>(n.value), out);
    }
  }
}

/// True for a pose (0, 0, t) without rotation, relative to `reference`.
bool plain_shift(const model::Pose& p, model::PoseReference reference) {
  return p.reference == reference && p.position[0].value == 0.0 && p.position[1].value == 0.0 &&
         p.rotation_deg[0].value == 0.0 && p.rotation_deg[1].value == 0.0 &&
         p.rotation_deg[2].value == 0.0;
}

/// The canonical form of ADR 0045, point 6 (on the resolved system).
bool canonical(const model::System& s) {
  if (!s.root.pose.is_identity() || s.object.reference) return false;
  if (s.paths.size() != 1 || !s.paths.front().automatic) return false;
  bool first = true;
  for (const model::Node& n : s.root.children) {
    const auto* e = std::get_if<Element>(&n.value);
    if (e == nullptr) return false;
    const auto& m = e->segment_materials;
    if (!m.empty() && std::all_of(m.begin(), m.end(), [&](const auto& x) { return x == m[0]; })) {
      return false;
    }
    if (first) {
      if (!plain_shift(e->pose, model::PoseReference::Absolute) ||
          e->pose.position[2].value != 0.0) {
        return false;
      }
    } else if (!plain_shift(e->pose, model::PoseReference::RelativeToPreceding)) {
      return false;
    }
    first = false;
    for (std::size_t i = 0; i < e->surfaces.size(); ++i) {
      const model::Pose& p = e->surfaces[i].pose;
      if (i == 0 ? !p.is_identity() : !plain_shift(p, model::PoseReference::RelativeToPreceding)) {
        return false;
      }
    }
  }
  return true;
}

/// The kind of the event the automatic path has at an element of this kind (compile()).
EventKind automatic_kind(ElementKind k) {
  switch (k) {
    case ElementKind::Lens:
    case ElementKind::Plate:
      return EventKind::Refract;
    case ElementKind::Mirror:
      return EventKind::Reflect;
    case ElementKind::ThinElement:
    case ElementKind::Stop:
    case ElementKind::Detector:
      return EventKind::Transmit;
  }
  return EventKind::Refract;
}

/// The shape seen from the global frame of a surface with the axis -e_z: z_global = -f(r), the
/// root of the conic holds only c^2, so -f has the curvature -c, the same conic and -A_i
/// (ADR 0045, point 2). + 0.0 keeps a zero coefficient +0.
void mirror_shape(model::ShapeStack& shape) {
  const auto negate = [](model::Param& p) { p.value = -p.value + 0.0; };
  if (auto* c = std::get_if<model::Conic>(&shape.base)) {
    negate(c->radius);
  } else if (auto* a = std::get_if<model::EvenAsphere>(&shape.base)) {
    negate(a->radius);
    for (model::Param& x : a->coefficients) negate(x);
  } else if (auto* o = std::get_if<model::OddAsphere>(&shape.base)) {
    negate(o->curvature);
    for (model::Param& x : o->coefficients) negate(x);
  }
  // other shapes have no export mapping (export.surface_unsupported in the writers)
}

/// Checks the interaction of one surface on the path and drops a coating (ADR 0045, point 4).
void check_interaction(ExportRow& row, const Element& e, std::vector<model::Diagnostic>& w) {
  const std::string& id = row.surface.id.str();
  if (!row.surface.phases.empty()) {
    export_fail("export.interaction_unsupported", id, "a phase layer or grating has no mapping");
  }
  if (e.crystal) {
    export_fail("export.interaction_unsupported", id, "a crystal has no mapping");
  }
  model::Interaction& i = row.surface.interaction;
  if (std::holds_alternative<model::CoatingRef>(i) ||
      std::holds_alternative<model::IdealAntiReflection>(i)) {
    w.push_back(export_warning("export.coating_dropped", row.location,
                               "coating of " + id + " not exported"));
    i = model::Fresnel{};
  } else if (std::holds_alternative<model::IdealMirror>(i)) {
    if (e.kind != ElementKind::Mirror) {
      export_fail("export.interaction_unsupported", id,
                  "an ideal mirror outside a mirror element has no mapping");
    }
  } else if (std::holds_alternative<model::IdealLens>(i)) {
    if (e.kind != ElementKind::ThinElement) {
      export_fail("export.interaction_unsupported", id,
                  "an ideal lens outside a thin element has no mapping");
    }
  } else if (!std::holds_alternative<model::Fresnel>(i)) {
    export_fail("export.interaction_unsupported", id, "the interaction has no mapping");
  }
}

/// Threshold of the frames of a folded path (ADR 0045, addendum #366): as kCoaxial, mm for
/// positions, dimensionless on the unit vectors.
constexpr double kFrame = 1e-12;

/// The sign of the propagation along z of the ZMX frame after `mirrors` mirrors: +1 after an
/// even number, -1 after an odd one (UG "Thickness"; KA-01588).
double parity(int mirrors) {
  return mirrors % 2 == 0 ? 1.0 : -1.0;
}

/// a^T b, element by element: GCC 13 with -O3 reports -Wnull-dereference inside the vectorized
/// Eigen 5 product of a transpose (vcpkg), so the 3x3 products are written out.
math::Mat3 tmul(const math::Mat3& a, const math::Mat3& b) {
  math::Mat3 out;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      out(i, j) = a(0, i) * b(0, j) + a(1, i) * b(1, j) + a(2, i) * b(2, j);
    }
  }
  return out;
}

/// a^T v, element by element (see tmul).
math::Vec3 tmul(const math::Mat3& a, const math::Vec3& v) {
  math::Vec3 out;
  for (int i = 0; i < 3; ++i) out(i) = a(0, i) * v(0) + a(1, i) * v(1) + a(2, i) * v(2);
  return out;
}

bool is_identity(const math::Mat3& r) {
  return (r - math::Mat3::Identity()).cwiseAbs().maxCoeff() <= kFrame;
}

/// The angle of the rotation r, rad.
double turn_angle(const math::Mat3& r) {
  return std::acos(std::clamp((r.trace() - 1.0) / 2.0, -1.0, 1.0));
}

/// A frame of the ZMX file in global coordinates: rotation (axes in the columns) and origin, mm.
struct ZmxFrame {
  math::Mat3 r = math::Mat3::Identity();
  math::Vec3 o = math::Vec3::Zero();
};

/// The step from the frame `a` of one row to the frame `b` of the next (ADR 0045, addendum
/// #366, point 3): the thickness at the row of `a` and the coordinate breaks after it. `s` is
/// the sign of the propagation on the segment, +1 or -1.
void frame_step(const ZmxFrame& a,
                const ZmxFrame& b,
                double s,
                double& thickness,
                std::vector<CoordinateBreak>& breaks) {
  const math::Mat3 turn = tmul(a.r, b.r);
  const math::Vec3 d = b.o - a.o;
  const math::Vec3 p = tmul(a.r, d);
  if (is_identity(turn)) {
    // cases 1 and 2: a thickness, a decenter beside it
    thickness = p.z();
    if (std::abs(p.x()) > kFrame || std::abs(p.y()) > kFrame) {
      breaks.push_back({p.x(), p.y(), {0.0, 0.0, 0.0}, 0.0});
    }
    return;
  }
  // case 3, the axis frame: break 1 turns z onto s d/|d| and carries s |d|, break 2 turns into b
  thickness = 0.0;
  math::Mat3 axis = a.r;
  const double length = d.norm();
  if (length > kFrame) {
    const math::Vec3 u = s * p / length;
    if (std::abs(u.x()) <= kFrame && std::abs(u.y()) <= kFrame) {
      thickness = s * length;
    } else {
      // Rx(ta) Ry(tb) e_z = (sin tb, -sin ta cos tb, cos ta cos tb) (derived in #366)
      const double tb = std::asin(std::clamp(u.x(), -1.0, 1.0));
      const double ta = std::atan2(-u.y(), u.z());
      breaks.push_back(
          {0.0, 0.0, {math::rad_to_deg(ta) + 0.0, math::rad_to_deg(tb) + 0.0, 0.0}, s * length});
      axis = a.r * math::rotation_xyz_intrinsic(math::Vec3(ta, tb, 0.0));
    }
  }
  const math::Mat3 rest = tmul(axis, b.r);
  if (!is_identity(rest)) breaks.push_back({0.0, 0.0, xyz_angles_deg(rest), 0.0});
}

/// The ZMX frames, mirrored shapes and coordinate breaks of a folded path (ADR 0045, addendum
/// #366, points 2 and 3). `rot` and `pos` are the global frames of the rows; the entrance frame
/// is the global one, in the axial object point (finite object) or in the first vertex.
void fold_path(ExportPath& out,
               const std::vector<math::Mat3>& rot,
               const std::vector<math::Vec3>& pos,
               double object_z) {
  const bool infinite = out.system.object.at_infinity;
  const std::size_t n = out.rows.size();
  const ZmxFrame entrance{math::Mat3::Identity(),
                          infinite ? pos.front() : math::Vec3(0.0, 0.0, object_z)};
  const math::Mat3 rz = math::Vec3(-1.0, -1.0, 1.0).asDiagonal();
  const math::Mat3 rx = math::Vec3(1.0, -1.0, -1.0).asDiagonal();
  const math::Mat3 ry = math::Vec3(-1.0, 1.0, -1.0).asDiagonal();
  std::vector<ZmxFrame> frames(n);
  int mirrors = 0;
  for (std::size_t k = 0; k < n; ++k) {
    ExportRow& row = out.rows[k];
    const std::string& id = row.surface.id.str();
    const int after = mirrors + (row.mirror ? 1 : 0);
    // true if the light on d runs along s z of the global frame of the row (point 2); (G^T d)_z
    // is the z axis of G dotted with d
    const math::Vec3 axis = rot.at(k).col(2);
    const auto along = [&](const math::Vec3& d, double s) {
      const double z = axis.dot(d);
      if (std::abs(z) <= kFrame * d.norm()) {
        export_fail("export.frame_ambiguous", id,
                    "the light runs in the plane of the surface, so its side is not defined");
      }
      return (z > 0.0) == (s > 0.0);
    };
    const math::Vec3 d_in =
        k > 0 ? math::Vec3(pos[k] - pos[k - 1])
              : (infinite ? math::Vec3(0.0, 0.0, 1.0) : math::Vec3(pos[k] - entrance.o));
    const math::Vec3 d_out = k + 1 < n ? math::Vec3(pos[k + 1] - pos[k]) : math::Vec3::Zero();
    bool keep = true;
    if (d_out.norm() > kFrame) {
      keep = along(d_out, parity(after));
      if (!row.mirror && d_in.norm() > kFrame && along(d_in, parity(mirrors)) != keep) {
        export_fail("export.frame_ambiguous", id,
                    "the light arrives and leaves on different sides of the surface");
      }
    } else if (d_in.norm() > kFrame) {
      keep = along(d_in, parity(mirrors));
    } else {
      export_fail("export.frame_ambiguous", id, "no direction of the light at the surface");
    }
    // the turn by 180 deg with the smaller rotation from the frame before; E or Rx on a tie
    const ZmxFrame& before = k > 0 ? frames[k - 1] : entrance;
    const math::Mat3 first = rot[k] * (keep ? math::Mat3(math::Mat3::Identity()) : rx);
    const math::Mat3 second = rot[k] * (keep ? rz : ry);
    const bool take_second =
        turn_angle(tmul(before.r, second)) < turn_angle(tmul(before.r, first)) - 1e-9;
    frames[k] = {take_second ? second : first, pos[k]};
    if (!keep) mirror_shape(row.surface.shape);
    if (k == 0) {
      // surface 0: its thickness along the global z, then one break with decenter and turn
      const math::Mat3 turn = frames[0].r;
      const math::Vec3 p = frames[0].o - entrance.o;
      out.object_distance = p.z();
      if (!is_identity(turn) || std::abs(p.x()) > kFrame || std::abs(p.y()) > kFrame) {
        out.object_breaks.push_back({p.x(), p.y(), xyz_angles_deg(turn), 0.0});
      }
    } else {
      ExportRow& last = out.rows[k - 1];
      frame_step(frames[k - 1], frames[k], parity(mirrors), last.thickness, last.breaks);
    }
    mirrors = after;
  }
}

}  // namespace

ExportPath walk_export_path(const model::System& system,
                            const material::MaterialLibrary& materials,
                            const ExportOptions& options,
                            bool folded) {
  ExportPath out;

  // Configuration (ADR 0045, point 3; F6).
  const std::size_t count = model::configuration_count(system);
  if (!options.configuration && count > 1) {
    export_fail("export.multi_config", "",
                "the system has " + std::to_string(count) + " configurations; choose one");
  }
  const std::size_t configuration = options.configuration.value_or(0);
  if (configuration >= count) {
    throw std::invalid_argument("configuration " + std::to_string(configuration) +
                                " is beyond the " + std::to_string(count) + " of the system");
  }

  // Compile without coatings: they are dropped anyway, and no coating library is needed.
  model::System bare = system;
  const auto strip = [](const auto& self, model::Assembly& a) -> void {
    for (model::Node& n : a.children) {
      if (auto* e = std::get_if<Element>(&n.value)) {
        for (model::Surface& s : e->surfaces) {
          if (std::holds_alternative<model::CoatingRef>(s.interaction)) {
            s.interaction = model::Fresnel{};
          }
        }
      } else {
        self(self, std::get<model::Assembly>(n.value));
      }
    }
  };
  strip(strip, bare.root);
  const compile::CompiledSystem cs = compile::compile(bare, materials, configuration);
  out.system = model::resolve_parameters(system, configuration);

  // Path (F6).
  compile::PathId path_id{0};
  if (options.path) {
    const auto found = cs.find_path(*options.path);
    if (!found) throw std::invalid_argument("the system has no path '" + *options.path + "'");
    path_id = *found;
  }
  if (cs.paths().size() > 1) {
    out.warnings.push_back(
        export_warning("export.paths_dropped", "/paths",
                       "only the path '" + cs.path(path_id).name + "' is exported, " +
                           std::to_string(cs.paths().size() - 1) + " other paths are dropped"));
  }
  const compile::CompiledPath& path = cs.path(path_id);

  std::vector<TreeSurface> tree;
  collect(out.system.root, tree);
  if (tree.size() != cs.surfaces().size()) {
    throw std::logic_error("export: the model tree and the compiled surfaces differ");
  }

  const std::uint32_t env = cs.environment_medium();
  const std::string& env_name = cs.media()[env].reference;

  // Object space before the path (ADR 0045, point 3).
  if (out.system.object.medium && *out.system.object.medium != env_name) {
    const compile::CompiledEvent* first = path.events.empty() ? nullptr : &path.events.front();
    const bool lens =
        first != nullptr && first->kind == EventKind::Refract && first->medium_after != env;
    if (!lens) {
      export_fail("export.object_medium_unsupported",
                  first != nullptr ? cs.surfaces()[first->surface].id.str() : std::string(),
                  "glass in object space needs a path that starts in a lens");
    }
    out.object_medium = out.system.object.medium;
  }

  const std::optional<std::uint32_t> stop_surface = cs.stop_surface();
  std::set<std::uint32_t> visited;
  std::vector<math::Mat3> rot;  // global frames of the rows
  std::vector<math::Vec3> pos;
  std::vector<bool> axis_down;  // the axis of the row is -e_z (coaxial path)
  bool coaxial_path = true;
  for (std::size_t k = 0; k < path.events.size(); ++k) {
    const compile::CompiledEvent& e = path.events[k];
    const compile::CompiledSurface& cs_surface = cs.surfaces()[e.surface];
    const TreeSurface& t = tree[e.surface];
    const Element& element = *t.element;
    const std::string& id = cs_surface.id.str();

    // refract, the default kind of an event, keeps the medium at a stop or detector (elements
    // without material) as transmit does (ADR 0045, addendum PR B)
    const bool same_as_transmit =
        e.kind == EventKind::Refract &&
        (element.kind == ElementKind::Stop || element.kind == ElementKind::Detector);
    if ((e.kind != automatic_kind(element.kind) && !same_as_transmit) || e.order != 0) {
      export_fail("export.event_unsupported", id,
                  "a ZMX surface holds only the event of the automatic path of its element");
    }
    if (!visited.insert(e.surface).second) {
      export_fail("export.event_unsupported", id, "the path visits the surface twice");
    }
    if (element.kind == ElementKind::Detector && k + 1 != path.events.size()) {
      export_fail("export.event_unsupported", id, "a detector before the last event");
    }

    // Coaxial (ADR 0045, point 2).
    const math::Vec3& p = cs_surface.to_global.translation();
    const math::Mat3& r = cs_surface.to_global.rotation();
    bool coaxial = std::abs(p.x()) <= kCoaxial && std::abs(p.y()) <= kCoaxial;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        if (i != j && std::abs(r(i, j)) > kCoaxial) coaxial = false;
      }
    }
    if (!coaxial && !folded) {
      export_fail("export.not_coaxial", id,
                  "only coaxial systems are exported: the vertex must lie on the z axis and the "
                  "surface axes along +-x, +-y, +-z");
    }
    coaxial_path = coaxial_path && coaxial;

    const bool glass_before = e.medium_before != env;
    const bool glass_after = e.medium_after != env;
    if (e.kind == EventKind::Reflect && glass_before) {
      export_fail("export.mirror_in_glass", id, "a mirror inside glass (Mangin)");
    }
    if ((element.kind == ElementKind::Stop || element.kind == ElementKind::ThinElement ||
         element.kind == ElementKind::Detector) &&
        (glass_before || glass_after)) {
      export_fail("export.element_in_glass", id,
                  "a stop, thin element or detector needs the environment on both sides");
    }

    ExportRow row;
    row.surface = *t.surface;
    row.element = element.name;
    row.kind = element.kind;
    row.location = cs_surface.location;
    row.z = p.z();
    row.mirror = e.kind == EventKind::Reflect;
    row.stop = stop_surface == std::optional<std::uint32_t>(e.surface);
    if (glass_after && !row.mirror) row.medium_after = cs.media()[e.medium_after].reference;
    row.medium_reference = cs.media()[e.medium_after].reference;
    rot.push_back(r);
    pos.push_back(p);
    axis_down.push_back(r(2, 2) < 0.0);
    check_interaction(row, element, out.warnings);
    if (element.kind == ElementKind::Mirror &&
        (element.material || !element.segment_materials.empty())) {
      out.warnings.push_back(export_warning("export.mirror_material_dropped", row.location,
                                            "substrate of " + element.name + " not exported"));
    }
    out.rows.push_back(std::move(row));
  }

  if (out.rows.empty() || out.rows.back().kind != ElementKind::Detector) {
    export_fail("export.image_missing",
                out.rows.empty() ? std::string() : out.rows.back().surface.id.str(),
                "the last event of the path must be a detector (the image surface)");
  }
  const auto stop =
      std::find_if(out.rows.begin(), out.rows.end(), [](const ExportRow& x) { return x.stop; });
  if (stop == out.rows.end()) {
    export_fail("export.stop_missing", "", "the exported path has no stop");
  }
  const std::optional<model::Aperture>& stop_aperture = stop->surface.aperture;
  if (stop_aperture.has_value() &&
      !std::holds_alternative<model::CircularAperture>(stop_aperture.value())) {
    export_fail("export.aperture_unsupported", stop->surface.id.str(),
                "a rectangular or elliptical stop has no mapping");
  }

  for (std::size_t i = 0; i < cs.surfaces().size(); ++i) {
    if (!visited.contains(static_cast<std::uint32_t>(i))) {
      out.warnings.push_back(
          export_warning("export.surface_off_path", cs.surfaces()[i].location,
                         "surface " + cs.surfaces()[i].id.str() + " is not on the exported path"));
    }
  }

  // An ideal lens after an odd number of mirrors (ADR 0045, addendum #366, point 5).
  int mirrors = 0;
  for (const ExportRow& row : out.rows) {
    if (std::holds_alternative<model::IdealLens>(row.surface.interaction) && mirrors % 2 != 0) {
      export_fail("export.surface_unsupported", row.surface.id.str(),
                  "an ideal lens after an odd number of mirrors: no source gives the sign of the "
                  "focal length of a PARAXIAL surface with the light along -z");
    }
    if (row.mirror) ++mirrors;
  }

  if (!coaxial_path) {
    // a folded path (ADR 0045, addendum #366): frames, mirrored shapes, coordinate breaks
    fold_path(out, rot, pos, cs.object_z());
    return out;
  }
  for (std::size_t k = 0; k < out.rows.size(); ++k) {
    if (axis_down[k]) mirror_shape(out.rows[k].surface.shape);
  }

  // Thicknesses (ADR 0045, point 4): the t of the poses if the system is canonical (then the
  // path is automatic and visits every surface in tree order), else differences of z.
  const bool exact = canonical(out.system);
  for (std::size_t k = 0; k + 1 < out.rows.size(); ++k) {
    out.rows[k].thickness = exact ? tree[k + 1].t : out.rows[k + 1].z - out.rows[k].z;
  }
  if (!out.system.object.at_infinity) {
    out.object_distance =
        exact ? out.system.object.distance.value : out.rows.front().z - cs.object_z();
  }
  return out;
}

}  // namespace detail
}  // namespace rtt::io
