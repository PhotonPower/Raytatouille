#pragma once

/// @file export_path.hpp
/// The walk along the exported path, shared by the ZMX export and the CSV surface table
/// (ADR 0045, points 2, 3 and 6): the coaxial subset, the events and elements a sequential file
/// can hold, the media after each event and the thicknesses along z.

#include <array>
#include <optional>
#include <string>
#include <vector>

#include "rtt/diagnostics/codes.hpp"
#include "rtt/io/zmx.hpp"
#include "rtt/model/element.hpp"
#include "rtt/model/surface.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/validate.hpp"

namespace rtt::material {
class MaterialLibrary;
}

namespace rtt::io::detail {

/// Throws ExportError with `code` for the surface with id `surface` (empty: the whole system).
[[noreturn]] void export_fail(diagnostics::DiagnosticCode code,
                              const std::string& surface,
                              const std::string& message);

/// A warning of the export at `location` (a JSON pointer).
[[nodiscard]] model::Diagnostic export_warning(diagnostics::DiagnosticCode code,
                                               std::string location,
                                               std::string message);

/// A ZMX coordinate break with order flag 0 (UG "Coordinate Break"): decenter (dx, dy), then
/// the tilts about x, the new y and the new z, then the thickness along the new z
/// (ADR 0045, addendum #366, point 3).
struct CoordinateBreak {
  double dx = 0.0;                   ///< decenter in x, mm
  double dy = 0.0;                   ///< decenter in y, mm
  std::array<double, 3> tilt_deg{};  ///< tilts about x, y, z, degree
  double thickness = 0.0;            ///< along the new z, mm, signed as UG "Thickness"
};

/// One event of the exported path, in the global frame.
struct ExportRow {
  /// The model surface with the parameters of the configuration. Its shape is given in the
  /// global frame: a surface with the axis -e_z is mirrored, R -> -R and A_i -> -A_i for conic
  /// and even asphere, c -> -c and A_i -> -A_i for the odd asphere (ADR 0045, point 2). A
  /// dropped coating (export.coating_dropped) is replaced by fresnel.
  model::Surface surface;
  std::string element;                                 ///< name of the element of the surface
  model::ElementKind kind = model::ElementKind::Lens;  ///< kind of that element
  std::string location;  ///< JSON pointer of the surface in the system file
  double z = 0.0;        ///< global z of the vertex, mm
  /// Material reference of the medium after the event; none for the environment and after a
  /// reflection.
  std::optional<std::string> medium_after;
  /// Material reference of the medium after the event, also when it is the environment (an
  /// inner segment of a lens may be of the environment medium, e.g. an air gap).
  std::string medium_reference;
  bool mirror = false;  ///< the event is a reflection
  bool stop = false;    ///< the surface is the system stop
  /// Distance to the vertex of the next row along the global z, mm (negative after an odd number
  /// of mirrors); 0 for the last row. For a canonical system (ADR 0045, point 6) exactly the t of
  /// the relative pose of the next surface, otherwise the difference of the global z. In a
  /// folded path the thickness along z of the ZMX frame of this row (ADR 0045, addendum #366).
  double thickness = 0.0;
  /// Coordinate breaks after this row, before the next; only in a folded path.
  std::vector<CoordinateBreak> breaks;
};

/// The exported path of one configuration.
struct ExportPath {
  model::System system;  ///< the system with the parameters of the configuration resolved
  std::vector<ExportRow> rows;
  std::optional<std::string> object_medium;  ///< reference; none for the environment
  /// Global z of the first vertex minus that of the axial object point, mm; finite object only.
  /// In a folded path the thickness of surface 0 along the global z.
  double object_distance = 0.0;
  /// Coordinate breaks after surface 0, before the first row; only in a folded path.
  std::vector<CoordinateBreak> object_breaks;
  std::vector<model::Diagnostic> warnings;  ///< export.* warnings in the order found
};

/// Compiles `system` for the configuration of `options` (coatings replaced by fresnel, so no
/// coating library is needed) and walks the chosen path. Checks the rules of ADR 0045, points 2
/// and 3, and the interactions a sequential file can hold: fresnel, ideal_mirror at a mirror,
/// ideal_lens at a thin element; coatings and ideal_anti_reflection are dropped with a warning.
/// With `folded`, a path that is not coaxial is accepted: its rows get the ZMX frames, mirrored
/// shapes and coordinate breaks of ADR 0045, addendum #366 (ZMX export); without it, such a path
/// is export.not_coaxial (surface table). A coaxial path is the same either way.
/// @throws ExportError with the codes of points 2 and 3, export.interaction_unsupported and, for
/// a folded path, export.frame_ambiguous and export.surface_unsupported
/// @throws std::invalid_argument for an unknown path name or configuration index
/// @throws compile::CompileError if the system does not compile
[[nodiscard]] ExportPath walk_export_path(const model::System& system,
                                          const material::MaterialLibrary& materials,
                                          const ExportOptions& options,
                                          bool folded = false);

}  // namespace rtt::io::detail
