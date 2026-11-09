#pragma once

/// @file export_path.hpp
/// The walk along the exported path, shared by the ZMX export and the CSV surface table
/// (ADR 0045, points 2, 3 and 6): the coaxial subset, the events and elements a sequential file
/// can hold, the media after each event and the thicknesses along z.

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
  bool mirror = false;  ///< the event is a reflection
  bool stop = false;    ///< the surface is the system stop
  /// Distance to the vertex of the next row along the global z, mm (negative after an odd number
  /// of mirrors); 0 for the last row. For a canonical system (ADR 0045, point 6) exactly the t of
  /// the relative pose of the next surface, otherwise the difference of the global z.
  double thickness = 0.0;
};

/// The exported path of one configuration.
struct ExportPath {
  model::System system;  ///< the system with the parameters of the configuration resolved
  std::vector<ExportRow> rows;
  std::optional<std::string> object_medium;  ///< reference; none for the environment
  /// Global z of the first vertex minus that of the axial object point, mm; finite object only.
  double object_distance = 0.0;
  std::vector<model::Diagnostic> warnings;  ///< export.* warnings in the order found
};

/// Compiles `system` for the configuration of `options` (coatings replaced by fresnel, so no
/// coating library is needed) and walks the chosen path. Checks the rules of ADR 0045, points 2
/// and 3, and the interactions a sequential file can hold: fresnel, ideal_mirror at a mirror,
/// ideal_lens at a thin element; coatings and ideal_anti_reflection are dropped with a warning.
/// @throws ExportError with the codes of points 2 and 3 and export.interaction_unsupported
/// @throws std::invalid_argument for an unknown path name or configuration index
/// @throws compile::CompileError if the system does not compile
[[nodiscard]] ExportPath walk_export_path(const model::System& system,
                                          const material::MaterialLibrary& materials,
                                          const ExportOptions& options);

}  // namespace rtt::io::detail
