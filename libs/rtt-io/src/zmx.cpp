#include "rtt/io/zmx.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "rtt/compile/compiled_system.hpp"
#include "rtt/diagnostics/codes.hpp"
#include "rtt/material/agf.hpp"
#include "rtt/material/air.hpp"
#include "rtt/material/material.hpp"
#include "rtt/math/isometry.hpp"
#include "rtt/math/units.hpp"
#include "rtt/model/element.hpp"
#include "rtt/model/surface.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/trace/sources.hpp"

// rtt-compile declares line_vacuum_wavelength_um in perturbation.hpp.
#include "numbers.hpp"
#include "rtt/compile/perturbation.hpp"
#include "zmx_names.hpp"

namespace rtt::io {

ZmxError::ZmxError(std::string code, std::size_t line, const std::string& message)
    : std::runtime_error("line " + std::to_string(line) + ": " + message + " [" + code + "]"),
      code_(std::move(code)),
      line_(line) {}

namespace {

using detail::mil_number;
using detail::upper;
using diagnostics::DiagnosticCode;
using model::Diagnostic;

[[noreturn]] void fail(DiagnosticCode code, std::size_t line, const std::string& message) {
  throw ZmxError(std::string(code.str()), line, message);
}

/// One line of the file: 1-based number, indentation and the whitespace-separated words.
struct Line {
  std::size_t number = 0;
  bool indented = false;
  std::vector<std::string> words;
  std::string rest;  ///< text after the keyword, trimmed (for NAME)
};

std::vector<Line> split_lines(const std::string& text) {
  std::vector<Line> out;
  std::size_t pos = 0;
  std::size_t number = 0;
  while (pos < text.size()) {
    std::size_t end = text.find('\n', pos);
    if (end == std::string::npos) end = text.size();
    std::string_view raw(text.data() + pos, end - pos);
    if (!raw.empty() && raw.back() == '\r') raw.remove_suffix(1);
    ++number;
    Line line;
    line.number = number;
    line.indented = !raw.empty() && (raw.front() == ' ' || raw.front() == '\t');
    std::size_t i = 0;
    while (i < raw.size()) {
      while (i < raw.size() && std::isspace(static_cast<unsigned char>(raw[i])) != 0) ++i;
      std::size_t j = i;
      while (j < raw.size() && std::isspace(static_cast<unsigned char>(raw[j])) == 0) ++j;
      if (j > i) {
        if (line.words.empty()) {
          std::size_t r = j;
          while (r < raw.size() && std::isspace(static_cast<unsigned char>(raw[r])) != 0) ++r;
          std::size_t e = raw.size();
          while (e > r && std::isspace(static_cast<unsigned char>(raw[e - 1])) != 0) --e;
          line.rest = std::string(raw.substr(r, e - r));
        }
        line.words.emplace_back(raw.substr(i, j - i));
      }
      i = j;
    }
    if (!line.words.empty()) out.push_back(std::move(line));
    pos = end + 1;
  }
  return out;
}

/// Number of word `i` of `line`, read independently of the locale (std::from_chars).
double number(const Line& line, std::size_t i) {
  if (i >= line.words.size()) {
    fail("zmx.parse", line.number,
         line.words[0] + ": expected at least " + std::to_string(i) + " values");
  }
  const std::string& w = line.words[i];
  const std::optional<double> v = detail::parse_number(w);
  if (!v) fail("zmx.parse", line.number, line.words[0] + ": '" + w + "' is not a finite number");
  return *v;
}

int integer(const Line& line, std::size_t i) {
  const double v = number(line, i);
  if (!(std::abs(v) <= 1e9) || v != static_cast<double>(static_cast<int>(v))) {
    fail("zmx.parse", line.number, line.words[0] + ": '" + line.words[i] + "' is not an integer");
  }
  return static_cast<int>(v);
}

/// Keywords that only carry display or note data, skipped without a warning (ADR 0044, point 2).
/// PR C adds the pickup and variable notes PPAR, PZUP, VPAR (the resolved value stands next to
/// them), TCED (expansion coefficient, thermal analyses only; reading) and SSID (surface id of
/// newer files).
bool silent(std::string_view k) {
  static const std::set<std::string_view> kSilent{
      "VERS", "NAME", "NOTE", "COMM", "HIDE", "SLAB", "POPS", "MIRR", "MAZH", "EDGE",
      "MEMA", "PUPD", "PICB", "ROPD", "COAT", "COFN", "POLS", "NSCD", "GFAC", "GLRS",
      "GSTD", "RAIM", "PUSH", "SDMA", "PPAR", "PZUP", "VPAR", "TCED", "SSID"};
  return kSilent.contains(k);
}

/// Angles (a, b, c) in degrees with R = Rx(a) Ry(b) Rz(c) (math::rotation_xyz_intrinsic).
/// Derived in #317 (ADR 0044, addendum PR C) by multiplying out: row 0 of R is
/// (cos b cos c, -cos b sin c, sin b) and column 2 is (sin b, -sin a cos b, cos a cos b), so
/// b = atan2(R02, hypot(R00, R01)) and a = atan2(-R12, R22); row 1 of Rx(a)^T R = Ry(b) Rz(c)
/// is (sin c, cos c, 0), so c = atan2(cos a R10 + sin a R20, cos a R11 + sin a R21). That form
/// needs no case at b = +-90 deg, where only a +- c is defined: c takes up the noise in a.
std::array<double, 3> xyz_angles_deg(const math::Mat3& r) {
  const double b = std::atan2(r(0, 2), std::hypot(r(0, 0), r(0, 1)));
  const double a = std::atan2(-r(1, 2), r(2, 2));
  const double ca = std::cos(a);
  const double sa = std::sin(a);
  const double c = std::atan2(ca * r(1, 0) + sa * r(2, 0), ca * r(1, 1) + sa * r(2, 1));
  // + 0.0 turns a -0.0 into 0.0 for the file
  return {math::rad_to_deg(a) + 0.0, math::rad_to_deg(b) + 0.0, math::rad_to_deg(c) + 0.0};
}

/// One step between two placed surfaces: a thickness along the current z, or a coordinate break
/// (UG "Coordinate Break": decenters dx, dy in mm, tilts a, b, c in degrees about x, y, z).
struct FrameStep {
  bool coordinate_break = false;
  double t = 0.0;
  double dx = 0.0;
  double dy = 0.0;
  std::array<double, 3> tilt{0.0, 0.0, 0.0};
  bool tilts_first = false;  ///< order flag != 0: tilt about z, y, x, then decenter
};

/// The way from the frame of the last placed surface to the next surface (ADR 0044, addendum
/// PR C, points 3 and 4).
struct FrameChain {
  std::optional<math::Mat3> untilt;  ///< R^T of the tilted plane placed last (TILTSURF)
  std::vector<FrameStep> steps;

  void thickness(double t) {
    if (!steps.empty() && !steps.back().coordinate_break) {
      steps.back().t += t;  // summed from front to back, as the dropped dummies of PR B
    } else {
      steps.push_back(FrameStep{false, t, 0.0, 0.0, {0.0, 0.0, 0.0}, false});
    }
  }
};

/// Rigid transform of a coordinate break (UG "Coordinate Break").
math::Isometry3 coordinate_break(const FrameStep& s) {
  using math::deg_to_rad;
  const math::Isometry3 shift(math::Mat3::Identity(), math::Vec3(s.dx, s.dy, 0.0));
  if (!s.tilts_first) {
    // decenter, then tilt about x, the new y, the new z: Rx Ry Rz
    const math::Vec3 angles(deg_to_rad(s.tilt[0]), deg_to_rad(s.tilt[1]), deg_to_rad(s.tilt[2]));
    return shift * math::Isometry3(math::rotation_xyz_intrinsic(angles), math::Vec3::Zero());
  }
  // tilt about z, the new y, the new x (Rz Ry Rx), then decenter along the turned axes
  const auto axis = [](int i, double deg) {
    math::Vec3 v = math::Vec3::Zero();
    v[i] = deg_to_rad(deg);
    return math::rotation_xyz_intrinsic(v);
  };
  const math::Mat3 r = axis(2, s.tilt[2]) * axis(1, s.tilt[1]) * axis(0, s.tilt[0]);
  return math::Isometry3(r, math::Vec3::Zero()) * shift;
}

/// The pose of the next surface relative to the frame of the last placed one (`own_tilt`: the
/// rotation of the next surface itself, for an inner lens surface that is a tilted plane).
/// Without a break and a tilt: position (0, 0, sum of thicknesses), as in PR B. One break with
/// flag 0, or with flag != 0 and at most one tilt: rotate_first with the tilts as in the file,
/// p_ref = Q + R (position + p - Q) (ADR 0028): flag 0 Q = (dx, dy, t1), flag != 0
/// Q = (0, 0, t1), position (dx, dy, t1 + t2). Only the untilt of a tilted plane: rotate_first
/// about the vertex with position (0, 0, t). Otherwise the multiplied-out transform, decomposed.
model::Pose chain_pose(const FrameChain& chain, const std::optional<math::Mat3>& own_tilt) {
  model::Pose pose;
  std::size_t breaks = 0;
  const FrameStep* cb = nullptr;
  for (const FrameStep& s : chain.steps) {
    if (s.coordinate_break) {
      ++breaks;
      cb = &s;
    }
  }
  const auto sum_t = [&](bool before_break) {
    double t = 0.0;
    bool seen = false;
    for (const FrameStep& s : chain.steps) {
      if (s.coordinate_break) {
        seen = true;
      } else if (seen != before_break) {
        t += s.t;
      }
    }
    return t;
  };
  if (!own_tilt && !chain.untilt && breaks == 0) {
    pose.position[2] = sum_t(true);
    return pose;
  }
  if (!own_tilt && !chain.untilt && breaks == 1 && cb != nullptr) {
    const int tilts = static_cast<int>(cb->tilt[0] != 0.0) + static_cast<int>(cb->tilt[1] != 0.0) +
                      static_cast<int>(cb->tilt[2] != 0.0);
    if (!cb->tilts_first || tilts <= 1) {
      const double t1 = sum_t(true);
      const double t2 = sum_t(false);
      pose.order = model::PoseOrder::RotateFirst;
      pose.pivot = cb->tilts_first ? std::array<double, 3>{0.0, 0.0, t1}
                                   : std::array<double, 3>{cb->dx, cb->dy, t1};
      pose.position = {cb->dx, cb->dy, t1 + t2};
      pose.rotation_deg = {cb->tilt[0], cb->tilt[1], cb->tilt[2]};
      return pose;
    }
  }
  if (!own_tilt && chain.untilt && breaks == 0) {
    const std::array<double, 3> a = xyz_angles_deg(*chain.untilt);
    pose.order = model::PoseOrder::RotateFirst;
    pose.position[2] = sum_t(true);
    pose.rotation_deg = {a[0], a[1], a[2]};
    return pose;
  }
  math::Isometry3 m(chain.untilt.value_or(math::Mat3::Identity()), math::Vec3::Zero());
  for (const FrameStep& s : chain.steps) {
    m = m * (s.coordinate_break
                 ? coordinate_break(s)
                 : math::Isometry3(math::Mat3::Identity(), math::Vec3(0.0, 0.0, s.t)));
  }
  if (own_tilt) m = m * math::Isometry3(*own_tilt, math::Vec3::Zero());
  const std::array<double, 3> a = xyz_angles_deg(m.rotation());
  pose.position = {m.translation().x(), m.translation().y(), m.translation().z()};
  pose.rotation_deg = {a[0], a[1], a[2]};
  return pose;
}

/// Tilt of a plane z = x tan_x + y tan_y (UG "Tilted") as R = Rx(alpha) Ry(beta) in degrees:
/// R e_z = (sin b, -sin a cos b, cos a cos b) is the normal (-tan_x, -tan_y, 1) / norm for
/// alpha = atan(tan_y), beta = -atan(tan_x cos alpha) (derived in #317; CEO computes the same
/// angles).
std::array<double, 3> plane_tilt_deg(double tan_x, double tan_y) {
  const double alpha = std::atan(tan_y);
  const double beta = -std::atan(tan_x * std::cos(alpha));
  return {math::rad_to_deg(alpha) + 0.0, math::rad_to_deg(beta) + 0.0, 0.0};
}

math::Mat3 rotation_of(const std::array<double, 3>& deg) {
  return math::rotation_xyz_intrinsic(
      math::Vec3(math::deg_to_rad(deg[0]), math::deg_to_rad(deg[1]), math::deg_to_rad(deg[2])));
}

struct Located {
  double value = 0.0;
  std::size_t line = 0;
};

/// The data of one SURF block.
struct ZSurface {
  int number = 0;
  std::size_t line = 0;
  std::string type = "STANDARD";
  std::size_t type_line = 0;
  Located curv;
  Located coni;
  std::map<int, Located> parm;
  bool infinite = false;
  Located disz;
  std::string glass;  ///< empty: air
  std::size_t glass_line = 0;
  std::vector<std::string> glass_words;
  Located diam;
  int diam_flag = 0;
  std::optional<std::pair<double, double>> clap;  ///< rmin, rmax
  std::size_t clap_line = 0;
  std::optional<std::pair<double, double>> flap;  ///< first field, rmax
  std::size_t flap_line = 0;
  std::optional<std::pair<double, double>> sqap;  ///< half widths x, y (lens units)
  std::size_t sqap_line = 0;
  std::optional<std::pair<double, double>> elap;  ///< half widths x, y (lens units)
  std::size_t elap_line = 0;
  bool stop = false;
  std::size_t stop_line = 0;
};

bool is_mirror(const ZSurface& s) {
  return s.glass == "MIRROR";
}

/// Tracks keywords that are neither mapped nor silent: one warning per keyword (ADR 0044,
/// point 2), in the order of their first appearance.
class Ignored {
 public:
  void add(const Line& line) {
    const std::string& k = line.words[0];
    auto it = counts_.find(k);
    if (it == counts_.end()) {
      order_.push_back(k);
      counts_.emplace(k, std::make_pair(std::size_t{1}, line.number));
    } else {
      ++it->second.first;
    }
  }
  void report(std::vector<Diagnostic>& out) const {
    for (const std::string& k : order_) {
      const auto& [count, first] = counts_.at(k);
      out.push_back(warning("zmx.keyword_ignored", first,
                            "keyword " + k + " ignored (" + std::to_string(count) +
                                (count == 1 ? " time)" : " times)")));
    }
  }
  static Diagnostic warning(DiagnosticCode code, std::size_t line, const std::string& message) {
    return {model::Severity::Warning, line == 0 ? "file" : "line " + std::to_string(line), message,
            std::string(code.str())};
  }

 private:
  std::vector<std::string> order_;
  std::map<std::string, std::pair<std::size_t, std::size_t>> counts_;
};

class Importer {
 public:
  Importer(const std::vector<Line>& lines, const material::MaterialLibrary& lib)
      : lines_(lines), lib_(lib) {}

  ZmxImport run();

 private:
  void read_header(const Line& l);
  void read_surface_line(ZSurface& s, const Line& l);
  void read_trailer(std::size_t from);
  [[nodiscard]] double length(double v) const { return v * scale_; }
  std::string resolve_glass(const ZSurface& s);
  model::ShapeStack shape(const ZSurface& s, bool& powered) const;
  std::optional<model::Aperture> aperture(const ZSurface& s);
  model::System build();

  const std::vector<Line>& lines_;
  const material::MaterialLibrary& lib_;
  std::vector<Diagnostic> warnings_;
  Ignored ignored_;

  std::string name_;
  double scale_ = 1.0;
  std::optional<std::string> aperture_kind_;  ///< ENPD, FNUM, OBNA, FLOA
  double aperture_value_ = 0.0;
  std::vector<std::string> gcat_;
  std::optional<Line> ftyp_;
  std::map<std::string, Line> field_lines_;  ///< XFLN/YFLN/FWGN or XFLD/YFLD/FWGT
  std::map<int, Line> wavm_;
  std::map<std::string, Line> wave_lists_;  ///< WAVN/WWGN or WAVL/WWGT
  std::vector<Line> vignetting_;
  int pwav_ = 1;
  std::vector<ZSurface> surfaces_;
  std::map<std::string, std::size_t> glass_first_line_;
  bool stop_from_diam_rule_ = false;  ///< the stop radius is 1.1 v (rule 3 of point 8)
  std::size_t stop_diam_line_ = 0;
};

void Importer::read_header(const Line& l) {
  const std::string& k = l.words[0];
  if (k == "NAME") {
    name_ = l.rest;
  } else if (k == "MODE") {
    if (l.words.size() < 2 || l.words[1] != "SEQ") {
      fail("zmx.mode_unsupported", l.number, "only sequential files (MODE SEQ) are imported");
    }
  } else if (k == "UNIT") {
    static const std::map<std::string, double> kUnits{
        {"MM", 1.0}, {"CM", 10.0}, {"IN", 25.4}, {"METER", 1000.0}};
    const auto it = l.words.size() > 1 ? kUnits.find(l.words[1]) : kUnits.end();
    if (it == kUnits.end()) {
      fail("zmx.unit_unsupported", l.number,
           "lens unit " + (l.words.size() > 1 ? l.words[1] : std::string()) + " is not supported");
    }
    scale_ = it->second;
  } else if (k == "ENPD" || k == "FNUM" || k == "OBNA" || k == "FLOA") {
    if (k != "FLOA") {
      aperture_value_ = number(l, 1);
      if ((k == "FNUM" || k == "OBNA") && l.words.size() > 2 && integer(l, 2) != 0) {
        fail("zmx.aperture_unsupported", l.number,
             k == "FNUM" ? "paraxial working F/# is not supported"
                         : "object cone angle is not supported");
      }
    }
    if (aperture_kind_) {
      fail("zmx.parse", l.number,
           "a second system aperture (" + k + " after " + *aperture_kind_ + ")");
    }
    aperture_kind_ = k;
  } else if (k == "ENVD") {
    if (number(l, 1) != 20.0 || number(l, 2) != 1.0) {
      fail("zmx.environment_unsupported", l.number,
           "only 20 degC and 1 atm are supported (wavelengths are relative to that air)");
    }
  } else if (k == "GCAT") {
    gcat_.assign(l.words.begin() + 1, l.words.end());
  } else if (k == "FTYP") {
    ftyp_ = l;
  } else if (k == "XFLN" || k == "YFLN" || k == "FWGN" || k == "XFLD" || k == "YFLD" ||
             k == "FWGT") {
    field_lines_[k] = l;
  } else if (k == "VDXN" || k == "VDYN" || k == "VCXN" || k == "VCYN" || k == "VANN" ||
             k == "ZVDX" || k == "ZVDY" || k == "ZVCX" || k == "ZVCY" || k == "ZVAN") {
    vignetting_.push_back(l);
  } else if (k == "WAVM") {
    wavm_[integer(l, 1)] = l;
  } else if (k == "WAVN" || k == "WWGN" || k == "WAVL" || k == "WWGT") {
    wave_lists_[k] = l;
  } else if (k == "PWAV") {
    pwav_ = integer(l, 1);
  } else if (!silent(k)) {
    ignored_.add(l);
  }
}

void Importer::read_surface_line(ZSurface& s, const Line& l) {
  const std::string& k = l.words[0];
  if (k == "TYPE") {
    s.type = l.words.size() > 1 ? l.words[1] : std::string();
    s.type_line = l.number;
  } else if (k == "CURV") {
    s.curv = {number(l, 1), l.number};
  } else if (k == "CONI") {
    s.coni = {number(l, 1), l.number};
  } else if (k == "PARM") {
    s.parm[integer(l, 1)] = {number(l, 2), l.number};
  } else if (k == "DISZ") {
    if (l.words.size() > 1 && l.words[1] == "INFINITY") {
      s.infinite = true;
      s.disz = {0.0, l.number};
    } else {
      s.disz = {number(l, 1), l.number};
    }
  } else if (k == "GLAS") {
    s.glass = l.words.size() > 1 ? l.words[1] : std::string();
    if (upper(s.glass) == "AIR") s.glass.clear();
    if (upper(s.glass) == "MIRROR") s.glass = "MIRROR";
    s.glass_line = l.number;
    s.glass_words = l.words;
  } else if (k == "DIAM") {
    s.diam = {number(l, 1), l.number};
    s.diam_flag = l.words.size() > 2 ? integer(l, 2) : 0;
  } else if (k == "CLAP") {
    s.clap = std::make_pair(number(l, 1), number(l, 2));
    s.clap_line = l.number;
  } else if (k == "FLAP") {
    s.flap = std::make_pair(number(l, 1), number(l, 2));
    s.flap_line = l.number;
  } else if (k == "STOP") {
    s.stop = true;
    s.stop_line = l.number;
  } else if (k == "OBDC") {
    if (number(l, 1) != 0.0 || number(l, 2) != 0.0) {
      fail("zmx.aperture_decenter", l.number, "decentred apertures are not supported");
    }
  } else if (k == "SQAP") {
    s.sqap = std::make_pair(number(l, 1), number(l, 2));
    s.sqap_line = l.number;
  } else if (k == "ELAP") {
    s.elap = std::make_pair(number(l, 1), number(l, 2));
    s.elap_line = l.number;
  } else if (k == "OBSC" || k == "SQOB" || k == "ELOB" || k == "SPID" || k == "USAP") {
    fail("zmx.aperture_unsupported", l.number, "aperture " + k + " is not supported");
  } else if (k == "VDSZ" || k == "FIMP") {
    // VDSZ (every surface of the OISC files, always zeros) and FIMP (newer files, no value) are
    // silent only in that form; their meaning is not documented (ADR 0044, addendum PR C).
    for (std::size_t i = 1; i < l.words.size(); ++i) {
      if (k == "FIMP" || number(l, i) != 0.0) {
        fail("zmx.keyword_unsupported", l.number,
             k + " with a value is not supported (its meaning is not documented)");
      }
    }
  } else if (!silent(k)) {
    // In a surface block every keyword may change the optics (e.g. the tilt and decenter tab,
    // whose keywords no free reader documents): an error, not a warning (addendum PR C).
    fail("zmx.keyword_unsupported", l.number, "unknown keyword " + k + " in a surface");
  }
}

/// Merit block, TOL lines and the configuration editor after the last surface (ADR 0044,
/// point 2 and 9).
void Importer::read_trailer(std::size_t from) {
  // A SURF block after a line that ended the surfaces would be skipped with the merit block.
  for (std::size_t k = from; k < lines_.size(); ++k) {
    if (lines_[k].words[0] == "SURF") {
      fail("zmx.parse", lines_[k].number,
           "SURF after the line " + std::to_string(lines_[from].number) +
               ", which ended the surfaces (a line that is not indented)");
    }
  }
  std::size_t i = from;
  while (i < lines_.size() && lines_[i].words[0] != "TOL" && lines_[i].words[0] != "MNUM") ++i;
  while (i < lines_.size() && lines_[i].words[0] != "MNUM") ++i;
  if (i == lines_.size()) return;
  const Line& mnum = lines_[i];
  if (integer(mnum, 1) != 1) {
    fail("zmx.multi_config", mnum.number, "more than one configuration is not supported");
  }
  for (++i; i < lines_.size(); ++i) {
    const Line& l = lines_[i];
    const std::string& k = l.words[0];
    if (k == "MOFF") continue;
    if (k == "THIC" || k == "CRVT") {
      const int surface = integer(l, 1);
      const double v = number(l, 3);
      const auto it = std::find_if(surfaces_.begin(), surfaces_.end(),
                                   [&](const ZSurface& s) { return s.number == surface; });
      const bool ok = it != surfaces_.end() &&
                      (k == "THIC" ? !it->infinite && it->disz.value == v : it->curv.value == v);
      if (!ok) {
        fail("zmx.config_mismatch", l.number,
             k + " of surface " + std::to_string(surface) + " differs from the surface data");
      }
      continue;
    }
    ignored_.add(l);
  }
}

std::string Importer::resolve_glass(const ZSurface& s) {
  const std::vector<std::string>& w = s.glass_words;
  if ((w.size() > 2 && w[2] == "1") || s.glass == "___BLANK" || mil_number(s.glass)) {
    fail("zmx.model_glass", s.glass_line,
         "model glass " + s.glass + " cannot be imported (its dispersion is not published)");
  }
  const std::vector<std::string> loaded = lib_.catalogs();
  std::string searched;
  for (const std::string& g : gcat_) {
    searched += (searched.empty() ? "" : ", ") + g;
    const auto cat = std::find_if(loaded.begin(), loaded.end(),
                                  [&](const std::string& c) { return upper(c) == upper(g); });
    if (cat == loaded.end()) continue;
    const std::string ref = *cat + ":" + s.glass;
    try {
      const auto material = lib_.resolve(ref);
      if (!glass_first_line_.contains(ref)) {
        glass_first_line_.emplace(ref, s.glass_line);
        const auto* catalog = dynamic_cast<const material::CatalogMaterial*>(material.get());
        const std::optional<material::SchottThermalCoefficients> thermal =
            catalog != nullptr ? catalog->thermal() : std::nullopt;
        if (thermal.has_value() && thermal->reference_temperature_c != 20.0) {
          warnings_.push_back(
              Ignored::warning("zmx.glass_reference_temperature", s.glass_line,
                               ref + " has its catalogue data at " +
                                   std::to_string(thermal->reference_temperature_c) +
                                   " degC, not at the 20 degC of the file"));
        }
      }
      return ref;
    } catch (const material::UnknownMaterial&) {
      continue;
    }
  }
  fail("zmx.glass_unknown", s.glass_line,
       "glass " + s.glass + " not found in the catalogues " +
           (searched.empty() ? std::string("(no GCAT)") : searched));
}

model::ShapeStack Importer::shape(const ZSurface& s, bool& powered) const {
  model::ShapeStack out;
  const double c = s.curv.value / scale_;
  powered = c != 0.0;
  if (s.type == "EVENASPH") {
    // alpha_i (r^(2i)) in L^(1 - 2i): times 1/s for (2i - 1) times, no pow (ADR 0044, point 3)
    std::vector<double> alpha(8, 0.0);
    for (const auto& [i, v] : s.parm) {
      if (i < 1 || i > 8) {
        if (v.value != 0.0) {
          fail("zmx.parse", v.line, "EVENASPH has the parameters 1 to 8 only");
        }
        continue;
      }
      double x = v.value;
      const double inv = 1.0 / scale_;
      for (int k = 0; k < 2 * i - 1; ++k) x *= inv;
      alpha[static_cast<std::size_t>(i - 1)] = x;
    }
    while (!alpha.empty() && alpha.back() == 0.0) alpha.pop_back();
    if (!alpha.empty()) {
      if (alpha[0] != 0.0) {
        fail("zmx.asphere_r2", s.parm.at(1).line,
             "EVENASPH with alpha1 (r^2) != 0 is not imported: its paraxial reading is not "
             "documented");
      }
      powered = true;
      if (c != 0.0) {
        model::EvenAsphere a{1.0 / c, s.coni.value, {}};
        for (std::size_t i = 1; i < alpha.size(); ++i) a.coefficients.emplace_back(alpha[i]);
        out.base = a;
      } else {
        // c = 0: odd_asphere with A_{2i} = alpha_i (ADR 0044, point 6); the conic has no effect
        model::OddAsphere a{0.0, 0.0, {0.0, 0.0, 0.0}};
        for (std::size_t i = 1; i < alpha.size(); ++i) {
          if (i > 1) a.coefficients.emplace_back(0.0);
          a.coefficients.emplace_back(alpha[i]);
        }
        out.base = a;
      }
      return out;
    }
  } else if (s.type == "TILTSURF" || s.type == "PARAXIAL") {
    // planes (UG "Tilted", "Paraxial"); the tilt goes into the pose, the lens into the interaction
    if (c != 0.0) {
      fail("zmx.surface_unsupported", s.curv.line, s.type + " is a plane, CURV must be 0");
    }
    powered = false;
    out.base = model::Plane{};
    return out;
  } else if (s.type != "STANDARD") {
    fail("zmx.surface_unsupported", s.type_line, "surface type " + s.type + " is not supported");
  }
  if (c != 0.0) {
    out.base = model::Conic{1.0 / c, s.coni.value};
  } else {
    out.base = model::Plane{};
  }
  return out;
}

std::optional<model::Aperture> Importer::aperture(const ZSurface& s) {
  const int count = static_cast<int>(s.clap.has_value()) + static_cast<int>(s.flap.has_value()) +
                    static_cast<int>(s.sqap.has_value()) + static_cast<int>(s.elap.has_value());
  if (count > 1) {
    fail("zmx.aperture_conflict", std::max({s.clap_line, s.flap_line, s.sqap_line, s.elap_line}),
         "a surface with more than one aperture (CLAP, FLAP, SQAP, ELAP)");
  }
  // SQAP and ELAP: half widths in x and y, local to the surface (UG "Aperture surface
  // properties")
  if (s.sqap) {
    const auto [x, y] = *s.sqap;
    if (!(x > 0.0) || !(y > 0.0))
      fail("zmx.aperture_invalid", s.sqap_line, "SQAP needs widths > 0");
    return model::RectangularAperture{length(x), length(y)};
  }
  if (s.elap) {
    const auto [x, y] = *s.elap;
    if (!(x > 0.0) || !(y > 0.0))
      fail("zmx.aperture_invalid", s.elap_line, "ELAP needs widths > 0");
    return model::EllipticalAperture{length(x), length(y)};
  }
  if (s.clap) {
    const auto [rmin, rmax] = *s.clap;
    if (!(rmax > 0.0) || !(rmax > rmin) || rmin < 0.0) {
      fail("zmx.aperture_invalid", s.clap_line, "CLAP needs 0 <= rmin < rmax and rmax > 0");
    }
    return model::CircularAperture{length(rmax), length(rmin)};
  }
  if (s.flap) {
    const auto [first, rmax] = *s.flap;
    if (first != 0.0) {
      fail("zmx.aperture_unsupported", s.flap_line,
           "FLAP with a first field != 0 (its meaning is not documented)");
    }
    if (!(rmax > 0.0)) fail("zmx.aperture_invalid", s.flap_line, "FLAP with radius 0");
    if (s.diam_flag == 1 || s.diam_flag == 2) return model::CircularAperture{length(rmax), 0.0};
    warnings_.push_back(Ignored::warning(
        "zmx.aperture_dropped", s.flap_line,
        "FLAP next to an automatic semi-diameter dropped (it would clip real marginal rays)"));
  }
  return std::nullopt;
}

model::System Importer::build() {
  model::System sys;
  sys.name = name_;
  sys.environment = model::Environment{20.0, 1.0, "AIR"};

  // Wavelengths and fields (ADR 0044, point 4).
  if (!ftyp_) fail("zmx.parse", 0, "no FTYP line");
  const Line& ft = *ftyp_;
  const int ftype = integer(ft, 1);
  if (ftype == 3 || ftype == 4) {
    fail("zmx.field_unsupported", ft.number,
         ftype == 3 ? "real image height is not supported" : "theodolite angle is not supported");
  }
  if (ftype < 0 || ftype > 4) fail("zmx.parse", ft.number, "unknown field type");
  if (ft.words.size() > 2 && integer(ft, 2) != 0) {
    fail("zmx.aperture_unsupported", ft.number, "object space telecentric is not supported");
  }
  if (ft.words.size() > 7 && integer(ft, 7) != 0) {
    warnings_.push_back(Ignored::warning(
        "zmx.afocal_image_space", ft.number,
        "afocal image space ignored (it only changes the units of analyses in OpticStudio)"));
  }
  const int n_fields = integer(ft, 3);
  const int n_waves = integer(ft, 4);
  if (n_fields < 1 || n_waves < 1) fail("zmx.parse", ft.number, "no fields or wavelengths");

  for (const Line& v : vignetting_) {
    const auto count = static_cast<std::size_t>(n_fields);
    for (std::size_t i = 1; i <= count && i < v.words.size(); ++i) {
      if (number(v, i) != 0.0) {
        fail("zmx.vignetting_unsupported", v.number, "vignetting factors are not supported");
      }
    }
  }

  const auto pick = [&](const char* a, const char* b) -> const Line* {
    if (field_lines_.contains(a)) return &field_lines_.at(a);
    if (field_lines_.contains(b)) return &field_lines_.at(b);
    return nullptr;
  };
  const Line* xs = pick("XFLN", "XFLD");
  const Line* ys = pick("YFLN", "YFLD");
  const Line* ws = pick("FWGN", "FWGT");
  if (xs == nullptr || ys == nullptr) fail("zmx.parse", 0, "no field coordinates");
  static constexpr std::array<model::FieldType, 3> kTypes{model::FieldType::AngleDeg,
                                                          model::FieldType::ObjectHeight,
                                                          model::FieldType::ParaxialImageHeight};
  sys.fields.type = kTypes[static_cast<std::size_t>(ftype)];
  const double field_scale = ftype == 0 ? 1.0 : scale_;
  for (std::size_t at = 1; at <= static_cast<std::size_t>(n_fields); ++at) {
    model::Field f;
    f.x = number(*xs, at) * field_scale;
    f.y = number(*ys, at) * field_scale;
    f.weight = ws != nullptr ? number(*ws, at) : 1.0;
    sys.fields.points.push_back(f);
  }

  for (int i = 1; i <= n_waves; ++i) {
    double air = 0.0;
    double weight = 1.0;
    if (wavm_.contains(i)) {
      air = number(wavm_.at(i), 2);
      weight = number(wavm_.at(i), 3);
    } else {
      const auto list = [&](const char* a, const char* b) -> const Line* {
        if (wave_lists_.contains(a)) return &wave_lists_.at(a);
        if (wave_lists_.contains(b)) return &wave_lists_.at(b);
        return nullptr;
      };
      const Line* wl = list("WAVN", "WAVL");
      const Line* ww = list("WWGN", "WWGT");
      if (wl == nullptr) fail("zmx.parse", 0, "no wavelengths");
      air = number(*wl, static_cast<std::size_t>(i));
      if (ww != nullptr) weight = number(*ww, static_cast<std::size_t>(i));
    }
    sys.wavelengths.push_back({compile::line_vacuum_wavelength_um(air), weight, i == pwav_});
  }
  if (pwav_ < 1 || pwav_ > n_waves) fail("zmx.parse", 0, "PWAV outside the wavelengths");
  const double ref_um = sys.wavelengths[static_cast<std::size_t>(pwav_ - 1)].um;
  const double n_air = material::ciddor_air_index(ref_um, 20.0, 1.0);

  // System aperture.
  if (!aperture_kind_) fail("zmx.parse", 0, "no system aperture (ENPD, FNUM, OBNA or FLOA)");
  const std::string& kind = *aperture_kind_;
  using AT = model::SystemApertureType;
  if (kind == "ENPD") {
    sys.aperture.type = AT::EntrancePupilDiameter;
    sys.aperture.value = length(aperture_value_);
  } else if (kind == "FNUM") {
    sys.aperture.type = AT::ImageSpaceFNumber;
    sys.aperture.value = aperture_value_ / n_air;
  } else if (kind == "OBNA") {
    sys.aperture.type = AT::ObjectSpaceNA;
    sys.aperture.value = n_air * aperture_value_;
    sys.aperture.reading = model::ApertureReading::Real;
  } else {
    sys.aperture.type = AT::StopSize;
  }

  // Object (ADR 0044, point 5).
  if (surfaces_.size() < 3)
    fail("zmx.parse", 0, "the file needs an object, a surface and an image");
  const ZSurface& obj = surfaces_.front();
  if (obj.type != "STANDARD" || obj.curv.value != 0.0) {
    fail("zmx.surface_unsupported", obj.curv.value != 0.0 ? obj.curv.line : obj.type_line,
         "the object surface must be a plane STANDARD surface");
  }
  std::optional<std::string> object_medium;
  if (obj.infinite) {
    if (!obj.glass.empty()) {
      fail("zmx.object_medium_infinite", obj.glass_line,
           "a glass in object space needs a finite object");
    }
  } else {
    sys.object.at_infinity = false;
    sys.object.distance = length(obj.disz.value);
    if (!obj.glass.empty()) {
      if (is_mirror(obj)) fail("zmx.object_medium_unsupported", obj.glass_line, "MIRROR at 0");
      object_medium = resolve_glass(obj);
      sys.object.medium = object_medium;
      // the first surface that is not a coordinate break (addendum PR C)
      std::size_t first_real = 1;
      while (first_real + 1 < surfaces_.size() && surfaces_[first_real].type == "COORDBRK") {
        ++first_real;
      }
      const ZSurface& s1 = surfaces_[first_real];
      if (s1.glass.empty() || is_mirror(s1)) {
        fail("zmx.object_medium_unsupported", obj.glass_line,
             "surface 1 after a glass object space must start a lens");
      }
    }
  }

  // Exactly one STOP, on a surface between the object and the image.
  const auto stops =
      std::count_if(surfaces_.begin(), surfaces_.end(), [](const ZSurface& z) { return z.stop; });
  if (surfaces_.front().stop || surfaces_.back().stop) {
    const ZSurface& z = surfaces_.front().stop ? surfaces_.front() : surfaces_.back();
    fail("zmx.parse", z.stop_line, "the stop cannot be the object or the image surface");
  }
  if (stops == 0) fail("zmx.parse", 0, "no STOP surface");

  // Surfaces -> tree.
  std::vector<model::Node> children;
  int lenses = 0;
  int mirrors = 0;
  int thins = 0;
  FrameChain chain;  // from the frame of the last placed surface to the next one
  bool first = true;
  model::Element lens;  // the lens being built, valid while lens_open
  bool lens_open = false;
  std::vector<std::string> lens_media;
  bool before_glass = object_medium.has_value();
  std::string medium_name = object_medium ? obj.glass : std::string();  // glass name in the file
  const auto place = [&](model::Element& e) {
    e.pose = chain_pose(chain, std::nullopt);
    if (!first) e.pose.reference = model::PoseReference::RelativeToPreceding;
    first = false;
  };
  // After a placed surface: its thickness, and the untilt of a tilted plane (it keeps the frame).
  const auto advance = [&](double disz, const std::optional<std::array<double, 3>>& tilt) {
    chain = FrameChain{};
    if (tilt) chain.untilt = rotation_of(*tilt).transpose();
    chain.thickness(disz);
  };
  const auto aperture_line = [](const ZSurface& z) {
    for (const std::size_t line : {z.clap_line, z.flap_line, z.sqap_line, z.elap_line}) {
      if (line != 0) return line;
    }
    return std::size_t{0};
  };
  const auto stop_radius = [&](const ZSurface& s) {
    const double v = length(s.diam.value);
    if (kind == "FLOA" || s.diam_flag == 1 || s.diam_flag == 2) return v;
    stop_from_diam_rule_ = true;
    stop_diam_line_ = s.diam.line;
    return 1.1 * v;
  };
  const auto make_surface = [&](const ZSurface& s, bool& powered) {
    model::Surface out;
    out.id = model::SurfaceId("S" + std::to_string(s.number));
    out.shape = shape(s, powered);
    out.aperture = aperture(s);
    if (s.stop) {
      if (out.aperture && !std::holds_alternative<model::CircularAperture>(*out.aperture)) {
        // the stop rules of point 8 and the pupils need a stop radius (addendum PR C)
        fail("zmx.aperture_unsupported", aperture_line(s),
             "a rectangular or elliptical aperture on the stop is not supported");
      }
      if (out.aperture) {
        const auto* c = std::get_if<model::CircularAperture>(&*out.aperture);
        if (kind == "FLOA" && c != nullptr && c->radius != length(s.diam.value)) {
          fail("zmx.stop_radius_conflict", s.clap ? s.clap_line : s.flap_line,
               "with FLOA the aperture of the stop must equal its semi-diameter");
        }
      } else {
        if (!(s.diam.value > 0.0)) {
          fail("zmx.stop_radius_missing", s.diam.line != 0 ? s.diam.line : s.stop_line,
               "the stop has no semi-diameter and no aperture");
        }
        out.aperture = model::CircularAperture{stop_radius(s), 0.0};
      }
    }
    return out;
  };
  const auto close_lens = [&] {
    if (std::all_of(lens_media.begin(), lens_media.end(),
                    [&](const std::string& m) { return m == lens_media.front(); })) {
      lens.material = lens_media.front();
    } else {
      lens.segment_materials = lens_media;
    }
    children.push_back(model::Node{lens});
    lens_open = false;
    lens_media.clear();
  };

  const auto parm = [](const ZSurface& z, int i) {
    const auto it = z.parm.find(i);
    return it == z.parm.end() ? 0.0 : it->second.value;
  };

  for (std::size_t k = 1; k < surfaces_.size(); ++k) {
    const ZSurface& s = surfaces_[k];
    const bool image = k + 1 == surfaces_.size();
    const double disz = s.infinite ? 0.0 : length(s.disz.value);
    if (s.infinite && !image) fail("zmx.parse", s.disz.line, "INFINITY only at the object");
    bool powered = false;
    const std::string id = "S" + std::to_string(s.number);

    // Tilted plane (UG "Tilted"): the tilt goes into the pose of this surface only. Its aperture
    // would be tested in the untilted frame there and in the tilted plane here (addendum PR C).
    std::optional<std::array<double, 3>> tilt;
    if (s.type == "TILTSURF") {
      if (s.stop) {
        fail("zmx.aperture_unsupported", s.stop_line,
             "a tilted plane as the stop is not supported");
      }
      if (aperture_line(s) != 0) {
        fail("zmx.aperture_unsupported", aperture_line(s),
             "an aperture on a tilted plane is not supported (OpticStudio tests it untilted)");
      }
      if (parm(s, 1) != 0.0 || parm(s, 2) != 0.0) tilt = plane_tilt_deg(parm(s, 1), parm(s, 2));
    }

    if (image) {
      if (s.type == "COORDBRK" || s.type == "PARAXIAL") {
        fail("zmx.surface_unsupported", s.type_line, s.type + " cannot be the image surface");
      }
      if (before_glass) {
        const ZSurface& prev = surfaces_[k - 1];
        fail("zmx.image_in_glass", prev.glass_line, "glass in front of the image surface");
      }
      model::Element det{"image", model::ElementKind::Detector, {}, std::nullopt, {}};
      det.surfaces.push_back(make_surface(s, powered));
      if (tilt) {
        det.surfaces.back().pose.rotation_deg = {(*tilt)[0], (*tilt)[1], (*tilt)[2]};
        warnings_.push_back(Ignored::warning(
            "zmx.tilted_image", s.type_line,
            "tilted image plane: image coordinates are in the tilted plane, not projected on the "
            "untilted frame as in OpticStudio"));
      }
      place(det);
      children.push_back(model::Node{det});
      break;
    }

    // Coordinate break (UG "Coordinate Break"): no surface, a step of the frame chain; it keeps
    // the medium of the surface before it.
    if (s.type == "COORDBRK") {
      if (s.stop) {
        fail("zmx.coordbrk_surface", s.stop_line, "a coordinate break cannot be the stop");
      }
      if (aperture_line(s) != 0) {
        fail("zmx.coordbrk_surface", aperture_line(s),
             "an aperture on a coordinate break is not supported");
      }
      if (s.glass_line != 0 && (is_mirror(s) || s.glass != medium_name)) {
        fail("zmx.parse", s.glass_line,
             "a coordinate break takes the medium of the surface before it");
      }
      const FrameStep step{true,
                           0.0,
                           length(parm(s, 1)),
                           length(parm(s, 2)),
                           {parm(s, 3), parm(s, 4), parm(s, 5)},
                           parm(s, 6) != 0.0};
      if (step.dx != 0.0 || step.dy != 0.0 || step.tilt != std::array<double, 3>{0.0, 0.0, 0.0}) {
        chain.steps.push_back(step);
      }
      chain.thickness(disz);
      continue;
    }

    if (is_mirror(s) && before_glass) {
      fail("zmx.mirror_in_glass", s.glass_line, "a mirror inside glass (Mangin) is not supported");
    }

    // Paraxial surface (UG "Paraxial"): an ideal lens on a thin element (ADR 0031). In air on
    // both sides only: in glass its focal distance would be n(lambda) f; f = 0 has no power and
    // is a plane like any other.
    if (s.type == "PARAXIAL" && !is_mirror(s) && parm(s, 1) != 0.0) {
      if (before_glass || !s.glass.empty()) {
        fail("zmx.paraxial_medium", s.glass.empty() ? s.type_line : s.glass_line,
             "a paraxial surface needs air on both sides");
      }
      const double mode = parm(s, 2);
      const std::size_t mode_line = s.parm.contains(2) ? s.parm.at(2).line : s.type_line;
      if (mode != 0.0 && mode != 1.0 && mode != 2.0 && mode != 3.0) {
        fail("zmx.parse", mode_line, "the OPD mode of a paraxial surface is 0, 1, 2 or 3");
      }
      if (mode != 2.0) {
        warnings_.push_back(Ignored::warning(
            "zmx.paraxial_opd", mode_line,
            "OPD mode of the paraxial surface imported for infinite conjugates: the rays are the "
            "same, the OPD at a finite conjugate is not"));
      }
      ++thins;
      model::Element t{
          "T" + std::to_string(thins), model::ElementKind::ThinElement, {}, std::nullopt, {}};
      model::Surface surface = make_surface(s, powered);
      surface.interaction = model::IdealLens{length(parm(s, 1)), std::nullopt};
      t.surfaces.push_back(std::move(surface));
      place(t);
      children.push_back(model::Node{t});
      if (s.stop) sys.aperture.stop = model::SurfaceId(id);
      before_glass = false;
      medium_name.clear();
      advance(disz, std::nullopt);
      continue;
    }
    if (s.type == "PARAXIAL" && is_mirror(s)) {
      fail("zmx.surface_unsupported", s.glass_line, "a paraxial surface cannot be a mirror");
    }

    const bool after_glass = !s.glass.empty() && !is_mirror(s);
    if (lens_open) {
      model::Surface surface = make_surface(s, powered);
      surface.pose =
          chain_pose(chain, tilt ? std::optional<math::Mat3>(rotation_of(*tilt)) : std::nullopt);
      surface.pose.reference = model::PoseReference::RelativeToPreceding;
      lens.surfaces.push_back(std::move(surface));
      if (s.stop) sys.aperture.stop = model::SurfaceId("S" + std::to_string(s.number));
      if (after_glass) {
        lens_media.push_back(resolve_glass(s));
      } else {
        close_lens();
      }
    } else if (is_mirror(s)) {
      ++mirrors;
      model::Element m{
          "M" + std::to_string(mirrors), model::ElementKind::Mirror, {}, std::nullopt, {}};
      model::Surface surface = make_surface(s, powered);
      surface.interaction = model::IdealMirror{};
      if (tilt) surface.pose.rotation_deg = {(*tilt)[0], (*tilt)[1], (*tilt)[2]};
      m.surfaces.push_back(std::move(surface));
      place(m);
      children.push_back(model::Node{m});
      if (s.stop) sys.aperture.stop = model::SurfaceId("S" + std::to_string(s.number));
    } else if (after_glass) {
      ++lenses;
      lens = model::Element{
          "L" + std::to_string(lenses), model::ElementKind::Lens, {}, std::nullopt, {}};
      lens_open = true;
      lens.surfaces.push_back(make_surface(s, powered));
      if (tilt) lens.surfaces.back().pose.rotation_deg = {(*tilt)[0], (*tilt)[1], (*tilt)[2]};
      lens_media.push_back(resolve_glass(s));
      place(lens);
      if (s.stop) sys.aperture.stop = model::SurfaceId("S" + std::to_string(s.number));
    } else if (s.stop) {
      model::Element st{"stop", model::ElementKind::Stop, {}, std::nullopt, {}};
      st.surfaces.push_back(make_surface(s, powered));
      place(st);
      children.push_back(model::Node{st});
    } else {
      model::Surface surface = make_surface(s, powered);
      if (surface.aperture) {
        ++thins;
        model::Element t{
            "T" + std::to_string(thins), model::ElementKind::ThinElement, {}, std::nullopt, {}};
        t.surfaces.push_back(std::move(surface));
        place(t);
        children.push_back(model::Node{t});
      } else {
        if (powered) {
          warnings_.push_back(Ignored::warning(
              "zmx.surface_dropped", s.curv.line != 0 ? s.curv.line : s.line,
              "surface " + std::to_string(s.number) +
                  " between air and air with power dropped (it refracts nothing)"));
        }
        chain.thickness(disz);
        continue;
      }
    }
    before_glass = after_glass ? true : (is_mirror(s) ? before_glass : false);
    if (after_glass) {
      medium_name = s.glass;
    } else if (!is_mirror(s)) {
      medium_name.clear();
    }
    advance(disz, tilt);
  }

  sys.root = model::Assembly{"system", {}, std::move(children)};
  sys.paths.push_back(model::Path{"main", true, {}});
  return sys;
}

ZmxImport Importer::run() {
  // Header up to the first SURF, then the SURF blocks, then the trailer (ADR 0044, point 2).
  std::size_t i = 0;
  for (; i < lines_.size() && lines_[i].words[0] != "SURF"; ++i) read_header(lines_[i]);
  while (i < lines_.size() && lines_[i].words[0] == "SURF") {
    ZSurface s;
    s.number = integer(lines_[i], 1);
    s.line = lines_[i].number;
    if (s.number != static_cast<int>(surfaces_.size())) {
      fail("zmx.parse", s.line, "surfaces must be numbered 0, 1, 2, ... in order");
    }
    for (++i; i < lines_.size() && lines_[i].indented; ++i) read_surface_line(s, lines_[i]);
    surfaces_.push_back(std::move(s));
  }
  read_trailer(i);

  model::System sys = build();
  ZmxImport out;
  ignored_.report(out.warnings);
  out.warnings.insert(out.warnings.end(), warnings_.begin(), warnings_.end());

  try {
    const compile::CompiledSystem cs = compile::compile(sys, lib_);
    for (const Diagnostic& d : cs.diagnostics()) out.warnings.push_back(d);
    std::uint16_t ref = 0;
    for (std::size_t w = 0; w < sys.wavelengths.size(); ++w) {
      if (sys.wavelengths[w].reference) ref = static_cast<std::uint16_t>(w);
    }
    // A folded, tilted or decentred path has no paraxial model (ADR 0030, addendum #256): the
    // system is imported, analyses with paraxial aiming reject it (#357).
    bool paraxial = true;
    try {
      (void)paraxial::first_order(cs, compile::PathId{0}, ref);
    } catch (const paraxial::ParaxialError& e) {
      paraxial = false;
      out.warnings.push_back(Ignored::warning(
          "zmx.paraxial_unavailable", 0,
          std::string("no paraxial model, analyses with paraxial aiming need given start rays: ") +
              e.what()));
    }
    if (stop_from_diam_rule_ && paraxial) {
      const double r_s = trace::paraxial_stop_radius(cs, compile::PathId{0}, ref);
      double radius = 0.0;
      const auto visit = [&](const auto& self, const model::Assembly& a) -> void {
        for (const model::Node& n : a.children) {
          if (const auto* e = std::get_if<model::Element>(&n.value)) {
            for (const model::Surface& s : e->surfaces) {
              const bool is_stop = sys.aperture.stop ? s.id == *sys.aperture.stop
                                                     : e->kind == model::ElementKind::Stop;
              if (is_stop && s.aperture) {
                radius = std::get<model::CircularAperture>(*s.aperture).radius;
              }
            }
          } else {
            self(self, std::get<model::Assembly>(n.value));
          }
        }
      };
      visit(visit, sys.root);
      if (r_s > radius) {
        out.warnings.push_back(Ignored::warning(
            "zmx.stop_radius_small", stop_diam_line_,
            "paraxial stop radius " + std::to_string(r_s) + " mm exceeds 1.1 DIAM = " +
                std::to_string(radius) + " mm; real marginal rays are vignetted at the stop"));
      }
    }
  } catch (const ZmxError&) {
    throw;
  } catch (const std::exception& e) {
    fail("zmx.invalid_result", 0, std::string("the imported system is invalid: ") + e.what());
  }
  out.system = std::move(sys);
  return out;
}

}  // namespace

ZmxImport import_zmx(std::string_view bytes, const material::MaterialLibrary& materials) {
  if (bytes.size() >= 2 && static_cast<unsigned char>(bytes[0]) == 0xFE &&
      static_cast<unsigned char>(bytes[1]) == 0xFF) {
    fail("zmx.encoding_unsupported", 0, "UTF-16 big endian is not supported");
  }
  std::string text;
  try {
    text = material::decode_agf_text(bytes, "zmx");
  } catch (const std::exception& e) {
    fail("zmx.encoding_unsupported", 0, e.what());
  }
  const std::vector<Line> lines = split_lines(text);
  return Importer(lines, materials).run();
}

ZmxImport load_zmx(const std::filesystem::path& file, const material::MaterialLibrary& materials) {
  const std::ifstream in(file, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read " + file.string());
  std::ostringstream buffer;
  buffer << in.rdbuf();
  return import_zmx(buffer.str(), materials);
}

}  // namespace rtt::io
