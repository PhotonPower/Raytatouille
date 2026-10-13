#pragma once

// A system with every kind of tolerance target (ADR 0035; #253): shared by the tests of the
// class enumeration and of validate.

#include <optional>
#include <string>
#include <utility>

#include "rtt/model/model.hpp"
#include "test_support.hpp"

namespace rtt::model::test {

inline Surface surface(std::string id, BaseShape base, double z) {
  Surface s;
  s.id = SurfaceId(std::move(id));
  s.shape.base = std::move(base);
  s.pose = Pose::along_z(z);
  return s;
}

inline ParameterRow row(std::string name, double v) {
  ParameterRow r;
  r.name = std::move(name);
  r.form = v;
  return r;
}

inline ToleranceTarget at(ToleranceAnchor anchor, std::string name) {
  ToleranceTarget t;
  t.anchor = anchor;
  t.name = std::move(name);
  return t;
}

inline ToleranceValueOverride value(ToleranceAnchor anchor, std::string name, ToleranceKind kind) {
  return {at(anchor, std::move(name)), kind, 0.01, ToleranceDistribution::Normal};
}

inline ToleranceItem item(std::string id, ToleranceTarget target, ToleranceKind kind) {
  return {std::move(id), std::move(target), kind, 0.01, ToleranceDistribution::Normal};
}

/// A stop, an assembly "objective" with a cemented lens L1 (two catalog glasses, a conic at
/// L1.S3) and a CONST singlet L2 (even asphere, a plane back surface placed by the row T2), a
/// metal fold mirror and a detector placed by the row BFD.
inline System toleranced() {
  System s = make_singlet();
  s.parameters = {row("T2", 4.0), row("BFD", 100.0)};
  const Element stop = element(s, 0);
  Element image = element(s, 2);
  image.pose.position[2] = Param::bound("BFD");

  Element l1{"L1",
             ElementKind::Lens,
             Pose::along_z(5.0),
             std::nullopt,
             {surface("L1.S1", Conic{Param(61.5), Param(0.0)}, 0.0),
              surface("L1.S2", Conic{Param(-35.2), Param(0.0)}, 6.0),
              surface("L1.S3", Conic{Param(-136.9), Param(-0.5)}, 9.0)}};
  l1.segment_materials = {"SCHOTT:N-BK7", "SCHOTT:F2"};
  Surface back = surface("L2.S2", Plane{}, 0.0);
  back.pose.position[2] = Param::bound("T2");
  Element l2{"L2",
             ElementKind::Lens,
             Pose::along_z(20.0),
             "CONST:1.5",
             {surface("L2.S1", EvenAsphere{Param(40.0), Param(0.0), {Param(1e-6)}}, 0.0), back}};
  Assembly objective;
  objective.name = "objective";
  objective.children = {Node{l1}, Node{l2}};
  Element fold{"fold",
               ElementKind::Mirror,
               Pose::along_z(40.0),
               "CONST:1.2,7.26",
               {surface("FOLD", Plane{}, 0.0)}};

  s.root.children = {Node{stop}, Node{objective}, Node{fold}, Node{image}};
  s.tolerances.tolerance_class = ToleranceClass::Normal;
  return s;
}

}  // namespace rtt::model::test
