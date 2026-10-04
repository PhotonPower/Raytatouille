#include "rtt/model/element.hpp"

namespace rtt::model {

std::optional<std::string> Element::segment_material(std::size_t segment) const {
  if (surfaces.size() < 2 || segment > surfaces.size() - 2) return std::nullopt;
  if (!segment_materials.empty()) {
    if (segment < segment_materials.size()) return segment_materials[segment];
    return std::nullopt;
  }
  return material;
}

bool Assembly::operator==(const Assembly& other) const {
  return name == other.name && pose == other.pose && children == other.children;
}

bool Node::operator==(const Node& other) const {
  return value == other.value;
}

}  // namespace rtt::model
