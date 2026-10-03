#include "rtt/model/element.hpp"

namespace rtt::model {

bool Assembly::operator==(const Assembly& other) const {
  return name == other.name && pose == other.pose && children == other.children;
}

bool Node::operator==(const Node& other) const {
  return value == other.value;
}

}  // namespace rtt::model
