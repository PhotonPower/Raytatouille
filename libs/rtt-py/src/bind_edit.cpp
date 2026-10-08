// Editing a system from Python (#82, ADR 0024): the edit form as Python objects, JSON pointers,
// stable anchors and JSON Patch. The Editor with undo and redo is pure Python
// (raytatouille/edit.py) on top of apply_patch_json.

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <variant>
#include <vector>

#include "bindings.hpp"
#include "rtt/io/edit.hpp"
#include "rtt/model/system.hpp"
#include "rtt/model/validate.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace rtt::py {
namespace {

/// RFC 6901 Sec. 3 and 4: "" or "/" followed by tokens; "~1" decodes to "/", then "~0" to "~".
/// @throws std::invalid_argument (ValueError) for invalid syntax
std::vector<std::string> pointer_tokens(std::string_view pointer) {
  std::vector<std::string> out;
  if (pointer.empty()) return out;
  if (pointer.front() != '/') {
    throw std::invalid_argument("JSON pointer '" + std::string(pointer) +
                                "' must be empty or start with '/'");
  }
  std::string token;
  for (std::size_t i = 1; i <= pointer.size(); ++i) {
    if (i == pointer.size() || pointer[i] == '/') {
      out.push_back(token);
      token.clear();
    } else if (pointer[i] == '~') {
      if (i + 1 < pointer.size() && (pointer[i + 1] == '0' || pointer[i + 1] == '1')) {
        token += pointer[i + 1] == '0' ? '~' : '/';
        ++i;
      } else {
        throw std::invalid_argument("JSON pointer '" + std::string(pointer) +
                                    "': '~' must be followed by '0' or '1'");
      }
    } else {
      token += pointer[i];
    }
  }
  return out;
}

/// Walks the tree of `node` and finds the pointer of a surface id or a node name.
struct Finder {
  std::string_view surface;
  std::string_view node;

  std::optional<std::string> assembly(const model::Assembly& a, const std::string& at) const {
    if (!node.empty() && a.name == node) return at;
    for (std::size_t i = 0; i < a.children.size(); ++i) {
      const std::string child = at + "/children/" + std::to_string(i);
      const std::optional<std::string> found = std::visit(
          [&](const auto& n) -> std::optional<std::string> {
            using T = std::decay_t<decltype(n)>;
            if constexpr (std::is_same_v<T, model::Assembly>) {
              return assembly(n, child);
            } else {
              return element(n, child);
            }
          },
          a.children[i].value);
      if (found) return found;
    }
    return std::nullopt;
  }

  std::optional<std::string> element(const model::Element& e, const std::string& at) const {
    if (!node.empty() && e.name == node) return at;
    if (surface.empty()) return std::nullopt;
    for (std::size_t k = 0; k < e.surfaces.size(); ++k) {
      if (e.surfaces[k].id.str() == surface) return at + "/surfaces/" + std::to_string(k);
    }
    return std::nullopt;
  }
};

}  // namespace

nb::object edit_dict(const model::System& s) {
  return nb::module_::import_("json").attr("loads")(io::to_edit_json(s));
}

nb::object json_at(const model::System& s, std::string_view pointer) {
  const std::vector<std::string> tokens = pointer_tokens(pointer);
  nb::object at = edit_dict(s);
  for (const std::string& t : tokens) {
    if (nb::isinstance<nb::dict>(at)) {
      const nb::dict d = nb::borrow<nb::dict>(at);
      const nb::str key(t.data(), t.size());  // with its size: a key may contain NUL
      if (!d.contains(key)) throw nb::key_error(std::string(pointer).c_str());
      at = d[key];
    } else if (nb::isinstance<nb::list>(at)) {
      const nb::list l = nb::borrow<nb::list>(at);
      const bool digits = !t.empty() && t.find_first_not_of("0123456789") == std::string::npos;
      if (!digits || (t.size() > 1 && t.front() == '0')) {
        throw std::invalid_argument("JSON pointer '" + std::string(pointer) + "': '" + t +
                                    "' is no array index (digits without a leading zero)");
      }
      // More than 18 digits is a valid index but beyond any list.
      if (t.size() > 18 || std::stoull(t) >= nb::len(l)) {
        throw nb::key_error(std::string(pointer).c_str());
      }
      at = l[std::stoull(t)];
    } else {
      throw nb::key_error(std::string(pointer).c_str());
    }
  }
  return at;
}

std::optional<std::string> locate_surface(const model::System& s, std::string_view id) {
  return Finder{id, {}}.assembly(s.root, "/root");
}

std::optional<std::string> locate_node(const model::System& s, std::string_view name) {
  return Finder{{}, name}.assembly(s.root, "/root");
}

std::optional<std::string> locate_path(const model::System& s, std::string_view name) {
  for (std::size_t i = 0; i < s.paths.size(); ++i) {
    if (s.paths[i].name == name) return "/paths/" + std::to_string(i);
  }
  return std::nullopt;
}

void bind_edit(nb::module_& m) {
  m.def(
      "apply_patch_json",
      [](const model::System& system, std::string_view patch_json, bool structure_only) {
        io::PatchResult r = io::apply_patch(
            system, patch_json,
            structure_only ? io::EditCheck::StructureOnly : io::EditCheck::NoNewErrors);
        return std::make_tuple(std::move(r.system), std::move(r.inverse), std::move(r.diagnostics));
      },
      "system"_a, "patch_json"_a, "structure_only"_a,
      "Applies an RFC 6902 patch document (JSON text) to the edit form of `system` (ADR 0024) "
      "and returns (new system, inverse patch as JSON text, all diagnostics of the new "
      "system). Use raytatouille.apply_patch or raytatouille.Editor instead.\n\nRaises "
      "EditError (code, location, op_index, diagnostics) if the patch cannot be applied or adds "
      "an error of validate() (unless structure_only).");
}

}  // namespace rtt::py
