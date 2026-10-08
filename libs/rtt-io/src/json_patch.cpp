#include "json_patch.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "rtt/json/strict.hpp"

namespace rtt::io::detail {

PatchFailure::PatchFailure(std::string code,
                           std::string pointer,
                           std::optional<std::size_t> op_index,
                           const std::string& message)
    : std::runtime_error((pointer.empty() ? std::string("/") : pointer) + ": " + message),
      code_(std::move(code)),
      pointer_(std::move(pointer)),
      op_index_(op_index),
      message_(message) {}

namespace {

using Json = nlohmann::json;
using Tokens = std::vector<std::string>;

/// One operation of the patch being applied; reports failures with its index.
class Operation {
 public:
  Operation(std::size_t index, const Json& op) : index_(index), op_(op) {
    if (!op.is_object()) fail("edit.patch_invalid", "", "an operation must be a JSON object");
    name_ = member_string("op", "");
    path_text_ = member_string("path", "");
    path_ = tokens(path_text_);
  }

  [[noreturn]] void fail(const char* code,
                         const std::string& pointer,
                         const std::string& message) const {
    throw PatchFailure(code, pointer, index_,
                       "operation " + std::to_string(index_) + ": " + message);
  }

  [[nodiscard]] const std::string& name() const noexcept { return name_; }
  [[nodiscard]] const Tokens& path() const noexcept { return path_; }
  [[nodiscard]] const std::string& path_text() const noexcept { return path_text_; }

  /// "value" (add, replace, test).
  [[nodiscard]] const Json& value() const {
    const auto it = op_.find("value");
    if (it == op_.end()) fail("edit.patch_invalid", path_text_, "'" + name_ + "' needs a value");
    return *it;
  }

  /// "from" (move, copy) as text.
  [[nodiscard]] std::string from_text() const { return member_string("from", path_text_); }

  /// RFC 6901 Sec. 3 and 4: a pointer is "" or "/" followed by reference tokens; "~1" decodes to
  /// "/", then "~0" to "~"; any other "~" is invalid.
  [[nodiscard]] Tokens tokens(const std::string& pointer) const {
    Tokens out;
    if (pointer.empty()) return out;
    if (pointer.front() != '/') {
      fail("edit.patch_invalid", pointer, "a JSON pointer must be empty or start with '/'");
    }
    std::string token;
    for (std::size_t i = 1; i <= pointer.size(); ++i) {
      if (i == pointer.size() || pointer[i] == '/') {
        out.push_back(token);
        token.clear();
      } else if (pointer[i] == '~') {
        if (i + 1 < pointer.size() && pointer[i + 1] == '0') {
          token += '~';
        } else if (i + 1 < pointer.size() && pointer[i + 1] == '1') {
          token += '/';
        } else {
          fail("edit.patch_invalid", pointer, "'~' must be followed by '0' or '1'");
        }
        ++i;
      } else {
        token += pointer[i];
      }
    }
    return out;
  }

 private:
  [[nodiscard]] std::string member_string(const char* key, const std::string& where) const {
    const auto it = op_.find(key);
    if (it == op_.end() || !it->is_string()) {
      fail("edit.patch_invalid", where, std::string("needs a string member \"") + key + "\"");
    }
    return it->get<std::string>();
  }

  std::size_t index_;
  const Json& op_;
  std::string name_;
  std::string path_text_;
  Tokens path_;
};

/// The tokens as a JSON pointer (RFC 6901, Sec. 3: "~" -> "~0", "/" -> "~1").
std::string pointer(const Tokens& tokens) {
  std::string out;
  for (const std::string& t : tokens) out += "/" + rtt::json::pointer_token(t);
  return out;
}

/// An array index token (RFC 6901, Sec. 4): "0" or digits without a leading zero, and "-" for
/// the end where `end_allowed` (the target of add, move and copy). A token of another form is an
/// invalid pointer (edit.patch_invalid, ADR 0024 point 4); an index above the size, or equal to
/// it where the end is not allowed, does not exist (edit.path_not_found).
std::size_t array_index(const std::string& token,
                        std::size_t size,
                        bool end_allowed,
                        const Operation& op,
                        const std::string& text) {
  if (token == "-") {
    if (!end_allowed) {
      op.fail("edit.patch_invalid", text, "'-' is only allowed as the last token of a target");
    }
    return size;
  }
  const bool digits = !token.empty() && token.find_first_not_of("0123456789") == std::string::npos;
  if (!digits || (token.size() > 1 && token.front() == '0')) {
    op.fail("edit.patch_invalid", text,
            "'" + token + "' is no array index (digits without a leading zero)");
  }
  std::size_t value = 0;
  if (token.size() <= 18) {
    for (const char c : token) value = value * 10 + static_cast<std::size_t>(c - '0');
  }
  if (token.size() > 18 || value > size || (value == size && !end_allowed)) {
    op.fail("edit.path_not_found", text, "array index " + token + " out of range");
  }
  return value;
}

/// The value at the first `count` tokens, or nullptr if an object member is missing or a
/// primitive value is in the way.
Json* find(Json& doc,
           const Tokens& tokens,
           std::size_t count,
           const Operation& op,
           const std::string& text) {
  Json* at = &doc;
  for (std::size_t i = 0; i < count; ++i) {
    if (at->is_object()) {
      const auto it = at->find(tokens[i]);
      if (it == at->end()) return nullptr;
      at = &*it;
    } else if (at->is_array()) {
      at = &(*at)[array_index(tokens[i], at->size(), false, op, text)];
    } else {
      return nullptr;
    }
  }
  return at;
}

Json& existing(Json& doc, const Tokens& tokens, const Operation& op, const std::string& text) {
  Json* v = find(doc, tokens, tokens.size(), op, text);
  if (v == nullptr) op.fail("edit.path_not_found", text, "no value at '" + text + "'");
  return *v;
}

Json& container(Json& doc, const Tokens& tokens, const Operation& op, const std::string& text) {
  Json* parent = find(doc, tokens, tokens.size() - 1, op, text);
  if (parent == nullptr || !(parent->is_object() || parent->is_array())) {
    op.fail("edit.path_not_found", text, "no object or array to hold '" + text + "'");
  }
  return *parent;
}

Json inverse_op(const char* name, const Tokens& path, const Json* value = nullptr) {
  Json o = {{"op", name}, {"path", pointer(path)}};
  if (value != nullptr) o["value"] = *value;
  return o;
}

/// RFC 6902 Sec. 4.1. Returns the inverse operations in application order.
std::vector<Json> add(
    Json& doc, const Tokens& path, Json value, const Operation& op, const std::string& text) {
  if (path.empty()) {  // the root: the value replaces the document
    const Json old = std::move(doc);
    doc = std::move(value);
    return {inverse_op("replace", path, &old)};
  }
  Json& parent = container(doc, path, op, text);
  if (parent.is_object()) {
    const std::string& key = path.back();
    std::vector<Json> inverse;
    const auto it = parent.find(key);
    if (it != parent.end()) {
      inverse.push_back(inverse_op("replace", path, &*it));
    } else {
      inverse.push_back(inverse_op("remove", path));
    }
    parent[key] = std::move(value);
    return inverse;
  }
  const std::size_t k = array_index(path.back(), parent.size(), true, op, text);
  parent.insert(parent.begin() + static_cast<std::ptrdiff_t>(k), std::move(value));
  Tokens resolved = path;
  resolved.back() = std::to_string(k);
  return {inverse_op("remove", resolved)};
}

/// RFC 6902 Sec. 4.2. Returns the removed value and its path with the index resolved.
std::pair<Json, Tokens> remove(Json& doc,
                               const Tokens& path,
                               const Operation& op,
                               const std::string& text) {
  if (path.empty()) op.fail("edit.patch_invalid", text, "the whole document cannot be removed");
  Json& parent = container(doc, path, op, text);
  if (parent.is_object()) {
    const auto it = parent.find(path.back());
    if (it == parent.end()) op.fail("edit.path_not_found", text, "no value at '" + text + "'");
    Json value = std::move(*it);
    parent.erase(it);
    return {std::move(value), path};
  }
  const std::size_t k = array_index(path.back(), parent.size(), false, op, text);
  Json value = std::move(parent[k]);
  parent.erase(k);
  Tokens resolved = path;
  resolved.back() = std::to_string(k);
  return {std::move(value), std::move(resolved)};
}

bool proper_prefix(const Tokens& a, const Tokens& b) {
  if (a.size() >= b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (a[i] != b[i]) return false;
  }
  return true;
}

/// Applies one operation to `doc`; returns its inverse operations in application order.
std::vector<Json> apply(Json& doc, const Operation& op) {
  const std::string& name = op.name();
  const Tokens& path = op.path();
  const std::string& text = op.path_text();
  if (name == "add") return add(doc, path, op.value(), op, text);
  if (name == "remove") {
    auto [value, resolved] = remove(doc, path, op, text);
    return {inverse_op("add", resolved, &value)};
  }
  if (name == "replace") {  // Sec. 4.3: the target must exist
    const Json& value = op.value();
    Json& target = existing(doc, path, op, text);
    const Json old = std::move(target);
    target = value;
    return {inverse_op("replace", path, &old)};
  }
  if (name == "move" || name == "copy") {
    const std::string from_text = op.from_text();
    const Tokens from = op.tokens(from_text);
    if (name == "copy") return add(doc, path, existing(doc, from, op, from_text), op, text);
    // Sec. 4.4
    if (proper_prefix(from, path)) {
      op.fail("edit.patch_invalid", text, "'" + from_text + "' cannot move into its own child");
    }
    if (from == path) {
      (void)existing(doc, from, op, from_text);
      return {};
    }
    auto [value, resolved_from] = remove(doc, from, op, from_text);
    const Json moved = value;
    std::vector<Json> inverse = add(doc, path, std::move(value), op, text);
    inverse.push_back(inverse_op("add", resolved_from, &moved));
    return inverse;
  }
  if (name == "test") {  // Sec. 4.6: equality by value, members in any order
    if (existing(doc, path, op, text) != op.value()) {
      op.fail("edit.test_failed", text, "test failed at '" + text + "'");
    }
    return {};
  }
  op.fail("edit.patch_invalid", text, "unknown operation '" + name + "'");
}

}  // namespace

PatchApplied apply_json_patch(const nlohmann::json& document, const nlohmann::json& patch) {
  if (!patch.is_array()) {
    throw PatchFailure("edit.patch_invalid", "", std::nullopt,
                       "a JSON Patch document must be an array of operations");
  }
  PatchApplied out{document, Json::array()};
  std::vector<std::vector<Json>> inverses;
  inverses.reserve(patch.size());
  for (std::size_t i = 0; i < patch.size(); ++i) {
    const Operation op(i, patch[i]);
    inverses.push_back(apply(out.document, op));
  }
  for (auto it = inverses.rbegin(); it != inverses.rend(); ++it) {
    for (Json& o : *it) out.inverse.push_back(std::move(o));
  }
  return out;
}

}  // namespace rtt::io::detail
