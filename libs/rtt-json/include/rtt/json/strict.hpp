#pragma once

/// @file strict.hpp
/// Strict JSON parsing shared by the file readers of rtt-io (system files) and rtt-coating
/// (coating catalogues), ADR 0008 and ADR 0020: like nlohmann::json::parse, but a key that
/// appears twice in one object is an error instead of silently keeping the last value, and every
/// parse failure (syntax, number overflow) is one exception type with a JSON pointer.

#include <cstddef>
#include <functional>
#include <nlohmann/json.hpp>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace rtt::json {

/// Parse failure with a JSON pointer (RFC 6901) to the offending value; the pointer is empty for
/// errors of the whole text (syntax, number overflow). `what()` is "<pointer>: <message>" with
/// "/" for an empty pointer, the same text as rtt::io::ParseError.
class StrictParseError : public std::runtime_error {
 public:
  /// @param pointer JSON pointer, e.g. "/coatings/0/layers/0/thickness_um", or "" for the text
  /// @param message description without the pointer
  StrictParseError(std::string pointer, std::string message)
      : std::runtime_error((pointer.empty() ? std::string("/") : pointer) + ": " + message),
        pointer_(std::move(pointer)),
        message_(std::move(message)) {}

  /// JSON pointer (RFC 6901) to the offending key or value; empty for the whole text.
  [[nodiscard]] const std::string& pointer() const noexcept { return pointer_; }
  /// Description without the pointer, e.g. "duplicate key 'a'" or "invalid JSON: ...".
  [[nodiscard]] const std::string& message() const noexcept { return message_; }

 private:
  std::string pointer_;
  std::string message_;
};

/// An object key as JSON pointer reference token (RFC 6901, Sec. 3: '~' -> "~0", '/' -> "~1").
[[nodiscard]] inline std::string pointer_token(std::string_view key) {
  std::string out;
  out.reserve(key.size());
  for (const char c : key) {
    if (c == '~') {
      out += "~0";
    } else if (c == '/') {
      out += "~1";
    } else {
      out += c;
    }
  }
  return out;
}

/// Parses JSON text and rejects duplicate keys at any depth.
/// @param text JSON text (UTF-8)
/// @return the parsed document
/// @throws StrictParseError with the pointer of the second occurrence for a duplicate key, and
///         with an empty pointer for a syntax error or a number that overflows double
///         (nlohmann reports it as out_of_range 406; JSON has no other non-finite numbers)
[[nodiscard]] inline nlohmann::json parse_strict(std::string_view text) {
  using Json = nlohmann::json;
  using Event = Json::parse_event_t;
  // One frame per open object or array: keys seen and current key (object), number of
  // finished elements (array), and the pointer of the container itself.
  struct Frame {
    bool object = true;
    std::set<std::string, std::less<>> keys;
    std::string key;
    std::size_t index = 0;
    std::string pointer;
  };
  std::vector<Frame> stack;
  const auto child_pointer = [&stack]() -> std::string {
    if (stack.empty()) return "";
    const Frame& f = stack.back();
    return f.pointer + "/" + (f.object ? pointer_token(f.key) : std::to_string(f.index));
  };
  const auto element_finished = [&stack]() {
    if (!stack.empty() && !stack.back().object) ++stack.back().index;
  };
  const Json::parser_callback_t callback = [&](int /*depth*/, Event event, Json& parsed) {
    switch (event) {
      case Event::object_start:
      case Event::array_start: {
        Frame frame;
        frame.object = event == Event::object_start;
        frame.pointer = child_pointer();
        stack.push_back(std::move(frame));
        break;
      }
      case Event::key: {
        Frame& f = stack.back();
        std::string key = parsed.get<std::string>();
        if (!f.keys.insert(key).second) {
          throw StrictParseError(f.pointer + "/" + pointer_token(key),
                                 "duplicate key '" + key + "'");
        }
        f.key = std::move(key);
        break;
      }
      case Event::object_end:
      case Event::array_end:
        stack.pop_back();
        element_finished();
        break;
      case Event::value:
        element_finished();
        break;
    }
    return true;  // keep every value
  };
  try {
    return Json::parse(text.begin(), text.end(), callback);
  } catch (const Json::exception& e) {
    throw StrictParseError("", std::string("invalid JSON: ") + e.what());
  }
}

}  // namespace rtt::json
