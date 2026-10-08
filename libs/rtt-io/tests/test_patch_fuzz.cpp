// Differential fuzz test of the inverse patches (ADR 0024 point 3), ported from the Python
// scaffold of the second review of ADR 0024: random JSON documents, random patches of 1 to 3
// operations; for every accepted patch the inverse must restore the original document.
// Deterministic (rule 7): std::mt19937 with a fixed seed and choices by rng() % n, so the
// sequence is the same with every standard library.

#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <nlohmann/json.hpp>
#include <random>
#include <string>
#include <vector>

#include "json_patch.hpp"

using nlohmann::json;
using rtt::io::detail::apply_json_patch;
using rtt::io::detail::PatchApplied;
using rtt::io::detail::PatchFailure;

namespace {

class Fuzz {
 public:
  explicit Fuzz(std::uint32_t seed) : rng_(seed) {}

  std::size_t below(std::size_t n) { return static_cast<std::size_t>(rng_() % n); }
  bool chance(unsigned percent) { return below(100) < percent; }

  json document(int depth) {
    if (depth > 2 || chance(30)) {
      switch (below(5)) {
        case 0:
          return 1;
        case 1:
          return 2.5;
        case 2:
          return "s";
        case 3:
          return true;
        default:
          return nullptr;
      }
    }
    if (chance(50)) {
      json o = json::object();
      std::string keys = "abc";
      const std::size_t count = below(4);
      for (std::size_t i = 0; i < count; ++i) {
        const std::size_t k = below(keys.size());
        o[std::string(1, keys[k])] = document(depth + 1);
        keys.erase(k, 1);
      }
      return o;
    }
    json a = json::array();
    const std::size_t count = below(4);
    for (std::size_t i = 0; i < count; ++i) a.push_back(document(depth + 1));
    return a;
  }

  json operation(const json& doc) {
    std::vector<std::string> pointers;
    collect(doc, "", pointers);
    const std::string p = pointers[below(pointers.size())];
    const std::string q = pointers[below(pointers.size())];
    static const std::vector<std::string> kOps = {"add", "remove", "replace", "move", "copy"};
    static const std::vector<std::string> kLast = {"a", "b", "0", "1", "-"};
    const std::string op = kOps[below(kOps.size())];
    const std::string target = chance(50) ? p + "/" + kLast[below(kLast.size())] : p;
    json o = {{"op", op}};
    o["path"] = (op == "add" || op == "move" || op == "copy") ? target : p;
    if (op == "add" || op == "replace") o["value"] = document(2);
    if (op == "move" || op == "copy") o["from"] = q;
    return o;
  }

 private:
  static void collect(const json& d, const std::string& pointer, std::vector<std::string>& out) {
    out.push_back(pointer);
    if (d.is_object()) {
      for (const auto& item : d.items()) collect(item.value(), pointer + "/" + item.key(), out);
    } else if (d.is_array()) {
      for (std::size_t i = 0; i < d.size(); ++i)
        collect(d[i], pointer + "/" + std::to_string(i), out);
    }
  }

  std::mt19937 rng_;
};

}  // namespace

TEST_CASE("every accepted patch is undone by its inverse (fuzz, fixed seed)", "[io][patch][fuzz]") {
  constexpr std::size_t kPatches = 30000;
  Fuzz fuzz(1);
  std::size_t accepted = 0;
  const auto start = std::chrono::steady_clock::now();
  for (std::size_t n = 0; n < kPatches; ++n) {
    const json doc = {{"root", fuzz.document(0)}};
    json patch = json::array();
    const std::size_t ops = 1 + fuzz.below(3);
    for (std::size_t i = 0; i < ops; ++i) patch.push_back(fuzz.operation(doc));
    PatchApplied r;
    try {
      r = apply_json_patch(doc, patch);
    } catch (const PatchFailure&) {
      continue;
    }
    ++accepted;
    INFO("document " << doc.dump() << "\npatch " << patch.dump() << "\ninverse "
                     << r.inverse.dump());
    REQUIRE(apply_json_patch(r.document, r.inverse).document == doc);
    // Counter-checks of the forward result (ADR 0024 point 3): nlohmann's own patch agrees on
    // every patch accepted here (the cases where it deviates are rejected above), and its diff
    // applied by this engine reaches the same document. One more deviation of nlohmann: a move
    // with the root "" as "from" (allowed only with the path "", RFC 6902 Sec. 4.4, and without
    // effect) throws "JSON pointer has no parent" there, so such patches are left out of the
    // first counter-check.
    bool move_from_root = false;
    for (const json& op : patch) {
      move_from_root = move_from_root || (op["op"] == "move" && op["from"] == "");
    }
    if (!move_from_root) REQUIRE(doc.patch(patch) == r.document);
    REQUIRE(apply_json_patch(doc, json::diff(doc, r.document)).document == r.document);
  }
  const double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  // The sequence is deterministic, so the number of accepted patches is fixed: a change of it
  // means a changed behaviour of the engine (13 670 of 30 000, 46 %). The run time is reported,
  // not checked, so that slow sanitizer builds do not fail on it.
  INFO(accepted << " of " << kPatches << " patches accepted in " << seconds << " s");
  CHECK(accepted == 13670);
}
