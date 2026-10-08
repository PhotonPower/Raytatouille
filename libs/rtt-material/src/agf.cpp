#include "rtt/material/agf.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "rtt/diagnostics/codes.hpp"
#include "rtt/material/air.hpp"

namespace rtt::material {
namespace {

/// Appends the UTF-8 encoding of a Unicode scalar value.
void append_utf8(std::string& out, char32_t c) {
  const auto byte = [](std::uint32_t v) {
    return static_cast<char>(static_cast<unsigned char>(v));
  };
  const auto v = static_cast<std::uint32_t>(c);
  if (v < 0x80) {
    out += byte(v);
  } else if (v < 0x800) {
    out += byte(0xC0 | (v >> 6));
    out += byte(0x80 | (v & 0x3F));
  } else if (v < 0x10000) {
    out += byte(0xE0 | (v >> 12));
    out += byte(0x80 | ((v >> 6) & 0x3F));
    out += byte(0x80 | (v & 0x3F));
  } else {
    out += byte(0xF0 | (v >> 18));
    out += byte(0x80 | ((v >> 12) & 0x3F));
    out += byte(0x80 | ((v >> 6) & 0x3F));
    out += byte(0x80 | (v & 0x3F));
  }
}

/// Splits a line at spaces and tabs.
std::vector<std::string_view> split(std::string_view line) {
  std::vector<std::string_view> items;
  std::size_t pos = 0;
  while (pos < line.size()) {
    while (pos < line.size() && (line[pos] == ' ' || line[pos] == '\t')) ++pos;
    const std::size_t start = pos;
    while (pos < line.size() && line[pos] != ' ' && line[pos] != '\t') ++pos;
    if (pos > start) items.push_back(line.substr(start, pos - start));
  }
  return items;
}

/// Reads one AGF file and turns it into an AgfCatalog.
class Parser {
 public:
  Parser(std::string name, std::string file)
      : catalog_{std::move(name), std::move(file), {}, {}, {}} {}

  AgfCatalog run(std::string_view text) {
    std::size_t line_no = 0;
    std::size_t pos = 0;
    while (pos <= text.size()) {
      const std::size_t end = std::min(text.find('\n', pos), text.size());
      std::string_view line = text.substr(pos, end - pos);
      if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
      ++line_no;
      parse_line(line, line_no);
      if (end == text.size()) break;
      pos = end + 1;
    }
    if (stray_) fail(stray_->line, "unknown record '" + stray_->word + "'");  // last line
    finish_glass();
    merge_duplicates();
    // Duplicates are found after the last line; keep all warnings in file order.
    std::stable_sort(catalog_.warnings.begin(), catalog_.warnings.end(),
                     [](const LoadWarning& x, const LoadWarning& y) { return x.line < y.line; });
    return std::move(catalog_);
  }

 private:
  [[noreturn]] void fail(std::size_t line, const std::string& message) const {
    throw AgfError(catalog_.file, line, message);
  }

  void warn(diagnostics::DiagnosticCode code, std::size_t line, std::string message) {
    catalog_.warnings.push_back({std::string(code.str()), catalog_.file, line, std::move(message)});
  }

  static bool is_record(std::string_view m) {
    static constexpr std::array<std::string_view, 11> kRecords = {
        "CC", "NM", "CD", "LD", "TD", "ED", "GC", "MD", "OD", "IT", "BD"};
    return std::find(kRecords.begin(), kRecords.end(), m) != kRecords.end();
  }

  static bool is_word(std::string_view m) {
    return !m.empty() && std::all_of(m.begin(), m.end(), [](char c) {
      return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
    });
  }

  /// R1/R2 (#71): a repeated glass with the same data (all records equal, line aside) is read
  /// once; with other data both blocks stay and the glass is ambiguous.
  void merge_duplicates() {
    std::vector<AgfGlass> kept;
    for (AgfGlass& g : catalog_.glasses) {
      const AgfGlass* same = nullptr;
      const AgfGlass* other = nullptr;
      for (const AgfGlass& k : kept) {
        if (k.name != g.name) continue;
        AgfGlass a = k;
        a.line = g.line;
        if (a == g) {
          same = &k;
        } else if (other == nullptr) {
          other = &k;
        }
      }
      if (same != nullptr) {
        warn("agf.duplicate_glass", g.line,
             "glass " + g.name + " repeats the block at line " + std::to_string(same->line) +
                 " with identical data; the repetition is ignored");
        continue;
      }
      if (other != nullptr) {
        warn("agf.duplicate_glass_conflict", g.line,
             "glass " + g.name + " is defined again with different data (records " +
                 differing_records(*other, g) + " differ from the block at line " +
                 std::to_string(other->line) + "); the glass is ambiguous and cannot be resolved");
      }
      kept.push_back(std::move(g));
    }
    catalog_.glasses = std::move(kept);
  }

  /// Names of the records whose data differ, e.g. "NM, CD".
  static std::string differing_records(const AgfGlass& a, const AgfGlass& b) {
    std::vector<std::string> names;
    if (a.formula != b.formula || a.nd != b.nd || a.vd != b.vd ||
        a.exclude_substitution != b.exclude_substitution || a.status != b.status ||
        a.melt_frequency != b.melt_frequency) {
      names.emplace_back("NM");
    }
    if (a.comment != b.comment) names.emplace_back("GC");
    if (a.extra != b.extra) names.emplace_back("ED");
    if (a.coefficients != b.coefficients) names.emplace_back("CD");
    if (a.thermal != b.thermal) names.emplace_back("TD");
    if (a.mechanical != b.mechanical) names.emplace_back("MD");
    if (a.other != b.other) names.emplace_back("OD");
    if (a.range != b.range) names.emplace_back("LD");
    if (a.transmission != b.transmission) names.emplace_back("IT");
    std::string text;
    for (const std::string& n : names) text += (text.empty() ? "" : ", ") + n;
    return text;
  }

  double number(std::string_view text, std::size_t line) const {
    double value = 0.0;
    const char* const first = text.data();
    const char* const last = first + text.size();
    const auto [ptr, ec] = std::from_chars(first, last, value);
    if (ec != std::errc{} || ptr != last || !std::isfinite(value)) {
      fail(line, "'" + std::string(text) + "' is not a finite number");
    }
    return value;
  }

  std::vector<double> numbers(const std::vector<std::string_view>& items, std::size_t line) const {
    std::vector<double> values;
    values.reserve(items.size() - 1);
    for (std::size_t i = 1; i < items.size(); ++i) values.push_back(number(items[i], line));
    return values;
  }

  AgfGlass& current(std::string_view mnemonic, std::size_t line) {
    if (!glass_) fail(line, "record " + std::string(mnemonic) + " before the first NM record");
    return *glass_;
  }

  /// Text after the mnemonic, without leading and trailing spaces and tabs (CC, GC).
  static std::string text_after(std::string_view line, std::string_view mnemonic) {
    const std::size_t first = line.find_first_not_of(" \t", line.find(mnemonic) + mnemonic.size());
    if (first == std::string_view::npos) return {};
    const std::size_t last = line.find_last_not_of(" \t");
    return std::string(line.substr(first, last - first + 1));
  }

  /// A record of exactly `count` values where "_" means "not available" (MD, OD; placeholder of
  /// the manufacturers' files, #85).
  std::vector<std::optional<double>> values_with_placeholders(
      const std::vector<std::string_view>& items, std::size_t count, std::size_t line) const {
    if (items.size() - 1 != count) {
      fail(line, std::string(items[0]) + " needs exactly " + std::to_string(count) +
                     " values (\"_\" for a missing one), got " + std::to_string(items.size() - 1));
    }
    std::vector<std::optional<double>> values;
    values.reserve(count);
    for (std::size_t i = 1; i < items.size(); ++i) {
      if (items[i] == "_") {
        values.emplace_back(std::nullopt);
      } else {
        values.emplace_back(number(items[i], line));
      }
    }
    return values;
  }

  /// A resistance class of OD: "_" or "-" (not available), a number (one class), or a range
  /// "a-b" of two numbers with 0 <= a <= b as SCHOTT writes it (e.g. "1-2", docs/quellen.md,
  /// #85).
  std::optional<AgfClassRange> class_range(std::string_view item, std::size_t line) const {
    if (item == "_" || item == "-") return std::nullopt;
    const auto parse = [](std::string_view text) -> std::optional<double> {
      double value = 0.0;
      const char* const first = text.data();
      const char* const last = first + text.size();
      const auto [ptr, ec] = std::from_chars(first, last, value);
      if (text.empty() || ec != std::errc{} || ptr != last || !std::isfinite(value)) {
        return std::nullopt;
      }
      return value;
    };
    // One class: >= 0, or -1 for "not available" (format description; kept as written).
    if (const std::optional<double> v = parse(item); v && (*v >= 0.0 || *v == -1.0)) {
      return AgfClassRange{*v, *v};
    }
    const std::size_t dash = item.find('-', 1);
    if (dash != std::string_view::npos) {
      const std::optional<double> low = parse(item.substr(0, dash));
      const std::optional<double> high = parse(item.substr(dash + 1));
      if (low && high && *low >= 0.0 && *low <= *high) return AgfClassRange{*low, *high};
    }
    fail(line, "OD: resistance class '" + std::string(item) +
                   "' is not a number, \"_\", \"-\" or a range \"a-b\" with 0 <= a <= b");
  }

  /// An NM extra (exclude sub, status, melt freq): "_" or "-" is "not available", otherwise an
  /// integer, in [lo, hi] if a range is given (#85).
  std::optional<int> nm_extra(std::string_view item,
                              const char* what,
                              std::optional<std::pair<int, int>> range,
                              std::size_t line) const {
    if (item == "_" || item == "-") return std::nullopt;
    const double v = number(item, line);
    constexpr double kIntLimit = 2147483647.0;
    if (v != std::floor(v) || std::abs(v) > kIntLimit) {
      fail(line, std::string("NM: ") + what + " '" + std::string(item) + "' is not an integer");
    }
    if (range && (v < range->first || v > range->second)) {
      fail(line, std::string("NM: ") + what + " '" + std::string(item) + "' is not an integer in " +
                     std::to_string(range->first) + ".." + std::to_string(range->second));
    }
    return static_cast<int>(v);
  }

  /// Record that a continuation line may extend: the CD or TD record directly above it.
  enum class Continuable : std::uint8_t { kNone, kCd, kTd };

  /// True if an item starts like a number (digit, '+', '-' or '.'): a continuation line (#42).
  static bool starts_with_number(std::string_view item) {
    return std::string_view("+-.0123456789").find(item.front()) != std::string_view::npos;
  }

  /// Appends a continuation line to the CD or TD record above it. NIKON-HIKARI_201911.AGF
  /// (nikon.com) wraps both records this way; the format description does not mention it
  /// (docs/quellen.md, #42).
  void continue_record(const std::vector<std::string_view>& items, std::size_t line_no) {
    std::vector<double>* values = nullptr;
    if (glass_ && continuable_ == Continuable::kCd) {
      values = &glass_->coefficients;
    } else if (glass_ && continuable_ == Continuable::kTd && glass_->thermal) {
      values = &*glass_->thermal;  // set by the TD line above
    }
    if (values == nullptr) {
      fail(line_no,
           "line starts with a number but is no continuation of a CD or TD record (a "
           "continuation line must directly follow CD, TD or another continuation line)");
    }
    for (const std::string_view item : items) values->push_back(number(item, line_no));
    if (continuable_ == Continuable::kCd && values->size() > 10) {
      fail(line_no, "CD has more than 10 coefficients");
    }
    if (continuable_ == Continuable::kTd && values->size() > 7) {
      fail(line_no, "TD has more than 7 values");
    }
  }

  void parse_line(std::string_view line, std::size_t line_no) {
    const std::vector<std::string_view> items = split(line);
    if (stray_) {
      // R4 (#71, NIKON-HIKARI_201911.AGF line 10102): the next line that is not empty must be NM
      // (the manufacturer separates the glass blocks with empty lines).
      if (items.empty()) {
        continuable_ = Continuable::kNone;
        return;
      }
      if (items[0] != "NM") {
        fail(stray_->line, "unknown record '" + stray_->word + "'");
      }
      warn("agf.stray_line", stray_->line,
           "line '" + stray_->word + "' before the next NM record skipped");
      stray_.reset();
    }
    if (items.empty() || items[0].starts_with("!")) {  // empty or comment line
      continuable_ = Continuable::kNone;
      return;
    }
    if (starts_with_number(items[0])) {
      continue_record(items, line_no);
      return;
    }
    continuable_ = Continuable::kNone;
    const std::string_view m = items[0];
    if (!is_record(m)) {
      if (!seen_record_) {  // R3 (#71): text before the first CC or NM, e.g. a header line
        warn("agf.preamble_skipped", line_no,
             "text before the first record skipped: '" + std::string(line) + "'");
        return;
      }
      if (items.size() == 1 && is_word(m)) {  // R4: decided by the next line
        stray_ = Stray{line_no, std::string(m)};
        return;
      }
      fail(line_no, "unknown record '" + std::string(m) + "'");
    }
    if (m == "CC" || m == "NM") seen_record_ = true;
    if (m == "CC") {
      catalog_.comment = text_after(line, "CC");
    } else if (m == "NM") {
      finish_glass();
      if (items.size() < 6) {
        fail(line_no, "NM needs <name> <formula #> <MIL#> <N(d)> <V(d)>");
      }
      AgfGlass g;
      g.name = std::string(items[1]);
      const double formula = number(items[2], line_no);
      if (formula != std::floor(formula) || formula < 0.0 || formula > 1000.0) {
        fail(line_no, "formula number '" + std::string(items[2]) + "' is not an integer");
      }
      g.formula = static_cast<int>(formula);
      g.nd = number(items[4], line_no);
      g.vd = number(items[5], line_no);
      // Optional extras (format description). The melt frequency is kept as written: the
      // description says 1..5, the manufacturers' files also write -1, 0 and 9
      // (docs/quellen.md, #85); it is not used by Raytatouille.
      if (items.size() > 9) {
        fail(line_no, "NM has more than <exclude sub> <status> <melt freq> after V(d)");
      }
      if (items.size() > 6) {
        g.exclude_substitution = nm_extra(items[6], "exclude sub", std::pair{0, 1}, line_no);
      }
      if (items.size() > 7) g.status = nm_extra(items[7], "status", std::pair{0, 4}, line_no);
      if (items.size() > 8) g.melt_frequency = nm_extra(items[8], "melt freq", {}, line_no);
      g.line = line_no;
      glass_ = std::move(g);
    } else if (m == "CD") {
      AgfGlass& g = current(m, line_no);
      if (has_cd_) fail(line_no, "second CD record for glass " + g.name);
      g.coefficients = numbers(items, line_no);
      if (g.coefficients.empty()) fail(line_no, "CD record without coefficients");
      if (g.coefficients.size() > 10) fail(line_no, "CD has more than 10 coefficients");
      has_cd_ = true;
      continuable_ = Continuable::kCd;
    } else if (m == "LD") {
      AgfGlass& g = current(m, line_no);
      if (g.range) fail(line_no, "second LD record for glass " + g.name);
      const std::vector<double> v = numbers(items, line_no);
      if (v.size() != 2 || !(v[0] > 0.0) || !(v[0] < v[1])) {
        fail(line_no, "LD needs <min lambda> < <max lambda> in um, both > 0");
      }
      g.range = WavelengthRange{v[0], v[1]};
    } else if (m == "TD") {
      AgfGlass& g = current(m, line_no);
      if (g.thermal) fail(line_no, "second TD record for glass " + g.name);
      g.thermal = numbers(items, line_no);
      continuable_ = Continuable::kTd;
    } else if (m == "ED") {
      AgfGlass& g = current(m, line_no);
      if (g.extra) fail(line_no, "second ED record for glass " + g.name);
      g.extra = numbers(items, line_no);
    } else if (m == "GC") {
      AgfGlass& g = current(m, line_no);
      if (has_gc_) fail(line_no, "second GC record for glass " + g.name);
      g.comment = text_after(line, "GC");
      has_gc_ = true;
    } else if (m == "MD") {
      AgfGlass& g = current(m, line_no);
      if (g.mechanical) fail(line_no, "second MD record for glass " + g.name);
      g.mechanical = values_with_placeholders(items, 5, line_no);
    } else if (m == "OD") {
      AgfGlass& g = current(m, line_no);
      if (g.other) fail(line_no, "second OD record for glass " + g.name);
      if (items.size() != 7) {
        fail(line_no, "OD needs exactly 6 values (\"_\" for a missing one), got " +
                          std::to_string(items.size() - 1));
      }
      AgfOtherData od;
      // "_" and "-" (SCHOTT, old glasses): not available (docs/quellen.md, #85).
      if (items[1] != "_" && items[1] != "-") od.relative_cost = number(items[1], line_no);
      od.cr = class_range(items[2], line_no);
      od.fr = class_range(items[3], line_no);
      od.sr = class_range(items[4], line_no);
      od.ar = class_range(items[5], line_no);
      od.pr = class_range(items[6], line_no);
      g.other = od;
    } else if (m == "IT") {
      AgfGlass& g = current(m, line_no);
      const std::vector<double> v = numbers(items, line_no);
      if (v.size() != 3) {
        fail(line_no, "IT needs exactly <lambda> <transmittance> <thickness>, got " +
                          std::to_string(v.size()) + " values");
      }
      g.transmission.push_back(AgfTransmission{v[0], v[1], v[2]});
    } else if (m == "BD") {
      (void)current(m, line_no);  // described in the format, not used by Raytatouille
    } else {
      // Every name of is_record() needs a branch above; this keeps a forgotten one an error.
      fail(line_no, "unknown record '" + std::string(m) + "'");
    }
  }

  void finish_glass() {
    if (!glass_) return;
    if (!has_cd_) fail(glass_->line, "glass " + glass_->name + " has no CD record");
    catalog_.glasses.push_back(std::move(*glass_));
    glass_.reset();
    has_cd_ = false;
    has_gc_ = false;
  }

  AgfCatalog catalog_;
  std::optional<AgfGlass> glass_;
  bool has_cd_ = false;
  bool has_gc_ = false;
  Continuable continuable_ = Continuable::kNone;
  bool seen_record_ = false;  ///< a CC or NM record was read (R3)
  struct Stray {
    std::size_t line = 0;
    std::string word;
  };
  std::optional<Stray> stray_;  ///< R4 candidate, decided by the next line
};

std::string upper(std::string text) {
  for (char& c : text) {
    if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
  }
  return text;
}

}  // namespace

AgfError::AgfError(std::string file, std::size_t line, const std::string& message)
    : std::runtime_error(file + (line > 0 ? ":" + std::to_string(line) : std::string()) + ": " +
                         message),
      file_(std::move(file)),
      line_(line) {}

std::string decode_agf_text(std::string_view bytes, const std::string& file) {
  const auto at = [&](std::size_t i) { return static_cast<unsigned char>(bytes[i]); };
  if (bytes.size() >= 2 && at(0) == 0xFE && at(1) == 0xFF) {
    throw AgfError(file, 0, "UTF-16 big endian is not supported (use UTF-16LE or ANSI/UTF-8)");
  }
  if (bytes.size() >= 3 && at(0) == 0xEF && at(1) == 0xBB && at(2) == 0xBF) {
    return std::string(bytes.substr(3));
  }
  if (!(bytes.size() >= 2 && at(0) == 0xFF && at(1) == 0xFE)) return std::string(bytes);

  // UTF-16LE with byte order mark.
  if (bytes.size() % 2 != 0) throw AgfError(file, 0, "UTF-16 file with an odd number of bytes");
  std::string out;
  out.reserve(bytes.size() / 2);
  for (std::size_t i = 2; i < bytes.size(); i += 2) {
    const auto unit = static_cast<char32_t>(at(i) | (at(i + 1) << 8));
    if (unit >= 0xD800 && unit <= 0xDBFF) {  // high surrogate, needs a low one
      if (i + 3 >= bytes.size()) throw AgfError(file, 0, "truncated UTF-16 surrogate pair");
      const auto low = static_cast<char32_t>(at(i + 2) | (at(i + 3) << 8));
      if (low < 0xDC00 || low > 0xDFFF) throw AgfError(file, 0, "invalid UTF-16 surrogate pair");
      append_utf8(out, 0x10000 + ((unit - 0xD800) << 10) + (low - 0xDC00));
      i += 2;
    } else if (unit >= 0xDC00 && unit <= 0xDFFF) {
      throw AgfError(file, 0, "unpaired UTF-16 low surrogate");
    } else {
      append_utf8(out, unit);
    }
  }
  return out;
}

AgfCatalog parse_agf(std::string_view text, std::string name, std::string file) {
  return Parser(std::move(name), std::move(file)).run(text);
}

AgfCatalog load_agf(const std::filesystem::path& file) {
  const std::string display = file.filename().string();
  std::ifstream in(file, std::ios::binary | std::ios::ate);
  if (!in) throw AgfError(display, 0, "cannot open " + file.string());
  // Size, then one read (istreambuf_iterator trips GCC 13's -Wnull-dereference in Release).
  const std::streamoff size = in.tellg();
  if (size < 0) throw AgfError(display, 0, "cannot read " + file.string());
  std::string bytes(static_cast<std::size_t>(size), '\0');
  in.seekg(0);
  if (!in.read(bytes.data(), static_cast<std::streamsize>(size))) {
    throw AgfError(display, 0, "cannot read " + file.string());
  }
  return parse_agf(decode_agf_text(bytes, display), upper(file.stem().string()), display);
}

DispersionFormula agf_formula(const AgfGlass& glass, const std::string& where) {
  // Up to 10 coefficients; missing trailing ones count as 0.
  const auto c = [&](std::size_t i) {
    return i < glass.coefficients.size() ? glass.coefficients[i] : 0.0;
  };
  // CD values beyond those of the formula must be 0: they are 0 in all verifying catalogues, and
  // a line without mnemonic after CD continues the CD record (#42), so a value there most likely
  // is a lost record (e.g. "LD").
  const auto require_unused_zero = [&](std::size_t used, const char* formula_name) {
    for (std::size_t i = used; i < glass.coefficients.size(); ++i) {
      if (glass.coefficients[i] != 0.0) {
        std::ostringstream text;
        text << "glass " << glass.name << " (" << where << "): CD value " << i + 1 << " = "
             << glass.coefficients[i] << " is not used by AGF dispersion formula " << glass.formula
             << " (" << formula_name << ", " << used
             << " values) and must be 0; a line without mnemonic after CD continues the CD "
                "record, see #42";
        throw std::invalid_argument(text.str());
      }
    }
  };
  switch (glass.formula) {
    case 1:  // Schott: CD = a0 a1 a2 a3 a4 a5 (order verified, docs/quellen.md, #24)
      require_unused_zero(6, "Schott");
      return SchottCoefficients{{c(0), c(1), c(2), c(3), c(4), c(5)}};
    case 2:  // Sellmeier 1: CD = K1 L1 K2 L2 K3 L3 (order verified, docs/quellen.md, #24)
      require_unused_zero(6, "Sellmeier 1");
      return Sellmeier1Coefficients{{c(0), c(2), c(4)}, {c(1), c(3), c(5)}};
    case 6:  // Sellmeier 3: CD = K1 L1 .. K4 L4 (order verified, docs/quellen.md, #42)
      require_unused_zero(8, "Sellmeier 3");
      return Sellmeier3Coefficients{{c(0), c(2), c(4), c(6)}, {c(1), c(3), c(5), c(7)}};
    case 12:  // Extended 2: CD = a0..a6 verified (docs/quellen.md, #42); a7 is 0 in every
              // catalogue glass, so its position is not verified and a7 != 0 is rejected.
      if (c(7) != 0.0) {
        std::ostringstream text;
        text << "glass " << glass.name << " (" << where
             << "): AGF dispersion formula 12 (Extended 2) with a7 = " << c(7)
             << " != 0 is not supported: the CD position of a7 is not verified, see #42";
        throw std::invalid_argument(text.str());
      }
      require_unused_zero(8, "Extended 2");
      return Extended2Coefficients{{c(0), c(1), c(2), c(3), c(4), c(5), c(6), 0.0}};
    case 13:  // Extended 3: CD = a0..a8 (order verified, docs/quellen.md, #42)
      require_unused_zero(9, "Extended 3");
      return Extended3Coefficients{{c(0), c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8)}};
    default:
      break;
  }
  static constexpr std::array<const char*, 14> kNames{"unknown",
                                                      "Schott",
                                                      "Sellmeier 1",
                                                      "Herzberger",
                                                      "Sellmeier 2",
                                                      "Conrady",
                                                      "Sellmeier 3",
                                                      "Handbook of Optics 1",
                                                      "Handbook of Optics 2",
                                                      "Sellmeier 4",
                                                      "Extended",
                                                      "Sellmeier 5",
                                                      "Extended 2",
                                                      "Extended 3"};
  const std::size_t index = glass.formula > 0 && glass.formula < static_cast<int>(kNames.size())
                                ? static_cast<std::size_t>(glass.formula)
                                : 0;
  std::ostringstream text;
  text << "glass " << glass.name << " (" << where << "): AGF dispersion formula " << glass.formula
       << " (" << kNames[index]
       << ") is not supported: its coefficient order is not verified; supported are 1 (Schott), 2 "
          "(Sellmeier 1), 6 (Sellmeier 3), 12 (Extended 2 with a7 = 0) and 13 (Extended 3), "
          "see #42";
  throw std::invalid_argument(text.str());
}

CatalogMaterial::CatalogMaterial(AgfGlass glass, const std::string& catalog)
    : glass_(std::move(glass)),
      formula_(
          agf_formula(glass_, "catalog " + catalog + ", line " + std::to_string(glass_.line))) {
  if (glass_.range) {
    const WavelengthRange& r = *glass_.range;
    if (!std::isfinite(r.min_um) || !std::isfinite(r.max_um) || !(r.min_um > 0.0) ||
        !(r.min_um < r.max_um)) {
      throw std::invalid_argument("CatalogMaterial: glass " + glass_.name + " of catalog " +
                                  catalog + ": the LD range must be finite with 0 < min < max");
    }
  }
  if (glass_.thermal) {
    const std::vector<double>& td = *glass_.thermal;
    if (td.size() != 7) {
      throw std::invalid_argument("CatalogMaterial: glass " + glass_.name + " of catalog " +
                                  catalog + ": TD needs D0 D1 D2 E0 E1 Ltk Temp (7 values)");
    }
    thermal_ = SchottThermalCoefficients{td[0], td[1], td[2], td[3], td[4], td[5], td[6]};
  }
}

math::Complex CatalogMaterial::index(double wavelength_um,
                                     double temperature_c,
                                     double /*pressure_atm*/) const {
  // Conversion decided for #25 (class comment): relative wavelength and index at T_ref, 1 atm.
  const double t_ref = thermal_ ? thermal_->reference_temperature_c : kSchottReferenceTemperatureC;
  const double n_air = ciddor_air_index(wavelength_um, t_ref, 1.0);
  const double n_rel = refractive_index(formula_, wavelength_um / n_air);
  double n_abs = n_rel * n_air;  // TIE-19 Eq. (5) at T_ref and 1 atm
  if (thermal_) {
    // TIE-19 Eq. (3) with n_rel for n(lambda, T0) and (4); lambda is the vacuum wavelength.
    n_abs += schott_delta_n_abs(n_rel, wavelength_um, temperature_c - t_ref, *thermal_);
  }
  return {n_abs, 0.0};
}

}  // namespace rtt::material
