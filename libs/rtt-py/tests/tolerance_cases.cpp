// Tolerance cases of rtt_py_reference (#261): the same calls as CASES in
// test_bitwise_tolerance.py, on the M7 reference achromat with the SCHOTT catalogue.

#include "tolerance_cases.hpp"

#include <oneapi/tbb/task_arena.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "npy_writer.hpp"
#include "rtt/io/json_io.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/system.hpp"
#include "rtt/tolerance/expand.hpp"

namespace rtt::py::reference {
namespace {

namespace fs = std::filesystem;

/// Writes the arrays of one case as <out>/<case>.<name>.npy.
class Writer {
 public:
  Writer(fs::path out, std::string name) : out_(std::move(out)), name_(std::move(name)) {}

  void f8(const std::string& array, const std::vector<double>& v) const {
    write_column<double>(file(array), "<f8", v);
  }
  void u8(const std::string& array, const std::vector<std::uint64_t>& v) const {
    write_column<std::uint64_t>(file(array), "<u8", v);
  }
  /// Lines joined with '\n' as UTF-8 bytes.
  void text(const std::string& array, const std::vector<std::string>& lines) const {
    std::string joined;
    for (std::size_t i = 0; i < lines.size(); ++i) {
      if (i > 0) joined += '\n';
      joined += lines[i];
    }
    write_column<std::uint8_t>(file(array), "|u1",
                               std::vector<std::uint8_t>(joined.begin(), joined.end()));
  }

 private:
  [[nodiscard]] fs::path file(const std::string& array) const {
    return out_ / (name_ + "." + array + ".npy");
  }

  fs::path out_;
  std::string name_;
};

template <typename E>
std::uint64_t code(E e) {
  return static_cast<std::uint64_t>(e);
}

/// 0 without, else the index + 1.
std::uint64_t optional_index(const std::optional<std::size_t>& i) {
  return i ? *i + 1 : 0;
}

// Flattening; the names match flatten_items() in test_bitwise_tolerance.py.

void write_items(const Writer& w, const tolerance::ToleranceItems& r) {
  std::vector<std::string> ids;
  std::vector<double> limit;
  std::vector<std::uint64_t> kind;
  std::vector<std::uint64_t> distribution;
  for (const model::ToleranceItem& i : r.items) {
    ids.push_back(i.id);
    limit.push_back(i.limit);
    kind.push_back(code(i.kind));
    distribution.push_back(code(i.distribution));
  }
  std::vector<std::uint64_t> origin;
  for (const tolerance::ItemOrigin& o : r.origins) {
    origin.push_back(code(o.source));
    origin.push_back(o.index);
  }
  std::vector<std::string> class_ids;
  std::vector<std::uint64_t> class_ints;
  for (const model::ClassTolerance& c : r.classes.items) {
    class_ids.push_back(c.id);
    class_ints.insert(class_ints.end(),
                      {code(c.kind), code(c.level), code(c.tolerance_class), code(c.skip),
                       optional_index(c.value_override), optional_index(c.replaced_by)});
  }
  std::vector<std::uint64_t> overrides;
  for (const model::OverrideStatus& v : r.classes.value_overrides) {
    overrides.insert(overrides.end(), {v.index, code(v.effect), code(v.skip)});
  }
  w.text("ids", ids);
  w.f8("limit", limit);
  w.u8("kind", kind);
  w.u8("distribution", distribution);
  w.u8("origin", origin);
  w.text("class_ids", class_ids);
  w.u8("class_ints", class_ints);
  w.u8("overrides", overrides);
}

}  // namespace

void run_tolerance_cases(const fs::path& reference_dir,
                         const fs::path& catalog_dir,
                         const fs::path& out,
                         int threads) {
  material::MaterialLibrary schott;
  schott.add_catalog(catalog_dir / "schott.agf");
  const model::System achromat = io::load_system(reference_dir / "m7/achromat_tolerances.rtt.json");
  oneapi::tbb::task_arena arena(threads);
  arena.execute(
      [&] { write_items(Writer(out, "achromat_expand"), tolerance::expand(achromat, schott)); });
}

}  // namespace rtt::py::reference
