// rtt - command line tool for Raytatouille system files.

#include <cstddef>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "rtt/io/json_io.hpp"
#include "rtt/io/zmx.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/validate.hpp"

namespace {

constexpr int kOk = 0;
constexpr int kInvalid = 1;
constexpr int kUsage = 2;

void usage(std::ostream& os) {
  os << "usage:\n"
        "  rtt --version\n"
        "  rtt validate <file.rtt.json>...      check structure and semantics\n"
        "  rtt format <file.rtt.json>...        rewrite files in canonical form\n"
        "  rtt format --check <file.rtt.json>...  fail if a file is not canonical\n"
        "  rtt import <file.zmx> <out.rtt.json> --catalog <file.agf>...\n"
        "                                       import the sequential part of a ZMX file\n";
}

std::string read_text(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open '" + path + "'");
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

/// Warnings of the reader (e.g. io.pickup_dropped when a file before schema 0.4 is migrated) go
/// to stderr, so that nothing a file contained disappears silently (#162).
void print_read_warnings(const std::string& file, const std::vector<rtt::model::Diagnostic>& w) {
  for (const auto& d : w) std::cerr << file << ": " << rtt::model::to_string(d) << "\n";
}

int validate(const std::vector<std::string>& files) {
  int result = kOk;
  for (const std::string& file : files) {
    try {
      std::vector<rtt::model::Diagnostic> read_warnings;
      const rtt::model::System system = rtt::io::load_system(file, read_warnings);
      print_read_warnings(file, read_warnings);
      const auto diagnostics = rtt::model::validate(system);
      for (const auto& d : diagnostics)
        std::cout << file << ": " << rtt::model::to_string(d) << "\n";
      if (rtt::model::has_errors(diagnostics)) {
        result = kInvalid;
      } else {
        std::cout << file << ": ok\n";
      }
    } catch (const std::exception& e) {
      std::cout << file << ": error " << e.what() << "\n";
      result = kInvalid;
    }
  }
  return result;
}

int format(const std::vector<std::string>& files, bool check_only) {
  int result = kOk;
  for (const std::string& file : files) {
    try {
      const std::string original = read_text(file);
      std::vector<rtt::model::Diagnostic> read_warnings;
      const std::string canonical =
          rtt::io::to_json(rtt::io::parse_system(original, read_warnings));
      print_read_warnings(file, read_warnings);
      if (canonical == original) continue;
      if (check_only) {
        std::cout << file << ": not in canonical form (run 'rtt format')\n";
        result = kInvalid;
      } else {
        std::ofstream out(file, std::ios::binary | std::ios::trunc);
        out << canonical;
        std::cout << file << ": formatted\n";
      }
    } catch (const std::exception& e) {
      std::cout << file << ": error " << e.what() << "\n";
      result = kInvalid;
    }
  }
  return result;
}

/// rtt import (ADR 0044): reads a ZMX file, resolves its glasses in the given AGF catalogues and
/// writes the system in canonical form. Warnings go to stderr; exit 1 if the file cannot be
/// imported (also a missing file), 0 otherwise.
int import_file(const std::string& zmx,
                const std::string& out,
                const std::vector<std::string>& catalogs) {
  try {
    rtt::material::MaterialLibrary lib;
    for (const std::string& c : catalogs) lib.add_catalog(c);
    const rtt::io::ZmxImport r = rtt::io::load_zmx(zmx, lib);
    for (const auto& d : r.warnings) std::cerr << zmx << ": " << rtt::model::to_string(d) << "\n";
    rtt::io::save_system(r.system, out);
    std::cout << zmx << ": imported to " << out << "\n";
    return kOk;
  } catch (const std::exception& e) {
    std::cout << zmx << ": error " << e.what() << "\n";
    return kInvalid;
  }
}

}  // namespace

int main(int argc, char** argv) {
  const std::vector<std::string> args(argv + 1, argv + argc);
  if (args.empty()) {
    usage(std::cerr);
    return kUsage;
  }
  const std::string& command = args[0];
  if (command == "--version") {
    std::cout << "rtt " << RTT_VERSION << " (schema " << rtt::model::kSchemaVersion << ")\n";
    return kOk;
  }
  if (command == "--help" || command == "-h") {
    usage(std::cout);
    return kOk;
  }
  std::vector<std::string> rest(args.begin() + 1, args.end());
  if (command == "validate" && !rest.empty()) return validate(rest);
  if (command == "import" && rest.size() >= 2) {
    std::vector<std::string> catalogs;
    bool ok = true;
    for (std::size_t i = 2; i < rest.size(); ++i) {
      if (rest[i] == "--catalog" && i + 1 < rest.size()) {
        catalogs.push_back(rest[++i]);
      } else {
        ok = false;
      }
    }
    if (ok) return import_file(rest[0], rest[1], catalogs);
  }
  if (command == "format") {
    const bool check = !rest.empty() && rest[0] == "--check";
    if (check) rest.erase(rest.begin());
    if (!rest.empty()) return format(rest, check);
  }
  usage(std::cerr);
  return kUsage;
}
