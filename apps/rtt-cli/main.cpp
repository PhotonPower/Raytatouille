// rtt - command line tool for Raytatouille system files.

#include <cctype>
#include <charconv>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include "rtt/io/json_io.hpp"
#include "rtt/io/surface_table.hpp"
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
        "                                       import the sequential part of a ZMX file\n"
        "  rtt import <file.csv> <out.rtt.json> --base <system.rtt.json> --catalog <file.agf>...\n"
        "                                       read a surface table into a base system\n"
        "  rtt export <file.rtt.json> <out.zmx|out.csv> --catalog <file.agf>... [--path <name>]\n"
        "             [--config <index>]        export one path of a coaxial system as ZMX\n"
        "                                       or as a surface table\n";
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

/// The extension of `file` in lower case, with the dot ("" if it has none).
std::string extension_of(const std::string& file) {
  const std::size_t dot = file.find_last_of('.');
  std::string e = dot == std::string::npos ? std::string() : file.substr(dot);
  for (char& c : e) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return e;
}

/// rtt import (ADR 0044; ADR 0045, point 5): reads a ZMX file, or a surface table into the base
/// system `base`, resolves the glasses in the given AGF catalogues and writes the system in
/// canonical form. Warnings go to stderr; exit 1 if the file cannot be imported (also a missing
/// file), 0 otherwise.
int import_file(const std::string& in,
                const std::string& out,
                const std::vector<std::string>& catalogs,
                const std::optional<std::string>& base) {
  try {
    rtt::material::MaterialLibrary lib;
    for (const std::string& c : catalogs) lib.add_catalog(c);
    rtt::model::System system;
    std::vector<rtt::model::Diagnostic> warnings;
    if (base) {
      std::vector<rtt::model::Diagnostic> read_warnings;
      rtt::model::System b = rtt::io::load_system(*base, read_warnings);
      print_read_warnings(*base, read_warnings);
      rtt::io::TableImport r = rtt::io::read_surface_table(read_text(in), std::move(b), lib);
      system = std::move(r.system);
      warnings = std::move(r.warnings);
    } else {
      rtt::io::ZmxImport r = rtt::io::load_zmx(in, lib);
      system = std::move(r.system);
      warnings = std::move(r.warnings);
    }
    for (const auto& d : warnings) std::cerr << in << ": " << rtt::model::to_string(d) << "\n";
    rtt::io::save_system(system, out);
    std::cout << in << ": imported to " << out << "\n";
    return kOk;
  } catch (const std::exception& e) {
    std::cout << in << ": error " << e.what() << "\n";
    return kInvalid;
  }
}

/// rtt export (ADR 0045, point 5): writes one path of a coaxial system as a ZMX file or, for an
/// output ending in .csv, as a surface table. Warnings go to stderr; exit 1 if the system cannot
/// be read or exported, 0 otherwise.
int export_file(const std::string& file,
                const std::string& out,
                const std::vector<std::string>& catalogs,
                const rtt::io::ExportOptions& options) {
  try {
    std::vector<rtt::model::Diagnostic> read_warnings;
    const rtt::model::System system = rtt::io::load_system(file, read_warnings);
    print_read_warnings(file, read_warnings);
    rtt::material::MaterialLibrary lib;
    for (const std::string& c : catalogs) lib.add_catalog(c);
    std::vector<rtt::model::Diagnostic> warnings;
    if (extension_of(out) == ".csv") {
      rtt::io::TableExport table = rtt::io::write_surface_table(system, lib, options);
      std::ofstream f(out, std::ios::binary | std::ios::trunc);
      f << table.text;
      if (!f) throw std::runtime_error("cannot write " + out);
      warnings = std::move(table.warnings);
    } else {
      warnings = rtt::io::save_zmx(out, system, lib, options);
    }
    for (const auto& d : warnings) std::cerr << file << ": " << rtt::model::to_string(d) << "\n";
    std::cout << file << ": exported to " << out << "\n";
    return kOk;
  } catch (const std::exception& e) {
    std::cout << file << ": error " << e.what() << "\n";
    return kInvalid;
  }
}

/// rtt export with the arguments after the command: the exit code, none for a wrong call.
std::optional<int> run_export(const std::vector<std::string>& rest) {
  if (rest.size() < 2 || (extension_of(rest[1]) != ".zmx" && extension_of(rest[1]) != ".csv")) {
    return std::nullopt;
  }
  std::vector<std::string> catalogs;
  rtt::io::ExportOptions options;
  for (std::size_t i = 2; i < rest.size(); i += 2) {
    if (i + 1 >= rest.size()) return std::nullopt;
    const std::string& value = rest[i + 1];
    if (rest[i] == "--catalog") {
      catalogs.push_back(value);
    } else if (rest[i] == "--path") {
      options.path = value;
    } else if (rest[i] == "--config") {
      std::size_t index = 0;
      const auto [end, ec] = std::from_chars(value.data(), value.data() + value.size(), index);
      if (ec != std::errc() || end != value.data() + value.size()) return std::nullopt;
      options.configuration = index;
    } else {
      return std::nullopt;
    }
  }
  return export_file(rest[0], rest[1], catalogs, options);
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
    std::optional<std::string> base;
    bool ok = true;
    for (std::size_t i = 2; i < rest.size(); ++i) {
      if (rest[i] == "--catalog" && i + 1 < rest.size()) {
        catalogs.push_back(rest[++i]);
      } else if (rest[i] == "--base" && i + 1 < rest.size() && !base) {
        base = rest[++i];
      } else {
        ok = false;
      }
    }
    // a table needs a base system, a ZMX file takes none
    if (ok && (extension_of(rest[0]) == ".csv") == base.has_value()) {
      return import_file(rest[0], rest[1], catalogs, base);
    }
  }
  if (command == "export") {
    if (const std::optional<int> code = run_export(rest)) return *code;
  }
  if (command == "format") {
    const bool check = !rest.empty() && rest[0] == "--check";
    if (check) rest.erase(rest.begin());
    if (!rest.empty()) return format(rest, check);
  }
  usage(std::cerr);
  return kUsage;
}
