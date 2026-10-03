// rtt - command line tool for Raytatouille system files.

#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "rtt/io/json_io.hpp"
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
        "  rtt format --check <file.rtt.json>...  fail if a file is not canonical\n";
}

std::string read_text(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open '" + path + "'");
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

int validate(const std::vector<std::string>& files) {
  int result = kOk;
  for (const std::string& file : files) {
    try {
      const auto diagnostics = rtt::model::validate(rtt::io::load_system(file));
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
      const std::string canonical = rtt::io::to_json(rtt::io::parse_system(original));
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
  if (command == "format") {
    const bool check = !rest.empty() && rest[0] == "--check";
    if (check) rest.erase(rest.begin());
    if (!rest.empty()) return format(rest, check);
  }
  usage(std::cerr);
  return kUsage;
}
