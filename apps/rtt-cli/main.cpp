// rtt - command line tool for Raytatouille system files (ADR 0044, 0045, 0052).

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "cli.hpp"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

int main(int argc, char** argv) {
#ifdef _WIN32
  // The data on stdout are the bytes of the result (ADR 0052, point 2): no \n -> \r\n.
  static_cast<void>(_setmode(_fileno(stdout), _O_BINARY));
#endif
  const std::vector<std::string> args(argv + 1, argv + argc);
  rtt::cli::install_interrupt_handler();  // Ctrl-C cancels trace, analyze and optimize
  return rtt::cli::run(args, std::cout, std::cerr, rtt::cli::interrupt_token());
}
