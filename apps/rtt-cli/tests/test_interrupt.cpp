// Ctrl-C at the real program (ADR 0052, points 6 and 8), POSIX only: rtt optimize runs as a child
// with --ftol 0 --xtol 0 (several seconds and some 26000 iterations until CONVERGED_STEP, so a
// missing cancellation shows as exit 0) and with --progress; after the first progress line on
// stderr (so the handler is installed and the run is under way) the test sends SIGINT. stderr is
// drained until EOF before waiting, so the child never blocks on a full pipe; after 60 s the child
// is killed and the test fails.
// rtt must end with 130, not write the system and write --result with the status CANCELLED.
// Under Windows the handler function is tested directly (test_cli.cpp).

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#ifndef _WIN32
#include <poll.h>
#include <sys/wait.h>
#include <unistd.h>

#include <csignal>
#include <ctime>
#endif

namespace fs = std::filesystem;

TEST_CASE("SIGINT cancels rtt optimize with exit 130 (ADR 0052, point 6)", "[cli][interrupt]") {
#ifdef _WIN32
  SUCCEED("POSIX only; the handler is tested in test_cli.cpp");
#else
  const fs::path dir = fs::temp_directory_path() / "rtt_cli_interrupt";
  fs::create_directories(dir);
  const fs::path out = dir / "long.rtt.json";
  const fs::path result = dir / "long.result.json";
  fs::remove(out);
  fs::remove(result);
  const std::string input = std::string(RTT_REFERENCE_DIR) + "/m5/singlet_optim.rtt.json";
  std::vector<std::string> args = {RTT_CLI_EXE, "optimize",     input, "--ftol",
                                   "0",         "--xtol",       "0",   "--max-iterations",
                                   "100000000", "--progress",   "-o",  out.string(),
                                   "--result",  result.string()};
  std::vector<char*> argv;
  argv.reserve(args.size() + 1);
  for (std::string& a : args) argv.push_back(a.data());
  argv.push_back(nullptr);
  std::array<int, 2> pipe_fds{};
  REQUIRE(pipe(pipe_fds.data()) == 0);
  // fork and execv (not posix_spawn: its header does not compile with GCC under macOS).
  const pid_t pid = fork();
  REQUIRE(pid >= 0);
  if (pid == 0) {
    dup2(pipe_fds[1], STDERR_FILENO);
    close(pipe_fds[0]);
    close(pipe_fds[1]);
    execv(RTT_CLI_EXE, argv.data());
    _exit(127);  // execv failed
  }
  close(pipe_fds[1]);
  // Wait (at most 60 s) for the first progress line.
  std::string err;
  bool ready = false;
  const std::time_t start = std::time(nullptr);
  while (!ready && std::time(nullptr) - start < 60) {
    pollfd p{pipe_fds[0], POLLIN, 0};
    if (poll(&p, 1, 1000) <= 0) continue;
    std::array<char, 256> buffer{};
    const ssize_t n = read(pipe_fds[0], buffer.data(), buffer.size());
    if (n <= 0) break;
    err.append(buffer.data(), static_cast<std::size_t>(n));
    ready = err.find("progress ") != std::string::npos;
  }
  INFO(err);
  REQUIRE(ready);
  REQUIRE(kill(pid, SIGINT) == 0);
  // Drain stderr until EOF (the child ended) or 60 s, then wait; kill the child after the limit.
  const std::time_t sent = std::time(nullptr);
  bool eof = false;
  while (!eof && std::time(nullptr) - sent < 60) {
    pollfd p{pipe_fds[0], POLLIN, 0};
    if (poll(&p, 1, 1000) <= 0) continue;
    std::array<char, 4096> buffer{};
    eof = read(pipe_fds[0], buffer.data(), buffer.size()) <= 0;
  }
  close(pipe_fds[0]);
  if (!eof) kill(pid, SIGKILL);
  int status = 0;
  REQUIRE(waitpid(pid, &status, 0) == pid);
  REQUIRE(eof);  // the child ended within 60 s of SIGINT
  REQUIRE(WIFEXITED(status));
  CHECK(WEXITSTATUS(status) == 130);
  CHECK_FALSE(fs::exists(out));
  const std::ifstream in(result, std::ios::binary);
  std::ostringstream text;
  text << in.rdbuf();
  CHECK(text.str().find(R"("status": "CANCELLED")") != std::string::npos);
#endif
}
