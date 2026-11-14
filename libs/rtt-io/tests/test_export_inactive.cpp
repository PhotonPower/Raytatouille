// The ZMX/CSV export with inactive nodes (#380, ADR 0049, point 6): refused until the export
// handles them, so that it never writes wrong thicknesses.

#include <catch2/catch_test_macros.hpp>
#include <stdexcept>
#include <string>

#include "rtt/io/surface_table.hpp"
#include "rtt/io/zmx.hpp"
#include "rtt/material/material.hpp"
#include "rtt/model/system.hpp"
#include "test_support.hpp"

TEST_CASE("the export refuses inactive nodes until it handles them (ADR 0049, point 6)",
          "[io][0.10][export]") {
  rtt::model::System s = rtt::model::test::make_singlet();
  rtt::material::MaterialLibrary lib;
  lib.add_catalog(std::string(RTT_CATALOG_DIR) + "/schott.agf");
  REQUIRE_NOTHROW(rtt::io::export_zmx(s, lib));  // content guard: the active system exports
  rtt::model::test::element(s, 1).active = false;
  CHECK_THROWS_AS(rtt::io::export_zmx(s, lib), std::invalid_argument);
  CHECK_THROWS_AS(rtt::io::write_surface_table(s, lib), std::invalid_argument);
}
