#pragma once

/// @file agf.hpp
/// Glass catalogues in the ANSI Glass Format (AGF) of Zemax OpticStudio.
///
/// Format (docs/quellen.md): Ansys Zemax OpticStudio User Guide, Release 2025 R1, "The AGF & BGF
/// File Formats". A header line "CC <comment>", optional lines starting with '!', then one
/// record per glass. Each line is a two-letter mnemonic followed by data items separated by
/// spaces:
///   NM <name> <formula #> <MIL#> <N(d)> <V(d)> [<exclude sub> <status> <melt freq>]
///   CD <dispersion coefficients 1 - 10>
///   LD <min lambda> <max lambda>                  (um)
///   TD <D0> <D1> <D2> <E0> <E1> <Ltk> <Temp>      (dn/dT model, reference temperature in degC)
///   ED <TCE -30..70> <TCE 100..300> <density> <dPgF> <ignore thermal expansion>
///   GC, MD, OD, IT, BD                            (read and ignored)
/// The order of the CD coefficients is not part of that description; it is verified for the
/// supported formulas against the free SCHOTT catalogue (docs/quellen.md, #24): formula 1
/// (Schott) a0..a5, formula 2 (Sellmeier 1) K1 L1 K2 L2 K3 L3. All other formula numbers are
/// rejected when the glass is resolved (#42).

#include <cstddef>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "rtt/material/dispersion.hpp"
#include "rtt/material/material.hpp"

namespace rtt::material {

/// Error in an AGF file. what() names the file and the line.
class AgfError : public std::runtime_error {
 public:
  AgfError(std::string file, std::size_t line, const std::string& message);

  [[nodiscard]] const std::string& file() const noexcept { return file_; }
  /// 1-based line number; 0 if the error concerns the whole file.
  [[nodiscard]] std::size_t line() const noexcept { return line_; }

 private:
  std::string file_;
  std::size_t line_;
};

/// One glass of an AGF catalogue as read from the file (no unit conversion).
struct AgfGlass {
  std::string name;
  int formula = 0;                             ///< AGF dispersion formula number
  double nd = 0.0;                             ///< N(d) from NM, for reference only
  double vd = 0.0;                             ///< V(d) from NM, for reference only
  std::vector<double> coefficients;            ///< CD, as in the file
  std::optional<WavelengthRange> range;        ///< LD in um
  std::optional<std::vector<double>> thermal;  ///< TD: D0 D1 D2 E0 E1 Ltk Temp (unused, #25)
  std::optional<std::vector<double>> extra;    ///< ED (unused)
  std::size_t line = 0;                        ///< line of the NM record
};

/// A parsed AGF catalogue.
struct AgfCatalog {
  std::string name;     ///< catalogue name (see MaterialLibrary::add_catalog)
  std::string file;     ///< file name for messages
  std::string comment;  ///< CC header text
  std::vector<AgfGlass> glasses;
};

/// Converts raw file bytes to UTF-8 text: UTF-16LE with byte order mark is converted, a UTF-8
/// byte order mark is removed, anything else (ANSI, UTF-8) is passed on unchanged. ANSI bytes
/// above 0x7F (code page characters in comments) are therefore not valid UTF-8 afterwards; the
/// parser only interprets ASCII (mnemonics, names, numbers), so this affects comments only.
/// @param file file name for error messages
/// @throws AgfError for UTF-16BE or malformed UTF-16
[[nodiscard]] std::string decode_agf_text(std::string_view bytes, const std::string& file);

/// Parses AGF text (UTF-8, LF or CRLF line ends).
/// @param text text of the file, see decode_agf_text
/// @param name catalogue name stored in the result
/// @param file file name for error messages
/// @throws AgfError with file and line for malformed or unknown records, records before the
///         first NM, missing or empty CD, duplicate CD/LD/TD/ED records of a glass, duplicate
///         glass names and invalid numbers
[[nodiscard]] AgfCatalog parse_agf(std::string_view text, std::string name, std::string file);

/// Reads, decodes and parses an AGF file; the catalogue name is the file name without
/// extension in upper case (#24).
/// @throws AgfError as parse_agf, or if the file cannot be read
[[nodiscard]] AgfCatalog load_agf(const std::filesystem::path& file);

/// Non-absorbing material of a catalogue glass with a supported formula. The index is the
/// formula value, absolute until the conversion relative to air (#25); temperature and pressure
/// are ignored until then. The AGF data (including TD and ED) stay available via glass().
class CatalogMaterial final : public Material {
 public:
  /// @param glass   glass data, e.g. from parse_agf
  /// @param catalog catalogue name for messages
  /// @throws std::invalid_argument if the formula number is not supported (see agf_formula) or
  ///         the LD range is not finite with 0 < min < max
  CatalogMaterial(AgfGlass glass, const std::string& catalog);

  /// n(lambda) + 0i of the glass formula; see Material::index.
  [[nodiscard]] math::Complex index(double wavelength_um,
                                    double temperature_c,
                                    double pressure_atm) const override;

  /// LD range of the glass; none if the catalogue has no LD record for it.
  [[nodiscard]] std::optional<WavelengthRange> wavelength_range_um() const override {
    return glass_.range;
  }

  /// The glass data as read from the catalogue.
  [[nodiscard]] const AgfGlass& glass() const noexcept { return glass_; }
  /// The dispersion formula built from the CD record.
  [[nodiscard]] const DispersionFormula& formula() const noexcept { return formula_; }

 private:
  AgfGlass glass_;
  DispersionFormula formula_;
};

/// Builds the dispersion formula of a glass. Supported: formula 1 (Schott, CD = a0..a5) and
/// formula 2 (Sellmeier 1, CD = K1 L1 K2 L2 K3 L3). Missing trailing coefficients count as 0
/// (the format allows "up to 10").
/// @param where text naming catalogue, file and line for messages
/// @throws std::invalid_argument for other formula numbers (see #42)
[[nodiscard]] DispersionFormula agf_formula(const AgfGlass& glass, const std::string& where);

}  // namespace rtt::material
