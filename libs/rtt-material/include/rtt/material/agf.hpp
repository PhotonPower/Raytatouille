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
///   GC <comment>
///   MD <E> <Poisson> <HK> <c_p> <heat conductivity> (exactly 5 values)
///   OD <rel. cost> <CR> <FR> <SR> <AR> <PR>       (exactly 6 values; classes may be "a-b")
///   IT <lambda> <internal transmittance> <thickness> (one line per point, exactly 3 values)
///   BD                                            (read and ignored)
/// Placeholders (#85, ADR 0008 addendum): the manufacturers write "_" for a missing value in
/// MD, OD and the NM extras and "-" in OD and the NM extras; both mean "not available". IT has
/// no placeholders. Units of ED, MD and IT: see AgfGlass.
/// A line whose first character is a digit, '+', '-' or '.' (a number) continues the CD or TD
/// record directly above it (or
/// above its previous continuation line); vendor files wrap long records this way (#42). The
/// format description does not mention it; the evidence is NIKON-HIKARI_201911.AGF
/// (docs/quellen.md).
/// The order of the CD coefficients is not part of that description; it is verified for the
/// supported formulas against N(d) and V(d) of free vendor catalogues (docs/quellen.md): formula
/// 1 (Schott) a0..a5 and formula 2 (Sellmeier 1) K1 L1 K2 L2 K3 L3 (SCHOTT, #24); formula 6
/// (Sellmeier 3) K1 L1 .. K4 L4, formula 12 (Extended 2) a0..a6 and formula 13 (Extended 3)
/// a0..a8 (NIKON-HIKARI, #42). All other formula numbers, and formula 12 with a7 != 0, are
/// rejected when the glass is resolved.

#include <cstddef>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "rtt/material/dispersion.hpp"
#include "rtt/material/material.hpp"
#include "rtt/material/thermal.hpp"

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

/// A resistance class of the OD record: one class (low == high) or a range "a-b" as SCHOTT
/// writes it (e.g. "1-2", docs/quellen.md, #85).
struct AgfClassRange {
  double low = 0.0;
  double high = 0.0;
  bool operator==(const AgfClassRange&) const = default;
};

/// The OD record ("other data"): relative cost and the resistance classes CR (climate), FR
/// (stain), SR (acid), AR (alkali), PR (phosphate). The format description writes -1 for "not
/// available", which is kept as written; "_" and "-" in the file are nullopt.
struct AgfOtherData {
  std::optional<double> relative_cost;
  std::optional<AgfClassRange> cr;
  std::optional<AgfClassRange> fr;
  std::optional<AgfClassRange> sr;
  std::optional<AgfClassRange> ar;
  std::optional<AgfClassRange> pr;
  bool operator==(const AgfOtherData&) const = default;
};

/// One IT record: internal transmittance of a sample of the given thickness.
struct AgfTransmission {
  double wavelength_um = 0.0;  ///< wavelength in um (as LD; checked against SCHOTT, quellen.md)
  double transmittance = 0.0;  ///< internal transmittance tau_i, 0..1
  double thickness_mm = 0.0;   ///< sample thickness in mm
};

/// One glass of an AGF catalogue as read from the file (no unit conversion). Units of the data
/// records follow the format description and, where it gives none, the SCHOTT N-BK7 data sheet
/// (docs/quellen.md, #85). A value written as the placeholder "_" (MD, OD, NM extras) or "-" (NM
/// extras) is "not available" (nullopt), as in the manufacturers' files.
struct AgfGlass {
  std::string name;
  int formula = 0;  ///< AGF dispersion formula number
  double nd = 0.0;  ///< N(d) from NM, for reference only
  double vd = 0.0;  ///< V(d) from NM, for reference only
  /// NM: exclude from glass substitution, 0 (no) or 1 (yes)
  std::optional<int> exclude_substitution;
  /// NM: 0 standard, 1 preferred, 2 obsolete, 3 special, 4 melt
  std::optional<int> status;
  /// NM: relative melt frequency 1..5. The format description gives 1..5; the manufacturers'
  /// files also write -1 and 0 for "not given", which are stored as nullopt.
  std::optional<int> melt_frequency;
  std::string comment;                         ///< GC text (empty without GC)
  std::vector<double> coefficients;            ///< CD, as in the file
  std::optional<WavelengthRange> range;        ///< LD in um
  std::optional<std::vector<double>> thermal;  ///< TD: D0 D1 D2 E0 E1 Ltk Temp (thermal.hpp)
  /// ED: TCE -30..70 degC and TCE 100..300 degC in 1e-6/K, density in g/cm^3, dPgF, ignore
  /// thermal expansion (0/1); as many values as the file gives (manufacturers write 4 to 6).
  /// SCHOTT writes alpha(+20/+300 degC) of its data sheet as the second value, which the format
  /// description calls TCE 100..300 degC ("currently not used").
  std::optional<std::vector<double>> extra;
  /// MD, exactly 5 values: Young's modulus in GPa, Poisson's ratio, Knoop hardness HK in
  /// kgf/mm^2,
  /// specific heat capacity in J/(kg K), heat conductivity in W/(m K)
  std::optional<std::vector<std::optional<double>>> mechanical;
  /// OD, exactly 6 values: relative cost (a number, "_" or "-"), CR, FR, SR, AR, PR (a number,
  /// "_", "-" or a class range "a-b" with a <= b)
  std::optional<AgfOtherData> other;
  std::vector<AgfTransmission> transmission;  ///< IT records in file order
  std::size_t line = 0;                       ///< line of the NM record
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
///         glass names, invalid numbers, a continuation line that does not follow CD or TD
///         (or follows an empty or comment line), more than 10 CD or (with continuation
///         lines) 7 TD values; MD, OD and IT with another number of values than 5, 6 and 3, a
///         second GC, MD or OD record, and NM extras that are not "_"/"-" or in range
///         (exclude sub 0/1, status 0..4, melt freq -1..5) or more than three of them
[[nodiscard]] AgfCatalog parse_agf(std::string_view text, std::string name, std::string file);

/// Reads, decodes and parses an AGF file; the catalogue name is the file name without
/// extension in upper case (#24).
/// @throws AgfError as parse_agf, or if the file cannot be read
[[nodiscard]] AgfCatalog load_agf(const std::filesystem::path& file);

/// Non-absorbing material of a catalogue glass with a supported formula.
///
/// Catalogue data are relative to air at the reference temperature T_ref of the glass and
/// 1 atm, as a function of the wavelength in that air (SCHOTT TIE-29 p. 1; Ansys OpticStudio
/// "Index of Refraction Computation"). index() returns the absolute index at temperature T
/// (decided for #25):
///   n_air = n_air(lambda_vac, T_ref, 1 atm)            (Ciddor, air.hpp)
///   n_rel = formula(lambda_vac / n_air)                 (relative wavelength)
///   n_abs = n_rel * n_air + Delta n_abs(T - T_ref)      (TIE-19 Eq. (5) and (3)/(4))
/// T_ref and the coefficients of Delta n_abs come from the TD record (thermal.hpp); without TD,
/// Delta n_abs = 0 and T_ref = 20 degC. The pressure of the medium is ignored (a solid).
/// The AGF data (including TD and ED) stay available via glass().
class CatalogMaterial final : public Material {
 public:
  /// @param glass   glass data, e.g. from parse_agf
  /// @param catalog catalogue name for messages
  /// @throws std::invalid_argument if the formula number is not supported (see agf_formula),
  ///         the LD range is not finite with 0 < min < max, or a TD record does not have the
  ///         seven values D0 D1 D2 E0 E1 Ltk Temp
  CatalogMaterial(AgfGlass glass, const std::string& catalog);

  /// Absolute index n + 0i at the vacuum wavelength and temperature, see the class comment.
  [[nodiscard]] math::Complex index(double wavelength_um,
                                    double temperature_c,
                                    double pressure_atm) const override;

  /// LD range of the glass; none if the catalogue has no LD record for it. LD is given in the
  /// relative wavelength; it is compared with vacuum wavelengths (difference < 0.03 %).
  [[nodiscard]] std::optional<WavelengthRange> wavelength_range_um() const override {
    return glass_.range;
  }

  /// The glass data as read from the catalogue.
  [[nodiscard]] const AgfGlass& glass() const noexcept { return glass_; }
  /// The dispersion formula built from the CD record (relative index over relative wavelength).
  [[nodiscard]] const DispersionFormula& formula() const noexcept { return formula_; }
  /// Thermal coefficients from the TD record; none without TD.
  [[nodiscard]] const std::optional<SchottThermalCoefficients>& thermal() const noexcept {
    return thermal_;
  }

 private:
  AgfGlass glass_;
  DispersionFormula formula_;
  std::optional<SchottThermalCoefficients> thermal_;
};

/// Builds the dispersion formula of a glass (lambda in um, coefficients as in dispersion.hpp).
/// Supported, with the verified CD order (docs/quellen.md): formula 1 (Schott, CD = a0..a5),
/// 2 (Sellmeier 1, K1 L1 K2 L2 K3 L3), 6 (Sellmeier 3, K1 L1 .. K4 L4), 12 (Extended 2,
/// a0..a6; a7 must be 0) and 13 (Extended 3, a0..a8). Missing trailing coefficients count as 0
/// (the format allows "up to 10"); CD values beyond those of the formula must be 0 (they are 0
/// in all glasses of the verifying catalogues, and a value there most likely is a record whose
/// mnemonic was lost and that continues CD, #42).
/// @param where text naming catalogue, file and line for messages
/// @throws std::invalid_argument for other formula numbers, for formula 12 with a7 != 0,
///         whose CD position is not verified, and for a CD value != 0 beyond those of the
///         formula (see #42)
[[nodiscard]] DispersionFormula agf_formula(const AgfGlass& glass, const std::string& where);

}  // namespace rtt::material
