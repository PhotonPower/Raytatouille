#pragma once

/// @file result_json.hpp
/// Writer of the result format raytatouille-result (ADR 0023) in C++ for the command line (ADR
/// 0052, point 1). The text is byte for byte what raytatouille.results.to_json(result) gives in
/// Python with indent=None for the same result: the envelope {"format", "schema_version",
/// "type", "data"}, the keys in the order of to_dict(), separators ", " and ": ", floats as
/// Python's repr, NaN and +-infinity as {"float": "NaN"} (scalars) or "NaN" (array elements),
/// arrays as {"dtype", "shape", "values"} in C order, UTF-8 without a final newline. What a
/// result contains follows the Python binding (libs/rtt-py/src/bind_*.cpp); a test compares both
/// sides for every type. Units are those of the result types.

#include <string>
#include <string_view>

#include "rtt/analysis/opd.hpp"
#include "rtt/analysis/psf.hpp"
#include "rtt/analysis/spot.hpp"
#include "rtt/analysis/zernike.hpp"
#include "rtt/optim/optimize.hpp"
#include "rtt/paraxial/paraxial.hpp"
#include "rtt/paraxial/seidel.hpp"
#include "rtt/trace/ray_batch.hpp"
#include "rtt/trace/ray_paths.hpp"
#include "rtt/trace/sequential.hpp"

namespace rtt::io {

/// The schema version written, the same as raytatouille.results.SCHEMA_VERSION.
inline constexpr std::string_view kResultSchemaVersion = "0.1.25";

/// A double as Python's repr: the shortest digits that give the same double, in fixed notation
/// for a decimal exponent d with -4 < d <= 16 (value = 0.digits * 10^d), with ".0" for an
/// integral value, else in scientific notation with a sign and at least two exponent digits
/// (1e-05, 1e+16); -0.0 keeps its sign. Only for finite values.
[[nodiscard]] std::string python_float_repr(double x);

/// The result as JSON text in the format raytatouille-result (type name as in Python).
[[nodiscard]] std::string to_result_json(const trace::RayBatch& result);
[[nodiscard]] std::string to_result_json(const trace::RayPaths& result);
[[nodiscard]] std::string to_result_json(const trace::TraceStats& result);
[[nodiscard]] std::string to_result_json(const paraxial::FirstOrder& result);
[[nodiscard]] std::string to_result_json(const paraxial::Seidel& result);
[[nodiscard]] std::string to_result_json(const analysis::SpotDiagram& result);
[[nodiscard]] std::string to_result_json(const analysis::OpdMap& result);
[[nodiscard]] std::string to_result_json(const analysis::ZernikeFit& result);
[[nodiscard]] std::string to_result_json(const analysis::Mtf& result);
[[nodiscard]] std::string to_result_json(const analysis::Psf& result);
[[nodiscard]] std::string to_result_json(const optim::OptimResult& result);

}  // namespace rtt::io
