/// @file instantiate.cpp
/// Explicit instantiations for double so that every template is compiled (and checked by
/// clang-tidy) once as part of the library.

#include "rtt/geom/conic.hpp"
#include "rtt/geom/intersect.hpp"
#include "rtt/geom/plane.hpp"
#include "rtt/geom/shape.hpp"

namespace rtt::geom {

template class Shape<double>;
template class Plane<double>;
template class Conic<double>;
template struct Intersection<double>;

template Intersection<double> intersect_conic<double>(
    double, double, const math::Vec3T<double>&, const math::Vec3T<double>&, double);
template Intersection<double> intersect<double>(const Conic<double>&,
                                                const math::Vec3T<double>&,
                                                const math::Vec3T<double>&,
                                                double);
template Intersection<double> intersect<double>(const Plane<double>&,
                                                const math::Vec3T<double>&,
                                                const math::Vec3T<double>&,
                                                double);

}  // namespace rtt::geom
