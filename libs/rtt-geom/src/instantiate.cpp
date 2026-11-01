/// @file instantiate.cpp
/// Explicit instantiations for double so that every template is compiled (and checked by
/// clang-tidy) once as part of the library.

#include "rtt/geom/asphere.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/geom/intersect.hpp"
#include "rtt/geom/phase.hpp"
#include "rtt/geom/plane.hpp"
#include "rtt/geom/shape.hpp"

namespace rtt::geom {

static_assert(IntersectableShape<Plane<double>, double>);
static_assert(IntersectableShape<Conic<double>, double>);
static_assert(IntersectableShape<EvenAsphere<double>, double>);

template class Plane<double>;
template class Conic<double>;
template class EvenAsphere<double>;
template struct Intersection<double>;
template class LinearGratingPhase<double>;
template class RadialPhasePolynomial<double>;
template double phase<double>(std::span<const PhaseFunction<double>>, double, double) noexcept;
template std::pair<double, double> phase_grad<double>(std::span<const PhaseFunction<double>>,
                                                      double,
                                                      double) noexcept;
template math::Vec3T<double> tangential_gradient<double>(std::pair<double, double>,
                                                         const math::Vec3T<double>&) noexcept;

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
template Intersection<double> intersect_newton<double, EvenAsphere<double>>(
    const EvenAsphere<double>&, const math::Vec3T<double>&, const math::Vec3T<double>&, double);
template Intersection<double> intersect<double, EvenAsphere<double>>(const EvenAsphere<double>&,
                                                                     const math::Vec3T<double>&,
                                                                     const math::Vec3T<double>&,
                                                                     double);

}  // namespace rtt::geom
