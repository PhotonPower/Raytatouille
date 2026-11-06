/// @file instantiate.cpp
/// Explicit instantiations for double so that every template is compiled (and checked by
/// clang-tidy) once as part of the library.

#include "rtt/geom/anamorphic.hpp"
#include "rtt/geom/asphere.hpp"
#include "rtt/geom/axicon.hpp"
#include "rtt/geom/conic.hpp"
#include "rtt/geom/forbes.hpp"
#include "rtt/geom/intersect.hpp"
#include "rtt/geom/odd_asphere.hpp"
#include "rtt/geom/phase.hpp"
#include "rtt/geom/plane.hpp"
#include "rtt/geom/shape.hpp"

namespace rtt::geom {

static_assert(IntersectableShape<Plane<double>, double>);
static_assert(IntersectableShape<Conic<double>, double>);
static_assert(IntersectableShape<EvenAsphere<double>, double>);
static_assert(ParaxialShape<Plane<double>, double>);
static_assert(ParaxialShape<Conic<double>, double>);
static_assert(ParaxialShape<EvenAsphere<double>, double>);
static_assert(DomainShape<Biconic<double>, double>);
static_assert(DomainShape<Toroid<double>, double>);
static_assert(DomainShape<Cylinder<double>, double>);
static_assert(ParaxialShape<Biconic<double>, double>);
static_assert(ParaxialShape<Toroid<double>, double>);
static_assert(ParaxialShape<Cylinder<double>, double>);
static_assert(IntersectableShape<OddAsphere<double>, double>);
static_assert(IntersectableShape<Axicon<double>, double>);
static_assert(ParaxialShape<OddAsphere<double>, double>);
static_assert(ParaxialShape<Axicon<double>, double>);
static_assert(IntersectableShape<QCon<double>, double>);
static_assert(DomainShape<QBfs<double>, double>);
static_assert(ParaxialShape<QCon<double>, double>);
static_assert(ParaxialShape<QBfs<double>, double>);

template class Plane<double>;
template class Conic<double>;
template class EvenAsphere<double>;
template class Biconic<double>;
template class Toroid<double>;
template class Cylinder<double>;
template class OddAsphere<double>;
template class Axicon<double>;
template class QCon<double>;
template class QBfs<double>;
template struct QBfsMatrix<double>;
template QBfsMatrix<double> qbfs_matrix<double>(std::size_t);
template std::vector<double> qbfs_auxiliary<double>(std::span<const double>);
template std::pair<double, double> qbfs_sum<double>(std::span<const double>, double);
template std::pair<double, double> qcon_sum<double>(std::span<const double>, double);
template struct Intersection<double>;
template struct VertexExpansion<double>;
template VertexExpansion<double> stack_expansion<double>(const VertexExpansion<double>&,
                                                         std::span<const VertexExpansion<double>>);
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
template Intersection<double> intersect<double>(const Axicon<double>&,
                                                const math::Vec3T<double>&,
                                                const math::Vec3T<double>&,
                                                double);
template Intersection<double> intersect_newton<double, EvenAsphere<double>>(
    const EvenAsphere<double>&, const math::Vec3T<double>&, const math::Vec3T<double>&, double);
template Intersection<double> intersect_newton<double, Biconic<double>>(const Biconic<double>&,
                                                                        const math::Vec3T<double>&,
                                                                        const math::Vec3T<double>&,
                                                                        double);
template Intersection<double> intersect_newton<double, Toroid<double>>(const Toroid<double>&,
                                                                       const math::Vec3T<double>&,
                                                                       const math::Vec3T<double>&,
                                                                       double);
template Intersection<double> intersect_newton<double, Cylinder<double>>(const Cylinder<double>&,
                                                                         const math::Vec3T<double>&,
                                                                         const math::Vec3T<double>&,
                                                                         double);
template Intersection<double> intersect_newton<double, OddAsphere<double>>(
    const OddAsphere<double>&, const math::Vec3T<double>&, const math::Vec3T<double>&, double);
template Intersection<double> intersect_newton<double, Axicon<double>>(const Axicon<double>&,
                                                                       const math::Vec3T<double>&,
                                                                       const math::Vec3T<double>&,
                                                                       double);
template Intersection<double> intersect_newton<double, QCon<double>>(const QCon<double>&,
                                                                     const math::Vec3T<double>&,
                                                                     const math::Vec3T<double>&,
                                                                     double);
template Intersection<double> intersect_newton<double, QBfs<double>>(const QBfs<double>&,
                                                                     const math::Vec3T<double>&,
                                                                     const math::Vec3T<double>&,
                                                                     double);
template Intersection<double> intersect<double, EvenAsphere<double>>(const EvenAsphere<double>&,
                                                                     const math::Vec3T<double>&,
                                                                     const math::Vec3T<double>&,
                                                                     double);

}  // namespace rtt::geom
