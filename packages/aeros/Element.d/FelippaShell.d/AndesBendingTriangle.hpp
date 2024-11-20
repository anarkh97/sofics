#ifndef _ANDESBENDINGTRIANGLE_HPP_
#define _ANDESBENDINGTRIANGLE_HPP_

#ifdef USE_EIGEN3
#include <Eigen/Core>

template<typename doublereal>
class AndesBendingTriangle
{
  public:
    static Eigen::Matrix<doublereal,9,3>
    L(doublereal x[3], doublereal y[3], doublereal clr, doublereal cqr);

    static Eigen::Matrix<doublereal,3,9>
    Bd(doublereal x[3], doublereal y[3], doublereal betab, doublereal zeta[3]);
};
#endif

#endif
