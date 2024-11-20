#ifndef _NOBENDINGTRIANGLE_HPP_
#define _NOBENDINGTRIANGLE_HPP_

#ifdef USE_EIGEN3
#include <Eigen/Core>

template<typename doublereal>
class NoBendingTriangle
{
  public:
    static Eigen::Matrix<doublereal,9,3>
    L(doublereal x[3], doublereal y[3], doublereal clr, doublereal cqr)
    {
      return Eigen::Matrix<doublereal,9,3>::Zero();
    }

    static Eigen::Matrix<doublereal,3,9>
    Bd(doublereal x[3], doublereal y[3], doublereal betab, doublereal zeta[3])
    {
      return Eigen::Matrix<doublereal,3,9>::Zero();
    }
};
#endif

#endif
