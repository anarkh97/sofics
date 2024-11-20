#ifdef USE_EIGEN3
#include <Element.d/FelippaShell.d/AndesBendingTriangleTemplate.cpp>

template
Eigen::Matrix<double,9,3>
AndesBendingTriangle<double>::L(double x[3], double y[3], double clr, double cqr);

template
Eigen::Matrix<double,3,9>
AndesBendingTriangle<double>::Bd(double x[3], double y[3], double f, double zeta[3]);
#endif
