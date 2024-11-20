#ifdef USE_EIGEN3
#include <Element.d/FelippaShell.d/EffMembraneTriangleTemplate.cpp>

template
Eigen::Matrix<double,9,3>
EffMembraneTriangle<double>::L(double x[3], double y[3], double alphab);

template
Eigen::Matrix<double,3,9>
EffMembraneTriangle<double>::Bd(double x[3], double y[3], double f, double zeta[3]);
#endif
