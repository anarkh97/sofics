#ifdef USE_EIGEN3
#ifndef _EFFMEMBRANETRIANGLETEMPLATE_CPP_
#define _EFFMEMBRANETRIANGLETEMPLATE_CPP_

#include <Element.d/FelippaShell.d/EffMembraneTriangle.hpp>
#include <cmath>
#include <stdexcept>
#include <Eigen/Core>

template<typename doublereal>
Eigen::Matrix<doublereal,9,3>
EffMembraneTriangle<doublereal>::L(doublereal x[3], doublereal y[3], doublereal alphab) 
{
    // Reference:
    // "Membrane triangles with corner drilling freedoms I. The EFF element"
    // Alvin, de la Fuente, Haugen & Felippa,
    // Finite Elements in Analysis and Design Vol. 12 (1992) pp. 163-187

    // Local variables 
    doublereal x21, x13, x32, y21, y32, y13, x12, x23, x31, y12, y23, y31;

    Eigen::Matrix<doublereal,9,3> Lm;

// .....GET THE DISTANCES BETWEEN NODAL POINT X- AND Y- COORDINATES 

    x21 = x[1] - x[0];
    x12 = -x21;
    x32 = x[2] - x[1];
    x23 = -x32;
    x13 = x[0] - x[2];
    x31 = -x13;
    y21 = y[1] - y[0];
    y12 = -y21;
    y32 = y[2] - y[1];
    y23 = -y32;
    y13 = y[0] - y[2];
    y31 = -y13;

// .....ASSEMBLE THE MATRIX [Lm] W/ SHAPE FUNCTION DERIVATIVES
//      see Alvin et al equation 15

    Lm(0, 0) = y23;
    Lm(1, 0) = 0;
    Lm(2, 0) = y31;
    Lm(3, 0) = 0;
    Lm(4, 0) = y12;
    Lm(5, 0) = 0;

    Lm(0, 1) = 0;
    Lm(1, 1) = x32;
    Lm(2, 1) = 0;
    Lm(3, 1) = x13;
    Lm(4, 1) = 0;
    Lm(5, 1) = x21;

    Lm(0, 2) = x32;
    Lm(1, 2) = y23;
    Lm(2, 2) = x13;
    Lm(3, 2) = y31;
    Lm(4, 2) = x21;
    Lm(5, 2) = y12;

    Lm(6, 0) = y23 * (y13 - y21) * alphab / 6.;
    Lm(6, 1) = x32 * (x31 - x12) * alphab / 6.;
    Lm(6, 2) = (x31 * y13 - x12 * y21) * alphab / 3.;

    Lm(7, 0) = y31 * (y21 - y32) * alphab / 6.;
    Lm(7, 1) = x13 * (x12 - x23) * alphab / 6.;
    Lm(7, 2) = (x12 * y21 - x23 * y32) * alphab / 3.;

    Lm(8, 0) = y12 * (y32 - y13) * alphab / 6.;
    Lm(8, 1) = x21 * (x23 - x31) * alphab / 6.;
    Lm(8, 2) = (x23 * y32 - x31 * y13) * alphab / 3.;

    return Lm/2;
}

template<typename doublereal>
Eigen::Matrix<doublereal,3,9>
EffMembraneTriangle<doublereal>::Bd(doublereal x[3], doublereal y[3], doublereal beta, doublereal zeta[3])
{
    // Reference:
    // "Membrane triangles with corner drilling freedoms I. The EFF element"
    // Alvin, de la Fuente, Haugen & Felippa,
    // Finite Elements in Analysis and Design Vol. 12 (1992) pp. 163-187
    // Note: (1-gamma) in the reference is equal to the function argument beta

    // Initialized data 
    doublereal zero = 0.;

    // Builtin functions 
    using std::sqrt;

    // Local variables 
    int i, j;
    doublereal twicearea, x0, y0,
      ca, x10, x20, x12, x21, x23, x32, x31, x13, y12, y21, y23,
      y32, y31, y13, x30, y10, y20, y30, aa12, aa31, aa23,
      ss12, ss31, ss23, caa12, caa31, caa23, area, 
      cax10, cax20, cax30, cay10, cay20, cay30, sum123, sum456;

    Eigen::Matrix<doublereal,3,6> Bhbar;
    Eigen::Matrix<doublereal,6,9> Hmv;
    Eigen::Matrix<doublereal,6,9> Hh;
    Eigen::Matrix<doublereal,6,1> Z;

// .....RETURN IF THE SCALING FACTOR [beta] IS ZERO 

    if (beta == zero) {
        return Eigen::Matrix<doublereal,3,9>::Zero();
    }

// .....CHECK IF THE SCALING FACTOR [beta] IS POSITIVE 

    if (beta < zero) {
        throw std::runtime_error("\n"
          "*** FATAL ERROR in EffMembraneTriangle::Bd ***\n"
          "*** The scaling factor [beta] is negative  ***\n"
          "*** Check the calling sequence:            ***\n"
          "*** Factor [beta] must be positive or zero ***\n");
    }

// .....GET THE DISTANCES BETWEEN NODAL POINT X- AND Y- COORDINATES 

    x21 = x[1] - x[0];
    x12 = -x21;
    x32 = x[2] - x[1];
    x23 = -x32;
    x13 = x[0] - x[2];
    x31 = -x13;
    y21 = y[1] - y[0];
    y12 = -y21;
    y32 = y[2] - y[1];
    y23 = -y32;
    y13 = y[0] - y[2];
    y31 = -y13;

// .....CALCULATE TWICE THE AREA OF THE TRIANGLE 

    twicearea = y21 * x13 - x21 * y13;

// .....CALCULATE THE AREA OF THE TRIANGLE 

    area = twicearea * .5;

// .....GET THE COORDINATES OF THE CENTROID OF THE TRIANGLE 

    x0 = (x[0] + x[1] + x[2]) / 3.;
    y0 = (y[0] + y[1] + y[2]) / 3.;

// .....GET THE DISTANCES BETWEEN NODES 1, 2 AND 3 AND THE CENTROID 

    x10 = x[0] - x0;
    x20 = x[1] - x0;
    x30 = x[2] - x0;

    y10 = y[0] - y0;
    y20 = y[1] - y0;
    y30 = y[2] - y0;

// .....CALCULATE BASIC COEFFICIENTS FOR THE ASSEMBLY OF MATRIX [HMV] 

    aa12 = (x30 * x30 + y30 * y30) * 2.25;
    aa23 = (x10 * x10 + y10 * y10) * 2.25;
    aa31 = (x20 * x20 + y20 * y20) * 2.25;

    caa12 = 15. / (aa12 * 128.);
    caa23 = 15. / (aa23 * 128.);
    caa31 = 15. / (aa31 * 128.);

    ss12 = x12 * x12 + y12 * y12;
    ss23 = x23 * x23 + y23 * y23;
    ss31 = x31 * x31 + y31 * y31;

    ca = 3. / (area * 16.);

    cax10 = ca * x10;
    cax20 = ca * x20;
    cax30 = ca * x30;

    cay10 = ca * y10;
    cay20 = ca * y20;
    cay30 = ca * y30;

// .....CONSTRUCT LOCAL MATRIX [HMV] W/ SHAPE FUNCTION DERIVATIVES
//      Although it's not given in Alvin et al in this form, presumably Hmv = Hms * Hst * Hrtheta * Hthetav
//      Most likely this has alphah = 5/4 built in, which is optimal for pure bending

    Hmv(0, 0) = cay30 * x32;
    Hmv(0, 1) = cay30 * y32;
    Hmv(0, 2) = cay30 * x13;
    Hmv(0, 3) = cay30 * y13;
    Hmv(0, 4) = cay30 * x21;
    Hmv(0, 5) = cay30 * y21;
    Hmv(0, 6) = (ss23 - ss31 + aa12 * 2.4) * y30 * caa12 + area * 4. * x30 * caa12;
    Hmv(0, 7) = y30 * .5625 - Hmv(0, 6);
    Hmv(0, 8) = y30 * .1875;

    Hmv(1, 0) = cay10 * x32;
    Hmv(1, 1) = cay10 * y32;
    Hmv(1, 2) = cay10 * x13;
    Hmv(1, 3) = cay10 * y13;
    Hmv(1, 4) = cay10 * x21;
    Hmv(1, 5) = cay10 * y21;
    Hmv(1, 6) = y10 * .1875;
    Hmv(1, 7) = (ss31 - ss12 + aa23 * 2.4) * y10 * caa23 + area * 4. * x10 * caa23;
    Hmv(1, 8) = y10 * .5625 - Hmv(1, 7);

    Hmv(2, 0) = cay20 * x32;
    Hmv(2, 1) = cay20 * y32;
    Hmv(2, 2) = cay20 * x13;
    Hmv(2, 3) = cay20 * y13;
    Hmv(2, 4) = cay20 * x21;
    Hmv(2, 5) = cay20 * y21;
    Hmv(2, 6) = (ss23 - ss12 + aa31 * 2.4) * y20 * caa31 - area * 4. * x20 * caa31;
    Hmv(2, 7) = y20 * .1875;
    Hmv(2, 8) = y20 * .5625 - Hmv(2, 6);

    Hmv(3, 0) = -cax30 * x32;
    Hmv(3, 1) = -cax30 * y32;
    Hmv(3, 2) = -cax30 * x13;
    Hmv(3, 3) = -cax30 * y13;
    Hmv(3, 4) = -cax30 * x21;
    Hmv(3, 5) = -cax30 * y21;
    Hmv(3, 6) = (ss31 - ss23 - aa12 * 2.4) * x30 * caa12 + area * 4. * y30 * caa12;
    Hmv(3, 7) = x30 * -.5625 - Hmv(3, 6);
    Hmv(3, 8) = x30 * -.1875;

    Hmv(4, 0) = -cax10 * x32;
    Hmv(4, 1) = -cax10 * y32;
    Hmv(4, 2) = -cax10 * x13;
    Hmv(4, 3) = -cax10 * y13;
    Hmv(4, 4) = -cax10 * x21;
    Hmv(4, 5) = -cax10 * y21;
    Hmv(4, 6) = x10 * -.1875;
    Hmv(4, 7) = (ss12 - ss31 - aa23 * 2.4) * x10 * caa23 + area * 4. * y10 * caa23;
    Hmv(4, 8) = x10 * -.5625 - Hmv(4, 7);

    Hmv(5, 0) = -cax20 * x32;
    Hmv(5, 1) = -cax20 * y32;
    Hmv(5, 2) = -cax20 * x13;
    Hmv(5, 3) = -cax20 * y13;
    Hmv(5, 4) = -cax20 * x21;
    Hmv(5, 5) = -cax20 * y21;
    Hmv(5, 6) = (ss12 - ss23 - aa31 * 2.4) * x20 * caa31 - area * 4. * y20 * caa31;
    Hmv(5, 7) = x20 * -.1875;
    Hmv(5, 8) = x20 * -.5625 - Hmv(5, 6);

// .....CONSTRUCT LOCAL MATRIX [Hh] FROM [Hmv]
//      This is an optimized implementaion of the matrix matrix product: Hh = Hqm * Hmv

    for (j = 0; j < 9; ++j) {

        sum123 = (Hmv(0, j) + Hmv(1, j) + Hmv(2, j)) * 2/9.;

        Hh(0, j) = sum123 - Hmv(0, j) * 12/9.;
        Hh(1, j) = sum123 - Hmv(1, j) * 12/9.;
        Hh(2, j) = sum123 - Hmv(2, j) * 12/9.;

        sum456 = (Hmv(3, j) + Hmv(4, j) + Hmv(5, j)) * 2/9.;

        Hh(3, j) = sum456 - Hmv(3, j) * 12/9.;
        Hh(4, j) = sum456 - Hmv(4, j) * 12/9.;
        Hh(5, j) = sum456 - Hmv(5, j) * 12/9.;

    }

// .....FORM THE LOCAL MATRIX [B] Reference: Alvin et al equation 26

    Bhbar(0, 0) = y30;
    Bhbar(1, 0) = zero;
    Bhbar(2, 0) = -x30;

    Bhbar(0, 1) = y10;
    Bhbar(1, 1) = zero;
    Bhbar(2, 1) = -x10;

    Bhbar(0, 2) = y20;
    Bhbar(1, 2) = zero;
    Bhbar(2, 2) = -x20;

    Bhbar(0, 3) = zero;
    Bhbar(1, 3) = -x30;
    Bhbar(2, 3) = y30;

    Bhbar(0, 4) = zero;
    Bhbar(1, 4) = -x10;
    Bhbar(2, 4) = y10;

    Bhbar(0, 5) = zero;
    Bhbar(1, 5) = -x20;
    Bhbar(2, 5) = y20;

    Bhbar *= 3/area;

// .....CALCULATE THE DIAGONAL TRIANGULAR COORDINATE MATRIX [Z]
// .....AT THE GAUSS INTEGRATION POINT (OR NODE) Reference: Alvin et al equation 26

    Z[0] = zeta[1] - zeta[0];
    Z[1] = zeta[2] - zeta[1];
    Z[2] = zeta[0] - zeta[2];
    Z[3] = zeta[1] - zeta[0];
    Z[4] = zeta[2] - zeta[1];
    Z[5] = zeta[0] - zeta[2];

// .....FORM THE LOCAL PRODUCT: 
// .....[Bdm] = -sqrt(beta) * [Bhbar] * [Z] * [Hh] 

    return -sqrt(beta) * (Bhbar * Z.asDiagonal() * Hh);
}
#endif
#endif
