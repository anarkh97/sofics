#ifdef USE_EIGEN3
#ifndef _ANDESBENDINGTRIANGLETEMPLATE_CPP_
#define _ANDESBENDINGTRIANGLETEMPLATE_CPP_

#include <Element.d/FelippaShell.d/AndesBendingTriangle.hpp>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <Eigen/Core>

template<typename doublereal>
Eigen::Matrix<doublereal,9,3>
AndesBendingTriangle<doublereal>::L(doublereal x[3], doublereal y[3], doublereal clr, doublereal cqr)
{
    // Reference:
    // "The first ANDES elements: 9-dof plate bending triangles",
    // Militello & Felippa, Comput. Methods Appl. Mech. Engrg,
    // Vol 93 (1991) pp 217-246

    // Builtin functions 
    using std::sqrt;

    // Local variables 
    int i, j;
    doublereal x0, y0, x1, y1, x2, y2, x3, y3, c12, c23, c31, s12, s31, s23, x21, 
      x13, x32, y21, y32, y13, x12, x23, x31, y12, y23, y31, cc12, cc31,
      cc23, cs12, cs31, cs23, ss12, ss23, ss31, l12, l31, l23, l21, l13, l32;
    Eigen::Matrix<doublereal,9,3> Ll;
    Eigen::Matrix<doublereal,9,3> Lq;

// .....CHECK IF [CLR] AND [CQR] SATISFY THE CONSTRAINT [CLR]+[CQR]=1 

    if (clr + cqr != 1) {
        throw std::runtime_error("\n"
          "*** FATAL ERROR in AndesBendingTriangle::L ***\n"
          "*** The factors [clr] and [cqr] violate    ***\n"
          "*** the constraint [clr]+[cqr]=1:          ***\n"
          "*** Check the calling sequence             ***\n");
    }

// .....GET THE COORDINATES OF THE CENTROID OF THE TRIANGLE 

    x0 = (x[0] + x[1] + x[2]) / 3.;
    y0 = (y[0] + y[1] + y[2]) / 3.;

// .....GET THE DISTANCES BETWEEN NODES 1, 2 AND 3 AND THE CENTROID 

    x1 = x[0] - x0; 
    x2 = x[1] - x0; 
    x3 = x[2] - x0; 

    y1 = y[0] - y0; 
    y2 = y[1] - y0; 
    y3 = y[2] - y0; 

// .....GET THE DISTANCES BETWEEN NODAL POINT X- AND Y- COORDINATES 

    x21 = x2 - x1;
    x12 = -x21;
    x32 = x3 - x2;
    x23 = -x32;
    x13 = x1 - x3;
    x31 = -x13;
    y21 = y2 - y1;
    y12 = -y21;
    y32 = y3 - y2;
    y23 = -y32;
    y13 = y1 - y3;
    y31 = -y13;

// .....GET THE DISTANCES BETWEEN NODES 1-2, 2-3 AND 3-1 

    l12 = l21 = sqrt(x12 * x12 + y12 * y12);
    l23 = l32 = sqrt(x23 * x23 + y23 * y23);
    l31 = l13 = sqrt(x31 * x31 + y31 * y31);

// .....ASSEMBLE THE LOCAL MATRIX [LLR] W/ SHAPE FUNCTION DERIVATIVES 

    if (clr != 0) { // linear normal rotation (eq. 52)

        Ll.setZero();

	Ll(2, 0) = y32 * .5;
	Ll(5, 0) = y13 * .5;
	Ll(8, 0) = y21 * .5;

	Ll(1, 1) = x32 * .5;
	Ll(4, 1) = x13 * .5;
	Ll(7, 1) = x21 * .5;

	Ll(1, 2) = y23 * .5;
	Ll(2, 2) = x23 * .5;
	Ll(4, 2) = y31 * .5;
	Ll(5, 2) = x31 * .5;
	Ll(7, 2) = y12 * .5;
	Ll(8, 2) = x12 * .5;

    }

// .....ASSEMBLE THE LOCAL MATRIX [LQR] W/ SHAPE FUNCTION DERIVATIVES 

    if (cqr != 0) { // quadratic normal rotation (eq. 56)
 
        c12 = x21 / l12;
        s12 = y21 / l12;
        c23 = x32 / l23;
        s23 = y32 / l23;
        c31 = x13 / l31;
        s31 = y13 / l31;

        ss12 = s12 * s12;
        ss23 = s23 * s23;
        ss31 = s31 * s31;
        cc12 = c12 * c12;
        cc23 = c23 * c23;
        cc31 = c31 * c31;
        cs12 = c12 * s12;
        cs23 = c23 * s23;
        cs31 = c31 * s31;

	Lq(0, 0) = -cs12 + cs31;
	Lq(0, 1) = -Lq(0, 0);
	Lq(0, 2) = ss31 - cc31 - (ss12 - cc12);

	Lq(1, 0) = (ss12 * x12 + ss31 * x31) * .5;
	Lq(1, 1) = (cc12 * x12 + cc31 * x31) * .5;
	Lq(1, 2) = cc12 * y21 + cc31 * y13;

	Lq(2, 0) = (ss12 * y21 + ss31 * y13) * -.5;
	Lq(2, 1) = Lq(1, 2) * -.5;
	Lq(2, 2) = Lq(1, 0) * -2.;

	Lq(3, 0) = -cs23 + cs12;
	Lq(3, 1) = -Lq(3, 0);
	Lq(3, 2) = ss12 - cc12 - (ss23 - cc23);

	Lq(4, 0) = (ss12 * x12 + ss23 * x23) * .5;
	Lq(4, 1) = (cc12 * x12 + cc23 * x23) * .5;
	Lq(4, 2) = cc12 * y21 + cc23 * y32;

	Lq(5, 0) = (ss12 * y21 + ss23 * y32) * -.5;
	Lq(5, 1) = Lq(4, 2) * -.5;
	Lq(5, 2) = Lq(4, 0) * -2.;

	Lq(6, 0) = -cs31 + cs23;
	Lq(6, 1) = -Lq(6, 0);
	Lq(6, 2) = ss23 - cc23 - (ss31 - cc31);

	Lq(7, 0) = (ss23 * x23 + ss31 * x31) * .5;
	Lq(7, 1) = (cc23 * x23 + cc31 * x31) * .5;
	Lq(7, 2) = cc23 * y32 + cc31 * y13;

	Lq(8, 0) = (ss23 * y32 + ss31 * y13) * -.5;
	Lq(8, 1) = Lq(7, 2) * -.5;
	Lq(8, 2) = Lq(7, 0) * -2.;

    }

// .....ASSEMBLE THE MATRIX [Lb] AS -([CLR]*[LLR] + [CQR]*[LQR])

    if (clr == 0) {
        return -Lq;
    }

    else if (cqr == 0) {
        return -Ll;
    }

    else { // clr != 0 && cqr != 0
        return -(clr*Ll + cqr*Lq);
    }

}

template<typename doublereal>
Eigen::Matrix<doublereal,3,9>
AndesBendingTriangle<doublereal>::Bd(doublereal x[3], doublereal y[3], doublereal beta, doublereal zeta[3])
{
    // Reference:
    // "The first ANDES elements: 9-dof plate bending triangles",
    // Militello & Felippa, Comput. Methods Appl. Mech. Engrg,
    // Vol 93 (1991) pp 217-246

    // Builtin functions 
    using std::sqrt;

    // Local variables 
    int i, j;
    doublereal twicearea, x0, y0, x1, x2, x3, y1, y2, y3, x21, x13, x32, x12, x31, x23,
      y21, y32, y13, y12, y23, y31, lambda13, lambda21, lambda32, lambda31, lambda12, lambda23,
      area, x2ap3, l12, l21, l13, l31, l23, l32, zeta10, zeta20, zeta30;

    Eigen::Matrix<doublereal,6,9> Q;
    Eigen::Matrix<doublereal,3,3> T;
    Eigen::Matrix<doublereal,3,6> A;

// .....RETURN IF THE SCALING FACTOR IS ZERO 

    if (beta == 0) {
        return Eigen::Matrix<doublereal,3,9>::Zero();
    }

// .....CHECK IF THE SCALING FACTOR [beta] IS POSITIVE 

    if (beta < 0) {
        throw std::runtime_error("\n"
          "*** FATAL ERROR in AndesBendingTriangle::Bd ***\n"
          "*** The scaling factor [beta] is negative   ***\n"
          "*** Check the calling sequence:             ***\n"
          "*** Factor [beta] must be positive or zero  ***\n");
    }

// .....GET THE COORDINATES OF THE CENTROID OF THE TRIANGLE 

    x0 = (x[0] + x[1] + x[2]) / 3.;
    y0 = (y[0] + y[1] + y[2]) / 3.;

// .....GET THE DISTANCES BETWEEN NODES 1, 2 AND 3 AND THE CENTROID 

    x1 = x[0] - x0;
    x2 = x[1] - x0;
    x3 = x[2] - x0;

    y1 = y[0] - y0;
    y2 = y[1] - y0;
    y3 = y[2] - y0;

// .....GET THE DISTANCES BETWEEN NODAL POINT X- AND Y- COORDINATES 

    x21 = x2 - x1;
    x12 = -x21;
    x32 = x3 - x2;
    x23 = -x32;
    x13 = x1 - x3;
    x31 = -x13;

    y21 = y2 - y1;
    y12 = -y21;
    y32 = y3 - y2;
    y23 = -y32;
    y13 = y1 - y3;
    y31 = -y13;

// .....CALCULATE TWICE THE AREA OF THE TRIANGLE 

    twicearea = y21 * x13 - x21 * y13;

// .....CALCULATE THE AREA OF THE TRIANGLE 

    area = twicearea * .5;

// .....GET THE DISTANCES BETWEEN NODES 1-2, 2-3 AND 3-1 

    l12 = l21 = sqrt(x21 * x21 + y21 * y21);
    l23 = l32 = sqrt(x32 * x32 + y32 * y32);
    l31 = l13 = sqrt(x13 * x13 + y13 * y13);

// .....GET THE SIDE PROJECTION RATIOS 

    lambda31 = (x31 * x32 + y13 * y23) / (l31 * l31);
    lambda12 = (x12 * x13 + y21 * y31) / (l12 * l12);
    lambda23 = (x23 * x21 + y32 * y12) / (l23 * l23);

    lambda13 = 1. - lambda31;
    lambda21 = 1. - lambda12;
    lambda32 = 1. - lambda23;

// .....FORM THE MATRIX [Q] W/ SHAPE FUNCTION DERIVATIVES 

    Q(0, 0) = -6.;
    Q(0, 1) = y21 * -4.;
    Q(0, 2) = x21 * 4.;
    Q(0, 3) = 6.;
    Q(0, 4) = y21 * -2.;
    Q(0, 5) = x21 * 2.;
    Q(0, 6) = 0.;
    Q(0, 7) = 0.;
    Q(0, 8) = 0.;

    Q(1, 0) = 6.;
    Q(1, 1) = y21 * 2.;
    Q(1, 2) = x21 * -2.;
    Q(1, 3) = -6.;
    Q(1, 4) = y21 * 4.;
    Q(1, 5) = x21 * -4.;
    Q(1, 6) = 0.;
    Q(1, 7) = 0.;
    Q(1, 8) = 0.;

    Q(2, 0) = 0.;
    Q(2, 1) = 0.;
    Q(2, 2) = 0.;
    Q(2, 3) = -6.;
    Q(2, 4) = y32 * -4.;
    Q(2, 5) = x32 * 4.;
    Q(2, 6) = 6.;
    Q(2, 7) = y32 * -2.;
    Q(2, 8) = x32 * 2.;

    Q(3, 0) = 0.;
    Q(3, 1) = 0.;
    Q(3, 2) = 0.;
    Q(3, 3) = 6.;
    Q(3, 4) = y32 * 2.;
    Q(3, 5) = x32 * -2.;
    Q(3, 6) = -6.;
    Q(3, 7) = y32 * 4.;
    Q(3, 8) = x32 * -4.;

    Q(4, 0) = 6.;
    Q(4, 1) = y13 * -2.;
    Q(4, 2) = x13 * 2.;
    Q(4, 3) = 0.;
    Q(4, 4) = 0.;
    Q(4, 5) = 0.;
    Q(4, 6) = -6.;
    Q(4, 7) = y13 * -4.;
    Q(4, 8) = x13 * 4.;

    Q(5, 0) = -6.;
    Q(5, 1) = y13 * 4.;
    Q(5, 2) = x13 * -4.;
    Q(5, 3) = 0.;
    Q(5, 4) = 0.;
    Q(5, 5) = 0.;
    Q(5, 6) = 6.;
    Q(5, 7) = y13 * 2.;
    Q(5, 8) = x13 * -2.;


// .....GET THE MATRIX [T] THAT REPRESENTS THE INVERSE OF THE MATRIX 
// .....RELATING INSIDE CURVATURES WITH BOUNDARY CURVATURES 

    x2ap3 = twicearea * twicearea * twicearea;

    T(0, 0) = (x13 * y13 * y32 * y32 - x32 * y13 * y13 * y32) / x2ap3;
    T(0, 1) = (x21 * y21 * y13 * y13 - x13 * y21 * y21 * y13) / x2ap3;
    T(0, 2) = (-x21 * y21 * y32 * y32 + x32 * y21 * y21 * y32) / x2ap3;

    T(1, 0) = (-x13 * x32 * x32 * y13 + x13 * x13 * x32 * y32) / x2ap3;
    T(1, 1) = (-x21 * x13 * x13 * y21 + x21 * x21 * x13 * y13) / x2ap3;
    T(1, 2) = (x21 * x32 * x32 * y21 - x21 * x21 * x32 * y32) / x2ap3;

    T(2, 0) = (-x13 * x13 * y32 * y32 + x32 * x32 * y13 * y13) / x2ap3;
    T(2, 1) = (-x21 * x21 * y13 * y13 + x13 * x13 * y21 * y21) / x2ap3;
    T(2, 2) = (x21 * x21 * y32 * y32 - x32 * x32 * y21 * y21) / x2ap3;

// .....ESTIMATE THE MATRIX [A] AT THE GAUSS INTEGRATION POINT (OR NODE) 
    zeta10 = zeta[0] - 1/3.;
    zeta20 = zeta[1] - 1/3.;
    zeta30 = zeta[2] - 1/3.;

    // NOTE: this doesn't quite match equation (47) in the reference (the lambda indices are reversed)
    // it remains to be seen which one is correct

    A(0, 0) = zeta10 + lambda21 * zeta30;
    A(0, 1) = zeta20 + lambda12 * zeta30;
    A(0, 2) = 0;
    A(0, 3) = 0;
    A(0, 4) = 0;
    A(0, 5) = 0;

    A(1, 0) = 0;
    A(1, 1) = 0;
    A(1, 2) = zeta20 + lambda32 * zeta10;
    A(1, 3) = zeta30 + lambda23 * zeta10;
    A(1, 4) = 0;
    A(1, 5) = 0;

    A(2, 0) = 0;
    A(2, 1) = 0;
    A(2, 2) = 0;
    A(2, 3) = 0;
    A(2, 4) = zeta30 + lambda13 * zeta20;
    A(2, 5) = zeta10 + lambda31 * zeta20;
//
//    A(0, 0) = zeta10 + lambda12 * zeta30;
//    A(0, 1) = zeta20 + lambda21 * zeta30;
//    A(0, 2) = 0;
//    A(0, 3) = 0;
//    A(0, 4) = 0;
//    A(0, 5) = 0;
//
//    A(1, 0) = 0;
//    A(1, 1) = 0;
//    A(1, 2) = zeta20 + lambda23 * zeta10;
//    A(1, 3) = zeta30 + lambda32 * zeta10;
//    A(1, 4) = 0;
//    A(1, 5) = 0;
//
//    A(2, 0) = 0;
//    A(2, 1) = 0;
//    A(2, 2) = 0;
//    A(2, 3) = 0;
//    A(2, 4) = zeta30 + lambda31 * zeta20;
//    A(2, 5) = zeta10 + lambda13 * zeta20;


// .....ASSEMBLE MATRIX [Bdb] FOR THE GAUSS INTEGRATION POINT (OR NODE) 
// .....[Bdb] = sqrt(beta) * [T] * [A] * [Q] 

    return sqrt(beta)*T*A*Q;
}
#endif
#endif
