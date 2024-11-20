/*
 * Copyright (C) 1997    F. Sebastian Grassia
 *
 * Permission to use and modify in any way, and for any purpose, this
 * software, is granted by the author.  Permission to redistribute
 * unmodified copies is also granted.  Modified copies may only be
 * redistributed with the express written consent of F. Sebastian Grassia,
 * F.Sebastian.Grassia@cs.cmu.edu
 *
 * This source code can be found at http://www.cs.cmu.edu/~spiff/exp-map
 */


#include <cmath>
#include <cassert>
#include "exp-map.h"
#include <iostream>

/* vector indices */
#define X 0
#define Y 1
#define Z 2
#define W 3

typedef double Quat[4];

/* crossover point to Taylor Series approximation.  Figuring 16
 * decimal digits of precision for doubles, the Taylor approximation
 * should be indistinguishable (to machine precision) for angles
 * of magnitude less than 1e-4. To be conservative, we add on three
 * additional orders of magnitude.  */
const double MIN_ANGLE = 1e-7;

/* Angle beyond which we perform dynamic reparameterization of a 3 DOF EM */
const double CUTOFF_ANGLE = M_PI;


double V3Magnitude(const double vec[3])
{
    return sqrt(vec[X]*vec[X] + vec[Y]*vec[Y] + vec[Z]*vec[Z]);
}

void V3Scale(const double v1[3], const double s1, double prod[3])
{
    prod[X] = v1[X] * s1;
    prod[Y] = v1[Y] * s1;
    prod[Z] = v1[Z] * s1;
}



/* -----------------------------------------------------------------
 * 'Q_To_Matrix' convert unit quaternion 'q' into rotation matrix 'R',
 * which transforms column vectors. Thus 'R' is the transpose of the
 * matrix found in Shoemake's 1985 paper.
 * -----------------------------------------------------------------*/
void Q_To_Matrix(Quat q, double R[3][3])
{
    double  xy, xz, yz;
    double  wx, wy, wz;
    double  xx, yy, zz;
    int     i;
    
    xy=2.*q[X]*q[Y]; xz=2.*q[X]*q[Z]; yz=2.*q[Y]*q[Z];
    wx=2.*q[W]*q[X]; wy=2.*q[W]*q[Y]; wz=2.*q[W]*q[Z];
    xx=2.*q[X]*q[X]; yy=2.*q[Y]*q[Y]; zz=2.*q[Z]*q[Z];

    R[X][X] = 1.-(yy+zz);   R[X][Y] = (xy-wz);      R[X][Z] = (xz+wy);
    R[Y][X] = (xy+wz);      R[Y][Y] = 1.-(xx+zz);   R[Y][Z] = (yz-wx);
    R[Z][X] = (xz-wy);      R[Z][Y] = (yz+wx);      R[Z][Z] = 1.-(xx+yy);
}


/* -----------------------------------------------------------------
 * 'Check_Parameterization' To escape the vanishing derivatives at
 * shells of 2PI rotations, we reparameterize to a rotation of (2PI -
 * theta) about the opposite axis when we get too close to 2PI
 * -----------------------------------------------------------------*/
int Check_Parameterization(double v[3], double *theta)
{
    int     rep = 0;
    *theta = V3Magnitude(v);

    //std::cerr << "theta = " << (*theta) << std::endl;
    if (*theta > CUTOFF_ANGLE){
	double scl = *theta;
	if (*theta > 2*M_PI){	/* first get theta into range 0..2PI */
	    *theta = fmod(*theta, 2*M_PI);
	    scl = *theta/scl;
	    V3Scale(v, scl, v);
	    rep = 1;
	}
	if (*theta > CUTOFF_ANGLE){
	    scl = *theta;
	    *theta = 2*M_PI - *theta;
	    scl = 1.0 - 2*M_PI/scl;
	    V3Scale(v, scl, v);
	    rep = 1;
	}
    }

    return rep;
}

/* -----------------------------------------------------------------
 * 'EM_To_Q' Convert a 3 DOF EM vector 'v' into its corresponding
 * quaternion 'q'. If 'reparam' is non-zero, perform dynamic
 * reparameterization, if necessary, storing the reparameterized EM in
 * 'v' and returning 1.  Returns 0 if no reparameterization was
 * performed.
 * -----------------------------------------------------------------*/
int EM_To_Q(double v[3], Quat q, int reparam)
{
    int      rep=0;
    double   cosp, sinp, theta;


    if (reparam)
      rep = Check_Parameterization(v, &theta);
    else
      theta = V3Magnitude(v);
    
    cosp = cos(.5*theta);
    sinp = sin(.5*theta);

    q[W] = cosp;
    if (theta < MIN_ANGLE)
      V3Scale(v, .5 - theta*theta/48.0, q);	/* Taylor Series for sinc */
    else
      V3Scale(v, sinp/theta, q);

    return rep;
}



/* -----------------------------------------------------------------
 * 'EM3_To_R' Convert a 3 DOF EM vector 'v' into a rotation matrix.
 * -----------------------------------------------------------------*/
void EM3_To_R(double v[3], double R[3][3])
{
    Quat     q;

    EM_To_Q(v, q, 0);
    Q_To_Matrix(q, R);
}


/* -----------------------------------------------------------------
 * 'EM2_To_EM3'Convert a 2 DOF EM into its corresponding EM 3-vector
 * -----------------------------------------------------------------*/
void EM2_To_EM3(double r[2], double s[3], double t[3], double v[3])
{
    v[X] = r[X]*s[X] + r[Y]*t[X];
    v[Y] = r[X]*s[Y] + r[Y]*t[Y];
    v[Z] = r[X]*s[Z] + r[Y]*t[Z];
}    

/* -----------------------------------------------------------------
 * 'EM2_To_R' Convert a 2 DOF EM vector 'r' with unit basis vectors
 * 's' and 't' into a rotation matrix.
 * -----------------------------------------------------------------*/
void EM2_To_R(double r[2], double s[3], double t[3], double R[3][3])
{
    double       v[3];

    EM2_To_EM3(r, s, t, v);
    EM3_To_R(v, R);
}


/* -----------------------------------------------------------------
 * 'Partial_R_Partial_Vi' Given a quaternion 'q' computed from the
 * current 2 or 3 degree of freedom EM vector 'v', and the partial
 * derivative of the quaternion with respect to the i'th element of
 * 'v' in 'dqdvi' (computed using 'Partial_Q_Partial_3V' or
 * 'Partial_Q_Partial_2V'), compute and store in 'dRdvi' the i'th
 * partial derivative of the rotation matrix 'R' with respect to the
 * i'th element of 'v'.
 * -----------------------------------------------------------------*/
void Partial_R_Partial_Vi(Quat q, Quat dqdvi, double dRdvi[3][3])
{
    double    prod[9];
    
    /* This efficient formulation is arrived at by writing out the
     * entire chain rule product dRdq * dqdv in terms of 'q' and 
     * noticing that all the entries are formed from sums of just
     * nine products of 'q' and 'dqdv' */
    prod[0] = -4*q[X]*dqdvi[X];
    prod[1] = -4*q[Y]*dqdvi[Y];
    prod[2] = -4*q[Z]*dqdvi[Z];
    prod[3] = 2*(q[Y]*dqdvi[X] + q[X]*dqdvi[Y]);
    prod[4] = 2*(q[W]*dqdvi[Z] + q[Z]*dqdvi[W]);
    prod[5] = 2*(q[Z]*dqdvi[X] + q[X]*dqdvi[Z]);
    prod[6] = 2*(q[W]*dqdvi[Y] + q[Y]*dqdvi[W]);
    prod[7] = 2*(q[Z]*dqdvi[Y] + q[Y]*dqdvi[Z]);
    prod[8] = 2*(q[W]*dqdvi[X] + q[X]*dqdvi[W]);

    /* first row, followed by second and third */
    dRdvi[0][0] = prod[1] + prod[2];
    dRdvi[0][1] = prod[3] - prod[4];
    dRdvi[0][2] = prod[5] + prod[6];

    dRdvi[1][0] = prod[3] + prod[4];
    dRdvi[1][1] = prod[0] + prod[2];
    dRdvi[1][2] = prod[7] - prod[8];

    dRdvi[2][0] = prod[5] - prod[6];
    dRdvi[2][1] = prod[7] + prod[8];
    dRdvi[2][2] = prod[0] + prod[1];

}

void Second_Partial_R_Partial_Vij(Quat q, Quat dqdvi, Quat dqdvj, Quat d2qdvidvj, double d2Rdvidvj[3][3])
{
/*
    double q1 = q[X], q2 = q[Y], q3 = q[Z], q4 = q[W];
    double dq1dvi = dqdvi[X], dq2dvi = dqdvi[Y], dq3dvi = dqdvi[Z], dq4dvi = dqdvi[W];
    double dq1dvj = dqdvj[X], dq2dvj = dqdvj[Y], dq3dvj = dqdvj[Z], dq4dvj = dqdvj[W];
    double d2q1dvidvj = d2qdvidvj[X], d2q2dvidvj = d2qdvidvj[Y], d2q3dvidvj = d2qdvidvj[Z], d2q4dvidvj = d2qdvidvj[W];

    d2Rdvidvj[0][0] = - 4*dq2dvi*dq2dvj - 4*dq3dvi*dq3dvj - 4*d2q2dvidvj*q2 - 4*d2q3dvidvj*q3;
    d2Rdvidvj[0][1] = 2*dq1dvi*dq2dvj + 2*dq1dvj*dq2dvi - 2*dq3dvi*dq4dvj - 2*dq3dvj*dq4dvi + 2*d2q2dvidvj*q1 + 2*d2q1dvidvj*q2 - 2*d2q4dvidvj*q3 - 2*d2q3dvidvj*q4;
    d2Rdvidvj[0][2] = 2*dq1dvi*dq3dvj + 2*dq1dvj*dq3dvi + 2*dq2dvi*dq4dvj + 2*dq2dvj*dq4dvi + 2*d2q3dvidvj*q1 + 2*d2q4dvidvj*q2 + 2*d2q1dvidvj*q3 + 2*d2q2dvidvj*q4;
    d2Rdvidvj[1][0] = 2*dq1dvi*dq2dvj + 2*dq1dvj*dq2dvi + 2*dq3dvi*dq4dvj + 2*dq3dvj*dq4dvi + 2*d2q2dvidvj*q1 + 2*d2q1dvidvj*q2 + 2*d2q4dvidvj*q3 + 2*d2q3dvidvj*q4;
    d2Rdvidvj[1][1] = - 4*dq1dvi*dq1dvj - 4*dq3dvi*dq3dvj - 4*d2q1dvidvj*q1 - 4*d2q3dvidvj*q3;
    d2Rdvidvj[1][2] = 2*dq2dvi*dq3dvj - 2*dq1dvj*dq4dvi - 2*dq1dvi*dq4dvj + 2*dq2dvj*dq3dvi - 2*d2q4dvidvj*q1 + 2*d2q3dvidvj*q2 + 2*d2q2dvidvj*q3 - 2*d2q1dvidvj*q4;
    d2Rdvidvj[2][0] = 2*dq1dvi*dq3dvj + 2*dq1dvj*dq3dvi - 2*dq2dvi*dq4dvj - 2*dq2dvj*dq4dvi + 2*d2q3dvidvj*q1 - 2*d2q4dvidvj*q2 + 2*d2q1dvidvj*q3 - 2*d2q2dvidvj*q4;
    d2Rdvidvj[2][1] = 2*dq1dvi*dq4dvj + 2*dq1dvj*dq4dvi + 2*dq2dvi*dq3dvj + 2*dq2dvj*dq3dvi + 2*d2q4dvidvj*q1 + 2*d2q3dvidvj*q2 + 2*d2q2dvidvj*q3 + 2*d2q1dvidvj*q4;
    d2Rdvidvj[2][2] = - 4*dq1dvi*dq1dvj - 4*dq2dvi*dq2dvj - 4*d2q1dvidvj*q1 - 4*d2q2dvidvj*q2;
*/
    double    prod[9];

    // efficient implementation
    prod[0]  = -4*(d2qdvidvj[X]*q[X] + dqdvi[X]*dqdvj[X]);
    prod[1]  = -4*(d2qdvidvj[Y]*q[Y] + dqdvi[Y]*dqdvj[Y]);
    prod[2]  = -4*(d2qdvidvj[Z]*q[Z] + dqdvi[Z]*dqdvj[Z]);
    prod[3]  = 2*(d2qdvidvj[Y]*q[X] + d2qdvidvj[X]*q[Y] + dqdvi[X]*dqdvj[Y] + dqdvj[X]*dqdvi[Y]);
    prod[4]  = 2*(d2qdvidvj[W]*q[Z] + d2qdvidvj[Z]*q[W] + dqdvi[Z]*dqdvj[W] + dqdvj[Z]*dqdvi[W]);
    prod[5]  = 2*(d2qdvidvj[Z]*q[X] + d2qdvidvj[X]*q[Z] + dqdvi[X]*dqdvj[Z] + dqdvj[X]*dqdvi[Z]);
    prod[6]  = 2*(d2qdvidvj[W]*q[Y] + d2qdvidvj[Y]*q[W] + dqdvi[Y]*dqdvj[W] + dqdvj[Y]*dqdvi[W]);
    prod[7]  = 2*(d2qdvidvj[Z]*q[Y] + d2qdvidvj[Y]*q[Z] + dqdvi[Y]*dqdvj[Z] + dqdvj[Y]*dqdvi[Z]);
    prod[8]  = 2*(d2qdvidvj[W]*q[X] + d2qdvidvj[X]*q[W] + dqdvi[X]*dqdvj[W] + dqdvj[X]*dqdvi[W]);

    d2Rdvidvj[0][0] = prod[1] + prod[2];
    d2Rdvidvj[0][1] = prod[3] - prod[4];
    d2Rdvidvj[0][2] = prod[5] + prod[6];

    d2Rdvidvj[1][0] = prod[3] + prod[4];
    d2Rdvidvj[1][1] = prod[0] + prod[2];
    d2Rdvidvj[1][2] = prod[7] - prod[8];

    d2Rdvidvj[2][0] = prod[5] - prod[6];
    d2Rdvidvj[2][1] = prod[7] + prod[8];
    d2Rdvidvj[2][2] = prod[0] + prod[1];

}


/* -----------------------------------------------------------------
 * 'Partial_Q_Partial_3V' Partial derivative of quaternion wrt i'th
 * component of EM vector 'v'
 * -----------------------------------------------------------------*/
void Partial_Q_Partial_3V(double v[3], int i, Quat dqdx)
{
    double   theta = V3Magnitude(v);
    double   cosp = cos(.5*theta), sinp = sin(.5*theta);
    
    assert(i>=0 && i<3);

    /* This is an efficient implementation of the derivatives given
     * in Appendix A of the paper with common subexpressions factored out */
    if (theta < MIN_ANGLE){
	const int i2 = (i+1)%3, i3 = (i+2)%3;
	double Tsinc = 0.5 - theta*theta/48.0;
	double vTerm = v[i] * (theta*theta/40.0 - 1.0) / 24.0;
	
	dqdx[W] = -.5*v[i]*Tsinc;
	dqdx[i]  = pow(theta,4)/3840.0 + v[i]* vTerm + Tsinc; // PJSA: added theta^4/3840 term
	dqdx[i2] = v[i2]*vTerm;
	dqdx[i3] = v[i3]*vTerm;
    }
    else{
	const int i2 = (i+1)%3, i3 = (i+2)%3;
	const double  ang = 1.0/theta, ang2 = ang*ang*v[i], sang = sinp*ang;
	const double  cterm = ang2*(.5*cosp - sang);
	
	dqdx[i]  = cterm*v[i] + sang;
	dqdx[i2] = cterm*v[i2];
	dqdx[i3] = cterm*v[i3];
	dqdx[W] = -.5*v[i]*sang;
    }
}


/* -----------------------------------------------------------------
 * 'Second_Partial_Q_Partial_3V' Second partial derivative of 
 * quaternion wrt i'th component of EM vector 'v'
 * -----------------------------------------------------------------*/
/*
void Second_Partial_Q_Partial_3V(double v[3], int i, Quat d2qdx2)
{
    double   theta = V3Magnitude(v);
    double   cosp = cos(.5*theta), sinp = sin(.5*theta);

    assert(i>=0 && i<3);

    // This is an efficient implementation of the derivatives with 
    // common subexpressions factored out
    if (theta < MIN_ANGLE){
        const int i2 = (i+1)%3, i3 = (i+2)%3;
        double Tsinc = v[i]*v[i]/32.0 + v[i2]*v[i2]/96.0 + v[i3]*v[i3]/96.0;
        double vTerm = Tsinc/10.0 - 1.0/24.0;

        d2qdx2[W] = Tsinc - 0.25;
        d2qdx2[i]  = v[i]*(vTerm + theta*theta/480.0 - 1.0/12.0);
        d2qdx2[i2] = v[i2]*vTerm;
        d2qdx2[i3] = v[i3]*vTerm;
    }
    else{
        const int i2 = (i+1)%3, i3 = (i+2)%3;
        const double  ang = 1.0/theta, ang2 = ang*ang*v[i], sang = sinp*ang;
        const double  cterm = ang2*(.5*cosp - sang);
        const double theta2 = theta*theta, ang3 = ang*ang*ang*v[i]*v[i];
        const double cterm2 = 0.5*cosp/theta2 - sang/theta2 - 1.5*ang3*cosp/theta - 0.25*ang3*sinp + 3.0*ang3*sinp/theta2;

        d2qdx2[i]  = v[i]*cterm2 + 2.0*cterm;
        d2qdx2[i2] = v[i2]*cterm2;
        d2qdx2[i3] = v[i3]*cterm2;
        d2qdx2[W] = -.5*v[i]*cterm - .5*sang;
   }
}
*/

void Second_Partial_Q_Partial_3V(double v[3], int i, int j, Quat d2qdxidxj)
{
    double   theta = V3Magnitude(v);
    double   cosp = cos(.5*theta), sinp = sin(.5*theta);
    double   theta2 = theta*theta, theta3 = theta2*theta, theta4 = theta3*theta, theta5 = theta4*theta;

    assert(i>=0 && i<3);
    assert(j>=0 && j<3);

    if (theta < MIN_ANGLE) {
      if(i == j) {
        const int i2 = (i+1)%3, i3 = (i+2)%3;
        d2qdxidxj[i ] = (v[i]*theta2)/480 - v[i]/12 + v[i]*(v[i]*v[i]/320 + v[i2]*v[i2]/960 + v[i3]*v[i3]/960 - 1/24.);
        d2qdxidxj[i2] = v[i2]*(v[i]*v[i]/320 + v[i2]*v[i2]/960 + v[i3]*v[i3]/960 - 1/24.);
        d2qdxidxj[i3] = v[i3]*(v[i]*v[i]/320 + v[i2]*v[i2]/960 + v[i3]*v[i3]/960 - 1/24.);
        d2qdxidxj[W ] = v[i]*v[i]/32 + v[i2]*v[i2]/96 + v[i3]*v[i3]/96 - .25;
      } 
      else { 
        const int k = 3-i-j;
        d2qdxidxj[i] = (v[i]*v[i]*v[j])/480 - v[j]/24 + (v[j]*theta2)/960;
        d2qdxidxj[j] = (v[i]*v[j]*v[j])/480 - v[i]/24 + (v[i]*theta2)/960;
        d2qdxidxj[k] = (v[i]*v[j]*v[k])/480;
        d2qdxidxj[W] = (v[i]*v[j])/48;
      }
    }
    else {
      if(i == j) {
        const int i2 = (i+1)%3, i3 = (i+2)%3;
        d2qdxidxj[i ] = (3*v[i]*cosp)/(2*theta2) - (3*v[i]*sinp)/theta3 - (3*v[i]*v[i]*v[i]*cosp)/(2*theta4) - (v[i]*v[i]*v[i]*sinp)/(4*theta3) + (3*v[i]*v[i]*v[i]*sinp)/theta5;
        d2qdxidxj[i2] = (v[i2]*cosp)/(2*theta2) - (v[i2]*sinp)/theta3 - (3*v[i]*v[i]*v[i2]*cosp)/(2*theta4) - (v[i]*v[i]*v[i2]*sinp)/(4*theta3) + (3*v[i]*v[i]*v[i2]*sinp)/theta5;
        d2qdxidxj[i3] = (v[i3]*cosp)/(2*theta2) - (v[i3]*sinp)/theta3 - (3*v[i]*v[i]*v[i3]*cosp)/(2*theta4) - (v[i]*v[i]*v[i3]*sinp)/(4*theta3) + (3*v[i]*v[i]*v[i3]*sinp)/theta5;
        d2qdxidxj[W ] = (v[i]*v[i]*sinp)/(2*theta3) - (v[i]*v[i]*cosp)/(4*theta2) - sinp/(2*theta);
      }
      else {
        const int k = 3-i-j;
        d2qdxidxj[i] = (v[j]*cosp)/(2*theta2) - (v[j]*sinp)/theta3 - (3*v[i]*v[i]*v[j]*cosp)/(2*theta4) - (v[i]*v[i]*v[j]*sinp)/(4*theta3) + (3*v[i]*v[i]*v[j]*sinp)/theta5;
        d2qdxidxj[j] = (v[i]*cosp)/(2*theta2) - (v[i]*sinp)/theta3 - (3*v[i]*v[j]*v[j]*cosp)/(2*theta4) - (v[i]*v[j]*v[j]*sinp)/(4*theta3) + (3*v[i]*v[j]*v[j]*sinp)/theta5;
        d2qdxidxj[k] = (3*v[i]*v[j]*v[k]*sinp)/theta5 - (v[i]*v[j]*v[k]*sinp)/(4*theta3) - (3*v[i]*v[j]*v[k]*cosp)/(2*theta4);
        d2qdxidxj[W] = (v[i]*v[j]*sinp)/(2*theta3) - (v[i]*v[j]*cosp)/(4*theta2);
      }
   }
/* work in progress
    // efficient implementation with common sub-expressions factored out
    if (theta < MIN_ANGLE) {
      if(i == j) {
        const int i2 = (i+1)%3, i3 = (i+2)%3;
        d2qdxidxj[i ] = 3*v[i]*theta2/960 - 3*v[i]/24 + v[i]*v[i]*v[i]/480;
        d2qdxidxj[i2] = v[i2]*(v[i]*v[i]/480 + theta2/960 - 1/24.);
        d2qdxidxj[i3] = v[i3]*(v[i]*v[i]/320 + v[i2]*v[i2]/960 + v[i3]*v[i3]/960 - 1/24.);
        d2qdxidxj[W ] = v[i]*v[i]/32 + v[i2]*v[i2]/96 + v[i3]*v[i3]/96 - .25;
      }
      else {
        const int k = 3-i-j;
        d2qdxidxj[i] = (v[i]*v[i]*v[j])/480 - v[j]/24 + (v[j]*theta2)/960;
        d2qdxidxj[j] = (v[i]*v[j]*v[j])/480 - v[i]/24 + (v[i]*theta2)/960;
        d2qdxidxj[k] = (v[i]*v[j]*v[k])/480;
        d2qdxidxj[W] = (v[i]*v[j])/48;
      }
    }
    else {
      if(i == j) {
        const int i2 = (i+1)%3, i3 = (i+2)%3;
        d2qdxidxj[i ] = (3*v[i]*cosp)/(2*theta2) - (3*v[i]*sinp)/theta3 - (3*v[i]*v[i]*v[i]*cosp)/(2*theta4) - (v[i]*v[i]*v[i]*sinp)/(4*theta3) + (3*v[i]*v[i]*v[i]*sinp)/theta5;
        d2qdxidxj[i2] = (v[i2]*cosp)/(2*theta2) - (v[i2]*sinp)/theta3 - (3*v[i]*v[i]*v[i2]*cosp)/(2*theta4) - (v[i]*v[i]*v[i2]*sinp)/(4*theta3) + (3*v[i]*v[i]*v[i2]*sinp)/theta5;
        d2qdxidxj[i3] = (v[i3]*cosp)/(2*theta2) - (v[i3]*sinp)/theta3 - (3*v[i]*v[i]*v[i3]*cosp)/(2*theta4) - (v[i]*v[i]*v[i3]*sinp)/(4*theta3) + (3*v[i]*v[i]*v[i3]*sinp)/theta5;
        d2qdxidxj[W ] = (v[i]*v[i]*sinp)/(2*theta3) - (v[i]*v[i]*cosp)/(4*theta2) - sinp/(2*theta);
      }
      else {
        const int k = 3-i-j;
        d2qdxidxj[i] = (v[j]*cosp)/(2*theta2) - (v[j]*sinp)/theta3 - (3*v[i]*v[i]*v[j]*cosp)/(2*theta4) - (v[i]*v[i]*v[j]*sinp)/(4*theta3) + (3*v[i]*v[i]*v[j]*sinp)/theta5;
        d2qdxidxj[j] = (v[i]*cosp)/(2*theta2) - (v[i]*sinp)/theta3 - (3*v[i]*v[j]*v[j]*cosp)/(2*theta4) - (v[i]*v[j]*v[j]*sinp)/(4*theta3) + (3*v[i]*v[j]*v[j]*sinp)/theta5;
        d2qdxidxj[k] = (3*v[i]*v[j]*v[k]*sinp)/theta5 - (v[i]*v[j]*v[k]*sinp)/(4*theta3) - (3*v[i]*v[j]*v[k]*cosp)/(2*theta4);
        d2qdxidxj[W] = (v[i]*v[j]*sinp)/(2*theta3) - (v[i]*v[j]*cosp)/(4*theta2);
      }
   }
*/
}


/* -----------------------------------------------------------------
 * 'Partial_Q_Partial_2V'Partial derivative of quaternion wrt i'th
 * component of 2 DOF EM vector 'r' with basis vectors 's' and 't'.
 * See Note about caching.
 * -----------------------------------------------------------------*/
void Partial_Q_Partial_2V(double r[2], double s[3], double t[3], int i, Quat dqdx)
{
    int      j;
    Quat     dQdV[3];
    double   v[3], *dvdr = i ? t : s;
    
    assert(i>=0 && i<2);

    /* Since each derivative of 'q' with respect to 'r' depends on all
     * derivs of 'q' with respect to 'v' (through the chain rule), we
     * SHOULD compute 'dQdV' just once and cache it somewhere with the
     * EM rotation object. */
    EM2_To_EM3(r, s, t, v);
    for (i=0; i<3; i++)
      Partial_Q_Partial_3V(v, i, dQdV[i]);

    /* this is just applying the chain rule: dQdri = dQdV * dvdri */
    for (j=0; j<4; j++)	
      dqdx[j] = ((double*)&dQdV[0])[i]*dvdr[X] + 
	        ((double*)&dQdV[1])[i]*dvdr[Y] +
		((double*)&dQdV[2])[i]*dvdr[Z];
}



/* -----------------------------------------------------------------
 * 'Partial_R_Partial_EM3'Compute the i'th partial derivative of
 * the rotation matrix with respect to EM parameter 'v', storing result
 * in 'dRdvi'.  If 'v' is near a singularity, it will be dynamically
 * reparameterized in place and the value 1 is returned; otherwise,
 * 0 is returned.
 * -----------------------------------------------------------------*/
int Partial_R_Partial_EM3(double v[3], int i, double dRdvi[3][3], int reparam)
{
    Quat  q, dqdvi;
    int   rep = EM_To_Q(v, q, reparam);

    Partial_Q_Partial_3V(v, i, dqdvi);
    Partial_R_Partial_Vi(q, dqdvi, dRdvi);

    return rep;
}

int Second_Partial_R_Partial_EM3(double v[3], int i, int j, double d2Rdvidvj[3][3], int reparam)
{
    Quat  q, dqdvi, dqdvj, d2qdvidvj;
    int   rep = EM_To_Q(v, q, reparam);

    Partial_Q_Partial_3V(v, i, dqdvi);
    Partial_Q_Partial_3V(v, j, dqdvj);
    Second_Partial_Q_Partial_3V(v, i, j, d2qdvidvj);
    Second_Partial_R_Partial_Vij(q, dqdvi, dqdvj, d2qdvidvj, d2Rdvidvj);

    return rep;
}
    
    
/* -----------------------------------------------------------------
 * 'Partial_R_Partial_EM2'Compute the i'th partial derivative of the
 * rotation matrix with respect to 2 DOF EM parameter 'r', storing
 * result in 'dRdvi'.  If used as intended in the paper, no dynamic
 * reparameterization should ever be necessary.  However, we do detect
 * and rectify the situation if it occurs, as for the 3 DOF version.
 * -----------------------------------------------------------------*/
int Partial_R_Partial_EM2(double r[3], double s[3], double t[3],
			  int i, double dRdvi[3][3])
{
    Quat      q, dqdvi;
    double    v[3];
    int       rep;

    EM2_To_EM3(r, s, t, v);
    rep = EM_To_Q(v, q, 1);
    if (rep){
	/* Since 's' and 't' are orthonormal basis vectors, we can
	 * properly reparameterize 'r' by temporarily considering it a
	 * regular EM vector in the XY plane. */
	double  tmp[3];
	double theta;
	tmp[X] = r[X]; tmp[Y] = r[Y]; tmp[Z] = 0;
	Check_Parameterization(tmp, &theta);
	r[X] = tmp[X]; r[Y] = tmp[Y];
    }
    Partial_Q_Partial_2V(r, s, t, i, dqdvi);
    Partial_R_Partial_Vi(q, dqdvi, dRdvi);

    return rep;
}


/* -----------------------------------------------------------------
 * 'Vdot' Compute the vdot necessary for dynamic simulation as a function
 * of the current EM orientation 'v' and the current angular velocity
 * 'omega'.  The results are undefined when 'v' represents a rotation
 * of 2*n*PI about any axis.
 * -----------------------------------------------------------------*/
void Vdot(double v[3], double omega[3], double vdot[3])
{
    double    theta = V3Magnitude(v);
    double    cosp = cos(.5*theta), sinp = sin(.5*theta), cotp;
    double    gamma, eta;

    if (theta < MIN_ANGLE){
	gamma = (12.0 - theta*theta) / 6.0;
	eta  = (v[X]*omega[X] + v[Y]*omega[Y] + v[Z]*omega[Z]) *
	        (60.0 + theta*theta) / 360.0;
    }
    else {
	cotp = cosp/sinp;
	gamma = theta*cotp;
	eta = (v[X]*omega[X]+v[Y]*omega[Y]+v[Z]*omega[Z])/theta * (cotp - 2.0/theta);
    }
    
    vdot[X] = .5*(gamma*omega[X] - eta*v[X] + (omega[Y]*v[Z] - omega[Z]*v[Y]));
    vdot[Y] = .5*(gamma*omega[Y] - eta*v[Y] + (omega[Z]*v[X] - omega[X]*v[Z]));
    vdot[Z] = .5*(gamma*omega[Z] - eta*v[Z] + (omega[X]*v[Y] - omega[Y]*v[X]));
}
