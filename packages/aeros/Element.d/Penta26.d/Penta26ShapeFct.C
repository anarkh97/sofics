// ---------------------------------------------------------------------
// HB - 05-22-05
// ---------------------------------------------------------------------
// Shape fcts & derivatives for 26 nodes wedge element 
// Serendipity finite element basis
// Iso-parametric formulation
// ---------------------------------------------------------------------
// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <Utils.d/linkfc.h>

// HB (05-21-05): for given reference shape fct derivatives dShape & nodes coordinates (X.Y,Z), 
//                compute the jacobian & shape fcts derivative w.r.t to the global coordinate system.
//                assuming an iso-parametric formulation.
double
computePenta26DShapeFct(double dShape[26][3], double X[26], double Y[26], double Z[26], double (*DShape)[3]=0)
{
  double xd1,xd2,xd3;
  double yd1,yd2,yd3;
  double zd1,zd2,zd3;

  xd1 = 0.0; xd2 = 0.0; xd3 = 0.0;
  yd1 = 0.0; yd2 = 0.0; yd3 = 0.0;
  zd1 = 0.0; zd2 = 0.0; zd3 = 0.0;

  for(int i=0; i<26; i+=2){
    xd1 += dShape[i][0]*X[i] + dShape[i+1][0]*X[i+1];
    xd2 += dShape[i][1]*X[i] + dShape[i+1][1]*X[i+1];
    xd3 += dShape[i][2]*X[i] + dShape[i+1][2]*X[i+1];
    
    yd1 += dShape[i][0]*Y[i] + dShape[i+1][0]*Y[i+1];
    yd2 += dShape[i][1]*Y[i] + dShape[i+1][1]*Y[i+1];
    yd3 += dShape[i][2]*Y[i] + dShape[i+1][2]*Y[i+1];
    
    zd1 += dShape[i][0]*Z[i] + dShape[i+1][0]*Z[i+1];
    zd2 += dShape[i][1]*Z[i] + dShape[i+1][1]*Z[i+1];
    zd3 += dShape[i][2]*Z[i] + dShape[i+1][2]*Z[i+1];
  }

  double a11 = yd2*zd3 - yd3*zd2; double a12 = yd3*zd1 - yd1*zd3; double a13 = yd1*zd2 - yd2*zd1;
  double a21 = xd3*zd2 - xd2*zd3; double a22 = xd1*zd3 - zd1*xd3; double a23 = xd2*zd1 - xd1*zd2;
  double a31 = xd2*yd3 - xd3*yd2; double a32 = yd1*xd3 - xd1*yd3; double a33 = xd1*yd2 - yd1*xd2;
 
  /* -------> DETERMINANT OF THE JACOBIAN MATRIX <--- */

  double J = xd1*a11 + yd1*a21 + zd1*a31 ;
  if(J == 0.0) {
    fprintf(stderr," *** WARNING: NULL JACOBIAN IN computePenta26DShapeFct.C ROUTINE.\n");
  }
  if(J < 0.0) {
    fprintf(stderr," *** WARNING: NEGATIVE JACOBIAN IN computePenta26DShapeFct.C ROUTINE.\n");
  }
  double cdet = 1.0/J;

  /* -------> DERIVATIVE OF THE SHAPE FCTS IN THE "REAL" ELEMENT <--- */
  if(DShape) {
    for(int i=0; i<26; i+=2) {
      DShape[i  ][0] = cdet * ( a11*dShape[i  ][0] + a12*dShape[i  ][1] + a13*dShape[i  ][2] );
      DShape[i  ][1] = cdet * ( a21*dShape[i  ][0] + a22*dShape[i  ][1] + a23*dShape[i  ][2] );
      DShape[i  ][2] = cdet * ( a31*dShape[i  ][0] + a32*dShape[i  ][1] + a33*dShape[i  ][2] );

      DShape[i+1][0] = cdet * ( a11*dShape[i+1][0] + a12*dShape[i+1][1] + a13*dShape[i+1][2] );
      DShape[i+1][1] = cdet * ( a21*dShape[i+1][0] + a22*dShape[i+1][1] + a23*dShape[i+1][2] );
      DShape[i+1][2] = cdet * ( a31*dShape[i+1][0] + a32*dShape[i+1][1] + a33*dShape[i+1][2] );
    } 
  }

  return(J);
}

// z in [-1,1]
void Penta26ShapeFct(double Shape[26], double m[3])
{
  double r = m[0];
  double s = m[1];
  double z = m[2];

  double t = 1.-r-s;
  
  double zzm = 1.-z*z;
  double zm  = 1.-z;
  double zp  = 1.+z;
  double z3m = 1.-3.*z;
  double z3p = 1.+3.*z;
  
  double c0  = 2./9.;
  double c1  = 9./16.;
  double c2  = 9./4.;
  double c3  = 27./2.;
  
  // Lower & upper wedge's corners  
  Shape[ 0] = c1*zm*t*(4.*(t*(t-1.)+c0) - zzm);
  Shape[ 1] = c1*zm*r*(4.*(r*(r-1.)+c0) - zzm);
  Shape[ 2] = c1*zm*s*(4.*(s*(s-1.)+c0) - zzm);
  Shape[ 3] = c1*zp*t*(4.*(t*(t-1.)+c0) - zzm);
  Shape[ 4] = c1*zp*r*(4.*(r*(r-1.)+c0) - zzm);
  Shape[ 5] = c1*zp*s*(4.*(s*(s-1.)+c0) - zzm);

  // Lower & upper wedge's triangular edges' nodes
  Shape[ 6] = c2*r*t*(3.*t-1.)*zm;
  Shape[ 7] = c2*r*t*(3.*r-1.)*zm;
  Shape[ 8] = c2*r*s*(3.*r-1.)*zm;
  Shape[ 9] = c2*r*s*(3.*s-1.)*zm;
  Shape[10] = c2*t*s*(3.*s-1.)*zm;
  Shape[11] = c2*t*s*(3.*t-1.)*zm;

  Shape[12] = c2*r*t*(3.*t-1.)*zp;
  Shape[13] = c2*r*t*(3.*r-1.)*zp;
  Shape[14] = c2*r*s*(3.*r-1.)*zp;
  Shape[15] = c2*r*s*(3.*s-1.)*zp;
  Shape[16] = c2*t*s*(3.*s-1.)*zp;
  Shape[17] = c2*t*s*(3.*t-1.)*zp;

  // Wedge vertical edges' nodes
  Shape[18] = c1*t*zzm*z3m;
  Shape[19] = c1*t*zzm*z3p;
  Shape[20] = c1*r*zzm*z3m;
  Shape[21] = c1*r*zzm*z3p;
  Shape[22] = c1*s*zzm*z3m;
  Shape[23] = c1*s*zzm*z3p;

  // Lower & upper wedge triangular bubble nodes
  Shape[24] = c3*r*s*t*zm;
  Shape[25] = c3*r*s*t*zp;
}

// z in [-1,1]
void Penta26dShapeFct(double dShape[26][3], double m[3])
{
  double r = m[0];
  double s = m[1];
  double z = m[2];

  double t = 1.-r-s;

  double zzm = 1.-z*z;
  double zm  = 1.-z;
  double zp  = 1.+z;
  double z3m = 1.-3.*z;
  double z3p = 1.+3.*z;

  double c0  = 2./9.;
  double c1  = 9./16.;
  double c2  = 9./4.;
  double c3  = 27./2.;

  // Lower & upper wedge's corners   
  dShape[ 0][0] = -c1*zm*(4.*(t*(t-1.)+c0) + 4.*t*(2.*t-1.) - zzm);
  dShape[ 0][1] = -c1*zm*(4.*(t*(t-1.)+c0) + 4.*t*(2.*t-1.) - zzm);
  dShape[ 0][2] = -c1*t*(4.*(t*(t-1.)+c0) - zzm - 2.*zm*z);

  dShape[ 1][0] =  c1*zm*(4.*(r*(r-1.)+c0) + 4.*r*(2.*r-1) - zzm);
  dShape[ 1][1] =  0.0;
  dShape[ 1][2] = -c1*r*(4.*(r*(r-1.)+c0) - zzm - 2.*zm*z);

  dShape[ 2][0] =  0.0;
  dShape[ 2][1] =  c1*zm*(4.*(s*(s-1.)+c0) + 4.*s*(2.*s-1.) - zzm);
  dShape[ 2][2] = -c1*s*(4.*(s*(s-1.)+c0) - zzm - 2.*zm*z);

  dShape[ 3][0] = -c1*zp*(4.*(t*(t-1.)+c0) + 4.*t*(2.*t-1.) - zzm);
  dShape[ 3][1] = -c1*zp*(4.*(t*(t-1.)+c0) + 4.*t*(2.*t-1.) - zzm);
  dShape[ 3][2] =  c1*t*(4.*(t*(t-1.)+c0) - zzm  + 2.*zp*z);

  dShape[ 4][0] =  c1*zp*(4.*(r*(r-1.)+c0) + 4.*r*(2.*r-1.) - zzm);
  dShape[ 4][1] =  0.0;
  dShape[ 4][2] =  c1*r*(4.*(r*(r-1.)+c0) - zzm + 2.*zp*z);

  dShape[ 5][0] =  0.0;
  dShape[ 5][1] =  c1*zp*(4.*(s*(s-1.)+c0) + 4.*s*(2.*s-1.) - zzm);
  dShape[ 5][2] =  c1*s*(4.*(s*(s-1.)+c0) - zzm + 2.*zp*z);

  // Lower & upper wedge's triangular edges' nodes  
  dShape[ 6][0] =  c2*zm*(t*(3.*t-1.) - r*(6.*t-1.));
  dShape[ 6][1] = -c2*zm*r*(6.*t-1.);
  dShape[ 6][2] = -c2*r*t*(3.*t-1.);

  dShape[ 7][0] =  c2*zm*(-r*(3.*r-1.) + t*(6.*r-1.));
  dShape[ 7][1] = -c2*zm*r*(3.*r-1.);
  dShape[ 7][2] = -c2*r*t*(3.*r-1.);

  dShape[ 8][0] =  c2*zm*s*(6.*r-1.);
  dShape[ 8][1] =  c2*zm*r*(3.*r-1.);
  dShape[ 8][2] = -c2*r*s*(3.*r-1.);

  dShape[ 9][0] =  c2*zm*s*(3.*s-1.);
  dShape[ 9][1] =  c2*zm*r*(6.*s-1.);
  dShape[ 9][2] = -c2*r*s*(3.*s-1.);

  dShape[10][0] = -c2*zm*s*(3.*s-1.);
  dShape[10][1] =  c2*zm*(-s*(3.*s-1.) + t*(6.*s-1.)); // PJSA 6/8/2014
  dShape[10][2] = -c2*t*s*(3.*s-1.);

  dShape[11][0] = -c2*zm*s*(6.*t-1.);
  dShape[11][1] =  c2*zm*(t*(3.*t-1.) - s*(6.*t-1.)); // PJSA 6/8/2014
  dShape[11][2] = -c2*t*s*(3.*t-1.);

  dShape[12][0] =  c2*zp*(t*(3.*t-1.) - r*(6.*t-1.));
  dShape[12][1] = -c2*zp*r*(6.*t-1.);
  dShape[12][2] =  c2*r*t*(3.*t-1.);

  dShape[13][0] =  c2*zp*(-r*(3.*r-1.) + t*(6.*r-1.));
  dShape[13][1] = -c2*zp*r*(3.*r-1.);
  dShape[13][2] =  c2*r*t*(3.*r-1.);

  dShape[14][0] =  c2*zp*s*(6.*r-1.);
  dShape[14][1] =  c2*zp*r*(3.*r-1.);
  dShape[14][2] =  c2*r*s*(3.*r-1.);

  dShape[15][0] =  c2*zp*s*(3.*s-1.);
  dShape[15][1] =  c2*zp*r*(6.*s-1.);
  dShape[15][2] =  c2*r*s*(3.*s-1.);

  dShape[16][0] = -c2*zp*s*(3.*s-1.);
  dShape[16][1] =  c2*zp*(-s*(3.*s-1.) + t*(6.*s-1.));
  dShape[16][2] =  c2*t*s*(3.*s-1.);

  dShape[17][0] = -c2*zp*s*(6.*t-1.);
  dShape[17][1] =  c2*zp*(t*(3.*t-1.) - s*(6.*t-1.));
  dShape[17][2] =  c2*t*s*(3.*t-1.);

  // Wedge vertical edges' nodes
  dShape[18][0] = -c1*zzm*z3m;
  dShape[18][1] = -c1*zzm*z3m;
  dShape[18][2] =  c1*t*(-2.*z*z3m - 3.*zzm);

  dShape[19][0] = -c1*zzm*z3p;
  dShape[19][1] = -c1*zzm*z3p;
  dShape[19][2] =  c1*t*(-2.*z*z3p + 3.*zzm);

  dShape[20][0] =  c1*zzm*z3m;
  dShape[20][1] =  0.0;
  dShape[20][2] =  c1*r*(-2.*z*z3m - 3.*zzm);

  dShape[21][0] =  c1*zzm*z3p;
  dShape[21][1] =  0.0;
  dShape[21][2] =  c1*r*(-2.*z*z3p + 3.*zzm);

  dShape[22][0] =  0.0;
  dShape[22][1] =  c1*zzm*z3m;
  dShape[22][2] =  c1*s*(-2.*z*z3m - 3.*zzm);

  dShape[23][0] =  0.0;
  dShape[23][1] =  c1*zzm*z3p;
  dShape[23][2] =  c1*s*(-2.*z*z3p + 3.*zzm);

  // Lower & upper wedge triangular bubble nodes      
  dShape[24][0] =  c3*s*zm*(t-r);
  dShape[24][1] =  c3*r*zm*(t-s);
  dShape[24][2] = -c3*r*s*t;

  dShape[25][0] =  c3*s*zp*(t-r);
  dShape[25][1] =  c3*r*zp*(t-s);
  dShape[25][2] =  c3*r*s*t;
}

// HB (05-21-05): shape fct & derivatives for the 26 nodes wedge element (Iso-parametric formulation)
void 
Penta26ShapeFct(double Shape[26], double dShape[26][3], double m[3])
{
  Penta26ShapeFct(Shape,m);
  Penta26dShapeFct(dShape,m);
}

// HB (05-21-05): shape fct & derivatives for the 26 nodes wedge element (Iso-parametric formulation)
double
Penta26ShapeFct(double Shape[26], double DShape[26][3], double m[3], double X[26], double Y[26], double Z[26])
{
  double dShape[26][3];
  Penta26ShapeFct(Shape, dShape, m);

  return(computePenta26DShapeFct(dShape, X, Y, Z, DShape));
}

