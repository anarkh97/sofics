// ------------------------------------------------------------
// HB - 06-19-05 
// ------------------------------------------------------------

// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>

// HB (06-19-05): for given reference shape fct derivatives dShape & nodes coordinates (X.Y,Z), 
//                compute the jacobian & shape fcts derivative w.r.t to the global coordinate system.
//                assuming an iso-parametric formulation.
double
computeTet10DShapeFct(double dShape[10][3], double X[10], double Y[10], double Z[10], double (*DShape)[3]=0)
{
  double xd1,xd2,xd3;
  double yd1,yd2,yd3;
  double zd1,zd2,zd3;

  xd1 = 0.0; xd2 = 0.0; xd3 = 0.0;
  yd1 = 0.0; yd2 = 0.0; yd3 = 0.0;
  zd1 = 0.0; zd2 = 0.0; zd3 = 0.0;

  for(int i=0; i<10; i+=5){
    xd1 += dShape[i][0]*X[i] + dShape[i+1][0]*X[i+1] + dShape[i+2][0]*X[i+2] + dShape[i+3][0]*X[i+3] + dShape[i+4][0]*X[i+4];
    xd2 += dShape[i][1]*X[i] + dShape[i+1][1]*X[i+1] + dShape[i+2][1]*X[i+2] + dShape[i+3][1]*X[i+3] + dShape[i+4][1]*X[i+4];
    xd3 += dShape[i][2]*X[i] + dShape[i+1][2]*X[i+1] + dShape[i+2][2]*X[i+2] + dShape[i+3][2]*X[i+3] + dShape[i+4][2]*X[i+4];
    yd1 += dShape[i][0]*Y[i] + dShape[i+1][0]*Y[i+1] + dShape[i+2][0]*Y[i+2] + dShape[i+3][0]*Y[i+3] + dShape[i+4][0]*Y[i+4];
    yd2 += dShape[i][1]*Y[i] + dShape[i+1][1]*Y[i+1] + dShape[i+2][1]*Y[i+2] + dShape[i+3][1]*Y[i+3] + dShape[i+4][1]*Y[i+4];
    yd3 += dShape[i][2]*Y[i] + dShape[i+1][2]*Y[i+1] + dShape[i+2][2]*Y[i+2] + dShape[i+3][2]*Y[i+3] + dShape[i+4][2]*Y[i+4];
    zd1 += dShape[i][0]*Z[i] + dShape[i+1][0]*Z[i+1] + dShape[i+2][0]*Z[i+2] + dShape[i+3][0]*Z[i+3] + dShape[i+4][0]*Z[i+4];
    zd2 += dShape[i][1]*Z[i] + dShape[i+1][1]*Z[i+1] + dShape[i+2][1]*Z[i+2] + dShape[i+3][1]*Z[i+3] + dShape[i+4][1]*Z[i+4];
    zd3 += dShape[i][2]*Z[i] + dShape[i+1][2]*Z[i+1] + dShape[i+2][2]*Z[i+2] + dShape[i+3][2]*Z[i+3] + dShape[i+4][2]*Z[i+4];
  }

  double a11 = yd2*zd3 - yd3*zd2; double a12 = yd3*zd1 - yd1*zd3; double a13 = yd1*zd2 - yd2*zd1;
  double a21 = xd3*zd2 - xd2*zd3; double a22 = xd1*zd3 - zd1*xd3; double a23 = xd2*zd1 - xd1*zd2;
  double a31 = xd2*yd3 - xd3*yd2; double a32 = yd1*xd3 - xd1*yd3; double a33 = xd1*yd2 - yd1*xd2;

  /* -------> DETERMINANT OF THE JACOBIAN MATRIX <--- */

  double J = xd1*a11 + yd1*a21 + zd1*a31 ;
  if(J == 0.0) {
   fprintf(stderr," *** WARNING: NULL JACOBIAN IN computeTet10DShapeFct ROUTINE.\n");
  }
  if(J < 0.0) {
    fprintf(stderr," *** WARNING: NEGATIVE JACOBIAN IN computeTet10DShapeFct ROUTINE.\n");
  }
  double cdet = 1.0/J;

  /* -------> DERIVATIVE OF THE SHAPE FCTS IN THE "REAL" ELEMENT <--- */
  if(DShape) {
    DShape[0][0] = cdet * ( a11*dShape[0][0] + a12*dShape[0][1] + a13*dShape[0][2] );
    DShape[1][0] = cdet * ( a11*dShape[1][0] + a12*dShape[1][1] + a13*dShape[1][2] );
    DShape[2][0] = cdet * ( a11*dShape[2][0] + a12*dShape[2][1] + a13*dShape[2][2] );
    DShape[3][0] = cdet * ( a11*dShape[3][0] + a12*dShape[3][1] + a13*dShape[3][2] );
    DShape[4][0] = cdet * ( a11*dShape[4][0] + a12*dShape[4][1] + a13*dShape[4][2] );
    DShape[5][0] = cdet * ( a11*dShape[5][0] + a12*dShape[5][1] + a13*dShape[5][2] );
    DShape[6][0] = cdet * ( a11*dShape[6][0] + a12*dShape[6][1] + a13*dShape[6][2] );
    DShape[7][0] = cdet * ( a11*dShape[7][0] + a12*dShape[7][1] + a13*dShape[7][2] );
    DShape[8][0] = cdet * ( a11*dShape[8][0] + a12*dShape[8][1] + a13*dShape[8][2] );
    DShape[9][0] = cdet * ( a11*dShape[9][0] + a12*dShape[9][1] + a13*dShape[9][2] );

    DShape[0][1] = cdet * ( a21*dShape[0][0] + a22*dShape[0][1] + a23*dShape[0][2] );
    DShape[1][1] = cdet * ( a21*dShape[1][0] + a22*dShape[1][1] + a23*dShape[1][2] );
    DShape[2][1] = cdet * ( a21*dShape[2][0] + a22*dShape[2][1] + a23*dShape[2][2] );
    DShape[3][1] = cdet * ( a21*dShape[3][0] + a22*dShape[3][1] + a23*dShape[3][2] );
    DShape[4][1] = cdet * ( a21*dShape[4][0] + a22*dShape[4][1] + a23*dShape[4][2] );
    DShape[5][1] = cdet * ( a21*dShape[5][0] + a22*dShape[5][1] + a23*dShape[5][2] );
    DShape[6][1] = cdet * ( a21*dShape[6][0] + a22*dShape[6][1] + a23*dShape[6][2] );
    DShape[7][1] = cdet * ( a21*dShape[7][0] + a22*dShape[7][1] + a23*dShape[7][2] );
    DShape[8][1] = cdet * ( a21*dShape[8][0] + a22*dShape[8][1] + a23*dShape[8][2] );
    DShape[9][1] = cdet * ( a21*dShape[9][0] + a22*dShape[9][1] + a23*dShape[9][2] );

    DShape[0][2] = cdet * ( a31*dShape[0][0] + a32*dShape[0][1] + a33*dShape[0][2] );
    DShape[1][2] = cdet * ( a31*dShape[1][0] + a32*dShape[1][1] + a33*dShape[1][2] );
    DShape[2][2] = cdet * ( a31*dShape[2][0] + a32*dShape[2][1] + a33*dShape[2][2] );
    DShape[3][2] = cdet * ( a31*dShape[3][0] + a32*dShape[3][1] + a33*dShape[3][2] );
    DShape[4][2] = cdet * ( a31*dShape[4][0] + a32*dShape[4][1] + a33*dShape[4][2] );
    DShape[5][2] = cdet * ( a31*dShape[5][0] + a32*dShape[5][1] + a33*dShape[5][2] );
    DShape[6][2] = cdet * ( a31*dShape[6][0] + a32*dShape[6][1] + a33*dShape[6][2] );
    DShape[7][2] = cdet * ( a31*dShape[7][0] + a32*dShape[7][1] + a33*dShape[7][2] );
    DShape[8][2] = cdet * ( a31*dShape[8][0] + a32*dShape[8][1] + a33*dShape[8][2] );
    DShape[9][2] = cdet * ( a31*dShape[9][0] + a32*dShape[9][1] + a33*dShape[9][2] );
  }
  
  return(J);
}

// HB (06-19-05): shape fct & derivatives for the 10 nodes tetrahedral element (Iso-parametric formulation)
void 
Tetra10ShapeFct(double Shape[10], double dShape[10][3], double m[3])
{
  double r = m[0]; // = x
  double s = m[1]; // = y
  double t = m[2]; // = z
  double u = 1.-r-s-t;

  // -------> VALUE OF THE SHAPE FCTS IN THE "REFERENCE" ELEMENT

  Shape[0] = u*(2.*u-1.);
  Shape[1] = r*(2.*r-1.);
  Shape[2] = s*(2.*s-1.);
  Shape[3] = t*(2.*t-1.);

  Shape[4] = 4.*r*u;
  Shape[5] = 4.*r*s;
  Shape[6] = 4.*s*u;

  Shape[7] = 4.*t*u;
  Shape[8] = 4.*t*r;
  Shape[9] = 4.*t*s;

  // -------> DERIVATIVE OF THE SHAPE FCTS IN THE "REFERENCE" ELEMENT

  dShape[0][0] = 1.-4.*u; dShape[0][1] = 1.-4.*u; dShape[0][2] = 1.-4.*u;
  dShape[1][0] = 4.*r-1.; dShape[1][1] =   0.0  ; dShape[1][2] =   0.0  ;
  dShape[2][0] =   0.0  ; dShape[2][1] = 4.*s-1.; dShape[2][2] =   0.0  ;
  dShape[3][0] =   0.0  ; dShape[3][1] =   0.0  ; dShape[3][2] = 4.*t-1.;

  dShape[4][0] =4.*(u-r); dShape[4][1] = -4.*r  ; dShape[4][2] = -4.*r  ;
  dShape[5][0] =  4.*s  ; dShape[5][1] =  4.*r  ; dShape[5][2] =  0.0   ;
  dShape[6][0] = -4.*s  ; dShape[6][1] =4.*(u-s); dShape[6][2] = -4.*s  ;

  dShape[7][0] = -4.*t  ; dShape[7][1] = -4.*t  ; dShape[7][2] =  4.*(u-t); // PJSA 6/8/2014
  dShape[8][0] =  4.*t  ; dShape[8][1] =  0.0   ; dShape[8][2] =  4.*r  ;
  dShape[9][0] =  0.0   ; dShape[9][1] =  4.*t  ; dShape[9][2] =  4.*s  ;
}

// HB (06-19-05): shape fct & derivatives for the 10 nodes tetrahedral element (Iso-parametric formulation)
double
Tetra10ShapeFct(double Shape[10], double DShape[10][3], double m[3], double X[10], double Y[10], double Z[10])
{
  double dShape[10][3];
  Tetra10ShapeFct(Shape, dShape, m);

  return(computeTet10DShapeFct(dShape, X, Y, Z, DShape));
}

