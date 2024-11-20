// ------------------------------------------------------------
// HB - 09-21-03 
// ------------------------------------------------------------
// !!! HERE z in [-1,1] !!!
// ------------------------------------------------------------

// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>

double
computePenta6DShapeFct(double dShape[6][3], double X[6], double Y[6], double Z[6], double (*DShape)[3]=0)
{ 
  double xd1,xd2,xd3;
  double yd1,yd2,yd3;
  double zd1,zd2,zd3;

  xd1 = 0.0; xd2 = 0.0; xd3 = 0.0;
  yd1 = 0.0; yd2 = 0.0; yd3 = 0.0;
  zd1 = 0.0; zd2 = 0.0; zd3 = 0.0;

  for(int i=0; i<6; i+=3){
    xd1 += dShape[i][0]*X[i] + dShape[i+1][0]*X[i+1] + dShape[i+2][0]*X[i+2];
    xd2 += dShape[i][1]*X[i] + dShape[i+1][1]*X[i+1] + dShape[i+2][1]*X[i+2];
    xd3 += dShape[i][2]*X[i] + dShape[i+1][2]*X[i+1] + dShape[i+2][2]*X[i+2];
    
    yd1 += dShape[i][0]*Y[i] + dShape[i+1][0]*Y[i+1] + dShape[i+2][0]*Y[i+2];
    yd2 += dShape[i][1]*Y[i] + dShape[i+1][1]*Y[i+1] + dShape[i+2][1]*Y[i+2];
    yd3 += dShape[i][2]*Y[i] + dShape[i+1][2]*Y[i+1] + dShape[i+2][2]*Y[i+2];
    
    zd1 += dShape[i][0]*Z[i] + dShape[i+1][0]*Z[i+1] + dShape[i+2][0]*Z[i+2];
    zd2 += dShape[i][1]*Z[i] + dShape[i+1][1]*Z[i+1] + dShape[i+2][1]*Z[i+2];
    zd3 += dShape[i][2]*Z[i] + dShape[i+1][2]*Z[i+1] + dShape[i+2][2]*Z[i+2];
  }

  double a11 = yd2*zd3 - yd3*zd2; double a12 = yd3*zd1 - yd1*zd3; double a13 = yd1*zd2 - yd2*zd1;
  double a21 = xd3*zd2 - xd2*zd3; double a22 = xd1*zd3 - zd1*xd3; double a23 = xd2*zd1 - xd1*zd2;
  double a31 = xd2*yd3 - xd3*yd2; double a32 = yd1*xd3 - xd1*yd3; double a33 = xd1*yd2 - yd1*xd2;

   /* -------> DETERMINANT OF THE JACOBIAN MATRIX <--- */
  double J = xd1*a11 + yd1*a21 + zd1*a31;  
  if(J == 0.0) {
    fprintf(stderr," *** WARNING: NULL JACOBIAN IN computePenta6DShapeFct ROUTINE.\n"); 
  }
  if(J < 0.0) { 
    fprintf(stderr," *** WARNING: NEGATIVE JACOBIAN IN computePenta6DShapeFct ROUTINE.\n"); 
  }
  double cdet = 1.0/J;

  /* -------> DERIVATIVE OF THE SHAPE FCTS IN THE "REAL" ELEMENT <--- */
  DShape[0][0] = cdet * ( a11*dShape[0][0] + a12*dShape[0][1] + a13*dShape[0][2] );
  DShape[1][0] = cdet * ( a11*dShape[1][0] + a12*dShape[1][1] + a13*dShape[1][2] );
  DShape[2][0] = cdet * ( a11*dShape[2][0] + a12*dShape[2][1] + a13*dShape[2][2] );
  DShape[3][0] = cdet * ( a11*dShape[3][0] + a12*dShape[3][1] + a13*dShape[3][2] );
  DShape[4][0] = cdet * ( a11*dShape[4][0] + a12*dShape[4][1] + a13*dShape[4][2] );
  DShape[5][0] = cdet * ( a11*dShape[5][0] + a12*dShape[5][1] + a13*dShape[5][2] );

  DShape[0][1] = cdet * ( a21*dShape[0][0] + a22*dShape[0][1] + a23*dShape[0][2] );
  DShape[1][1] = cdet * ( a21*dShape[1][0] + a22*dShape[1][1] + a23*dShape[1][2] );
  DShape[2][1] = cdet * ( a21*dShape[2][0] + a22*dShape[2][1] + a23*dShape[2][2] );
  DShape[3][1] = cdet * ( a21*dShape[3][0] + a22*dShape[3][1] + a23*dShape[3][2] );
  DShape[4][1] = cdet * ( a21*dShape[4][0] + a22*dShape[4][1] + a23*dShape[4][2] );
  DShape[5][1] = cdet * ( a21*dShape[5][0] + a22*dShape[5][1] + a23*dShape[5][2] );

  DShape[0][2] = cdet * ( a31*dShape[0][0] + a32*dShape[0][1] + a33*dShape[0][2] );
  DShape[1][2] = cdet * ( a31*dShape[1][0] + a32*dShape[1][1] + a33*dShape[1][2] );
  DShape[2][2] = cdet * ( a31*dShape[2][0] + a32*dShape[2][1] + a33*dShape[2][2] );
  DShape[3][2] = cdet * ( a31*dShape[3][0] + a32*dShape[3][1] + a33*dShape[3][2] );
  DShape[4][2] = cdet * ( a31*dShape[4][0] + a32*dShape[4][1] + a33*dShape[4][2] );
  DShape[5][2] = cdet * ( a31*dShape[5][0] + a32*dShape[5][1] + a33*dShape[5][2] );
  
  return(J);
}

// HB (09-23-03): shape fct & derivative for the 6 nodes pentahedral element
void
Penta6ShapeFct(double Shape[6], double dShape[6][3], double m[3])
{
  double r = m[0]; // = x
  double s = m[1]; // = y
  double z = m[2];
  double t = 1.-r-s;
  double zp = 0.5*(1.0+z); // z in [-1,1]
  double zm = 0.5*(1.0-z);

  // -------> VALUE OF THE SHAPE FCTS IN THE "REFERENCE" ELEMENT <--- 
  
  Shape[0] =  t*zm;
  Shape[1] =  r*zm;
  Shape[2] =  s*zm;
  Shape[3] =  t*zp;
  Shape[4] =  r*zp;
  Shape[5] =  s*zp;

  // -------> DERIVATIVE OF THE SHAPE FCTS IN THE "REFERENCE" ELEMENT <--- 
  
  dShape[0][0] = -zm ; dShape[0][1] = -zm ; dShape[0][2] = -0.5*t;
  dShape[1][0] =  zm ; dShape[1][1] =  0.0; dShape[1][2] = -0.5*r; 
  dShape[2][0] =  0.0; dShape[2][1] =  zm ; dShape[2][2] = -0.5*s;
  dShape[3][0] = -zp ; dShape[3][1] = -zp ; dShape[3][2] =  0.5*t;
  dShape[4][0] =  zp ; dShape[4][1] =  0.0; dShape[4][2] =  0.5*r;
  dShape[5][0] =  0.0; dShape[5][1] =  zp ; dShape[5][2] =  0.5*s;
}

double
Penta6ShapeFct(double Shape[6], double DShape[6][3], double m[3], double X[6], double Y[6], double Z[6])
{ 
  double dShape[6][3];
  Penta6ShapeFct(Shape, dShape, m);
#ifndef PENT6_SHAPE_OPT
  return(computePenta6DShapeFct(dShape, X, Y, Z, DShape));
#else
  fprintf(stderr," *** ERROR: Penta6ShapeFct.C: optimization NOT implemented. Abort.\n");
  exit(-1);
  double r = m[0]; // = x
  double s = m[1]; // = y
  double z = m[2];
  double t = 1.-r-s;
  double zp = z;
  double zm = 1.0-z;
  double A = zm;
  double B = zp;
  double J, cdet;

  double xd1,xd2,xd3;
  double yd1,yd2,yd3;
  double zd1,zd2,zd3;

  xd1 = 0.0; xd2 = 0.0; xd3 = 0.0;
  yd1 = 0.0; yd2 = 0.0; yd3 = 0.0;
  zd1 = 0.0; zd2 = 0.0; zd3 = 0.0;
 
  xd1 = (X[1]-X[0])*Zm + (X[4]-X[3])*Zp;
  yd1 = (Y[1]-Y[0])*Zm + (Y[4]-Y[3])*Zp;
  zd1 = (Z[1]-Z[0])*Zm + (Z[4]-Z[3])*Zp;

  xd2 = (X[2]-X[0])*Zm + (X[5]-X[3])*Zp;
  yd2 = (Y[2]-Y[0])*Zm + (Y[5]-Y[3])*Zp;
  zd2 = (Z[2]-Z[0])*Zm + (Z[5]-Z[3])*Zp;

  xd3 = (X[3]-X[0])*t + (X[4]-X[1])*r + (X[5]-X[2])*s; xd3 *= 0.5;
  yd3 = (Y[3]-Y[0])*t + (Y[4]-Y[1])*r + (Y[5]-Y[2])*s; yd3 *= 0.5;
  zd3 = (Z[3]-Z[0])*t + (Z[4]-Z[1])*r + (Z[5]-Z[2])*s; zd3 *= 0.5;
 
  double a11 = yd2*zd3 - yd3*zd2; double a12 = yd3*zd1 - yd1*zd3; double a13 = yd1*zd2 - yd2*zd1;
  double a21 = xd3*zd2 - xd2*zd3; double a22 = xd1*zd3 - zd1*xd3; double a23 = xd2*zd1 - xd1*zd2;
  double a31 = xd2*yd3 - xd3*yd2; double a32 = yd1*xd3 - xd1*yd3; double a33 = xd1*yd2 - yd1*xd2;
  
  // -------> DETERMINANT OF THE JACOBIAN MATRIX <--- 

  double J = xd1*a11 + yd1*a21 + zd1*a31;
  if(J == 0.0) { 
   fprintf(stderr," *** WARNING: NULL JACOBIAN IN Penta6ShapeFct.C ROUTINE.\n"); 
  }
  if(J < 0.0) { 
    fprintf(stderr," *** WARNING: NEGATIVE JACOBIAN IN Penta6ShapeFct.C ROUTINE.\n"); 
  }
  cdet = 1.0/J;

  // -------> DERIVATIVE OF THE SHAPE FCTS IN THE "REAL" ELEMENT <--- 
  
  DShape[0][0] = cdet * ( a11*dShape[0][0] + a12*dShape[0][1] + a13*dShape[0][2] );
  DShape[1][0] = cdet * ( a11*dShape[1][0] + a12*dShape[1][1] + a13*dShape[1][2] );
  DShape[2][0] = cdet * ( a11*dShape[2][0] + a12*dShape[2][1] + a13*dShape[2][2] );
  DShape[3][0] = cdet * ( a11*dShape[3][0] + a12*dShape[3][1] + a13*dShape[3][2] );
  DShape[4][0] = cdet * ( a11*dShape[4][0] + a12*dShape[4][1] + a13*dShape[4][2] );
  DShape[5][0] = cdet * ( a11*dShape[5][0] + a12*dShape[5][1] + a13*dShape[5][2] );

  DShape[0][1] = cdet * ( a21*dShape[0][0] + a22*dShape[0][1] + a23*dShape[0][2] );
  DShape[1][1] = cdet * ( a21*dShape[1][0] + a22*dShape[1][1] + a23*dShape[1][2] );
  DShape[2][1] = cdet * ( a21*dShape[2][0] + a22*dShape[2][1] + a23*dShape[2][2] );
  DShape[3][1] = cdet * ( a21*dShape[3][0] + a22*dShape[3][1] + a23*dShape[3][2] );
  DShape[4][1] = cdet * ( a21*dShape[4][0] + a22*dShape[4][1] + a23*dShape[4][2] );
  DShape[5][1] = cdet * ( a21*dShape[5][0] + a22*dShape[5][1] + a23*dShape[5][2] );

  DShape[0][2] = cdet * ( a31*dShape[0][0] + a32*dShape[0][1] + a33*dShape[0][2] );
  DShape[1][2] = cdet * ( a31*dShape[1][0] + a32*dShape[1][1] + a33*dShape[1][2] );
  DShape[2][2] = cdet * ( a31*dShape[2][0] + a32*dShape[2][1] + a33*dShape[2][2] );
  DShape[3][2] = cdet * ( a31*dShape[3][0] + a32*dShape[3][1] + a33*dShape[3][2] );
  DShape[4][2] = cdet * ( a31*dShape[4][0] + a32*dShape[4][1] + a33*dShape[4][2] );
  DShape[5][2] = cdet * ( a31*dShape[5][0] + a32*dShape[5][1] + a33*dShape[5][2] );
  
  return(J);
#endif
}

