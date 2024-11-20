// ------------------------------------------------------------
// HB - 04-15-05 
// ------------------------------------------------------------

// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>

double
computeTetra4Jacobian(double X[4], double Y[4], double Z[4], double a[3][3])
{
  double xd1,xd2,xd3;
  double yd1,yd2,yd3;
  double zd1,zd2,zd3;

  xd1 = (X[1]-X[0]); xd2 = (X[2]-X[0]);  xd3 = (X[3]-X[0]);
  yd1 = (Y[1]-Y[0]); yd2 = (Y[2]-Y[0]);  yd3 = (Y[3]-Y[0]);
  zd1 = (Z[1]-Z[0]); zd2 = (Z[2]-Z[0]);  zd3 = (Z[3]-Z[0]);

  a[0][0] = yd2*zd3 - yd3*zd2; a[0][1] = yd3*zd1 - yd1*zd3; a[0][2] = yd1*zd2 - yd2*zd1;
  a[1][0] = xd3*zd2 - xd2*zd3; a[1][1] = xd1*zd3 - zd1*xd3; a[1][2] = xd2*zd1 - xd1*zd2;
  a[2][0] = xd2*yd3 - xd3*yd2; a[2][1] = yd1*xd3 - xd1*yd3; a[2][2] = xd1*yd2 - yd1*xd2;
   // -------> DETERMINANT OF THE JACOBIAN MATRIX <--- 

  double J = xd1*a[0][0] + yd1*a[1][0] + zd1*a[2][0];
  if(J == 0.0) { 
    fprintf(stderr," *** WARNING: NULL JACOBIAN IN computeTetra4Jacobian ROUTINE.\n"); 
  }
  if(J < 0.0) {
    fprintf(stderr," *** WARNING: NEGATIVE JACOBIAN IN computeTetra4Jacobian ROUTINE.\n"); 
  }
  return J;
}

void
computeTetra4dadx(double X[4], double Y[4], double Z[4], double dadx[3][3][12])
{
  dadx[0][0][1] = Z[2]-Z[3];  dadx[0][0][2] = Y[3]-Y[2];  dadx[0][0][7] = Z[3]-Z[0];  dadx[0][0][0] = 0.0;  dadx[0][0][3] = 0.0;  dadx[0][0][4] = 0.0;
  dadx[0][0][10] = Z[0]-Z[2]; dadx[0][0][8] = Y[0]-Y[3];  dadx[0][0][11]= Y[2]-Y[0];  dadx[0][0][5] = 0.0;  dadx[0][0][6] = 0.0;  dadx[0][0][9] = 0.0;
  dadx[0][1][1] = Z[3]-Z[1];  dadx[0][1][2] = Y[1]-Y[3];  dadx[0][1][4] = Z[0]-Z[3];  dadx[0][1][0] = 0.0;  dadx[0][1][3] = 0.0;  dadx[0][1][6] = 0.0;
  dadx[0][1][10]= Z[1]-Z[0];  dadx[0][1][5] = Y[3]-Y[0];  dadx[0][1][11]= Y[0]-Y[1];  dadx[0][1][7] = 0.0;  dadx[0][1][8] = 0.0;  dadx[0][1][9] = 0.0;
  dadx[0][2][1] = Z[1]-Z[2];  dadx[0][2][2] = Y[2]-Y[1];  dadx[0][2][4] = Z[2]-Z[0];  dadx[0][2][0] = 0.0;  dadx[0][2][3] = 0.0;  dadx[0][2][6] = 0.0;
  dadx[0][2][7] = Z[0]-Z[1];  dadx[0][2][5] = Y[0]-Y[2];  dadx[0][2][8] = Y[1]-Y[0];  dadx[0][2][9] = 0.0;  dadx[0][2][10] = 0.0; dadx[0][2][11] = 0.0;

  dadx[1][0][0] = Z[3]-Z[2];  dadx[1][0][2] = X[2]-X[3];  dadx[1][0][6] = Z[0]-Z[3];  dadx[1][0][1] = 0.0;  dadx[1][0][3] = 0.0;  dadx[1][0][4] = 0.0;
  dadx[1][0][9] = Z[2]-Z[0];  dadx[1][0][8] = X[3]-X[0];  dadx[1][0][11]= X[0]-X[2];  dadx[1][0][5] = 0.0;  dadx[1][0][7] = 0.0;  dadx[1][0][10] = 0.0;
  dadx[1][1][0] = Z[1]-Z[3];  dadx[1][1][2] = X[3]-X[1];  dadx[1][1][3] = Z[3]-Z[0];  dadx[1][1][1] = 0.0;  dadx[1][1][4] = 0.0;  dadx[1][1][6] = 0.0;
  dadx[1][1][9] = Z[0]-Z[1];  dadx[1][1][5] = X[0]-X[3];  dadx[1][1][11]= X[1]-X[0];  dadx[1][1][7] = 0.0;  dadx[1][1][8] = 0.0;  dadx[1][1][10] = 0.0;
  dadx[1][2][0] = Z[2]-Z[1];  dadx[1][2][2] = X[1]-X[2];  dadx[1][2][3] = Z[0]-Z[2];  dadx[1][2][1] = 0.0;  dadx[1][2][4] = 0.0;  dadx[1][2][7] = 0.0;
  dadx[1][2][6] = Z[1]-Z[0];  dadx[1][2][5] = X[2]-X[0];  dadx[1][2][8] = X[0]-X[1];  dadx[1][2][9] = 0.0;  dadx[1][2][10] = 0.0;  dadx[1][2][11] = 0.0;

  dadx[2][0][0] = Y[2]-Y[3];  dadx[2][0][1] = X[3]-X[2];  dadx[2][0][6] = Y[3]-Y[0];  dadx[2][0][2] = 0.0;  dadx[2][0][3] = 0.0;  dadx[2][0][4] = 0.0;
  dadx[2][0][9] = Y[0]-Y[2];  dadx[2][0][7] = X[0]-X[3];  dadx[2][0][10]= X[2]-X[0];  dadx[2][0][5] = 0.0;  dadx[2][0][8] = 0.0;  dadx[2][0][11] = 0.0;
  dadx[2][1][0] = Y[3]-Y[1];  dadx[2][1][1] = X[1]-X[3];  dadx[2][1][3] = Y[0]-Y[3];  dadx[2][1][2] = 0.0;  dadx[2][1][5] = 0.0;  dadx[2][1][6] = 0.0;
  dadx[2][1][9] = Y[1]-Y[0];  dadx[2][1][4] = X[3]-X[0];  dadx[2][1][10]= X[0]-X[1];  dadx[2][1][7] = 0.0;  dadx[2][1][8] = 0.0;  dadx[2][1][11] = 0.0;
  dadx[2][2][0] = Y[1]-Y[2];  dadx[2][2][1] = X[2]-X[1];  dadx[2][2][3] = Y[2]-Y[0];  dadx[2][2][2] = 0.0;  dadx[2][2][5] = 0.0;  dadx[2][2][8] = 0.0;
  dadx[2][2][6] = Y[0]-Y[1];  dadx[2][2][4] = X[0]-X[2];  dadx[2][2][7] = X[1]-X[0];  dadx[2][2][9] = 0.0;  dadx[2][2][10] = 0.0; dadx[2][2][11] = 0.0;
}

void
computeTetra4dadxTimesdShapeFct(double dShape[4][3], double J, double dadx[3][3][12], double DDShape[4][3][12])
{
  double cdet = 1.0/J;

  for(int i=0; i<4; ++i) {

    DDShape[i][0][0] = 0.0;
    DDShape[i][0][1] = cdet * ( dadx[0][0][1]*dShape[i][0] + dadx[0][1][1]*dShape[i][1] + dadx[0][2][1]*dShape[i][2] );
    DDShape[i][0][2] = cdet * ( dadx[0][0][2]*dShape[i][0] + dadx[0][1][2]*dShape[i][1] + dadx[0][2][2]*dShape[i][2] );
    DDShape[i][0][3] = 0.0;
    DDShape[i][0][4] = cdet * ( dadx[0][1][4]*dShape[i][1] + dadx[0][2][4]*dShape[i][2] );
    DDShape[i][0][5] = cdet * ( dadx[0][1][5]*dShape[i][1] + dadx[0][2][5]*dShape[i][2] );
    DDShape[i][0][6] = 0.0;
    DDShape[i][0][7] = cdet * ( dadx[0][0][7]*dShape[i][0] + dadx[0][2][7]*dShape[i][2] );
    DDShape[i][0][8] = cdet * ( dadx[0][0][8]*dShape[i][0] + dadx[0][2][8]*dShape[i][2] );
    DDShape[i][0][9] = 0.0;
    DDShape[i][0][10]= cdet * ( dadx[0][0][10]*dShape[i][0] + dadx[0][1][10]*dShape[i][1] );
    DDShape[i][0][11]= cdet * ( dadx[0][0][11]*dShape[i][0] + dadx[0][1][11]*dShape[i][1] );
    
    DDShape[i][1][0] = cdet * ( dadx[1][0][0]*dShape[i][0] + dadx[1][1][0]*dShape[i][1] + dadx[1][2][0]*dShape[i][2] );
    DDShape[i][1][1] = 0.0;
    DDShape[i][1][2] = cdet * ( dadx[1][0][2]*dShape[i][0] + dadx[1][1][2]*dShape[i][1] + dadx[1][2][2]*dShape[i][2] );
    DDShape[i][1][3] = cdet * ( dadx[1][1][3]*dShape[i][1] + dadx[1][2][3]*dShape[i][2] );
    DDShape[i][1][4] = 0.0;
    DDShape[i][1][5] = cdet * ( dadx[1][1][5]*dShape[i][1] + dadx[1][2][5]*dShape[i][2] );
    DDShape[i][1][6] = cdet * ( dadx[1][0][6]*dShape[i][0] + dadx[1][2][6]*dShape[i][2] );
    DDShape[i][1][7] = 0.0;
    DDShape[i][1][8] = cdet * ( dadx[1][0][8]*dShape[i][0] + dadx[1][2][8]*dShape[i][2] );
    DDShape[i][1][9] = cdet * ( dadx[1][0][9]*dShape[i][0] + dadx[1][1][9]*dShape[i][1] );
    DDShape[i][1][10] = 0.0;
    DDShape[i][1][11] = cdet * ( dadx[1][0][11]*dShape[i][0] + dadx[1][1][11]*dShape[i][1] );

    DDShape[i][2][0] = cdet * ( dadx[2][0][0]*dShape[i][0] + dadx[2][1][0]*dShape[i][1] + dadx[2][2][0]*dShape[i][2] );
    DDShape[i][2][1] = cdet * ( dadx[2][0][1]*dShape[i][0] + dadx[2][1][1]*dShape[i][1] + dadx[2][2][1]*dShape[i][2] );
    DDShape[i][1][2] = 0.0;
    DDShape[i][2][3] = cdet * ( dadx[2][1][3]*dShape[i][1] + dadx[2][2][3]*dShape[i][2] );
    DDShape[i][2][4] = cdet * ( dadx[2][1][4]*dShape[i][1] + dadx[2][2][4]*dShape[i][2] );
    DDShape[i][1][5] = 0.0;
    DDShape[i][2][6] = cdet * ( dadx[2][0][6]*dShape[i][0] + dadx[2][2][6]*dShape[i][2] );
    DDShape[i][2][7] = cdet * ( dadx[2][0][7]*dShape[i][0] + dadx[2][2][7]*dShape[i][2] );
    DDShape[i][1][8] = 0.0;
    DDShape[i][2][9] = cdet * ( dadx[2][0][9]*dShape[i][0] + dadx[2][1][9]*dShape[i][1] );
    DDShape[i][2][10]= cdet * ( dadx[2][0][10]*dShape[i][0] + dadx[2][1][10]*dShape[i][1] );
    DDShape[i][1][11] = 0.0;

  }

}

void
computeTetra4DShapeFct(double dShape[4][3], double J, double a[3][3], double (*DShape)[3]=0)
{
  double cdet = 1.0/J;

  // -------> DERIVATIVE OF THE SHAPE FCTS IN THE "REAL" ELEMENT <--- 
 
  DShape[0][0] = cdet * ( a[0][0]*dShape[0][0] + a[0][1]*dShape[0][1] + a[0][2]*dShape[0][2] );
  DShape[1][0] = cdet * ( a[0][0]*dShape[1][0] + a[0][1]*dShape[1][1] + a[0][2]*dShape[1][2] );
  DShape[2][0] = cdet * ( a[0][0]*dShape[2][0] + a[0][1]*dShape[2][1] + a[0][2]*dShape[2][2] );
  DShape[3][0] = cdet * ( a[0][0]*dShape[3][0] + a[0][1]*dShape[3][1] + a[0][2]*dShape[3][2] );

  DShape[0][1] = cdet * ( a[1][0]*dShape[0][0] + a[1][1]*dShape[0][1] + a[1][2]*dShape[0][2] );
  DShape[1][1] = cdet * ( a[1][0]*dShape[1][0] + a[1][1]*dShape[1][1] + a[1][2]*dShape[1][2] );
  DShape[2][1] = cdet * ( a[1][0]*dShape[2][0] + a[1][1]*dShape[2][1] + a[1][2]*dShape[2][2] );
  DShape[3][1] = cdet * ( a[1][0]*dShape[3][0] + a[1][1]*dShape[3][1] + a[1][2]*dShape[3][2] );

  DShape[0][2] = cdet * ( a[2][0]*dShape[0][0] + a[2][1]*dShape[0][1] + a[2][2]*dShape[0][2] );
  DShape[1][2] = cdet * ( a[2][0]*dShape[1][0] + a[2][1]*dShape[1][1] + a[2][2]*dShape[1][2] );
  DShape[2][2] = cdet * ( a[2][0]*dShape[2][0] + a[2][1]*dShape[2][1] + a[2][2]*dShape[2][2] );
  DShape[3][2] = cdet * ( a[2][0]*dShape[3][0] + a[2][1]*dShape[3][1] + a[2][2]*dShape[3][2] ); 
}

// HB (04-15-05): shape fct & derivative for the 4 nodes tetrahedral element (Iso-parametric formulation)
// This implementation can be a bit more optimized ...
void
Tetra4ShapeFct(double Shape[4], double dShape[4][3], double m[3])
{
  double r = m[0]; // = x
  double s = m[1]; // = y
  double t = m[2]; // = z
  double u = 1.-r-s-t;

  /* -------> VALUE OF THE SHAPE FCTS IN THE "REFERENCE" ELEMENT <--- */
 
  Shape[0] = u;
  Shape[1] = r;
  Shape[2] = s;
  Shape[3] = t;

  /* -------> DERIVATIVE OF THE SHAPE FCTS IN THE "REFERENCE" ELEMENT <--- */
 
  dShape[0][0] = -1.0; dShape[0][1] = -1.0; dShape[0][2] = -1.0;
  dShape[1][0] =  1.0; dShape[1][1] =  0.0; dShape[1][2] =  0.0;
  dShape[2][0] =  0.0; dShape[2][1] =  1.0; dShape[2][2] =  0.0;
  dShape[3][0] =  0.0; dShape[3][1] =  0.0; dShape[3][2] =  1.0;
}

void
Tetra4ShapeFct(double Shape[4], double DShape[4][3], double m[3], double J, double a[3][3])
{
  double dShape[4][3];
  Tetra4ShapeFct(Shape, dShape, m);
  computeTetra4DShapeFct(dShape, J, a, DShape);
}

void computedShape(double dShape[4][3])
{
  dShape[0][0] = -1.0; dShape[0][1] = -1.0; dShape[0][2] = -1.0;
  dShape[1][0] =  1.0; dShape[1][1] =  0.0; dShape[1][2] =  0.0;
  dShape[2][0] =  0.0; dShape[2][1] =  1.0; dShape[2][2] =  0.0;
  dShape[3][0] =  0.0; dShape[3][1] =  0.0; dShape[3][2] =  1.0;
}
