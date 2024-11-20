// ---------------------------------------------------------------------
// HB - 01-15-06
// ---------------------------------------------------------------------
// Shape fcts & derivatives for 15 nodes wedge element 
// Serendipity finite element basis
// Iso-parametric formulation
#include <cstdio>
#include <iostream>

// HB (01-15-06): for given reference shape fct derivatives dShape & nodes coordinates (X.Y,Z), 
//                compute the jacobian & shape fcts derivative w.r.t to the global coordinate system.
//                assuming an iso-parametric formulation.
double
computePenta15DShapeFct(double dShape[15][3], double X[15], double Y[15], double Z[15], double (*DShape)[3]=0)
{
  double xd1,xd2,xd3;
  double yd1,yd2,yd3;
  double zd1,zd2,zd3;

  xd1 = 0.0 ; xd2 = 0.0 ; xd3 = 0.0 ;
  yd1 = 0.0 ; yd2 = 0.0 ; yd3 = 0.0 ;
  zd1 = 0.0 ; zd2 = 0.0 ; zd3 = 0.0 ;

  for(int i=0; i<15; i+=3){
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

  double a11 = yd2*zd3 - yd3*zd2 ; double a12 = yd3*zd1 - yd1*zd3 ; double a13 = yd1*zd2 - yd2*zd1 ;
  double a21 = xd3*zd2 - xd2*zd3 ; double a22 = xd1*zd3 - zd1*xd3 ; double a23 = xd2*zd1 - xd1*zd2 ;
  double a31 = xd2*yd3 - xd3*yd2 ; double a32 = yd1*xd3 - xd1*yd3 ; double a33 = xd1*yd2 - yd1*xd2 ;

  /* -------> DETERMINANT OF THE JACOBIAN MATRIX <--- */
  double J = xd1*a11 + yd1*a21 + zd1*a31 ;
  if(J == 0.0) {
    fprintf(stderr," *** WARNING: NULL JACOBIAN IN computePenta15DShapeFct.C ROUTINE.\n");
  }
  if(J < 0.0) {
    fprintf(stderr," *** WARNING: NEGATIVE JACOBIAN IN computePenta15DShapeFct.C ROUTINE.\n");
  }
  double cdet = 1.0/J;

  /* -------> DERIVATIVE OF THE SHAPE FCTS IN THE "REAL" ELEMENT <--- */
  if(DShape){
    for(int i=0; i<15; i+=3){
      DShape[i  ][0] = cdet * ( a11*dShape[i  ][0] + a12*dShape[i  ][1] + a13*dShape[i  ][2] );
      DShape[i  ][1] = cdet * ( a21*dShape[i  ][0] + a22*dShape[i  ][1] + a23*dShape[i  ][2] );
      DShape[i  ][2] = cdet * ( a31*dShape[i  ][0] + a32*dShape[i  ][1] + a33*dShape[i  ][2] );
     
      DShape[i+1][0] = cdet * ( a11*dShape[i+1][0] + a12*dShape[i+1][1] + a13*dShape[i+1][2] );
      DShape[i+1][1] = cdet * ( a21*dShape[i+1][0] + a22*dShape[i+1][1] + a23*dShape[i+1][2] );
      DShape[i+1][2] = cdet * ( a31*dShape[i+1][0] + a32*dShape[i+1][1] + a33*dShape[i+1][2] );

      DShape[i+2][0] = cdet * ( a11*dShape[i+2][0] + a12*dShape[i+2][1] + a13*dShape[i+2][2] );
      DShape[i+2][1] = cdet * ( a21*dShape[i+2][0] + a22*dShape[i+2][1] + a23*dShape[i+2][2] );
      DShape[i+2][2] = cdet * ( a31*dShape[i+2][0] + a32*dShape[i+2][1] + a33*dShape[i+2][2] );
    } 
  }
  
  return(J);
}

void Penta15ShapeFct(double Shape[15], double m[3])
{
  double xi = m[0];
  double et = m[1];
  double ze = m[2];
  double a = 1.0-xi-et;

  Shape[0] = -0.5*a*(1.0-ze)*(2.0*xi+2.0*et+ze);
  Shape[1] = 0.5*xi*(1.0-ze)*(2.0*xi-2.0-ze);
  Shape[2] = 0.5*et*(1.0-ze)*(2.0*et-2.0-ze);
  Shape[3] = -0.5*a*(1.0+ze)*(2.0*xi+2.0*et-ze);
  Shape[4] = 0.5*xi*(1.0+ze)*(2.0*xi-2.0+ze);
  Shape[5] = 0.5*et*(1.0+ze)*(2.0*et-2.0+ze);
  Shape[6] = 2.0*xi*a*(1.0-ze);
  Shape[7] = 2.0*xi*et*(1.0-ze);
  Shape[8] = 2.0*et*a*(1.0-ze);
  Shape[9] = 2.0*xi*a*(1.0+ze);
  Shape[10] = 2.0*xi*et*(1.0+ze);
  Shape[11] = 2.0*et*a*(1.0+ze);
  Shape[12] = a*(1.0-ze*ze);
  Shape[13] = xi*(1.0-ze*ze);
  Shape[14] = et*(1.0-ze*ze);

}

void Penta15dShapeFct(double dShape[15][3], double m[3])
{
  double xi = m[0];
  double et = m[1];
  double ze = m[2];
  double a = 1.0-xi-et;
//
// local derivatives of the shape functions: xi-derivative
//
  dShape[0][0] = 0.5*(1.0-ze)*(4.0*xi+4.0*et+ze-2.0);
  dShape[1][0] = 0.5*(1.0-ze)*(4.0*xi-ze-2.0);
  dShape[2][0] = 0.0;
  dShape[3][0] = 0.5*(1.0+ze)*(4.0*xi+4.0*et-ze-2.0);
  dShape[4][0] = 0.5*(1.0+ze)*(4.0*xi+ze-2.0);
  dShape[5][0] = 0.0;
  dShape[6][0] = 2.0*(1.0-ze)*(1.0-2.0*xi-et);
  dShape[7][0] = 2.0*et*(1.0-ze);
  dShape[8][0] = -2.0*et*(1.0-ze);
  dShape[9][0] = 2.0*(1.0+ze)*(1.0-2.0*xi-et);
  dShape[10][0] = 2.0*et*(1.0+ze);
  dShape[11][0] = -2.0*et*(1.0+ze);
  dShape[12][0] = -(1.0-ze*ze);
  dShape[13][0] = (1.0-ze*ze);
  dShape[14][0] = 0.0;
//
// local derivatives of the shape functions: eta-derivative
//
  dShape[0][1] = 0.5*(1.0-ze)*(4.0*xi+4.0*et+ze-2.0);
  dShape[1][1] = 0.0;
  dShape[2][1] = 0.5*(1.0-ze)*(4.0*et-ze-2.0);
  dShape[3][1] = 0.5*(1.0+ze)*(4.0*xi+4.0*et-ze-2.0);
  dShape[4][1] = 0.0;
  dShape[5][1] = 0.5*(1.0+ze)*(4.0*et+ze-2.0);
  dShape[6][1] =-2.0*xi*(1.0-ze);
  dShape[7][1] = 2.0*xi*(1.0-ze);
  dShape[8][1] = 2.0*(1.0-ze)*(1.0-xi-2.0*et);
  dShape[9][1] =-2.0*xi*(1.0+ze);
  dShape[10][1] = 2.0*xi*(1.0+ze);
  dShape[11][1] = 2.0*(1.0+ze)*(1.0-xi-2.0*et);
  dShape[12][1] =-(1.0-ze*ze);
  dShape[13][1] = 0.0;
  dShape[14][1] = (1.0-ze*ze);
//
// local derivatives of the shape functions: zeta-derivative
//
  dShape[0][2] = a*(xi+et+ze-0.5);
  dShape[1][2] = xi*(-xi+ze+0.5);
  dShape[2][2] = et*(-et+ze+0.5);
  dShape[3][2] = a*(-xi-et+ze+0.5);
  dShape[4][2] = xi*(xi+ze-0.5);
  dShape[5][2] = et*(et+ze-0.5);
  dShape[6][2] = -2.0*xi*a;
  dShape[7][2] = -2.0*xi*et;
  dShape[8][2] = -2.0*et*a;
  dShape[9][2] = 2.0*xi*a;
  dShape[10][2] = 2.0*xi*et;
  dShape[11][2] = 2.0*et*a;
  dShape[12][2] = -2.0*a*ze;
  dShape[13][2] = -2.0*xi*ze;
  dShape[14][2] = -2.0*et*ze;
}

// HB (01-15-06): shape fct & derivatives for the 15 nodes wedge element (Iso-parametric formulation)
void 
Penta15ShapeFct(double Shape[15], double dShape[15][3], double m[3])
{
  Penta15ShapeFct(Shape,m);
  Penta15dShapeFct(dShape,m);
}

// HB (01-15-06): shape fct & derivatives for the 15 nodes wedge element(Iso-parametric formulation)
double
Penta15ShapeFct(double Shape[15], double DShape[15][3], double m[3], double X[15], double Y[15], double Z[15])
{
  double dShape[15][3];
  Penta15ShapeFct(Shape, dShape, m);

  return(computePenta15DShapeFct(dShape, X, Y, Z, DShape));
}

