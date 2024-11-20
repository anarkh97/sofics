// ------------------------------------------------------------
// HB - 05-16-05
// ------------------------------------------------------------
// Hexa8 shape fct & derivative. 
// Wrapper for Fortran routine h8shpe.f
// ------------------------------------------------------------
                                                                                                                                       
// Std C/C++ lib
#include <cstdlib>
#include <Utils.d/linkfc.h>

extern "C" {
void _FORTRAN(h8shpe)(double &, double &, double &, double *, double *,
                      double *, double *, double *, double *, double *,
                      double &);
}
                                                                                                                                       
double 
Hexa8ShapeFct(double Shape[8], double DShape[8][3],double m[3], double X[8], double Y[8], double Z[8])
{
  double xi = m[0]; 
  double eta= m[1]; 
  double emu= m[2]; 
  double qx[8], qy[8], qz[8];
  double J;
  
  _FORTRAN(h8shpe)(xi,eta,emu,X,Y,Z,Shape,qx,qy,qz,J);

  DShape[0][0] = qx[0]; DShape[0][1] = qy[0]; DShape[0][2] = qz[0];
  DShape[1][0] = qx[1]; DShape[1][1] = qy[1]; DShape[1][2] = qz[1];
  DShape[2][0] = qx[2]; DShape[2][1] = qy[2]; DShape[2][2] = qz[2];
  DShape[3][0] = qx[3]; DShape[3][1] = qy[3]; DShape[3][2] = qz[3];
  DShape[4][0] = qx[4]; DShape[4][1] = qy[4]; DShape[4][2] = qz[4];
  DShape[5][0] = qx[5]; DShape[5][1] = qy[5]; DShape[5][2] = qz[5];
  DShape[6][0] = qx[6]; DShape[6][1] = qy[6]; DShape[6][2] = qz[6];
  DShape[7][0] = qx[7]; DShape[7][1] = qy[7]; DShape[7][2] = qz[7];
  
  return(J);
}

