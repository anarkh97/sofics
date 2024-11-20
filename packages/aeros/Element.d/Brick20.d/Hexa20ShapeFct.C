// ------------------------------------------------------------
// HB - 05-19-05
// ------------------------------------------------------------
// Hexa20 shape fct & derivative. 
// Wrapper for Fortran routine h20shpe.f
// ------------------------------------------------------------
                                                                                                                                       
// Std C/C++ lib
#include <cstdlib>
#include <Utils.d/linkfc.h>

extern "C" {
void _FORTRAN(h20shpe)(double &, double &, double &, double *, double *,
                      double *, double *, double *, double *, double *,
                      double &);
}
                                                                                                                                       
double
Hexa20ShapeFct(double Shape[20], double DShape[20][3], double m[3], double X[20], double Y[20], double Z[20])
{
  double xi = m[0]; 
  double eta= m[1]; 
  double emu= m[2]; 
  double qx[20], qy[20], qz[20];
  double J = 0.0;
  
  _FORTRAN(h20shpe)(xi,eta,emu,X,Y,Z,Shape,qx,qy,qz,J);

  for(int i=0; i<20; i+=5) {
    DShape[i+0][0] = qx[i+0]; DShape[i+0][1] = qy[i+0]; DShape[i+0][2] = qz[i+0];
    DShape[i+1][0] = qx[i+1]; DShape[i+1][1] = qy[i+1]; DShape[i+1][2] = qz[i+1];
    DShape[i+2][0] = qx[i+2]; DShape[i+2][1] = qy[i+2]; DShape[i+2][2] = qz[i+2];
    DShape[i+3][0] = qx[i+3]; DShape[i+3][1] = qy[i+3]; DShape[i+3][2] = qz[i+3];
    DShape[i+4][0] = qx[i+4]; DShape[i+4][1] = qy[i+4]; DShape[i+4][2] = qz[i+4];
  }
  return(J);
}

