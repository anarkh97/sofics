/*
  FEEL - Finite Elements Educational Library
  Antonini Puppin-Macedo
  File: q4shape.C
  Purpose: Shape functions routines for the 4 node quad
  Version: 1.0 - last update Dec. 11th 1996
*/

#include<cmath>

void Q4shape2
(double jdet[4],double N[4][4],double dNX[4][4],double dNY[4][4],int nint,int element2, double *xx2, double *yy2, int **connection2)

/**************************  void Q4SHAPE  *********************************/
/*
#   This routine calculates the shape functions and its derivatives
#   for the four node bilinear quadrilateral element Q4.
#
#   N[i][j]    --> shape function
#   dNX[i][j]  --> global N,x shape function derivative dN/dx
#   dNY[i][j]  --> global N,y shape function derivative dN/dy
#
#   eta,xi     --> local element coordinates ("xi" and "eta")
#
#   dNx[i][j]  --> local "xi" derivative shape function
#   dNe[i][j]  --> local "eta" derivative of shape function
#   dxx[j]     --> local x,xi derivative: dx/dxi
#   dyx[j]     --> local y,xi derivative: dy/dxi
#   dxe[j]     --> local x,eta derivative: dx/deta
#   dye[j]     --> local y,eta derivative: dy/deta
#   i          --> local node number
#   j          --> integration point number
#   nint       --> number of integration points, equal to 1 or 4
#   jdet[j]    --> jacobian matrix determinant
*/
{
  int          i, j;
  double        aux1;
  double        x, y; // node coordinates

  // Memory allocation
  double xi[4];
  double eta[4];
  double xi_int[4];
  double eta_int[4];
  double dNx[4][4];
  double dNe[4][4];
  double dxx[4];
  double dyx[4];
  double dxe[4];
  double dye[4];

  for (i=0;i<4;i++) {
    xi[i] = 0.0; eta[i] = 0.0;
    for (j=0;j<nint;j++) {
      xi_int[j] = 0.0; eta_int[j] = 0.0;
      dxx[j] = 0.0; dyx[j] = 0.0; dxe[j] = 0.0; dye[j] = 0.0;
      dNx[i][j] = 0.0; dNe[i][j] = 0.0;
    }
  }

  // Local "xi" coordinates of the serendipity 4 node bilinear element
  xi[0]=-1.0; xi[1]=1.0; xi[2]=1.0; xi[3]=-1.0;
  // Local "eta" coordinates of the serendipity 4 node bilinear element
  eta[0]=-1.0; eta[1]=-1.0; eta[2]=1.0; eta[3]=1.0;

  // The following coordinates for 1 point integration
  // Local "xi" coordinates of integration point
  xi_int[0]=0.0;
  // Local "eta" coordinates of integration point
  eta_int[0]=0.0;

  // The following coordinates for 4 points integration
  if (nint == 4)
    {
      aux1=1.0/sqrt(3.0);
      // Local "xi" coordinates of integration points
      xi_int[0]=-aux1; xi_int[1]=aux1; xi_int[2]=aux1; xi_int[3]=-aux1;
      // Local "eta" coordinates of integration points
      eta_int[0]=-aux1; eta_int[1]=-aux1; eta_int[2]=aux1; eta_int[3]=aux1;
    }

  /*
    # This loop calculates the shape functions N[i][j],
    # the derivatives of the shape functions to xi  dNx[i][j] and
    # the derivatives of the shape functions to eta dNe[i][j].
    # The index i indicates the shape function (and local node) number.
    # The index j indicates the integration point number.
  */
  for (i=0;i<4;i++)
    for (j=0;j<nint;j++)
      {
        N[i][j] = (0.25)*(1.0+xi[i]*xi_int[j])*(1.0+eta[i]*eta_int[j]);
        dNx[i][j] = (0.25)*(xi[i]+xi[i]*eta[i]*eta_int[j]);
        dNe[i][j] = (0.25)*(eta[i]+xi[i]*eta[i]*xi_int[j]);
      }


  /*
    # This loop calculates the derivatives: dx/dxi, dx/deta, dy/dxi, dy/deta,
    # (or x,xi  x,eta  y,xi  y,eta) and the jacobian determinant jdet.
    # The index i indicates the shape function number.
    # The index j indicates the integration point number.
  */

  for (j=0;j<nint;j++)
    {
      for (i=0;i<4;i++)
        {
          x = xx2[connection2[element2][i]-1];
          y = yy2[connection2[element2][i]-1];
          dxx[j] = dxx[j] + dNx[i][j]*x;
          dxe[j] = dxe[j] + dNe[i][j]*x;
          dyx[j] = dyx[j] + dNx[i][j]*y;
          dye[j] = dye[j] + dNe[i][j]*y;
        }
      jdet[j] = dxx[j]*dye[j] - dxe[j]*dyx[j];
    }


  /*
    # This loop calculates the derivatives of the shape functions
    # to the x and y for the integration points.
    # The index i indicates the shape function number.
    # The index j indicates the integration point number.
  */

  for(i=0;i<4;i++)
    for(j=0;j<nint;j++)
      {
        dNX[i][j] =  (dNx[i][j]*dye[j] - dNe[i][j]*dyx[j])/jdet[j];
        dNY[i][j] = -(dNx[i][j]*dxe[j] - dNe[i][j]*dxx[j])/jdet[j];
      }

} // void Q4shape
