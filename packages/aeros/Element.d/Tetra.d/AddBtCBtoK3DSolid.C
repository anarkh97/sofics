// ------------------------------------------------------------
// HB - 04-15-05 
// ------------------------------------------------------------
#include <Math.d/FullSquareMatrix.h>

// HB (04-15-05): compute Bt.C.B & add to K for 3D solid el. (TET...,BRICK...,PENTA...,PYRA...)
//                K <- K + alpha.(Bt.C.B)
//                Assume symmetry of K.
//                Assume constitutive matric C maps 
//                  [e11, e22, e33, 2.e12, 2.e13, 2.e23] to [s11, s22, s33, s12, s13, s23]
//
// Inputs: K     : elementary stiffness matrix
//         DShape: derivative of shape functions with respect to global coordinates
//                 DShape[i][j]: derivative of shape fct i w.t.r direction j  (0=X, 1=Y, 2=Z)
//         C     : 6x6 constitutive matrix (assume C symmetric)
//         alpha : multiplier coefficient (i.e. K <- K + alpha.(Bt.C.B))
//         nnodes: number of nodes (i.e. shape functions) 
//         ls    : array (of length 3*nnodes, see 3 dofs/nodes) defining the arrangement of rows 
//                 and columns of K:
//                 for node i, ls[i]          = row/column position of nodal dof X in K
//                             ls[i+nnodes]   = row/column position of nodal dof Y in K
//                             ls[i+2*nnodes] = row/column position of nodal dof Y in K
void
AddBtCBtoK3DSolid(FullSquareMatrix &K, double (*DShape)[3], double C[6][6], double alpha, int nnodes, int* ls)
{
  for(int j=0; j<nnodes; j++){
    int jx = ls[j];
    int jy = ls[j+nnodes];
    int jz = ls[j+2*nnodes];
    double c1x = (C[0][0]*DShape[j][0] + C[0][3]*DShape[j][1] + C[0][5]*DShape[j][2])*alpha;
    double c1y = (C[0][3]*DShape[j][0] + C[0][1]*DShape[j][1] + C[0][4]*DShape[j][2])*alpha;
    double c1z = (C[0][5]*DShape[j][0] + C[0][4]*DShape[j][1] + C[0][2]*DShape[j][2])*alpha;
    double c2x = (C[1][0]*DShape[j][0] + C[1][3]*DShape[j][1] + C[1][5]*DShape[j][2])*alpha;
    double c2y = (C[1][3]*DShape[j][0] + C[1][1]*DShape[j][1] + C[1][4]*DShape[j][2])*alpha;
    double c2z = (C[1][5]*DShape[j][0] + C[1][4]*DShape[j][1] + C[1][2]*DShape[j][2])*alpha;
    double c3x = (C[2][0]*DShape[j][0] + C[2][3]*DShape[j][1] + C[2][5]*DShape[j][2])*alpha;
    double c3y = (C[2][3]*DShape[j][0] + C[2][1]*DShape[j][1] + C[2][4]*DShape[j][2])*alpha;
    double c3z = (C[2][5]*DShape[j][0] + C[2][4]*DShape[j][1] + C[2][2]*DShape[j][2])*alpha;
    double c4x = (C[3][0]*DShape[j][0] + C[3][3]*DShape[j][1] + C[3][5]*DShape[j][2])*alpha;
    double c4y = (C[3][3]*DShape[j][0] + C[3][1]*DShape[j][1] + C[3][4]*DShape[j][2])*alpha;
    double c4z = (C[3][5]*DShape[j][0] + C[3][4]*DShape[j][1] + C[3][2]*DShape[j][2])*alpha;
    double c5x = (C[4][0]*DShape[j][0] + C[4][3]*DShape[j][1] + C[4][5]*DShape[j][2])*alpha;
    double c5y = (C[4][3]*DShape[j][0] + C[4][1]*DShape[j][1] + C[4][4]*DShape[j][2])*alpha;
    double c5z = (C[4][5]*DShape[j][0] + C[4][4]*DShape[j][1] + C[4][2]*DShape[j][2])*alpha;
    double c6x = (C[5][0]*DShape[j][0] + C[5][3]*DShape[j][1] + C[5][5]*DShape[j][2])*alpha;
    double c6y = (C[5][3]*DShape[j][0] + C[5][1]*DShape[j][1] + C[5][4]*DShape[j][2])*alpha;
    double c6z = (C[5][5]*DShape[j][0] + C[5][4]*DShape[j][1] + C[5][2]*DShape[j][2])*alpha;
    for(int i=j; i<nnodes; i++){
      int ix = ls[i];
      int iy = ls[i+nnodes];
      int iz = ls[i+2*nnodes];		    
      K[ix][jx] += DShape[i][0]*c1x + DShape[i][1]*c4x + DShape[i][2]*c6x;
      K[jx][ix]  = K[ix][jx];
      K[iy][jy] += DShape[i][0]*c4y + DShape[i][1]*c2y + DShape[i][2]*c5y;
      K[jy][iy]  = K[iy][jy];
      K[iz][jz] += DShape[i][0]*c6z + DShape[i][1]*c5z + DShape[i][2]*c3z;
      K[jz][iz]  = K[iz][jz];
      K[ix][jy] += DShape[i][0]*c1y + DShape[i][1]*c4y + DShape[i][2]*c6y;
      K[iy][jx] += DShape[i][0]*c4x + DShape[i][1]*c2x + DShape[i][2]*c5x;
      K[jy][ix]  = K[ix][jy];
      K[jx][iy]  = K[iy][jx];
      K[ix][jz] += DShape[i][0]*c1z + DShape[i][1]*c4z + DShape[i][2]*c6z;
      K[iz][jx] += DShape[i][0]*c6x + DShape[i][1]*c5x + DShape[i][2]*c3x;
      K[jz][ix]  = K[ix][jz];
      K[jx][iz]  = K[iz][jx];
      K[iy][jz] += DShape[i][0]*c4z + DShape[i][1]*c2z + DShape[i][2]*c5z;
      K[iz][jy] += DShape[i][0]*c6y + DShape[i][1]*c5y + DShape[i][2]*c3y;
      K[jz][iy]  = K[iy][jz];
      K[jy][iz]  = K[iz][jy];
    }
  } 
}  


/*for(int j=0; j<nnodes; j++){
  int jx = ls[j];
  int jy = ls[j+nnodes];
  int jz = ls[j+2*nnodes];
  double c1x = (C[0][0]*nGrad[j][0] + C[0][3]*nGrad[j][1] + C[0][5]*nGrad[j][2])*w;
  double c1y = (C[0][3]*nGrad[j][0] + C[0][1]*nGrad[j][1] + C[0][4]*nGrad[j][2])*w;
  double c1z = (C[0][5]*nGrad[j][0] + C[0][4]*nGrad[j][1] + C[0][2]*nGrad[j][2])*w;
  double c2x = (C[1][0]*nGrad[j][0] + C[1][3]*nGrad[j][1] + C[1][5]*nGrad[j][2])*w;
  double c2y = (C[1][3]*nGrad[j][0] + C[1][1]*nGrad[j][1] + C[1][4]*nGrad[j][2])*w;
  double c2z = (C[1][5]*nGrad[j][0] + C[1][4]*nGrad[j][1] + C[1][2]*nGrad[j][2])*w;
  double c3x = (C[2][0]*nGrad[j][0] + C[2][3]*nGrad[j][1] + C[2][5]*nGrad[j][2])*w;
  double c3y = (C[2][3]*nGrad[j][0] + C[2][1]*nGrad[j][1] + C[2][4]*nGrad[j][2])*w;
  double c3z = (C[2][5]*nGrad[j][0] + C[2][4]*nGrad[j][1] + C[2][2]*nGrad[j][2])*w;
  double c4x = (C[3][0]*nGrad[j][0] + C[3][3]*nGrad[j][1] + C[3][5]*nGrad[j][2])*w;
  double c4y = (C[3][3]*nGrad[j][0] + C[3][1]*nGrad[j][1] + C[3][4]*nGrad[j][2])*w;
  double c4z = (C[3][5]*nGrad[j][0] + C[3][4]*nGrad[j][1] + C[3][2]*nGrad[j][2])*w;
  double c5x = (C[4][0]*nGrad[j][0] + C[4][3]*nGrad[j][1] + C[4][5]*nGrad[j][2])*w;
  double c5y = (C[4][3]*nGrad[j][0] + C[4][1]*nGrad[j][1] + C[4][4]*nGrad[j][2])*w;
  double c5z = (C[4][5]*nGrad[j][0] + C[4][4]*nGrad[j][1] + C[4][2]*nGrad[j][2])*w;
  double c6x = (C[5][0]*nGrad[j][0] + C[5][3]*nGrad[j][1] + C[5][5]*nGrad[j][2])*w;
  double c6y = (C[5][3]*nGrad[j][0] + C[5][1]*nGrad[j][1] + C[5][4]*nGrad[j][2])*w;
  double c6z = (C[5][5]*nGrad[j][0] + C[5][4]*nGrad[j][1] + C[5][2]*nGrad[j][2])*w;
  for(int i=j; i<nnodes; i++){
    int ix = ls[i];
    int iy = ls[i+nnodes];
    int iz = ls[i+2*nnodes];		    
    K[ix][jx] = K[ix][jx]+nGrad[i][0]*c1x + nGrad[i][1]*c4x + nGrad[i][2]*c6x;
    K[jx][ix] = K[ix][jx];
    K[iy][jy] = K[iy][jy]+nGrad[i][0]*c4y + nGrad[i][1]*c2y + nGrad[i][2]*c5y;
    K[jy][iy] = K[iy][jy];
    K[iz][jz] = K[iz][jz]+nGrad[i][0]*c6z + nGrad[i][1]*c5z + nGrad[i][2]*c3z;
    K[jz][iz] = K[iz][jz];
    K[ix][jy] = K[ix][jy]+nGrad[i][0]*c1y + nGrad[i][1]*c4y + nGrad[i][2]*c6y;
    K[iy][jx] = K[iy][jx]+nGrad[i][0]*c4x + nGrad[i][1]*c2x + nGrad[i][2]*c5x;
    K[jy][ix] = K[ix][jy];
    K[jx][iy] = K[iy][jx];
    K[ix][jz] = K[ix][jz]+nGrad[i][0]*c1z + nGrad[i][1]*c4z + nGrad[i][2]*c6z;
    K[iz][jx] = K[iz][jx]+nGrad[i][0]*c6x + nGrad[i][1]*c5x + nGrad[i][2]*c3x;
    K[jz][ix] = K[ix][jz];
    K[jx][iz] = K[iz][jx];
    K[iy][jz] = K[iy][jz]+nGrad[i][0]*c4z + nGrad[i][1]*c2z + nGrad[i][2]*c5z;
    K[iz][jy] = K[iz][jy]+nGrad[i][0]*c6y + nGrad[i][1]*c5y + nGrad[i][2]*c3y;
    K[jz][iy] = K[iy][jz];
    K[jy][iz] = K[iz][jy];
  }
} */
