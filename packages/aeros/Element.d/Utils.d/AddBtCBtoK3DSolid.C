// ------------------------------------------------------------
// HB - 04-15-05 
// ------------------------------------------------------------
#include <Math.d/FullSquareMatrix.h>

// HB (04-15-05): compute Bt.C.B & add to K for 3D solid el. (TET...,BRICK...,PENTA...,PYRA...)
//                K <- K + alpha.(Bt.C.B)
//                Assume symmetry of K.
//                Assume 6x6 constitutive matrix C maps 
//                  [e11, e22, e33, 2.e12, 2.e13, 2.e23] to [s11, s22, s33, s12, s13, s23]
//                B is the operator that maps the dofs U to [e11, e22, e33, 2.e12, 2.e13, 2.e23] 
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
//                             ls[i+2*nnodes] = row/column position of nodal dof Z in K
//
// Outputs: K    : elementary stiffness matrix with added contribution alpha.(Bt.C.B)
//
void
addBtCBtoK3DSolid(FullSquareMatrix &K, double (*DShape)[3], double C[6][6], double alpha, int nnodes, int* ls)
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

// HB (05-05-05): compute Nt.D.N & add to M for 3D solid el. (TET...,BRICK...,PENTA...,PYRA...)
//                M <- M + alpha.(Nt.D.N)
//                Assume symmetry of M.
//                Assume D is a 3x3 SPD operator. Usualy, D = rho.I 
//                N are the shape functions
//
// Inputs: M     : elementary mass matrix
//         Shape : shape functions (Shape[i] = shape fct on node i) 
//         alpha : multiplier coefficient (i.e. M <- M + alpha.(Nt.D.N))
//         nnodes: number of nodes (i.e. shape functions) 
//         ls    : array (of length 3*nnodes, see 3 dofs/nodes) defining the arrangement of rows 
//                 and columns of M:
//                 for node i, ls[i]          = row/column position of nodal dof X in M
//                             ls[i+nnodes]   = row/column position of nodal dof Y in M
//                             ls[i+2*nnodes] = row/column position of nodal dof Z in M
//         D     : 3x3 "mass constitutive operator" (assume D SPD, stored row-wise)
//                 If NOT provided assume D = Identity
//
// Outputs:M     : elementary mass matrix with added contribution alpha.(Nt.D.N)
//
void
addNtDNtoM3DSolid(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls, double (*D)[3] = 0)
{
  if(D) {
    for(int j=0; j<nnodes; j++){
      int jx = ls[j];
      int jy = ls[j+nnodes];
      int jz = ls[j+2*nnodes];
      for(int i=j; i<nnodes; i++){
        int ix = ls[i];
        int iy = ls[i+nnodes];
        int iz = ls[i+2*nnodes];		    
        M[ix][jx] += alpha*(Shape[i]*D[0][0]*Shape[j]);
        M[jx][ix]  = M[ix][jx];

        M[iy][jy] += alpha*(Shape[i]*D[1][1]*Shape[j]);
        M[jy][iy]  = M[iy][jy];
        
        M[iz][jz] += alpha*(Shape[i]*D[2][2]*Shape[j]);
        M[jz][iz]  = M[ix][jx];
        
        M[ix][jy] += alpha*(Shape[i]*D[0][1]*Shape[j]);
        M[jy][ix]  = M[ix][jy];

        M[ix][jz] += alpha*(Shape[i]*D[0][2]*Shape[j]);
        M[jz][ix]  = M[ix][jz];
        
        M[iy][jz] += alpha*(Shape[i]*D[1][2]*Shape[j]);
        M[jz][iy]  = M[iy][jz];
      }
    } 
  } else {
    for(int j=0; j<nnodes; j++){
      int jx = ls[j];
      int jy = ls[j+nnodes];
      int jz = ls[j+2*nnodes];
      for(int i=j; i<nnodes; i++){
        int ix = ls[i];
        int iy = ls[i+nnodes];
        int iz = ls[i+2*nnodes];		    
        M[ix][jx] += alpha*Shape[i]*Shape[j];
        M[jx][ix]  = M[ix][jx];

        M[iy][jy]  = M[ix][jx];
        M[jy][iy]  = M[iy][jy];
        
        M[iz][jz]  = M[ix][jx];
        M[jz][iz]  = M[ix][jx];
      }
    }   
  }
}  

/*
void
addBtSigmaToFint3DSolid(Vector& F, double Sigma[6], double (*DShape)[3], double alpha, int nnodes, int* ls)
{
   for(int i=0; i<nnodes; i++){
     int ix = ls[i];
     int iy = ls[i+nnodes];
     int iz = ls[i+2*nnodes];
     F[ix] += alpha*(DShape[i][0]*Sigma[0] + DShape[i][1]*Sigma[3] + DShape[i][2]*Sigma[5]);
     F[iy] += alpha*(DShape[i][0]*Sigma[3] + DShape[i][1]*Sigma[1] + DShape[i][2]*Sigma[4]);
     F[iz] += alpha*(DShape[i][0]*Sigma[5] + DShape[i][1]*Sigma[4] + DShape[i][2]*Sigma[2]);
  }
}
*/

// HB (05-05-05): compute Bt.B & add to K for 3D Helmhotz el. (TET...,BRICK...,PENTA...,PYRA...)
//                K <- K + alpha.(Bt.B)
//                Assume symmetry of K.
//                B is the operator that maps the dofs p to [dp/dX dp/dY dp/dZ] 
//
// Inputs: K     : elementary stiffness matrix
//         DShape: derivative of shape functions with respect to global coordinates
//                 DShape[i][j]: derivative of shape fct i w.t.r direction j  (0=X, 1=Y, 2=Z)
//         alpha : multiplier coefficient (i.e. K <- K + alpha.(Bt.B))
//         nnodes: number of nodes (i.e. shape functions) 
//         ls    : array (of length nnodes, see 1 dofs/nodes) defining the arrangement of rows 
//                 and columns of K:
//                 for node i, ls[i] = row/column position of nodal dof p in K
//
// Outputs: K    : elementary stiffness matrix with added contribution alpha.(Bt.B)
//
void
addBtBtoK3DHelm(FullSquareMatrix &K, double (*DShape)[3], double alpha, int nnodes, int* ls)
{
  for(int j=0; j<nnodes; j++){
    int jx = ls[j];
    for(int i=j; i<nnodes; i++){
      int ix = ls[i];
      K[ix][jx] += alpha*(DShape[i][0]*DShape[j][0] + DShape[i][1]*DShape[j][1] + DShape[i][2]*DShape[j][2]);
      K[jx][ix]  = K[ix][jx];
    }
  }
}

// HB (05-05-05): compute Nt.N & add to M for 3D Helmhotz el. (TET...,BRICK...,PENTA...,PYRA...)
//                M <- M + alpha.(Nt.N)
//                Assume symmetry of M.
//                N are the shape fucntion p = sum[Ni.pi]
//
// Inputs: M     : elementary mass matrix
//         Shape : shape functions (Shape[i] = shape fct on node i) 
//         alpha : multiplier coefficient (i.e. M <- M + alpha.(Nt.N))
//         nnodes: number of nodes (i.e. shape functions) 
//         ls    : array (of length nnodes, see 1 dofs/nodes) defining the arrangement of rows 
//                 and columns of M:
//                 for node i, ls[i] = row/column position of nodal dof p in M
//
// Outputs: M    : elementary mass matrix with added contribution alpha.(Nt.N)
//
void
addNtNtoM3DHelm(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls)
{
  for(int j=0; j<nnodes; j++){
    int jx = ls[j];
    for(int i=j; i<nnodes; i++){
      int ix = ls[i];
      M[ix][jx] += alpha*Shape[i]*Shape[j];
      M[jx][ix]  = M[ix][jx];
    }
  }
}


void
addDBtCDBtodKdx3DSolid(FullSquareMatrix *&dKdx, double (*DShape)[3], double DDShape[4][3][12], double C[6][6], double alpha, int nnodes, int* ls)
{

  for(int k=0; k<12; ++k) {
    for(int j=0; j<nnodes; j++) {
      int jx = ls[j];
      int jy = ls[j+nnodes];
      int jz = ls[j+2*nnodes];
      double c1x = (C[0][0]*DDShape[j][0][k] + C[0][3]*DDShape[j][1][k] + C[0][5]*DDShape[j][2][k])*alpha;
      double c1y = (C[0][3]*DDShape[j][0][k] + C[0][1]*DDShape[j][1][k] + C[0][4]*DDShape[j][2][k])*alpha;
      double c1z = (C[0][5]*DDShape[j][0][k] + C[0][4]*DDShape[j][1][k] + C[0][2]*DDShape[j][2][k])*alpha;
      double c2x = (C[1][0]*DDShape[j][0][k] + C[1][3]*DDShape[j][1][k] + C[1][5]*DDShape[j][2][k])*alpha;
      double c2y = (C[1][3]*DDShape[j][0][k] + C[1][1]*DDShape[j][1][k] + C[1][4]*DDShape[j][2][k])*alpha;
      double c2z = (C[1][5]*DDShape[j][0][k] + C[1][4]*DDShape[j][1][k] + C[1][2]*DDShape[j][2][k])*alpha;
      double c3x = (C[2][0]*DDShape[j][0][k] + C[2][3]*DDShape[j][1][k] + C[2][5]*DDShape[j][2][k])*alpha;
      double c3y = (C[2][3]*DDShape[j][0][k] + C[2][1]*DDShape[j][1][k] + C[2][4]*DDShape[j][2][k])*alpha;
      double c3z = (C[2][5]*DDShape[j][0][k] + C[2][4]*DDShape[j][1][k] + C[2][2]*DDShape[j][2][k])*alpha;
      double c4x = (C[3][0]*DDShape[j][0][k] + C[3][3]*DDShape[j][1][k] + C[3][5]*DDShape[j][2][k])*alpha;
      double c4y = (C[3][3]*DDShape[j][0][k] + C[3][1]*DDShape[j][1][k] + C[3][4]*DDShape[j][2][k])*alpha;
      double c4z = (C[3][5]*DDShape[j][0][k] + C[3][4]*DDShape[j][1][k] + C[3][2]*DDShape[j][2][k])*alpha;
      double c5x = (C[4][0]*DDShape[j][0][k] + C[4][3]*DDShape[j][1][k] + C[4][5]*DDShape[j][2][k])*alpha;
      double c5y = (C[4][3]*DDShape[j][0][k] + C[4][1]*DDShape[j][1][k] + C[4][4]*DDShape[j][2][k])*alpha;
      double c5z = (C[4][5]*DDShape[j][0][k] + C[4][4]*DDShape[j][1][k] + C[4][2]*DDShape[j][2][k])*alpha;
      double c6x = (C[5][0]*DDShape[j][0][k] + C[5][3]*DDShape[j][1][k] + C[5][5]*DDShape[j][2][k])*alpha;
      double c6y = (C[5][3]*DDShape[j][0][k] + C[5][1]*DDShape[j][1][k] + C[5][4]*DDShape[j][2][k])*alpha;
      double c6z = (C[5][5]*DDShape[j][0][k] + C[5][4]*DDShape[j][1][k] + C[5][2]*DDShape[j][2][k])*alpha;
      for(int i=j; i<nnodes; i++) {
        int ix = ls[i];
        int iy = ls[i+nnodes];
        int iz = ls[i+2*nnodes];		    
        dKdx[k][ix][jx] += DShape[i][0]*c1x + DShape[i][1]*c4x + DShape[i][2]*c6x;
        dKdx[k][iy][jy] += DShape[i][0]*c4y + DShape[i][1]*c2y + DShape[i][2]*c5y;
        dKdx[k][iz][jz] += DShape[i][0]*c6z + DShape[i][1]*c5z + DShape[i][2]*c3z;
        dKdx[k][ix][jy] += DShape[i][0]*c1y + DShape[i][1]*c4y + DShape[i][2]*c6y;
        dKdx[k][iy][jx] += DShape[i][0]*c4x + DShape[i][1]*c2x + DShape[i][2]*c5x;
        dKdx[k][ix][jz] += DShape[i][0]*c1z + DShape[i][1]*c4z + DShape[i][2]*c6z;
        dKdx[k][iz][jx] += DShape[i][0]*c6x + DShape[i][1]*c5x + DShape[i][2]*c3x;
        dKdx[k][iy][jz] += DShape[i][0]*c4z + DShape[i][1]*c2z + DShape[i][2]*c5z;
        dKdx[k][iz][jy] += DShape[i][0]*c6y + DShape[i][1]*c5y + DShape[i][2]*c3y;
      }
    }
  }
 
  for(int j=0; j<nnodes; j++) {
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
    for(int i=j; i<nnodes; i++) {
      int ix = ls[i];
      int iy = ls[i+nnodes];
      int iz = ls[i+2*nnodes];		    
      for(int k=0; k<12; ++k) {
        dKdx[k][ix][jx] += DDShape[i][0][k]*c1x + DDShape[i][1][k]*c4x + DDShape[i][2][k]*c6x;
        dKdx[k][jx][ix]  = dKdx[k][ix][jx];
        dKdx[k][iy][jy] += DDShape[i][0][k]*c4y + DDShape[i][1][k]*c2y + DDShape[i][2][k]*c5y;
        dKdx[k][jy][iy]  = dKdx[k][iy][jy];
        dKdx[k][iz][jz] += DDShape[i][0][k]*c6z + DDShape[i][1][k]*c5z + DDShape[i][2][k]*c3z;
        dKdx[k][jz][iz]  = dKdx[k][iz][jz];
        dKdx[k][ix][jy] += DDShape[i][0][k]*c1y + DDShape[i][1][k]*c4y + DDShape[i][2][k]*c6y;
        dKdx[k][iy][jx] += DDShape[i][0][k]*c4x + DDShape[i][1][k]*c2x + DDShape[i][2][k]*c5x;
        dKdx[k][jy][ix]  = dKdx[k][ix][jy];
        dKdx[k][jx][iy]  = dKdx[k][iy][jx];
        dKdx[k][ix][jz] += DDShape[i][0][k]*c1z + DDShape[i][1][k]*c4z + DDShape[i][2][k]*c6z;
        dKdx[k][iz][jx] += DDShape[i][0][k]*c6x + DDShape[i][1][k]*c5x + DDShape[i][2][k]*c3x;
        dKdx[k][jz][ix]  = dKdx[k][ix][jz];
        dKdx[k][jx][iz]  = dKdx[k][iz][jx];
        dKdx[k][iy][jz] += DDShape[i][0][k]*c4z + DDShape[i][1][k]*c2z + DDShape[i][2][k]*c5z;
        dKdx[k][iz][jy] += DDShape[i][0][k]*c6y + DDShape[i][1][k]*c5y + DDShape[i][2][k]*c3y;
        dKdx[k][jz][iy]  = dKdx[k][iy][jz];
        dKdx[k][jy][iz]  = dKdx[k][iz][jy];
      }
    }
  } 
}  
