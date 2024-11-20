// ------------------------------------------------------------
// HB - 05-19-05
// ------------------------------------------------------------
#include <cmath>

//         ls    : array (of length 3*nnodes, see 3 dofs/nodes) defining the position of the
//                 dofs of node i in the element displacement array U: 
//                 U[ls[i       ]]  = Ux of node i
//                 U[ls[i+nnodes]]  = Uy of node i
//                 U[ls[i+2*nnodes]]= Uz of node i
//                 If not provided (i.e. ls = 0), assume following ordering of the given
//                 element displacement array U: [Ux1,Uy1,Uz1,Ux2,Uy2,Uz2, ... ]
void
computeEngStrain3DSolid(double Strain[6], double (*DShape)[3], double* U, int nnodes, int* ls=0)
{
  Strain[0] = Strain[1] = Strain[2] = Strain[3] = Strain[4] = Strain[5] = 0.0;
  int ix, iy, iz;
  for(int i=0; i<nnodes; i++){
    if(ls) { ix = ls[i]; iy = ls[i+nnodes]; iy = ls[i+nnodes]; }
    else { ix = 3*i; iy = ix+1; iz = iy+1; }
    Strain[0] += DShape[i][0]*U[ix]; // e11
    Strain[1] += DShape[i][1]*U[iy]; // e22
    Strain[2] += DShape[i][2]*U[iz]; // e33
                                                                                                                                                                     
    Strain[3] += DShape[i][1]*U[ix] + DShape[i][0]*U[iy]; // g12 = 2.e12
    Strain[4] += DShape[i][2]*U[iy] + DShape[i][1]*U[iz]; // g23 = 2.e23
    Strain[5] += DShape[i][2]*U[ix] + DShape[i][0]*U[iz]; // g13 = 2.e13
  }
}
                                                                                                                                        
void
computeStress3DSolid(double Stress[6], double Strain[6], double C[6][6])
{
  Stress[0] = C[0][0]*Strain[0] + C[0][1]*Strain[1] + C[0][2]*Strain[2] + C[0][3]*Strain[3] + C[0][4]*Strain[4] + C[0][5]*Strain[5];
  Stress[1] = C[1][0]*Strain[0] + C[1][1]*Strain[1] + C[1][2]*Strain[2] + C[1][3]*Strain[3] + C[1][4]*Strain[4] + C[1][5]*Strain[5];
  Stress[2] = C[2][0]*Strain[0] + C[2][1]*Strain[1] + C[2][2]*Strain[2] + C[2][3]*Strain[3] + C[2][4]*Strain[4] + C[2][5]*Strain[5];
  Stress[3] = C[3][0]*Strain[0] + C[3][1]*Strain[1] + C[3][2]*Strain[2] + C[3][3]*Strain[3] + C[3][4]*Strain[4] + C[3][5]*Strain[5];
  Stress[4] = C[4][0]*Strain[0] + C[4][1]*Strain[1] + C[4][2]*Strain[2] + C[4][3]*Strain[3] + C[4][4]*Strain[4] + C[4][5]*Strain[5];
  Stress[5] = C[5][0]*Strain[0] + C[5][1]*Strain[1] + C[5][2]*Strain[2] + C[5][3]*Strain[3] + C[5][4]*Strain[4] + C[5][5]*Strain[5];
}

//         ls    : array (of length 3*nnodes, see 3 dofs/nodes) defining the position of the
//                 dofs of node i in the element displacement array U:
//                 U[ls[i       ]]  = Ux of node i
//                 U[ls[i+nnodes]]  = Uy of node i
//                 U[ls[i+2*nnodes]]= Uz of node i
//                 If not provided (i.e. ls = 0), assume following ordering of the given
//                 element displacement array U: [Ux1,Uy1,Uz1,Ux2,Uy2,Uz2, ...]
void
computeStressAndEngStrain3DSolid(double Stress[6], double Strain[6], double C[6][6], double (*DShape)[3], double* U, int nnodes, int* ls=0)
{
  computeEngStrain3DSolid(Strain, DShape, U, nnodes, ls);
  computeStress3DSolid(Stress, Strain, C);
}

double
computeVonMisesStress(double Stress[6])
{
  double p = (Stress[0]+Stress[1]+Stress[2])/3.;
  double dsxx = Stress[0] - p;
  double dsyy = Stress[1] - p;
  double dszz = Stress[2] - p;
  double dsxy = Stress[3];
  double dsyz = Stress[4];
  double dsxz = Stress[5];

  double j2 = ((dsxx*dsxx)+(dsyy*dsyy)+(dszz*dszz))/2. +
	       (dsxy*dsxy)+(dsyz*dsyz)+(dsxz*dsxz);
	       
  return(sqrt(3.*j2));
}

double
computeVonMisesStrain(double Strain[6])
{
  double ev = (Strain[0]+Strain[1]+Strain[2])/3.;
  double dsxx = Strain[0] - ev;
  double dsyy = Strain[1] - ev;
  double dszz = Strain[2] - ev;
  double dsxy = Strain[3];
  double dsyz = Strain[4];
  double dsxz = Strain[5];

  double j2 = ((dsxx*dsxx)+(dsyy*dsyy)+(dszz*dszz))/2. +
	       (dsxy*dsxy)+(dsyz*dsyz)+(dsxz*dsxz);
	       
  return(sqrt(3.*j2));
}
