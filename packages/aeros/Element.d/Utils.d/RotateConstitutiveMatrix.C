#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif

//#define SYMMETRIZE 
//#define FILTER_SMALLTERMS 
#define LOOP_UNROLL
void
printMat(double *A, int n, int m, char* mssg=0)
{
  if(mssg) fprintf(stderr,"%s\n",mssg);
  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++)
      fprintf(stderr," %e ",A[i*m+j]);
    fprintf(stderr,"\n");
  }
}

//HB: transform/rotate a given 6x6 constitutive operator Cin according to the 
//    3x3 transformation/rotation operator T33:
//       Cout = T'.Cin.T
//    Assume constitutive matric Cin maps
//      [e11, e22, e33, 2.e12, 2.e13, 2.e23] to [s11, s22, s33, s12, s13, s23]
void
rotateConstitutiveMatrix(double *_Cin, double *T33, double Cout[6][6])
{
  double Cin[6][6]; for(int i=0; i<6; ++i) for(int j=i; j<6; ++j) { Cin[i][j] = Cin[j][i] = _Cin[6*i+j]; }
  double l1 = T33[0]; double m1 = T33[1]; double n1 = T33[2];
  double l2 = T33[3]; double m2 = T33[4]; double n2 = T33[5];
  double l3 = T33[6]; double m3 = T33[7]; double n3 = T33[8];

  double TNN[3][3] = {{l1*l1,m1*m1,n1*n1},{l2*l2,m2*m2,n2*n2},{l3*l3,m3*m3,n3*n3}};
  double TNS[3][3] = {{l1*m1,l1*n1,m1*n1},{l2*m2,l2*n2,m2*n2},{l3*m3,l3*n3,m3*n3}};
  double TSN[3][3] = {{l1*l2,m1*m2,n1*n2},{l1*l3,m1*m3,n1*n3},{l2*l3,m2*m3,n2*n3}};
  double TSS[3][3] = {{l1*m2+l2*m1,l1*n2+l2*n1,m1*n2+m2*n1},
                      {l1*m3+l3*m1,l1*n3+l3*n1,m1*n3+m3*n1},
                      {l2*m3+l3*m2,l2*n3+l3*n2,m2*n3+m3*n2}};

  // need to check for the presence or not of the sqrt(2) in the stress/strain vectors
  double sqrt2 = sqrt(2.0);
  for(int i=0; i<3; i++)
#ifdef LOOP_UNROLL
    { TNS[i][0] *= sqrt2; TNS[i][1] *= sqrt2; TNS[i][2] *= sqrt2; 
      TSN[i][0] *= sqrt2; TSN[i][1] *= sqrt2; TSN[i][2] *= sqrt2; }
#else
    for(int j=0; j<3; j++){ TNS[i][j] *= sqrt2;  TSN[i][j] *= sqrt2; }
#endif

  double T66[6][6] = {{TNN[0][0],TNN[0][1],TNN[0][2],TNS[0][0],TNS[0][1],TNS[0][2]},
                      {TNN[1][0],TNN[1][1],TNN[1][2],TNS[1][0],TNS[1][1],TNS[1][2]},  
                      {TNN[2][0],TNN[2][1],TNN[2][2],TNS[2][0],TNS[2][1],TNS[2][2]},  
                      {TSN[0][0],TSN[0][1],TSN[0][2],TSS[0][0],TSS[0][1],TSS[0][2]},
                      {TSN[1][0],TSN[1][1],TSN[1][2],TSS[1][0],TSS[1][1],TSS[1][2]},  
                      {TSN[2][0],TSN[2][1],TSN[2][2],TSS[2][0],TSS[2][1],TSS[2][2]}};

  // Chat = Rhat C Rhat
#ifdef LOOP_UNROLL
  for(int i=0; i<6; i++) { Cin[i][3] *= sqrt2; Cin[i][4] *= sqrt2; Cin[i][5] *= sqrt2; }
  for(int j=0; j<6; j++) { Cin[3][j] *= sqrt2; Cin[4][j] *= sqrt2; Cin[5][j] *= sqrt2; }
#else
  for(int i=0; i<6; i++) for(int j=3; j<6; j++) Cin[i][j] *= sqrt2;
  for(int i=3; i<6; i++) for(int j=0; j<6; j++) Cin[i][j] *= sqrt2;
#endif

  // Compute T'.(C.T) (could call a level 3 BLAS routine for optimization)
  double CC[6][6] = {0.};

  for(int i=0; i<6; i++)
    for(int j=0; j<6; j++)
#ifdef LOOP_UNROLL
      CC[i][j] = Cin[i][0]*T66[0][j] + Cin[i][1]*T66[1][j] + Cin[i][2]*T66[2][j]
               + Cin[i][3]*T66[3][j] + Cin[i][4]*T66[4][j] + Cin[i][5]*T66[5][j];
#else
      for(int k=0; k<6; k++) CC[i][j] += Cin[i][k]*T66[k][j];
#endif

  for(int i=0; i<6; i++)
#ifdef LOOP_UNROLL
    for(int j=i; j<6; j++) // take advantage of symmetry
        Cout[i][j] = Cout[j][i] = T66[0][i]*CC[0][j] + T66[1][i]*CC[1][j] + T66[2][i]*CC[2][j]
                                + T66[3][i]*CC[3][j] + T66[4][i]*CC[4][j] + T66[5][i]*CC[5][j];
#else
    for(int j=0; j<6; j++){
      Cout[i][j] = 0.0;
      for(int k=0; k<6; k++) Cout[i][j] += T66[k][i]*CC[k][j];
    } 
#endif

  // Chat = Rhat C Rhat
  double isqrt2 = 1./sqrt2;
#ifdef LOOP_UNROLL
  for(int i=0; i<6; i++) { Cout[i][3] *= isqrt2; Cout[i][4] *= isqrt2; Cout[i][5] *= isqrt2; }
  for(int j=0; j<6; j++) { Cout[3][j] *= isqrt2; Cout[4][j] *= isqrt2; Cout[5][j] *= isqrt2; }
#else
  for(int i=0; i<6; i++) for(int j=3; j<6; j++) Cout[i][j] *= isqrt2;
  for(int i=3; i<6; i++) for(int j=0; j<6; j++) Cout[i][j] *= isqrt2;
#endif

#ifdef SYMMETRIZE // enforce symmetry
  for(int i=0; i<6; i++)
    for(int j=i+1; j<6; j++){
      Cout[i][j] = 0.5*(Cout[i][j]+Cout[j][i]);
      Cout[j][i] = Cout[i][j];
    }
#endif
#ifdef FILTER_SMALLTERMS
  double maxdiag = fabs(Cout[0][0]);
  for(int i=1; i<6; i++) if(fabs(Cout[i][i])>maxdiag) maxdiag = fabs(Cout[i][i]);
  for(int i=0; i<6; i++)
    for(int j=i+1; j<6; j++)
      if(fabs(Cout[i][j])/maxdiag < 1.E-14) {
         Cout[i][j] = Cout[j][i] = 0.0;
      }
#endif
}

void rotateVector(double *_w, double *cFrame, double *_alpha)
{
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(cFrame);
  Eigen::Map<Eigen::Matrix<double,6,1> > alpha(_alpha);
  Eigen::Map<Eigen::Matrix<double,6,1> > w(_w);
  alpha.head<3>() = T.transpose()*w.head<3>();
  alpha.tail<3>() = T.transpose()*w.tail<3>();
#else
  std::cerr << " *** WARNING: USE_EIGEN3 is not defined in Element.d/Utils.d/RotateConstitutiveMatrix.C\n";
#endif
}

/*
int main()
{
  double E = 1.0;
  double nu = 0.3;
  double lam = nu*E/((1.0+nu)*(1.0-2.0*nu));
  double mu  = E/(2.0*(1.0+nu));
  double FEMCiso[6][6] = {{lam+2*mu,   lam  ,   lam  ,0.0,0.0,0.0},
                          {   lam  ,lam+2*mu,   lam  ,0.0,0.0,0.0},
                          {   lam  ,   lam  ,lam+2*mu,0.0,0.0,0.0},
                          {   0.0  ,   0.0  ,   0.0  , mu,0.0,0.0},
                          {   0.0  ,   0.0  ,   0.0  , 0.0,mu,0.0},
                          {   0.0  ,   0.0  ,   0.0  , 0.0,0.0,mu}};
  
  printMat(&FEMCiso[0][0],6,6,"FEMCiso = ");
  double theta = 1.0;
  double T[3][3] = {{cos(theta),-sin(theta),0.0},
                    {sin(theta), cos(theta),0.0},
                    {    0.0   ,     0.0   ,1.0}};
  printMat(&T[0][0],3,3,"T = ");

  double NewFEMCiso[6][6];
  rotateConstitutiveMatrix(FEMCiso,T,NewFEMCiso);
  printMat(&NewFEMCiso[0][0],6,6,"NewFEMCiso = ");
 
  return(1);
}
*/
