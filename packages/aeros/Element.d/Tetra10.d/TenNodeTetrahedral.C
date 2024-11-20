#include <cstdio>
#include <iostream>
#include <cmath>

#include <Element.d/Tetra10.d/TenNodeTetrahedral.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Corotational.d/Tet10Corotator.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/NLTetrahedral.h>
#include <Element.d/Utils.d/SolidElemUtils.h>
#include <Corotational.d/MatNLCorotator.h>
#include <Driver.d/PolygonSet.h>
#include <Element.d/Helm.d/ARubberF.h>

#define CHECK_JACOBIAN // force check nullity & constant sign of jacobian over el.

extern "C" {
void _FORTRAN(sands25)(const int&, double*, double*, double*,
                       double&, double&, double*, double*,
                       double*, const int&, const int&, const int&,
                       const int&, const int&, const int&);

void _FORTRAN(brkcmt)(double&, double&, double*);
}

double Tetra10ShapeFct(double Shape[10], double DShape[10][3], double m[3], double X[10], double Y[10], double Z[10]);
double computeTet10DShapeFct(double dShape[10][3], double X[10], double Y[10], double Z[10], double (*DShape)[3] = 0);

extern bool useFull;

double weight3d5[15] = { 1.975308731198311E-02, 1.198951396316977E-02,
                         1.198951396316977E-02, 1.198951396316977E-02,
                         1.198951396316977E-02, 1.151136787104540E-02,
                         1.151136787104540E-02, 1.151136787104540E-02,
                         1.151136787104540E-02, 8.818342350423336E-03,
                         8.818342350423336E-03, 8.818342350423336E-03,
                         8.818342350423336E-03, 8.818342350423336E-03,
                         8.818342350423336E-03 };

double gauss3d5[15][3] = { {0.250000000000000000000, 0.250000000000000000000, 0.250000000000000000000},
                           {0.724086765841830901630, 0.091971078052723032789, 0.091971078052723032789},
                           {0.091971078052723032789, 0.724086765841830901630, 0.091971078052723032789},
                           {0.091971078052723032789, 0.091971078052723032789, 0.724086765841830901630},
                           {0.091971078052723032789, 0.091971078052723032789, 0.091971078052723032789},
                           {0.040619116511110274837, 0.319793627829629908390, 0.319793627829629908390},
                           {0.319793627829629908390, 0.040619116511110274837, 0.319793627829629908390},
                           {0.319793627829629908390, 0.319793627829629908390, 0.040619116511110274837},
                           {0.319793627829629908390, 0.319793627829629908390, 0.319793627829629908390},
                           {0.443649167310370844260, 0.443649167310370844260, 0.056350832689629155741},
                           {0.443649167310370844260, 0.056350832689629155741, 0.443649167310370844260},
                           {0.443649167310370844260, 0.056350832689629155741, 0.056350832689629155741},
                           {0.056350832689629155741, 0.443649167310370844260, 0.443649167310370844260},
                           {0.056350832689629155741, 0.443649167310370844260, 0.056350832689629155741},
                           {0.056350832689629155741, 0.056350832689629155741, 0.443649167310370844260} };

TenNodeTetrahedral::TenNodeTetrahedral(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
  nn[2] = nodenums[2];
  nn[3] = nodenums[3];
  nn[4] = nodenums[4];
  nn[5] = nodenums[5];
  nn[6] = nodenums[6];
  nn[7] = nodenums[7];
  nn[8] = nodenums[8];
  nn[9] = nodenums[9];

  cFrame = 0;
  cCoefs = 0;
  mat = 0;
}

TenNodeTetrahedral::~TenNodeTetrahedral()
{
  if(cCoefs && mat) delete mat;
}

Element *
TenNodeTetrahedral::clone()
{
  return new TenNodeTetrahedral(*this);
}

void
TenNodeTetrahedral::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
  nn[5] = table[nn[5]];
  nn[6] = table[nn[6]];
  nn[7] = table[nn[7]];
  nn[8] = table[nn[8]];
  nn[9] = table[nn[9]];
}

void
TenNodeTetrahedral::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
  nn[5] = table[nn[5]];
  nn[6] = table[nn[6]];
  nn[7] = table[nn[7]];
  nn[8] = table[nn[8]];
  nn[9] = table[nn[9]];
}

void
TenNodeTetrahedral::getVonMises(Vector& stress, Vector& weight, CoordSet &cs,
                                Vector& elDisp, int strInd, int surface, double *ndTemps,
                                double ylayer, double zlayer, int avgnum)
{
  if(cCoefs) {
    getVonMisesAniso(stress, weight, cs,
                     elDisp, strInd, surface, ndTemps,
                     ylayer, zlayer, avgnum);
    return;
  }

  weight = 1.0;

  double x[10], y[10], z[10];
  cs.getCoordinates(nn, numNodes(), x, y, z);

  int vmflg = 0, strainFlg = 0;
  // Flags sands25 to calculate Von Mises stress
  if(strInd == 6) vmflg = 1;
  // Flags sands25 to calculate Von Mises strain
  if(strInd == 13) strainFlg = 1;

  const int maxgus = 10;
  const int maxstr = 7;
  const int elm    = 1;
  const int msize  = 1;
  const int outerr = 6;

  double elStress[maxgus][maxstr], elStrain[maxgus][maxstr];

  _FORTRAN(sands25)(elm, x, y, z, prop->E, prop->nu, elDisp.data(),
                    (double*)elStress, (double*)elStrain,
                    maxgus, maxstr, msize, outerr, vmflg, strainFlg);

  if(strInd < 7) {
    double thermalStress[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    if(strInd == 0 || strInd == 1 || strInd == 2) {
      double &Tref  = prop->Ta;
      double &alpha = prop->W;
      for(int i=0; i<10; ++i)
        thermalStress[i] = prop->E*alpha*(ndTemps[i] - Tref)/(1.0 - 2.0*prop->nu);
    }
    for(int i=0; i<10; ++i)
      stress[i] = elStress[i][strInd] - thermalStress[i];
  }
  else if(strInd < 14) {
    for(int i=0; i<10; ++i)
      stress[i] = elStrain[i][strInd-7];
  }
  else {
    for(int i=0; i<10; ++i)
      stress[i] = 0;
  }
}

void
TenNodeTetrahedral::getAllStress(FullM& stress, Vector& weight, CoordSet &cs,
                                 Vector& elDisp, int strInd, int surface, double *ndTemps)
{
  if(cCoefs) {
    getAllStressAniso(stress, weight, cs,
                      elDisp, strInd, surface, ndTemps);
    return;
  }

  weight = 1.0;

  double x[10], y[10], z[10];
  cs.getCoordinates(nn, numNodes(), x, y, z);

  int vmflg, strainFlg;
  vmflg = 0;
  strainFlg = 0;

  const int maxgus = 10;
  const int maxstr = 7;
  const int elm    = 1;
  const int msize  = 1;
  const int outerr = 6;

  double elStress[maxgus][maxstr], elStrain[maxgus][maxstr];

  _FORTRAN(sands25)(elm, x, y, z, prop->E, prop->nu, elDisp.data(),
                    (double*)elStress, (double*)elStrain,
                    maxgus, maxstr, msize, outerr, vmflg, strainFlg);

  // Store all Stress or all Strain as defined by strInd
  int i,j;
  if(strInd == 0) {
    double thermalStress[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    if(ndTemps) {
      double &Tref  = prop->Ta;
      double &alpha = prop->W;      
      for(i=0; i<10; ++i)
        thermalStress[i] = prop->E*alpha*(ndTemps[i] - Tref)/(1.0 - 2.0*prop->nu);
    }
    for (i=0; i<10; ++i) {
      for (j=0; j<3; ++j) {
        stress[i][j] = elStress[i][j] - thermalStress[i];
      }
      for (j=3; j<6; ++j) {
        stress[i][j] = elStress[i][j];
      }
    }
  }
  else {
    for (i=0; i<10; ++i) {
      for (j=0; j<6; ++j) {
        stress[i][j] = elStrain[i][j];
      }
    }
  }

  // Get Element Principals for each node without averaging
  double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double pvec[3] = {0.0,0.0,0.0};

  for (i=0; i<10; ++i) {
    for (j=0; j<6; ++j) {
      svec[j] = stress[i][j];
    }
    // Convert Engineering to Tensor Strains
    if(strInd != 0) {
      svec[3] /= 2;
      svec[4] /= 2;
      svec[5] /= 2;
    }
    pstress(svec,pvec);
    for (j=0; j<3; ++j) {
      stress[i][j+6] = pvec[j];
    }
  }
}

double
TenNodeTetrahedral::getMass(const CoordSet& cs) const
{
  double x[10], y[10], z[10];
  cs.getCoordinates(nn, numNodes(), x, y, z);

  // integration: loop over Gauss pts
  // reuse the 15 pts integration rule (order 5) -> reuse arrays vp1 & dp
  const int numgauss = 15;
  extern double dp[15][10][3]; // contains the values of the Tet10 shape fct at the 15 integration pts
  double dOmega; // det of jacobian
  double volume = 0.0;
  for(int i = 0; i < numgauss; i++) {
    dOmega = computeTet10DShapeFct(dp[i], x, y, z);
    volume += fabs(dOmega)*weight3d5[i];
  }

  return volume*prop->rho;
}

void
TenNodeTetrahedral::getGravityForce(CoordSet& cs, double *gravityAcceleration,
                                    Vector& gravityForce, int gravflg, GeomState *geomState)
{
  int nnodes = 10;

  // Lumped
  if (gravflg != 2) {

    double totmas = getMass(cs);

    double m[30*30];
    FullSquareMatrix result = massMatrix(cs, m, 1);
    std::vector<double> factors;
    lumpMatrix(result, factors);

    // divvy up the total body force using same ratio as the corresponding diagonal of the lumped mass matrix to the total mass
    for(int i = 0; i < nnodes; ++i) {
      gravityForce[3*i+0] = totmas*gravityAcceleration[0]*(3.0*factors[3*i+0]);
      gravityForce[3*i+1] = totmas*gravityAcceleration[1]*(3.0*factors[3*i+1]);
      gravityForce[3*i+2] = totmas*gravityAcceleration[2]*(3.0*factors[3*i+2]);
    }
  }
  // Consistent
  else {

    double x[10], y[10], z[10];
    cs.getCoordinates(nn, numNodes(), x, y, z);

    double lforce[10];
    for(int i = 0; i < nnodes; ++i) lforce[i] = 0.0;

    // integration: loop over Gauss pts
    // reuse the 15 pts integration rule (order 5) -> reuse arrays vp1 & dp
    const int numgauss = 15;
    extern double dp[15][10][3]; // arrays vp1 & dp contain the values of the Tet10 shape fct  & their 
    extern double vp1[15][10];   // derivatives w.r.t reference coordinate system at the 15 integration pts
    double dOmega; // det of jacobian

    for(int i = 0; i < numgauss; i++) {
      dOmega = computeTet10DShapeFct(dp[i], x, y, z);
      double w = fabs(dOmega)*weight3d5[i]*prop->rho;

      for(int n = 0; n < nnodes; ++n)
        lforce[n] += w*vp1[i][n];
    }

    for(int i = 0; i < nnodes; ++i) {
      gravityForce[3*i+0] = lforce[i]*gravityAcceleration[0];
      gravityForce[3*i+1] = lforce[i]*gravityAcceleration[1];
      gravityForce[3*i+2] = lforce[i]*gravityAcceleration[2];
    }
  }
}

void
TenNodeTetrahedral::getThermalForce(CoordSet &cs, Vector &ndTemps,
                                    Vector &elementThermalForce, int glflag,
                                    GeomState *geomState)
{
  // ASSUME CONSTANT THERMAL EXPANSION COEFF. & REFERENCE TEMPERATURE OVER THE ELEMENT 
  const int nnodes = 10;
  const int ndofs = 30;

  // initialize nodal thermal forces
  for(int i=0; i<ndofs; i++) elementThermalForce[i] = 0.0;

  // for nonlinear analyses, the thermal load for this element is now computed in getStiffAndForce
  if(geomState) return;

  double X[10], Y[10], Z[10];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  // get material props & constitutive matrix
  double &Tref = prop->Ta;
  double alpha[6];
  double C[6][6];
  if(cCoefs) { // anisotropic material
    // transform local constitutive matrix to global frame
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
    // transform local coefficients of thermal expansion to global frame
    rotateVector(cCoefs+36, cFrame, alpha);
  }
  else { // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
    alpha[0] = alpha[1] = alpha[2] = prop->W;
    alpha[3] = alpha[4] = alpha[5] = 0;
  }

  // Integate over the element: F = Int[Bt.ThermaStress]
  // with ThermalStress = C.ThermalStrain, with ThermalStrain = alpha.theta.[1, 1, 1, 0, 0, 0]'
  // where theta = T(M)-Tref = Sum[inode][N[inode]*(ndTemps[inode] - Tref)]
  // N[inode] is the shape fct at node inode
  // M is the position in the real frame, m its associated position in the reference
  // element frame
  // NUMERICAL INTEGRATION BY GAUSS PTS
  const int numgauss = 15;
  extern double dp[15][10][3]; // arrays vp1 & dp contain the values of the Tet10 shape fct & their
  extern double vp1[15][10];   // derivatives w.r.t reference coordinate system at the 15 integration pts
  extern double gauss3d5[15][3];
  double w, J;
  double DShape[10][3];

  for(int i=0; i<numgauss; i++) {
    J = computeTet10DShapeFct(dp[i],X,Y,Z,DShape);
    double* Shape = &vp1[i][0];
    w = fabs(J)*weight3d5[i];
    // compute thermal stresses
    double eT = 0.0;
    for(int inode=0; inode<nnodes; inode++) eT += Shape[inode]*(ndTemps[inode] - Tref);
    double thermalStrain[6];
    for(int l=0; l<6; ++l) thermalStrain[l] = alpha[l]*eT;
    double thermalStress[6] = {0.0,0.0,0.0,0.0,0.0,0.0}; 
    computeStress3DSolid(thermalStress, thermalStrain, C); // thermalStress <- C.thermalStrain
    // sum contribution
    for(int inode=0; inode<nnodes; inode++) {
      elementThermalForce[3*inode  ] += w*(DShape[inode][0]*thermalStress[0] + DShape[inode][1]*thermalStress[3] + DShape[inode][2]*thermalStress[5]);
      elementThermalForce[3*inode+1] += w*(DShape[inode][0]*thermalStress[3] + DShape[inode][1]*thermalStress[1] + DShape[inode][2]*thermalStress[4]);
      elementThermalForce[3*inode+2] += w*(DShape[inode][0]*thermalStress[5] + DShape[inode][1]*thermalStress[4] + DShape[inode][2]*thermalStress[2]);
    }
  }
}

FullSquareMatrix
TenNodeTetrahedral::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
  const int nnodes = 10;
  const int ndofs = 30;

  double X[10], Y[10], Z[10];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  FullSquareMatrix M(ndofs,mel);

  if(cmflg) { // consistent mass matrix
    M.zero();
    int ls[30] = {0, 3, 6, 9,12,15,18,21,24,27,
                  1, 4, 7,10,13,16,19,22,25,28,
                  2, 5, 8,11,14,17,20,23,26,29};

    // integration: loop over Gauss pts
    // reuse the 15 pts integration rule (order 5) -> reuse arrays vp1 & dp
    const int numgauss = 15;     // use 15 pts integration rule (order 5)
    extern double dp[15][10][3]; // arrays vp1 & dp contain the values of the Tet10 shape fct & their
    extern double vp1[15][10];   // derivatives w.r.t reference coordinate system at the 15 integration pts
    double dOmega; // det of jacobian
    int jSign = 0;

    for(int i = 0; i < numgauss; i++) {
      dOmega = computeTet10DShapeFct(dp[i], X, Y, Z);
#ifdef CHECK_JACOBIAN
      checkJacobian(&dOmega, &jSign, getGlNum()+1, "TenNodeTetrahedral::massMatrix");
#endif
      double w = fabs(dOmega)*weight3d5[i]*prop->rho;
      addNtDNtoM3DSolid(M, vp1[i], w, nnodes, ls);
    }
  }
  else { // Lumped mass matrix
    fprintf(stderr," *** In TenNodeTetrahedral::massMatrix: Lumped mass matrix NOT implemented. Abort.\n");
    exit(-1);
  }

  return M;
}

//HB (06/19/05)  new implementation of the Tetra4 stiffness matrix to deal
//               with anisotropic constitutive matrix
FullSquareMatrix
TenNodeTetrahedral::stiffness(const CoordSet &cs, double *d, int flg) const
{
  const int nnodes = 10;
  const int ndofs = 30;

  double X[10], Y[10], Z[10];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  int ls[30] = {0, 3, 6, 9,12,15,18,21,24,27,
                1, 4, 7,10,13,16,19,22,25,28,
                2, 5, 8,11,14,17,20,23,26,29};
  FullSquareMatrix K(ndofs,d);
  K.zero();

  // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic material
    // transform local constitutive matrix to global frame
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);

  // integration: loop over Gauss pts
  // reuse the 15 pts integration rule (order 5) -> reuse arrays vp1 & dp
  const int numgauss = 15;     // use 15 pts integration rule (order 5) (order 2 is exact for stiffness if linear mapping)
  extern double dp[15][10][3]; // arrays vp1 & dp contain the values of the Tet10 shape fct & their
  extern double vp1[15][10];   // derivatives w.r.t reference coordinate system at the 15 integration pts
  double DShape[10][3];
  double dOmega; // det of jacobian
  int jSign = 0;

  for(int i=0; i<numgauss; i++) {
    dOmega = computeTet10DShapeFct(dp[i],X,Y,Z,DShape);
#ifdef CHECK_JACOBIAN
    checkJacobian(&dOmega, &jSign, getGlNum()+1, "TenNodeTetrahedral::stiffness");
#endif
    double w = fabs(dOmega)*weight3d5[i];
    addBtCBtoK3DSolid(K, DShape, C, w, nnodes, ls);
  }

  return K;
}

void
TenNodeTetrahedral::aRubberStiffnessDerivs(CoordSet & cs, complex<double> *d,
                                            int n, double omega) {
  const int nnodes = 10;
  const int ndofs = 30;

  double X[10], Y[10], Z[10];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  int ls[30] = {0, 3, 6, 9,12,15,18,21,24,27,
                1, 4, 7,10,13,16,19,22,25,28,
                2, 5, 8,11,14,17,20,23,26,29};

  // integration: loop over Gauss pts
  // reuse the 15 pts integration rule (order 5) -> reuse arrays vp1 & dp
  const int numgauss = 15;     // use 15 pts integration rule (order 5) (order 2 is exact for stiffness if linear mapping)
  extern double dp[15][10][3]; // arrays vp1 & dp contain the values of the Tet10 shape fct & their
  extern double vp1[15][10];   // derivatives w.r.t reference coordinate system at the 15 integration pts
  double DShape[10][3];
  double dOmega; // det of jacobian
  int jSign = 0;

  double Cm[6][6], Cl[6][6];
  for(int i=0;i<6;i++) for(int j=0;j<6;j++) Cm[i][j] = 0.0;
  for(int i=0;i<6;i++) for(int j=0;j<6;j++) Cl[i][j] = 0.0;
  double *Km = new double[ndofs*ndofs];
  double *Kl = new double[ndofs*ndofs];
  FullSquareMatrix Kmr(ndofs,Km);
  FullSquareMatrix Klr(ndofs,Kl);
  Kmr.zero();
  Klr.zero();
  Cm[0][0] = Cm[1][1] = Cm[2][2] = 2.0;
  Cm[3][3] = Cm[4][4] = Cm[5][5] = 1.0;
  Cl[0][0] = Cl[1][1] = Cl[2][2] = 1.0;
  Cl[0][1] = Cl[1][0] = Cl[0][2] = Cl[2][0] = Cl[1][2] = Cl[2][1] = 1.0;
  Cl[3][3] = Cl[4][4] = Cl[5][5] = 0.0;

  for(int i=0; i<numgauss; i++) {
    dOmega = computeTet10DShapeFct(dp[i],X,Y,Z,DShape);
#ifdef CHECK_JACOBIAN
    checkJacobian(&dOmega, &jSign, getGlNum()+1, "TenNodeTetrahedral::stiffness");
#endif
    double w = fabs(dOmega)*weight3d5[i];
    addBtCBtoK3DSolid(Kmr, DShape, Cm, w, nnodes, ls);
    addBtCBtoK3DSolid(Klr, DShape, Cl, w, nnodes, ls);
  }


  ARubberF ar(n,omega,
              prop->E0,prop->dE,prop->mu0,prop->dmu,
              prop->eta_E,prop->deta_E,prop->eta_mu,prop->deta_mu);

  for(int i=0;i<ndofs*ndofs;i++)  d[i+(0)*ndofs*ndofs] = Km[i];
  for(int i=0;i<ndofs*ndofs;i++)  d[i+(1)*ndofs*ndofs] = Kl[i];
  for(int j=0;j<=n;j++)
    for(int i=0;i<ndofs*ndofs;i++)
        d[i+(j+2)*ndofs*ndofs] = ar.d_lambda(j)*Kl[i]+ ar.d_mu(j)*Km[i];

 delete[] Kl;
 delete[] Km;
}

int
TenNodeTetrahedral::numNodes() const
{
  if(useFull) return 10;
  return 4;
}

int
TenNodeTetrahedral::numDofs() const
{
  return 30;
}

int
TenNodeTetrahedral::getTopNumber() const
{
  return 125;
}

int*
TenNodeTetrahedral::nodes(int *p) const
{
  if(useFull) {
    if(!p) p = new int[10];
    p[0] = nn[0];
    p[1] = nn[1];
    p[2] = nn[2];
    p[3] = nn[3];
    p[4] = nn[4];
    p[5] = nn[5];
    p[6] = nn[6];
    p[7] = nn[7];
    p[8] = nn[8];
    p[9] = nn[9];
  }
  else {
    if(!p) p = new int[4];
    p[0] = nn[0];
    p[1] = nn[1];
    p[2] = nn[2];
    p[3] = nn[3];
  }

  return p;
}

int*
TenNodeTetrahedral::dofs(DofSetArray &dsa, int *p) const
{
  if(!p) p = new int[30];

  dsa.number(nn[0], DofSet::XYZdisp, p  );
  dsa.number(nn[1], DofSet::XYZdisp, p+3);
  dsa.number(nn[2], DofSet::XYZdisp, p+6);
  dsa.number(nn[3], DofSet::XYZdisp, p+9);
  dsa.number(nn[4], DofSet::XYZdisp, p+12);
  dsa.number(nn[5], DofSet::XYZdisp, p+15);
  dsa.number(nn[6], DofSet::XYZdisp, p+18);
  dsa.number(nn[7], DofSet::XYZdisp, p+21);
  dsa.number(nn[8], DofSet::XYZdisp, p+24);
  dsa.number(nn[9], DofSet::XYZdisp, p+27);

  return p;
}

void
TenNodeTetrahedral::markDofs(DofSetArray &dsa) const
{
  dsa.mark(nn, numNodes(), DofSet::XYZdisp);
}

/*int
TenNodeTetrahedral::facelist(PolygonSet &ps,int *list)
{
 if(list != 0) {
   list[0] = ps.addTri6(nn[0],nn[2],nn[1],nn[6],nn[5],nn[4]);
   list[1] = ps.addTri6(nn[0],nn[3],nn[2],nn[7],nn[9],nn[6]);
   list[2] = ps.addTri6(nn[0],nn[1],nn[3],nn[4],nn[8],nn[7]);
   list[3] = ps.addTri6(nn[1],nn[2],nn[3],nn[5],nn[9],nn[8]);
 }
 return 4;
}
*/

// Stress evaluation in case of anisotropic elastic constitutive matrix
void
TenNodeTetrahedral::getVonMisesAniso(Vector &stress, Vector &weight, CoordSet &cs,
                                     Vector &elDisp, int strInd, int surface, double *ndTemps,
                                     double ylayer, double zlayer, int avgnum)
{
  const int nnodes = 10;
  weight = 1.0;
  
  double X[10], Y[10], Z[10];
  cs.getCoordinates(nn, nnodes, X, Y, Z);
     
  // Flags sands17 to calculate Von Mises stress and/or Von Mises strain
  bool vmflg     = (strInd==6) ? true : false;
  bool strainFlg = (strInd==13)? true : false;
  bool meanVms   = false; // This can be used to force averaging of the Von Mises stress & strain.

  double elStress[10][7];
  double elStrain[10][7];
 
  // get constitutive matrix and coefficients of thermal expansion
  double C[6][6], alpha[6];
  // transform local constitutive matrix to global frame
  rotateConstitutiveMatrix(cCoefs, cFrame, C);
  // transform local coefficients of thermal expansion to global frame
  if(ndTemps) rotateVector(cCoefs+36, cFrame, alpha);

  // Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[10][3] = {{0.0,0.0,0.0},{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0},
                                {0.5,0.0,0.0},{0.5,0.5,0.0},{0.0,0.5,0.0},
                                {0.0,0.0,0.5},{0.5,0.0,0.5},{0.0,0.5,0.5}};

  double Shape[10], DShape[10][3];
  for(int inode=0; inode<nnodes; inode++) {
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Tetra10ShapeFct(Shape, DShape, m, X, Y, Z); 
    computeStressAndEngStrain3DSolid(elStress[inode], elStrain[inode], C, DShape, elDisp.data(), nnodes);

    if(ndTemps) {
      double &Tref = prop->Ta;
      double thermalStrain[6];
      for(int i=0; i<6; ++i) thermalStrain[i] = alpha[i]*(ndTemps[inode]-Tref);
      double thermalStress[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
      computeStress3DSolid(thermalStress, thermalStrain, C);
      elStress[inode][0] -= thermalStress[0];
      elStress[inode][1] -= thermalStress[1];
      elStress[inode][2] -= thermalStress[2];
      elStress[inode][3] -= thermalStress[3];
      elStress[inode][4] -= thermalStress[4];
      elStress[inode][5] -= thermalStress[5];
    }
    
    if(vmflg) elStress[inode][6] = computeVonMisesStress(elStress[inode]);
    else elStress[inode][6] = 0.0;
 
    if(strainFlg) elStrain[inode][6] = computeVonMisesStrain(elStrain[inode]);
    else elStrain[inode][6] = 0.0;
  }

  // compute average Von Mises stress and/or Von Mises strain: to match old Fortran code 
  if(vmflg && meanVms) {
    double vms = 0.0;
    for(int inode=0; inode<nnodes; inode++) vms += elStress[inode][6];
    vms /= nnodes;
    for(int inode=0; inode<nnodes; inode++) elStress[inode][6] = vms;
  }
  if(strainFlg && meanVms) {
    double vms = 0.0;
    for(int inode=0; inode<nnodes; inode++) vms += elStrain[inode][6];
    vms /= nnodes;
    for(int inode=0; inode<nnodes; inode++) elStrain[inode][6] = vms;
  }
  
  // fill the output array stress with the requested stress or strain component 
  for(int inode=0; inode<nnodes; inode++) {
    if(strInd < 7) 
      stress[inode] = elStress[inode][strInd];
    else if(strInd < 14)
      stress[inode] = elStrain[inode][strInd-7];
    else
      stress[inode] = 0;
  }
}

void
TenNodeTetrahedral::getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
                                      Vector &elDisp, int strInd, int surface, double *ndTemps)
{
  const int nnodes = 10;
  weight = 1.0;
  
  double X[10], Y[10], Z[10];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  double elStress[10][6];
  double elStrain[10][6];
 
  // get constitutive matrix and coefficients of thermal expansion
  double C[6][6], alpha[6];
  // transform local constitutive matrix to global frame
  rotateConstitutiveMatrix(cCoefs, cFrame, C);
  // transform local coefficients of thermal expansion to global frame
  if(ndTemps) rotateVector(cCoefs+36, cFrame, alpha);
 
  // Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[10][3] = {{0.0,0.0,0.0},{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0},
                                {0.5,0.0,0.0},{0.5,0.5,0.0},{0.0,0.5,0.0},
                                {0.0,0.0,0.5},{0.5,0.0,0.5},{0.0,0.5,0.5}};

  double Shape[10], DShape[10][3];
  for(int inode=0; inode<nnodes; inode++) {
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Tetra10ShapeFct(Shape, DShape, m, X, Y, Z); 
    computeStressAndEngStrain3DSolid(elStress[inode], elStrain[inode], C, DShape, elDisp.data(), nnodes);

    if(ndTemps) {
      double &Tref = prop->Ta;
      double thermalStrain[6];
      for(int i=0; i<6; ++i) thermalStrain[i] = alpha[i]*(ndTemps[inode]-Tref);
      double thermalStress[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
      computeStress3DSolid(thermalStress, thermalStrain, C);
      elStress[inode][0] -= thermalStress[0];
      elStress[inode][1] -= thermalStress[1];
      elStress[inode][2] -= thermalStress[2];
      elStress[inode][3] -= thermalStress[3];
      elStress[inode][4] -= thermalStress[4];
      elStress[inode][5] -= thermalStress[5];
    }
  }

  // Store all Stress or all Strain as defined by strInd
  if(strInd == 0) {
    for(int i=0; i<nnodes; ++i)
      for(int j=0; j<6; ++j) 
        stress[i][j] = elStress[i][j];
  } else {
    for(int i=0; i<nnodes; ++i)
      for(int j=0; j<6; ++j) 
        stress[i][j] = elStrain[i][j];
  }       
  
  // Get Element Principals without averaging
  double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double pvec[3] = {0.0,0.0,0.0};
  for(int i=0; i<nnodes; ++i) {
    for(int j=0; j<6; ++j) 
      svec[j] = stress[i][j];

    // Convert Engineering to Tensor Strains
    if(strInd != 0) { svec[3] /= 2; svec[4] /= 2; svec[5] /= 2; }
    pstress(svec,pvec); // compute principal stress (or strain) & direction
    for(int j=0; j<3; ++j) 
      stress[i][j+6] = pvec[j];
  }
}

void
TenNodeTetrahedral::setMaterial(NLMaterial *_mat)
{
  if(cCoefs) { // anisotropic material
    mat = _mat->clone();
    if(mat) {
      double C[6][6], alpha[6];
      // transform local constitutive matrix to global frame
      rotateConstitutiveMatrix(cCoefs, cFrame, C);
      mat->setTangentMaterial(C);
      // transform local coefficients of thermal expansion to global frame
      rotateVector(cCoefs+36, cFrame, alpha);
      mat->setThermalExpansionCoef(alpha);
    }
  }
  else {
    mat = _mat;
  }
}

int
TenNodeTetrahedral::numStates()
{
#ifdef USE_EIGEN3
  int numGaussPoints = NLTetrahedral10::numGaussPoints;
  return (mat) ? numGaussPoints*mat->getNumStates(): 0;
#else
  return 0;
#endif
}

void
TenNodeTetrahedral::initStates(double *st)
{
#ifdef USE_EIGEN3
  if(mat) {
    int ninterns = mat->getNumStates();
    int numGaussPoints = NLTetrahedral10::numGaussPoints;

    for(int i = 0; i < numGaussPoints; ++i)
      mat->initStates(st+i*ninterns);
  }
#endif
}

Corotator *
TenNodeTetrahedral::getCorotator(CoordSet &cs, double *kel, int, int)
{
  if(cCoefs && !mat) {
    double C[6][6], alpha[6];
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
    rotateVector(cCoefs+36, cFrame, alpha);
    mat = new StVenantKirchhoffMat(prop->rho, C, prop->Ta, alpha);
  }
  if(mat) {
#ifdef USE_EIGEN3
    mat->setTDProps(prop->ymtt, prop->ctett);
    MatNLElement *ele = new NLTetrahedral10(nn);
    ele->setMaterial(mat);
    ele->setGlNum(glNum);
    ele->setProp(prop);
    return new MatNLCorotator(ele);
#endif
  }
  else {
    return new Tet10Corotator(nn, prop->E, prop->nu, cs, prop->Ta, prop->W, prop->ymtt, prop->ctett);
  }
  printf("WARNING: Corotator not implemented for element %d\n", glNum+1); return 0;
}

int
TenNodeTetrahedral::getDecFace(int iFace, int *fn)
{
  switch(iFace) {
    case 0: fn[0] = nn[0]; fn[1] = nn[2]; fn[2] = nn[1]; break;
    case 1: fn[0] = nn[0]; fn[1] = nn[1]; fn[2] = nn[3]; break;
    case 2: fn[0] = nn[0]; fn[1] = nn[3]; fn[2] = nn[2]; break;
    default:
    case 3: fn[0] = nn[2]; fn[1] = nn[3]; fn[2] = nn[1]; break;
  }
  return 3;
}

int
TenNodeTetrahedral::getFace(int iFace, int *fn)
{
  switch(iFace) {
    case 0: fn[0] = nn[0]; fn[1] = nn[2]; fn[2] = nn[1];
            fn[3] = nn[6]; fn[4] = nn[5]; fn[5] = nn[4]; break;
    case 1: fn[0] = nn[0]; fn[1] = nn[1]; fn[2] = nn[3];
            fn[3] = nn[4]; fn[4] = nn[8]; fn[5] = nn[7]; break;
    case 2: fn[0] = nn[0]; fn[1] = nn[3]; fn[2] = nn[2];
            fn[3] = nn[7]; fn[4] = nn[9]; fn[5] = nn[6]; break;
    default:
    case 3: fn[0] = nn[2]; fn[1] = nn[3]; fn[2] = nn[1];
            fn[3] = nn[9]; fn[4] = nn[8]; fn[5] = nn[5]; break;
  }
  return 6;
}

void
TenNodeTetrahedral::getCFrame(CoordSet &cs, double cFrame[3][3]) const
{
  if(TenNodeTetrahedral::cFrame) {
    cFrame[0][0] = TenNodeTetrahedral::cFrame[0]; cFrame[0][1] = TenNodeTetrahedral::cFrame[1]; cFrame[0][2] = TenNodeTetrahedral::cFrame[2];
    cFrame[1][0] = TenNodeTetrahedral::cFrame[3]; cFrame[1][1] = TenNodeTetrahedral::cFrame[4]; cFrame[1][2] = TenNodeTetrahedral::cFrame[5];
    cFrame[2][0] = TenNodeTetrahedral::cFrame[6]; cFrame[2][1] = TenNodeTetrahedral::cFrame[7]; cFrame[2][2] = TenNodeTetrahedral::cFrame[8];
  }
  else {
    cFrame[0][0] = cFrame[1][1] = cFrame[2][2] = 1.;
    cFrame[0][1] = cFrame[0][2] = cFrame[1][0] = cFrame[1][2] = cFrame[2][0] = cFrame[2][1] = 0.;
  }
}

