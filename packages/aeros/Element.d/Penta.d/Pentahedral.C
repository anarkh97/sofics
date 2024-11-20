#include <cstdio>
#include <iostream>
#include <cmath>

#include <Element.d/Penta.d/Pentahedral.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Corotational.d/PentaCorotator.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/NLPentahedral.h>
#include <Element.d/Utils.d/SolidElemUtils.h>
#include <Corotational.d/MatNLCorotator.h>

#define CHECK_JACOBIAN // force check nullity & constant sign of jacobian over el.

extern "C" {
void _FORTRAN(mstf24)(double*, double*, double*, double&, double&, double*, 
                      const int&, int&);
void _FORTRAN(sands24)(const int&, double*, double*, double*, double&, 
                       double&, double*, double*, double*, const int&, const int&,
                       const int&, const int&, const int&,const int&);
void  _FORTRAN(lgauss)(int &, int &, double *, double *); 
void _FORTRAN(brkcmt)(double&, double&, double*);
}

double Penta6ShapeFct(double Shape[6], double DShape[6][3], double m[3], double X[6], double Y[6], double Z[6]);

Pentahedral::Pentahedral(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
  nn[2] = nodenums[2];
  nn[3] = nodenums[3];
  nn[4] = nodenums[4];
  nn[5] = nodenums[5];

  cFrame = 0;
  cCoefs = 0;
  mat = 0;
}

Pentahedral::~Pentahedral()
{
  if(cCoefs && mat) delete mat;
}

Element *
Pentahedral::clone()
{
  return new Pentahedral(*this);
}

void
Pentahedral::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
  nn[5] = table[nn[5]];
}

void
Pentahedral::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
  nn[5] = table[nn[5]];
}

void
Pentahedral::getVonMises(Vector& stress, Vector& weight, CoordSet &cs,
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

  double x[6], y[6], z[6];
  cs.getCoordinates(nn, numNodes(), x, y, z);

  int vmflg = 0, strainFlg = 0;
  // Flags sands24 to calculate Von Mises stress
  if(strInd == 6) vmflg = 1;
  // Flags sands24 to calculate Von Mises strain
  if(strInd == 13) strainFlg = 1;

  const int maxgus = 6;
  const int maxstr = 7;
  const int elm    = 1;
  const int msize  = 1;
  const int outerr = 6;

  double elStress[maxgus][maxstr], elStrain[maxgus][maxstr];

  _FORTRAN(sands24)(elm, x, y, z, prop->E, prop->nu, elDisp.data(),
                    (double*)elStress, (double*)elStrain,
                    maxgus, maxstr, msize, outerr, vmflg, strainFlg);

  if(strInd < 7) {
    double thermalStress[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    if(strInd == 0 || strInd == 1 || strInd == 2) {
      double &Tref  = prop->Ta;
      double &alpha = prop->W;
      double coef = (prop->E*alpha)/(1.0 - 2.0*prop->nu);
      thermalStress[0] = coef*(ndTemps[0]-Tref);
      thermalStress[1] = coef*(ndTemps[1]-Tref);
      thermalStress[2] = coef*(ndTemps[2]-Tref);
      thermalStress[3] = coef*(ndTemps[3]-Tref);
      thermalStress[4] = coef*(ndTemps[4]-Tref);
      thermalStress[5] = coef*(ndTemps[5]-Tref);
    }
    stress[0] = elStress[0][strInd] - thermalStress[0];
    stress[1] = elStress[1][strInd] - thermalStress[1];
    stress[2] = elStress[2][strInd] - thermalStress[2];
    stress[3] = elStress[3][strInd] - thermalStress[3];
    stress[4] = elStress[4][strInd] - thermalStress[4];
    stress[5] = elStress[5][strInd] - thermalStress[5];
  }
  else if(strInd < 14) {
    stress[0] = elStrain[0][strInd-7];
    stress[1] = elStrain[1][strInd-7];
    stress[2] = elStrain[2][strInd-7];
    stress[3] = elStrain[3][strInd-7];
    stress[4] = elStrain[4][strInd-7];
    stress[5] = elStrain[5][strInd-7];
  }
  else {
    stress[0] = 0;
    stress[1] = 0;
    stress[2] = 0;
    stress[3] = 0;
    stress[4] = 0;
    stress[5] = 0;
  }
}

void
Pentahedral::getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
                          Vector &elDisp, int strInd, int surface, double *ndTemps)
{
  if(cCoefs) {
    getAllStressAniso(stress, weight, cs,
                      elDisp, strInd, surface, ndTemps);
    return;
  }

  weight = 1.0;

  double x[6], y[6], z[6];
  cs.getCoordinates(nn, numNodes(), x, y, z);

  int vmflg, strainFlg;
  vmflg = 0;
  strainFlg = 0;

  const int maxgus = 6;
  const int maxstr = 7;
  const int elm    = 1;
  const int msize  = 1;
  const int outerr = 6;

  double elStress[maxgus][maxstr], elStrain[maxgus][maxstr];

  _FORTRAN(sands24)(elm, x, y, z, prop->E, prop->nu, elDisp.data(),
                    (double*)elStress, (double*)elStrain,
                    maxgus, maxstr, msize, outerr, vmflg, strainFlg);

  // Store all Stress or all Strain as defined by strInd
  int i,j;
  if(strInd == 0) {
    double thermalStress[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    if(ndTemps) {
      double &Tref  = prop->Ta;
      double &alpha = prop->W;
      double coef = (prop->E*alpha)/(1.0 - 2.0*prop->nu);
      thermalStress[0] = coef*(ndTemps[0]-Tref);
      thermalStress[1] = coef*(ndTemps[1]-Tref);
      thermalStress[2] = coef*(ndTemps[2]-Tref);
      thermalStress[3] = coef*(ndTemps[3]-Tref);
      thermalStress[4] = coef*(ndTemps[4]-Tref);
      thermalStress[5] = coef*(ndTemps[5]-Tref);
    }
    for (i=0; i<6; ++i) {
      for (j=0; j<3; ++j) {
        stress[i][j] = elStress[i][j] - thermalStress[i];
      }
      for (j=3; j<6; ++j) {
        stress[i][j] = elStress[i][j];
      }
    }
  }
  else {
    for (i=0; i<6; ++i) {
      for (j=0; j<6; ++j) {
        stress[i][j] = elStrain[i][j];
      }
    }
  }

  // Get Element Principals for each node without averaging
  double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double pvec[3] = {0.0,0.0,0.0};

  for (i=0; i<6; ++i) {
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
Pentahedral::getMass(const CoordSet& cs) const
{
  const int nnodes = 6;
  double X[6], Y[6], Z[6];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  // integration: loop over Gauss pts  
  // hard coded order 2 triangle quadrature rule: {r,s,t(=1-r-s),w}
  double TriGPt3[3][4] = {{1./6.,1./6.,1./6.,1./6.},
                          {2./3.,1./6.,1./6.,1./6.},
                          {1./6.,2./3.,1./6.,1./6.}};
  double wxy, wz;
  int ngpz = 2;  // number of (linear) Gauss pts in the (local) z direction
  int ngpxy = 3; // number of (triangular) integration pts (in the local x-y plane)
  double m[3], Shape[6], DShape[6][3];
  double dOmega;
  double volume = 0.0;

  for(int iz = 1; iz <= ngpz; iz++) {
    // get z position & weight of the Gauss pt
    _FORTRAN(lgauss)(ngpz, iz, &m[2], &wz);
    for(int ixy = 0; ixy < ngpxy; ixy++) { // triangle Gauss pts
      // get x, y position & weight of the Gauss pt
      m[0] = TriGPt3[ixy][0]; m[1] = TriGPt3[ixy][1]; wxy = TriGPt3[ixy][3];
      dOmega = Penta6ShapeFct(Shape, DShape, m, X, Y, Z);
      volume += fabs(dOmega)*wxy*wz;
    }
  }

  return volume*prop->rho;
}

void
Pentahedral::getGravityForce(CoordSet& cs, double *gravityAcceleration,
                             Vector& gravityForce, int gravflg, GeomState *geomState)
{
  const int nnodes = 6;

  // Lumped
  if (gravflg != 2) {

    double totmas = getMass(cs);

    double m[18*18];
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

    const int ndofs = 18;

    double lforce[20];
    for(int i=0; i<nnodes; ++i) lforce[i] = 0.0;

    double X[6], Y[6], Z[6];
    cs.getCoordinates(nn, nnodes, X, Y, Z);

    // integration: loop over Gauss pts
    // hard coded order 2 triangle quadrature rule: {r,s,t(=1-r-s),w}
    double TriGPt3[3][4] = {{1./6.,1./6.,1./6.,1./6.},
                            {2./3.,1./6.,1./6.,1./6.},
                            {1./6.,2./3.,1./6.,1./6.}};
    double wxy, wz, w;
    int ngpz = 2; // number of (linear) Gauss pts in the (local) z direction
    int ngpxy = 3; // number of (triangular) integration pts (in the local x-y plane)
    double m[3], Shape[6], DShape[6][3];
    double dOmega;

    for(int iz = 1; iz <= ngpz; iz++) {
      // get z position & weight of the Gauss pt
      _FORTRAN(lgauss)(ngpz, iz, &m[2], &wz);
      for(int ixy = 0; ixy < ngpxy; ixy++) { // triangle Gauss pts
        // get x, y  position & weight of the Gauss pt
        m[0] = TriGPt3[ixy][0]; m[1] = TriGPt3[ixy][1]; wxy = TriGPt3[ixy][3];
        dOmega = Penta6ShapeFct(Shape, DShape, m, X, Y, Z);
        w = fabs(dOmega)*wxy*wz*prop->rho;

        for(int n = 0; n < nnodes; ++n)
          lforce[n] += w*Shape[n];
      }
    }

    for(int i = 0; i < nnodes; ++i) {
      gravityForce[3*i+0] = lforce[i]*gravityAcceleration[0];
      gravityForce[3*i+1] = lforce[i]*gravityAcceleration[1];
      gravityForce[3*i+2] = lforce[i]*gravityAcceleration[2];
    }
  }
}

void
Pentahedral::getThermalForce(CoordSet &cs, Vector &ndTemps,
                             Vector &elementThermalForce, int glflag,
                             GeomState *geomState)
{
  // ASSUME CONSTANT THERMAL EXPANSION COEFF. & REFERENCE TEMPERATURE OVER THE ELEMENT 
  const int nnodes = 6;
  const int ndofs = 18;

  // initialize nodal thermal forces
  for(int i=0; i<ndofs; i++) elementThermalForce[i] = 0.0;

  // for nonlinear analyses, the thermal load for this element is now computed in getStiffAndForce
  if(geomState) return;

  // extract nodes coordinates
  double X[6], Y[6], Z[6];
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
  // USE NUMERICAL INTEGRATION BY GAUSS PTS
  double Shape[6], DShape[6][3], m[3];
  double wxy,wz,w,J;
  int ngpz  = 2; // number of (linear) Gauss pts in the (local) z direction
  int ngpxy = 3; // number of (triangular) integration pts (in the local x-y plane)
  // hard coded order 2 triangle quadrature rule: {r,s,t(=1-r-s),w}
  double TriGPt3[3][4] = {{1./6.,1./6.,1./6.,1./6.},
                          {2./3.,1./6.,1./6.,1./6.},
                          {1./6.,2./3.,1./6.,1./6.}};

  // integration: loop over Gauss pts
  for(int iz=1; iz<=ngpz; iz++) { // z Gauss pts
    // get z position & weight of the Gauss pt
    _FORTRAN(lgauss)(ngpz,iz,&m[2],&wz);
    for(int ixy=0; ixy<ngpxy; ixy++) { // triangle Gauss pts
      // get x, y  position & weight of the Gauss pt
      m[0] = TriGPt3[ixy][0]; m[1] = TriGPt3[ixy][1]; wxy = TriGPt3[ixy][3];
      // compute shape fcts & their derivatives at the Gauss pt
      J = Penta6ShapeFct(Shape, DShape, m, X, Y, Z);
      w = wxy*wz*fabs(J);
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
}

FullSquareMatrix
Pentahedral::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
  const int nnodes = 6;
  const int ndofs = 18;

  double X[6], Y[6], Z[6];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  FullSquareMatrix M(ndofs,mel);

  if(cmflg) { // consistent mass matrix
    M.zero();
    int ls[18] = {0,3,6, 9,12,15,
                  1,4,7,10,13,16,
                  2,5,8,11,14,17};

    // integration: loop over Gauss pts
    // hard coded order 2 triangle quadrature rule: {r,s,t(=1-r-s),w}
    double TriGPt3[3][4] = {{1./6.,1./6.,1./6.,1./6.},
                            {2./3.,1./6.,1./6.,1./6.},
                            {1./6.,2./3.,1./6.,1./6.}};
    double wxy, wz, w;
    int ngpz = 2; // number of (linear) Gauss pts in the (local) z direction
    int ngpxy= 3; // number of (triangular) integration pts (in the local x-y plane)
    double m[3], Shape[6], DShape[6][3];
    double dOmega;
    int jSign = 0;

    for(int iz = 1; iz <= ngpz; iz++) {
      // get z position & weight of the Gauss pt
      _FORTRAN(lgauss)(ngpz, iz, &m[2], &wz);
      for(int ixy = 0; ixy < ngpxy; ixy++) { // triangle Gauss pts
        // get x, y position & weight of the Gauss pt
        m[0] = TriGPt3[ixy][0]; m[1] = TriGPt3[ixy][1]; wxy = TriGPt3[ixy][3];
        dOmega = Penta6ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
        checkJacobian(&dOmega, &jSign, getGlNum()+1, "Pentahedral::massMatrix");
#endif
        w = fabs(dOmega)*wxy*wz*prop->rho;
        addNtDNtoM3DSolid(M, Shape, w, nnodes, ls);
      }
    }
  }
  else { // Lumped mass matrix
    fprintf(stderr," *** In Pentahedral::massMatrix: Lumped mass matrix NOT implemented. Abort.\n");
    exit(-1);
  }

  return M;
}

//HB (04/15/05)  new implementation of the Penta6 stiffness matrix to deal
//               with anisotropic constitutive matrix
FullSquareMatrix
Pentahedral::stiffness(const CoordSet &cs,double *d, int flg) const
{
  const int nnodes = 6;

  double X[6], Y[6], Z[6];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  int ls[18] = {0,3,6, 9,12,15,
                1,4,7,10,13,16,
                2,5,8,11,14,17};
  FullSquareMatrix K(18,d);
  K.zero();

  // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic material
    // transform local constitutive matrix to global frame
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);

  // integration: loop over Gauss pts
  // hard coded order 2 triangle quadrature rule: {r,s,t(=1-r-s),w}
  double TriGPt3[3][4] = {{1./6.,1./6.,1./6.,1./6.},
                          {2./3.,1./6.,1./6.,1./6.},
                          {1./6.,2./3.,1./6.,1./6.}};
  double wxy,wz,w;
  int ngpz = 2; // number of (linear) Gauss pts in the (local) z direction
  int ngpxy= 3; // number of (triangular) integration pts (in the local x-y plane)
  double m[3], Shape[6], DShape[6][3];
  double dOmega;
  int jSign = 0;

  for(int iz=1;iz<=ngpz;iz++) { 
    // get z position & weight of the Gauss pt
   _FORTRAN(lgauss)(ngpz,iz,&m[2],&wz);
    for(int ixy=0; ixy<ngpxy; ixy++) { // triangle Gauss pts
      // get x, y position & weight of the Gauss pt
      m[0] = TriGPt3[ixy][0]; m[1] = TriGPt3[ixy][1]; wxy = TriGPt3[ixy][3];
      dOmega = Penta6ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
      checkJacobian(&dOmega, &jSign, getGlNum()+1, "Pentahedral::stiffness");
#endif
      w = fabs(dOmega)*wxy*wz;
      addBtCBtoK3DSolid(K, DShape, C, w, nnodes, ls);
    }
  }

  return K;
}

int
Pentahedral::numNodes() const
{
  return 6;
}

int
Pentahedral::numDofs() const
{
  return 18;
}

int
Pentahedral::getTopNumber() const
{
  return 124;
}

int*
Pentahedral::nodes(int *p) const
{
  if(!p) p = new int[numNodes()];
  p[0] = nn[0];
  p[1] = nn[1];
  p[2] = nn[2];
  p[3] = nn[3];
  p[4] = nn[4];
  p[5] = nn[5];
  return p;
}

int*
Pentahedral::dofs(DofSetArray &dsa, int *p) const
{
  if(!p) p = new int[numDofs()];

  dsa.number(nn[0], DofSet::XYZdisp, p  );
  dsa.number(nn[1], DofSet::XYZdisp, p+3);
  dsa.number(nn[2], DofSet::XYZdisp, p+6);
  dsa.number(nn[3], DofSet::XYZdisp, p+9);
  dsa.number(nn[4], DofSet::XYZdisp, p+12);
  dsa.number(nn[5], DofSet::XYZdisp, p+15);

  return p;
}

void
Pentahedral::markDofs(DofSetArray &dsa) const
{
  dsa.mark(nn, numNodes(), DofSet::XYZdisp);
}

// Stress evaluation in case of anisotropic elastic constitutive matrix
void
Pentahedral::getVonMisesAniso(Vector &stress, Vector &weight, CoordSet &cs,
                              Vector &elDisp, int strInd, int surface, double *ndTemps,
                              double ylayer, double zlayer, int avgnum)
{
  const int nnodes = 6;
  weight = 1.0;
  
  double X[6], Y[6], Z[6];
  cs.getCoordinates(nn, nnodes, X, Y, Z);
     
  // Flags sands17 to calculate Von Mises stress and/or Von Mises strain
  bool vmflg     = (strInd==6) ? true : false;
  bool strainFlg = (strInd==13)? true : false;
  bool meanVms   = false; // This can be used to force averaging of the Von Mises stress & strain.

  double elStress[6][7];
  double elStrain[6][7];
 
  // get constitutive matrix and coefficients of thermal expansion
  double C[6][6], alpha[6];
  // transform local constitutive matrix to global frame
  rotateConstitutiveMatrix(cCoefs, cFrame, C);
  // transform local coefficients of thermal expansion to global frame
  if(ndTemps) rotateVector(cCoefs+36, cFrame, alpha);
 
  // Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[6][3] = {{0.0,0.0,-1.0},{1.0,0.0,-1.0},{0.0,1.0,-1.0},
                               {0.0,0.0, 1.0},{1.0,0.0, 1.0},{0.0,1.0, 1.0}};

  double Shape[6], DShape[6][3];
  for(int inode=0; inode<nnodes; inode++) {
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Penta6ShapeFct(Shape, DShape, m, X, Y, Z); 
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
Pentahedral::getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
                               Vector &elDisp, int strInd, int surface, double *ndTemps)
{
  const int nnodes = 6;
  weight = 1.0;
  
  double X[6], Y[6], Z[6];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  double elStress[6][6];
  double elStrain[6][6];
 
  // get constitutive matrix and coefficients of thermal expansion
  double C[6][6], alpha[6];
  // transform local constitutive matrix to global frame
  rotateConstitutiveMatrix(cCoefs, cFrame, C);
  // transform local coefficients of thermal expansion to global frame
  if(ndTemps) rotateVector(cCoefs+36, cFrame, alpha);
 
  // Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[6][3] = {{0.0,0.0,-1.0},{1.0,0.0,-1.0},{0.0,1.0,-1.0},
                               {0.0,0.0, 1.0},{1.0,0.0, 1.0},{0.0,1.0, 1.0}};

  double Shape[6], DShape[6][3];
  for(int inode=0; inode<nnodes; inode++) {
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Penta6ShapeFct(Shape, DShape, m, X, Y, Z); 
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
Pentahedral::setMaterial(NLMaterial *_mat)
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
Pentahedral::numStates()
{
#ifdef USE_EIGEN3
  int numGaussPoints = NLPentahedral6::numGaussPoints;
  return (mat) ? numGaussPoints*mat->getNumStates(): 0;
#else
  return 0;
#endif
}

void
Pentahedral::initStates(double *st)
{
#ifdef USE_EIGEN3
  if(mat) {
    int ninterns = mat->getNumStates();
    int numGaussPoints = NLPentahedral6::numGaussPoints;

    for(int i = 0; i < numGaussPoints; ++i)
      mat->initStates(st+i*ninterns);
  }
#endif
}

Corotator *
Pentahedral::getCorotator(CoordSet &cs, double *kel, int, int)
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
    MatNLElement *ele = new NLPentahedral6(nn);
    ele->setMaterial(mat);
    ele->setGlNum(glNum);
    ele->setProp(prop);
    return new MatNLCorotator(ele);
#endif
  }
  else {
    return new PentaCorotator(nn, prop->E, prop->nu, cs, prop->Ta, prop->W, prop->ymtt, prop->ctett);
  }
  printf("WARNING: Corotator not implemented for element %d\n", glNum+1); return 0;
}

int
Pentahedral::getDecFace(int iFace, int *fn)
{
  int count;
  switch(iFace) {
    case 0: fn[0] = nn[0]; fn[1] = nn[2]; fn[2] = nn[1]; count = 3; break;
    case 1: fn[0] = nn[3]; fn[1] = nn[4]; fn[2] = nn[5]; count = 3; break;
    case 2: fn[0] = nn[0]; fn[1] = nn[1]; fn[2] = nn[4]; fn[3] = nn[3]; count = 4; break;
    case 3: fn[0] = nn[1]; fn[1] = nn[2]; fn[2] = nn[5]; fn[3] = nn[4]; count = 4; break;
    default:
    case 4: fn[0] = nn[2]; fn[1] = nn[0]; fn[2] = nn[3]; fn[3] = nn[5]; count = 4; break;
  }
  return count;
}

void
Pentahedral::getCFrame(CoordSet &cs, double cFrame[3][3]) const
{
  if(Pentahedral::cFrame) {
    cFrame[0][0] = Pentahedral::cFrame[0]; cFrame[0][1] = Pentahedral::cFrame[1]; cFrame[0][2] = Pentahedral::cFrame[2];
    cFrame[1][0] = Pentahedral::cFrame[3]; cFrame[1][1] = Pentahedral::cFrame[4]; cFrame[1][2] = Pentahedral::cFrame[5];
    cFrame[2][0] = Pentahedral::cFrame[6]; cFrame[2][1] = Pentahedral::cFrame[7]; cFrame[2][2] = Pentahedral::cFrame[8];
  }
  else {
    cFrame[0][0] = cFrame[1][1] = cFrame[2][2] = 1.;
    cFrame[0][1] = cFrame[0][2] = cFrame[1][0] = cFrame[1][2] = cFrame[2][0] = cFrame[2][1] = 0.;
  }
}
