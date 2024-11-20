// ---------------------------------------------------------------------
// HB - 05-22-05
// ---------------------------------------------------------------------
// 32 nodes brick element 
// Serendipity finite element basis
// Iso-parametric formulation
// ---------------------------------------------------------------------
#include <cstdio>
#include <iostream>
#include <Element.d/Brick32.d/Brick32.h>
#include <Element.d/Brick20.d/Brick20.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Corotational.d/Brick32Corotator.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/NLHexahedral.h>
#include <Element.d/Utils.d/SolidElemUtils.h>
#include <Corotational.d/MatNLCorotator.h>
#include <Element.d/Helm.d/ARubberF.h>

#define CHECK_JACOBIAN // force check nullity & constant sign of jacobian over el.

extern "C" {
void _FORTRAN(brkcmt)(double&, double&, double*);

void _FORTRAN(lgauss)(const int&, int&, double*, double*);
}

double Hexa32ShapeFct(double Shape[32], double DShape[32][3], double m[3], double X[32], double Y[32], double Z[32]);
double computeHexa32DShapeFct(double dShape[32][3], double X[32], double Y[32], double Z[32], double (*DShape)[3] = 0);

extern bool useFull;

//-----------------------------------------------------------------------
//    Brick32 node numbering
//                     24     23
//                5+----+-----+----+8
//                /|              /|
//             17+ |           22+ |
//              /  +26          /  +32
//           18+   |         21+   |
//            /  19|    20    /    |
//          6+----++25--+----+7    +31
//           |     |         |     |
//           |     |         |     |
//         28+    1+----+--30++----+4
//           |    /    16    |15  /
//           |   +9          |   +14
//         27+  /          29+  /
//           | +10           | +13
//           |/              |/
//          2+----+-----+----+3
//               11     12
//-----------------------------------------------------------------------

Brick32::Brick32(int* nodenums)
{
  for(int i=0; i<32; i++)
    nn[i] = nodenums[i];

  cFrame = 0;
  cCoefs = 0;
  mat = 0;
}

Brick32::~Brick32()
{
  if(cCoefs && mat) delete mat;
}

Element *
Brick32::clone()
{
  return new Brick32(*this);
}

void
Brick32::renum(const int *table)
{
  for(int i=0; i<32; i++)
    nn[i] = table[nn[i]];
}

void
Brick32::renum(EleRenumMap& table)
{
  for(int i=0; i<32; i++)
    nn[i] = table[nn[i]];
}

void
Brick32::getVonMises(Vector& stress, Vector& weight, CoordSet &cs,
                     Vector& elDisp, int strInd, int surface, double *ndTemps,
                     double ylayer, double zlayer, int avgnum)
{
  const int nnodes = 32;
  weight = 1.0;

  double X[32], Y[32], Z[32];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  // Flags to calculate Von Mises stress and/or Von Mises strain
  bool vmflg     = (strInd==6) ? true : false;
  bool strainFlg = (strInd==13)? true : false;
  bool meanVms   = false; // This can be used to force averaging of the Von Mises stress & strain.

  double elStress[32][7];
  double elStrain[32][7];

  // get constitutive matrix and coefficients of thermal expansion
  double C[6][6], alpha[6];
  if(cCoefs) { // anisotropic material
    // transform local constitutive matrix to global frame
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
    // transform local coefficients of thermal expansion to global frame
    if(ndTemps) rotateVector(cCoefs+36, cFrame, alpha);
  }
  else { // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
    alpha[0] = alpha[1] = alpha[2] = prop->W;
    alpha[3] = alpha[4] = alpha[5] = 0;
  }
 
  // Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[32][3] = {{-1.   ,-1.,-1.},{ 1.   ,-1.,-1.},{ 1., 1.   ,-1.},{-1., 1.   ,-1.},
                                {-1.   ,-1., 1.},{ 1.   ,-1., 1.},{ 1., 1.   , 1.},{-1., 1.   , 1.},
                                {-1./3.,-1.,-1.},{ 1./3.,-1.,-1.},{ 1.,-1./3.,-1.},{ 1., 1./3.,-1.},
                                { 1./3., 1.,-1.},{-1./3., 1.,-1.},{-1., 1./3.,-1.},{-1.,-1./3.,-1.},
                                {-1./3.,-1., 1.},{ 1./3.,-1., 1.},{ 1.,-1./3., 1.},{ 1., 1./3., 1.},
                                { 1./3., 1., 1.},{-1./3., 1., 1.},{-1., 1./3., 1.},{-1.,-1./3., 1.},
                                {-1.,-1.,-1./3.},{-1.,-1., 1./3.},{ 1.,-1.,-1./3.},{ 1.,-1., 1./3.},
                                { 1., 1.,-1./3.},{ 1., 1., 1./3.},{-1., 1.,-1./3.},{-1., 1., 1./3.}};

  double Shape[32], DShape[32][3];
  for(int inode=0; inode<nnodes; inode++) {
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Hexa32ShapeFct(Shape, DShape, m, X, Y, Z); 
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
Brick32::getAllStress(FullM& stress, Vector& weight, CoordSet &cs,
                      Vector& elDisp, int strInd,int surface, double *ndTemps)
{
  const int nnodes = 32;
  weight = 1.0;
  
  double X[32], Y[32], Z[32];
  cs.getCoordinates(nn, nnodes, X, Y, Z);
      
  double elStress[32][7];
  double elStrain[32][7];

  // get constitutive matrix and coefficients of thermal expansion
  double C[6][6], alpha[6];
  if(cCoefs) { // anisotropic material
    // transform local constitutive matrix to global frame
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
    // transform local coefficients of thermal expansion to global frame
    if(ndTemps) rotateVector(cCoefs+36, cFrame, alpha);
  }
  else { // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);
    alpha[0] = alpha[1] = alpha[2] = prop->W;
    alpha[3] = alpha[4] = alpha[5] = 0;
  }

  // Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[32][3] = {{-1.   ,-1.,-1.},{ 1.   ,-1.,-1.},{ 1., 1.   ,-1.},{-1., 1.   ,-1.},
                                {-1.   ,-1., 1.},{ 1.   ,-1., 1.},{ 1., 1.   , 1.},{-1., 1.   , 1.},
                                {-1./3.,-1.,-1.},{ 1./3.,-1.,-1.},{ 1.,-1./3.,-1.},{ 1., 1./3.,-1.},
                                { 1./3., 1.,-1.},{-1./3., 1.,-1.},{-1., 1./3.,-1.},{-1.,-1./3.,-1.},
                                {-1./3.,-1., 1.},{ 1./3.,-1., 1.},{ 1.,-1./3., 1.},{ 1., 1./3., 1.},
                                { 1./3., 1., 1.},{-1./3., 1., 1.},{-1., 1./3., 1.},{-1.,-1./3., 1.},
                                {-1.,-1.,-1./3.},{-1.,-1., 1./3.},{ 1.,-1.,-1./3.},{ 1.,-1., 1./3.},
                                { 1., 1.,-1./3.},{ 1., 1., 1./3.},{-1., 1.,-1./3.},{-1., 1., 1./3.}};

  double Shape[32], DShape[32][3];
  for(int inode=0; inode<nnodes; inode++) {
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Hexa32ShapeFct(Shape, DShape, m, X, Y, Z); 
    computeStressAndEngStrain3DSolid(elStress[inode], elStrain[inode], C, DShape, elDisp.data(), nnodes);

    if(ndTemps) {
      double Tref = prop->Ta;
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
  
  // Get Element Principals for each node without averaging
  double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double pvec[3] = {0.0,0.0,0.0};

  for(int i=0; i<nnodes; ++i) {
    for(int j=0; j<6; ++j) {
      svec[j] = stress[i][j];
    }
    // Convert Engineering to Tensor Strains
    if(strInd != 0) {
      svec[3] /= 2;
      svec[4] /= 2;
      svec[5] /= 2;
    }
    pstress(svec,pvec);
    for(int j=0; j<3; ++j) {
      stress[i][j+6] = pvec[j];
    }
  }
}

double
Brick32::getMass(const CoordSet& cs) const
{
  const int nnodes = 32;
  double x[32], y[32], z[32];
  cs.getCoordinates(nn, nnodes, x, y, z);

  // integration: loop over gauss pts
  const int numgauss = 4;
  double wx, wy, wz;
  double m[3], Shape[32], DShape[32][3];
  double dOmega;
  double volume = 0.0;

  for(int i = 1; i <= numgauss; i++) {
    _FORTRAN(lgauss)(numgauss, i, &m[0], &wx);
    for(int j = 1; j <= numgauss; j++) {
      _FORTRAN(lgauss)(numgauss, j, &m[1], &wy);
      for(int k = 1; k <= numgauss; k++) {
        _FORTRAN(lgauss)(numgauss, k, &m[2], &wz);
        dOmega = Hexa32ShapeFct(Shape, DShape, m, x, y, z);
        volume += fabs(dOmega)*wx*wy*wz;
      }
    }
  }

  return volume*prop->rho;
}

void
Brick32::getGravityForce(CoordSet& cs, double *gravityAcceleration, 
                         Vector& gravityForce, int gravflg, GeomState *geomState)
{
  int nnodes = 32;

  // Lumped
  if (gravflg != 2) {

    double totmas = getMass(cs);

    double m[96*96];
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
    double x[32], y[32], z[32];
    cs.getCoordinates(nn, numNodes(), x, y, z);

    double lforce[32];
    for(int i = 0; i < nnodes; ++i) lforce[i] = 0.0;

    // integration: loop over Gauss pts
    int numgauss = 4;
    double wx, wy, wz, w;
    double m[3], Shape[32], DShape[32][3];
    double dOmega; // det of jacobian

    for(int i = 1; i <= numgauss; i++) {
      _FORTRAN(lgauss)(numgauss, i, &m[0], &wx);
      for(int j = 1; j <= numgauss; j++) {
        _FORTRAN(lgauss)(numgauss, j, &m[1], &wy);
        for(int k = 1; k <= numgauss; k++) {
          _FORTRAN(lgauss)(numgauss, k, &m[2], &wz);

          dOmega = Hexa32ShapeFct(Shape, DShape, m, x, y, z);
          w = fabs(dOmega)*wx*wy*wz*prop->rho;

          for(int n = 0; n < nnodes; ++n)
            lforce[n] += w*Shape[n];
        }
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
Brick32::getThermalForce(CoordSet &cs, Vector &ndTemps,
                         Vector &elementThermalForce, int glflag,
                         GeomState *geomState)
{
  // ASSUME CONSTANT THERMAL EXPANSION COEFF. & REFERENCE TEMPERATURE OVER THE ELEMENT 
  const int nnodes = 32;
  const int ndofs = 96;

  // initialize nodal thermal forces
  for(int i=0; i<ndofs; i++) elementThermalForce[i] = 0.0;

  // for nonlinear analyses, the thermal load for this element is now computed in getStiffAndForce
  if(geomState) return;

  double X[32], Y[32], Z[32];
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
  const int numgauss = 4;
  double Shape[32], DShape[32][3], m[3];
  double wx,wy,wz,w,J;

  for(int i=1; i<=numgauss; i++) {
    _FORTRAN(lgauss)(numgauss,i,&m[0],&wx);
    for(int j=1; j<=numgauss; j++) {
      _FORTRAN(lgauss)(numgauss,j,&m[1],&wy);
      for(int k=1; k<=numgauss; k++) {
        _FORTRAN(lgauss)(numgauss,k,&m[2],&wz);
        J = Hexa32ShapeFct(Shape, DShape, m, X, Y, Z);
        w = fabs(J)*wx*wy*wz;
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
}

FullSquareMatrix
Brick32::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
  const int nnodes = 32;
  const int ndofs = 96;

  FullSquareMatrix M(ndofs, mel);

  if(cmflg) { // consistent mass matrix

    double X[32], Y[32], Z[32];
    cs.getCoordinates(nn, nnodes, X, Y, Z);

    M.zero();
    int ls[96] = {0,3,6, 9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,
                  1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67,70,73,76,79,82,85,88,91,94,
                  2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92,95};

    // integration: loop over Gauss pts
    const int numgauss = 4;
    double wx,wy,wz,w;
    double m[3], Shape[32], DShape[32][3];
    double dOmega;
    int jSign = 0;
    for(int i=1; i<=numgauss; i++) {
      _FORTRAN(lgauss)(numgauss,i,&m[0],&wx);
      for(int j=1; j<=numgauss; j++) {
        _FORTRAN(lgauss)(numgauss,j,&m[1],&wy);
        for(int k=1; k<=numgauss; k++) {
          _FORTRAN(lgauss)(numgauss,k,&m[2],&wz);
          dOmega = Hexa32ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
          checkJacobian(&dOmega, &jSign, getGlNum()+1, "Brick32::massMatrix");
#endif
          w = fabs(dOmega)*wx*wy*wz*prop->rho;
          addNtDNtoM3DSolid(M, Shape, w, nnodes, ls);
        }
      }
    }
  }
  else { // lumped mass matrix
    fprintf(stderr," *** In Brick32::massMatrix: Lumped mass matrix NOT implemented. Abort.\n");
    exit(-1);
  }

  return M;
}

FullSquareMatrix
Brick32::stiffness(const CoordSet &cs, double *d, int flg) const
{
  const int nnodes = 32;
  const int ndofs = 96;

  double X[32], Y[32], Z[32];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic material
    // transform local constitutive matrix to global frame
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);

  FullSquareMatrix K(ndofs, d);
  K.zero();

  int ls[96] = {0,3,6, 9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,
                1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67,70,73,76,79,82,85,88,91,94,
                2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92,95};

  // integration: loop over Gauss pts
  const int numgauss= 4;
  double wx,wy,wz,w;
  double m[3], Shape[32], DShape[32][3];
  double dOmega; // det of jacobian
  int jSign = 0;

  for(int i=1; i<=numgauss; i++) {
    _FORTRAN(lgauss)(numgauss,i,&m[0],&wx);
    for(int j=1; j<=numgauss; j++) {
      _FORTRAN(lgauss)(numgauss,j,&m[1],&wy);
      for(int k=1; k<=numgauss; k++) {
        _FORTRAN(lgauss)(numgauss,k,&m[2],&wz);
        dOmega = Hexa32ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
        checkJacobian(&dOmega, &jSign, getGlNum()+1, "Brick32::stiffness");
#endif
        w = fabs(dOmega)*wx*wy*wz;
        addBtCBtoK3DSolid(K, DShape, C, w, nnodes, ls);
      }
    }
  }

  return K;
}


void Brick32::aRubberStiffnessDerivs(CoordSet & cs, complex<double> *d,
                                            int n, double omega) {
  const int nnodes = 32;
  const int ndofs = 96;

  double X[32], Y[32], Z[32];
  cs.getCoordinates(nn, nnodes, X, Y, Z);


  int ls[96] = {0,3,6, 9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,
                1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67,70,73,76,79,82,85,88,91,94,
                2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92,95};

  // integration: loop over Gauss pts
  const int numgauss= 4;
  double wx,wy,wz,w;
  double m[3], Shape[32], DShape[32][3];
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

  for(int i=1; i<=numgauss; i++) {
    _FORTRAN(lgauss)(numgauss,i,&m[0],&wx);
    for(int j=1; j<=numgauss; j++) {
      _FORTRAN(lgauss)(numgauss,j,&m[1],&wy);
      for(int k=1; k<=numgauss; k++) {
        _FORTRAN(lgauss)(numgauss,k,&m[2],&wz);
        dOmega = Hexa32ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
        checkJacobian(&dOmega, &jSign, getGlNum()+1, "Brick32::stiffness");
#endif
        w = fabs(dOmega)*wx*wy*wz;
        addBtCBtoK3DSolid(Kmr, DShape, Cm, w, nnodes, ls);
        addBtCBtoK3DSolid(Klr, DShape, Cl, w, nnodes, ls);
      }
    }
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
Brick32::numNodes() const
{ 
  if(useFull)
    return 32; 
  else
    return 8; // to ignore effect of mid-side nodes in dec
}

int
Brick32::numDofs() const
{
  return 96;
}

int
Brick32::getTopNumber() const
{
  return 191;
}

int*
Brick32::nodes(int *p) const
{
  if(useFull) {
    if(!p) p = new int[32];
    for(int i=0; i<32; i+=2) {
      p[i] = nn[i]; p[i+1] = nn[i+1];
    }
  }
  else {
    if(!p) p = new int[8];
    for(int i=0; i<8; ++i)
      p[i] = nn[i];
  }

  return p;
}

int*
Brick32::dofs(DofSetArray &dsa, int *p) const
{
  if(!p) p = new int[96];

  for(int i=0; i<32; i++)
    dsa.number(nn[i], DofSet::XYZdisp, p+3*i);

  return p;
}

void
Brick32::markDofs(DofSetArray &dsa) const
{
  dsa.mark(nn, numNodes(), DofSet::XYZdisp);
}

void
Brick32::setMaterial(NLMaterial *_mat)
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
Brick32::numStates()
{
#ifdef USE_EIGEN3
  int numGaussPoints = NLHexahedral32::numGaussPoints;
  return (mat) ? numGaussPoints*mat->getNumStates() : 0;
#else
  return 0;
#endif
}

void
Brick32::initStates(double *st)
{
#ifdef USE_EIGEN3
  if(mat) {
    int ninterns = mat->getNumStates();
    int numGaussPoints = NLHexahedral32::numGaussPoints;

    for(int i = 0; i < numGaussPoints; ++i)
      mat->initStates(st+i*ninterns);
  }
#endif
}

Corotator *
Brick32::getCorotator(CoordSet &cs, double *kel, int, int)
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
    MatNLElement *ele = new NLHexahedral32(nn);
    ele->setMaterial(mat);
    ele->setGlNum(glNum);
    ele->setProp(prop);
    return new MatNLCorotator(ele);
#endif
  }
  else {
    return new Brick32Corotator(nn, prop->E, prop->nu, cs, prop->Ta, prop->W, prop->ymtt, prop->ctett);
  }
  printf("WARNING: Corotator not implemented for element %d\n", glNum+1); return 0;
}

int
Brick32::getDecFace(int iFace, int *fn)
{
  switch(iFace) {
    case 0: fn[0] = nn[0]; fn[1] = nn[1]; fn[2] = nn[2]; fn[3] = nn[3]; break;
    case 1: fn[0] = nn[4]; fn[1] = nn[5]; fn[2] = nn[6]; fn[3] = nn[7]; break;
    case 2: fn[0] = nn[3]; fn[1] = nn[0]; fn[2] = nn[4]; fn[3] = nn[7]; break;
    case 3: fn[0] = nn[0]; fn[1] = nn[1]; fn[2] = nn[5]; fn[3] = nn[4]; break;
    case 4: fn[0] = nn[2]; fn[1] = nn[1]; fn[2] = nn[5]; fn[3] = nn[6]; break;
    default:
    case 5: fn[0] = nn[3]; fn[1] = nn[2]; fn[2] = nn[6]; fn[3] = nn[7]; break;
  }
  return 4;
}

int
Brick32::getFace(int iFace, int *fn)
{
  // note: all face normals are outward pointing
  switch(iFace) {
    case 0: fn[0] = nn[3];  fn[1] = nn[2];  fn[2] = nn[1];  fn[3] = nn[0];
            fn[4] = nn[13]; fn[5] = nn[12]; fn[6] = nn[11]; fn[7] = nn[10];
            fn[8] = nn[9];  fn[9] = nn[8];  fn[10]= nn[15]; fn[11]= nn[14]; break;
    case 1: fn[0] = nn[4];  fn[1] = nn[5];  fn[2] = nn[6];  fn[3] = nn[7];
            fn[4] = nn[16]; fn[5] = nn[17]; fn[6] = nn[18]; fn[7] = nn[19];
            fn[8] = nn[20]; fn[9] = nn[21]; fn[10]= nn[22]; fn[11]= nn[23]; break;
    case 2: fn[0] = nn[3];  fn[1] = nn[0];  fn[2] = nn[4];  fn[3] = nn[7];
            fn[4] = nn[14]; fn[5] = nn[15]; fn[6] = nn[24]; fn[7] = nn[25];
            fn[8] = nn[23]; fn[9] = nn[22]; fn[10]= nn[31]; fn[11]= nn[30]; break;
    case 3: fn[0] = nn[0];  fn[1] = nn[1];  fn[2] = nn[5];  fn[3] = nn[4];
            fn[4] = nn[8];  fn[5] = nn[9];  fn[6] = nn[26]; fn[7] = nn[27];
            fn[8] = nn[17]; fn[9] = nn[16]; fn[10]= nn[25]; fn[11]= nn[24]; break;
    case 4: fn[0] = nn[6];  fn[1] = nn[5];  fn[2] = nn[1];  fn[3] = nn[2];
            fn[4] = nn[19]; fn[5] = nn[18]; fn[6] = nn[27]; fn[7] = nn[26];
            fn[8] = nn[10]; fn[9] = nn[11]; fn[10]= nn[28]; fn[11]= nn[29]; break;
    default:
    case 5: fn[0] = nn[7];  fn[1] = nn[6];  fn[2] = nn[2];  fn[3] = nn[3];
            fn[4] = nn[21]; fn[5] = nn[20]; fn[6] = nn[29]; fn[7] = nn[28];
            fn[8] = nn[12]; fn[9] = nn[13]; fn[10]= nn[30]; fn[11]= nn[31]; break;
  }
  return 12;
}

void
Brick32::getCFrame(CoordSet &cs, double cFrame[3][3]) const
{
  if(Brick32::cFrame) {
    cFrame[0][0] = Brick32::cFrame[0]; cFrame[0][1] = Brick32::cFrame[1]; cFrame[0][2] = Brick32::cFrame[2];
    cFrame[1][0] = Brick32::cFrame[3]; cFrame[1][1] = Brick32::cFrame[4]; cFrame[1][2] = Brick32::cFrame[5];
    cFrame[2][0] = Brick32::cFrame[6]; cFrame[2][1] = Brick32::cFrame[7]; cFrame[2][2] = Brick32::cFrame[8];
  }
  else {
    cFrame[0][0] = cFrame[1][1] = cFrame[2][2] = 1.;
    cFrame[0][1] = cFrame[0][2] = cFrame[1][0] = cFrame[1][2] = cFrame[2][0] = cFrame[2][1] = 0.;
  }
}

