// ---------------------------------------------------------------------
// HB - 05-22-05
// ---------------------------------------------------------------------
// 15 nodes wedge element 
// Serendipity finite element basis
// Iso-parametric formulation
// ---------------------------------------------------------------------
#include <cstdio>
#include <iostream>
#include <Element.d/Penta15.d/Penta15.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Corotational.d/Penta15Corotator.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/NLPentahedral.h>
#include <Element.d/Utils.d/SolidElemUtils.h>
#include <Corotational.d/MatNLCorotator.h>

extern "C" {
void _FORTRAN(brkcmt)(double&, double&, double*);
}

double Penta15ShapeFct(double Shape[15], double DShape[15][3], double m[3], double X[15], double Y[15], double Z[15]);
void Penta15ShapeFct(double Shape[15], double m[3]);

extern bool useFull;

double weight3d8[9] = { 0.092592592592593,0.092592592592593,0.092592592592593,
                        0.148148148148148,0.148148148148148,0.148148148148148,
                        0.092592592592593,0.092592592592593,0.092592592592593};
double gauss3d8[9][3] = { {0.166666666666667,0.166666666666667,-0.774596669241483},
                          {0.666666666666667,0.166666666666667,-0.774596669241483},
                          {0.166666666666667,0.666666666666667,-0.774596669241483},
                          {0.166666666666667,0.166666666666667,0.},
                          {0.666666666666667,0.166666666666667,0.},
                          {0.166666666666667,0.666666666666667,0.},
                          {0.166666666666667,0.166666666666667,0.774596669241483},
                          {0.666666666666667,0.166666666666667,0.774596669241483},
                          {0.166666666666667,0.666666666666667,0.774596669241483} };

Penta15::Penta15(int* nodenums)
{
  for(int i=0; i<15; i++)
    nn[i] = nodenums[i];

  cFrame = 0;
  cCoefs = 0;
  mat = 0;
}

Penta15::~Penta15()
{
  if(cCoefs && mat) delete mat;
}

Element *
Penta15::clone()
{
  return new Penta15(*this);
}

void
Penta15::renum(const int *table)
{
  for(int i=0; i<15; i++)
    nn[i] = table[nn[i]];
}

void
Penta15::renum(EleRenumMap& table)
{
  for(int i=0; i<15; i++)
    nn[i] = table[nn[i]];
}

void
Penta15::getVonMises(Vector& stress, Vector& weight, CoordSet &cs,
                     Vector& elDisp, int strInd, int surface, double *ndTemps,
                     double ylayer, double zlayer, int avgnum)
{
  const int nnodes = 15;
  weight = 1.0;

  double X[15], Y[15], Z[15];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  // Flags to calculate Von Mises stress and/or Von Mises strain
  bool vmflg     = (strInd==6) ? true : false;
  bool strainFlg = (strInd==13)? true : false;
  bool meanVms   = false; // This can be used to force averaging of the Von Mises stress & strain.

  double elStress[15][7];
  double elStrain[15][7];

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
  double nodeRefCoord[15][3] = {{  0., 0., -1.},{  1.,  0., -1.},{ 0.,  1., -1.},
                                {  0., 0.,  1.},{  1.,  0.,  1.},{ 0.,  1.,  1.},
                                { 0.5, 0., -1.},{ 0.5, 0.5, -1.},{ 0., 0.5, -1.},
                                { 0.5, 0.,  1.},{ 0.5, 0.5,  1.},{ 0., 0.5,  1.},
                                {  0., 0.,  0.},{  1.,  0.,  0.},{ 0.,  1.,  0.}};

  double Shape[15], DShape[15][3];
  for(int inode=0; inode<nnodes; inode++) {
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Penta15ShapeFct(Shape, DShape, m, X, Y, Z); 
    computeStressAndEngStrain3DSolid(elStress[inode],elStrain[inode], C, DShape, elDisp.data(), nnodes);

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
Penta15::getAllStress(FullM& stress, Vector& weight, CoordSet &cs,
                      Vector& elDisp, int strInd,int surface, double *ndTemps)
{
  const int nnodes = 15;
  weight = 1.0;
  
  double X[15], Y[15], Z[15];
  cs.getCoordinates(nn, nnodes, X, Y, Z);
      
  double elStress[15][7];
  double elStrain[15][7];

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
  double nodeRefCoord[15][3] = {{  0., 0., -1.},{  1.,  0., -1.},{ 0.,  1., -1.},
                                {  0., 0.,  1.},{  1.,  0.,  1.},{ 0.,  1.,  1.},
                                { 0.5, 0., -1.},{ 0.5, 0.5, -1.},{ 0., 0.5, -1.},
                                { 0.5, 0.,  1.},{ 0.5, 0.5,  1.},{ 0., 0.5,  1.},
                                {  0., 0.,  0.},{  1.,  0.,  0.},{ 0.,  1.,  0.}};

  double Shape[15], DShape[15][3];
  for(int inode=0; inode<nnodes; inode++) {
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Penta15ShapeFct(Shape, DShape, m, X, Y, Z); 
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
Penta15::getMass(const CoordSet& cs) const
{
  const int nnodes = 15;
  double X[15], Y[15], Z[15];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  // integration: loop over gauss pts
  const int ngauss = 9;
  double shape[nnodes];
  double dShape[nnodes][3];
  double volume = 0.0;

  for(int i = 0; i < ngauss; i++) {
    double dOmega = Penta15ShapeFct(shape, dShape, gauss3d8[i], X, Y, Z);
    volume += fabs(dOmega)*weight3d8[i];
  }

  return volume*prop->rho;
}

void
Penta15::getGravityForce(CoordSet& cs, double *gravityAcceleration, 
                         Vector& gravityForce, int gravflg, GeomState *geomState)
{
  int nnodes = 15;

  // Lumped
  if (gravflg != 2) {

    double totmas = getMass(cs);

    double m[45*45];
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
    double x[15], y[15], z[15];
    cs.getCoordinates(nn, numNodes(), x, y, z);

    double lforce[15];
    for(int i = 0; i < nnodes; ++i) lforce[i] = 0.0;

    // integration: loop over Gauss pts
    const int ngauss = 9;
    double shape[15], dShape[15][3];

    for(int i = 0; i < ngauss; i++) {
      double dOmega = Penta15ShapeFct(shape, dShape, gauss3d8[i], x, y, z);
      for(int j = 0; j < nnodes; ++j) lforce[j] += fabs(dOmega)*weight3d8[i]*shape[j];
    }

    for(int i = 0; i < nnodes; ++i)
      for(int j = 0; j < 3; ++j)
        gravityForce[3*i+j] = lforce[i]*gravityAcceleration[j]*prop->rho;

  }
}

void
Penta15::getThermalForce(CoordSet &cs, Vector &ndTemps,
                         Vector &elementThermalForce, int glflag,
                         GeomState *geomState)
{
  // ASSUME CONSTANT THERMAL EXPANSION COEFF. & REFERENCE TEMPERATURE OVER THE ELEMENT 
  const int nnodes = 15;
  const int ndofs = 45;

  // initialize nodal thermal forces
  for(int i=0; i<ndofs; i++) elementThermalForce[i] = 0.0;

  // for nonlinear analyses, the thermal load for this element is now computed in getStiffAndForce
  if(geomState) return;

  double X[15], Y[15], Z[15];
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
  int ngauss = 9;
  double Shape[15], DShape[15][3];
  double w, J;

  for(int i = 0; i < ngauss; i++) {
    // compute shape fcts & their derivatives at the Gauss pt
    J = Penta15ShapeFct(Shape, DShape, gauss3d8[i], X, Y, Z);
    w = weight3d8[i]*fabs(J);
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
Penta15::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
  const int nnodes = 15;
  const int ndofs = 45;

  FullSquareMatrix M(ndofs, mel);

  if(cmflg) { // consistent mass matrix

    double x[15], y[15], z[15];
    cs.getCoordinates(nn, nnodes, x, y, z);

    M.zero();
    int ls[45] = {0,3,6, 9,12,15,18,21,24,27,30,33,36,39,42,
                  1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,
                  2,5,8,11,14,17,20,23,26,29,32,35,38,41,44};

    // integration: loop over Gauss pts
    const int ngauss = 9;
    double shape[nnodes];
    double dShape[nnodes][3];

    for(int i = 0; i < ngauss; i++) {
      double dOmega = Penta15ShapeFct(shape, dShape, gauss3d8[i], x, y, z);
      double w = fabs(dOmega)*weight3d8[i]*prop->rho;
      addNtDNtoM3DSolid(M, shape, w, nnodes, ls);
    }
  }
  else { // lumped mass matrix
    fprintf(stderr," *** In Penta15::massMatrix: Lumped mass matrix NOT implemented. Abort.\n");
    exit(-1);
  }

  return M;
}

FullSquareMatrix
Penta15::stiffness(const CoordSet &cs, double *d, int flg) const
{
  const int nnodes = 15;
  const int ndofs = 45;

  double X[15], Y[15], Z[15];
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

  int ls[45] = {0,3,6, 9,12,15,18,21,24,27,30,33,36,39,42,
                1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,
                2,5,8,11,14,17,20,23,26,29,32,35,38,41,44};

  // integration: loop over Gauss pts
  const int ngauss = 9;
  double shape[nnodes];
  double dShape[nnodes][3];

  for(int i = 0; i < ngauss; i++) {
    double dOmega = Penta15ShapeFct(shape, dShape, gauss3d8[i], X, Y, Z);
    double w = fabs(dOmega)*weight3d8[i];
    addBtCBtoK3DSolid(K, dShape, C, w, nnodes, ls);
  }

  return K;
}

int
Penta15::numNodes() const
{ 
  if(useFull)
    return 15; 
  else
    return 6;
}

int
Penta15::numDofs() const
{
  return 45;
}

// treat as a 6-node penta
// this is because xpost does not have a 15-node penta
int
Penta15::getTopNumber() const
{
  return 124;
}

int
Penta15::numTopNodes() const
{
  return 6;
}

int*
Penta15::nodes(int *p) const
{
  if(useFull) {
    if(!p) p = new int[15];
    for(int i=0; i<15; i+=3) {
      p[i] = nn[i]; p[i+1] = nn[i+1]; p[i+2] = nn[i+2];
    }
  }
  else {
    if(!p) p = new int[6];
    for(int i=0; i<6; i++)
      p[i] = nn[i];
  }

  return p;
}

int*
Penta15::dofs(DofSetArray &dsa, int *p) const
{
  if(!p) p = new int[45];

  for(int i=0; i<15; i++)
    dsa.number(nn[i], DofSet::XYZdisp, p+3*i);

  return p;
}

void
Penta15::markDofs(DofSetArray &dsa) const
{
  dsa.mark(nn, numNodes(), DofSet::XYZdisp);
}

void
Penta15::setMaterial(NLMaterial *_mat)
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
Penta15::numStates()
{
#ifdef USE_EIGEN3
  int numGaussPoints = NLPentahedral15::numGaussPoints;
  return (mat) ? numGaussPoints*mat->getNumStates() : 0;
#else
  return 0;
#endif
}

void
Penta15::initStates(double *st)
{
#ifdef USE_EIGEN3
  if(mat) {
    int ninterns = mat->getNumStates();
    int numGaussPoints = NLPentahedral15::numGaussPoints;

    for(int i = 0; i < numGaussPoints; ++i)
      mat->initStates(st+i*ninterns);
  }
#endif
}

Corotator *
Penta15::getCorotator(CoordSet &cs, double *kel, int, int)
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
    MatNLElement *ele = new NLPentahedral15(nn);
    ele->setMaterial(mat);
    ele->setGlNum(glNum);
    ele->setProp(prop);
    return new MatNLCorotator(ele);
#endif
  }
  else {
    return new Penta15Corotator(nn, prop->E, prop->nu, cs, prop->Ta, prop->W, prop->ymtt, prop->ctett);
  }
  printf("WARNING: Corotator not implemented for element %d\n", glNum+1); return 0;
}

int
Penta15::getDecFace(int iFace, int *fn)
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

int
Penta15::getFace(int iFace, int *fn)
{
  int count;
  switch(iFace) {
    case 0: fn[0] = nn[0]; fn[1] = nn[2];  fn[2] = nn[1];
            fn[3] = nn[8]; fn[4] = nn[7];  fn[5] = nn[6]; count = 6; break;
    case 1: fn[0] = nn[3]; fn[1] = nn[4];  fn[2] = nn[5];
            fn[3] = nn[9]; fn[4] = nn[10]; fn[5] = nn[11]; count = 6; break;
    case 2: fn[0] = nn[0]; fn[1] = nn[1];  fn[2] = nn[4];  fn[3] = nn[3];
            fn[4] = nn[6]; fn[5] = nn[13]; fn[6] = nn[9];  fn[7] = nn[12]; count = 8; break;
    case 3: fn[0] = nn[1]; fn[1] = nn[2];  fn[2] = nn[5];  fn[3] = nn[4];
            fn[4] = nn[7]; fn[5] = nn[14]; fn[6] = nn[10]; fn[7] = nn[13]; count = 8; break;
    default:
    case 4: fn[0] = nn[2]; fn[1] = nn[0];  fn[2] = nn[3];  fn[3] = nn[5];
            fn[4] = nn[8]; fn[5] = nn[12]; fn[6] = nn[11]; fn[7] = nn[14]; count = 8; break;
  }
  return count;
}

void
Penta15::getCFrame(CoordSet &cs, double cFrame[3][3]) const
{
  if(Penta15::cFrame) {
    cFrame[0][0] = Penta15::cFrame[0]; cFrame[0][1] = Penta15::cFrame[1]; cFrame[0][2] = Penta15::cFrame[2];
    cFrame[1][0] = Penta15::cFrame[3]; cFrame[1][1] = Penta15::cFrame[4]; cFrame[1][2] = Penta15::cFrame[5];
    cFrame[2][0] = Penta15::cFrame[6]; cFrame[2][1] = Penta15::cFrame[7]; cFrame[2][2] = Penta15::cFrame[8];
  }
  else {
    cFrame[0][0] = cFrame[1][1] = cFrame[2][2] = 1.;
    cFrame[0][1] = cFrame[0][2] = cFrame[1][0] = cFrame[1][2] = cFrame[2][0] = cFrame[2][1] = 0.;
  }
}
