#include <cstdio>
#include <iostream>
#include <cmath>

#include <Element.d/Tetra.d/Tetrahedral.h>
#include <Element.d/Tetra.d/TetraElementTemplate.cpp>
#include <Element.d/Tetra.d/TetraElementStressWRTNodalCoordinateSensitivity.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Corotational.d/TetCorotator.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/NLTetrahedral.h>
#include <Element.d/Utils.d/SolidElemUtils.h>
#include <Element.d/Helm.d/ARubberF.h>
#include <Corotational.d/MatNLCorotator.h>

#define CHECK_JACOBIAN // force check nullity & constant sign of jacobian over el.

extern int verboseFlag;

extern "C" {
void _FORTRAN(mass23)(double*, double*, double*, double&, double*,
                      const int&, double*, double*, const int&, 
                      double&, const int&);

void _FORTRAN(sands23)(const int&, double*, double*, double*, double&, 
                       double&, double*, double*, double*, const int&,
                       const int&,
                       const int&, const int&, const int&, const int&, double*);

void _FORTRAN(brkcmt)(double&, double&, double*);
}

void Tetra4ShapeFct(double Shape[4], double DShape[4][3], double m[3], double J, double a[3][3]);
double computeTetra4Jacobian(double X[4], double Y[4], double Z[4], double a[3][3]);
void computeTetra4dadx(double X[4], double Y[4], double Z[4], double dadx[3][3][12]);
void computedShape(double dShape[4][3]);
void computeTetra4dadxTimesdShapeFct(double dShape[4][3], double J, double dadx[3][3][12], double [4][3][12]);

Tetrahedral::Tetrahedral(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
  nn[2] = nodenums[2];
  nn[3] = nodenums[3];

  cFrame = 0;
  cCoefs = 0;
  mat = 0;
}

Tetrahedral::~Tetrahedral()
{
  if(cCoefs && mat) delete mat;
}

Element *
Tetrahedral::clone()
{
  return new Tetrahedral(*this);
}

void
Tetrahedral::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

void
Tetrahedral::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

void
Tetrahedral::getVonMises(Vector& stress, Vector& weight, CoordSet &cs,
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

  double x[4], y[4], z[4];
  cs.getCoordinates(nn, numNodes(), x, y, z);

  int vmflg = 0, strainFlg = 0;
  // Flags sands23 to calculate Von Mises stress
  if(strInd == 6) vmflg = 1;
  // Flags sands23 to calculate Von Mises strain
  if(strInd == 13) strainFlg = 1;

  const int maxgus = 4;
  const int maxstr = 7;
  const int elm    = 1;
  const int msize  = 1;
  const int outerr = 6;

  double elStress[maxgus][maxstr], elStrain[maxgus][maxstr], elStressSen[4][12][6];

  _FORTRAN(sands23)(elm, x, y, z, prop->E, prop->nu, elDisp.data(),
                    (double*)elStress, (double*)elStrain,
                    maxgus, maxstr, msize, outerr, vmflg, strainFlg,(double*)elStressSen);

  if(strInd < 7) {
    double thermalStress[4] = {0.0,0.0,0.0,0.0};
    if(strInd == 0 || strInd == 1 || strInd == 2) {
      double &Tref  = prop->Ta;
      double &alpha = prop->W;
      double coef = (prop->E*alpha)/(1.0 - 2.0*prop->nu);
      thermalStress[0] = coef*(ndTemps[0]-Tref);
      thermalStress[1] = coef*(ndTemps[1]-Tref);
      thermalStress[2] = coef*(ndTemps[2]-Tref);
      thermalStress[3] = coef*(ndTemps[3]-Tref);
    }
    stress[0] = elStress[0][strInd] - thermalStress[0];
    stress[1] = elStress[1][strInd] - thermalStress[1];
    stress[2] = elStress[2][strInd] - thermalStress[2];
    stress[3] = elStress[3][strInd] - thermalStress[3];
  }
  else if(strInd < 14) {
    stress[0] = elStrain[0][strInd-7];
    stress[1] = elStrain[1][strInd-7];
    stress[2] = elStrain[2][strInd-7];
    stress[3] = elStrain[3][strInd-7];
  }
  else {
    stress[0] = 0;
    stress[1] = 0;
    stress[2] = 0;
    stress[3] = 0;
  }
}

void
Tetrahedral::getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *dDispDisp,
                                                CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                                double *ndTemps, int avgnum, double ylayer, double zlayer)
{
   if(strInd != 6) {
     std::cerr << " ... Error: strInd must be 6 in Tetrahedral::getVonMisesDisplacementSensitivity\n";
     exit(-1);
   }
   if(dStdDisp.numRow() != 12 || dStdDisp.numCol() != 4) {
     std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";
     exit(-1);
   }
   weight = 1.0;

    double x[4], y[4], z[4];
    cs.getCoordinates(nn, numNodes(), x, y, z);

    int vmflg = 0, strainFlg = 0;
    // Flags sands23 to calculate Von Mises stress
    if(strInd == 6) vmflg = 1;
    // Flags sands23 to calculate Von Mises strain
    if(strInd == 13) strainFlg = 1;

    const int maxgus = 4;
    const int maxstr = 7;
    const int elm    = 1;
    const int msize  = 1;
    const int outerr = 6;

    double elStress[maxgus][maxstr], elStrain[maxgus][maxstr], elStressSen[4][12][6];

    _FORTRAN(sands23)(elm, x, y, z, prop->E, prop->nu, elDisp.data(),
                      (double*)elStress, (double*)elStrain,
                      maxgus, maxstr, msize, outerr, vmflg, strainFlg, (double*)elStressSen);

    Vector stress(4);
    if(strInd < 7) {
      stress[0] = elStress[0][strInd];
      stress[1] = elStress[1][strInd];
      stress[2] = elStress[2][strInd];
      stress[3] = elStress[3][strInd];
    }

    double dvmsdStress[4][6];
    for(int i=0; i<4; ++i) {
      dvmsdStress[i][0] = (2.*elStress[i][0]-elStress[i][1]-elStress[i][2])/(2.*elStress[i][6]);    
      dvmsdStress[i][1] = (2.*elStress[i][1]-elStress[i][0]-elStress[i][2])/(2.*elStress[i][6]);    
      dvmsdStress[i][2] = (2.*elStress[i][2]-elStress[i][1]-elStress[i][0])/(2.*elStress[i][6]);    
      dvmsdStress[i][3] = (3.*elStress[i][3])/elStress[i][6];    
      dvmsdStress[i][4] = (3.*elStress[i][4])/elStress[i][6];    
      dvmsdStress[i][5] = (3.*elStress[i][5])/elStress[i][6];
    }

    for(int i=0; i<12; ++i) 
      for(int j=0; j<4; ++j) 
        for(int k=0; k<6; ++k)
          dStdDisp[i][j] += dvmsdStress[j][k]*elStressSen[j][i][k];

    if(dDispDisp) dStdDisp ^= (*dDispDisp);
}

void
Tetrahedral::getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                                   double* ndTemps, int avgnum, double ylayer, double zlayer)
{
#ifdef USE_EIGEN3
  if(strInd != 6) {
    std::cerr << " ... Error: strInd must be 6 in Tetrahedral::getVonMisesNodalCoordinateSensitivity\n";
    exit(-1);
  }
  if(dStdx.numRow() != 12 || dStdx.numCol() != 4) {
    std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";
    exit(-1);
  }
  if(ndTemps != 0) {
    std::cerr << " ... Error: thermal stress should not be passed in sensitivity computation\n";
    exit(-1);
  }
  if(cCoefs) {
    std::cerr << " ... Error: Anisotropic material not supported in sensitivity computation\n";
    exit(-1);                 
  }
  weight = 1.0;

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);
  auto &nd4 = cs.getNode(nn[3]);

  Eigen::Array<double,14,1> dconst;
  dconst.segment<12>(0) = Eigen::Map<Eigen::Matrix<double,12,1> >(elDisp.data()).segment(0,12); // displacements
  dconst[12] = prop->E;
  dconst[13] = prop->nu;

  // integer parameters
  Eigen::Array<int,1,1> iconst;
  iconst[0] = surface;
 
  // inputs
  Eigen::Matrix<double,12,1> q;
  q << nd1.x, nd1.y, nd1.z, nd2.x, nd2.y, nd2.z, nd3.x, nd3.y, nd3.z, nd4.x, nd4.y, nd4.z;
  Eigen::Matrix<double,4,12> dStressdx;
 
  Simo::Jacobian<double,TetraElementStressWRTNodalCoordinateSensitivity> dSdx(dconst,iconst);
  dStressdx = dSdx(q, 0);
  dStdx.copy(dStressdx.data());

#else
  std::cerr << " ... Error! Tetrahedral::getVonMisesNodalCoordinateSensitivity needs Eigen library.\n";
  exit(-1);
#endif
}

void
Tetrahedral::getAllStress(FullM& stress, Vector& weight, CoordSet &cs,
                          Vector& elDisp, int strInd, int surface, double *ndTemps)
{
  if(cCoefs) {
    getAllStressAniso(stress, weight, cs,
                      elDisp, strInd, surface, ndTemps);
    return;
  }

  weight = 1.0;

  double x[4], y[4], z[4];
  cs.getCoordinates(nn, numNodes(), x, y, z);

  int vmflg, strainFlg;
  vmflg = 0;
  strainFlg = 0;

  const int maxgus = 4;
  const int maxstr = 7;
  const int elm    = 1;
  const int msize  = 1;
  const int outerr = 6;

  double elStress[maxgus][maxstr], elStrain[maxgus][maxstr], elStressSen[4][12][6];

  _FORTRAN(sands23)(elm, x, y, z, prop->E, prop->nu, elDisp.data(),
                    (double*)elStress, (double*)elStrain,
                    maxgus, maxstr, msize, outerr, vmflg, strainFlg, (double*)elStressSen);

  // Store all Stress or all Strain as defined by strInd
  int i,j;
  if(strInd == 0) {
    double thermalStress[4] = {0.0,0.0,0.0,0.0};
    if(ndTemps) {
      double &Tref  = prop->Ta;
      double &alpha = prop->W;
      double coef = (prop->E*alpha)/(1.0 - 2.0*prop->nu);
      thermalStress[0] = coef*(ndTemps[0]-Tref);
      thermalStress[1] = coef*(ndTemps[1]-Tref);
      thermalStress[2] = coef*(ndTemps[2]-Tref);
      thermalStress[3] = coef*(ndTemps[3]-Tref);
    }
    for (i=0; i<4; ++i) {
      for (j=0; j<3; ++j) {
        stress[i][j] = elStress[i][j] - thermalStress[i];
      }
      for (j=3; j<6; ++j) {
        stress[i][j] = elStress[i][j];
      }
    }
  }
  else {
    for (i=0; i<4; ++i) {
      for (j=0; j<6; ++j) {
        stress[i][j] = elStrain[i][j];
      }
    }
  }

  // Get Element Principals for each node without averaging
  double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double pvec[3] = {0.0,0.0,0.0};

  for (i=0; i<4; ++i) {
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

void 
Tetrahedral::computeDjDx(double x[4], double y[4], double z[4], double J, double djdx[12])
{
  for(int i=0; i<12; ++i) djdx[i] = 0.0;
  double xd1,xd2,xd3, yd1,yd2,yd3, zd1,zd2,zd3;
  xd1 = (x[1]-x[0]); xd2 = (x[2]-x[0]);  xd3 = (x[3]-x[0]);
  yd1 = (y[1]-y[0]); yd2 = (y[2]-y[0]);  yd3 = (y[3]-y[0]);
  zd1 = (z[1]-z[0]); zd2 = (z[2]-z[0]);  zd3 = (z[3]-z[0]);
  djdx[0] = - yd2*zd3 - yd3*zd1 - yd1*zd2 + yd2*zd1 + yd1*zd3 + yd3*zd2;
  djdx[1] = - xd1*zd3 - xd2*zd1 - xd3*zd2 + xd3*zd1 + xd2*zd3 + xd1*zd2;
  djdx[2] = - xd1*yd2 - xd2*yd3 - xd3*yd1 + xd3*yd2 + xd2*yd1 + xd1*yd3;
  djdx[3] =   yd2*zd3 - yd3*zd2;
  djdx[4] =   xd3*zd2 - xd2*zd3;
  djdx[5] =   xd2*yd3 - xd3*yd2;
  djdx[6] =   yd3*zd1 - yd1*zd3;
  djdx[7] =   xd1*zd3 - xd3*zd1;
  djdx[8] =   xd3*yd1 - xd1*yd3;
  djdx[9] =   yd1*zd2 - yd2*zd1;
  djdx[10] =  xd2*zd1 - xd1*zd2;
  djdx[11] =  xd1*yd2 - xd2*yd1;
  if(J < 0) 
    for(int i=0; i<12; ++i) djdx[i] *= -1;
}

void
Tetrahedral::getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet& cs, double *gravityAcceleration)
{
  int nnodes = 4;
  double x[4], y[4], z[4];
  cs.getCoordinates(nn, numNodes(), x, y, z);

  for(int k = 0; k < 12; ++k) dwdx[k] = 0.0;

  int ngp = 4;
  double w1 = 4.166666666666666e-02;
  double r1 = 1.0;
  double s1 = 0.0;
  double t1 = 0.0;
  double u1 = 0.0;
  double TetGPt2[4][5] = {{r1, s1, t1, u1, w1},
                          {s1, r1, u1, u1, w1},
                          {t1, u1, r1, s1, w1},
                          {u1, u1, s1, r1, w1}};
  double djdx[12];
  double m[3], Shape[8], DShape[8][3], w, a[3][3];
  double dOmega = computeTetra4Jacobian(x,y,z,a); //det of jacobian
  computeDjDx(x, y, z, dOmega, djdx);
//  for(int i=0; i<12; ++i) fprintf(stderr,"::weightDerivativeWRTNodalCoordinate   djdx[%d] = %e\n", i, djdx[i]);
//  for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) fprintf(stderr,"::weightDerivativeWRTNodalCoordinate   a[%d][%d] = %e\n", i, j, a[i][j]);

  for(int igp = 0; igp < ngp; igp++) {
    // get x, y, z  position & weight of the integration pt
    m[0] = TetGPt2[igp][0]; m[1] = TetGPt2[igp][1]; m[2] = TetGPt2[igp][2]; w = TetGPt2[igp][4];
    Tetra4ShapeFct(Shape, DShape, m, dOmega, a);
    w *= prop->rho;

    for(int k = 0; k < 12; ++k) dwdx[k] += w*djdx[k];
  }

  double gravAccNorm = sqrt(gravityAcceleration[0]*gravityAcceleration[0] +
                            gravityAcceleration[1]*gravityAcceleration[1] +
                            gravityAcceleration[2]*gravityAcceleration[2]);
  dwdx *= gravAccNorm;

}

double
Tetrahedral::getMass(const CoordSet& cs) const
{
  double x[4], y[4], z[4]; 
  cs.getCoordinates(nn, numNodes(), x, y, z);

  double ElementMassMatrix[12][12];

  double *gravityAcceleration = 0, *grvfor = 0, totmas = 0.0;

  int grvflg = 0, masflg = 1;

  const int numdof = 12;
  _FORTRAN(mass23)(x, y, z, prop->rho, (double*)ElementMassMatrix, numdof,
                   gravityAcceleration, grvfor, grvflg, totmas, masflg);

  return totmas;
}

void
Tetrahedral::getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration,
                                                       GenFullM<double> &dGfdx, int gravflg, GeomState *geomState)
{
  int nnodes = 4;

  double x[4], y[4], z[4];
  cs.getCoordinates(nn, numNodes(), x, y, z);

  // Lumped
  if (gravflg != 2) {
    std::cerr << "Tetrahedral::getGravityForceNodalCoordinateSensitivity is not implemented for lumped\n";
  }
  // Consistent
  else {
   
    double lforce[12][nnodes];
    for(int k = 0; k < 12; ++k) for(int i = 0; i < nnodes; ++i) lforce[k][i] = 0.0;

    // integration: loop over Gauss pts
    // hard coded order 2 tetrahedral quadrature rule: {r,s,t,u(=1-r-s-t),w}
    int ngp = 4;
    double w1 = 1./24.;
    double r1 = (5.+3.*sqrt(5.))/20.;
    double s1 = (5.-  sqrt(5.))/20.;
    double t1 = s1;
    double u1 = 1.-r1-s1-t1;
    double TetGPt2[4][5] = {{r1, s1, t1, u1, w1},
                            {s1, t1, u1, r1, w1},
                            {t1, u1, r1, s1, w1},
                            {u1, r1, s1, t1, w1}};
    double w;
    double m[3], Shape[8], DShape[8][3], a[3][3];
    double dOmega = computeTetra4Jacobian(x,y,z,a); //det of jacobian

    double djdx[12];
    computeDjDx(x,y,z,dOmega, djdx);

    for(int igp = 0; igp < ngp; igp++) {
      // get x, y, z  position & weight of the integration pt
      m[0] = TetGPt2[igp][0]; m[1] = TetGPt2[igp][1]; m[2] = TetGPt2[igp][2]; w = TetGPt2[igp][4];
      Tetra4ShapeFct(Shape, DShape, m, dOmega, a);
      w *= prop->rho;

      for(int k = 0; k < 12; ++k)
        for(int n = 0; n < nnodes; ++n)
          lforce[k][n] += w*djdx[k]*Shape[n];
    }

    for(int i = 0; i < nnodes; ++i) {
      for(int k = 0; k < 12; ++k) {
        dGfdx[k][3*i+0] = lforce[k][i]*gravityAcceleration[0];
        dGfdx[k][3*i+1] = lforce[k][i]*gravityAcceleration[1];
        dGfdx[k][3*i+2] = lforce[k][i]*gravityAcceleration[2];
      }
    }
  }
}

void
Tetrahedral::getGravityForce(CoordSet& cs, double *gravityAcceleration,
                             Vector& gravityForce, int gravflg, GeomState *geomState)
{
  int nnodes = 4;

  double x[4], y[4], z[4];
  cs.getCoordinates(nn, numNodes(), x, y, z);

  // Lumped
  if (gravflg != 2) {

    double ElementMassMatrix[12][12];
    double grvfor[3], totmas = 0.0;
    int grvflg = 1, masflg = 0;
    const int numdof = 12;
    _FORTRAN(mass23)(x, y, z, prop->rho, (double*)ElementMassMatrix, numdof,
                     gravityAcceleration, grvfor, grvflg, totmas, masflg);

    grvfor[0] /= double(nnodes);
    grvfor[1] /= double(nnodes);
    grvfor[2] /= double(nnodes);

    for(int i = 0; i < nnodes; ++i) {
      gravityForce[3*i+0] = grvfor[0];
      gravityForce[3*i+1] = grvfor[1];
      gravityForce[3*i+2] = grvfor[2];
    }
  }
  // Consistent
  else {

    double lforce[4];
    for(int i = 0; i < nnodes; ++i) lforce[i] = 0.0;

    // integration: loop over Gauss pts
    // hard coded order 2 tetrahedral quadrature rule: {r,s,t,u(=1-r-s-t),w}
    int ngp = 4;
    double w1 = 1./24.;
    double r1 = (5.+3.*sqrt(5.))/20.;
    double s1 = (5.-  sqrt(5.))/20.;
    double t1 = s1;
    double u1 = 1.-r1-s1-t1;
    double TetGPt2[4][5] = {{r1, s1, t1, u1, w1},
                            {s1, t1, u1, r1, w1},
                            {t1, u1, r1, s1, w1},
                            {u1, r1, s1, t1, w1}};
    double w;
    double m[3], Shape[8], DShape[8][3], a[3][3];
    double dOmega = computeTetra4Jacobian(x,y,z,a); //det of jacobian

    for(int igp = 0; igp < ngp; igp++) {
      // get x, y, z  position & weight of the integration pt
      m[0] = TetGPt2[igp][0]; m[1] = TetGPt2[igp][1]; m[2] = TetGPt2[igp][2]; w = TetGPt2[igp][4];
      Tetra4ShapeFct(Shape, DShape, m, dOmega, a);
      w *= prop->rho*fabs(dOmega);

      for(int n = 0; n < nnodes; ++n)
        lforce[n] += w*Shape[n];
    }

    for(int i = 0; i < nnodes; ++i) {
      gravityForce[3*i+0] = lforce[i]*gravityAcceleration[0];
      gravityForce[3*i+1] = lforce[i]*gravityAcceleration[1];
      gravityForce[3*i+2] = lforce[i]*gravityAcceleration[2];
    }
  }
}

void
Tetrahedral::getThermalForce(CoordSet &cs, Vector &ndTemps,
                             Vector &elementThermalForce, int glflag,
                             GeomState *geomState)
{
  // ASSUME CONSTANT THERMAL EXPANSION COEFF. & REFERENCE TEMPERATURE OVER THE ELEMENT 
  const int nnodes = 4;
  const int ndofs = 12;

  // initialize nodal thermal forces
  for(int i=0; i<ndofs; i++) elementThermalForce[i] = 0.0;

  // for nonlinear analyses, the thermal load for this element is now computed in getStiffAndForce
  if(geomState) return;

  double X[4], Y[4], Z[4];
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
 
  // NUMERICAL INTEGRATION BY GAUSS PTS
  // integration: loop over Gauss pts
  // hard coded order 1 tetrahedral quadrature rule: {r,s,t,u(=1-r-s-t),w}
  int ngp = 1;
  double TetGPt1[1][5] = {{1./4.,1./4.,1./4.,1./4.,1./6.}};
  double m[3], Shape[4], DShape[4][3];
  double w, J, a[3][3];
  int jSign = 0;

  J = computeTetra4Jacobian(X,Y,Z, a);
  for(int igp=0; igp<ngp; igp++) {
    // get x, y, z position & weight of the integration pt
    m[0] = TetGPt1[igp][0]; m[1] = TetGPt1[igp][1]; m[2] = TetGPt1[igp][2]; w = TetGPt1[igp][4];
    Tetra4ShapeFct(Shape, DShape, m, J, a);
#ifdef CHECK_JACOBIAN
    checkJacobian(&J, &jSign, getGlNum()+1, "Tetrahedral::getThermalForce");
#endif
    w *= fabs(J);
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
Tetrahedral::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
  const int nnodes = 4;
  const int ndofs = 12;

  double X[4], Y[4], Z[4];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  FullSquareMatrix M(ndofs,mel);

  if(cmflg) { // consistent mass matrix
    M.zero();
    int ls[12] = {0,3,6,9,1,4,7,10,2,5,8,11};

    // integration: loop over Gauss pts
    // hard coded order 2 tetrahedral quadrature rule: {r,s,t,u(=1-r-s-t),w}
    int ngp = 4;
    double w1 = 1./24.;
    double r1 = (5.+3.*sqrt(5.))/20.;
    double s1 = (5.-  sqrt(5.))/20.;
    double t1 = s1;
    double u1 = 1.-r1-s1-t1;
    double TetGPt2[4][5] = {{r1, s1, t1, u1, w1},
                            {s1, t1, u1, r1, w1},
                            {t1, u1, r1, s1, w1},
                            {u1, r1, s1, t1, w1}};
    double w;
    double m[3], Shape[8], DShape[8][3], a[3][3];
    double dOmega = computeTetra4Jacobian(X,Y,Z, a); // det of jacobian
    int jSign = 0;

    for(int igp = 0; igp < ngp; igp++) {
      // get x, y, z  position & weight of the integration pt
      m[0] = TetGPt2[igp][0]; m[1] = TetGPt2[igp][1]; m[2] = TetGPt2[igp][2]; w = TetGPt2[igp][4];
      Tetra4ShapeFct(Shape, DShape, m, dOmega, a);
#ifdef CHECK_JACOBIAN
      checkJacobian(&dOmega, &jSign, getGlNum()+1, "Tetrahedral::massMatrix");
#endif
      w *= prop->rho*fabs(dOmega);
      addNtDNtoM3DSolid(M, Shape, w, nnodes, ls);
    }
  }
  else { // Lumped mass matrix
    double *gravityAcceleration = 0, *grvfor = 0, totmas = 0.0;
    int grvflg = 0, masflg = 0;
    _FORTRAN(mass23)(X, Y, Z, prop->rho, (double*)mel, ndofs,
                     gravityAcceleration, grvfor, grvflg, totmas, masflg);
  }

  return M;
}

void
Tetrahedral::getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs)
{
  for(int i=0; i<12; ++i) { 
    if(dStiffdx[i].dim() != 12) { 
      std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";   exit(-1); }
    dStiffdx[i].zero();
  }

  const int nnodes = 4;
  const int ndofs = 12;

  double X[4], Y[4], Z[4], djdx[12], dadx[3][3][12], DDShape[4][3][12], dShape[4][3];
  cs.getCoordinates(nn, nnodes, X, Y, Z);
  double m[3], Shape[4], DShape[4][3], a[3][3];
  double w, dOmega;

  dOmega = computeTetra4Jacobian(X, Y, Z, a);
  computeTetra4dadx(X,Y,Z,dadx);

  computeDjDx(X,Y,Z,dOmega,djdx);
  computedShape(dShape);

  int ls[12] = {0,3,6,9,1,4,7,10,2,5,8,11};

  // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic material
    // transform local constitutive matrix to global frame
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);

  // integration: loop over Gauss pts
  // hard coded order 1 tetrahedral quadrature rule: {r,s,t,u(=1-r-s-t),w}

  int ngp = 1;
  double TetGPt1[1][5] = {{1./4.,1./4.,1./4.,1./4.,1./6.}};
  for(int igp=0; igp<ngp; igp++) {
    // get x, y, z  position & weight of the integration pt
    m[0] = TetGPt1[igp][0]; m[1] = TetGPt1[igp][1]; m[2] = TetGPt1[igp][2]; w = TetGPt1[igp][4];
    Tetra4ShapeFct(Shape, DShape, m, dOmega, a);
    computeTetra4dadxTimesdShapeFct(dShape, dOmega, dadx, DDShape);
    for(int k=0; k<12; ++k) {
      addBtCBtoK3DSolid(dStiffdx[k], DShape, C, -w*djdx[k], nnodes, ls);
    }
    addDBtCDBtodKdx3DSolid(dStiffdx, DShape, DDShape, C, w*fabs(dOmega), nnodes, ls); 
  }
}

//HB (04/15/05)  new implementation of the Tetra4 stiffness matrix to deal
//               with anisotropic constitutive matrix
FullSquareMatrix
Tetrahedral::stiffness(const CoordSet &cs, double *d, int flg) const
{
  const int nnodes = 4;
  const int ndofs = 12;

  double X[4], Y[4], Z[4];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  int ls[12] = {0,3,6,9,1,4,7,10,2,5,8,11};
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
  // hard coded order 1 tetrahedral quadrature rule: {r,s,t,u(=1-r-s-t),w}
  int ngp = 1;
  double TetGPt1[1][5] = {{1./4.,1./4.,1./4.,1./4.,1./6.}};
  double m[3], Shape[4], DShape[4][3], a[3][3];
  double w, dOmega;
  int jSign = 0;

  dOmega = computeTetra4Jacobian(X, Y, Z, a);
  for(int igp=0; igp<ngp; igp++) {
    // get x, y, z  position & weight of the integration pt
    m[0] = TetGPt1[igp][0]; m[1] = TetGPt1[igp][1]; m[2] = TetGPt1[igp][2]; w = TetGPt1[igp][4];
    Tetra4ShapeFct(Shape, DShape, m, dOmega, a);
#ifdef CHECK_JACOBIAN
    checkJacobian(&dOmega, &jSign, getGlNum()+1, "Tetrahedral::stiffness");
#endif
    w *= fabs(dOmega);
    addBtCBtoK3DSolid(K, DShape, C, w, nnodes, ls);
  }

  return K;
}


void Tetrahedral::aRubberStiffnessDerivs(CoordSet & cs, complex<double> *d,
                                            int n, double omega) {
  const int nnodes = 4;
  const int ndofs = 12;

  double X[4], Y[4], Z[4];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  int ls[12] = {0,3,6,9,1,4,7,10,2,5,8,11};

  double *Kl = new double[ndofs*ndofs];
  FullSquareMatrix mKl(ndofs,Kl);
  mKl.zero();

  double Cl[6][6];
  for(int i=0;i<6;i++) for(int j=0;j<6;j++) Cl[i][j] = 0.0;
  Cl[0][0] = Cl[1][1] = Cl[2][2] = 1.0;
  Cl[0][1] = Cl[1][0] = Cl[0][2] = Cl[2][0] = Cl[1][2] = Cl[2][1] = 1.0;
  Cl[3][3] = Cl[4][4] = Cl[5][5] = 0.0;

  double *Km = new double[ndofs*ndofs];
  FullSquareMatrix mKm(ndofs,Km);
  mKm.zero();
  double Cm[6][6];
  for(int i=0;i<6;i++) for(int j=0;j<6;j++) Cm[i][j] = 0.0;
  Cm[0][0] = Cm[1][1] = Cm[2][2] = 2.0;
  Cm[3][3] = Cm[4][4] = Cm[5][5] = 1.0;

  // integration: loop over Gauss pts
  // hard coded order 1 tetrahedral quadrature rule: {r,s,t,u(=1-r-s-t),w}
  int ngp = 1;
  double TetGPt1[1][5] = {{1./4.,1./4.,1./4.,1./4.,1./6.}};
  double m[3], Shape[4], DShape[4][3], a[3][3];
  double w, dOmega;
  int jSign = 0;

  dOmega = computeTetra4Jacobian(X, Y, Z, a);
  for(int igp=0; igp<ngp; igp++) {
    // get x, y, z  position & weight of the integration pt
    m[0] = TetGPt1[igp][0]; m[1] = TetGPt1[igp][1]; m[2] = TetGPt1[igp][2]; w = TetGPt1[igp][4];
    Tetra4ShapeFct(Shape, DShape, m, dOmega, a);
#ifdef CHECK_JACOBIAN
    checkJacobian(&dOmega, &jSign, getGlNum()+1, "Tetrahedral::stiffness");
#endif
    w *= fabs(dOmega);
    addBtCBtoK3DSolid(mKl, DShape, Cl, w, nnodes, ls);
    addBtCBtoK3DSolid(mKm, DShape, Cm, w, nnodes, ls);
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
Tetrahedral::numNodes() const
{
  return 4;
}

int
Tetrahedral::numDofs() const
{
  return 12;
}

int
Tetrahedral::getTopNumber() const
{
  return 123;
}

int*
Tetrahedral::nodes(int *p) const
{
  if(!p) p = new int[numNodes()];
  p[0] = nn[0];
  p[1] = nn[1];
  p[2] = nn[2];
  p[3] = nn[3];
  return p;
}

int*
Tetrahedral::dofs(DofSetArray &dsa, int *p) const
{
  if(!p) p = new int[numDofs()];

  dsa.number(nn[0], DofSet::XYZdisp, p  );
  dsa.number(nn[1], DofSet::XYZdisp, p+3);
  dsa.number(nn[2], DofSet::XYZdisp, p+6);
  dsa.number(nn[3], DofSet::XYZdisp, p+9);

  return p;
}

void
Tetrahedral::markDofs(DofSetArray &dsa) const
{
  dsa.mark(nn, numNodes(), DofSet::XYZdisp);
}

// Stress evaluation in case of anisotropic elastic constitutive matrix
void
Tetrahedral::getVonMisesAniso(Vector &stress, Vector &weight, CoordSet &cs,
                              Vector &elDisp, int strInd, int surface, double *ndTemps,
                              double ylayer, double zlayer, int avgnum)
{
  const int nnodes = 4;
  weight = 1.0;
  
  double X[4], Y[4], Z[4];
  cs.getCoordinates(nn, nnodes, X, Y, Z);
     
  // Flags sands17 to calculate Von Mises stress and/or Von Mises strain
  bool vmflg     = (strInd==6) ? true : false;
  bool strainFlg = (strInd==13)? true : false;
  bool meanVms   = false; // This can be used to force averaging of the Von Mises stress & strain.

  double elStress[4][7];
  double elStrain[4][7];
 
  // get constitutive matrix and coefficients of thermal expansion
  double C[6][6], alpha[6];
  // transform local constitutive matrix to global frame
  rotateConstitutiveMatrix(cCoefs, cFrame, C);
  // transform local coefficients of thermal expansion to global frame
  if(ndTemps) rotateVector(cCoefs+36, cFrame, alpha);
 
  // Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[4][3] = {{0.0,0.0,0.0},{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

  double Shape[4], DShape[4][3], a[3][3];
  double dOmega = computeTetra4Jacobian(X,Y,Z,a);
  for(int inode=0; inode<nnodes; inode++) {
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Tetra4ShapeFct(Shape, DShape, m, dOmega, a); 
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
Tetrahedral::getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
                               Vector &elDisp, int strInd, int surface, double *ndTemps)
{
  const int nnodes = 4;
  weight = 1.0;
  
  double X[4], Y[4], Z[4];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  double elStress[4][6];
  double elStrain[4][6];
 
  // get constitutive matrix and coefficients of thermal expansion
  double C[6][6], alpha[6];
  // transform local constitutive matrix to global frame
  rotateConstitutiveMatrix(cCoefs, cFrame, C);
  // transform local coefficients of thermal expansion to global frame
  if(ndTemps) rotateVector(cCoefs+36, cFrame, alpha);
 
  // Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[4][3] = {{0.0,0.0,0.0},{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

  double Shape[4], DShape[4][3], a[3][3];
  double dOmega = computeTetra4Jacobian(X,Y,Z,a);
  for(int inode=0; inode<nnodes; inode++) {
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Tetra4ShapeFct(Shape, DShape, m, dOmega, a); 
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
Tetrahedral::setMaterial(NLMaterial *_mat)
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
Tetrahedral::numStates()
{
#ifdef USE_EIGEN3
  int numGaussPoints = NLTetrahedral4::numGaussPoints;
  return (mat) ? numGaussPoints*mat->getNumStates(): 0;
#else
  return 0;
#endif
}

void
Tetrahedral::initStates(double *st)
{
#ifdef USE_EIGEN3
  if(mat) {
    int ninterns = mat->getNumStates();
    int numGaussPoints = NLTetrahedral4::numGaussPoints;

    for(int i = 0; i < numGaussPoints; ++i)
      mat->initStates(st+i*ninterns);
  }
#endif
}

Corotator *
Tetrahedral::getCorotator(CoordSet &cs, double *kel, int, int)
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
    MatNLElement *ele = new NLTetrahedral4(nn);
    ele->setMaterial(mat);
    ele->setGlNum(glNum);
    ele->setProp(prop);
    return new MatNLCorotator(ele);
#endif
  }
  else {
    return new TetCorotator(nn, prop->E, prop->nu, cs, prop->Ta, prop->W, prop->ymtt, prop->ctett);
  }
  printf("WARNING: Corotator not implemented for element %d\n", glNum+1); return 0;
}

int
Tetrahedral::getDecFace(int iFace, int *fn)
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

void
Tetrahedral::getCFrame(CoordSet &cs, double cFrame[3][3]) const
{
  if(Tetrahedral::cFrame) {
    cFrame[0][0] = Tetrahedral::cFrame[0]; cFrame[0][1] = Tetrahedral::cFrame[1]; cFrame[0][2] = Tetrahedral::cFrame[2];
    cFrame[1][0] = Tetrahedral::cFrame[3]; cFrame[1][1] = Tetrahedral::cFrame[4]; cFrame[1][2] = Tetrahedral::cFrame[5];
    cFrame[2][0] = Tetrahedral::cFrame[6]; cFrame[2][1] = Tetrahedral::cFrame[7]; cFrame[2][2] = Tetrahedral::cFrame[8];
  }
  else {
    cFrame[0][0] = cFrame[1][1] = cFrame[2][2] = 1.;
    cFrame[0][1] = cFrame[0][2] = cFrame[1][0] = cFrame[1][2] = cFrame[2][0] = cFrame[2][1] = 0.;
  }
}
