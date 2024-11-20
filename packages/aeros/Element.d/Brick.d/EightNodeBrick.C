#include <cstdio>
#include <iostream>
#include <cmath>

#include <Element.d/Brick.d/EightNodeBrick.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Corotational.d/BrickCorotator.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/NLHexahedral.h>
#include <Element.d/Utils.d/SolidElemUtils.h>
#include <Corotational.d/MatNLCorotator.h>
#include <Element.d/Helm.d/ARubberF.h>

extern "C" {
void _FORTRAN(brkcmt)(double&, double&, double*);

void _FORTRAN(brik8v)(double*, double*, double*, double*, const int&,
                      double*, int &);

void _FORTRAN(sands17)(const int&, double*, double*, double*, double*,
                       double*, double*, double*, const int&, const int&,
                       const int&, const int&, const int&);

void _FORTRAN(hxgaus)(int&, int&, int&, int&, int&, int&,
                      double&, double&, double&, double&);

void _FORTRAN(h8shpe)(double&, double&, double&, double*, double*,
                      double*, double*, double*, double*, double*,
                      double&);

void _FORTRAN(lgauss)(const int&, int&, double*, double*);
}

double Hexa8ShapeFct(double Shape[8], double DShape[8][3], double m[3], double X[8], double Y[8], double Z[8]);

EightNodeBrick::EightNodeBrick(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
  nn[2] = nodenums[2];
  nn[3] = nodenums[3];
  nn[4] = nodenums[4];
  nn[5] = nodenums[5];
  nn[6] = nodenums[6];
  nn[7] = nodenums[7];

  cFrame = 0;
  cCoefs = 0;
  mat = 0;
}

EightNodeBrick::~EightNodeBrick()
{
  if(cCoefs && mat) delete mat;
}

Element *
EightNodeBrick::clone()
{
  return new EightNodeBrick(*this);
}

void
EightNodeBrick::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
  nn[5] = table[nn[5]];
  nn[6] = table[nn[6]];
  nn[7] = table[nn[7]];
}

void
EightNodeBrick::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
  nn[5] = table[nn[5]];
  nn[6] = table[nn[6]];
  nn[7] = table[nn[7]];
}

void
EightNodeBrick::getVonMises(Vector& stress, Vector& weight, CoordSet &cs,
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

  double x[8], y[8], z[8];
  cs.getCoordinates(nn, numNodes(), x, y, z);

  double c[6][6];
  _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)c);

  int vmflg = 0, strainFlg = 0;
  // Flags sands17 to calculate Von Mises stress
  if(strInd == 6) vmflg = 1;
  // Flags sands17 to calculate Von Mises strain
  if(strInd == 13) strainFlg = 1;

  const int maxgus = 8;
  const int maxstr = 7;
  const int elm    = 1;
  const int msize  = 1;

  double elStress[maxgus][maxstr], elStrain[maxgus][maxstr];

  _FORTRAN(sands17)(elm, x, y, z, (double*)c, elDisp.data(),
                    (double*)elStress, (double*)elStrain,
                    maxgus, maxstr, msize, vmflg, strainFlg);

  if(strInd < 7) {
    double thermalStress[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
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
      thermalStress[6] = coef*(ndTemps[6]-Tref);
      thermalStress[7] = coef*(ndTemps[7]-Tref);
    }
    stress[0] = elStress[0][strInd] - thermalStress[0];
    stress[1] = elStress[1][strInd] - thermalStress[1];
    stress[2] = elStress[2][strInd] - thermalStress[2];
    stress[3] = elStress[3][strInd] - thermalStress[3];
    stress[4] = elStress[4][strInd] - thermalStress[4];
    stress[5] = elStress[5][strInd] - thermalStress[5];
    stress[6] = elStress[6][strInd] - thermalStress[6];
    stress[7] = elStress[7][strInd] - thermalStress[7];
  }
  else if(strInd < 14) {
    stress[0] = elStrain[0][strInd-7];
    stress[1] = elStrain[1][strInd-7];
    stress[2] = elStrain[2][strInd-7];
    stress[3] = elStrain[3][strInd-7];
    stress[4] = elStrain[4][strInd-7];
    stress[5] = elStrain[5][strInd-7];
    stress[6] = elStrain[6][strInd-7];
    stress[7] = elStrain[7][strInd-7];
  }
  else {
    stress[0] = 0;
    stress[1] = 0;
    stress[2] = 0;
    stress[3] = 0;
    stress[4] = 0;
    stress[5] = 0;
    stress[6] = 0;
    stress[7] = 0;
  }
}

void
EightNodeBrick::getAllStress(FullM& stress, Vector& weight, CoordSet &cs,
                             Vector& elDisp, int strInd, int surface, double *ndTemps)
{
  if(cCoefs) {
    getAllStressAniso(stress, weight, cs,
                      elDisp, strInd, surface, ndTemps);
    return;
  }

  weight = 1.0;

  double x[8], y[8], z[8];
  cs.getCoordinates(nn, numNodes(), x, y, z);

  double c[6][6];
  _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)c);

  int vmflg, strainFlg;
  vmflg = 0;
  strainFlg = 0;

  const int maxgus = 8;
  const int maxstr = 7;
  const int elm    = 1;
  const int msize  = 1;

  double elStress[maxgus][maxstr], elStrain[maxgus][maxstr];

  _FORTRAN(sands17)(elm, x, y, z, (double*)c, elDisp.data(),
                    (double*)elStress, (double*)elStrain,
                    maxgus, maxstr, msize, vmflg, strainFlg);

  // Store all Stress or all Strain as defined by strInd
  int i,j;
  if(strInd == 0) {
    double thermalStress[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
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
      thermalStress[6] = coef*(ndTemps[6]-Tref);
      thermalStress[7] = coef*(ndTemps[7]-Tref);
    }
    for (i=0; i<8; ++i) {
      for (j=0; j<3; ++j) {
        stress[i][j] = elStress[i][j] - thermalStress[i];
      }
      for (j=3; j<6; ++j) {
        stress[i][j] = elStress[i][j];
      }
    }
  }
  else {
    for (i=0; i<8; ++i) {
      for (j=0; j<6; ++j) {
        stress[i][j] = elStrain[i][j];
      }
    }
  }

  // Get Element Principals for each node without averaging
  double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double pvec[3] = {0.0,0.0,0.0};

  for (i=0; i<8; ++i) {
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
EightNodeBrick::getMass(const CoordSet& cs) const
{
  const int nnodes = 8;
  const int numgauss = 2;
 
  double X[8], Y[8], Z[8]; 
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  double m[3], Shape[8], DShape[8][3];
  double wx, wy, wz, dOmega;
  double volume = 0.0;
  for(int i = 1; i <= numgauss; i++) {
    _FORTRAN(lgauss)(numgauss, i, &m[0], &wx);
    for(int j = 1; j <= numgauss; j++) {
      _FORTRAN(lgauss)(numgauss, j, &m[1], &wy);
      for(int k = 1; k <= numgauss; k++) {
        _FORTRAN(lgauss)(numgauss, k, &m[2], &wz);
        dOmega = Hexa8ShapeFct(Shape, DShape, m, X, Y, Z);
        volume += fabs(dOmega)*wx*wy*wz;
      }
    }
  }

  return prop->rho*volume;
}

double
EightNodeBrick::weight(CoordSet& cs, double *gravityAcceleration, int altitude_direction)
{
  if (prop == NULL) {
    return 0.0;
  }

  double _mass = getMass(cs);
  return _mass*gravityAcceleration[altitude_direction];
}

void
EightNodeBrick::getGravityForce(CoordSet& cs, double *gravityAcceleration, 
                                Vector& gravityForce, int gravflg, GeomState *geomState)
{
  const int nnodes = 8;

  // Lumped
  if (gravflg != 2) {

    double totmas = getMass(cs);

    double m[24*24];
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

    double X[8], Y[8], Z[8];
    cs.getCoordinates(nn, nnodes, X, Y, Z);

    double lforce[8];
    for(int i = 0; i < nnodes; ++i) lforce[i] = 0.0;

    // integration: loop over Gauss pts
    int numgauss = 2;

    for (int pt1 = 1; pt1 <= numgauss; pt1++) {
      for (int pt2 = 1; pt2 <= numgauss; pt2++) {
        for (int pt3 = 1; pt3 <= numgauss; pt3++) {
          // get gauss point
          double xi, eta, mu, wt;
          _FORTRAN(hxgaus)(numgauss, pt1, numgauss, pt2, numgauss, pt3, xi, eta, mu, wt);

          //compute shape functions
          double shapeFunc[8], shapeGradX[8], shapeGradY[8], shapeGradZ[8];
          double detJ; // det of jacobian

          _FORTRAN(h8shpe)(xi, eta, mu, X, Y, Z,
                           shapeFunc, shapeGradX, shapeGradY, shapeGradZ, detJ);

          for(int n = 0; n < 8; ++n)
            lforce[n] += wt*shapeFunc[n]*detJ;
        }
      }
    }

    for(int i = 0; i < nnodes; ++i) {
      gravityForce[3*i+0] = lforce[i]*prop->rho*gravityAcceleration[0];
      gravityForce[3*i+1] = lforce[i]*prop->rho*gravityAcceleration[1];
      gravityForce[3*i+2] = lforce[i]*prop->rho*gravityAcceleration[2];
    }
  }
}

void
EightNodeBrick::getThermalForce(CoordSet &cs, Vector &ndTemps,
                                Vector &elementThermalForce, int glflag,
                                GeomState *geomState)
{
  // ASSUME CONSTANT THERMAL EXPANSION COEFF. & REFERENCE TEMPERATURE OVER THE ELEMENT 
  const int nnodes = 8;
  const int ndofs = 24;

  // initialize nodal thermal forces
  for(int i=0; i<ndofs; i++) elementThermalForce[i] = 0.0;

  // for nonlinear analyses, the thermal load for this element is now computed in getStiffAndForce
  if(geomState) return;

  double X[8], Y[8], Z[8];
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
  const int numgauss = 2;
  double Shape[8], DShape[8][3], m[3];
  double wx,wy,wz,w,J; 

  for(int i=1; i<=numgauss; i++) {
    _FORTRAN(lgauss)(numgauss,i,&m[0],&wx);
    for(int j=1; j<=numgauss; j++) {
      _FORTRAN(lgauss)(numgauss,j,&m[1],&wy);
      for(int k=1; k<=numgauss; k++) {
        _FORTRAN(lgauss)(numgauss,k,&m[2],&wz);
        J = Hexa8ShapeFct(Shape, DShape, m, X, Y, Z);
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
EightNodeBrick::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
  const int nnodes = 8;
  const int ndofs = 24;

  double X[8], Y[8], Z[8];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  FullSquareMatrix M(ndofs,mel);

  if(cmflg) { // consistent mass matrix
    M.zero();
    int ls[24] = {0,3,6,9,12,15,18,21,
                  1,4,7,10,13,16,19,22,
                  2,5,8,11,14,17,20,23};

    // integration: loop over Gauss pts
    const int numgauss = 2;
    double m[3], Shape[8], DShape[8][3];
    double wx, wy, wz, w;
    double dOmega; // det of jacobian
    int jSign = 0;

    for(int i = 1; i <= numgauss; i++) {
      _FORTRAN(lgauss)(numgauss, i, &m[0], &wx);
      for(int j = 1; j <= numgauss; j++) {
        _FORTRAN(lgauss)(numgauss, j, &m[1], &wy);
        for(int k = 1; k <= numgauss; k++) {
          _FORTRAN(lgauss)(numgauss, k, &m[2], &wz);
          dOmega = Hexa8ShapeFct(Shape, DShape, m, X, Y, Z);
#ifdef CHECK_JACOBIAN
          checkJacobian(&dOmega, &jSign, "EightNodeBrick::massMatrix");
#endif
          w = fabs(dOmega)*wx*wy*wz*prop->rho;
          addNtDNtoM3DSolid(M, Shape, w, nnodes, ls);
        }
      }
    }
  }
  else { // Lumped mass matrix
    fprintf(stderr," *** In EightNodeBrick::massMatrix: Lumped mass matrix NOT implemented. Abort.\n");
    exit(-1);
  }

  return M;
}

FullSquareMatrix
EightNodeBrick::stiffness(const CoordSet &cs, double *d, int flg) const
{
  const int nnodes = 8;
  const int ndofs = 24;

  double X[8], Y[8], Z[8];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  FullSquareMatrix K(ndofs,d);

  // get constitutive matrix
  double C[6][6];
  if(cCoefs) { // anisotropic material
    // transform local constitutive matrix to global frame
    rotateConstitutiveMatrix(cCoefs, cFrame, C);
  } else // isotropic material
    _FORTRAN(brkcmt)(prop->E, prop->nu, (double*)C);

  const int numgauss = 2;
  int status = 0;

  _FORTRAN(brik8v)(X, Y, Z, (double *)C, numgauss, (double *)d, status);

  return K;
}

void EightNodeBrick::aRubberStiffnessDerivs(CoordSet & cs, complex<double> *d,
                                            int n, double omega) {
  const int nnodes = 8;
  const int ndofs = 24;

  double X[8], Y[8], Z[8];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  double C[6][6];
  for(int i=0;i<6;i++) for(int j=0;j<6;j++) C[i][j] = 0.0;

  double *Km = new double[ndofs*ndofs];
  C[0][0] = C[1][1] = C[2][2] = 2.0;
  C[3][3] = C[4][4] = C[5][5] = 1.0;
  const int numgauss = 2;
  int status = 0;
  _FORTRAN(brik8v)(X, Y, Z, (double *)C, numgauss, (double *)Km, status);

  double *Kl = new double[ndofs*ndofs];
  C[0][0] = C[1][1] = C[2][2] = 1.0;
  C[0][1] = C[1][0] = C[0][2] = C[2][0] = C[1][2] = C[2][1] = 1.0;
  C[3][3] = C[4][4] = C[5][5] = 0.0;
  _FORTRAN(brik8v)(X, Y, Z, (double *)C, numgauss, (double *)Kl, status);

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
EightNodeBrick::numNodes() const
{
  return 8;
}

int
EightNodeBrick::numDofs() const
{
  return 24;
}

int
EightNodeBrick::getTopNumber() const
{
  return 117;
}

int*
EightNodeBrick::nodes(int *p) const
{
  if(!p) p = new int[numNodes()];
  p[0] = nn[0];
  p[1] = nn[1];
  p[2] = nn[2];
  p[3] = nn[3];
  p[4] = nn[4];
  p[5] = nn[5];
  p[6] = nn[6];
  p[7] = nn[7];
  return p;
}

int*
EightNodeBrick::dofs(DofSetArray &dsa, int *p) const
{
  if(!p) p = new int[numDofs()];

  dsa.number(nn[0], DofSet::XYZdisp, p  );
  dsa.number(nn[1], DofSet::XYZdisp, p+3);
  dsa.number(nn[2], DofSet::XYZdisp, p+6);
  dsa.number(nn[3], DofSet::XYZdisp, p+9);
  dsa.number(nn[4], DofSet::XYZdisp, p+12);
  dsa.number(nn[5], DofSet::XYZdisp, p+15);
  dsa.number(nn[6], DofSet::XYZdisp, p+18);
  dsa.number(nn[7], DofSet::XYZdisp, p+21);

  return p;
}

void
EightNodeBrick::markDofs(DofSetArray &dsa) const
{
  dsa.mark(nn, numNodes(), DofSet::XYZdisp);
}

// Stress evaluation in case of anisotropic elastic constitutive matrix
void
EightNodeBrick::getVonMisesAniso(Vector &stress, Vector &weight, CoordSet &cs,
                                 Vector &elDisp, int strInd, int surface, double *ndTemps,
                                 double ylayer, double zlayer, int avgnum)
{
  const int nnodes = 8;
  weight = 1.0;
  
  double X[8], Y[8], Z[8];
  cs.getCoordinates(nn, nnodes, X, Y, Z);
     
  // Flags sands17 to calculate Von Mises stress and/or Von Mises strain
  bool vmflg     = (strInd==6) ? true : false;
  bool strainFlg = (strInd==13)? true : false;
  bool meanVms   = false; // This can be used to force averaging of the Von Mises stress & strain.

  double elStress[8][7];
  double elStrain[8][7];
 
  // get constitutive matrix and coefficients of thermal expansion
  double C[6][6], alpha[6];
  // transform local constitutive matrix to global frame
  rotateConstitutiveMatrix(cCoefs, cFrame, C);
  // transform local coefficients of thermal expansion to global frame
  if(ndTemps) rotateVector(cCoefs+36, cFrame, alpha);
 
  // Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[8][3] = {{-1.0,-1.0,-1.0},{1.0,-1.0,-1.0},{ 1.0,1.0,-1.0},{-1.0,1.0,-1.0},
                               {-1.0,-1.0, 1.0},{1.0,-1.0, 1.0},{ 1.0,1.0, 1.0},{-1.0,1.0, 1.0}};

  double Shape[8], DShape[8][3];
  for(int inode=0; inode<nnodes; inode++) {
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Hexa8ShapeFct(Shape, DShape, m, X, Y, Z); 
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
EightNodeBrick::getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
                                  Vector &elDisp, int strInd, int surface, double *ndTemps)
{
  const int nnodes = 8;
  weight = 1.0;
  
  double X[8], Y[8], Z[8];
  cs.getCoordinates(nn, nnodes, X, Y, Z);

  double elStress[8][6];
  double elStrain[8][6];
 
  // get constitutive matrix and coefficients of thermal expansion
  double C[6][6], alpha[6];
  // transform local constitutive matrix to global frame
  rotateConstitutiveMatrix(cCoefs, cFrame, C);
  // transform local coefficients of thermal expansion to global frame
  if(ndTemps) rotateVector(cCoefs+36, cFrame, alpha);
 
  // Loop over nodes -> compute nodal strains & stresses
  double nodeRefCoord[8][3] = {{-1.0,-1.0,-1.0},{1.0,-1.0,-1.0},{ 1.0,1.0,-1.0},{-1.0,1.0,-1.0},
                               {-1.0,-1.0, 1.0},{1.0,-1.0, 1.0},{ 1.0,1.0, 1.0},{-1.0,1.0, 1.0}};

  double Shape[8], DShape[8][3];
  for(int inode=0; inode<nnodes; inode++) {
    // compute shape fcts & their derivatives at node
    double* m = &nodeRefCoord[inode][0];
    Hexa8ShapeFct(Shape, DShape, m, X, Y, Z); 
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
EightNodeBrick::setMaterial(NLMaterial *_mat)
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
EightNodeBrick::numStates()
{
#ifdef USE_EIGEN3
  int numGaussPoints = NLHexahedral8::numGaussPoints;
  return (mat) ? numGaussPoints*mat->getNumStates(): 0;
#else
  return 0;
#endif
}

void
EightNodeBrick::initStates(double *st)
{
#ifdef USE_EIGEN3
  if(mat) {
    int ninterns = mat->getNumStates();
    int numGaussPoints = NLHexahedral8::numGaussPoints;

    for(int i = 0; i < numGaussPoints; ++i)
      mat->initStates(st+i*ninterns);
  }
#endif
}

Corotator *
EightNodeBrick::getCorotator(CoordSet &cs, double *kel, int, int)
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
    MatNLElement *ele = new NLHexahedral8(nn);
    ele->setMaterial(mat);
    ele->setGlNum(glNum);
    ele->setProp(prop);
    return new MatNLCorotator(ele);
#endif
  }
  else {
    return new BrickCorotator(nn, prop->E, prop->nu, cs, prop->Ta, prop->W, prop->ymtt, prop->ctett);
  }
  printf("WARNING: Corotator not implemented for element %d\n", glNum+1); return 0;
}

int
EightNodeBrick::getDecFace(int iFace, int *fn)
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
EightNodeBrick::getFace(int iFace, int *fn)
{
  // note: all face normals are outward pointing.
  switch(iFace) {
    case 0: fn[0] = nn[3]; fn[1] = nn[2]; fn[2] = nn[1]; fn[3] = nn[0]; break;
    case 1: fn[0] = nn[4]; fn[1] = nn[5]; fn[2] = nn[6]; fn[3] = nn[7]; break;
    case 2: fn[0] = nn[3]; fn[1] = nn[0]; fn[2] = nn[4]; fn[3] = nn[7]; break;
    case 3: fn[0] = nn[0]; fn[1] = nn[1]; fn[2] = nn[5]; fn[3] = nn[4]; break;
    case 4: fn[0] = nn[6]; fn[1] = nn[5]; fn[2] = nn[1]; fn[3] = nn[2]; break;
    default:
    case 5: fn[0] = nn[7]; fn[1] = nn[6]; fn[2] = nn[2]; fn[3] = nn[3]; break;
  }
  return 4;
}

void
EightNodeBrick::getCFrame(CoordSet &cs, double cFrame[3][3]) const
{
  if(EightNodeBrick::cFrame) {
    cFrame[0][0] = EightNodeBrick::cFrame[0]; cFrame[0][1] = EightNodeBrick::cFrame[1]; cFrame[0][2] = EightNodeBrick::cFrame[2];
    cFrame[1][0] = EightNodeBrick::cFrame[3]; cFrame[1][1] = EightNodeBrick::cFrame[4]; cFrame[1][2] = EightNodeBrick::cFrame[5];
    cFrame[2][0] = EightNodeBrick::cFrame[6]; cFrame[2][1] = EightNodeBrick::cFrame[7]; cFrame[2][2] = EightNodeBrick::cFrame[8];
  }
  else {
    cFrame[0][0] = cFrame[1][1] = cFrame[2][2] = 1.;
    cFrame[0][1] = cFrame[0][2] = cFrame[1][0] = cFrame[1][2] = cFrame[2][0] = cFrame[2][1] = 0.;
  }
}
