#include <cmath>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Utils.d/MFTT.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Element.d/Element.h>
#include <Corotational.d/TetCorotator.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>

extern "C" {
  void  _FORTRAN(vmelmv)(double*, int &, int &, int &, int &, int &);
  void  _FORTRAN(strainvm)(double*, int &, int &, int &, int &);
};

TetCorotator::TetCorotator(int nodeNumbers[4], double _em, double _nu,
                           CoordSet& cs, double _Tref, double _alpha,
                           MFTTData *_ymtt, MFTTData *_ctett)
{
  nodeNum[0] = nodeNumbers[0];
  nodeNum[1] = nodeNumbers[1];
  nodeNum[2] = nodeNumbers[2];
  nodeNum[3] = nodeNumbers[3];

  em    = _em;    // Elastic modulus
  nu    = _nu;    // Poisson's ratio
  Tref  = _Tref;  // Ambient temperature
  alpha = _alpha; // Thermal expansion coefficient
  ymtt = _ymtt;
  ctett = _ctett;
}

// geomState -> contains the updated nodal coordinates
// cs        -> contains the original nodal coordinates
void
TetCorotator::getStiffAndForce(GeomState &geomState, CoordSet &cs,
                               FullSquareMatrix &K, double *f, double dt, double t)
{
  int i,j,k;
  double nGrad[4][3];

  // get the nodal temperatures
  Vector ndTemps(4);
  geomState.get_temperature(4, nodeNum, ndTemps, Tref);
  
  // compute dN_i/dX_j, also obtain dOmega
  double dOmega;

  // dOmega is here off by a factor of 1/4
  // this divides all by a factor of 4 later
  dOmega = computeShapeGrad(cs, nGrad)/24;
 
  // now get F_ij = dPhi_i/dX_j = x^k_i dN_k/dX_j
  double F[3][3];
  
  for(j = 0; j < 3; ++j) {
    F[0][j] = geomState[nodeNum[0]].x*nGrad[0][j]
            + geomState[nodeNum[1]].x*nGrad[1][j]
            + geomState[nodeNum[2]].x*nGrad[2][j]
            + geomState[nodeNum[3]].x*nGrad[3][j];

    F[1][j] = geomState[nodeNum[0]].y*nGrad[0][j]
            + geomState[nodeNum[1]].y*nGrad[1][j]
            + geomState[nodeNum[2]].y*nGrad[2][j]
            + geomState[nodeNum[3]].y*nGrad[3][j];
    
    F[2][j] = geomState[nodeNum[0]].z*nGrad[0][j]
            + geomState[nodeNum[1]].z*nGrad[1][j]
            + geomState[nodeNum[2]].z*nGrad[2][j]
            + geomState[nodeNum[3]].z*nGrad[3][j];
  }
    
  // compute e_ij = 0.5*(F_ki Fkj - delta_ij)
  // here these are off by factor of 2
  double e_11 = (F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][0]-1.0);
  double e_22 = (F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1]-1.0);
  double e_33 = (F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2]-1.0);
  double e_12 = (F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1]);
  double e_13 = (F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2]);
  double e_23 = (F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2]);

  // Subtract thermal strain (off by factor of 2)
  double theta = (ndTemps[0]+ndTemps[1]+ndTemps[2]+ndTemps[3])/4 - Tref;
  double alpha = (ctett) ? ctett->getValAlt(theta) : TetCorotator::alpha;
  e_11 -= 2*alpha*theta;
  e_22 -= 2*alpha*theta;
  e_33 -= 2*alpha*theta;

  double sigma[6];
  double em = (ymtt) ? ymtt->getValAlt(theta) : TetCorotator::em;
  double E  = em;
  double E2 = E*nu/((1+nu)*(1-2*nu));
  double E1 = E2+E/(1+nu);
  // no factor of 1/2 on G2 due to using tensor strain
  double G2 = E/(1+nu);
  // these here are off by a factor of 2
  sigma[0] = E1*e_11+E2*(e_22+e_33);
  sigma[1] = E1*e_22+E2*(e_11+e_33);
  sigma[2] = E1*e_33+E2*(e_11+e_22);
  // these here are off by a factor of 4
  sigma[3] = 2*G2*e_12;
  sigma[4] = 2*G2*e_13;
  sigma[5] = 2*G2*e_23;

  // Compute de_ij/dUl for the symmetric part
  // First we get dF_ij/dUl in a very compact form.
  // dF_ij/dUl = dN_p/dX_j delta_iq; with p = int(l/3)+1 and l-1=q-1 mod(3)
  // this means that dF_ij/dUl is already contained in dN_k/dX_j
  double dedU[12][6];
  for(i = 0; i < 4; ++i)
    for(j = 0; j < 3; ++j) {
      dedU[3*i+j][0] = 2*nGrad[i][0]*F[j][0];
      dedU[3*i+j][1] = 2*nGrad[i][1]*F[j][1];
      dedU[3*i+j][2] = 2*nGrad[i][2]*F[j][2];
      dedU[3*i+j][3] = (nGrad[i][0]*F[j][1]+nGrad[i][1]*F[j][0]);
      dedU[3*i+j][4] = (nGrad[i][0]*F[j][2]+nGrad[i][2]*F[j][0]);
      dedU[3*i+j][5] = (nGrad[i][1]*F[j][2]+nGrad[i][2]*F[j][1]);
    }
  // all dedU terms are here off by a factor of 2

  // Get the force:
  for(i = 0; i < 12; ++i) {
    f[i] = dOmega*( dedU[i][0]*sigma[0] +
                    dedU[i][1]*sigma[1] +
                    dedU[i][2]*sigma[2] +
                    dedU[i][3]*sigma[3] +
                    dedU[i][4]*sigma[4] +
                    dedU[i][5]*sigma[5]);
  }
  // normal terms of f are
  // (1/4)*(2*2) = 1 => right level
  // shear terms of f are
  // (1/4)*(2*4)*(1/2 for off-diag mult) = 1 => right level

  // now get ds_ij/dUl
  double dsdU[12][6];
  for(i = 0; i < 12; ++i) {
    dsdU[i][0] = E1*dedU[i][0]+E2*(dedU[i][1]+dedU[i][2]);
    dsdU[i][1] = E1*dedU[i][1]+E2*(dedU[i][0]+dedU[i][2]);
    dsdU[i][2] = E1*dedU[i][2]+E2*(dedU[i][0]+dedU[i][1]);
    dsdU[i][3] = 2*G2*dedU[i][3];
    dsdU[i][4] = 2*G2*dedU[i][4];
    dsdU[i][5] = 2*G2*dedU[i][5];
  }
  // normal terms are here off by a factor of 2
  // shear terms are here off by a factor of 2*2 = 4

  // multiply modified dsdU by dedU. Only do the symmetric part
  for(i = 0; i < 12; ++i)
    for(j = 0; j <= i; ++j)
      K[j][i] = dsdU[i][0]*dedU[j][0] + 
                dsdU[i][1]*dedU[j][1] +
                dsdU[i][2]*dedU[j][2] +
                dsdU[i][3]*dedU[j][3] +
                dsdU[i][4]*dedU[j][4] +
                dsdU[i][5]*dedU[j][5];
  // normal terms are here off by a factor of 2*2=4
  // shear terms are here off by a factor of 4*2 = 8
  //                        *(1/2 for off-diag mult) = 4

  // add s*d2e/dU_idUj (symmetric part only)
  for(i = 0; i < 4; ++i)
    for(k = 0; k <= i; ++k)
      for(j = 0; j < 3; ++j)
        K[3*k+j][3*i+j] +=
                sigma[0]*(2*nGrad[i][0]*nGrad[k][0]) +
                sigma[1]*(2*nGrad[i][1]*nGrad[k][1]) +
                sigma[2]*(2*nGrad[i][2]*nGrad[k][2]) +
                sigma[3]*(nGrad[i][0]*nGrad[k][1]+nGrad[i][1]*nGrad[k][0]) +
                sigma[4]*(nGrad[i][0]*nGrad[k][2]+nGrad[i][2]*nGrad[k][0]) +
                sigma[5]*(nGrad[i][1]*nGrad[k][2]+nGrad[i][2]*nGrad[k][1]);
  // normal terms are here off by a factor of 2*2=4
  // shear terms are here off by a factor of 4*2 = 8
  //                        *(1/2 for off-diag mult) = 4

  for(i = 0; i < 12; ++i)
    for(j = 0; j <= i; ++j)
      K[i][j] = (K[j][i] *= dOmega);
// all terms are 1/4*4 = 1 => right level
}

void
TetCorotator::getInternalForce(GeomState &geomState, CoordSet &cs,
                               FullSquareMatrix &, double *f, double dt, double t)
{
  int i,j,k;
  double nGrad[4][3];

  // get the nodal temperatures
  Vector ndTemps(4);
  geomState.get_temperature(4, nodeNum, ndTemps, Tref);
  
  // compute dN_i/dX_j, also obtain dOmega 
  double dOmega;

  // dOmega is here off by a factor of 1/4
  // this divides all by a factor of 4 later
  dOmega = computeShapeGrad(cs, nGrad)/24;
 
  // now get F_ij = dPhi_i/dX_j = x^k_i dN_k/dX_j
  double F[3][3];
  
  for(j = 0; j < 3; ++j) {
    F[0][j] = geomState[nodeNum[0]].x*nGrad[0][j]
            + geomState[nodeNum[1]].x*nGrad[1][j]
            + geomState[nodeNum[2]].x*nGrad[2][j]
            + geomState[nodeNum[3]].x*nGrad[3][j];

    F[1][j] = geomState[nodeNum[0]].y*nGrad[0][j]
            + geomState[nodeNum[1]].y*nGrad[1][j]
            + geomState[nodeNum[2]].y*nGrad[2][j]
            + geomState[nodeNum[3]].y*nGrad[3][j];
    
    F[2][j] = geomState[nodeNum[0]].z*nGrad[0][j]
            + geomState[nodeNum[1]].z*nGrad[1][j]
            + geomState[nodeNum[2]].z*nGrad[2][j]
            + geomState[nodeNum[3]].z*nGrad[3][j];
  }
    
  // compute e_ij = 0.5*(F_ki Fkj - delta_ij)
  // here these are off by factor of 2
  double e_11 = (F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][0]-1.0);
  double e_22 = (F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1]-1.0);
  double e_33 = (F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2]-1.0);
  double e_12 = (F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1]);
  double e_13 = (F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2]);
  double e_23 = (F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2]);

  // Subtract thermal strain (off by factor of 2)
  double theta = (ndTemps[0]+ndTemps[1]+ndTemps[2]+ndTemps[3])/4 - Tref;
  double alpha = (ctett) ? ctett->getValAlt(theta) : TetCorotator::alpha;
  e_11 -= 2*alpha*theta;
  e_22 -= 2*alpha*theta;
  e_33 -= 2*alpha*theta;

  double sigma[6];
  double em = (ymtt) ? ymtt->getValAlt(theta) : TetCorotator::em;
  double E  = em;
  double E2 = E*nu/((1+nu)*(1-2*nu));
  double E1 = E2+E/(1+nu);
  // no factor of 1/2 on G2 due to using tensor strain
  double G2 = E/(1+nu);
  // these here are off by a factor of 2
  sigma[0] = E1*e_11+E2*(e_22+e_33);
  sigma[1] = E1*e_22+E2*(e_11+e_33);
  sigma[2] = E1*e_33+E2*(e_11+e_22);
  // these here are off by a factor of 4
  sigma[3] = 2*G2*e_12;
  sigma[4] = 2*G2*e_13;
  sigma[5] = 2*G2*e_23;

  // Compute de_ij/dUl for the symmetric part
  // First we get dF_ij/dUl in a very compact form.
  // dF_ij/dUl = dN_p/dX_j delta_iq; with p = int(l/3)+1 and l-1=q-1 mod(3)
  // this means that dF_ij/dUl is already contained in dN_k/dX_j
  double dedU[12][6];
  for(i = 0; i < 4; ++i)
    for(j = 0; j < 3; ++j) {
      dedU[3*i+j][0] = 2*nGrad[i][0]*F[j][0];
      dedU[3*i+j][1] = 2*nGrad[i][1]*F[j][1];
      dedU[3*i+j][2] = 2*nGrad[i][2]*F[j][2];
      dedU[3*i+j][3] = (nGrad[i][0]*F[j][1]+nGrad[i][1]*F[j][0]);
      dedU[3*i+j][4] = (nGrad[i][0]*F[j][2]+nGrad[i][2]*F[j][0]);
      dedU[3*i+j][5] = (nGrad[i][1]*F[j][2]+nGrad[i][2]*F[j][1]);
    }
  // all dedU terms are here off by a factor of 2

  // Get the force:
  for(i = 0; i < 12; ++i) {
    f[i] = dOmega*( dedU[i][0]*sigma[0] +
                    dedU[i][1]*sigma[1] +
                    dedU[i][2]*sigma[2] +
                    dedU[i][3]*sigma[3] +
                    dedU[i][4]*sigma[4] +
                    dedU[i][5]*sigma[5]);
  }
  // normal terms of f are
  // (1/4)*(2*2) = 1 => right level
  // shear terms of f are
  // (1/4)*(2*4)*(1/2 for off-diag mult) = 1 => right level
}

//-----------------------------------------------------------------------------------

void
TetCorotator::computeStrainGrad(GeomState &geomState, double nGrad[4][3], 
                                double dedU[12][6])
{
  double F[3][3];

  int i, j;
  for(j = 0; j < 3; ++j) {
    F[0][j] = geomState[nodeNum[0]].x*nGrad[0][j]
            + geomState[nodeNum[1]].x*nGrad[1][j]
            + geomState[nodeNum[2]].x*nGrad[2][j]
            + geomState[nodeNum[3]].x*nGrad[3][j];

    F[1][j] = geomState[nodeNum[0]].y*nGrad[0][j]
            + geomState[nodeNum[1]].y*nGrad[1][j]
            + geomState[nodeNum[2]].y*nGrad[2][j]
            + geomState[nodeNum[3]].y*nGrad[3][j];

    F[2][j] = geomState[nodeNum[0]].z*nGrad[0][j]
            + geomState[nodeNum[1]].z*nGrad[1][j]
            + geomState[nodeNum[2]].z*nGrad[2][j]
            + geomState[nodeNum[3]].z*nGrad[3][j];
  }

  for(i = 0; i < 4; ++i)
    for(j = 0; j < 3; ++j) {
      dedU[3*i+j][0] = nGrad[i][0]*F[j][0];
      dedU[3*i+j][1] = nGrad[i][1]*F[j][1];
      dedU[3*i+j][2] = nGrad[i][2]*F[j][2];
      dedU[3*i+j][3] = 0.5*(nGrad[i][0]*F[j][1]+nGrad[i][1]*F[j][0]);
      dedU[3*i+j][4] = 0.5*(nGrad[i][0]*F[j][2]+nGrad[i][2]*F[j][0]);
      dedU[3*i+j][5] = 0.5*(nGrad[i][1]*F[j][2]+nGrad[i][2]*F[j][1]);
    }

}

//----------------------------------------------------------------------------------- 

void
TetCorotator::getNLVonMises(Vector& stress, Vector& weight, GeomState &geomState,
                            GeomState *, CoordSet& cs, int strInd, int,
                            double, double, int, int)
{
  weight = 1.0;

  int maxstr = 7;
  int elm    = 1;
  int nno    = 4;

  double elStress[4][7];
  double elStrain[4][7];

  // Compute NL Stress/Strain
  computePiolaStress(geomState, cs, elStress, elStrain);

  // Compute Von Mises
  if(strInd == 6)
    _FORTRAN(vmelmv)((double*)elStress,nno,maxstr,elm,elm,nno);
  if(strInd == 13)
    _FORTRAN(strainvm)((double*)elStrain,nno,maxstr,elm,nno);

  // Store Stress or Strain as defined by strInd
  if(strInd < 7) {
    stress[0] = elStress[0][strInd];
    stress[1] = elStress[1][strInd];
    stress[2] = elStress[2][strInd];
    stress[3] = elStress[3][strInd];
  } else if(strInd < 14) {
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
TetCorotator::getNLAllStress(FullM &stress, Vector &weight, GeomState &geomState,
                             GeomState *, CoordSet &cs, int strInd, int, int)
{
  weight = 1.0;

  int i,j;

  double elStress[4][7];
  double elStrain[4][7];

  // Compute NL Stress/Strain
  computePiolaStress(geomState, cs, elStress, elStrain);

  // Store all Stress or all Strain as defined by strInd
  if(strInd == 0) {
    for (i=0; i<4; ++i) {
      for (j=0; j<6; ++j) {
        stress[i][j] = elStress[i][j];
      }
    }
  } else {
    for (i=0; i<4; ++i) {
      for (j=0; j<6; ++j) {
        stress[i][j] = elStrain[i][j];
      }
    }
  }

  // Get Element Principals without averaging
  double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double pvec[3] = {0.0,0.0,0.0};
  for (i=0; i<4; ++i) {
    for (j=0; j<6; ++j)
      svec[j] = stress[i][j];

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
TetCorotator::computePiolaStress(GeomState &geomState, CoordSet &cs,
                                 double stress[4][7], double strain[4][7])
{
  int i,j;
  double nGrad[4][3];

  // get the nodal temperatures
  Vector ndTemps(4);
  geomState.get_temperature(4, nodeNum, ndTemps, Tref);
  
  double dOmega = computeShapeGrad(cs, nGrad)/6;
 
  // now get F_ij = dPhi_i/dX_j = x^k_i dN_k/dX_j
  double F[3][3];
  
  for(j = 0; j < 3; ++j) {
    F[0][j] = geomState[nodeNum[0]].x*nGrad[0][j]
            + geomState[nodeNum[1]].x*nGrad[1][j]
            + geomState[nodeNum[2]].x*nGrad[2][j]
            + geomState[nodeNum[3]].x*nGrad[3][j];

    F[1][j] = geomState[nodeNum[0]].y*nGrad[0][j]
            + geomState[nodeNum[1]].y*nGrad[1][j]
            + geomState[nodeNum[2]].y*nGrad[2][j]
            + geomState[nodeNum[3]].y*nGrad[3][j];
    
    F[2][j] = geomState[nodeNum[0]].z*nGrad[0][j]
            + geomState[nodeNum[1]].z*nGrad[1][j]
            + geomState[nodeNum[2]].z*nGrad[2][j]
            + geomState[nodeNum[3]].z*nGrad[3][j];
  }

  // compute e_ij = 1/2(F_ki Fkj - delta_ij)
  double e_11 = 0.5*(F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][0]-1.0);
  double e_22 = 0.5*(F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1]-1.0);
  double e_33 = 0.5*(F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2]-1.0);
  // here engineering strain is computed
  double e_12 = (F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1]);
  double e_13 = (F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2]);
  double e_23 = (F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2]);

  // Reorder strain
  for(i=0; i<4; ++i) {
    strain[i][0] = e_11;
    strain[i][1] = e_22;
    strain[i][2] = e_33;
    strain[i][3] = e_12;
    strain[i][4] = e_23;
    strain[i][5] = e_13;
  }

  for(i=0; i<4; ++i) {
    // Subtract thermal strain
    double alpha = (ctett) ? ctett->getValAlt(ndTemps[i]) : TetCorotator::alpha;
    double e_11_m = e_11 - alpha*(ndTemps[i]-Tref);
    double e_22_m = e_22 - alpha*(ndTemps[i]-Tref);
    double e_33_m = e_33 - alpha*(ndTemps[i]-Tref);

    double sigma[6];
    double em = (ymtt) ? ymtt->getValAlt(ndTemps[i]) : TetCorotator::em;
    double E2 = em*nu/((1+nu)*(1-2*nu));
    double G2 = em/(2*(1+nu));
    double E1 = E2+em/(1+nu);

    sigma[0] = E1*e_11_m+E2*(e_22_m+e_33_m);
    sigma[1] = E1*e_22_m+E2*(e_11_m+e_33_m);
    sigma[2] = E1*e_33_m+E2*(e_11_m+e_22_m);
    sigma[3] = G2*e_12;
    sigma[4] = G2*e_23;
    sigma[5] = G2*e_13;

    // Reorder stress
    stress[i][0] = sigma[0];
    stress[i][1] = sigma[1];
    stress[i][2] = sigma[2];
    stress[i][3] = sigma[3];
    stress[i][4] = sigma[4];
    stress[i][5] = sigma[5];
  }
}

double TetCorotator::computeShapeGrad(CoordSet &nodes, double nGrad[4][3])
{
  double jac[3][3];

  //Jacobian
  // J_ij = dx_i/dxi_j
  jac[0][0] = nodes[ nodeNum[1] ]->x - nodes[ nodeNum[0] ]->x;
  jac[0][1] = nodes[ nodeNum[2] ]->x - nodes[ nodeNum[0] ]->x;
  jac[0][2] = nodes[ nodeNum[3] ]->x - nodes[ nodeNum[0] ]->x;
  jac[1][0] = nodes[ nodeNum[1] ]->y - nodes[ nodeNum[0] ]->y;
  jac[1][1] = nodes[ nodeNum[2] ]->y - nodes[ nodeNum[0] ]->y;
  jac[1][2] = nodes[ nodeNum[3] ]->y - nodes[ nodeNum[0] ]->y;
  jac[2][0] = nodes[ nodeNum[1] ]->z - nodes[ nodeNum[0] ]->z;
  jac[2][1] = nodes[ nodeNum[2] ]->z - nodes[ nodeNum[0] ]->z;
  jac[2][2] = nodes[ nodeNum[3] ]->z - nodes[ nodeNum[0] ]->z;

  // compute determinant of jac
  double dOmega = jac[0][0] * (jac[1][1] * jac[2][2]
                - jac[1][2] * jac[2][1])
                + jac[1][0] * (jac[0][2] * jac[2][1]
                - jac[0][1] * jac[2][2])
                + jac[2][0] * (jac[0][1] * jac[1][2]
                - jac[0][2] * jac[1][1]);

  // compute inverse matrix of jac
  // Maple code used
  double t17 = -1.0/dOmega;

  //compute shape function gradients
  nGrad[1][0] =  (-jac[1][1] * jac[2][2] + jac[1][2] * jac[2][1] ) * t17;
  nGrad[1][1] =  ( jac[0][1] * jac[2][2] - jac[0][2] * jac[2][1] ) * t17;
  nGrad[1][2] = -( jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1] ) * t17;

  // change by tl
  nGrad[2][0] = -(-jac[1][0] * jac[2][2] + jac[1][2] * jac[2][0] ) * t17;
  //nGrad[2][0] =  (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0] ) * t17;

  nGrad[2][1] = -( jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0] ) * t17;
  nGrad[2][2] =  ( jac[0][0] * jac[1][2] - jac[0][2] * jac[1][0] ) * t17;
  nGrad[3][0] = -( jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0] ) * t17;
  nGrad[3][1] =  ( jac[0][0] * jac[2][1] - jac[0][1] * jac[2][0] ) * t17;
  nGrad[3][2] = -( jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0] ) * t17;

  // Shape function gradients dN_i/dx_i = dN/dxi * transpose(jInv)
  // Note: 1st index = shape function #
  // 2nd index = direction (0=x, 1=y, 2=z)

  nGrad[0][0] = -( nGrad[1][0] + nGrad[2][0] + nGrad[3][0] );
  nGrad[0][1] = -( nGrad[1][1] + nGrad[2][1] + nGrad[3][1] );
  nGrad[0][2] = -( nGrad[1][2] + nGrad[2][2] + nGrad[3][2] );

  return dOmega;
}

void
TetCorotator::extractDeformations(GeomState &geomState, CoordSet &cs, 
                                  double *vld, int &nlflag)
{
 // Set Flag to Use Non-Linear Routines for Stress
 nlflag = 2;
}

void
TetCorotator::extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                     double *vlr)
{
}

double
TetCorotator::getElementEnergy(GeomState &geomState, CoordSet &cs)
{
// Computes Internal Energy of Element in Given State

  int j;
  double nGrad[4][3];

  // get the nodal temperatures
  Vector ndTemps(4);
  geomState.get_temperature(4, nodeNum, ndTemps, Tref);
 
  // compute dN_i/dX_j, also obtain dOmega
  double dOmega;
  dOmega = computeShapeGrad(cs, nGrad)/6;

  // now get F_ij = dPhi_i/dX_j = x^k_i dN_k/dX_j
  double F[3][3];

  for(j = 0; j < 3; ++j) {
    F[0][j] = geomState[nodeNum[0]].x*nGrad[0][j]
            + geomState[nodeNum[1]].x*nGrad[1][j]
            + geomState[nodeNum[2]].x*nGrad[2][j]
            + geomState[nodeNum[3]].x*nGrad[3][j];

    F[1][j] = geomState[nodeNum[0]].y*nGrad[0][j]
            + geomState[nodeNum[1]].y*nGrad[1][j]
            + geomState[nodeNum[2]].y*nGrad[2][j]
            + geomState[nodeNum[3]].y*nGrad[3][j];

    F[2][j] = geomState[nodeNum[0]].z*nGrad[0][j]
            + geomState[nodeNum[1]].z*nGrad[1][j]
            + geomState[nodeNum[2]].z*nGrad[2][j]
            + geomState[nodeNum[3]].z*nGrad[3][j];
  }

  // compute e_ij = 0.5*(F_ki Fkj - delta_ij)
  double e_11 = 0.5*(F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][0]-1.0);
  double e_22 = 0.5*(F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1]-1.0);
  double e_33 = 0.5*(F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2]-1.0);
  // here engineering strain is computed
  double e_12 = (F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1]);
  double e_13 = (F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2]);
  double e_23 = (F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2]);

  // Subtract thermal strain
  double theta = (ndTemps[0]+ndTemps[1]+ndTemps[2]+ndTemps[3])/4 - Tref;
  double alpha = (ctett) ? ctett->getValAlt(theta) : TetCorotator::alpha;
  e_11 -= alpha*theta;
  e_22 -= alpha*theta;
  e_33 -= alpha*theta;

  double em = (ymtt) ? ymtt->getValAlt(theta) : TetCorotator::em;
  double E2 = em*nu/((1+nu)*(1-2*nu));
  double G2 = em/(2*(1+nu));
  double E1 = E2+em/(1+nu);
  double s_11 = E1*e_11+E2*(e_22+e_33);
  double s_22 = E1*e_22+E2*(e_11+e_33);
  double s_33 = E1*e_33+E2*(e_11+e_22);
  double s_12 = G2*e_12;
  double s_13 = G2*e_13;
  double s_23 = G2*e_23;

  double Energy = dOmega*((e_11*s_11) +
                          (e_22*s_22) +
                          (e_33*s_33) +
                          (e_12*s_12) +
                          (e_13*s_13) +
                          (e_23*s_23));

  Energy *= 0.5;
  return Energy;
}
