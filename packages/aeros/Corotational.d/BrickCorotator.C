#include <cmath>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Utils.d/MFTT.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Element.d/Element.h>
#include <Corotational.d/BrickCorotator.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>

extern "C" {
  void _FORTRAN(hxgaus)(int &, int &, int &, int &, int &, int &,
                        double &, double &, double &, double &);
  void _FORTRAN(h8shpe)(double &, double &, double &, double *, double *,
                        double *, double *, double *, double *, double *,
                        double &);
  void _FORTRAN(vmelmv)(double*, int &, int &, int &, int &, int &);
  void _FORTRAN(strainvm)(double*, int &, int &, int &, int &);
};

BrickCorotator::BrickCorotator(int nodeNumbers[8], double _em, double _nu,
                               CoordSet& cs, double _Tref, double _alpha,
                               MFTTData *_ymtt, MFTTData *_ctett)
{
  nodeNum[0] = nodeNumbers[0];
  nodeNum[1] = nodeNumbers[1];
  nodeNum[2] = nodeNumbers[2];
  nodeNum[3] = nodeNumbers[3];
  nodeNum[4] = nodeNumbers[4];
  nodeNum[5] = nodeNumbers[5];
  nodeNum[6] = nodeNumbers[6];
  nodeNum[7] = nodeNumbers[7];

  em = _em;       // Elastic modulus
  nu = _nu;       // Poisson's ratio
  Tref = _Tref;   // Ambient temperature
  alpha = _alpha; // Thermal expansion coefficient
  ymtt = _ymtt; 
  ctett = _ctett;
}

// geomState -> contains the updated nodal coordinates
// cs        -> contains the original nodal coordinates
void
BrickCorotator::getStiffAndForce(GeomState &geomState, CoordSet &cs,
                                 FullSquareMatrix &K, double *f, double dt, double t)
{
  int i,j,k;
  int numLinGaussPts = 2;
  double nGrad[8][3];

  // initialize forces and stiffness matrix
  for (i = 0; i < 24; i++) {
    f[i] = 0;
    for (j = 0; j < 24; j++)
      K[i][j] = 0;
  }

  // reformat cs to accommodate fortran routine
  double xNodes[8], yNodes[8], zNodes[8];
  for (i = 0; i < 8; i++) {
    xNodes[i] = cs[nodeNum[i]]->x;
    yNodes[i] = cs[nodeNum[i]]->y;
    zNodes[i] = cs[nodeNum[i]]->z;
  }

  // get the nodal temperatures
  Vector ndTemps(8);
  geomState.get_temperature(8, nodeNum, ndTemps, Tref);

  int fortran = 1; // fortran routines start from index 1
  int pt1, pt2, pt3;
  for (pt1 = 0 + fortran; pt1 < numLinGaussPts + fortran; pt1++) {
    for (pt2 = 0 + fortran; pt2 < numLinGaussPts + fortran; pt2++) {
      for (pt3 = 0 + fortran; pt3 < numLinGaussPts + fortran; pt3++) {
        // get gauss point
        double xi, eta, mu, wt;
        _FORTRAN(hxgaus)(numLinGaussPts, pt1, numLinGaussPts, pt2,
                         numLinGaussPts, pt3, xi, eta, mu, wt);

        // compute shape functions
        double shapeFunc[8], shapeGradX[8], shapeGradY[8], shapeGradZ[8];
        double dOmega; // det of jacobian

        _FORTRAN(h8shpe)(xi, eta, mu, xNodes, yNodes, zNodes,
                         shapeFunc, shapeGradX, shapeGradY, shapeGradZ, dOmega);

        // get volume
        // dOmega is here off by a factor of 1/4
        // this divides all by a factor of 4 later
        dOmega *= wt;
        dOmega /= 4;

        for (i = 0; i < 8; ++i)  {
          nGrad[i][0] = shapeGradX[i];
          nGrad[i][1] = shapeGradY[i];
          nGrad[i][2] = shapeGradZ[i];
        }

        // now get F_ij = dPhi_i/dX_j = x^k_i dN_k/dX_j
        double F[3][3];

        for(j = 0; j < 3; ++j)
          F[0][j] = geomState[nodeNum[0]].x * nGrad[0][j]
                  + geomState[nodeNum[1]].x * nGrad[1][j]
                  + geomState[nodeNum[2]].x * nGrad[2][j]
                  + geomState[nodeNum[3]].x * nGrad[3][j]
                  + geomState[nodeNum[4]].x * nGrad[4][j]
                  + geomState[nodeNum[5]].x * nGrad[5][j]
                  + geomState[nodeNum[6]].x * nGrad[6][j]
                  + geomState[nodeNum[7]].x * nGrad[7][j];

        for(j = 0; j < 3; ++j)
          F[1][j] = geomState[nodeNum[0]].y * nGrad[0][j]
                  + geomState[nodeNum[1]].y * nGrad[1][j]
                  + geomState[nodeNum[2]].y * nGrad[2][j]
                  + geomState[nodeNum[3]].y * nGrad[3][j]
                  + geomState[nodeNum[4]].y * nGrad[4][j]
                  + geomState[nodeNum[5]].y * nGrad[5][j]
                  + geomState[nodeNum[6]].y * nGrad[6][j]
                  + geomState[nodeNum[7]].y * nGrad[7][j];

        for(j = 0; j < 3; ++j)
          F[2][j] = geomState[nodeNum[0]].z * nGrad[0][j]
                  + geomState[nodeNum[1]].z * nGrad[1][j]
                  + geomState[nodeNum[2]].z * nGrad[2][j]
                  + geomState[nodeNum[3]].z * nGrad[3][j]
                  + geomState[nodeNum[4]].z * nGrad[4][j]
                  + geomState[nodeNum[5]].z * nGrad[5][j]
                  + geomState[nodeNum[6]].z * nGrad[6][j]
                  + geomState[nodeNum[7]].z * nGrad[7][j];

        //double detF = F[0][0]*F[1][1]*F[2][2]+F[0][1]*F[1][2]*F[2][0]+F[0][2]*F[1][0]*F[2][1]
        //             -F[0][0]*F[1][2]*F[2][1]-F[0][1]*F[1][0]*F[2][2]-F[0][2]*F[1][1]*F[2][0];
        //if(detF < 0) cerr << " *** WARNING: in BrickCorotator::getStiffAndForce |F| = " << detF << endl;

        // compute e_ij = 0.5*(F_ki Fkj - delta_ij)
        // here these are off by factor of 2
        double e_11 = (F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][0]-1.0);
        double e_22 = (F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1]-1.0);
        double e_33 = (F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2]-1.0);
        double e_12 = (F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1]);
        double e_13 = (F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2]);
        double e_23 = (F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2]);

        // Subtract thermal strain (off by factor of 2)
        double theta = 0.0;
        for(j = 0; j < 8; j++) theta += shapeFunc[j]*(ndTemps[j] - Tref);
        double alpha = (ctett) ? ctett->getValAlt(theta) : BrickCorotator::alpha;
        e_11 -= 2*alpha*theta;
        e_22 -= 2*alpha*theta;
        e_33 -= 2*alpha*theta;

        double sigma[6];
        double em = (ymtt) ? ymtt->getValAlt(theta) : BrickCorotator::em;
        double E2 = em*nu/((1+nu)*(1-2*nu));
        double E1 = E2+em/(1+nu);
        // no factor of 1/2 on G2 due to using tensor strain
        double G2 = em/(1+nu);
        // these here are off by a factor of 2
        sigma[0] = E1*e_11+E2*(e_22+e_33);
        sigma[1] = E1*e_22+E2*(e_11+e_33);
        sigma[2] = E1*e_33+E2*(e_11+e_22);
        // these here are off by a factor of 4
        sigma[3] = 2*G2*e_12;
        sigma[4] = 2*G2*e_13;
        sigma[5] = 2*G2*e_23;

        // Compute de_ij/dUl for the symmetric part
        // First we get df_ij/dUl in a very compact form.
        // df_ij/dUl=dN_p/dX_j delta_iq; with p = int(l/3)+1 and l-1=q-1 mod(3)
        // this means that df_ij/dUl is already contained in dN_k/dX_j
        double dedU[24][6];
        for(i = 0; i < 8; ++i)
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
        for(i = 0; i < 24; ++i)
          f[i] += dOmega*( dedU[i][0]*sigma[0] +
                           dedU[i][1]*sigma[1] +
                           dedU[i][2]*sigma[2] +
                           dedU[i][3]*sigma[3] +
                           dedU[i][4]*sigma[4] +
                           dedU[i][5]*sigma[5]);
        // normal terms of f are
        // (1/4)*(2*2) = 1 => right level
        // shear terms of f are
        // (1/4)*(2*4)*(1/2 for off-diag mult) = 1 => right level

        // now get ds_ij/dUl
        double dsdU[24][6];
        for(i = 0; i < 24; ++i) {
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
        for(i = 0; i < 24; ++i)
          for(j = 0; j <= i; ++j)
            K[j][i] += dOmega*(dsdU[i][0]*dedU[j][0] +
                               dsdU[i][1]*dedU[j][1] +
                               dsdU[i][2]*dedU[j][2] +
                               dsdU[i][3]*dedU[j][3] +
                               dsdU[i][4]*dedU[j][4] +
                               dsdU[i][5]*dedU[j][5]);
        // normal terms are (1/4)*2*2 = 1 => right level
        // shear terms are  (1/4)*4*2*(1/2 for off-diag mult) = 1 => right level

        // add s*d2e/dU_idUj (symmetric part only)
        for(i = 0; i < 8; ++i)
          for(k = 0; k <= i; ++k)
            for(j = 0; j < 3; ++j) {
              K[3*k+j][3*i+j] +=
              dOmega*(sigma[0]*(2*nGrad[i][0]*nGrad[k][0]) +
                      sigma[1]*(2*nGrad[i][1]*nGrad[k][1]) +
                      sigma[2]*(2*nGrad[i][2]*nGrad[k][2]) +
                      sigma[3]*(nGrad[i][0]*nGrad[k][1]+nGrad[i][1]*nGrad[k][0]) +
                      sigma[4]*(nGrad[i][0]*nGrad[k][2]+nGrad[i][2]*nGrad[k][0]) +
                      sigma[5]*(nGrad[i][1]*nGrad[k][2]+nGrad[i][2]*nGrad[k][1]));
            }
        // normal terms are (1/4)*2*2 = 1 => right level
        // shear terms are  (1/4)*4*2*(1/2 for off-diag mult) = 1 => right level
      }
    }
  }

  // Symmetrize
  for(i = 0; i < 24; ++i)
    for(j = 0; j <= i; ++j)
      K[i][j] = K[j][i];
}

void
BrickCorotator::getInternalForce(GeomState &geomState, CoordSet &cs,
                                 FullSquareMatrix &, double *f, double dt, double t)
{
  int i,j,k;
  int numLinGaussPts = 2;
  double nGrad[8][3];

  // initialize forces
  for (i = 0; i < 24; i++) {
    f[i] = 0;
  }

  // reformat cs to accommodate fortran routine
  double xNodes[8], yNodes[8], zNodes[8];
  for (i = 0; i < 8; i++) {
    xNodes[i] = cs[nodeNum[i]]->x;
    yNodes[i] = cs[nodeNum[i]]->y;
    zNodes[i] = cs[nodeNum[i]]->z;
  }

  // get the nodal temperatures
  Vector ndTemps(8);
  geomState.get_temperature(8, nodeNum, ndTemps, Tref);

  int fortran = 1; // fortran routines start from index 1
  int pt1, pt2, pt3;
  for (pt1 = 0 + fortran; pt1 < numLinGaussPts + fortran; pt1++) {
    for (pt2 = 0 + fortran; pt2 < numLinGaussPts + fortran; pt2++) {
      for (pt3 = 0 + fortran; pt3 < numLinGaussPts + fortran; pt3++) {
        // get gauss point
        double xi, eta, mu, wt;
        _FORTRAN(hxgaus)(numLinGaussPts, pt1, numLinGaussPts, pt2,
                         numLinGaussPts, pt3, xi, eta, mu, wt);

        // compute shape functions
        double shapeFunc[8], shapeGradX[8], shapeGradY[8], shapeGradZ[8];
        double dOmega; // det of jacobian

        _FORTRAN(h8shpe)(xi, eta, mu, xNodes, yNodes, zNodes,
                         shapeFunc, shapeGradX, shapeGradY, shapeGradZ, dOmega);

        // get volume
        // dOmega is here off by a factor of 1/4
        // this divides all by a factor of 4 later
        dOmega *= wt;
        dOmega /= 4;

        for (i = 0; i < 8; ++i) {
          nGrad[i][0] = shapeGradX[i];
          nGrad[i][1] = shapeGradY[i];
          nGrad[i][2] = shapeGradZ[i];
        }

        // now get F_ij = dPhi_i/dX_j = x^k_i dN_k/dX_j
        double F[3][3];

        for(j = 0; j < 3; ++j)
          F[0][j] = geomState[nodeNum[0]].x * nGrad[0][j]
                  + geomState[nodeNum[1]].x * nGrad[1][j]
                  + geomState[nodeNum[2]].x * nGrad[2][j]
                  + geomState[nodeNum[3]].x * nGrad[3][j]
                  + geomState[nodeNum[4]].x * nGrad[4][j]
                  + geomState[nodeNum[5]].x * nGrad[5][j]
                  + geomState[nodeNum[6]].x * nGrad[6][j]
                  + geomState[nodeNum[7]].x * nGrad[7][j];

        for(j = 0; j < 3; ++j)
          F[1][j] = geomState[nodeNum[0]].y * nGrad[0][j]
                  + geomState[nodeNum[1]].y * nGrad[1][j]
                  + geomState[nodeNum[2]].y * nGrad[2][j]
                  + geomState[nodeNum[3]].y * nGrad[3][j]
                  + geomState[nodeNum[4]].y * nGrad[4][j]
                  + geomState[nodeNum[5]].y * nGrad[5][j]
                  + geomState[nodeNum[6]].y * nGrad[6][j]
                  + geomState[nodeNum[7]].y * nGrad[7][j];

        for(j = 0; j < 3; ++j)
          F[2][j] = geomState[nodeNum[0]].z * nGrad[0][j]
                  + geomState[nodeNum[1]].z * nGrad[1][j]
                  + geomState[nodeNum[2]].z * nGrad[2][j]
                  + geomState[nodeNum[3]].z * nGrad[3][j]
                  + geomState[nodeNum[4]].z * nGrad[4][j]
                  + geomState[nodeNum[5]].z * nGrad[5][j]
                  + geomState[nodeNum[6]].z * nGrad[6][j]
                  + geomState[nodeNum[7]].z * nGrad[7][j];

        //double detF = F[0][0]*F[1][1]*F[2][2]+F[0][1]*F[1][2]*F[2][0]+F[0][2]*F[1][0]*F[2][1]
        //             -F[0][0]*F[1][2]*F[2][1]-F[0][1]*F[1][0]*F[2][2]-F[0][2]*F[1][1]*F[2][0];
        //if(detF < 0) cerr << " *** WARNING: in BrickCorotator::getInternalForce |F| = " << detF << endl;

        // compute e_ij = 0.5*(F_ki Fkj - delta_ij)
        // here these are off by factor of 2
        double e_11 = (F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][0]-1.0);
        double e_22 = (F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1]-1.0);
        double e_33 = (F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2]-1.0);
        double e_12 = (F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1]);
        double e_13 = (F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2]);
        double e_23 = (F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2]);

        // Subtract thermal strain (off by factor of 2)
        double theta = 0.0;
        for(j = 0; j < 8; j++) theta += shapeFunc[j]*(ndTemps[j] - Tref);
        double alpha = (ctett) ? ctett->getValAlt(theta) : BrickCorotator::alpha;
        e_11 -= 2*alpha*theta;
        e_22 -= 2*alpha*theta;
        e_33 -= 2*alpha*theta;

        double sigma[6];
        double em = (ymtt) ? ymtt->getValAlt(theta) : BrickCorotator::em;
        double E2 = em*nu/((1+nu)*(1-2*nu));
        double E1 = E2+em/(1+nu);
        // no factor of 1/2 on G2 due to using tensor strain
        double G2 = em/(1+nu);
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
        double dedU[24][6];
        for(i = 0; i < 8; ++i)
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
        for(i = 0; i < 24; ++i)
          f[i] += dOmega*( dedU[i][0]*sigma[0] +
                           dedU[i][1]*sigma[1] +
                           dedU[i][2]*sigma[2] +
                           dedU[i][3]*sigma[3] +
                           dedU[i][4]*sigma[4] +
                           dedU[i][5]*sigma[5]);
        // normal terms of f are
        // (1/4)*(2*2) = 1 => right level
        // shear terms of f are
        // (1/4)*(2*4)*(1/2 for off-diag mult) = 1 => right level
      }
    }
  }

}

//-------------------------------------------------------------------------------

double
BrickCorotator::computeStrainGrad(GeomState &geomState, CoordSet &cs,
                                  double dedU[24][6], int pt1, int pt2, int pt3)
{
  int i, j;
  double dOmega;
  int numLinGaussPts = 2;
  double nGrad[8][3];

  // reformat cs to accommodate fortran routine
  double xNodes[8], yNodes[8], zNodes[8];
  for (i = 0; i < 8; i++) {
    xNodes[i] = cs[nodeNum[i]]->x;
    yNodes[i] = cs[nodeNum[i]]->y;
    zNodes[i] = cs[nodeNum[i]]->z;
  }

  // get gauss point
  double xi, eta, mu, wt;
  _FORTRAN(hxgaus)(numLinGaussPts, pt1, numLinGaussPts, pt2,
                   numLinGaussPts, pt3, xi, eta, mu, wt);

  // compute shape functions
  double shapeFunc[8], shapeGradX[8], shapeGradY[8], shapeGradZ[8];

  _FORTRAN(h8shpe)(xi, eta, mu, xNodes, yNodes, zNodes,
                   shapeFunc, shapeGradX, shapeGradY, shapeGradZ, dOmega);

  for (i = 0; i < 8; ++i) {
    nGrad[i][0] = shapeGradX[i];
    nGrad[i][1] = shapeGradY[i];
    nGrad[i][2] = shapeGradZ[i];
  }

  // now get F_ij = dPhi_i/dX_j = x^k_i dN_k/dX_j
  double F[3][3];

  for (j = 0; j < 3; ++j)
    F[0][j] = geomState[nodeNum[0]].x * nGrad[0][j]
            + geomState[nodeNum[1]].x * nGrad[1][j]
            + geomState[nodeNum[2]].x * nGrad[2][j]
            + geomState[nodeNum[3]].x * nGrad[3][j]
            + geomState[nodeNum[4]].x * nGrad[4][j]
            + geomState[nodeNum[5]].x * nGrad[5][j]
            + geomState[nodeNum[6]].x * nGrad[6][j]
            + geomState[nodeNum[7]].x * nGrad[7][j];

  for (j = 0; j < 3; ++j)
    F[1][j] = geomState[nodeNum[0]].y * nGrad[0][j]
            + geomState[nodeNum[1]].y * nGrad[1][j]
            + geomState[nodeNum[2]].y * nGrad[2][j]
            + geomState[nodeNum[3]].y * nGrad[3][j]
            + geomState[nodeNum[4]].y * nGrad[4][j]
            + geomState[nodeNum[5]].y * nGrad[5][j]
            + geomState[nodeNum[6]].y * nGrad[6][j]
            + geomState[nodeNum[7]].y * nGrad[7][j];

  for (j = 0; j < 3; ++j)
    F[2][j] = geomState[nodeNum[0]].z * nGrad[0][j]
            + geomState[nodeNum[1]].z * nGrad[1][j]
            + geomState[nodeNum[2]].z * nGrad[2][j]
            + geomState[nodeNum[3]].z * nGrad[3][j]
            + geomState[nodeNum[4]].z * nGrad[4][j]
            + geomState[nodeNum[5]].z * nGrad[5][j]
            + geomState[nodeNum[6]].z * nGrad[6][j]
            + geomState[nodeNum[7]].z * nGrad[7][j];

  for (i = 0; i < 8; ++i)
    for (j = 0; j < 3; ++j) {
      dedU[3*i+j][0] = nGrad[i][0]*F[j][0];
      dedU[3*i+j][1] = nGrad[i][1]*F[j][1];
      dedU[3*i+j][2] = nGrad[i][2]*F[j][2];
      dedU[3*i+j][3] = 0.5*(nGrad[i][0]*F[j][1]+nGrad[i][1]*F[j][0]);
      dedU[3*i+j][4] = 0.5*(nGrad[i][0]*F[j][2]+nGrad[i][2]*F[j][0]);
      dedU[3*i+j][5] = 0.5*(nGrad[i][1]*F[j][2]+nGrad[i][2]*F[j][1]);
    }

  return dOmega;
}

//-------------------------------------------------------------------------------

void
BrickCorotator::getNLVonMises(Vector& stress, Vector& weight, GeomState &geomState,
                              GeomState *, CoordSet& cs, int strInd, int,
                              double, double, int, int)
{
  weight = 1.0;

  int maxstr = 7;
  int elm    = 1;
  int nno    = 8;

  double elStress[8][7];
  double elStrain[8][7];

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
    stress[4] = elStress[4][strInd];
    stress[5] = elStress[5][strInd];
    stress[6] = elStress[6][strInd];
    stress[7] = elStress[7][strInd];
  } else if(strInd < 14) {
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
BrickCorotator::getNLAllStress(FullM &stress, Vector &weight, GeomState &geomState,
                               GeomState *, CoordSet &cs, int strInd, int, int)
{
  weight = 1.0;

  int i,j;

  double elStress[8][7];
  double elStrain[8][7];

  // Compute NL Stress/Strain
  computePiolaStress(geomState, cs, elStress, elStrain);

  // Store all Stress or all Strain as defined by strInd
  if(strInd == 0) {
    for (i=0; i<8; ++i) {
      for (j=0; j<6; ++j) {
        stress[i][j] = elStress[i][j];
      }
    }
  } else {
    for (i=0; i<8; ++i) {
      for (j=0; j<6; ++j) {
        stress[i][j] = elStrain[i][j];
      }
    }
  }

  // Get Element Principals without averaging
  double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double pvec[3] = {0.0,0.0,0.0};
  for (i=0; i<8; ++i) {
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
BrickCorotator::computePiolaStress(GeomState &geomState, CoordSet &cs,
                                   double stress[8][7], double strain[8][7])
{
  int i,j,n;
  double nGrad[8][3];

  double xinod[8] = {-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0};
  double etanod[8] = {-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0};
  double munod[8] = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0};
  double xi,eta,mu;

  // reformat cs to accommodate fortran routine
  double xNodes[8], yNodes[8], zNodes[8];
  for (i = 0; i < 8; i++) {
    xNodes[i] = cs[nodeNum[i]]->x;
    yNodes[i] = cs[nodeNum[i]]->y;
    zNodes[i] = cs[nodeNum[i]]->z;
  }

  // get the nodal temperatures
  Vector ndTemps(8);
  geomState.get_temperature(8, nodeNum, ndTemps, Tref);

  for (n = 0; n < 8; n++) {
    xi  = xinod[n];
    eta = etanod[n];
    mu  = munod[n];

    // compute shape functions
    double shapeFunc[8], shapeGradX[8], shapeGradY[8], shapeGradZ[8];
    double dOmega; // det of jacobian

    _FORTRAN(h8shpe)(xi, eta, mu, xNodes, yNodes, zNodes,
                     shapeFunc, shapeGradX, shapeGradY, shapeGradZ, dOmega);

    for (i = 0; i < 8; ++i)  {
      nGrad[i][0] = shapeGradX[i];
      nGrad[i][1] = shapeGradY[i];
      nGrad[i][2] = shapeGradZ[i];
    }

    // now get F_ij = dPhi_i/dX_j = x^k_i dN_k/dX_j
    double F[3][3];

    for(j = 0; j < 3; ++j)
      F[0][j] = geomState[nodeNum[0]].x * nGrad[0][j]
              + geomState[nodeNum[1]].x * nGrad[1][j]
              + geomState[nodeNum[2]].x * nGrad[2][j]
              + geomState[nodeNum[3]].x * nGrad[3][j]
              + geomState[nodeNum[4]].x * nGrad[4][j]
              + geomState[nodeNum[5]].x * nGrad[5][j]
              + geomState[nodeNum[6]].x * nGrad[6][j]
              + geomState[nodeNum[7]].x * nGrad[7][j];

    for(j = 0; j < 3; ++j)
      F[1][j] = geomState[nodeNum[0]].y * nGrad[0][j]
              + geomState[nodeNum[1]].y * nGrad[1][j]
              + geomState[nodeNum[2]].y * nGrad[2][j]
              + geomState[nodeNum[3]].y * nGrad[3][j]
              + geomState[nodeNum[4]].y * nGrad[4][j]
              + geomState[nodeNum[5]].y * nGrad[5][j]
              + geomState[nodeNum[6]].y * nGrad[6][j]
              + geomState[nodeNum[7]].y * nGrad[7][j];

    for(j = 0; j < 3; ++j)
      F[2][j] = geomState[nodeNum[0]].z * nGrad[0][j]
              + geomState[nodeNum[1]].z * nGrad[1][j]
              + geomState[nodeNum[2]].z * nGrad[2][j]
              + geomState[nodeNum[3]].z * nGrad[3][j]
              + geomState[nodeNum[4]].z * nGrad[4][j]
              + geomState[nodeNum[5]].z * nGrad[5][j]
              + geomState[nodeNum[6]].z * nGrad[6][j]
              + geomState[nodeNum[7]].z * nGrad[7][j];

    // compute e_ij = 0.5*(F_ki Fkj - delta_ij)
    double e_11 = 0.5*(F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][0]-1.0);
    double e_22 = 0.5*(F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1]-1.0);
    double e_33 = 0.5*(F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2]-1.0);
    // here engineering strain is computed
    double e_12 = (F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1]);
    double e_13 = (F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2]);
    double e_23 = (F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2]);

    // Reorder strain
    strain[n][0] = e_11;
    strain[n][1] = e_22;
    strain[n][2] = e_33;
    strain[n][3] = e_12;
    strain[n][4] = e_23;
    strain[n][5] = e_13;

    // Subtract thermal strain
    double alpha = (ctett) ? ctett->getValAlt(ndTemps[n]) : BrickCorotator::alpha;
    e_11 -= alpha*(ndTemps[n]-Tref);
    e_22 -= alpha*(ndTemps[n]-Tref);
    e_33 -= alpha*(ndTemps[n]-Tref);

    double sigma[6];
    double em = (ymtt) ? ymtt->getValAlt(ndTemps[n]) : BrickCorotator::em;
    double E2 = em*nu/((1+nu)*(1-2*nu));
    double G2 = em/(2*(1+nu));
    double E1 = E2+em/(1+nu);
    sigma[0] = E1*e_11+E2*(e_22+e_33);
    sigma[1] = E1*e_22+E2*(e_11+e_33);
    sigma[2] = E1*e_33+E2*(e_11+e_22);
    sigma[3] = G2*e_12;
    sigma[4] = G2*e_23;
    sigma[5] = G2*e_13;

    // Reorder stress
    stress[n][0] = sigma[0];
    stress[n][1] = sigma[1];
    stress[n][2] = sigma[2];
    stress[n][3] = sigma[3];
    stress[n][4] = sigma[4];
    stress[n][5] = sigma[5];
  }
}

void
BrickCorotator::extractDeformations(GeomState &geomState, CoordSet &cs,
                                    double *vld, int &nlflag)
{
  // Set Flag to Use Non-Linear Routines for Stress
  nlflag = 2;
}

void
BrickCorotator::extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                       double *vlr)
{
}

double
BrickCorotator::getElementEnergy(GeomState &geomState, CoordSet &cs)
{
  // Computes Internal Energy of Element in Given State
  int i,j;
  int numLinGaussPts = 2;
  double nGrad[8][3];

  // initialize energy to zero
  double Energy = 0;

  // reformat cs to accommodate fortran routines
  double xNodes[8], yNodes[8], zNodes[8];
  for (i = 0; i < 8; i++) {
    xNodes[i] = cs[nodeNum[i]]->x;
    yNodes[i] = cs[nodeNum[i]]->y;
    zNodes[i] = cs[nodeNum[i]]->z;
  }

  // get the nodal temperatures
  Vector ndTemps(8);
  geomState.get_temperature(8, nodeNum, ndTemps, Tref);

  int fortran = 1; // fortran routines start from index 1
  int pt1, pt2, pt3;
  for (pt1 = 0 + fortran; pt1 < numLinGaussPts + fortran; pt1++) {
    for (pt2 = 0 + fortran; pt2 < numLinGaussPts + fortran; pt2++) {
      for (pt3 = 0 + fortran; pt3 < numLinGaussPts + fortran; pt3++) {
        // get gauss point
        double xi, eta, mu, wt;
        _FORTRAN(hxgaus)(numLinGaussPts, pt1, numLinGaussPts, pt2,
                         numLinGaussPts, pt3, xi, eta, mu, wt);

        // compute shape functions
        double shapeFunc[8], shapeGradX[8], shapeGradY[8], shapeGradZ[8];
        double dOmega; // det of jacobian

        _FORTRAN(h8shpe)(xi, eta, mu, xNodes, yNodes, zNodes,
                       shapeFunc, shapeGradX, shapeGradY, shapeGradZ, dOmega);

        // get volume
        dOmega *= wt;

        for (i = 0; i < 8; ++i)  {
          nGrad[i][0] = shapeGradX[i];
          nGrad[i][1] = shapeGradY[i];
          nGrad[i][2] = shapeGradZ[i];
        }

        // now get F_ij = dPhi_i/dX_j = x^k_i dN_k/dX_j
        double F[3][3];

        for(j = 0; j < 3; ++j)
          F[0][j] = geomState[nodeNum[0]].x * nGrad[0][j]
                  + geomState[nodeNum[1]].x * nGrad[1][j]
                  + geomState[nodeNum[2]].x * nGrad[2][j]
                  + geomState[nodeNum[3]].x * nGrad[3][j]
                  + geomState[nodeNum[4]].x * nGrad[4][j]
                  + geomState[nodeNum[5]].x * nGrad[5][j]
                  + geomState[nodeNum[6]].x * nGrad[6][j]
                  + geomState[nodeNum[7]].x * nGrad[7][j];

        for(j = 0; j < 3; ++j)
          F[1][j] = geomState[nodeNum[0]].y * nGrad[0][j]
                  + geomState[nodeNum[1]].y * nGrad[1][j]
                  + geomState[nodeNum[2]].y * nGrad[2][j]
                  + geomState[nodeNum[3]].y * nGrad[3][j]
                  + geomState[nodeNum[4]].y * nGrad[4][j]
                  + geomState[nodeNum[5]].y * nGrad[5][j]
                  + geomState[nodeNum[6]].y * nGrad[6][j]
                  + geomState[nodeNum[7]].y * nGrad[7][j];

        for(j = 0; j < 3; ++j)
          F[2][j] = geomState[nodeNum[0]].z * nGrad[0][j]
                  + geomState[nodeNum[1]].z * nGrad[1][j]
                  + geomState[nodeNum[2]].z * nGrad[2][j]
                  + geomState[nodeNum[3]].z * nGrad[3][j]
                  + geomState[nodeNum[4]].z * nGrad[4][j]
                  + geomState[nodeNum[5]].z * nGrad[5][j]
                  + geomState[nodeNum[6]].z * nGrad[6][j]
                  + geomState[nodeNum[7]].z * nGrad[7][j];

        // compute e_ij = 0.5*(F_ki Fkj - delta_ij)
        double e_11 = 0.5*(F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][0]-1.0);
        double e_22 = 0.5*(F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1]-1.0);
        double e_33 = 0.5*(F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2]-1.0);
        // here engineering strain is computed
        double e_12 = (F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1]);
        double e_13 = (F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2]);
        double e_23 = (F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2]);

        // Subtract thermal strain
        double theta = 0.0;
        for(j = 0; j < 8; j++) theta += shapeFunc[j]*(ndTemps[j] - Tref);
        double alpha = (ctett) ? ctett->getValAlt(theta) : BrickCorotator::alpha;
        e_11 -= alpha*theta;
        e_22 -= alpha*theta;
        e_33 -= alpha*theta;

        double em = (ymtt) ? ymtt->getValAlt(theta) : BrickCorotator::em;
        double E2 = em*nu/((1+nu)*(1-2*nu));
        double G2 = em/(2*(1+nu));
        double E1 = E2+em/(1+nu);
        double s_11 = E1*e_11+E2*(e_22+e_33);
        double s_22 = E1*e_22+E2*(e_11+e_33);
        double s_33 = E1*e_33+E2*(e_11+e_22);
        double s_12 = G2*e_12;
        double s_13 = G2*e_13;
        double s_23 = G2*e_23;

        Energy += dOmega*((e_11*s_11) +
                          (e_22*s_22) +
                          (e_33*s_33) +
                          (e_12*s_12) +
                          (e_13*s_13) +
                          (e_23*s_23));

      }
    }
  }
  Energy *= 0.5;
  return Energy;
}
