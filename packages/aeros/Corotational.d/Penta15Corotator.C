#include <cmath>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Utils.d/MFTT.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Element.d/Element.h>
#include <Corotational.d/Penta15Corotator.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>

extern "C" {
  void _FORTRAN(vmelmv)(double*, int &, int &, int &, int &, int &);
  void _FORTRAN(strainvm)(double*, int &, int &, int &, int &);
};

extern double Penta15ShapeFct(double Shape[15], double DShape[15][3], double m[3], double X[15], double Y[15], double Z[15]);

extern double weight3d8[9];
extern double gauss3d8[9][3]; 

Penta15Corotator::Penta15Corotator(int nodeNumbers[15], double _em, double _nu, CoordSet& cs, double _Tref, double _alpha,
                                   MFTTData *_ymtt, MFTTData *_ctett)
{
  for(int i = 0; i < 15; ++i)
    nodeNum[i] = nodeNumbers[i];

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
Penta15Corotator::getStiffAndForce(GeomState &geomState, CoordSet &cs,
                                   FullSquareMatrix &K, double *f, double dt, double t)
{
  int i,j,k;

  // initialize forces and stiffness matrix
  for (i = 0; i < 45; i++) f[i] = 0;
  K.zero();

  // reformat cs to accommodate shape function routine
  double X[15], Y[15], Z[15];
  for (i = 0; i < 15; i++) {
    X[i] = cs[nodeNum[i]]->x;
    Y[i] = cs[nodeNum[i]]->y;
    Z[i] = cs[nodeNum[i]]->z;
  }

  // get the nodal temperatures
  Vector ndTemps(15);
  geomState.get_temperature(15, nodeNum, ndTemps, Tref);

  // integration: loop over Gauss pts
  double Shape[15], nGrad[15][3];
  double dOmega; // det of jacobian

  for(int n = 0; n < 9; n++) {

    // compute shape functions
    dOmega = Penta15ShapeFct(Shape, nGrad, gauss3d8[n], X, Y, Z);

    // get volume
    // dOmega is here off by a factor of 1/4
    // this divides all by a factor of 4 later
    dOmega *= weight3d8[n];
    dOmega /= 4;

    // now get F_ij = dPhi_i/dX_j = x^k_i dN_k/dX_j
    double F[3][3];

    for(j = 0; j < 3; ++j)
      F[0][j] = geomState[nodeNum[ 0]].x * nGrad[ 0][j]
              + geomState[nodeNum[ 1]].x * nGrad[ 1][j]
              + geomState[nodeNum[ 2]].x * nGrad[ 2][j]
              + geomState[nodeNum[ 3]].x * nGrad[ 3][j]
              + geomState[nodeNum[ 4]].x * nGrad[ 4][j]
              + geomState[nodeNum[ 5]].x * nGrad[ 5][j]
              + geomState[nodeNum[ 6]].x * nGrad[ 6][j]
              + geomState[nodeNum[ 7]].x * nGrad[ 7][j]
              + geomState[nodeNum[ 8]].x * nGrad[ 8][j]
              + geomState[nodeNum[ 9]].x * nGrad[ 9][j]
              + geomState[nodeNum[10]].x * nGrad[10][j]
              + geomState[nodeNum[11]].x * nGrad[11][j]
              + geomState[nodeNum[12]].x * nGrad[12][j]
              + geomState[nodeNum[13]].x * nGrad[13][j]
              + geomState[nodeNum[14]].x * nGrad[14][j];

    for(j = 0; j < 3; ++j)
      F[1][j] = geomState[nodeNum[ 0]].y * nGrad[ 0][j]
              + geomState[nodeNum[ 1]].y * nGrad[ 1][j]
              + geomState[nodeNum[ 2]].y * nGrad[ 2][j]
              + geomState[nodeNum[ 3]].y * nGrad[ 3][j]
              + geomState[nodeNum[ 4]].y * nGrad[ 4][j]
              + geomState[nodeNum[ 5]].y * nGrad[ 5][j]
              + geomState[nodeNum[ 6]].y * nGrad[ 6][j]
              + geomState[nodeNum[ 7]].y * nGrad[ 7][j]
              + geomState[nodeNum[ 8]].y * nGrad[ 8][j]
              + geomState[nodeNum[ 9]].y * nGrad[ 9][j]
              + geomState[nodeNum[10]].y * nGrad[10][j]
              + geomState[nodeNum[11]].y * nGrad[11][j]
              + geomState[nodeNum[12]].y * nGrad[12][j]
              + geomState[nodeNum[13]].y * nGrad[13][j]
              + geomState[nodeNum[14]].y * nGrad[14][j];

    for(j = 0; j < 3; ++j)
      F[2][j] = geomState[nodeNum[ 0]].z * nGrad[ 0][j]
              + geomState[nodeNum[ 1]].z * nGrad[ 1][j]
              + geomState[nodeNum[ 2]].z * nGrad[ 2][j]
              + geomState[nodeNum[ 3]].z * nGrad[ 3][j]
              + geomState[nodeNum[ 4]].z * nGrad[ 4][j]
              + geomState[nodeNum[ 5]].z * nGrad[ 5][j]
              + geomState[nodeNum[ 6]].z * nGrad[ 6][j]
              + geomState[nodeNum[ 7]].z * nGrad[ 7][j]
              + geomState[nodeNum[ 8]].z * nGrad[ 8][j]
              + geomState[nodeNum[ 9]].z * nGrad[ 9][j]
              + geomState[nodeNum[10]].z * nGrad[10][j]
              + geomState[nodeNum[11]].z * nGrad[11][j]
              + geomState[nodeNum[12]].z * nGrad[12][j]
              + geomState[nodeNum[13]].z * nGrad[13][j]
              + geomState[nodeNum[14]].z * nGrad[14][j];

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
    for(j = 0; j < 15; j++) theta += Shape[j]*(ndTemps[j] - Tref);
    double alpha = (ctett) ? ctett->getValAlt(theta) : Penta15Corotator::alpha;
    e_11 -= 2*alpha*theta;
    e_22 -= 2*alpha*theta;
    e_33 -= 2*alpha*theta;

    double sigma[6];
    double em = (ymtt) ? ymtt->getValAlt(theta) : Penta15Corotator::em;
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
    double dedU[45][6];
    for(i = 0; i < 15; ++i)
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
    for(i = 0; i < 45; ++i)
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
    double dsdU[45][6];
    for(i = 0; i < 45; ++i) {
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
    for(i = 0; i < 45; ++i)
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
    for(i = 0; i < 15; ++i)
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

  // Symmetrize
  for(i = 0; i < 45; ++i)
    for(j = 0; j <= i; ++j)
      K[i][j] = K[j][i];
}

void
Penta15Corotator::getInternalForce(GeomState &geomState, CoordSet &cs,
                                   FullSquareMatrix &, double *f, double dt, double t)
{
  int i,j,k;

  // initialize forces
  for (i = 0; i < 45; i++) f[i] = 0;

  // reformat cs to accommodate shape function routine
  double X[15], Y[15], Z[15];
  for (i = 0; i < 15; i++) {
    X[i] = cs[nodeNum[i]]->x;
    Y[i] = cs[nodeNum[i]]->y;
    Z[i] = cs[nodeNum[i]]->z;
  }

  // get the nodal temperatures
  Vector ndTemps(15);
  geomState.get_temperature(15, nodeNum, ndTemps, Tref);

  // integration: loop over Gauss pts
  double Shape[15], nGrad[15][3];
  double dOmega; // det of jacobian

  for(int n = 0; n < 9; n++) { 

    // compute shape functions
    dOmega = Penta15ShapeFct(Shape, nGrad, gauss3d8[n], X, Y, Z);

    // get volume
    // dOmega is here off by a factor of 1/4
    // this divides all by a factor of 4 later
    dOmega *= weight3d8[n];
    dOmega /= 4;

    // now get F_ij = dPhi_i/dX_j = x^k_i dN_k/dX_j
    double F[3][3];
    for(j = 0; j < 3; ++j)
      F[0][j] = geomState[nodeNum[ 0]].x * nGrad[ 0][j]
              + geomState[nodeNum[ 1]].x * nGrad[ 1][j]
              + geomState[nodeNum[ 2]].x * nGrad[ 2][j]
              + geomState[nodeNum[ 3]].x * nGrad[ 3][j]
              + geomState[nodeNum[ 4]].x * nGrad[ 4][j]
              + geomState[nodeNum[ 5]].x * nGrad[ 5][j]
              + geomState[nodeNum[ 6]].x * nGrad[ 6][j]
              + geomState[nodeNum[ 7]].x * nGrad[ 7][j]
              + geomState[nodeNum[ 8]].x * nGrad[ 8][j]
              + geomState[nodeNum[ 9]].x * nGrad[ 9][j]
              + geomState[nodeNum[10]].x * nGrad[10][j]
              + geomState[nodeNum[11]].x * nGrad[11][j]
              + geomState[nodeNum[12]].x * nGrad[12][j]
              + geomState[nodeNum[13]].x * nGrad[13][j]
              + geomState[nodeNum[14]].x * nGrad[14][j];

    for(j = 0; j < 3; ++j)
      F[1][j] = geomState[nodeNum[ 0]].y * nGrad[ 0][j]
              + geomState[nodeNum[ 1]].y * nGrad[ 1][j]
              + geomState[nodeNum[ 2]].y * nGrad[ 2][j]
              + geomState[nodeNum[ 3]].y * nGrad[ 3][j]
              + geomState[nodeNum[ 4]].y * nGrad[ 4][j]
              + geomState[nodeNum[ 5]].y * nGrad[ 5][j]
              + geomState[nodeNum[ 6]].y * nGrad[ 6][j]
              + geomState[nodeNum[ 7]].y * nGrad[ 7][j]
              + geomState[nodeNum[ 8]].y * nGrad[ 8][j]
              + geomState[nodeNum[ 9]].y * nGrad[ 9][j]
              + geomState[nodeNum[10]].y * nGrad[10][j]
              + geomState[nodeNum[11]].y * nGrad[11][j]
              + geomState[nodeNum[12]].y * nGrad[12][j]
              + geomState[nodeNum[13]].y * nGrad[13][j]
              + geomState[nodeNum[14]].y * nGrad[14][j];

    for(j = 0; j < 3; ++j)
      F[2][j] = geomState[nodeNum[ 0]].z * nGrad[ 0][j]
              + geomState[nodeNum[ 1]].z * nGrad[ 1][j]
              + geomState[nodeNum[ 2]].z * nGrad[ 2][j]
              + geomState[nodeNum[ 3]].z * nGrad[ 3][j]
              + geomState[nodeNum[ 4]].z * nGrad[ 4][j]
              + geomState[nodeNum[ 5]].z * nGrad[ 5][j]
              + geomState[nodeNum[ 6]].z * nGrad[ 6][j]
              + geomState[nodeNum[ 7]].z * nGrad[ 7][j]
              + geomState[nodeNum[ 8]].z * nGrad[ 8][j]
              + geomState[nodeNum[ 9]].z * nGrad[ 9][j]
              + geomState[nodeNum[10]].z * nGrad[10][j]
              + geomState[nodeNum[11]].z * nGrad[11][j]
              + geomState[nodeNum[12]].z * nGrad[12][j]
              + geomState[nodeNum[13]].z * nGrad[13][j]
              + geomState[nodeNum[14]].z * nGrad[14][j];

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
    for(j = 0; j < 15; j++) theta += Shape[j]*(ndTemps[j] - Tref);
    double alpha = (ctett) ? ctett->getValAlt(theta) : Penta15Corotator::alpha;
    e_11 -= 2*alpha*theta;
    e_22 -= 2*alpha*theta;
    e_33 -= 2*alpha*theta;

    double sigma[6];
    double em = (ymtt) ? ymtt->getValAlt(theta) : Penta15Corotator::em;
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
    double dedU[45][6];
    for(i = 0; i < 15; ++i)
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
    for(i = 0; i < 45; ++i)
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

//-------------------------------------------------------------------------------

double
Penta15Corotator::computeStrainGrad(GeomState &geomState, CoordSet &cs,
                                    double dedU[45][6], double m[3])
{
  int i, j;
  double dOmega; // det of jacobian
  double nGrad[15][3];

  // reformat cs to accommodate shape function routine
  double X[15], Y[15], Z[15];
  for (i = 0; i < 15; i++) {
    X[i] = cs[nodeNum[i]]->x;
    Y[i] = cs[nodeNum[i]]->y;
    Z[i] = cs[nodeNum[i]]->z;
  }

  // compute shape function & derivative
  double Shape[15];
  dOmega = Penta15ShapeFct(Shape, nGrad, m, X, Y, Z);

  // now get F_ij = dPhi_i/dX_j = x^k_i dN_k/dX_j
  double F[3][3];

  for (j = 0; j < 3; ++j)
    F[0][j] = geomState[nodeNum[ 0]].x * nGrad[ 0][j]
            + geomState[nodeNum[ 1]].x * nGrad[ 1][j]
            + geomState[nodeNum[ 2]].x * nGrad[ 2][j]
            + geomState[nodeNum[ 3]].x * nGrad[ 3][j]
            + geomState[nodeNum[ 4]].x * nGrad[ 4][j]
            + geomState[nodeNum[ 5]].x * nGrad[ 5][j]
            + geomState[nodeNum[ 6]].x * nGrad[ 6][j]
            + geomState[nodeNum[ 7]].x * nGrad[ 7][j]
            + geomState[nodeNum[ 8]].x * nGrad[ 8][j]
            + geomState[nodeNum[ 9]].x * nGrad[ 9][j]
            + geomState[nodeNum[10]].x * nGrad[10][j]
            + geomState[nodeNum[11]].x * nGrad[11][j]
            + geomState[nodeNum[12]].x * nGrad[12][j]
            + geomState[nodeNum[13]].x * nGrad[13][j]
            + geomState[nodeNum[14]].x * nGrad[14][j];

  for (j = 0; j < 3; ++j)
    F[1][j] = geomState[nodeNum[ 0]].y * nGrad[ 0][j]
            + geomState[nodeNum[ 1]].y * nGrad[ 1][j]
            + geomState[nodeNum[ 2]].y * nGrad[ 2][j]
            + geomState[nodeNum[ 3]].y * nGrad[ 3][j]
            + geomState[nodeNum[ 4]].y * nGrad[ 4][j]
            + geomState[nodeNum[ 5]].y * nGrad[ 5][j]
            + geomState[nodeNum[ 6]].y * nGrad[ 6][j]
            + geomState[nodeNum[ 7]].y * nGrad[ 7][j]
            + geomState[nodeNum[ 8]].y * nGrad[ 8][j]
            + geomState[nodeNum[ 9]].y * nGrad[ 9][j]
            + geomState[nodeNum[10]].y * nGrad[10][j]
            + geomState[nodeNum[11]].y * nGrad[11][j]
            + geomState[nodeNum[12]].y * nGrad[12][j]
            + geomState[nodeNum[13]].y * nGrad[13][j]
            + geomState[nodeNum[14]].y * nGrad[14][j];

  for (j = 0; j < 3; ++j)
    F[2][j] = geomState[nodeNum[ 0]].z * nGrad[ 0][j]
            + geomState[nodeNum[ 1]].z * nGrad[ 1][j]
            + geomState[nodeNum[ 2]].z * nGrad[ 2][j]
            + geomState[nodeNum[ 3]].z * nGrad[ 3][j]
            + geomState[nodeNum[ 4]].z * nGrad[ 4][j]
            + geomState[nodeNum[ 5]].z * nGrad[ 5][j]
            + geomState[nodeNum[ 6]].z * nGrad[ 6][j]
            + geomState[nodeNum[ 7]].z * nGrad[ 7][j]
            + geomState[nodeNum[ 8]].z * nGrad[ 8][j]
            + geomState[nodeNum[ 9]].z * nGrad[ 9][j]
            + geomState[nodeNum[10]].z * nGrad[10][j]
            + geomState[nodeNum[11]].z * nGrad[11][j]
            + geomState[nodeNum[12]].z * nGrad[12][j]
            + geomState[nodeNum[13]].z * nGrad[13][j]
            + geomState[nodeNum[14]].z * nGrad[14][j];

  for (i = 0; i < 15; ++i)
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
Penta15Corotator::getNLVonMises(Vector& stress, Vector& weight, GeomState &geomState,
                                GeomState *, CoordSet& cs, int strInd, int,
                                double, double, int, int)
{
  weight = 1.0;

  int maxstr = 7;
  int elm    = 1;
  int nno    = 15;

  double elStress[15][7];
  double elStrain[15][7];

  // Compute NL Stress/Strain
  computePiolaStress(geomState, cs, elStress, elStrain);

  // Compute Von Mises
  if(strInd == 6)
    _FORTRAN(vmelmv)((double*)elStress,nno,maxstr,elm,elm,nno);
  if(strInd == 13)
    _FORTRAN(strainvm)((double*)elStrain,nno,maxstr,elm,nno);

  // Store Stress or Strain as defined by strInd
  if(strInd < 7) {
    for(int i = 0; i < 15; ++i)
      stress[i] = elStress[i][strInd];
  }
  else if(strInd < 14) {
    for(int i = 0; i < 15; ++i)
      stress[i] = elStrain[i][strInd-7];
  }
  else {
    for(int i = 0; i < 15; ++i)
      stress[i] = 0;
  }
}

void
Penta15Corotator::getNLAllStress(FullM &stress, Vector &weight, GeomState &geomState,
                                 GeomState *, CoordSet &cs, int strInd, int, int)
{
  weight = 1.0;

  int i,j;

  double elStress[15][7];
  double elStrain[15][7];

  // Compute NL Stress/Strain
  computePiolaStress(geomState, cs, elStress, elStrain);

  // Store all Stress or all Strain as defined by strInd
  if(strInd == 0) {
    for (i=0; i<15; ++i) {
      for (j=0; j<6; ++j) {
        stress[i][j] = elStress[i][j];
      }
    }
  } else {
    for (i=0; i<15; ++i) {
      for (j=0; j<6; ++j) {
        stress[i][j] = elStrain[i][j];
      }
    }
  }

  // Get Element Principals without averaging
  double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double pvec[3] = {0.0,0.0,0.0};
  for (i=0; i<15; ++i) {
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
Penta15Corotator::computePiolaStress(GeomState &geomState, CoordSet &cs,
                                     double stress[15][7], double strain[15][7])
{
  int i,j,n;
  double nGrad[15][3];
  double Shape[15];

  // coordinates of the nodes in the reference element
  double nodeRefCoord[15][3] = {{  0., 0., -1.},{  1.,  0., -1.},{ 0.,  1., -1.},
                                {  0., 0.,  1.},{  1.,  0.,  1.},{ 0.,  1.,  1.},
                                { 0.5, 0., -1.},{ 0.5, 0.5, -1.},{ 0., 0.5, -1.},
                                { 0.5, 0.,  1.},{ 0.5, 0.5,  1.},{ 0., 0.5,  1.},
                                {  0., 0.,  0.},{  1.,  0.,  0.},{ 0.,  1.,  0.}};

  // reformat cs to accommodate shape function routine
  double X[15], Y[15], Z[15];
  for (i = 0; i < 15; i++) {
    X[i] = cs[nodeNum[i]]->x;
    Y[i] = cs[nodeNum[i]]->y;
    Z[i] = cs[nodeNum[i]]->z;
  }

  // get the nodal temperatures
  Vector ndTemps(15);
  geomState.get_temperature(15, nodeNum, ndTemps, Tref);

  for (n = 0; n < 15; n++) { // loop over nodes

    // compute shape functions
    Penta15ShapeFct(Shape, nGrad, nodeRefCoord[n], X, Y, Z);

    // now get F_ij = dPhi_i/dX_j = x^k_i dN_k/dX_j
    double F[3][3];

    for(j = 0; j < 3; ++j)
      F[0][j] = geomState[nodeNum[ 0]].x * nGrad[ 0][j]
              + geomState[nodeNum[ 1]].x * nGrad[ 1][j]
              + geomState[nodeNum[ 2]].x * nGrad[ 2][j]
              + geomState[nodeNum[ 3]].x * nGrad[ 3][j]
              + geomState[nodeNum[ 4]].x * nGrad[ 4][j]
              + geomState[nodeNum[ 5]].x * nGrad[ 5][j]
              + geomState[nodeNum[ 6]].x * nGrad[ 6][j]
              + geomState[nodeNum[ 7]].x * nGrad[ 7][j]
              + geomState[nodeNum[ 8]].x * nGrad[ 8][j]
              + geomState[nodeNum[ 9]].x * nGrad[ 9][j]
              + geomState[nodeNum[10]].x * nGrad[10][j]
              + geomState[nodeNum[11]].x * nGrad[11][j]
              + geomState[nodeNum[12]].x * nGrad[12][j]
              + geomState[nodeNum[13]].x * nGrad[13][j]
              + geomState[nodeNum[14]].x * nGrad[14][j];

    for(j = 0; j < 3; ++j)
      F[1][j] = geomState[nodeNum[ 0]].y * nGrad[ 0][j]
              + geomState[nodeNum[ 1]].y * nGrad[ 1][j]
              + geomState[nodeNum[ 2]].y * nGrad[ 2][j]
              + geomState[nodeNum[ 3]].y * nGrad[ 3][j]
              + geomState[nodeNum[ 4]].y * nGrad[ 4][j]
              + geomState[nodeNum[ 5]].y * nGrad[ 5][j]
              + geomState[nodeNum[ 6]].y * nGrad[ 6][j]
              + geomState[nodeNum[ 7]].y * nGrad[ 7][j]
              + geomState[nodeNum[ 8]].y * nGrad[ 8][j]
              + geomState[nodeNum[ 9]].y * nGrad[ 9][j]
              + geomState[nodeNum[10]].y * nGrad[10][j]
              + geomState[nodeNum[11]].y * nGrad[11][j]
              + geomState[nodeNum[12]].y * nGrad[12][j]
              + geomState[nodeNum[13]].y * nGrad[13][j]
              + geomState[nodeNum[14]].y * nGrad[14][j];

    for(j = 0; j < 3; ++j)
      F[2][j] = geomState[nodeNum[ 0]].z * nGrad[ 0][j]
              + geomState[nodeNum[ 1]].z * nGrad[ 1][j]
              + geomState[nodeNum[ 2]].z * nGrad[ 2][j]
              + geomState[nodeNum[ 3]].z * nGrad[ 3][j]
              + geomState[nodeNum[ 4]].z * nGrad[ 4][j]
              + geomState[nodeNum[ 5]].z * nGrad[ 5][j]
              + geomState[nodeNum[ 6]].z * nGrad[ 6][j]
              + geomState[nodeNum[ 7]].z * nGrad[ 7][j]
              + geomState[nodeNum[ 8]].z * nGrad[ 8][j]
              + geomState[nodeNum[ 9]].z * nGrad[ 9][j]
              + geomState[nodeNum[10]].z * nGrad[10][j]
              + geomState[nodeNum[11]].z * nGrad[11][j]
              + geomState[nodeNum[12]].z * nGrad[12][j]
              + geomState[nodeNum[13]].z * nGrad[13][j]
              + geomState[nodeNum[14]].z * nGrad[14][j];

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
    double alpha = (ctett) ? ctett->getValAlt(ndTemps[n]) : Penta15Corotator::alpha;
    e_11 -= alpha*(ndTemps[n]-Tref);
    e_22 -= alpha*(ndTemps[n]-Tref);
    e_33 -= alpha*(ndTemps[n]-Tref);

    double sigma[6];
    double em = (ymtt) ? ymtt->getValAlt(ndTemps[n]) : Penta15Corotator::em;
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
Penta15Corotator::extractDeformations(GeomState &geomState, CoordSet &cs,
                                      double *vld, int &nlflag)
{
  // Set Flag to Use Non-Linear Routines for Stress
  nlflag = 2;
}

void
Penta15Corotator::extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                         double *vlr)
{
}

double
Penta15Corotator::getElementEnergy(GeomState &geomState, CoordSet &cs)
{
  // Computes Internal Energy of Element in Given State
  int i,j;

  // initialize energy to zero
  double Energy = 0;

  // reformat cs to accommodate shape function routine
  double X[15], Y[15], Z[15];
  for (i = 0; i < 15; i++) {
    X[i] = cs[nodeNum[i]]->x;
    Y[i] = cs[nodeNum[i]]->y;
    Z[i] = cs[nodeNum[i]]->z;
  }

  // get the nodal temperatures
  Vector ndTemps(15);
  geomState.get_temperature(15, nodeNum, ndTemps, Tref);

  // integration: loop over Gauss pts
  double Shape[15], nGrad[15][3];
  double dOmega; // det of jacobian

  for(int n = 0; n < 9; n++) {

    // compute shape functions
    dOmega = Penta15ShapeFct(Shape, nGrad, gauss3d8[n], X, Y, Z);

    // get volume
    dOmega *= weight3d8[n];

    // now get F_ij = dPhi_i/dX_j = x^k_i dN_k/dX_j
    double F[3][3];
    for(j = 0; j < 3; ++j)
      F[0][j] = geomState[nodeNum[ 0]].x * nGrad[ 0][j]
              + geomState[nodeNum[ 1]].x * nGrad[ 1][j]
              + geomState[nodeNum[ 2]].x * nGrad[ 2][j]
              + geomState[nodeNum[ 3]].x * nGrad[ 3][j]
              + geomState[nodeNum[ 4]].x * nGrad[ 4][j]
              + geomState[nodeNum[ 5]].x * nGrad[ 5][j]
              + geomState[nodeNum[ 6]].x * nGrad[ 6][j]
              + geomState[nodeNum[ 7]].x * nGrad[ 7][j]
              + geomState[nodeNum[ 8]].x * nGrad[ 8][j]
              + geomState[nodeNum[ 9]].x * nGrad[ 9][j]
              + geomState[nodeNum[10]].x * nGrad[10][j]
              + geomState[nodeNum[11]].x * nGrad[11][j]
              + geomState[nodeNum[12]].x * nGrad[12][j]
              + geomState[nodeNum[13]].x * nGrad[13][j]
              + geomState[nodeNum[14]].x * nGrad[14][j];

    for(j = 0; j < 3; ++j)
      F[1][j] = geomState[nodeNum[ 0]].y * nGrad[ 0][j]
              + geomState[nodeNum[ 1]].y * nGrad[ 1][j]
              + geomState[nodeNum[ 2]].y * nGrad[ 2][j]
              + geomState[nodeNum[ 3]].y * nGrad[ 3][j]
              + geomState[nodeNum[ 4]].y * nGrad[ 4][j]
              + geomState[nodeNum[ 5]].y * nGrad[ 5][j]
              + geomState[nodeNum[ 6]].y * nGrad[ 6][j]
              + geomState[nodeNum[ 7]].y * nGrad[ 7][j]
              + geomState[nodeNum[ 8]].y * nGrad[ 8][j]
              + geomState[nodeNum[ 9]].y * nGrad[ 9][j]
              + geomState[nodeNum[10]].y * nGrad[10][j]
              + geomState[nodeNum[11]].y * nGrad[11][j]
              + geomState[nodeNum[12]].y * nGrad[12][j]
              + geomState[nodeNum[13]].y * nGrad[13][j]
              + geomState[nodeNum[14]].y * nGrad[14][j];

    for(j = 0; j < 3; ++j)
      F[2][j] = geomState[nodeNum[ 0]].z * nGrad[ 0][j]
              + geomState[nodeNum[ 1]].z * nGrad[ 1][j]
              + geomState[nodeNum[ 2]].z * nGrad[ 2][j]
              + geomState[nodeNum[ 3]].z * nGrad[ 3][j]
              + geomState[nodeNum[ 4]].z * nGrad[ 4][j]
              + geomState[nodeNum[ 5]].z * nGrad[ 5][j]
              + geomState[nodeNum[ 6]].z * nGrad[ 6][j]
              + geomState[nodeNum[ 7]].z * nGrad[ 7][j]
              + geomState[nodeNum[ 8]].z * nGrad[ 8][j]
              + geomState[nodeNum[ 9]].z * nGrad[ 9][j]
              + geomState[nodeNum[10]].z * nGrad[10][j]
              + geomState[nodeNum[11]].z * nGrad[11][j]
              + geomState[nodeNum[12]].z * nGrad[12][j]
              + geomState[nodeNum[13]].z * nGrad[13][j]
              + geomState[nodeNum[14]].z * nGrad[14][j];

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
    for(j = 0; j < 15; j++) theta += Shape[j]*(ndTemps[j] - Tref);
    double alpha = (ctett) ? ctett->getValAlt(theta) : Penta15Corotator::alpha;
    e_11 -= alpha*theta;
    e_22 -= alpha*theta;
    e_33 -= alpha*theta;

    double em = (ymtt) ? ymtt->getValAlt(theta) : Penta15Corotator::em;
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
  Energy *= 0.5;
  return Energy;
}
