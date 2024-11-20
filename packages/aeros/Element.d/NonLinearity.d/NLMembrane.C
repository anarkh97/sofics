#include <Element.d/NonLinearity.d/NLMembrane.h>
#include <Element.d/NonLinearity.d/2DMat.h>
#include <Element.d/State.h>
#include <Element.d/Utils.d/SolidElemUtils.h>
#include <Element.d/NonLinearity.d/PlaneStressMat.h>
#include <Math.d/matrix.h>
#include <Math.d/TTensor.h>
#include <Parser.d/AuxDefs.h>
#include <Corotational.d/utilities.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/PhantomCorotator.h>
#include <Corotational.d/MatNLCorotator.h>
#include <Utils.d/dofset.h>
#include <Utils.d/SolverInfo.h>
#include <Element.d/FelippaShell.d/ShellElementTemplate.hpp>
#include <Element.d/FelippaShell.d/NoBendingTriangle.hpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangle.hpp>
#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif
#include <cmath>
#include "ElaLinIsoMat.h"
#include <Hetero.d/FlExchange.h>

extern SolverInfo &solInfo;

template <int n>
void
LinearStrain2D<n>::getE(typename TwoDTensorTypes<n>::StrainTensor  &e,
                        typename TwoDTensorTypes<n>::GradUTensor &gradU)
{
  e[0] = gradU[0][0];
  e[1] = gradU[1][1];
  e[2] = 0.5*(gradU[0][1]+gradU[1][0]);
}

template <int n>
void
LinearStrain2D<n>::getEandB(typename TwoDTensorTypes<n>::StrainTensor &e,
                            typename TwoDTensorTypes<n>::BTensor &B,
                            typename TwoDTensorTypes<n>::GradUTensor &gradU,
                            typename TwoDTensorTypes<n>::GradUDerivTensor &dgradUdqk)
{
  e[0] = gradU[0][0];
  e[1] = gradU[1][1];
  e[2] = 0.5*(gradU[0][1]+gradU[1][0]);
  int k;
  for(k = 0; k < n; ++k) {
    B[k][0] = dgradUdqk[k][0][0];
    B[k][1] = dgradUdqk[k][1][1];
    B[k][2] = 0.5*(dgradUdqk[k][0][1]+dgradUdqk[k][1][0]);
  }
}

template <int n>
void
LinearStrain2D<n>::getEBandDB(typename TwoDTensorTypes<n>::StrainTensor &e, 
                              typename TwoDTensorTypes<n>::BTensor &B, 
			      typename TwoDTensorTypes<n>::DBTensor &DB,
			      typename TwoDTensorTypes<n>::GradUTensor &gradU, 
			      typename TwoDTensorTypes<n>::GradUDerivTensor &dgradUdqk)
{
  e[0] = gradU[0][0];
  e[1] = gradU[1][1];
  e[2] = 0.5*(gradU[0][1]+gradU[1][0]);
  int k;
  for(k = 0; k < n; ++k) {
    B[k][0] = dgradUdqk[k][0][0];
    B[k][1] = dgradUdqk[k][1][1];
    B[k][2] = 0.5*(dgradUdqk[k][0][1]+dgradUdqk[k][1][0]);
  }
}

template <int n>
void
GLStrain2D<n>::getE(typename TwoDTensorTypes<n>::StrainTensor &e,
                    typename TwoDTensorTypes<n>::GradUTensor &gradU)
{
  // 1/2*((I+gradU)^T(I+gradU)-I)
  e[0] = gradU[0][0] +
         0.5*(gradU[0][0]*gradU[0][0]+gradU[0][1]*gradU[0][1]+gradU[0][2]*gradU[0][2]);
  e[1] = gradU[1][1] +
         0.5*(gradU[1][0]*gradU[1][0]+gradU[1][1]*gradU[1][1]+gradU[1][2]*gradU[1][2]);
  e[2] = 0.5*(gradU[0][1]+gradU[1][0] +
         gradU[1][0]*gradU[0][0]+gradU[1][1]*gradU[0][1]+gradU[1][2]*gradU[0][2]);
}

template <int n>
void
GLStrain2D<n>::getEandB(typename TwoDTensorTypes<n>::StrainTensor &e, 
                        typename TwoDTensorTypes<n>::BTensor &B, 
                        typename TwoDTensorTypes<n>::GradUTensor &gradU, 
                        typename TwoDTensorTypes<n>::GradUDerivTensor &dgradUdqk)
{
  // 1/2*((I+gradU)^T(I+gradU)-I)
  e[0] = gradU[0][0] +
         0.5*(gradU[0][0]*gradU[0][0]+gradU[0][1]*gradU[0][1]+gradU[0][2]*gradU[0][2]);
  e[1] = gradU[1][1] +
         0.5*(gradU[1][0]*gradU[1][0]+gradU[1][1]*gradU[1][1]+gradU[1][2]*gradU[1][2]);
  e[2] = 0.5*(gradU[0][1]+gradU[1][0] +
           gradU[1][0]*gradU[0][0]+gradU[1][1]*gradU[0][1]+gradU[1][2]*gradU[0][2]);

  int k;
  for(k = 0; k < n; ++k) {
    B[k][0] = dgradUdqk[k][0][0]
        + gradU[0][0]*dgradUdqk[k][0][0]
	+ gradU[0][1]*dgradUdqk[k][0][1]
	+ gradU[0][2]*dgradUdqk[k][0][2];
    B[k][1] = dgradUdqk[k][1][1]
        + gradU[1][0]*dgradUdqk[k][1][0]
	+ gradU[1][1]*dgradUdqk[k][1][1]
	+ gradU[1][2]*dgradUdqk[k][1][2];
    B[k][2] = 0.5*(dgradUdqk[k][0][1]+dgradUdqk[k][1][0]
        + gradU[0][0]*dgradUdqk[k][1][0]
	+ gradU[0][1]*dgradUdqk[k][1][1]
	+ gradU[0][2]*dgradUdqk[k][1][2]
	+ gradU[1][0]*dgradUdqk[k][0][0]
	+ gradU[1][1]*dgradUdqk[k][0][1]
	+ gradU[1][2]*dgradUdqk[k][0][2]);
  }

}

template <int n>
void
GLStrain2D<n>::getEBandDB(typename TwoDTensorTypes<n>::StrainTensor &e, 
                          typename TwoDTensorTypes<n>::BTensor &B, 
			  typename TwoDTensorTypes<n>::DBTensor &DB,
			  typename TwoDTensorTypes<n>::GradUTensor &gradU, 
			  typename TwoDTensorTypes<n>::GradUDerivTensor &dgradUdqk)
{
  // 1/2*((I+gradU)^T(I+gradU)-I)
  e[0] = gradU[0][0] +
         0.5*(gradU[0][0]*gradU[0][0]+gradU[0][1]*gradU[0][1]+gradU[0][2]*gradU[0][2]);
  e[1] = gradU[1][1] +
         0.5*(gradU[1][0]*gradU[1][0]+gradU[1][1]*gradU[1][1]+gradU[1][2]*gradU[1][2]);
  e[2] = 0.5*(gradU[0][1]+gradU[1][0] +
           gradU[1][0]*gradU[0][0]+gradU[1][1]*gradU[0][1]+gradU[1][2]*gradU[0][2]);

  int k;
  for(k = 0; k < n; ++k) {
    B[k][0] = dgradUdqk[k][0][0]
        + gradU[0][0]*dgradUdqk[k][0][0]
	+ gradU[0][1]*dgradUdqk[k][0][1]
	+ gradU[0][2]*dgradUdqk[k][0][2];
    B[k][1] = dgradUdqk[k][1][1]
        + gradU[1][0]*dgradUdqk[k][1][0]
	+ gradU[1][1]*dgradUdqk[k][1][1]
	+ gradU[1][2]*dgradUdqk[k][1][2];
    B[k][2] = 0.5*(dgradUdqk[k][0][1]+dgradUdqk[k][1][0]
        + gradU[0][0]*dgradUdqk[k][1][0]
	+ gradU[0][1]*dgradUdqk[k][1][1]
	+ gradU[0][2]*dgradUdqk[k][1][2]
	+ gradU[1][0]*dgradUdqk[k][0][0]
	+ gradU[1][1]*dgradUdqk[k][0][1]
	+ gradU[1][2]*dgradUdqk[k][0][2]);
  }
  int l;
  for(k = 0; k < n; ++k) {
    for(l=0; l < n; ++l) {
      DB[k][l][0] = dgradUdqk[k][0][0]*dgradUdqk[l][0][0]
                  + dgradUdqk[k][0][1]*dgradUdqk[l][0][1]
                  + dgradUdqk[k][0][2]*dgradUdqk[l][0][2];
      DB[k][l][1] = dgradUdqk[k][1][0]*dgradUdqk[l][1][0]
                  + dgradUdqk[k][1][1]*dgradUdqk[l][1][1]
                  + dgradUdqk[k][1][2]*dgradUdqk[l][1][2];
      DB[k][l][2] = 0.5*( dgradUdqk[k][0][0]*dgradUdqk[l][1][0]
                  + dgradUdqk[k][0][1]*dgradUdqk[l][1][1]
                  + dgradUdqk[k][0][2]*dgradUdqk[l][1][2]
		  + dgradUdqk[l][0][0]*dgradUdqk[k][1][0]
                  + dgradUdqk[l][0][1]*dgradUdqk[k][1][1]
                  + dgradUdqk[l][0][2]*dgradUdqk[k][1][2]);
    }
  }

}

void
TriMembraneShapeFunct::getGlobalGrads(Grad2D &gradU, Grad2DDeriv9 &dGradUdqk,
                                      double *_jac,
                                      Node *nodes, double xi[3], Vector &disp)
{
  int i,j,k;
  // First obtain the frame X, in which xi[i] runs along X[i].
  double d[2][3] = { { nodes[1].x-nodes[0].x, 
                       nodes[1].y-nodes[0].y,
                       nodes[1].z-nodes[0].z },
                     { nodes[2].x-nodes[0].x, 
                       nodes[2].y-nodes[0].y,
                       nodes[2].z-nodes[0].z },
		   };

  double X[3][3];
  X[2][0] = d[0][1]*d[1][2] - d[0][2]*d[1][1];
  X[2][1] = d[0][2]*d[1][0] - d[0][0]*d[1][2];
  X[2][2] = d[0][0]*d[1][1] - d[0][1]*d[1][0];
 
  double l1inv = 1.0/sqrt(d[0][0]*d[0][0]+d[0][1]*d[0][1]+d[0][2]*d[0][2]);
  double l3inv = 1.0/sqrt(X[2][0]*X[2][0]+X[2][1]*X[2][1]+X[2][2]*X[2][2]);
  for(i = 0; i < 3; ++i) {
    X[0][i] = l1inv*d[0][i];
    X[2][i] *= l3inv;
  }
  X[1][0] = X[2][1]*X[0][2] - X[2][2]*X[0][1];
  X[1][1] = X[2][2]*X[0][0] - X[2][0]*X[0][2];
  X[1][2] = X[2][0]*X[0][1] - X[2][1]*X[0][0];

  // coordinates of point P in global frame M = (1-ξ-η)*nodes[0] + ξ*nodes[1] + η*nodes[2]
  // ∂M/∂ξ = nodes[1]-nodes[0] = d[0]
  // ∂M/∂η = nodes[2]-nodes[0] = d[1]
  // coordinates of point P in element frame M' = X*M
  // ∂M'/∂ξ = X*d[0]
  // ∂M'/∂η = X*d[1]

  // This is: [ ∂M'₀/∂ξ  ∂M'₁/∂ξ ]
  //          [ ∂M'₀/∂η  ∂M'₁/∂η ]
  double dXdxi[2][2] = { { d[0][0]*X[0][0]+d[0][1]*X[0][1]+d[0][2]*X[0][2],
                           d[0][0]*X[1][0]+d[0][1]*X[1][1]+d[0][2]*X[1][2] },
                         { d[1][0]*X[0][0]+d[1][1]*X[0][1]+d[1][2]*X[0][2],
                           d[1][0]*X[1][0]+d[1][1]*X[1][1]+d[1][2]*X[1][2] }
		       };
  double &jac = *_jac;
  jac = dXdxi[0][0]*dXdxi[1][1]-dXdxi[0][1]*dXdxi[1][0];
  double jacInv = 1.0/jac;
  // This is the transpose of the inverse of dXdxi, i.e: [ ∂ξ/∂M'₀ ∂ξ/∂M'₁ ]  
  //                                                     [ ∂η/∂M'₀ ∂η/∂M'₁ ]
  double dxidX[2][2] = { { jacInv*dXdxi[1][1], -jacInv*dXdxi[1][0] },
                         { -jacInv*dXdxi[0][1], jacInv*dXdxi[0][0] }
                       };
  // Create the local displacement vectors.
  double Utilde[3][3] = 
    { { X[0][0]*disp[0] + X[0][1]*disp[1] + X[0][2]*disp[2],
        X[1][0]*disp[0] + X[1][1]*disp[1] + X[1][2]*disp[2],
	X[2][0]*disp[0] + X[2][1]*disp[1] + X[2][2]*disp[2]
      },
      { X[0][0]*disp[3] + X[0][1]*disp[4] + X[0][2]*disp[5],
        X[1][0]*disp[3] + X[1][1]*disp[4] + X[1][2]*disp[5],
	X[2][0]*disp[3] + X[2][1]*disp[4] + X[2][2]*disp[5]
      },
      { X[0][0]*disp[6] + X[0][1]*disp[7] + X[0][2]*disp[8],
        X[1][0]*disp[6] + X[1][1]*disp[7] + X[1][2]*disp[8],
	X[2][0]*disp[6] + X[2][1]*disp[7] + X[2][2]*disp[8]
      } };
  // Derivatives of the scalar shape functions:
  // N = [ N₀ ] = [ 1-ξ-η ]
  //     [ N₁ ]   [   ξ   ]
  //     [ N₂ ]   [   η   ]
  // dNdx = [ ∂N₀/∂ξ  ∂N₀/∂η ] = [ -1 -1 ]
  //        [ ∂N₁/∂ξ  ∂N₁/∂η ]   [  1  0 ]
  //        [ ∂N₂/∂ξ  ∂N₂/∂η ]   [  0  1 ]
  // dNdX = dNdxi*dxidX
  //      = [ ∂N₀/∂ξ  ∂N₀/∂η ] * [ ∂ξ/∂M'₀ ∂ξ/∂M'₁ ]
  //        [ ∂N₁/∂ξ  ∂N₁/∂η ]   [ ∂η/∂M'₀ ∂η/∂M'₁ ]
  //        [ ∂N₂/∂ξ  ∂N₂/∂η ]
  //      = [ -1 -1 ]   [ dxidX_00 dxidX_01 ]
  //        [  1  0 ] * [ dxidX_10 dxidX_11 ]
  //        [  0  1 ]
  double dNdX[3][2] = {  { -(dxidX[0][0] + dxidX[1][0]), 
                           -(dxidX[0][1] + dxidX[1][1]) },
                         { dxidX[0][0], dxidX[0][1] },
                         { dxidX[1][0], dxidX[1][1] }
		      };
  for(i = 0; i < 2; ++i) {
    gradU[i][0] = dNdX[0][i]*Utilde[0][0]
                + dNdX[1][i]*Utilde[1][0]
                + dNdX[2][i]*Utilde[2][0];
    gradU[i][1] = dNdX[0][i]*Utilde[0][1]
                + dNdX[1][i]*Utilde[1][1]
                + dNdX[2][i]*Utilde[2][1];
    gradU[i][2] = dNdX[0][i]*Utilde[0][2]
                + dNdX[1][i]*Utilde[1][2]
                + dNdX[2][i]*Utilde[2][2];
  }

  for(k = 0; k < 9; ++k) {
    int dir, nd;
    dir = k%3;
    nd = k/3;
    for(i = 0; i < 2; ++i)
      for(j = 0; j < 3; ++j) {
        dGradUdqk[k][i][j] = dNdX[nd][i]*X[j][dir];
      }
  }
  // The volume of the element is only half of the jacobian.
  jac *=0.5;
}

void
TriMembraneShapeFunct::getGradU(Grad2D &gradU, Node *nodes, double xi[3], Vector &disp)
{
  int i;
  // First obtain the frame X, in which xi[i] runs along X[i].
  double d[2][3] = { { nodes[1].x-nodes[0].x, 
                       nodes[1].y-nodes[0].y,
                       nodes[1].z-nodes[0].z },
                     { nodes[2].x-nodes[0].x, 
                       nodes[2].y-nodes[0].y,
                       nodes[2].z-nodes[0].z },
		   };
  double X[3][3];
  X[2][0] = d[0][1]*d[1][2] - d[0][2]*d[1][1];
  X[2][1] = d[0][2]*d[1][0] - d[0][0]*d[1][2];
  X[2][2] = d[0][0]*d[1][1] - d[0][1]*d[1][0];
 
  double l1inv = 1.0/sqrt(d[0][0]*d[0][0]+d[0][1]*d[0][1]+d[0][2]*d[0][2]);
  double l3inv = 1.0/sqrt(X[2][0]*X[2][0]+X[2][1]*X[2][1]+X[2][2]*X[2][2]);
  for(i = 0; i < 3; ++i) {
    X[0][i] = l1inv*d[0][i];
    X[2][i] *= l3inv;
  }
  X[1][0] = X[2][1]*X[0][2] - X[2][2]*X[0][1];
  X[1][1] = X[2][2]*X[0][0] - X[2][0]*X[0][2];
  X[1][2] = X[2][0]*X[0][1] - X[2][1]*X[0][0];
 
  double dXdxi[2][2] = { { d[0][0]*X[0][0]+d[0][1]*X[0][1]+d[0][2]*X[0][2],
                           d[0][0]*X[1][0]+d[0][1]*X[1][1]+d[0][2]*X[1][2] },
                         { d[1][0]*X[0][0]+d[1][1]*X[0][1]+d[1][2]*X[0][2],
                           d[1][0]*X[1][0]+d[1][1]*X[1][1]+d[1][2]*X[1][2] }
		       };
  double jac;
  jac = dXdxi[0][0]*dXdxi[1][1]-dXdxi[0][1]*dXdxi[1][0];
  double jacInv = 1.0/jac;
  double dxidX[2][2] = { { jacInv*dXdxi[1][1], -jacInv*dXdxi[1][0] },
                         { -jacInv*dXdxi[0][1], jacInv*dXdxi[0][0] }
                       };
  // Create the local displacement vectors.
  double Utilde[3][3] = 
    { { X[0][0]*disp[0] + X[0][1]*disp[1] + X[0][2]*disp[2],
        X[1][0]*disp[0] + X[1][1]*disp[1] + X[1][2]*disp[2],
	X[2][0]*disp[0] + X[2][1]*disp[1] + X[2][2]*disp[2]
      },
      { X[0][0]*disp[3] + X[0][1]*disp[4] + X[0][2]*disp[5],
        X[1][0]*disp[3] + X[1][1]*disp[4] + X[1][2]*disp[5],
	X[2][0]*disp[3] + X[2][1]*disp[4] + X[2][2]*disp[5]
      },
      { X[0][0]*disp[6] + X[0][1]*disp[7] + X[0][2]*disp[8],
        X[1][0]*disp[6] + X[1][1]*disp[7] + X[1][2]*disp[8],
	X[2][0]*disp[6] + X[2][1]*disp[7] + X[2][2]*disp[8]
      } };
  // Derivatives of the scalar shape functions:
  double dNdX[3][2] = {  { -(dxidX[0][0] + dxidX[1][0]), 
                           -(dxidX[0][1] + dxidX[1][1]) },
                         { dxidX[0][0], dxidX[0][1] },
                         { dxidX[1][0], dxidX[1][1] }
		      };
  for(i = 0; i < 2; ++i) {
    gradU[i][0] = dNdX[0][i]*Utilde[0][0]
                + dNdX[1][i]*Utilde[1][0]
                + dNdX[2][i]*Utilde[2][0];
    gradU[i][1] = dNdX[0][i]*Utilde[0][1]
                + dNdX[1][i]*Utilde[1][1]
                + dNdX[2][i]*Utilde[2][1];
    gradU[i][2] = dNdX[0][i]*Utilde[0][2]
                + dNdX[1][i]*Utilde[1][2]
                + dNdX[2][i]*Utilde[2][2];
  }
}

double
TriMembraneShapeFunct::interpolateScalar(double *q, double xi[3])
{
  return xi[0]*q[1] + xi[1]*q[2] + (1-xi[0]-xi[1])*q[0];
}

static TriMembraneShapeFunct shpFct;

NLMembrane::NLMembrane(int *nd)
 : material(NULL), linearMaterial(NULL), useDefaultMaterial(false), cFrame(NULL), cCoefs(NULL)
{
  for(int i = 0; i < 3; ++i)
    n[i] = nd[i];
  pbc = 0;
}

NLMembrane::~NLMembrane()
{
  if(material && (useDefaultMaterial || cCoefs)) delete material;
  if(linearMaterial) delete linearMaterial;
}

int
NLMembrane::getNumGaussPoints() const
{
  return 1;
  //return 3;
}

void
NLMembrane::getGaussPointAndWeight(int n, double *point, double &weight) const
{
  const double third = 1.0/3.0;
  point[0] = third;
  point[1] = third;
  point[2] = 0.0;
  weight = 1.0;
/*
  double w_save[3] = {
    0.33333333333333333333,
    0.33333333333333333333,
    0.33333333333333333333 };
  double xy_save[2*3] = {
    0.66666666666666666667,  0.16666666666666666667,
    0.16666666666666666667,  0.66666666666666666667,
    0.16666666666666666667,  0.16666666666666666667 };

  point[0] = xy_save[2*n+0];
  point[1] = xy_save[2*n+1];
  point[2] = 0.0;
  weight = w_save[n];
*/
}

void
NLMembrane::getLocalNodalCoords(int n, double *point)
{
  switch(n) {
    case 0: {
      point[0] = 0;
      point[1] = 0;
    } break;
    case 1: {
      point[0] = 1;
      point[1] = 0;
    } break;
    case 2: {
      point[0] = 0;
      point[1] = 1;
    } break;
  }
}

void
NLMembrane::renum(const int *table)
{
  n[0] = table[n[0]];
  n[1] = table[n[1]];
  n[2] = table[n[2]];
}

void
NLMembrane::renum(EleRenumMap& table)
{
  n[0] = table[n[0]];
  n[1] = table[n[1]];
  n[2] = table[n[2]];
}

int*
NLMembrane::nodes(int *nd) const
{
  if(nd == 0) nd = new int[3];
  nd[0] = n[0];
  nd[1] = n[1];
  nd[2] = n[2];
  return nd;
}

int *
NLMembrane::dofs(DofSetArray &dsa, int*df) const
{
  if(df == 0) df = new int[9];
  dsa.number(n[0], DofSet::XYZdisp, df);
  dsa.number(n[1], DofSet::XYZdisp, df+3);
  dsa.number(n[2], DofSet::XYZdisp, df+6);

  return df;
}

LinearStrain2D<9> linStrain2D;
GLStrain2D<9> glStrain2D;

GenStrainEvaluator<TwoDTensorTypes<9> > *
NLMembrane::getGenStrainEvaluator() const
{
  return material->getGenStrainEvaluator();
}

const NLMaterial *
NLMembrane::getMaterial() const
{
  return material;
}

NLMaterial *
NLMembrane::getLinearMaterial() const
{
  return linearMaterial;
}

GenShapeFunction< TwoDTensorTypes<9> > *
NLMembrane::getShapeFunction() const
{
 return &shpFct;
}

void
NLMembrane::markDofs(DofSetArray &dsa) const
{
  dsa.mark(n, 3, DofSet::XYZdisp);
}

void
NLMembrane::setProp(StructProp *p, bool _myProp)
{
  Element::setProp(p, _myProp);
}

void
NLMembrane::buildFrame(CoordSet &cs) {
#ifdef USE_EIGEN3
  int elm = glNum+1;
  Node &nd1 = cs.getNode(n[0]);
  Node &nd2 = cs.getNode(n[1]);
  Node &nd3 = cs.getNode(n[2]);
  double x[3], y[3], z[3];
  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  if(!eframe) eframe = new double[9];
  double xlp[3], ylp[3], zlp[3];
  double area;

  typedef ShellElementTemplate<double,EffMembraneTriangle,NoBendingTriangle> Impl;
  Impl::andescrd(elm, x, y, z, eframe, xlp, ylp, zlp, area);
#endif
}

void
NLMembrane::rotateCFrame(const CoordSet &cs, double *_T, double *_Tinv) const
{
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor>> T(_T), Tinv(_Tinv);
  Tinv = ShellMaterial<double>::andesinvt(eframe, cFrame, 0.);
  T = Tinv.inverse();
#else
  std::cerr << " *** ERROR: NLMembrane::rotateCFrame requires AERO-S to be configured with the Eigen library. Exiting...\n";
  exit(-1);
#endif
}

void
NLMembrane::setCompositeData(int type, int, double *, double *coefs, double *frame)
{
  if(type != 1) {
    std::cerr << " *** ERROR: only COEF-type composite is supported for NLMembrane element. Exiting...\n";
    exit(-1);
  }
  cCoefs = coefs;
  cFrame = frame;
}

double *
NLMembrane::setCompositeData2(int type, int nlays, double *lData,
                              double *coefs, CoordSet &cs, double theta)
{
  // cframe is not pre-defined but calculated from nodal coordinates and angle theta
  // theta is the angle in degrees between node1-node2 and the material x axis

  if(type != 1) {
    std::cerr << " *** ERROR: only COEF-type composite is supported for NLMembrane element. Exiting...\n";
    exit(-1);
  }

  // compute cFrame
  cFrame = new double[9];

  auto &nd1 = cs.getNode(n[0]);
  auto &nd2 = cs.getNode(n[1]);
  auto &nd3 = cs.getNode(n[2]);

  double ab[3], ac[3], x[3], y[3], z[3];
  ab[0] = nd2.x - nd1.x; ab[1] = nd2.y - nd1.y; ab[2] = nd2.z - nd1.z;
  ac[0] = nd3.x - nd1.x; ac[1] = nd3.y - nd1.y; ac[2] = nd3.z - nd1.z;
  // x = AB
  x[0] = ab[0]; x[1] = ab[1]; x[2] = ab[2];
  // z = AB cross AC
  z[0] = ab[1]*ac[2]-ab[2]*ac[1];
  z[1] = ab[2]*ac[0]-ab[0]*ac[2];
  z[2] = ab[0]*ac[1]-ab[1]*ac[0];
  // y = z cross x
  y[0] = z[1]*x[2]-z[2]*x[1];
  y[1] = z[2]*x[0]-z[0]*x[2];
  y[2] = z[0]*x[1]-z[1]*x[0];
  // rotate x and y about z
  theta *= M_PI/180.; // convert to radians
  double c = cos(theta), s = sin(theta);

  // use Rodrigues' Rotation Formula to rotation x and y about z by an angle theta
  double R[3][3];
  normalize(x); normalize(y); normalize(z); double wx = z[0], wy = z[1], wz = z[2];
  R[0][0] = c + wx*wx*(1-c);
  R[0][1] = wx*wy*(1-c)-wz*s;
  R[0][2] = wy*s+wx*wz*(1-c);
  R[1][0] = wz*s+wx*wy*(1-c);
  R[1][1] = c+wy*wy*(1-c);
  R[1][2] = -wx*s+wy*wz*(1-c);
  R[2][0] = -wy*s+wx*wz*(1-c);
  R[2][1] = wx*s+wy*wz*(1-c);
  R[2][2] = c+wz*wz*(1-c);

  cFrame[0] = R[0][0]*x[0] + R[0][1]*x[1] + R[0][2]*x[2];
  cFrame[1] = R[1][0]*x[0] + R[1][1]*x[1] + R[1][2]*x[2];
  cFrame[2] = R[2][0]*x[0] + R[2][1]*x[1] + R[2][2]*x[2];
  cFrame[3] = R[0][0]*y[0] + R[0][1]*y[1] + R[0][2]*y[2];
  cFrame[4] = R[1][0]*y[0] + R[1][1]*y[1] + R[1][2]*y[2];
  cFrame[5] = R[2][0]*y[0] + R[2][1]*y[1] + R[2][2]*y[2];
  cFrame[6] = z[0];
  cFrame[7] = z[1];
  cFrame[8] = z[2];

  setCompositeData(type, nlays, lData, coefs, cFrame);
  return cFrame;
}

void
NLMembrane::setMaterial(NLMaterial *m)
{
  if(cCoefs) { // anisotropic material
    material = m->clone();
  }
  else {
    material = m;
  }
}

void
NLMembrane::computePressureForce(CoordSet& cs, Vector& force,
                                 GeomState *geomState, int cflg, double)
{
  Node nodes[3];
  double gs[9];
  for(int i = 0; i < 3; ++i) { 
    nodes[i] = *cs[n[i]];
    if(geomState) {
      gs[3*i  ] = (*geomState)[n[i]].x;
      gs[3*i+1] = (*geomState)[n[i]].y;
      gs[3*i+2] = (*geomState)[n[i]].z;
    }
    else {
      gs[3*i] = gs[3*i+1] = gs[3*i+2] = 0;
    }
  }

  double d[2][3] = { { nodes[1].x+gs[3]-nodes[0].x-gs[0], 
                       nodes[1].y+gs[4]-nodes[0].y-gs[1],
                       nodes[1].z+gs[5]-nodes[0].z-gs[2] },
                     { nodes[2].x+gs[6]-nodes[0].x-gs[0], 
                       nodes[2].y+gs[7]-nodes[0].y-gs[1],
                       nodes[2].z+gs[8]-nodes[0].z-gs[2] },
                   };
  double n[3];
  n[0] = d[0][1]*d[1][2] - d[0][2]*d[1][1];
  n[1] = d[0][2]*d[1][0] - d[0][0]*d[1][2];
  n[2] = d[0][0]*d[1][1] - d[0][1]*d[1][0];
  double p = getPressure()->val;

  if(solInfo.nlmembrane_pressure_type == 1) { // area computed using reference/undeformed configuration
    double d0[2][3] = { { nodes[1].x-nodes[0].x,
                          nodes[1].y-nodes[0].y,
                          nodes[1].z-nodes[0].z },
                        { nodes[2].x-nodes[0].x,
                          nodes[2].y-nodes[0].y,
                          nodes[2].z-nodes[0].z },
                      };
    double n0[3];
    n0[0] = d0[0][1]*d0[1][2] - d0[0][2]*d0[1][1];
    n0[1] = d0[0][2]*d0[1][0] - d0[0][0]*d0[1][2];
    n0[2] = d0[0][0]*d0[1][1] - d0[0][1]*d0[1][0];
    double A0 = 0.5*sqrt(n0[0]*n0[0]+n0[1]*n0[1]+n0[2]*n0[2]);
    double A = 0.5*sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
    for(int i = 0; i < 3; ++i)
      force[i] = force[i+3] = force[i+6] = 1.0/6.0*p*n[i]*A0/A;
  }
  else { // area computed using current/deformed configuration
    for(int i = 0; i < 3; ++i)
      force[i] = force[i+3] = force[i+6] = 1.0/6.0*p*n[i];
  }
}

Corotator*
NLMembrane::getCorotator(CoordSet &cs, double *, int , int)
{
  if(prop == NULL) return new PhantomCorotator();
  else {
    if(!material) { // if no nonlinear material has been assigned, construct a default one
      useDefaultMaterial = true;
      if(cCoefs)
        material = new PlaneStressMat<StVenantKirchhoffMat>(prop->rho, prop->E, prop->nu, prop->eh);
      else
        material = new StVenantKirchhoffMat2D(prop);
    }
    if(cCoefs) {
      material->setTangentMaterial(reinterpret_cast<double(*)[6]>(cCoefs));
      material->setThermalExpansionCoef(cCoefs+36);
    }
    material->setTDProps(prop->ymtt, prop->ctett);
    if(cFrame && (this->tframe == NULL)) {
      this->tframe = new double[9];
      this->tframe_inv = new double[9];
      rotateCFrame(cs, this->tframe, this->tframe_inv);
    }
    return new MatNLCorotator(this, false);
  }
}

FullSquareMatrix
NLMembrane::stiffness(const CoordSet& cs, double *k, int flg) const
{
	if(prop == NULL) {
		FullSquareMatrix ret(9,k);
		ret.zero();
		return ret;
	}
	else {
		// TODO Restore constness!!!!
		NLMembrane *nConstThis = const_cast<NLMembrane *>(this);

		if(!linearMaterial) {
			if(cCoefs)
				nConstThis->linearMaterial = new PlaneStressMat<ElaLinIsoMat>(prop->rho, prop->E, prop->nu, prop->eh);
			else
				nConstThis->linearMaterial = new ElaLinIsoMat2D(prop);
		}
		if(cCoefs) {
      linearMaterial->setTangentMaterial(reinterpret_cast<double(*)[6]>(cCoefs));
      linearMaterial->setThermalExpansionCoef(cCoefs+36);
		}
		linearMaterial->setTDProps(prop->ymtt, prop->ctett);
    if(cFrame && (this->tframe == NULL)) {
      // TODO Get rid of this mutablility!
      this->tframe = new double[9];
      this->tframe_inv = new double[9];
      rotateCFrame(cs, this->tframe, this->tframe_inv);
    }
		return GenGaussIntgElement<TwoDTensorTypes<9> >::stiffness(cs,k,flg);
	}
}

FullSquareMatrix
NLMembrane::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
  FullSquareMatrix ret(9,mel);
  if(prop == NULL) { ret.zero(); return ret; }

  if(cmflg) { // consistent mass matrix
#ifdef USE_EIGEN3
    double mass = getMass(cs);
    Eigen::Map<Eigen::Matrix<double,9,9> > M(mel);
    M << 2, 0, 0, 1, 0, 0, 1, 0, 0,
         0, 2, 0, 0, 1, 0, 0, 1, 0,
         0, 0, 2, 0, 0, 1, 0, 0, 1,
         1, 0, 0, 2, 0, 0, 1, 0, 0,
         0, 1, 0, 0, 2, 0, 0, 1, 0,
         0, 0, 1, 0, 0, 2, 0, 0, 1,
         1, 0, 0, 1, 0, 0, 2, 0, 0,
         0, 1, 0, 0, 1, 0, 0, 2, 0,
         0, 0, 1, 0, 0, 1, 0, 0, 2;
     M *= mass/12;
#else
     std::cerr << " *** ERROR: Consistent mass matrix for NLMembrane element requires AERO-S configured with Eigen library. Exiting...\n";
     exit(-1);
#endif
  }
  else { // lumped mass matrix

    double mass = getMass(cs);
    double massPerNode = mass/3.0;

    ret.zero();
    for(int i = 0; i < 9; ++i)
      ret[i][i] = massPerNode;
  }

  return ret;
}

double
NLMembrane::getMass(const CoordSet& cs) const
{
  if(prop == NULL) return 0;

  auto &nd1 = cs.getNode(n[0]);
  auto &nd2 = cs.getNode(n[1]);
  auto &nd3 = cs.getNode(n[2]);

  double r1[3], r2[3], r3[3], v1[3], v2[3], v3[3];

  r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
  r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
  r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;

  v1[0] = r3[0] - r1[0];
  v1[1] = r3[1] - r1[1];
  v1[2] = r3[2] - r1[2];

  v2[0] = r2[0] - r1[0];
  v2[1] = r2[1] - r1[1];
  v2[2] = r2[2] - r1[2];

  crossprod(v1, v2, v3);

  double area = 0.5*sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
  double density = prop->rho;
  double t       = prop->eh;

  double mass = area*t*density;

  return mass;
}

void
NLMembrane::computeDisp(CoordSet &cs, State &state, const InterpPoint &ip,
                        double *res, GeomState *gs)
{
  const double *gp = ip.xy;
  double xyz[3][6];
  state.getDV(n[0], xyz[0], xyz[0]+3);
  state.getDV(n[1], xyz[1], xyz[1]+3);
  state.getDV(n[2], xyz[2], xyz[2]+3);

  for(int j = 0; j < 6; ++j)
    res[j] = (1.0-gp[0]-gp[1]) * xyz[0][j] + gp[0]*xyz[1][j] + gp[1]*xyz[2][j];
}

void
NLMembrane::getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF,
                      double *resF, GeomState *gs)
{
  const double *gp = ip.xy;
  for(int i = 0; i < 3; ++i) {
    resF[i]   = (1.0-gp[0]-gp[1]) * flF[i];
    resF[3+i] = gp[0] * flF[i];
    resF[6+i] = gp[1] * flF[i];
  }
}

// Four node membrane comprising two three node membranes, 3 dof per node
NLMembrane4::NLMembrane4(int *nodenums)
{
  int i,j,k;
  nn = new int[4];
  for(i=0; i<4; ++i) nn[i] = nodenums[i];

  nSubElems = 2;
  subElems = new Element * [2];
  subElemNodes = new int * [2];
  subElemDofs = new int * [2];

  subElemNodes[0] = new int[3];
  subElemNodes[0][0] = 0; subElemNodes[0][1] = 1; subElemNodes[0][2] = 3;

  subElemNodes[1] = new int[3];
  subElemNodes[1][0] = 2; subElemNodes[1][1] = 3; subElemNodes[1][2] = 1;

  for(i=0; i<2; ++i) {
    int tmp[3];
    subElemDofs[i] = new int[9];
    for(j=0; j<3; ++j) {
      int nij = subElemNodes[i][j];
      tmp[j] = nodenums[nij]; // global node numbers
      for(k=0;k<3;++k) {
        subElemDofs[i][3*j+k] = 3*nij+k;
      }
    }
    subElems[i] = new NLMembrane(tmp);
  }
  nnodes = 4;
  ndofs = 12;
}

int
NLMembrane4::getTopNumber() const
{
  return 102;
}

void
NLMembrane4::computeDisp(CoordSet &cs, State &state, const InterpPoint &ip, double *res, GeomState *gs)
{
  const double *gp = ip.xy;
  double gpsum = gp[0] + gp[1];
  int i;
  InterpPoint subip;
  if(gpsum <= 1.) {
    i = 0;
    subip.xy[0] = gp[0];
    subip.xy[1] = gp[1];
  }
  else {
    i = 1;
    subip.xy[0] = 1.0 - gp[0];
    subip.xy[1] = 1.0 - gp[1];
  }

  subElems[i]->computeDisp(cs, state, subip, res);
}

void
NLMembrane4::getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF,
                       double *res, GeomState *gs)
{
  const double *gp = ip.xy;
  double gpsum = gp[0] + gp[1];
  int i;
  InterpPoint subip;
  if(gpsum <= 1.) {
    i = 0;
    subip.xy[0] = gp[0];
    subip.xy[1] = gp[1];
  }
  else {
    i = 1;
    subip.xy[0] = 1.0 - gp[0];
    subip.xy[1] = 1.0 - gp[1];
  }

  double subres[9];
  subElems[i]->getFlLoad(cs, subip, flF, subres);

  int j;
  for(j=0; j<12; ++j) res[j] = 0.0;
  for(j=0; j<9; ++j) res[subElemDofs[i][j]] = subres[j];
}

