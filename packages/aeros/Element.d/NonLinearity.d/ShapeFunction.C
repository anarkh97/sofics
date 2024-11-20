#include <cstdio>
#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/ShapeFunction.h>

Tensor *
ShapeFunction::getGradUInstance()
{
  return (new Tensor_d0s2());
}

Tensor *
ShapeFunction::getDgradUDqkInstance()
{
  return (new Tensor_d1s2_sparse(numdofs));
}

Tensor *
ShapeFunction::getNodesCoordinatesInstance()
{
  return (new Tensor_d1s0(numdofs));
}

Tensor *
ShapeFunction::getDisplacementsInstance()
{
  return (new Tensor_d1s0(numdofs));
}

Tensor *
ShapeFunction::getLocalDerivativesInstance()
{
  return (new Tensor_d1s2_sparse(numdofs));
}

void
ShapeFunction::getGlobalGrads(Tensor *_gradU, Tensor *_dgradUdqk, double *jac,
                              Node *nodes, double xi[3], Vector &disp)
{
  Tensor_d1s2_sparse localderivatives(numdofs);
  Tensor_d1s0 nodescoordinates(numdofs);
  Tensor_d0s2 localGrad;
  Tensor_d0s2 jacobian;
  Tensor_d0s2 invjacobian;

  // localDerivatives(i,j,k) = at node i, the derivative of j^{th} shape function w.r.t xi[k]
  getLocalDerivatives(&localderivatives, xi);

  for (int j = 0; j < numdofs/3; j++) {
    Node &nd = nodes[j]; 
    nodescoordinates[3*j] = nd.x;
    nodescoordinates[3*j+1] = nd.y;
    nodescoordinates[3*j+2] = nd.z;
  }

  //isoparametric elements
  jacobian = localderivatives%nodescoordinates; // dof contraction

  *jac = jacobian.getInverseAndDeterminant(invjacobian);

  Tensor_d1s0 displacements(numdofs);
  disp.vectorToTensor(displacements);

  localGrad = localderivatives%displacements; // dof contraction

  Tensor_d1s2_sparse * dgradUdqk = static_cast<Tensor_d1s2_sparse *>(_dgradUdqk);
  Tensor_d0s2 * gradU = static_cast<Tensor_d0s2 *>(_gradU);

  (*gradU) = localGrad|invjacobian; // space contraction
  (*dgradUdqk) = localderivatives|invjacobian;
}

void
ShapeFunction::getNodesCoordinates(Node *nodes, Tensor *_nodescoordinates)
{
  Tensor_d1s0 & nodescoordinates = static_cast<Tensor_d1s0 &>(*_nodescoordinates);

  for (int j = 0; j < numdofs/3; j++) {
    Node &nd = nodes[j];  
    nodescoordinates[3*j] = nd.x;
    nodescoordinates[3*j+1] = nd.y;
    nodescoordinates[3*j+2] = nd.z;
  }
}

void
ShapeFunction::getDisplacements(Vector &disp, Tensor *_displacements)
{
  Tensor_d1s0 & displacements = static_cast<Tensor_d1s0 &>(*_displacements);
  disp.vectorToTensor(displacements);
}

void
ShapeFunction::getGlobalGrads(Tensor *_gradU, Tensor *_dgradUdqk, double *jac,
                              Tensor *_nodescoordinates, double xi[3], Tensor *_displacements,
                              Tensor *_localderivatives)
{
  // alternative version in which nodescoordinates and displacements are pre-computed, and localderivatives is pre-allocated
  Tensor_d1s2_sparse & localderivatives = static_cast<Tensor_d1s2_sparse &>(*_localderivatives);
  Tensor_d1s0 & nodescoordinates = static_cast<Tensor_d1s0 &>(*_nodescoordinates);
  Tensor_d0s2 localGrad;
  Tensor_d0s2 jacobian;
  Tensor_d0s2 invjacobian;

  // localDerivatives(i,j,k) = at node i, the derivative of j^{th} shape function w.r.t xi[k]
  getLocalDerivatives(&localderivatives, xi);

  //isoparametric elements
  jacobian = localderivatives%nodescoordinates; // dof contraction

  *jac = jacobian.getInverseAndDeterminant(invjacobian);
  
  Tensor_d1s0 & displacements = static_cast<Tensor_d1s0 &>(*_displacements);
  
  localGrad = localderivatives%displacements; // dof contraction

  Tensor_d1s2_sparse * dgradUdqk = static_cast<Tensor_d1s2_sparse *>(_dgradUdqk);
  Tensor_d0s2 * gradU = static_cast<Tensor_d0s2 *>(_gradU);

  (*gradU) = localGrad|invjacobian; // space contraction
  (*dgradUdqk) = localderivatives|invjacobian;
}

void
ShapeFunction::getJacobianDeterminant(double *jac, Node *nodes, double xi[3])
{
  Tensor_d1s2_sparse localderivatives(numdofs);
  Tensor_d1s0 nodescoordinates(numdofs);
  Tensor_d0s2 jacobian;

  // localDerivatives(i,j,k) = at node i, the derivative of j^{th} shape function w.r.t xi[k]
  getLocalDerivatives(&localderivatives, xi);

  for (int j = 0; j < numdofs/3; j++) {
    Node &nd = nodes[j];
    nodescoordinates[3*j] = nd.x;
    nodescoordinates[3*j+1] = nd.y;
    nodescoordinates[3*j+2] = nd.z;
  }

  //isoparametric elements
  jacobian = localderivatives%nodescoordinates; // dof contraction

  jacobian.getDeterminant(*jac);
}

void
ShapeFunction::getGradU(Tensor *_gradU, Node *nodes, double xi[3], Vector &disp)
{
  Tensor_d0s2 * gradU = static_cast<Tensor_d0s2 *>(_gradU);
  Tensor_d0s2 jacobian;
  Tensor_d0s2 invjacobian;
  Tensor_d1s0 nodescoordinates(numdofs);
  Tensor_d1s2_sparse localderivatives(numdofs);

  getLocalDerivatives(&localderivatives, xi);

  for (int i = 0; i < numdofs; i++) {
    int j = i/3;
    Node &nd = nodes[j]; 
    nodescoordinates[3*j] = nd.x;
    nodescoordinates[3*j+1] = nd.y;
    nodescoordinates[3*j+2] = nd.z;
  }

  //isoparametric elements
  jacobian = localderivatives%nodescoordinates;//dof contraction

  jacobian.getInverse(invjacobian);

  Tensor_d0s2 localGrad;
  Tensor_d1s0 displacements(numdofs);
  // Can be done differently XFL
  disp.vectorToTensor(displacements);

  localGrad = (localderivatives%displacements);//dof contraction

  (*gradU) = (localGrad|invjacobian);//space contraction
}

void
ShapeFunction::getGradU(Tensor *_gradU, double *jac, Node *nodes, double xi[3], Vector &disp)
{
  Tensor_d0s2 * gradU = static_cast<Tensor_d0s2 *>(_gradU);
  Tensor_d0s2 jacobian;
  Tensor_d0s2 invjacobian;
  Tensor_d1s0 nodescoordinates(numdofs);
  Tensor_d1s2_sparse localderivatives(numdofs);

  getLocalDerivatives(&localderivatives, xi);

  for (int i = 0; i < numdofs; i++) {
    int j = i/3;
    Node &nd = nodes[j];
    nodescoordinates[3*j] = nd.x;
    nodescoordinates[3*j+1] = nd.y;
    nodescoordinates[3*j+2] = nd.z;
  }

  //isoparametric elements
  jacobian = localderivatives%nodescoordinates;//dof contraction

  jacobian.getInverse(invjacobian);

  Tensor_d0s2 localGrad;
  Tensor_d1s0 displacements(numdofs);
  // Can be done differently XFL
  disp.vectorToTensor(displacements);

  localGrad = (localderivatives%displacements);//dof contraction

  (*gradU) = (localGrad|invjacobian);//space contraction

  jacobian.getDeterminant(*jac);
}

double
ShapeFunction::interpolateScalar(double *_q, double _xi[3])
{
  std::cerr << "Error: ShapeFunction::interpolateScalar is not implemented\n";
  return 0;
}
