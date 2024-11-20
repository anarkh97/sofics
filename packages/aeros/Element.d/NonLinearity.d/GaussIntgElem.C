#include <cstdio>
#include <Utils.d/dbg_alloca.h>
#include <cmath>
#include <iostream>
#include <Element.d/NonLinearity.d/GaussIntgElem.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#include <Math.d/FullSquareMatrix.h>
#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif

FullSquareMatrix  
GaussIntgElement::stiffness(const CoordSet& cs, double *k, int) const
{
  int i;
  FullSquareMatrix kTan(numDofs(), k);

  Node *nodes = (Node *) dbg_alloca(numNodes()*sizeof(Node));
  int nnd = numNodes();
  int *ndn = (int *)dbg_alloca(nnd*sizeof(int)); 
  this->nodes(ndn);
  for(i = 0; i < nnd; ++i)
    nodes[i] = *cs[ndn[i]];
  int ndofs = numDofs();
  double *disp = (double *) dbg_alloca(ndofs*sizeof(double));
  for(i = 0; i < ndofs; ++i)
    disp[i] = 0;

  ShapeFunction *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  // Obtain the storage for gradU ( 3x3 )
  Tensor &gradU = *shapeF->getGradUInstance();

  // Obtain the storage for dgradUdqk ( ndof x3x3 )
  Tensor &dgradUdqk = *shapeF->getDgradUDqkInstance();
 
  // NDofsx3x3x-> 6xNDofs
  Tensor &B = *strainEvaluator->getBInstance(ndofs);

  // NdofsxNdofsx3x3x -> 6xNdofsxNdofs
  Tensor &DB = *strainEvaluator->getDBInstance(ndofs);

  Tensor &D = *strainEvaluator->getTMInstance();
  Tensor &e = *strainEvaluator->getStrainInstance();
  Tensor &s = *strainEvaluator->getStressInstance();
  Tensor &temp1 = *strainEvaluator->getBInstance(ndofs);
  Tensor_d2s0 temp2(ndofs,false);
  Tensor_d2s0 temp3(ndofs,false);

  // Obtain the storage for cached quantities (if any)
  Tensor *cache = strainEvaluator->getCacheInstance();
  
  //fprintf(stderr,"Je suis dans stiffness\n");

  int ngp = getNumGaussPoints();

  kTan.zero();
  for(i = 0; i < ngp; i++) {
    double point[3], weight, jac;
    getGaussPointAndWeight(i, point, weight);

    StackVector dispVec(disp,ndofs);
    shapeF->getGlobalGrads(&gradU, &dgradUdqk, &jac, nodes, point, dispVec);
    strainEvaluator->getEBandDB(e, B, DB, gradU, dgradUdqk, cache);

    //material->getStress(&s,e,0);          
    //material->getTangentMaterial(&D, e, 0);
    //material->getStressAndTangentMaterial(&s, &D, e, 0);

    material->integrate(&s, &D, e, e, 0, 0, 0, cache);

    temp1 =  D || B;
    temp2 =   B||temp1;
    temp3 = DB || s;
    temp3 = temp3 + temp2;
    temp3 = (weight * fabs(jac))*temp3;

    kTan += temp3;
  }
  delete &temp1;
  delete &s;
  delete &gradU;
  delete &dgradUdqk;
  delete &B;
  delete &DB; 
  delete &e;
  delete &D;
  if(cache) delete cache;
  return kTan;
}

FullSquareMatrix  
GaussIntgElement::massMatrix(const CoordSet& cs, double* k, int) const
{
  FullSquareMatrix m(numDofs(), k);
  std::cerr << "GaussIntgElement::massMatrix not implemented\n"; exit(-1);
  return m;
}
 
void
GaussIntgElement::getStiffAndForce(Node *nodes, double *disp,
                                   double *state, FullSquareMatrix &kTan,
                                   double *force)
{
  int ndofs = numDofs();
  ShapeFunction *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  // Obtain the storage for gradU ( 3x3 )
  Tensor &gradU = *shapeF->getGradUInstance();

  // Obtain the storage for dgradUdqk ( ndof x3x3 )
  Tensor &dgradUdqk = *shapeF->getDgradUDqkInstance();
  
  // NDofsx3x3x-> 6xNDofs
  Tensor &B = *strainEvaluator->getBInstance(ndofs);

  // NdofsxNdofsx3x3x -> 6xNdofsxNdofs but diagonal in dofs
  Tensor &DB = *strainEvaluator->getDBInstance(ndofs);

  Tensor &D = *strainEvaluator->getTMInstance();
  Tensor &e = *strainEvaluator->getStrainInstance();
  Tensor &s = *strainEvaluator->getStressInstance();

  Tensor_d1s0 nodeforce(ndofs);
  Tensor_d1s0 temp0(ndofs);
  Tensor &temp1 = *strainEvaluator->getBInstance(ndofs);
  Tensor_d2s0 temp2(ndofs,false);
  Tensor_d2s0 temp3(ndofs,false);

  // Obtain the storage for cached quantities (if any)
  Tensor *cache = strainEvaluator->getCacheInstance();

  //fprintf(stderr,"Je suis dans getStiffAndForce\n");

  int i,j;
  int ngp = getNumGaussPoints();

  kTan.zero();
  
  for(i = 0; i < ngp; i++) {
    double point[3], weight, jac;
    getGaussPointAndWeight(i, point, weight);

    StackVector dispVec(disp,ndofs);
        
    shapeF->getGlobalGrads(&gradU, &dgradUdqk,  &jac, nodes, point, dispVec);
    strainEvaluator->getEBandDB(e, B, DB, gradU, dgradUdqk, cache);

    //material->getStress(&s, e, 0);        
    //material->getStressAndTangentMaterial(&s, &D, e, 0);
        
    material->integrate(&s, &D, e, e, 0, 0, 0, cache);

    temp0 = s || B;
    temp0 = (weight*jac)*temp0;
    nodeforce = nodeforce + temp0;
    temp1 =  D || B;
    temp2 =   B||temp1;
    temp3 = DB || s;
    temp3 = temp3 + temp2;
    temp3 = (weight * fabs(jac))*temp3;
    kTan += temp3;
    temp3 = DB || s;
  }

  for(j = 0; j < ndofs; ++j) {
    force[j] = - nodeforce[j];
  }

  ////////////////////////////////////////////////////////
  //Check convergence with approximate tangent stiffness//
  ////////////////////////////////////////////////////////

  Tensor &eleft = *strainEvaluator->getStrainInstance();
  Tensor &eright = *strainEvaluator->getStrainInstance();
  Tensor &DBleft = *strainEvaluator->getDBInstance(ndofs);
  Tensor &DBright = *strainEvaluator->getDBInstance(ndofs);
  Tensor &gradUleft = *shapeF->getGradUInstance();
  Tensor &gradUright = *shapeF->getGradUInstance();
  Tensor &dgradUdqkleft = *shapeF->getDgradUDqkInstance();
  Tensor &dgradUdqkright = *shapeF->getDgradUDqkInstance();
  Tensor &sleft = *strainEvaluator->getStressInstance();
  Tensor &sright = *strainEvaluator->getStressInstance();
  Tensor &Bleft = *strainEvaluator->getBInstance(ndofs);
  Tensor &Bright = *strainEvaluator->getBInstance(ndofs);
  Tensor *cacheleft = strainEvaluator->getCacheInstance();
  Tensor *cacheright = strainEvaluator->getCacheInstance();

  double *kt = new double[ndofs*ndofs];
  FullSquareMatrix kTant(ndofs, kt);
 
  int a,b;
  double dispsqNorm =0.0;
  for (b=0; b < ndofs; ++b)
    dispsqNorm+=disp[b]*disp[b];
  double epsilon = 1e-8; //(1e-8)*sqrt(dispsqNorm);
  //if(dispsqNorm == 0.0) epsilon = 1e-8;
  kTant.zero();
  fprintf(stderr, "epsilon = %e\n" , epsilon);
  double* dispLeft = new double[ndofs];
  double* dispRight = new double[ndofs];
  for (b=0; b < ndofs; ++b) {
    dispLeft[b] = disp[b];
    dispRight[b] = disp[b]; 
  }

  for (b=0; b < ndofs; ++b) {

    if(disp[b] > 1 || disp[b] < -1)
      epsilon = 1e-8*disp[b];
    else
      epsilon = 1e-8;
    dispLeft[b] = disp[b]+epsilon;
    dispRight[b] = disp[b]-epsilon;
    Tensor_d1s0 nodeforceLeft(ndofs);
    Tensor_d1s0 nodeforceRight(ndofs);
    StackVector dispVecTestLeft(dispLeft,ndofs);
    StackVector dispVecTestRight(dispRight,ndofs);
     
    for(i = 0; i < ngp; i++) {
      double point[3], weight, jac;
      getGaussPointAndWeight(i, point, weight);

      shapeF->getGlobalGrads(&gradUleft, &dgradUdqkleft,  &jac, nodes, point, dispVecTestLeft);
      strainEvaluator->getEBandDB(eleft, Bleft, DBleft, gradUleft, dgradUdqkleft, cacheleft);
      material->integrate(&sleft, &D,  eleft, eleft, 0, 0, 0, cacheleft);
      temp0 = sleft || Bleft;
      temp0 = (weight*jac)*temp0;
      nodeforceLeft = nodeforceLeft + temp0;

      shapeF->getGlobalGrads(&gradUright, &dgradUdqkright,  &jac, nodes, point, dispVecTestRight);
      strainEvaluator->getEBandDB(eright, Bright, DBright, gradUright, dgradUdqkright, cacheright);
      material->integrate(&sright, &D,  eright, eright, 0, 0, 0, cacheright);
      temp0 = sright || Bright;
      temp0 = (weight*jac)*temp0;
      nodeforceRight = nodeforceRight + temp0;
    }

    for (a=0; a < ndofs; ++a) {              
      kTant[a][b] = (nodeforceLeft[a] - nodeforceRight[a])/(2*epsilon);
    }

    /* if(b == 1) fprintf(stderr, "%e vs %e or %e or %e\n", kTan[b][b], kTant[b][b],
                          (-nodeforceLeft[b]-force[b])/epsilon, (force[b]+nodeforceRight[b])/epsilon);*/

    dispLeft[b] = disp[b];
    dispRight[b] = disp[b];              
  }
  //double diffcol[ndofs];
  //double col[ndofs];
  double diffsqNorm = 0.0;
  double sqNorm = 0.0;
  double relativeNorm;

  for (a=0; a < ndofs; ++a)         
    for (b=0; b < ndofs; ++b) {
      //diffcol[a]+=fabs(kTan[a][b]-kTant[a][b]);
      //col[a]+=fabs(kTant[a][b])

      diffsqNorm+=(kTan[a][b]-kTant[a][b])*(kTan[a][b]-kTant[a][b]);
      sqNorm+=kTant[a][b]*kTant[a][b];
    }
  relativeNorm = sqrt(diffsqNorm/sqNorm);
  fprintf(stderr, "Relative Norm = %e\n", relativeNorm);
        
  delete &Bleft;
  delete &Bright;        
  delete &sleft;
  delete &sright;
  delete &gradUleft;
  delete &gradUright;
  delete &dgradUdqkleft; 
  delete &dgradUdqkright;
  delete &eleft;
  delete &eright;
  delete &DBleft;
  delete &DBright;
  if(cacheleft) delete cacheleft;
  if(cacheright) delete cacheright;

  /////////////////////////////////////////////////////////////////////////////

  delete &temp1;
  delete &gradU;
  delete &dgradUdqk;
  delete &B;
  delete &DB;
  delete &e;
  delete &s;
  delete &D;
  if(cache) delete cache;
}

void 
GaussIntgElement::integrate(Node *nodes, double *dispn, double *staten,
                            double *dispnp, double *statenp,
                            FullSquareMatrix &kTan, double *force,
                            double dt, double *temps)
{
  int ndofs = numDofs();
  ShapeFunction *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  // Obtain the storage for gradU ( 3x3 )
  Tensor &gradUn = *shapeF->getGradUInstance();
  Tensor &gradUnp = *shapeF->getGradUInstance();
  // Obtain the storage for dgradUdqk ( ndof x3x3 )
  Tensor &dgradUdqknp = *shapeF->getDgradUDqkInstance();

  // NDofsx3x3x-> 6xNDofs
  Tensor &Bnp = *strainEvaluator->getBInstance(ndofs);

  // NdofsxNdofsx3x3x -> 6xNdofsxNdofs but sparse
  Tensor &DBnp = *strainEvaluator->getDBInstance(ndofs);

  Tensor &Dnp = *strainEvaluator->getTMInstance();
  Tensor &en = *strainEvaluator->getStrainInstance();
  Tensor &enp = *strainEvaluator->getStrainInstance();
  Tensor &s = *strainEvaluator->getStressInstance();
  
  Tensor_d1s0 nodeforce(ndofs);
  Tensor_d1s0 temp0(ndofs);
  Tensor &temp1 = *strainEvaluator->getBInstance(ndofs);
  Tensor_d2s0 temp2(ndofs, false);
  Tensor_d2s0 temp3(ndofs, false);

  // Obtain storage for cached quantities (if any)
  //Tensor *cachen = strainEvaluator->getCacheInstance();
  Tensor *cachenp = strainEvaluator->getCacheInstance();

  int i,j;
  int ngp = getNumGaussPoints();
  int nstatepgp = material->getNumStates();
  
  kTan.zero();

  //fprintf(stderr,"Je suis dans integrate\n");

  Tensor &nodescoordinates = *shapeF->getNodesCoordinatesInstance();
  shapeF->getNodesCoordinates(nodes, &nodescoordinates);

  Tensor &displacements = *shapeF->getDisplacementsInstance();
  StackVector dispVecnp(dispnp, ndofs);
  StackVector dispVecn(dispn, ndofs);
  shapeF->getDisplacements(dispVecnp, &displacements);

  Tensor &localderivatives = *shapeF->getLocalDerivativesInstance();

  for(i = 0; i < ngp; i++) {

    double point[3], weight, jacn, jacnp, tempnp;
 
    getGaussPointAndWeight(i, point, weight);

    shapeF->getGradU(&gradUn, nodes, point, dispVecn);
    shapeF->getGlobalGrads(&gradUnp, &dgradUdqknp, &jacnp, &nodescoordinates, point, &displacements, &localderivatives);

    //strainEvaluator->getE(en, gradUn, cachen);
    strainEvaluator->getEBandDB(enp, Bnp, DBnp, gradUnp, dgradUdqknp, cachenp, staten + nstatepgp*i);
    strainEvaluator->getL(gradUnp, gradUn, cachenp, dt); // adds velocity gradient tensor to cache for GreenLag strain type

    // compute temperature at integration point
    tempnp = (temps) ? shapeF->interpolateScalar(temps, point) : 0;

    material->integrate(&s, &Dnp, en, enp,
                        staten + nstatepgp*i, statenp + nstatepgp*i, tempnp, cachenp, dt);
    temp0 = s || Bnp;
    temp0 = (weight*jacnp)*temp0;
    nodeforce = nodeforce + temp0;
    temp1 =  Dnp || Bnp;
    temp2 =   Bnp||temp1;
    temp3 = DBnp || s;
    temp3 = temp3 + temp2;
    if(jacnp < 0) std::cerr << "warning: jacnp < 0\n";
    temp3 = (weight * fabs(jacnp))*temp3;
    kTan += temp3;

  }

#ifdef USE_EIGEN3
  if(material->getPosdefifyTol() >= 0) {
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > K(kTan.data(),ndofs,ndofs);
    /*if(!K.isApprox(0.5*(K+K.transpose()))) {
      std::cerr << "Matrix K is not symmetric\n";
    }*/
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > es(K);
    if(es.eigenvalues().minCoeff() < -Eigen::NumTraits<double>::dummy_precision()) {
      Eigen::VectorXd d(ndofs);
      for(j=0; j<ndofs; ++j) 
        if(es.eigenvalues()[j] > material->getPosdefifyTol()) d[j] = es.eigenvalues()[j]; 
        else d[j] = 0;
      K = es.eigenvectors()*d.asDiagonal()*es.eigenvectors().transpose();
    }
  }
#endif

  for(j = 0; j < ndofs; ++j) {
    force[j] = - nodeforce[j];}

  delete &temp1;
  delete &gradUn;
  delete &gradUnp;
  delete &dgradUdqknp;
  delete &Bnp;
  delete &DBnp;
  delete &en;
  delete &enp;
  delete &s;
  delete &Dnp;
  //if(cachen) delete cachen;
  if(cachenp) delete cachenp;
  delete &nodescoordinates;
  delete &displacements;
  delete &localderivatives;
}

void 
GaussIntgElement::integrate(Node *nodes, double *dispn, double *staten,
                            double *dispnp, double *statenp,
                            double *force, double dt, double *temps)
{
  int ndofs = numDofs();
  ShapeFunction *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  // Obtain the storage for gradU ( 3x3 )
  Tensor &gradUn = *shapeF->getGradUInstance();
  Tensor &gradUnp = *shapeF->getGradUInstance();
  // Obtain the storage for dgradUdqk ( ndof x3x3 )
  Tensor &dgradUdqknp = *shapeF->getDgradUDqkInstance();

  // NDofsx3x3x-> 6xNDofs
  Tensor &Bnp = *strainEvaluator->getBInstance(ndofs);

  Tensor &en = *strainEvaluator->getStrainInstance();
  Tensor &enp = *strainEvaluator->getStrainInstance();
  Tensor &s = *strainEvaluator->getStressInstance();
  
  Tensor_d1s0 nodeforce(ndofs);
  Tensor_d1s0 temp0(ndofs);

  // Obtain storage for cached quantities (if any)
  //Tensor *cachen = strainEvaluator->getCacheInstance();
  Tensor *cachenp = strainEvaluator->getCacheInstance();

  int i,j;
  int ngp = getNumGaussPoints();
  int nstatepgp = material->getNumStates();
  
  //fprintf(stderr,"Je suis dans integrate\n");

  Tensor &nodescoordinates = *shapeF->getNodesCoordinatesInstance();
  shapeF->getNodesCoordinates(nodes, &nodescoordinates);

  Tensor &displacements = *shapeF->getDisplacementsInstance();
  StackVector dispVecnp(dispnp, ndofs);
  StackVector dispVecn(dispn, ndofs);
  shapeF->getDisplacements(dispVecnp, &displacements);

  Tensor &localderivatives = *shapeF->getLocalDerivativesInstance();
  
  for(i = 0; i < ngp; i++) {

    double point[3], weight, jacn, jacnp, tempnp;
 
    getGaussPointAndWeight(i, point, weight);

    shapeF->getGradU(&gradUn, nodes, point, dispVecn);
    shapeF->getGlobalGrads(&gradUnp, &dgradUdqknp, &jacnp, &nodescoordinates, point, &displacements, &localderivatives);

    //strainEvaluator->getE(en, gradUn, cachen);
    strainEvaluator->getEandB(enp, Bnp, gradUnp, dgradUdqknp, cachenp, staten + nstatepgp*i);
    strainEvaluator->getL(gradUnp, gradUn, cachenp, dt); // adds velocity gradient tensor to cache

    // compute temperature at integration point
    tempnp = (temps) ? shapeF->interpolateScalar(temps, point) : 0;

    material->integrate(&s, en, enp,
                        staten + nstatepgp*i, statenp + nstatepgp*i, tempnp, cachenp, dt);

    temp0 = s || Bnp;
    temp0 = (weight*jacnp)*temp0;
    nodeforce = nodeforce + temp0;
  }

  for(j = 0; j < ndofs; ++j) {
    force[j] = - nodeforce[j];}

  delete &gradUn;
  delete &gradUnp;
  delete &dgradUdqknp;
  delete &Bnp;
  delete &en;
  delete &enp;
  delete &s;
  //if(cachen) delete cachen;
  if(cachenp) delete cachenp;
  delete &nodescoordinates;
  delete &displacements;
  delete &localderivatives;
}

void
GaussIntgElement::initStates(double *st)
{
  NLMaterial *material = getMaterial();
  int ninterns = material->getNumStates();
  int ngp = getNumGaussPoints();

  for(int i = 0; i < ngp; ++i)
    material->initStates(st+i*ninterns);
}

void
GaussIntgElement::updateStates(Node *nodes, double *state, double *dispn, double *dispnp, double *temps, double dt)
{
  int ndofs = numDofs();
  ShapeFunction *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  // Obtain the storage for gradU ( 3x3 )
  Tensor &gradUn = *shapeF->getGradUInstance();
  Tensor &gradUnp = *shapeF->getGradUInstance();

  Tensor &Dnp = *strainEvaluator->getTMInstance();
  Tensor &en = *strainEvaluator->getStrainInstance();
  Tensor &enp = *strainEvaluator->getStrainInstance();
  Tensor &s = *strainEvaluator->getStressInstance();

  StackVector dispVecn(dispn, ndofs);
  StackVector dispVecnp(dispnp, ndofs);

  // Obtain the storage for cached quantities (if any)
  Tensor *cache = strainEvaluator->getCacheInstance();

  int nstatepgp = material->getNumStates();

  for(int i = 0; i < getNumGaussPoints(); i++) {

    double point[3], weight, jacn, jacnp, tempnp;

    getGaussPointAndWeight(i, point, weight);

    shapeF->getGradU(&gradUn, nodes, point, dispVecn);
    shapeF->getGradU(&gradUnp, nodes, point, dispVecnp);
    //strainEvaluator->getE(en, gradUn);
    strainEvaluator->getE(enp, gradUnp, cache, state + nstatepgp*i);
    strainEvaluator->getL(gradUnp, gradUn, cache, dt); // adds velocity gradient tensor to cache

    //material->updateStates(en, enp, state + nstatepgp*i);
    //material->getStress(&s, e, 0);
    //material->getStressAndTangentMaterial(&s, &D, enp, 0);

    // compute temperature at integration point
    tempnp = (temps) ? shapeF->interpolateScalar(temps, point) : 0;

    double *state_copy = new double[nstatepgp];
    for(int j = 0; j < nstatepgp; ++j) state_copy[j] = state[nstatepgp*i+j];
    
    material->integrate(&s, &Dnp, en, enp,
                        state_copy, state + nstatepgp*i, tempnp, cache, dt);

    delete [] state_copy;
  }

  delete &gradUn;
  delete &gradUnp;
  delete &en;
  delete &enp;
  delete &s;
  delete &Dnp;
  if(cache) delete cache;
}

void
copyTens(Tensor *stens, double *svec)
{
  Tensor_d0s2_Ss12 *symmetricTensor = dynamic_cast<Tensor_d0s2_Ss12 *>(stens);
  if(symmetricTensor) {
    for(int j = 0; j < 3; ++j)
      for(int k = 0; k < 3; ++k)
        svec[3*j+k] = (*symmetricTensor)(j,k);
  }
  else {
    Tensor_d0s2 *nonsymmetricTensor = dynamic_cast<Tensor_d0s2 *>(stens);
    if(nonsymmetricTensor) {
      for(int j = 0; j < 3; ++j)
        for(int k = 0; k < 3; ++k)
          svec[3*j+k] = (*nonsymmetricTensor)(j,k);
    }
    else { std::cerr << "ERROR: unrecognized tensor type"; exit(-1); }
  }
}

#ifdef USE_EIGEN3
void
copyTens(Tensor *stens, Eigen::Matrix3d &smat)
{
  Tensor_d0s2_Ss12 *symmetricTensor = dynamic_cast<Tensor_d0s2_Ss12 *>(stens);
  if(symmetricTensor) {
    for(int j = 0; j < 3; ++j)
      for(int k = 0; k < 3; ++k)
        smat(j,k) = (*symmetricTensor)(j,k);
  }
  else {
    Tensor_d0s2 *nonsymmetricTensor = dynamic_cast<Tensor_d0s2 *>(stens);
    if(nonsymmetricTensor) {
      for(int j = 0; j < 3; ++j)
        for(int k = 0; k < 3; ++k)
          smat(j,k) = (*nonsymmetricTensor)(j,k);
    }
    else { std::cerr << "ERROR: unrecognized tensor type"; exit(-1); }
  }
}
#endif

void
GaussIntgElement::getStrainTens(Node *nodes, double *dispnp, double (*result)[9], int avgnum)
{
  // compute strain tensor at element nodes for post processing
  int ndofs = numDofs();
  ShapeFunction *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the storage for gradU ( 3x3 )
  Tensor &gradUnp = *shapeF->getGradUInstance();

  Tensor &enp = *strainEvaluator->getStrainInstance();

  // Obtain the storage for cached quantities (if any)
  Tensor *cache = strainEvaluator->getCacheInstance();

  int numPoints = (avgnum == -1) ? getNumGaussPoints() : numNodes();
  for(int i = 0; i < numPoints; ++i) {

    StackVector dispVecnp(dispnp, ndofs);

    double point[3], weight;
    if(avgnum == -1) getGaussPointAndWeight(i, point, weight);
    else getLocalNodalCoords(i, point);

    shapeF->getGradU(&gradUnp, nodes, point, dispVecnp);

    strainEvaluator->getE(enp, gradUnp, cache);

    copyTens(&enp, result[i]);
  }
 
  delete &gradUnp;
  delete &enp;
  if(cache) delete cache;
}

void
GaussIntgElement::getVonMisesStrain(Node *nodes, double *dispnp, double *result, int avgnum)
{
  // compute von mises strain at element nodes for post processing
  int ndofs = numDofs();
  ShapeFunction *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the storage for gradU ( 3x3 )
  Tensor &gradUnp = *shapeF->getGradUInstance();

  Tensor &enp = *strainEvaluator->getStrainInstance();

  // Obtain the storage for cached quantities (if any)
  Tensor *cache = strainEvaluator->getCacheInstance();

  int numPoints = (avgnum == -1) ? getNumGaussPoints() : numNodes();
  for(int i = 0; i < numPoints; ++i) {

    StackVector dispVecnp(dispnp, ndofs);

    double point[3] = { 0, 0, 0 }, weight;
    if(avgnum == -1) getGaussPointAndWeight(i, point, weight);
    else getLocalNodalCoords(i, point);

    shapeF->getGradU(&gradUnp, nodes, point, dispVecnp);

    strainEvaluator->getE(enp, gradUnp, cache);

#ifdef USE_EIGEN3
    Eigen::Matrix3d M;
    copyTens(&enp,M);
    // compute the deviatoric stress/strain tensor and it's second invariant
    Eigen::Matrix3d dev = M - (M.trace()/3)*Eigen::Matrix3d::Identity();
    double J2 = 0.5*(dev*dev).trace();
    // compute the effective stress/strain
    result[i] = sqrt(3*J2);
#else
    result[i] = 0;
#endif
  }

  delete &gradUnp;
  delete &enp;
  if(cache) delete cache;
}

void
GaussIntgElement::getStressTens(Node *nodes, double *dispn, double *staten,
                                double *dispnp, double *statenp, double (*result)[9],
                                int avgnum, double *temps)
{
  // stress tensor at element nodes for post processing

  int ndofs = numDofs();
  ShapeFunction *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  // Obtain the storage for gradU ( 3x3 )
  Tensor &gradUnp = *shapeF->getGradUInstance();

  Tensor &enp = *strainEvaluator->getStrainInstance();
  Tensor &s = *strainEvaluator->getStressInstance();

  // Obtain the storage for cached quantities (if any)
  Tensor *cache = strainEvaluator->getCacheInstance();

  Tensor_d0s2_Ss12 S; // second Piola-Kirchhoff stress tensor

  int nstatepgp = material->getNumStates();
  if(nstatepgp == 0 && avgnum != -1) { // evaluate the stress at the nodes
    for(int i = 0; i < numNodes(); ++i) {

      StackVector dispVecnp(dispnp, ndofs);

      double point[3];
      getLocalNodalCoords(i, point);

      shapeF->getGradU(&gradUnp, nodes, point, dispVecnp);
      strainEvaluator->getE(enp, gradUnp, cache, statenp);

      material->getStress(&s, enp, statenp, (temps ? temps[i] : 0));
      strainEvaluator->transformStress(s, cache, S);
      copyTens(&S, result[i]);
    }
  }
  else { // evaluate at the gauss points and (if avgnum != -1) extrapolate to the nodes

    double (*gpstress)[9] = new double[getNumGaussPoints()][9];
    //Tensor &gradUn = *shapeF->getGradUInstance();
    Tensor &en = *strainEvaluator->getStrainInstance();

    for(int i = 0; i < getNumGaussPoints(); i++) {

      double point[3], weight, tempnp;
      StackVector dispVecnp(dispnp, ndofs);

      getGaussPointAndWeight(i, point, weight);

      shapeF->getGradU(&gradUnp, nodes, point, dispVecnp);
      //strainEvaluator->getE(en, gradUn);
      strainEvaluator->getE(enp, gradUnp, cache, staten + nstatepgp*i);

      // compute temperature at integration point
      tempnp = (temps) ? shapeF->interpolateScalar(temps, point) : 0;

      material->integrate(&s, en, enp,
                          staten + nstatepgp*i, statenp + nstatepgp*i, tempnp, cache);

      strainEvaluator->transformStress(s, cache, S);

      if(avgnum == -1) copyTens(&S, result[i]);
      else copyTens(&S, gpstress[i]);
    }

    if(avgnum != -1) {
      //TODO extrapolate from gauss points to nodes
      //temporary fix implemented here is to return the average of all the gauss points for each node
      double average[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
      for(int i = 0; i < getNumGaussPoints(); i++)
        for(int j = 0; j < 9; ++j)
          average[j] += gpstress[i][j];
      for(int i = 0; i < 9; ++i)
        average[i] /= getNumGaussPoints();
      for(int i = 0; i < numNodes(); ++i)
        for(int j = 0; j < 9; ++j)
          result[i][j] = average[j];
    }

    delete [] gpstress;
    //delete &gradUn;
    delete &en;
  }

  delete &gradUnp;
  delete &enp;
  delete &s;
  if(cache) delete cache;
}

void
GaussIntgElement::getVonMisesStress(Node *nodes, double *dispn, double *staten,
                                    double *dispnp, double *statenp, double *result,
                                    int avgnum, double *temps)
{
  // von mises stress at element nodes for post processing

  int ndofs = numDofs();
  ShapeFunction *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  // Obtain the storage for gradU ( 3x3 )
  Tensor &gradUnp = *shapeF->getGradUInstance();

  Tensor &enp = *strainEvaluator->getStrainInstance();
  Tensor &s = *strainEvaluator->getStressInstance();

  // Obtain the storage for cached quantities (if any)
  Tensor *cache = strainEvaluator->getCacheInstance();

  Tensor_d0s2_Ss12 S; // second Piola-Kirchhoff stress tensor

  int nstatepgp = material->getNumStates();

  if(nstatepgp == 0 && avgnum != -1) { // evaluate the stress at the nodes
    for(int i = 0; i < numNodes(); ++i) {

      StackVector dispVecnp(dispnp, ndofs);

      double point[3];
      getLocalNodalCoords(i, point);

      shapeF->getGradU(&gradUnp, nodes, point, dispVecnp);
      strainEvaluator->getE(enp, gradUnp, cache, statenp);

      material->getStress(&s, enp, statenp, (temps ? temps[i] : 0));
      strainEvaluator->transformStress(s, cache, S);
#ifdef USE_EIGEN3
      Eigen::Matrix3d M;
      copyTens(&S,M);
      // compute the deviatoric stress/strain tensor and it's second invariant
      Eigen::Matrix3d dev = M - (M.trace()/3)*Eigen::Matrix3d::Identity();
      double J2 = 0.5*(dev*dev).trace();
      // compute the effective stress/strain
      result[i] = sqrt(3*J2);
#else
      result[i] = 0;
#endif
    }
  }
  else { // evaluate at the gauss points and (if avgnum != -1) extrapolate to the nodes

    double *gpstress = new double[getNumGaussPoints()];
    //Tensor &gradUn = *shapeF->getGradUInstance();
    Tensor &en = *strainEvaluator->getStrainInstance();

    for(int i = 0; i < getNumGaussPoints(); i++) {

      double point[3], weight, tempnp;
      //StackVector dispVecn(dispn, ndofs);
      StackVector dispVecnp(dispnp, ndofs);

      getGaussPointAndWeight(i, point, weight);

      //shapeF->getGradU(&gradUn, nodes, point, dispVecn);
      shapeF->getGradU(&gradUnp, nodes, point, dispVecnp);
      //strainEvaluator->getE(en, gradUn);
      strainEvaluator->getE(enp, gradUnp, cache, staten + nstatepgp*i);

      // compute temperature at integration point
      tempnp = (temps) ? shapeF->interpolateScalar(temps, point) : 0;

      material->integrate(&s, en, enp,
                          staten + nstatepgp*i, statenp + nstatepgp*i, tempnp, cache);

      strainEvaluator->transformStress(s, cache, S);
#ifdef USE_EIGEN3
      Eigen::Matrix3d M;
      copyTens(&S,M);
      // compute the deviatoric stress/strain tensor and it's second invariant
      Eigen::Matrix3d dev = M - (M.trace()/3)*Eigen::Matrix3d::Identity();
      double J2 = 0.5*(dev*dev).trace();
      // compute the effective stress/strain
      if(avgnum == -1) result[i] = sqrt(3*J2);
      else gpstress[i] = sqrt(3*J2);
#else
      if(avgnum == -1) result[i] = 0;
      else gpstress[i] = 0;
#endif
    }

    if(avgnum != -1) {
      //TODO extrapolate from gauss points to nodes
      //temporary fix implemented here is to return the average of all the gauss points for each node
      double average = 0;
      for(int i = 0; i < getNumGaussPoints(); i++)
        average += gpstress[i];
      average /= getNumGaussPoints();
      for(int i = 0; i < numNodes(); ++i)
        result[i] = average;
    }

    delete [] gpstress;
    //delete &gradUn;
    delete &en;
  }

  delete &gradUnp;
  delete &enp;
  delete &s;
  if(cache) delete cache;
}

void
GaussIntgElement::getEquivPlasticStrain(double *statenp, double *result, int avgnum)
{
  // Obtain the material model
  NLMaterial *material = getMaterial();
  if(material->getNumStates() == 0) {
    int numPoints = (avgnum == -1) ? getNumGaussPoints() : numNodes();
    for(int i = 0; i < numPoints; ++i) {
      result[i] = 0;
    }
    return;
  }

  // TODO: replace averaging below with extrapolation using shape functions
  double average = 0;
  for(int i = 0; i < getNumGaussPoints(); i++) {
     double eps = material->getEquivPlasticStrain(statenp + i*material->getNumStates());
     if(avgnum == -1) result[i] = eps;
     else average += eps;
  }

  if(avgnum != -1) {
    average /= getNumGaussPoints();

    for(int i = 0; i <  numNodes(); ++i) {
      result[i] = average;
    }
  }
}

void
GaussIntgElement::getBackStressTens(double *statenp, double (*result)[9], int avgnum)
{
  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  // Storage for backstress
  Tensor &s = *strainEvaluator->getStressInstance();

  if(material->getNumStates() == 0) {
    int numPoints = (avgnum == -1) ? getNumGaussPoints() : numNodes();
    for(int i = 0; i < numPoints; ++i) {
      for(int j = 0; j < 9; ++j) {
        result[i][j] = 0;
      }
    }
    return;
  }

  // TODO: replace averaging below with extrapolation using shape functions
  double average[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  double gpbackstress[9];
  for(int i = 0; i < getNumGaussPoints(); i++) {
    if(material->getBackStress(statenp + i*material->getNumStates(), &s)) {
      if(avgnum == -1) copyTens(&s, result[i]);
      else {
        copyTens(&s, gpbackstress);
        for(int j = 0; j < 9; ++j) {
          average[j] += gpbackstress[j];
        }
      }
    }
    else {
      if(avgnum == -1) {
        for(int j = 0; j < 9; ++j) result[i][j] = 0;
      }
    }
  }

  if(avgnum != -1) {
    for(int j = 0; j < 9; ++j) {
      average[j] /= getNumGaussPoints();
    }

    for(int i = 0; i <  numNodes(); ++i) {
      for(int j = 0; j < 9; ++j) {
        result[i][j] = average[j];
      }
    }
  }
}

void
GaussIntgElement::getPlasticStrainTens(double *statenp, double (*result)[9], int avgnum)
{
  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  // Storage for plastic strain
  Tensor &e = *strainEvaluator->getStrainInstance();

  if(material->getNumStates() == 0) {
    int numPoints = (avgnum == -1) ? getNumGaussPoints() : numNodes();
    for(int i = 0; i < numPoints; ++i) {
      for(int j = 0; j < 9; ++j) {
        result[i][j] = 0;
      }
    }
    return;
  }

  // TODO: replace averaging below with extrapolation using shape functions
  double average[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  double gpplasticstrain[9];
  for(int i = 0; i < getNumGaussPoints(); i++) {
    if(material->getPlasticStrain(statenp + i*material->getNumStates(), &e)) {
      if(avgnum == -1) copyTens(&e, result[i]);
      else {
        copyTens(&e, gpplasticstrain);
        for(int j = 0; j < 9; ++j) {
          average[j] += gpplasticstrain[j];
        }
      }
    }
    else {
      if(avgnum == -1) {
        for(int j = 0; j < 9; ++j) result[i][j] = 0;
      }
    }
  }

  if(avgnum != -1) {
    for(int j = 0; j < 9; ++j) {
      average[j] /= getNumGaussPoints();
    }

    for(int i = 0; i <  numNodes(); ++i) {
      for(int j = 0; j < 9; ++j) {
        result[i][j] = average[j];
      }
    }
  }
}

void
GaussIntgElement::getDamage(double *statenp, double *result, int avgnum)
{
  // Obtain the material model
  NLMaterial *material = getMaterial();
  if(material->getNumStates() == 0) {
    int numPoints = (avgnum == -1) ? getNumGaussPoints() : numNodes();
    for(int i = 0; i < numPoints; ++i) {
      result[i] = 0;
    }
    return;
  }

  // TODO: replace averaging below with extrapolation using shape functions
  double average = 0;
  for(int i = 0; i < getNumGaussPoints(); i++) {
     double d = material->getDamage(statenp + i*material->getNumStates());
     if(avgnum == -1) result[i] = d;
     else average += d;
  }

  if(avgnum != -1) {
    average /= getNumGaussPoints();

    for(int i = 0; i <  numNodes(); ++i) {
      result[i] = average;
    }
  }
}

double
GaussIntgElement::getStrainEnergy(Node *nodes, double *dispnp, double *state, double *temps)
{
  int ndofs = numDofs();
  ShapeFunction *shapeF = getShapeFunction();

  // Obtain the strain function. It can be linear or non-linear
  StrainEvaluator *strainEvaluator = getStrainEvaluator();

  // Obtain the storage for gradU ( 3x3 )
  Tensor &gradUnp = *shapeF->getGradUInstance();

  Tensor &enp = *strainEvaluator->getStrainInstance();

  // Obtain the storage for cached quantities (if any)
  Tensor *cache = strainEvaluator->getCacheInstance();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  int i,j;
  int ngp = getNumGaussPoints();
  int nstatepgp = material->getNumStates();

  double W = 0;

  for(i = 0; i < ngp; i++) {

    StackVector dispVecnp(dispnp, ndofs);

    double point[3], weight, jac, tempnp;

    getGaussPointAndWeight(i, point, weight);

    shapeF->getGradU(&gradUnp, &jac, nodes, point, dispVecnp);

    strainEvaluator->getE(enp, gradUnp, cache, state + i*material->getNumStates());

    // compute temperature at integration point
    tempnp = (temps) ? shapeF->interpolateScalar(temps, point) : 0;

    W += (weight * fabs(jac))*material->getStrainEnergyDensity(enp, state + i*material->getNumStates(), tempnp);
  }

  if(cache) delete cache;

  return W;
}

double
GaussIntgElement::getDissipatedEnergy(Node *nodes, double *state)
{
  // fprintf(stdout, "In GaussIntgElement::getDissipatedEnergy()\n");
  int ndofs = numDofs();
  ShapeFunction *shapeF = getShapeFunction();

  // Obtain the material model
  NLMaterial *material = getMaterial();

  int i,j;
  int ngp = getNumGaussPoints();
  int nstatepgp = material->getNumStates();
  
  double D = 0;

  for(i = 0; i < ngp; i++) {

    double point[3], weight, jac;
 
    getGaussPointAndWeight(i, point, weight);

    shapeF->getJacobianDeterminant(&jac, nodes, point);

    D += (weight * fabs(jac))*(material->getDissipatedEnergy(state + i*material->getNumStates()));

  }

  return D;
}

bool
GaussIntgElement::checkFailure(double *statenp)
{
  NLMaterial *material = getMaterial();
  for(int i = 0; i < getNumGaussPoints(); i++) {
    if(material->getDamage(statenp + i*material->getNumStates()) < 1) return false;
  }
  return true;
}
