#include <Corotational.d/SuperCorotator.h>
#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dbg_alloca.h>

SuperCorotator::SuperCorotator(SuperElement *_superElem)
{
  superElem = _superElem;
  nSubElems = superElem->getNumSubElems();
  subElemCorotators = new Corotator * [nSubElems];
  origK = 0;
  sub_vld = 0;
  sub_dvld = 0;
  sub_vlr = 0;
}

SuperCorotator::~SuperCorotator() 
{ 
  for(int i = 0; i < nSubElems; ++i)
    if(!dynamic_cast<Element*>(subElemCorotators[i])) delete subElemCorotators[i];
  delete [] subElemCorotators; 
  if(origK) delete origK;
  if(sub_dvld) {
    int i;
    for(i=0; i<nSubElems; ++i)
      if(sub_dvld[i]) delete [] sub_dvld[i];
    delete [] sub_dvld;
  }
  if(sub_vld) {
    int i;
    for(i=0; i<nSubElems; ++i)
      if(sub_vld[i]) delete [] sub_vld[i];
    delete [] sub_vld;
  }
  if(sub_vlr) {
    int i;
    for(i=0; i<nSubElems; ++i)
      if(sub_vlr[i]) delete [] sub_vlr[i];
    delete [] sub_vlr;
  }
}

void
SuperCorotator::getStiffAndForce(GeomState &geomState, CoordSet &cs,
                                 FullSquareMatrix &elK, double *f, double dt, double t)
{
  int i, j;
  elK.zero();
  for(i=0; i<elK.dim(); ++i) f[i] = 0.0;

  for(i=0; i<nSubElems; ++i) {
    int ndofs = superElem->getSubElemNumDofs(i);
    FullSquareMatrix subK(ndofs);
    subK.zero();  // this is necessary because getStiffAndForce may not be implemented for all subelems
    double *subf = new double[ndofs];
    for(j=0; j<ndofs; ++j) subf[j] = 0.0;
    subElemCorotators[i]->getStiffAndForce(geomState, cs, subK, subf, dt, t);
    int *subElemDofs = superElem->getSubElemDofs(i);
    elK.add(subK, subElemDofs);
    for(j=0; j<ndofs; ++j) f[subElemDofs[j]] += subf[j];
    delete [] subf;
  }
}

void 
SuperCorotator::getStiffAndForce(GeomState *refState, GeomState &geomState, CoordSet &cs,
                                 FullSquareMatrix &elK, double *f, double dt, double t) 
{
  int i, j;
  elK.zero();
  for(i=0; i<elK.dim(); ++i) f[i] = 0.0;

  for(i=0; i<nSubElems; ++i) {
    int ndofs = superElem->getSubElemNumDofs(i);
    FullSquareMatrix subK(ndofs);
    subK.zero();  // this is necessary because getStiffAndForce may not be implemented for all subelems
    double *subf = new double[ndofs];
    for(j=0; j<ndofs; ++j) subf[j] = 0.0;
    subElemCorotators[i]->getStiffAndForce(refState, geomState, cs, subK, subf, dt, t);
    int *subElemDofs = superElem->getSubElemDofs(i);
    elK.add(subK, subElemDofs);
    for(j=0; j<ndofs; ++j) f[subElemDofs[j]] += subf[j];
    delete [] subf;
  }
}

void
SuperCorotator::getDExternalForceDu(GeomState &geomState, CoordSet &cs,
                                    FullSquareMatrix &elK, double *f)
{
  int i, j;

  for(i=0; i<nSubElems; ++i) {
    int ndofs = superElem->getSubElemNumDofs(i);
    FullSquareMatrix subK(ndofs);
    subK.zero();
    int *subElemDofs = superElem->getSubElemDofs(i);
    double *subf = superElem->getPreviouslyComputedSubExternalForce(i); 
    subElemCorotators[i]->getDExternalForceDu(geomState, cs, subK, subf);
    elK.add(subK, subElemDofs);
  }
}

void
SuperCorotator::getInternalForce(GeomState &geomState, CoordSet &cs,
                                 FullSquareMatrix &elK, double *f, double dt, double t)
{
  int i, j;                            
  for(i=0; i<elK.dim(); ++i) f[i] = 0.0;
  double *subf = (double *) dbg_alloca(sizeof(double)*elK.dim());

  for(i=0; i<nSubElems; ++i) {
    int ndofs = superElem->getSubElemNumDofs(i);
    for(j=0; j<ndofs; ++j) subf[j] = 0.0;
    subElemCorotators[i]->getInternalForce(geomState, cs, elK, subf, dt, t);
    int *subElemDofs = superElem->getSubElemDofs(i);
    for(j=0; j<ndofs; ++j) f[subElemDofs[j]] += subf[j];
  }
}

void
SuperCorotator::getInternalForce(GeomState *refState, GeomState &geomState, CoordSet &cs,
                                 FullSquareMatrix &elK, double *f, double dt, double t)
{
  int i, j;
  for(i=0; i<elK.dim(); ++i) f[i] = 0.0;
  double *subf = (double *) dbg_alloca(sizeof(double)*elK.dim());

  for(i=0; i<nSubElems; ++i) {
    int ndofs = superElem->getSubElemNumDofs(i);
    for(j=0; j<ndofs; ++j) subf[j] = 0.0;
    subElemCorotators[i]->getInternalForce(refState, geomState, cs, elK, subf, dt, t);
    int *subElemDofs = superElem->getSubElemDofs(i);
    for(j=0; j<ndofs; ++j) f[subElemDofs[j]] += subf[j];
  }
}


void
SuperCorotator::getExternalForce(GeomState &geomState, CoordSet &cs,
                                 double *f)
{
  int i, j;                            
  double *fg = new double[superElem->numDofs()];
  for(i=0; i<superElem->numDofs(); ++i) fg[i] = 0.0;

  for(i=0; i<nSubElems; ++i) {
    double *subf = superElem->getPreviouslyComputedSubExternalForce(i);
    if(subf) {
      int ndofs = superElem->getSubElemNumDofs(i);
      int *subElemDofs = superElem->getSubElemDofs(i);
      subElemCorotators[i]->getExternalForce(geomState, cs, subf);
      for(j=0; j<ndofs; ++j) fg[subElemDofs[j]] += subf[j];
    }
  }

  for(i=0; i<superElem->numDofs(); ++i) f[i] = fg[i];
  delete [] fg;
}

void 
SuperCorotator::formGeometricStiffness(GeomState &geomState, CoordSet &cs, 
                                       FullSquareMatrix &elK, double *f)
{
  int i, j;
  elK.zero();
  for(i=0; i<superElem->numDofs(); ++i) f[i] = 0.0;

  for(i=0; i<nSubElems; ++i) {
    int ndofs = superElem->getSubElemNumDofs(i);
    FullSquareMatrix subK(ndofs);
    subK.zero();  // this is necessary because formGeometricStiffness may not be implemented for all subelems
    double *subf = new double[ndofs];
    for(j=0; j<ndofs; ++j) subf[j] = 0.0;
    subElemCorotators[i]->formGeometricStiffness(geomState, cs, subK, subf);
    int *subElemDofs = superElem->getSubElemDofs(i);
    elK.add(subK, subElemDofs);
    for(j=0; j<ndofs; ++j) f[subElemDofs[j]] += subf[j];
    delete [] subf;
  }
}

double *
SuperCorotator::getOriginalStiffness() 
{
  if(!origK) origK = new FullM(superElem->numDofs());
  origK->zero();
  int i;
  for(i=0; i<nSubElems; ++i) {
    double *subKdata = subElemCorotators[i]->getOriginalStiffness();
    if(subKdata) {
      int ndofs = superElem->getSubElemNumDofs(i);
      FullM subK(subKdata, ndofs, ndofs);
      origK->add(subK, superElem->getSubElemDofs(i));
    }
  }
  return origK->data();
}

void 
SuperCorotator::extractDeformations(GeomState &geomState, CoordSet &cs, double *vld, int &nlflag) 
{
  int i, j;
  if(!sub_vld) { 
    sub_vld = new double * [nSubElems];
    for(i=0; i<nSubElems; ++i) sub_vld[i] = 0;
  }
  for(i=0; i<superElem->numDofs(); ++i) vld[i] = 0.0; // this shouldn't be used, should always use sub_vld 
                                                      // for non-linear because the same dof can have 2
                                                      //different vld displacements in different sub elements
  nlflag = 999;

  for(i=0; i<nSubElems; ++i) {
    int ndofs = superElem->getSubElemNumDofs(i);
    if(!sub_vld[i]) sub_vld[i] = new double[ndofs];
    for(j=0; j<ndofs; ++j) sub_vld[i][j] = 0.0;
    int subnlflag;
    subElemCorotators[i]->extractDeformations(geomState, cs, sub_vld[i], subnlflag);
    if(subnlflag < nlflag) nlflag = subnlflag;  // always take the lowest
  }
}

void 
SuperCorotator::getNLVonMises(Vector &stress, Vector &weight, GeomState &geomState, GeomState *refState,
                              CoordSet &cs, int strInd, int surface, double ylayer, double zlayer,
                              int avgnum, int measure)
{
  int i;
  stress.zero();
  weight.zero();

  for(i=0; i<nSubElems; ++i) {
    int nnodes = superElem->getSubElemNumNodes(i);
    Vector subStress(nnodes);
    subStress.zero();
    Vector subWeight(nnodes);
    subWeight.zero();
    int subAvgnum = (avgnum == 0) ? 1 : avgnum;
    subElemCorotators[i]->getNLVonMises(subStress, subWeight, geomState, refState, cs, strInd, surface,
                                        ylayer, zlayer, avgnum, measure);
    int *subElemNodes = superElem->getSubElemNodes(i);
    stress.add(subStress, subElemNodes);
    weight.add(subWeight, subElemNodes);
  }

  // Average stress/strain value at each node by the number of sub-elements attached to the node
  for(i = 0; i < superElem->numNodes(); ++i) {
    if(weight[i] == 0) {
      stress[i] = 0;
    }
    else {
      stress[i] /= weight[i];
      weight[i] = 1;
    }
  }
}

void 
SuperCorotator::getNLAllStress(FullM &stress, Vector &weight, GeomState &geomState, GeomState *refState,
                               CoordSet &cs, int strInd, int surface, int measure) 
{
  int i;
  stress.zero();
  weight.zero();

  for(i=0; i<nSubElems; ++i) {
    int nnodes = superElem->getSubElemNumNodes(i);
    FullM subStress(nnodes, 9);
    subStress.zero();
    Vector subWeight(nnodes);
    subWeight.zero();
    subElemCorotators[i]->getNLAllStress(subStress, subWeight, geomState, refState, cs, strInd, surface,
                                         measure);
    int *subElemNodes = superElem->getSubElemNodes(i);
    stress.addrows(subStress, subElemNodes);
    weight.add(subWeight, subElemNodes);
  }

  // Average stress/strain values at each node by the number of sub-elements attached to the node
  for(i = 0; i < superElem->numNodes(); ++i) {
    if(weight[i] == 0) {
      for(int j = 0; j < stress.numCol(); ++j) stress[i][j] = 0;
    }
    else {
      for(int j = 0; j < stress.numCol(); ++j) stress[i][j] /= weight[i];
      weight[i] = 1;
    }
  }
}

double 
SuperCorotator::getElementEnergy(GeomState &geomState, CoordSet &cs) 
{
  int i;
  double ret = 0.0;
  for(i=0; i<nSubElems; ++i)
    ret += subElemCorotators[i]->getElementEnergy(geomState, cs);
  return ret;
}

double
SuperCorotator::getDissipatedEnergy(GeomState &geomState, CoordSet &cs)
{
  int i;
  double ret = 0.0;
  for(i=0; i<nSubElems; ++i)
    ret += subElemCorotators[i]->getDissipatedEnergy(geomState, cs);
  return ret;
}

void 
SuperCorotator::extractRigidBodyMotion(GeomState &geomState, CoordSet &cs, double *vlr) 
{
  int i,j;
  if(!sub_vlr) {
    sub_vlr = new double * [nSubElems];
    for(i=0; i<nSubElems; ++i) sub_vlr[i] = 0;
  }
  for(i=0; i<superElem->numDofs(); ++i) vlr[i] = 0.0; // this shouldn't be used, should always use sub_vlr 
                                                      // for non-linear because the same dof can have 
                                                      // 2 different vlr displacements in different sub elements
  for(i=0; i<nSubElems; ++i) {
    int ndofs = superElem->getSubElemNumDofs(i);
    if(!sub_vlr[i]) sub_vlr[i] = new double[ndofs];
    for(j=0; j<ndofs; ++j) sub_vlr[i][j] = 0.0;
    subElemCorotators[i]->extractRigidBodyMotion(geomState, cs, sub_vlr[i]);
  }
}

void
SuperCorotator::updateStates(GeomState *refState, GeomState &curState, CoordSet &C0, double dt)
{
  int i;
  for(i=0; i<nSubElems; ++i)
    subElemCorotators[i]->updateStates(refState, curState, C0, dt);
}

bool
SuperCorotator::checkElementDeletion(GeomState &curState)
{
  int i;
  for(i=0; i<nSubElems; ++i)
    if(!subElemCorotators[i]->checkElementDeletion(curState)) return false;
  return true;
}

void
SuperCorotator::getResidualCorrection(GeomState &gs, double *r)
{
  int i, j;
  for(i=0; i<nSubElems; ++i) {
    int ndofs = superElem->getSubElemNumDofs(i);
    double *subr = new double[ndofs];
    for(j=0; j<ndofs; ++j) subr[j] = 0.0;
    subElemCorotators[i]->getResidualCorrection(gs, subr);
    int *subElemDofs = superElem->getSubElemDofs(i);
    for(j=0; j<ndofs; ++j) r[subElemDofs[j]] += subr[j];
    delete [] subr;
  }
}

void
SuperCorotator::initMultipliers(GeomState& c1)
{
  int i;
  for(i=0; i<nSubElems; ++i)
    subElemCorotators[i]->initMultipliers(c1);
}

void
SuperCorotator::updateMultipliers(GeomState& c1)
{
  int i;
  for(i=0; i<nSubElems; ++i)
    subElemCorotators[i]->updateMultipliers(c1);
}

double
SuperCorotator::getError(GeomState& c1)
{
  double err = 0;
  int i;
  for(i=0; i<nSubElems; ++i)
    err = std::max(err,subElemCorotators[i]->getError(c1));
  return err;
}

bool
SuperCorotator::useDefaultInertialStiffAndForce()
{
  int i;
  for(i=0; i<nSubElems; ++i) {
    if(!subElemCorotators[i]->useDefaultInertialStiffAndForce()) return false;
  }
  return true;
}

void
SuperCorotator::getInertialStiffAndForce(GeomState *refState, GeomState &curState, CoordSet &c0,
                                         FullSquareMatrix &elK, double *f, double dt, double t,
                                         double beta, double gamma, double alphaf, double alpham)
{
  int i, j;
  elK.zero();
  for(i=0; i<elK.dim(); ++i) f[i] = 0.0;

  for(i=0; i<nSubElems; ++i) {
    int ndofs = superElem->getSubElemNumDofs(i);
    FullSquareMatrix subK(ndofs);
    subK.zero();
    double *subf = new double[ndofs];
    for(j=0; j<ndofs; ++j) subf[j] = 0.0;
    int *subElemDofs = superElem->getSubElemDofs(i);
    subElemCorotators[i]->getInertialStiffAndForce(refState, curState, c0, subK, subf, dt, t, 
                                                   beta, gamma, alphaf, alpham);
    elK.add(subK, subElemDofs);
    for(j=0; j<ndofs; ++j) f[subElemDofs[j]] += subf[j];
    delete [] subf;
  }
}

void
SuperCorotator::getInternalForceThicknessSensitivity(GeomState *refState, GeomState &geomState, CoordSet &cs,
                                                     Vector &dFintdThick, double dt, double t)
{
  int i, j;
  for(i=0; i<dFintdThick.size(); ++i) dFintdThick[i] = 0.0;

  for(i=0; i<nSubElems; ++i) {
    int ndofs = superElem->getSubElemNumDofs(i);
    Vector subf(ndofs);
    for(j=0; j<ndofs; ++j) subf[j] = 0.0;
    subElemCorotators[i]->getInternalForceThicknessSensitivity(refState, geomState, cs, subf, dt, t);
    int *subElemDofs = superElem->getSubElemDofs(i);
    for(j=0; j<ndofs; ++j) dFintdThick[subElemDofs[j]] += subf[j];
  }
}

void
SuperCorotator::getInternalForceNodalCoordinateSensitivity(GeomState *refState, GeomState &geomState, CoordSet &cs,
                                                           Vector *&dFintdx, double dt, double t)
{
  int i, j, k;
  for(i=0; i<superElem->numNodes()*3; ++i)
    for(j=0; j<dFintdx[i].size(); ++i) dFintdx[i][j] = 0.0;

  for(i=0; i<nSubElems; ++i) {
    int nnodes = superElem->getSubElemNumNodes(i);
    int ndofs = superElem->getSubElemNumDofs(i);
    Vector *subf = new Vector[nnodes*3];
    for(j=0; j<nnodes*3; ++j) {
      subf[i].resize(ndofs);
      for(k=0; k<ndofs; ++k) subf[j][k] = 0.0;
    }
    subElemCorotators[i]->getInternalForceNodalCoordinateSensitivity(refState, geomState, cs, subf, dt, t);
    int *subElemDofs = superElem->getSubElemDofs(i);
    int *subElemNodes = superElem->getSubElemNodes(i);
    for(j=0; j<nnodes*3; ++j)
      for(k=0; k<ndofs; ++k) dFintdx[3*subElemNodes[j/3]+j%3][subElemDofs[k]] += subf[j][k];
    delete [] subf;
  }
}

void
SuperCorotator::extractDeformationsDisplacementSensitivity(GeomState &geomState, CoordSet &cs, double *dvld)
{
  int i, j, k;
  if(!sub_dvld) { 
    sub_dvld = new double * [nSubElems];
    for(i=0; i<nSubElems; ++i) sub_dvld[i] = 0;
  }
  int N = superElem->numDofs();
  for(i=0; i<N*N; ++i) dvld[i] = 0;  // this shouldn't be used, should always use sub_dvld for non-linear
                                     // because the same dof can have 2 different dvld displacement
                                     // sensitivities in different sub elements

 for(i=0; i<nSubElems; ++i) {
    int n = superElem->getSubElemNumDofs(i);
    if(!sub_dvld[i]) sub_dvld[i] = new double[n*n];
    for(j=0; j<n*n; ++j) sub_dvld[i][j] = 0.0;
    subElemCorotators[i]->extractDeformationsDisplacementSensitivity(geomState, cs, sub_dvld[i]);
    int *subElemDofs = superElem->getSubElemDofs(i);
    for(j=0; j<n; ++j)
      for(k=0; k<n; ++k) dvld[subElemDofs[j]*N+subElemDofs[k]] += sub_dvld[i][j*n+k];
  }
}
