#include "NlDynamOps.h"
#include "PitaNonLinDynam.h"
#include <Corotational.d/GeomState.h>

namespace Pita {

NlDynamOps::NlDynamOps(PitaNonLinDynamic * probDesc) :
  probDesc_(probDesc),
  refState_(probDesc->createGeomState()),
  geomState_(probDesc->createGeomState()),
  internalForceDummy_(probDesc->elemVecInfo(), 0.0),
  residualDummy_(probDesc->solVecInfo(), 0.0) 
{}

NlDynamOps::~NlDynamOps() {
  delete geomState_;
  delete refState_;
}

void
NlDynamOps::displacementIs(const GenVector<double> & disp, Seconds time) {
  *geomState_ = *refState_;
  geomState_->update(disp);
  probDesc_->getStiffAndForce(*geomState_, residualDummy_, internalForceDummy_, time.value());
  probDesc_->reBuildKonly();
}

const GenSparseMatrix<double> *
NlDynamOps::massMatrix() const {
  return probDesc_->getMassMatrix();
}

const GenSparseMatrix<double> *
NlDynamOps::stiffnessMatrix() const {
  return probDesc_->getStiffMatrix();
}

} // end namespace Pita
