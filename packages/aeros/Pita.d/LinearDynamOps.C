#include "LinearDynamOps.h"
#include <Problems.d/DynamDescr.h>

namespace Pita {

// LinearDynamOps implementation

LinearDynamOps::LinearDynamOps(LinearDynamOps::DynOpsType * allOps) :
  allOps_(allOps)
{}

LinearDynamOps::~LinearDynamOps() {
  delete allOps_;
}

const LinearDynamOps::MatrixType *
LinearDynamOps::massMatrix() const {
  return allOps_->M;
}

const LinearDynamOps::MatrixType *
LinearDynamOps::stiffnessMatrix() const {
  return allOps_->K;
}

const LinearDynamOps::MatrixType *
LinearDynamOps::dampingMatrix() const {
  return allOps_->C;
}

const LinearDynamOps::SolverType *
LinearDynamOps::dynamicMassSolver() const {
  return allOps_->dynMat;
}

LinearDynamOps::MatrixType *
LinearDynamOps::massMatrix() {
  return allOps_->M;
}

LinearDynamOps::MatrixType *
LinearDynamOps::stiffnessMatrix() {
  return allOps_->K;
}

LinearDynamOps::MatrixType *
LinearDynamOps::dampingMatrix() {
  return allOps_->C;
}

LinearDynamOps::SolverType *
LinearDynamOps::dynamicMassSolver() {
  return allOps_->dynMat;
}

// GeneralizedAlphaParameter implementation

GeneralizedAlphaParameter::GeneralizedAlphaParameter(Seconds stepSize, double rhoInf) :
  timeStepSize_(stepSize)
{
  rhoInfinityIs(rhoInf);
}

void
GeneralizedAlphaParameter::rhoInfinityIs(double rhoInf) {
  alpham_ = (2.0 * rhoInf - 1.0) / (rhoInf + 1.0);
  alphaf_ = rhoInf / (rhoInf + 1.0);
  beta_   = 1.0 - (alpham_ - alphaf_);
  beta_  *= (0.25 * beta_);
  gamma_  = 0.5 - (alpham_ - alphaf_);
  rhoInfinity_ = rhoInf;
}

void 
GeneralizedAlphaParameter::timeStepSizeIs(Seconds dt) {
  timeStepSize_ = dt;
}

bool
GeneralizedAlphaParameter::operator==(const GeneralizedAlphaParameter & other) const {
  return (timeStepSize() == other.timeStepSize()) && (rhoInfinity() == other.rhoInfinity());
}

bool
GeneralizedAlphaParameter::operator<(const GeneralizedAlphaParameter & other) const {
  if (timeStepSize() < other.timeStepSize()) {
    return true;
  } else if (timeStepSize() == other.timeStepSize()) {
    return (rhoInfinity() < other.rhoInfinity());
  }
  return false;
}


// LinearDynamOps::Manager implementation

LinearDynamOps::Manager::Manager(LinearDynamOps::Manager::ProblemDescriptor * pbDesc) :
  probDesc_(pbDesc),
  noPhysicalDamping_(pbDesc ? pbDesc->alphaDamp() == 0.0 && pbDesc->betaDamp() == 0.0 : false)
{}

LinearDynamOps::Ptr
LinearDynamOps::Manager::dynOpsNew(const GeneralizedAlphaParameter & param) {
  // Find insertion point
  DynOpsMap::iterator it = dynOpsMap_.lower_bound(param);
  if (it != dynOpsMap_.end() && it->first == param)
    throw Fwk::NameInUseException();

  LinearDynamOps::Ptr newDynOps;

  // Possible optimization when no physical damping
  if (noPhysicalDamping_) {
    GeneralizedAlphaParameter alternateParam(-param.timeStepSize(), param.rhoInfinity());
    newDynOps = this->dynOps(alternateParam); 
  }

  if (!newDynOps) {
    // Build new DynOps
    double dt = param.timeStepSize().value();
    double coefM = (1.0 - param.alpham()) / (1.0 - param.alphaf());
    double coefC = dt * param.gamma();
    double coefK = (dt * dt) * param.beta();
    LinearDynamOps::DynOpsType * newAllOps = probDesc_->buildOps(coefM, coefC, coefK); // Assemble + Factor
    newDynOps = new LinearDynamOps(newAllOps);
  }

  // Add new DynOps
  dynOpsMap_.insert(it, std::make_pair(param, newDynOps));
  return newDynOps; 
}

LinearDynamOps::Ptr
LinearDynamOps::Manager::dynOps(const GeneralizedAlphaParameter & param) const {
  DynOpsMap::const_iterator it = dynOpsMap_.find(param);
  return (it != dynOpsMap_.end()) ? it->second.ptr() : NULL;
}

void
LinearDynamOps::Manager::dynOpsDel(const GeneralizedAlphaParameter & param) {
  dynOpsMap_.erase(param);
}

size_t
LinearDynamOps::Manager::dynOpsCount() const {
  return dynOpsMap_.size();
}

} // end namespace Pita
