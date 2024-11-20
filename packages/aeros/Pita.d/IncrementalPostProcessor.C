#include "IncrementalPostProcessor.h"
#include "LinearGenAlphaIntegrator.h"
#include <Problems.d/DynamDescr.h>
#include <Driver.d/SysState.h>

namespace Pita {

IncrementalPostProcessor::IncrementalPostProcessor(GeoSource * gs, int lfc, const int * lfi, SDDynamPostProcessor * bpp) :
  GenPostProcessor<LinearGenAlphaIntegrator>(gs, lfc, lfi),
  basePostProcessor_(bpp),
  constantTermMap_()
{}

void
IncrementalPostProcessor::outputNew(FileSetId fileSetId, const LinearGenAlphaIntegrator * oi) {
  this->fileStatusIs(fileSetId, OPEN);

  FileSetMap::iterator it = constantTermMap_.lower_bound(fileSetId);
  if (it == constantTermMap_.end() || it->first != fileSetId) {
    it = constantTermMap_.insert(it, std::make_pair(fileSetId, DirectionMap()));
  }

  bool forwardIntegration = (oi->timeStepSize().value() > 0.0);
  DirectionMap::iterator jt = it->second.lower_bound(forwardIntegration);
  if (jt == it->second.end() || jt->first != forwardIntegration) {
    jt = it->second.insert(jt, std::make_pair(forwardIntegration, ConstantTermMap()));
  }

  ConstantTermMap::iterator kt = jt->second.lower_bound(oi->timeStepCount());
  if (kt != jt->second.end() && kt->first == oi->timeStepCount()) {
    // Add reference state 
    kt->second += oi->currentState();
  } else {
    // New reference state
    kt = jt->second.insert(kt, std::make_pair(oi->timeStepCount(), oi->currentState()));
  }

  const DynamState & affineState = kt->second;

  SysState<Vector> sysState(const_cast<Vector &>(affineState.displacement()),
                            const_cast<Vector &>(affineState.velocity()),
                            const_cast<Vector &>(oi->currentAcceleration()),
                            const_cast<Vector &>(oi->previousVelocity()));
  
  GenDynamMat<double> * dynamMat = const_cast<GenDynamMat<double> *>(oi->dynamOps()->dynamMat());

  this->basePostProcessor()->pitaDynamOutput(oi->timeStepCount().value(),
                                             *dynamMat,
                                             const_cast<Vector &>(oi->externalForce()),
                                             const_cast<Vector *>(&oi->aeroForce()),
                                             sysState,
                                             fileSetId.value(),
                                             oi->currentTime().value());
}

} // end namespace Pita
