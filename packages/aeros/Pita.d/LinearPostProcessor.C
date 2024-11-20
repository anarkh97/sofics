#include "LinearPostProcessor.h"
#include "LinearGenAlphaIntegrator.h"
#include <Problems.d/DynamDescr.h>
#include <Driver.d/SysState.h>

namespace Pita {

LinearPostProcessor::LinearPostProcessor(GeoSource * gs, int lfc, const int * lfi, SDDynamPostProcessor * bpp) :
  GenPostProcessor<LinearGenAlphaIntegrator>(gs, lfc, lfi),
  basePostProcessor_(bpp)
{}

void
LinearPostProcessor::outputNew(FileSetId fileSetId, const LinearGenAlphaIntegrator * oi) {
  this->fileStatusIs(fileSetId, OPEN);

  SysState<Vector> sysState(const_cast<Vector &>(oi->currentState().displacement()),
                            const_cast<Vector &>(oi->currentState().velocity()),
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
