#include "NlPostProcessor.h"

namespace Pita {

NlPostProcessor::NlPostProcessor(GeoSource * geoSource,
                                 int localFileCount,
                                 const int * localFileId,
                                 PitaNonLinDynamic * basePostProcessor) :
  GenPostProcessor<NlDynamTimeIntegrator>(geoSource, localFileCount, localFileId),
  basePostProcessor_(basePostProcessor),
  dummy_(basePostProcessor ? basePostProcessor->solVecInfo() : 0)
{}

void
NlPostProcessor::outputNew(NlPostProcessor::FileSetId fileSetId, const NlDynamTimeIntegrator * integrator) {
  fileStatusIs(fileSetId, OPEN);
 
  basePostProcessor()->pitaDynamOutput(
      fileSetId.value(),
      const_cast<GeomState*>(integrator->geomState()),
      const_cast<Vector&>(integrator->currentState().velocity()),
      dummy_,
      integrator->currentTime().value(),
      integrator->timeStepCount().value(),
      const_cast<Vector&>(integrator->externalForce()),
      dummy_,
      const_cast<Vector&>(integrator->acceleration()));
}

} /* end namespace Pita */
