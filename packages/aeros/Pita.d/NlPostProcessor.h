#ifndef PITA_NLPOSTPROCESSOR_H
#define PITA_NLPOSTPROCESSOR_H

#include "PostProcessor.h"

#include "NlDynamTimeIntegrator.h"
#include "PitaNonLinDynam.h"

namespace Pita {

class NlPostProcessor : public GenPostProcessor<NlDynamTimeIntegrator> {
public:
  EXPORT_PTRINTERFACE_TYPES(NlPostProcessor);

  virtual void outputNew(FileSetId fileSetId, const NlDynamTimeIntegrator * integrator);

  PitaNonLinDynamic *  basePostProcessor() const { return basePostProcessor_; }

  static Ptr New(GeoSource * geoSource,
                 int localFileCount,
                 const int * localFileId,
                 PitaNonLinDynamic * basePostProcessor) {
    return new NlPostProcessor(geoSource, localFileCount, localFileId, basePostProcessor);
  }

protected:
  NlPostProcessor(GeoSource * geoSource,
                  int localFileCount,
                  const int * localFileId,
                  PitaNonLinDynamic * basePostProcessor);

private:
  PitaNonLinDynamic * basePostProcessor_;
  Vector dummy_;
};

} /* end namespace Pita */

#endif /* PITA_NLPOSTPROCESSOR_H */
