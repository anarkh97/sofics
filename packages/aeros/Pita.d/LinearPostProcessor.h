#ifndef PITA_LINEARPOSTPROCESSOR_H
#define PITA_LINEARPOSTPROCESSOR_H

#include "Fwk.h"
#include "PostProcessor.h"

class SDDynamPostProcessor;

namespace Pita {

class LinearGenAlphaIntegrator;

class LinearPostProcessor : public GenPostProcessor<LinearGenAlphaIntegrator> {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearPostProcessor);

  virtual void outputNew(FileSetId fileSetId, const LinearGenAlphaIntegrator * oi);

  SDDynamPostProcessor * basePostProcessor() const { return basePostProcessor_; }

  static Ptr New(GeoSource * geoSource,
                 int localFileCount,
                 const int * localFileId,
                 SDDynamPostProcessor * basePostProcessor) {
    return new LinearPostProcessor(geoSource, localFileCount, localFileId, basePostProcessor);
  }

protected:
  LinearPostProcessor(GeoSource * gs, int lfc, const int * lfi, SDDynamPostProcessor * bpp);

private:
  SDDynamPostProcessor * basePostProcessor_;
};

} // namespace Pita

#endif /* PITA_LINEARPOSTPROCESSOR_H */
