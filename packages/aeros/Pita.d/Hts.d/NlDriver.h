#ifndef PITA_HTS_NLDRIVER_H
#define PITA_HTS_NLDRIVER_H

#include "Fwk.h"
#include "Types.h"

#include "../NlDriverImpl.h"

#include "../PitaNonLinDynam.h"
#include <Utils.d/SolverInfo.h>
#include <Comm.d/Communicator.h>

#include "SliceMapping.h"

#include "../PostProcessingManager.h"

class GeoSource;
class Domain;
class SolverInfo;
class Communicator;

namespace Pita { namespace Hts {

class NlDriver : public NlDriverImpl {
public:
  EXPORT_PTRINTERFACE_TYPES(NlDriver);

  virtual void solve();

  GeoSource * geoSource() const { return geoSource_; }
  Domain * domain() const { return domain_; }
  SolverInfo * solverInfo() const { return solverInfo_; }
  Communicator * baseComm() const { return baseComm_; }
  
  static Ptr New(PitaNonLinDynamic * pbDesc,
                 GeoSource * geoSource,
                 Domain * domain,
                 SolverInfo * solverInfo,
                 Communicator * baseComm) {
    return new NlDriver(pbDesc, geoSource, solverInfo, baseComm);
  }

protected:
  explicit NlDriver(PitaNonLinDynamic *, GeoSource *, SolverInfo *, Communicator *);

  void preprocess();
  void summarizeParameters() const;
  void solveParallel();
  PostProcessing::Manager::Ptr buildPostProcessor(); 

private:
  /* Primary sources */
  GeoSource * geoSource_;
  Domain * domain_;
  SolverInfo * solverInfo_;
  Communicator * baseComm_;
  
  /* Space-domain */
  size_t vectorSize_;
  
  /* Time-domain */ 
  Seconds fineTimeStep_;
  TimeStepCount halfSliceRatio_;
  TimeStepCount sliceRatio_;
  Seconds coarseTimeStep_;
  Seconds initialTime_;
  Seconds finalTime_;

  /* Load balancing */ 
  SliceMapping::Ptr mapping_;
  CpuRank localCpu_;

  /* Other parameters */
  IterationRank lastIteration_;
  double jumpCvgRatio_;
  double projectorTolerance_;
  bool userProvidedSeeds_;
  int globalBasisEnrichment_;

  /* PITA-specific output */
  bool jumpMagnOutput_;
};

} /* end namespace Hts */ } /* end namespace Pita */

Pita::NlDriver::Ptr nlReversiblePitaDriverNew(Pita::PitaNonLinDynamic * problemDescriptor);

#endif /* PITA_HTS_NLDRIVER_H */
