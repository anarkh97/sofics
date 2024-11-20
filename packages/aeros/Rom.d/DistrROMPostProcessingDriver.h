#ifndef ROM_DISTRROMPOSTPROCESSINGDRIVER_H
#define ROM_DISTRROMPOSTPROCESSINGDRIVER_H

#include "DriverInterface.h"
#include "Paral.d/MDDynam.h"
#include "Driver.d/DynamProbType.h"

#include "DistrVecBasis.h"

namespace Rom {

class DistrROMPostProcessingDriver : public MultiDomainDynam, public DriverInterface {

public:
  DistrROMPostProcessingDriver(Domain *);
  ~DistrROMPostProcessingDriver();

  virtual void solve();

private:
  void preProcess();
  void bufferReducedFiles();
  void setPODsize();
  void subUpdateStates(int i, double time);

  MultiDomDynPostProcessor *mddPostPro;

  DistrVecBasis  normalizedBasis_;

  SysState<GenDistrVector<double> > * curState;

  GenDistrVector<double> * fullDispBuffer;
  GenDistrVector<double> * fullVelBuffer;
  GenDistrVector<double> * fullAccBuffer;
  GenDistrVector<double> * fullVel2Buffer;
  GenDistrVector<double> * fullDummyBuffer;
  GenDistrVector<double> * fullConstForceBuffer;
  GenDistrVector<double> * fullExtForceBuffer;
  GenDistrVector<double> * fullInertForceBuffer;
  GenDistrVector<double> * fullResBuffer;

  std::vector< std::vector<double> >  TimeStamps;
  std::vector<double>                 reducedAccBuffer;
  std::vector<double>                 reducedDispBuffer;
  std::vector<double>                 reducedVelBuffer;
  std::vector<std::pair<int,int> >    DataType;

  MDDynamMat * dummyDynOps;

  int numConversionFiles;
  int projectionSubspaceSize;
  int VECsize;
};

} //end namespace ROM

Rom::DriverInterface *distrROMPostProcessingDriverNew(MultiDomainDynam *);

#endif /* ROM_DISTRROMPOSTPROCESSINGDRIVER_H*/
