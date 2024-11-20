#ifndef ROM_ROMPOSTPROCESSINGDRIVER_H
#define ROM_ROMPOSTPROCESSINGDRIVER_H

#include "DriverInterface.h"
#include "Problems.d/NonLinDynam.h"
#include "Driver.d/DynamProbType.h"

#include "VecBasis.h"

namespace Rom {

class ROMPostProcessingDriver : public NonLinDynamic, public DriverInterface {

public:
  ROMPostProcessingDriver(Domain *);
  ~ROMPostProcessingDriver();

  virtual void solve();

private:
  void preProcess();
  void bufferReducedFiles();
  void setPODsize();

  VecBasis normalizedBasis_;
  VecBasis adjointBasis_;

  SysState<GenVector<double> > * curState;

  GenVector<double> * fullDispBuffer;
  GenVector<double> * fullVelBuffer;
  GenVector<double> * fullAccBuffer;
  GenVector<double> * fullVel2Buffer;
  GenVector<double> * fullDummyBuffer;

  std::vector< std::vector<double> >  TimeStamps;
  std::vector<double>                 reducedAccBuffer;
  std::vector<double>                 reducedDispBuffer;
  std::vector<double>                 reducedVelBuffer;
  std::vector<std::pair<int,int> >    DataType;

  int numConversionFiles;
  int projectionSubspaceSize;
  int VECsize;
};

} //end namespace ROM

Rom::DriverInterface *ROMPostProcessingDriverNew(NonLinDynamic *);

#endif /* ROM_ROMPOSTPROCESSINGDRIVER_H*/
