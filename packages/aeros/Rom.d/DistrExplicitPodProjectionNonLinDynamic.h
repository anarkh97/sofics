#ifndef ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMIC_H
#define ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMIC_H

#include "DistrExplicitPodProjectionNonLinDynamicBase.h"

#include <Problems.d/DynamProbTraits.h>

#include <memory>

namespace Rom {

class DistrExplicitPodProjectionNonLinDynamic : public DistrExplicitPodProjectionNonLinDynamicBase {
public:
  explicit DistrExplicitPodProjectionNonLinDynamic(Domain *);

  // Overriding via hiding
  void preProcess(); // Additional pre-processing

  // Intercept passing information (2 overloads)
  void computeExtForce2(SysState<DistrVector> &distState,
                        DistrVector &f, DistrVector &cnst_f,
                        int tIndex, double t,
                        DistrVector *aero_f,
                        double gamma, double alphaf);
  void computeExtForce2(SysState<DistrVector> &distState,
                        DistrVector &f, DistrVector &cnst_f,
                        int tIndex, double t,
                        DistrVector *aero_f);

  // Added functionality
  void currentTimeIs(double t);
  void stateSnapshotAdd(const DistrVector &s);
  void accelerationSnapshotAdd(const DistrVector &s);
  void velocitySnapshotAdd(const DistrVector &s);
  void forceSnapshotAdd(const DistrVector &s);

  ~DistrExplicitPodProjectionNonLinDynamic();

  bool collectState;
  bool collectAccel;
  bool collectVeloc;
  bool collectForce;

protected:
  class SnapshotHandler;

private:
  friend class SnapshotHandler;
  std::unique_ptr<SnapshotHandler> snapshotHandler_;
};

inline
void
DistrExplicitPodProjectionNonLinDynamic::computeExtForce2(SysState<DistrVector> &distState,
                                                     DistrVector &f, DistrVector &cnst_f,
                                                     int tIndex, double t,
                                                     DistrVector *aero_f,
                                                     double gamma, double alphaf) {
  currentTimeIs(t);
  DistrExplicitPodProjectionNonLinDynamicBase::computeExtForce2(distState, f, cnst_f, tIndex, t, aero_f, gamma, alphaf);
} 

inline
void
DistrExplicitPodProjectionNonLinDynamic::computeExtForce2(SysState<DistrVector> &distState,
                                                     DistrVector &f, DistrVector &cnst_f,
                                                     int tIndex, double t,
                                                     DistrVector *aero_f) {
  currentTimeIs(t);
  DistrExplicitPodProjectionNonLinDynamicBase::computeExtForce2(distState, f, cnst_f, tIndex, t, aero_f);
}

} // end namespace Rom

inline
void
handleDisplacement(Rom::DistrExplicitPodProjectionNonLinDynamic &probDesc, DistrVector &d) {
 if(probDesc.collectState)
   probDesc.stateSnapshotAdd(d);
}

inline
void
handleAcceleration(Rom::DistrExplicitPodProjectionNonLinDynamic &probDesc, DistrVector &d) {
  if(probDesc.collectAccel)
    probDesc.accelerationSnapshotAdd(d);
}

inline
void
handleVelocity(Rom::DistrExplicitPodProjectionNonLinDynamic &probDesc, DistrVector &d) {
  if(probDesc.collectVeloc)
    probDesc.velocitySnapshotAdd(d);
}

inline
void
handleForce(Rom::DistrExplicitPodProjectionNonLinDynamic &probDesc, DistrVector &f) {
  if(probDesc.collectForce)
    probDesc.forceSnapshotAdd(f);
}

#endif /* ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMIC_H */
