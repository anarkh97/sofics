#ifndef _MD_MODAL_BASE_H_
#define _MD_MODAL_BASE_H_

#include <Math.d/Vector.h>
#include <Feti.d/DistrVector.h>
#include <Feti.d/DistrVectorSet.h>

class Domain;
template <class Scalar> class GenDecDomain;
typedef GenDecDomain<double> DecDomain;
template <class V> class SysState;
struct ModalOps;

//------------------------------------------------------------------------------
//******************************************************************************
//------------------------------------------------------------------------------

class MDModalBase
{
/* base class for ModalDescr
   NOTE to self: if modalizing dsp/vel/acc is desired,
     will need to store scale as data memeber
*/
protected:

  Domain *domain;
  DecDomain *decDomain;
  GenDistrVectorSet<double> *modesFl;   // eigenmodes of non-zero freq
  GenDistrVectorSet<double> *modesRB;   // rigid body modes
  double *freqs;     // circular frequencies (rad/s) of modesFl
  int numModes, numRBM, numFlex;

  DistrVector *fullTmpF, *fullTmpGrav, *fullAeroF;
  DistrVector *fullDsp, *fullVel, *fullAcc, *fullPrevVel;
  DistrVector *prevFrc, *prevFrcBackup;

public:

  virtual ~MDModalBase() {/* TODO */};  
  MDModalBase(){}
  MDModalBase(Domain *d){ domain = d; }

  virtual void preProcess() = 0;
  void preProcessBase();
  void populateRBModes();
  void populateFlexModes(double scale = 1.0, bool readAll = 0);

  void initStateBase(Vector& dsp, Vector& vel, Vector& acc,
    Vector& vel_p, int idxOffset = 0);

  void outputModal(SysState<Vector> &state, Vector& extF, int tIndex, ModalOps &ops);
};

#endif
