#ifndef PITA_OLD_NLDISTRTIMEDECOMPSOLVER_H
#define PITA_OLD_NLDISTRTIMEDECOMPSOLVER_H

#include "NLTimeSlice.h"
#include "TimeSliceMapping.h"
#include "PitaNonLinDynam.h"
#include "NLDynamTimeIntegrator.h"
#include "LinearizedTimeIntegrator.h"

#include <Pita.d/SimpleBuffer.h>
#include <Math.d/Vector.h>

#include <vector>

class Communicator;
class Connectivity;

namespace Pita { namespace Old {

template <typename Scalar> class DynamState;
template <typename Scalar> class DynamStateSet;

class NLDistrTimeDecompSolver
{
public:
  typedef GenVector<double> VecType;
  typedef VecType::DataType DataType;
  typedef DynamState<DataType> State;
  typedef DynamStateSet<DataType> StateSet;
    
  explicit NLDistrTimeDecompSolver(PitaNonLinDynamic * pbDesc);
  ~NLDistrTimeDecompSolver();
  void solve();

private:
  typedef std::vector<NLTimeSlice>::iterator sliceIterator;

  int myCPU;
  PitaNonLinDynamic * probDesc;
  NLDynamTimeIntegrator seqIntegrator; 
  LinearizedTimeIntegrator linIntegrator;
  Communicator * timeCom;
  TimeSliceMapping * sliceMap;
  std::vector<NLTimeSlice> timeSliceSet;
  sliceIterator firstActive, firstInactive; 
  SimpleBuffer<DataType> baseBuffer, seedBuffer;
 
  // Initialization
  void initialize();
  void buildConnectivities();
  void buildTimeSlices();

  // Algorithm steps
  void getInitialSeeds();
  void fineGridComputation();
  void improveBases();
  void clearBases();
  void computeCorrection();

  // Subroutines
  void improveBasesWithAllSeeds();
  void improveBasesWithLocalIncrements();

  // Timeslice routines 
  void coarseIntegrator(NLTimeSlice &);
  void fineIntegrator(NLTimeSlice &);
  void initializeSeqIntegrator(NLTimeSlice &, double);
  void computeSliceUpdate(NLTimeSlice &);
  void trivialSliceUpdate(NLTimeSlice &);

  // Base Management
  void reBuildLocalK(const NLTimeSlice &);
  void performOG(NLTimeSlice &);
  void addStateSetToBase(NLTimeSlice &, StateSet &);
  void addRawDataToBase(NLTimeSlice &, DataType *, int);

  // Helper functions
  void getIntegratorState(State &);
  void setIntegratorState(const State &);

  // Get number of unconstrained dofs (problem size)
  int getProbSize() const { return probDesc->solVecInfo(); }
};

inline void NLDistrTimeDecompSolver::initializeSeqIntegrator(NLTimeSlice & timeSlice, double timeStep)
{
  setIntegratorState(timeSlice.seedState);
  seqIntegrator.currentTimeIs(timeSlice.initialTime);
  seqIntegrator.timeStepIs(timeStep);
  seqIntegrator.currentTimeStepNumberIs(0);
}

inline void NLDistrTimeDecompSolver::computeSliceUpdate(NLTimeSlice & timeSlice)
{
  timeSlice.jumpState = timeSlice.seedState;
  timeSlice.jumpState -= timeSlice.oldSeedState;
  
  // Convergence criterion to be figured out
  timeSlice.converged = false;
  timeSlice.projector.projection(timeSlice.propBase, timeSlice.jumpState, timeSlice.nextSeedState);
  timeSlice.nextSeedState += timeSlice.propState;
}

inline void NLDistrTimeDecompSolver::trivialSliceUpdate(NLTimeSlice & timeSlice)
{
  timeSlice.converged = true;
  timeSlice.jumpState = 0.0;
  timeSlice.nextSeedState = timeSlice.propState;
}

inline void NLDistrTimeDecompSolver::addRawDataToBase(NLTimeSlice & timeSlice, NLDistrTimeDecompSolver::DataType * dataPtr, int numStates)
{
  timeSlice.projector.lastStateSetIs(numStates, dataPtr);
}

inline void NLDistrTimeDecompSolver::addStateSetToBase(NLTimeSlice & timeSlice, NLDistrTimeDecompSolver::StateSet & stateSet)
{
  timeSlice.projector.lastStateSetIs(stateSet);
}

inline void NLDistrTimeDecompSolver::performOG(NLTimeSlice & timeSlice)
{
  timeSlice.projector.metricIs(probDesc->getStiffMatrix(), probDesc->getMassMatrix(), timeSlice.propBase);
}

inline void NLDistrTimeDecompSolver::getIntegratorState(NLDistrTimeDecompSolver::State & state)
{
  seqIntegrator.getCurrentVelocity(state.vel());
  seqIntegrator.getCurrentDisplacement(state.disp());
}

inline void NLDistrTimeDecompSolver::setIntegratorState(const NLDistrTimeDecompSolver::State & state)
{
  seqIntegrator.setCurrentVelocity(state.vel());
  seqIntegrator.setCurrentDisplacement(state.disp());
}

} /* end namespace Old */ } /* end namespace Pita */

#endif
