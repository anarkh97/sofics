#include "NLDistrTimeDecompSolver.h"
#include "DynamStateSet.h"

#include <Utils.d/Connectivity.h>
#include <Comm.d/Communicator.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/Domain.h>

extern Communicator * structCom;
extern int verboseFlag;
extern GeoSource * geoSource;
extern Domain * domain;

namespace Pita { namespace Old {

NLDistrTimeDecompSolver::NLDistrTimeDecompSolver(PitaNonLinDynamic * pbdesc) :
  probDesc(pbdesc),
  seqIntegrator(*pbdesc),
  linIntegrator(*pbdesc),
  timeCom(NULL),
  sliceMap(NULL),
  firstActive(timeSliceSet.end()),
  firstInactive(timeSliceSet.end())
{
  probDesc->pitaTimers.start();
}

NLDistrTimeDecompSolver::~NLDistrTimeDecompSolver()
{
  delete sliceMap;
}

void NLDistrTimeDecompSolver::initialize()
{
  probDesc->pitaTimers.start("Driver Initialization");
  timeCom = structCom;
  myCPU  = structCom->myID(); 
  buildConnectivities(); 
  buildTimeSlices();
  probDesc->pitaTimers.stop();
}

// Main driver function
void NLDistrTimeDecompSolver::solve()
{
  initialize();
  probDesc->pitaTimers.start("Driver Solve");  
  // 1) (Sequential) Initial seeds on coarse grid + (Parallel) Local base initialization
  getInitialSeeds();
  int maxMainIter = probDesc->getKiter();
  if (maxMainIter > 0)
  {
    probDesc->pitaTimers.start("Main Iterations");
    probDesc->pitaTimers.newIteration();
    if (maxMainIter > 1)
    {
      // 2) Main Iteration loop
      int mainIter = 1;
      while(true)
      {
        //   a) (Parallel) Fine grid computation
        fineGridComputation();
        //   b) (Sequential) Correction
        computeCorrection();
        probDesc->pitaTimers.newIteration();
        if (++mainIter >= maxMainIter)
          break;
        //   c) (Sequential/Parallel) Local base improvement
        improveBases();
      }
    }
    clearBases(); // Last iteration : No correction needed
    fineGridComputation();
    probDesc->pitaTimers.stop();
  }
  probDesc->pitaTimers.stopAll();
  probDesc->printNLPitaTimerFile(myCPU); // Output Pita timing file
}

void NLDistrTimeDecompSolver::buildTimeSlices()
{
  int myNumTS = sliceMap->numSlicesOnCPU(myCPU);
  int* slicesOnMyCPU = sliceMap->slicesOnCPU(myCPU); 

  geoSource->duplicateFilesForPita(myNumTS, slicesOnMyCPU);
  timeSliceSet.resize(myNumTS);

  for (int i = 0; i < myNumTS; ++i)
  {
    timeSliceSet[i].initialize(*probDesc, slicesOnMyCPU[i]);
  }

  firstActive = timeSliceSet.begin();
  firstInactive = timeSliceSet.begin();
  int i = 0;
  while (firstInactive != timeSliceSet.end() && i < sliceMap->numMaxActiveSlices())
  {
    ++i;
    ++firstInactive;
  }
}

void NLDistrTimeDecompSolver::buildConnectivities()
{
  sliceMap = new TimeSliceMapping(probDesc->getNumTS(), structCom->numCPUs(), probDesc->getNumTSonCPU());
}

void NLDistrTimeDecompSolver::getInitialSeeds()
{
  probDesc->pitaTimers.start("Initial Seeds");
  int totalSeeds = sliceMap->numActiveSlices();
  StateSet localSeedBase(getProbSize(), totalSeeds); 
  int seedsFromFile = 0;
  if (domain->solInfo().pitaReadInitSeed) {
    seedsFromFile = std::min(geoSource->getUserProvidedSeedCount(), totalSeeds);
  }
  localSeedBase.addState();
  int seedRank = 0;
  if (seedsFromFile == 0) {
    probDesc->getInitState(localSeedBase[0]);
    seedRank = 1;
  }
  else
  {
    probDesc->pitaTimers.start("User-provided seeds");
    do
    {
      localSeedBase.addState();
      probDesc->getInitSeed(localSeedBase[seedRank], seedRank);
      ++seedRank;
    }
    while (seedRank < seedsFromFile);
    probDesc->pitaTimers.stop();
  }
  if (seedRank < totalSeeds)
  {
    probDesc->pitaTimers.start("Coarse Time-Grid Integration"); 
    setIntegratorState(localSeedBase[seedRank - 1]);
    seqIntegrator.timeStepIs(probDesc->getCoarseDt()); 
    seqIntegrator.currentTimeIs((seedRank - 1) * probDesc->getCoarseDt());
    seqIntegrator.currentTimeStepNumberIs(seedRank - 1);
    do
    {
      probDesc->pitaTimers.start("Time Loop");
      seqIntegrator.integrate(1);
      probDesc->pitaTimers.stop(); 
      localSeedBase.addState();
      getIntegratorState(localSeedBase[seedRank]);
      ++seedRank; 
    }
    while (seedRank < totalSeeds);
    probDesc->pitaTimers.stop();
  }
  for (sliceIterator it = firstActive; it != firstInactive; ++it)
  {
    it->active = true;
    it->seedState = localSeedBase[it->sliceRank];
    if (it->sliceRank < totalSeeds - 1)
    {
      it->nextSeedState = localSeedBase[it->sliceRank + 1];
    }
    it->locDispOG = it->seedState.disp();
    it->locTimeOG = it->initialTime;  
    probDesc->pitaTimers.start("Initial Base Building");
    probDesc->pitaTimers.start("Rebuild Stiffness Matrix");
    reBuildLocalK(*it);
    probDesc->pitaTimers.swap("Orthogonalization");
    addStateSetToBase(*it, localSeedBase);
    performOG(*it);
    probDesc->pitaTimers.stop();
    probDesc->pitaTimers.stop();
  }
  probDesc->pitaTimers.stop(); 
}

void NLDistrTimeDecompSolver::improveBases()
{
  if (probDesc->getBasisImprovementMethod() == 2)
    improveBasesWithLocalIncrements();
  else
    improveBasesWithAllSeeds();
}

void NLDistrTimeDecompSolver::improveBasesWithAllSeeds()
{
  probDesc->pitaTimers.start("Global Base Improvement");
  probDesc->pitaTimers.start("Global Data Exchange");

  bool usePropagatedSeeds = (probDesc->getBasisImprovementMethod() == 1);
  
  // Determine number of seeds to be exchanged -> buffer size
  int seedPerSlice = usePropagatedSeeds ? 2 : 1;
  int numSeeds = seedPerSlice * sliceMap->numActiveSlices(); // Main seeds + Propagated seeds
  int dataPerSeed = 2 * getProbSize(); 
  
  int bufferSize = numSeeds * dataPerSeed;
  baseBuffer.sizeIs(bufferSize);
 
  // Setup parameters for global MPI communication
  int numCPUs = sliceMap->numCPUs();
  int dataPerSlice = seedPerSlice * dataPerSeed;
  
  int * recv_counts = new int[numCPUs];
  int * displacements = new int[numCPUs];
  for (int i = 0; i < numCPUs; ++i)
  {
    recv_counts[i] = dataPerSlice * sliceMap->numActiveSlicesOnCPU(i);
  }
  displacements[0] = 0;
  for (int i = 1; i < numCPUs; ++i)
  {
    displacements[i] = displacements[i-1] + recv_counts[i-1];
  }
  
  // Fill the buffer
  int sliceShift = displacements[myCPU];
  for (sliceIterator it = firstActive; it != firstInactive; ++it)
  {
    it->seedState.getRaw(baseBuffer.array() + sliceShift);
    sliceShift += dataPerSeed;
    if (usePropagatedSeeds) {
      it->propState.getRaw(baseBuffer.array() + sliceShift);
      sliceShift += dataPerSeed;
    }
  }
 
  // Perform global communication
  timeCom->allGatherv<DataType>(baseBuffer.array(), recv_counts, displacements);
  
  // Garbage collection
  delete[] displacements;
  delete[] recv_counts;
 
  probDesc->pitaTimers.stop();
 
  for (sliceIterator it = firstActive; it != firstInactive; ++it)
  {
    if (it->sliceRank == sliceMap->numCvgSlices())
    {
      it->clearAllBases(); // Not necessary to propagate correction
    }
    else
    {
      probDesc->pitaTimers.start("Rebuild Stiffness Matrix");
      reBuildLocalK(*it);
      probDesc->pitaTimers.swap("Orthogonalization");
      addRawDataToBase(*it, baseBuffer.array(), numSeeds);
      performOG(*it);
      probDesc->pitaTimers.stop();
    }
  }
  probDesc->pitaTimers.stop();
}

void NLDistrTimeDecompSolver::improveBasesWithLocalIncrements()
{
  probDesc->pitaTimers.start("Local Base Improvement");
 
  int numStatesToExchange = probDesc->getJratio() + 1;
  int halfBufferSize = 2 * getProbSize() * numStatesToExchange; 
  int bufferSize = 2 * halfBufferSize;
  baseBuffer.sizeIs(bufferSize);
  
  RecInfo recInfo;
  sliceIterator it;
  int currSliceRank;
  int firstActiveSlice = sliceMap->numCvgSlices();
  
  for (sliceIterator it = firstActive; it != firstInactive; ++it)
  { 
    currSliceRank = it->sliceRank;
    if (currSliceRank != firstActiveSlice)
    {
      probDesc->pitaTimers.start("Rebuild Stiffness matrix");
      reBuildLocalK(*it);
      probDesc->pitaTimers.stop();
      if (currSliceRank + 1 < sliceMap->firstInactiveSlice())
      {
        probDesc->pitaTimers.start("Local Data Exchange");
        it->localBase.addState(it->seedState);
        it->localBase.getRaw(baseBuffer.array() +  halfBufferSize);
        probDesc->pitaTimers.stop();
      }
      probDesc->pitaTimers.start("Orthogonalization");
      addStateSetToBase(*it, it->localBase);
      probDesc->pitaTimers.stop();
      if (sliceMap->numCPU(currSliceRank - 1) != myCPU)
      {
        probDesc->pitaTimers.start("Local Data Exchange");
        timeCom->waitForAllReq();
        recInfo = timeCom->recFrom<DataType>(currSliceRank, baseBuffer.array(), halfBufferSize);
        probDesc->pitaTimers.stop();
      }
      probDesc->pitaTimers.start("Orthogonalization");
      addRawDataToBase(*it, baseBuffer.array(), numStatesToExchange);
      performOG(*it);
      probDesc->pitaTimers.stop();
    }

    ++currSliceRank;
    if (currSliceRank < sliceMap->firstInactiveSlice())
    {
      probDesc->pitaTimers.start("Data Exchange");
      if (sliceMap->numCPU(currSliceRank) != myCPU)
      {
        timeCom->sendTo<DataType>(sliceMap->numCPU(currSliceRank), currSliceRank, baseBuffer.array() + halfBufferSize, halfBufferSize);
        timeCom->waitForAllReq();
      }
      probDesc->pitaTimers.stop();
    }
 
  }

  timeCom->waitForAllReq();  
  probDesc->pitaTimers.stop();
}

void NLDistrTimeDecompSolver::clearBases()
{
  for (sliceIterator it = firstActive; it != firstInactive; ++it)
  {
    it->clearAllBases();
  } 
}

void NLDistrTimeDecompSolver::fineGridComputation()
{
  for (sliceIterator it = firstActive; it != firstInactive; ++it)
  {
    fineIntegrator(*it);
  }
}

void NLDistrTimeDecompSolver::computeCorrection()
{
  probDesc->pitaTimers.start("Correction");

  int bufferSize = 2 * getProbSize();
  seedBuffer.sizeIs(bufferSize);
  
  DynamState<DataType> currState(getProbSize());
  RecInfo recInfo;
  sliceIterator it;
  int currSliceRank;
  int firstActiveSlice = sliceMap->numCvgSlices();
  bool convergenceFlag = false;

  for (it = firstActive; it != firstInactive; ++it)
  {
    currSliceRank = it->sliceRank;
   
    if (currSliceRank > firstActiveSlice)
    {
      it->oldSeedState = it->seedState;
      if (sliceMap->numCPU(currSliceRank - 1) != myCPU)
      {
        // Receive initial value in buffer -> currState
        probDesc->pitaTimers.start("Waiting for Data");
        timeCom->waitForAllReq();
        recInfo = timeCom->recFrom<DataType>(currSliceRank, seedBuffer.array(), bufferSize);
        probDesc->pitaTimers.stop();
        currState.setRaw(seedBuffer.array());
      }
      it->seedState = currState;
      if (it->active) // TS is already active : Compute correction
      {
        probDesc->pitaTimers.start("Update");
        computeSliceUpdate(*it);
        probDesc->pitaTimers.stop();
      }
      else // TS is not currently active : Initialize and activate
      {
        it->locDispOG = currState.disp();
        it->locTimeOG = it->initialTime;
        it->active = true;
        coarseIntegrator(*it);
      }
    }
    else
    {
      probDesc->pitaTimers.start("Update");
      trivialSliceUpdate(*it);
      probDesc->pitaTimers.stop();
      it->active = false;
      convergenceFlag = true;
    }

    currState = it->nextSeedState;

    ++currSliceRank;
    if (currSliceRank < sliceMap->firstInactiveSlice() && sliceMap->numCPU(currSliceRank) != myCPU)
    {
      // Send final value
      currState.getRaw(seedBuffer.array());
      timeCom->sendTo<DataType>(sliceMap->numCPU(currSliceRank), currSliceRank, seedBuffer.array(), bufferSize);
    }

    if (convergenceFlag)
    {
      ++firstActive;
      if (firstInactive != timeSliceSet.end()) ++firstInactive;
      convergenceFlag = false;
    }
  }

  timeCom->waitForAllReq();
  sliceMap->incrementNumCvgSlices();
  probDesc->pitaTimers.stop();
}

void NLDistrTimeDecompSolver::reBuildLocalK(const NLTimeSlice & ts)
{
  GeomState * geomState = probDesc->createGeomState();
  geomState->update(ts.locDispOG);
  VecType elementInternalForce(probDesc->elemVecInfo(), 0.0);
  VecType dummy_residual(probDesc->solVecInfo(), 0.0);
  probDesc->getStiffAndForce(*geomState, dummy_residual, elementInternalForce, ts.locTimeOG);
  probDesc->reBuildKonly();
  delete geomState;
}

void NLDistrTimeDecompSolver::coarseIntegrator(NLTimeSlice & timeSlice)
{
  probDesc->pitaTimers.start("Coarse Time-Grid Integration");
  if (verboseFlag)
  {
    fprintf(stderr, "\n**** Begin coarse time-grid integration for TS %d - Initial time %e ****\n", timeSlice.sliceRank, timeSlice.initialTime);
  }
  probDesc->defaultPostProcessor().sliceRank(timeSlice.sliceRank);
  initializeSeqIntegrator(timeSlice, probDesc->getCoarseDt());
  probDesc->pitaTimers.start("Time Loop");
  seqIntegrator.integrate(1);
  probDesc->pitaTimers.stop();
  probDesc->defaultPostProcessor().sliceRank(-1);
  getIntegratorState(timeSlice.nextSeedState);
  probDesc->pitaTimers.stop(); // End Coarse Grid Time Integration
}

void NLDistrTimeDecompSolver::fineIntegrator(NLTimeSlice & timeSlice)
{
  probDesc->pitaTimers.start("Fine Time-Grid Integration");
  if (verboseFlag)
  {
    fprintf(stderr, "\n**** Begin fine time-grid integration for TS %d - Initial time %e ****\n", timeSlice.sliceRank, timeSlice.initialTime);
  }

  timeSlice.localBase.clear();
  probDesc->defaultPostProcessor().sliceRank(timeSlice.sliceRank);
  initializeSeqIntegrator(timeSlice, probDesc->getDt());
  for (int timeStepNum = 0; timeStepNum < probDesc->getJratio(); ++timeStepNum)
  {
    probDesc->pitaTimers.start("Time Loop");
    seqIntegrator.integrate(1);
    probDesc->pitaTimers.swap("Base Propagation").collapseIterations(true);
    linIntegrator.stepLinearizedIntegrate(timeSlice.propBase);
    if (probDesc->getBasisImprovementMethod() == 2)
    {
      probDesc->pitaTimers.swap("Save Increments");
      timeSlice.localBase.addState();
      getIntegratorState(timeSlice.localBase[timeStepNum]);
    }
    probDesc->pitaTimers.stop();
  }
  probDesc->defaultPostProcessor().sliceRank(-1);
  getIntegratorState(timeSlice.propState);
  probDesc->pitaTimers.stop(); // End Fine Grid Integration
}

} /* end namespace Old */ } /* end namespace Pita */
