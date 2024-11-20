#include "DistrElementSamplingDriver.h"
#include "SubElementSamplingDriver.h"
#include "RenumberingUtils.h"
#include "MeshDesc.h"
#include "DistrDomainUtils.h"

#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include "DistrMasterMapping.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrVecNodeDof6Conversion.h"
#include "PtrPtrIterAdapter.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include "BasisFileStream.h"
#include "VecBasisFile.h"
#include "BasisOps.h"
#include <Driver.d/SysState.h>

#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <memory>
#include <cassert>
#include <iostream>

extern GeoSource *geoSource;
extern int verboseFlag;

namespace Rom {

// Forward declarations
// ====================
int snapSize(BasisId::Type type, std::vector<int> &snapshotCounts, int j=-1, int podBasisSize = 0);

// Non-member functions
// ====================
std::string getMeshFilename(const FileNameInfo &fileInfo);

void readAndProjectSnapshots(BasisId::Type type, const DistrInfo &vectorSize, DistrNodeDof6Buffer &buffer,
                             DistrVecNodeDof6Conversion &converter, DistrVecBasis &podBasis,
                             std::vector<int> &snapshotCounts, std::vector<double> &timeStamps, DistrVecBasis *&config)
{
#ifdef PRINT_ESTIMERS
  double t1 = getTime();
#endif
  const int snapshotCount = snapSize(type, snapshotCounts, -1, podBasis.vectorCount());
  filePrint(stderr, " ... Reading in and Projecting %d %s Snapshots ...\n", snapshotCount, toString(type).c_str());

  config = new DistrVecBasis(snapshotCount, vectorSize);
  timeStamps.clear();
  timeStamps.reserve(snapshotCount);

  const int skipFactor = std::max(domain->solInfo().skipPodRom, 1); // skipFactor must be >= 1
  const int skipOffSet = std::max(domain->solInfo().skipOffSet, 0); // skipOffSet must be >= 0
  const int podVectorCount = podBasis.vectorCount();
  DistrVector snapshot(vectorSize);
  Vector podComponents(podVectorCount);
  FileNameInfo fileInfo;

  //this ain't pretty
  //allocate vector with length equal to total number of vectors
  std::vector<int> randRead;
  {
    std::string fileName = BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, 0);
    DistrBasisInputFile in(fileName);
    int numberOfVectorsInFile = in.stateCount();
    randRead.resize(numberOfVectorsInFile);
  }
  std::fill(randRead.begin(),randRead.end(),0);
  // if reading in N randomly selected vectors, call this portion and shuffle 
  if(domain->solInfo().randomVecSampling){
   // start at the offset, then fill in
   int numberOfVectors = std::max(domain->solInfo().skipPodRom,podVectorCount);// need at least as many vectors as the size of the subspace because reasons
   for(int setRead = 0+skipOffSet; setRead < numberOfVectors; setRead++)
    randRead[setRead] = 1;

   std::srand(podVectorCount); // use same seed for each file to get consistent random shuffles
   std::random_shuffle(randRead.begin()+skipOffSet, randRead.end());
  }


  int offset = 0;
  for(int i = 0; i < FileNameInfo::size(type, BasisId::SNAPSHOTS); i++) {
    std::string fileName = BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, i);
    filePrint(stderr, " ... Processing File: %s ...\n", fileName.c_str());
    DistrBasisInputFile in(fileName);

    long int yolo = 0;
    int count = 0;
    int skipCounter = 0;
    if(!domain->solInfo().randomVecSampling)
      skipCounter = skipFactor - skipOffSet;

    while(count < snapshotCounts[i]) {
      if((skipCounter == skipFactor) || (domain->solInfo().randomVecSampling && randRead[yolo+skipOffSet])) {
        assert(in.validCurrentState());
        in.currentStateBuffer(buffer);
        converter.vector(buffer, snapshot);
        expand(podBasis, reduce(podBasis, snapshot, podComponents), (*config)[offset+count]); // do projection
        timeStamps.push_back(in.currentStateHeaderValue());
        if(!domain->solInfo().randomVecSampling) skipCounter = 1;
        ++count;
        filePrint(stderr, "\r ... timeStamp = %8.2e, %3d%% done ...", in.currentStateHeaderValue(), (count*100)/snapshotCounts[i]);
      }
      else if(!domain->solInfo().randomVecSampling){
        ++skipCounter;
      }
      in.currentStateIndexInc();
      yolo++;
    }
    filePrint(stderr, "\r ... timeStamp = %8.2e, %3d%% done... \n", in.currentStateHeaderValue(), 100);

    offset += snapshotCounts[i];
  }

  assert(timeStamps.size() == snapshotCount);
#ifdef PRINT_ESTIMERS
  filePrint(stderr, "time for readAndProjectSnapshots = %f\n", (getTime()-t1)/1000.0);
#endif
}

// Member functions
// ================
const DistrInfo&
DistrElementSamplingDriver::vectorSize() const
{
  return decDomain->masterSolVecInfo();
}

DistrElementSamplingDriver::DistrElementSamplingDriver(Domain *domain, Communicator *comm) :
  MultiDomainDynam(domain),
  comm_(comm),
  solver_(NULL)
{}

void
DistrElementSamplingDriver::solve()
{
  MultiDomainDynam::preProcess();
  if(domain->solInfo().newmarkBeta == 0) {
    decDomain->assembleNodalInertiaTensors(melArray);
  }

  FileNameInfo fileInfo;

  DistrVecBasis podBasis;
  std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
  fileName.append(".orthonormalized");
  DistrBasisInputFile podBasisFile(fileName);

  const int projectionSubspaceSize = domain->solInfo().maxSizePodRom ?
                                     std::min(domain->solInfo().maxSizePodRom, podBasisFile.stateCount()) :
                                     podBasisFile.stateCount();

  podBasis.dimensionIs(projectionSubspaceSize, vectorSize());

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());

  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
  DistrNodeDof6Buffer buffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());

  for (DistrVecBasis::iterator it = podBasis.begin(),
                               it_end = podBasis.end();
                               it != it_end; ++it) {
    assert(podBasisFile.validCurrentState());

    podBasisFile.currentStateBuffer(buffer);
    converter.vector(buffer, *it);

    podBasisFile.currentStateIndexInc();
  }

  const int podVectorCount = podBasis.vectorCount();

  // Read some displacement snapshots from one or more files and project them on to the basis
  DistrVecBasis *displac; 
  std::vector<double> timeStamps;
  std::vector<int> snapshotCounts;
  readAndProjectSnapshots(BasisId::STATE, vectorSize(), buffer, converter, podBasis,
                          snapshotCounts, timeStamps, displac);

  // Optionally, read some velocity snapshots and project them on to the reduced order basis
  DistrVecBasis *veloc = 0;
  if(!domain->solInfo().velocPodRomFile.empty()) {
    std::vector<double> velTimeStamps;
    std::vector<int> velSnapshotCounts;
    readAndProjectSnapshots(BasisId::VELOCITY, vectorSize(), buffer, converter, podBasis,
                            velSnapshotCounts, velTimeStamps, veloc);
    if(velSnapshotCounts != snapshotCounts) std::cerr << " *** WARNING: inconsistent velocity snapshots\n";
  }

  // Optionally, read some acceleration snapshots and project them on to the reduced order basis
  DistrVecBasis *accel = 0;
  if(!domain->solInfo().accelPodRomFile.empty()) {
    std::vector<double> accTimeStamps;
    std::vector<int> accSnapshotCounts;
    readAndProjectSnapshots(BasisId::ACCELERATION, vectorSize(), buffer, converter, podBasis,
                            accSnapshotCounts, accTimeStamps, accel);
    if(accSnapshotCounts != snapshotCounts) std::cerr << " *** WARNING: inconsistent acceleration snapshots\n";
  }

  const int snapshotCount = std::accumulate(snapshotCounts.begin(), snapshotCounts.end(), 0);

  // read in mass-normalized basis
  if(domain->solInfo().newmarkBeta == 0 || domain->solInfo().useMassNormalizedBasis) {
    std::string normalizedBasisFileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
    normalizedBasisFileName.append(".massorthonormalized");
    DistrBasisInputFile normalizedBasisFile(normalizedBasisFileName);
    for(DistrVecBasis::iterator it = podBasis.begin(), it_end = podBasis.end(); it != it_end; ++it) {
      assert(normalizedBasisFile.validCurrentState());
      normalizedBasisFile.currentStateBuffer(buffer);
      converter.vector(buffer, *it);
      normalizedBasisFile.currentStateIndexInc();
    }
  }

  // read in ndscfg coordinates
  DistrVecBasis *ndscfgCoords = 0;
  if(domain->solInfo().NodeTrainingFiles.size() > 0) {
    DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub(), false);
    ndscfgCoords = new DistrVecBasis[domain->solInfo().NodeTrainingFiles.size()];
    DistrInfo coordVecInfo;
    decDomain->makeBasicDistrInfo(coordVecInfo, &Domain::coordVecSize);
    for(int i = 0; i < domain->solInfo().NodeTrainingFiles.size(); ++i) {
      DistrBasisInputFile in(domain->solInfo().NodeTrainingFiles[i]);
      ndscfgCoords[i].dimensionIs(1, coordVecInfo);
      in.currentStateBuffer(buffer);
      converter.vector(buffer, *ndscfgCoords[i].begin());
    }
  }

  // read in ndscfg mass orthogonal bases
  DistrVecBasis *ndscfgMassOrthogonalBases = 0;
  if(domain->solInfo().MassOrthogonalBasisFiles.size() > 0) {
    ndscfgMassOrthogonalBases = new DistrVecBasis[domain->solInfo().MassOrthogonalBasisFiles.size()];
    for(int i = 0; i < domain->solInfo().MassOrthogonalBasisFiles.size(); ++i) {
      DistrBasisInputFile in(domain->solInfo().MassOrthogonalBasisFiles[i]);
      ndscfgMassOrthogonalBases[i].dimensionIs(projectionSubspaceSize, vectorSize());
      for(DistrVecBasis::iterator it = ndscfgMassOrthogonalBases[i].begin(), it_end = ndscfgMassOrthogonalBases[i].end(); it != it_end; ++it) {
        assert(in.validCurrentState());
        in.currentStateBuffer(buffer);
        converter.vector(buffer, *it);
        in.currentStateIndexInc();
      }
    }
  }
 
  int glNumSubs = decDomain->getNumSub();
  structCom->globalSum(1,&glNumSubs);
 
  SubElementSamplingDriver **subDrivers = new SubElementSamplingDriver * [decDomain->getNumSub()];
  SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t> **subSolvers = new SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t> * [decDomain->getNumSub()];
  Vector *solutions = new Vector[decDomain->getNumSub()];
  int numCPUs = (structCom) ? structCom->numCPUs() : 1;
  int myID = (structCom) ? structCom->myID() : 0;
  solver_ = new ParallelSparseNonNegativeLeastSquaresSolver(decDomain->getNumSub(), subSolvers);

  if(domain->solInfo().stackedElementSampling){
    filePrint(stderr," ... ALLOCATING EXTRA MEMORY FOR STACKED TRAINING ...\n");
    solver_->problemSizeIs(2*podVectorCount*snapshotCount, domain->numElements());
  } else {
    solver_->problemSizeIs(podVectorCount*snapshotCount, domain->numElements());
  }

#ifdef PRINT_ESTIMERS
  double t2 = getTime();
#endif
  StackVector glTrainingTarget(solver_->rhsBuffer(), solver_->equationCount());
  glTrainingTarget.zero();
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
  for(int i = 0; i < decDomain->getNumSub(); ++i) {
    subDrivers[i] = new SubElementSamplingDriver(decDomain->getAllSubDomains()[i]);
    subSolvers[i] = &subDrivers[i]->solver();

    std::vector<StackVector> subPodBasis(podVectorCount);
    for(int j = 0; j < podVectorCount; ++j) {
      subPodBasis[j].setData(podBasis[j].subData(i), podBasis[j].subLen(i));
    }

    std::vector<StackVector> subDisplac(snapshotCount);
    for(int j = 0; j < snapshotCount; ++j) {
      subDisplac[j].setData((*displac)[j].subData(i), (*displac)[j].subLen(i));
    }

    subDrivers[i]->timeStampsIs(timeStamps);
    subDrivers[i]->snapshotCountsIs(snapshotCounts);

    std::vector<StackVector> *subVeloc = 0;
    if(veloc) {
      subVeloc = new std::vector<StackVector>(snapshotCount);
      for(int j = 0; j < snapshotCount; ++j) {
        (*subVeloc)[j].setData((*veloc)[j].subData(i), (*veloc)[j].subLen(i));
      }
    }

    std::vector<StackVector> *subAccel = 0;
    if(accel) {
      subAccel = new std::vector<StackVector>(snapshotCount); 
      for(int j = 0; j < snapshotCount; ++j) {
        (*subAccel)[j].setData((*accel)[j].subData(i), (*accel)[j].subLen(i));
      }
    }

    std::vector<StackVector> *subNdscfgCoords = 0;
    if(ndscfgCoords) {
      subNdscfgCoords = new std::vector<StackVector>[domain->solInfo().NodeTrainingFiles.size()];
      for(int j = 0; j < domain->solInfo().NodeTrainingFiles.size(); ++j) {
        subNdscfgCoords[j].resize(1);
        subNdscfgCoords[j][0].setData(ndscfgCoords[j][0].subData(i), ndscfgCoords[j][0].subLen(i));
      }
    }

    std::vector<StackVector> *subNdscfgMassOrthogonalBases = 0;
    if(ndscfgMassOrthogonalBases) {
      subNdscfgMassOrthogonalBases = new std::vector<StackVector>[domain->solInfo().MassOrthogonalBasisFiles.size()];
      for(int j = 0; j < domain->solInfo().MassOrthogonalBasisFiles.size(); ++j) {
        subNdscfgMassOrthogonalBases[j].resize(podVectorCount);
        for(int k = 0; k < podVectorCount; ++k) {
          subNdscfgMassOrthogonalBases[j][k].setData(ndscfgMassOrthogonalBases[j][k].subData(i), ndscfgMassOrthogonalBases[j][k].subLen(i));
        }
      }
    }

    subDrivers[i]->preProcess();
    if(domain->solInfo().stackedElementSampling)
      subDrivers[i]->solver().problemSizeIs(2*podVectorCount*snapshotCount, subDrivers[i]->elementCount());
    else
      subDrivers[i]->solver().problemSizeIs(podVectorCount*snapshotCount, subDrivers[i]->elementCount());

    for(int j=0; j<subDrivers[i]->solver().equationCount(); ++j) subDrivers[i]->solver().rhsBuffer()[j] = 0.0;
    subDrivers[i]->assembleTrainingData(subPodBasis, podVectorCount, subDisplac, subVeloc, subAccel,
                                        subNdscfgCoords, subNdscfgMassOrthogonalBases);
    subDrivers[i]->clean();
    if(subVeloc) delete subVeloc;
    if(subAccel) delete subAccel;
    if(subNdscfgCoords) delete [] subNdscfgCoords;
    if(subNdscfgMassOrthogonalBases) delete [] subNdscfgMassOrthogonalBases;

    StackVector *trainingTarget;
    if(domain->solInfo().stackedElementSampling)
      trainingTarget = new StackVector(subDrivers[i]->solver().rhsBuffer(), 2*podVectorCount*snapshotCount);
    else
      trainingTarget = new StackVector(subDrivers[i]->solver().rhsBuffer(), podVectorCount*snapshotCount);
#if defined(_OPENMP)
    #pragma omp critical
#endif
    glTrainingTarget += *trainingTarget;
  }
  delete displac;
  if(veloc) delete veloc;
  if(accel) delete accel;
  if(ndscfgCoords) delete [] ndscfgCoords;
  if(ndscfgMassOrthogonalBases) delete [] ndscfgMassOrthogonalBases;

  if(structCom) 
    structCom->globalSum(glTrainingTarget.size(), glTrainingTarget.data()); 
#ifdef PRINT_ESTIMERS
  filePrint(stderr, "time for assembleTrainingData = %f\n", (getTime()-t2)/1000.0);
#endif

  computeSolution(solutions, domain->solInfo().tolPodRom);

  std::vector<double> lweights; 
  std::vector<int> lelemIds;
  for(int i = 0; i < decDomain->getNumSub(); i++) {
    subDrivers[i]->getGlobalWeights(solutions[i], lweights, lelemIds, verboseFlag);
    delete subDrivers[i];
  }
  delete [] solutions;
  delete [] subSolvers;
  delete [] subDrivers;
  
  std::vector<double> gweights(domain->numElements());
  std::vector<int> gelemIds(domain->numElements());
  
  // Gather weights and IDs from all processors
  int numLocalElems = lweights.size(); 
  if(structCom) {
    int recvcnts[numCPUs];
    int displacements[numCPUs];
    structCom->allGather(&numLocalElems, 1, &recvcnts[0], 1);
    int location = 0;
    for(int i = 0; i < numCPUs; i++) {
      displacements[i] = location;
      location += recvcnts[i];
    }
    structCom->gatherv(&lweights[0], lweights.size(), &gweights[0], &recvcnts[0], &displacements[0], 0);
    structCom->gatherv(&lelemIds[0], lelemIds.size(), &gelemIds[0], &recvcnts[0], &displacements[0], 0);
  }
  else {
    gweights = lweights;
    gelemIds = lelemIds;
  }

  // Compute the reduced forces
  DistrVector forceFull(decDomain->masterSolVecInfo());
  GenAssembler<double> * assembler = decDomain->getSolVecAssembler();
  // 1) gravity
  Vector gravForceRed(podBasis.vectorCount());
  if(domain->gravityFlag()) {
    MultiDomainDynam::getGravityForce(forceFull);
    assembler->assemble(forceFull);
    reduce(podBasis, forceFull, gravForceRed);
  }
  // 2) constant force or constant part of time-dependent forces (default loadset only) TODO add support for multiple loadsets
  MultiDomainDynam::getUnamplifiedExtForce(forceFull, 0);
  assembler->assemble(forceFull);
  Vector constForceRed(podBasis.vectorCount());
  bool reduce_f = (forceFull.norm() != 0);
  if(reduce_f) reduce(podBasis, forceFull, constForceRed);

  // compute the reduced initial conditions
  DistrVector d0Full(MultiDomainDynam::solVecInfo()),
              v0Full(MultiDomainDynam::solVecInfo());
  DistrVector tmp(decDomain->masterSolVecInfo());
  SysState<DistrVector> inState(d0Full, v0Full, tmp, tmp);
  MultiDomainDynam::getInitState(inState);
  Vector d0Red(podBasis.vectorCount()),
         v0Red(podBasis.vectorCount());
  bool reduce_idis = (d0Full.norm() != 0),
       reduce_ivel = (v0Full.norm() != 0);
  if(domain->solInfo().useMassNormalizedBasis || domain->solInfo().newmarkBeta == 0) {
    SubDOp *M = NULL;
    if(reduce_idis || reduce_ivel) {
      SparseMatrix **subM = new SparseMatrix * [decDomain->getNumSub()];
      execParal(decDomain->getNumSub(), this, &DistrElementSamplingDriver::subMakeMass, subM);
      M = new SubDOp(decDomain->getNumSub(), subM);
    }
    if(reduce_idis) {
      M->mult(d0Full, tmp);
      assembler->assemble(tmp);
      reduce(podBasis, tmp, d0Red);
    }
    if(reduce_ivel) {
      M->mult(v0Full, tmp);
      assembler->assemble(tmp);
      reduce(podBasis, tmp, v0Red);
    }
    if(M) delete M;
  }
  else {
    if(reduce_idis) reduce(podBasis, d0Full, d0Red);
    if(reduce_ivel) reduce(podBasis, v0Full, v0Red);
  }

  if(myID == 0) {

    // Weights output file generation
    outputFullWeights(gweights, gelemIds);

    // Mesh output file generation
    std::map<int,double> weightsMap;
    std::vector<int> reducedelemIds;
    for(int i = 0; i < gweights.size(); i++) {
      if(gweights[i] > 0) {
	reducedelemIds.push_back(gelemIds[i]);
        weightsMap.insert(std::pair<int,double>(gelemIds[i], gweights[i]));
      }
    }

    addContactElems(reducedelemIds, weightsMap);

    const FileNameInfo fileInfo;
    for(int i = 0; i < decDomain->getNumSub(); i++) {
      decDomain->getSubDomain(i)->renumberElementsGlobal();
    }

    // Read in the truncated basis into a (non-distributed) VecBasis XXX
    domain->preProcessing();
    buildDomainCdsa();
    std::string fileName2 = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
    if(domain->solInfo().newmarkBeta == 0 || domain->solInfo().useMassNormalizedBasis) fileName2.append(".massorthonormalized");
    else fileName2.append(".orthonormalized");
    const VecNodeDof6Conversion vecDofConversion(*domain->getCDSA());
    BasisInputStream<6> in(fileName2, vecDofConversion);
    VecBasis podBasis;
    const int podSizeMax = domain->solInfo().maxSizePodRom;
    if(podSizeMax != 0) {
      readVectors(in, podBasis, podSizeMax);
    } else {
      readVectors(in, podBasis);
    }

    // Output the reduced mesh
    Elemset &inputElemSet = *(geoSource->getElemSet());
    std::unique_ptr<Connectivity> elemToNode(new Connectivity(inputElemSet.asSet()));
    const MeshRenumbering meshRenumbering(reducedelemIds.begin(), reducedelemIds.end(), *elemToNode, verboseFlag);
    const MeshDesc reducedMesh(domain, geoSource, meshRenumbering, weightsMap); 
    outputMeshFile(fileInfo, reducedMesh, podBasis.vectorCount());

    // Output the reduced forces
    std::ofstream meshOut(getMeshFilename(fileInfo).c_str(), std::ios_base::app);
    if(domain->solInfo().reduceFollower) meshOut << "*\nEXTFOL\n";
    if(domain->gravityFlag()) {
      meshOut << "*\nFORCES -1\nMODAL\n"; // note: gravity forces are put in loadset -1 so that MFTT (if present) will not be applied
      meshOut.precision(std::numeric_limits<double>::digits10+1);
      for(int i = 0; i < podBasis.vectorCount(); ++i)
        meshOut << i+1 << " " << gravForceRed[i] << std::endl;
    }
    if(reduce_f) {
      meshOut << "*\nFORCES\nMODAL\n";
      meshOut.precision(std::numeric_limits<double>::digits10+1);
      for(int i = 0; i < podBasis.vectorCount(); ++i)
        meshOut << i+1 << " " << constForceRed[i] << std::endl;
    }

    // output the reduced initial conditions
    if(reduce_idis) {
      meshOut << "*\nIDISPLACEMENTS\nMODAL\n";
      meshOut.precision(std::numeric_limits<double>::digits10+1);
      for(int i=0; i<podBasis.vectorCount(); ++i)
        meshOut << i+1 << " " << d0Red[i] << std::endl;
    }
    if(reduce_ivel) {
      meshOut << "*\nIVELOCITIES\nMODAL\n";
      meshOut.precision(std::numeric_limits<double>::digits10+1);
      for(int i=0; i<podBasis.vectorCount(); ++i)
        meshOut << i+1 << " " << v0Red[i] << std::endl;
    }

    meshOut.close();

#ifdef USE_EIGEN3
    // Build and output compressed basis
    DofSetArray reduced_dsa(reducedMesh.nodes().size(), const_cast<Elemset&>(reducedMesh.elements()));
    int num_bc = reducedMesh.dirichletBConds().size();
    BCond *bc = (num_bc > 0) ? const_cast<BCond*>(&reducedMesh.dirichletBConds()[0]) : NULL;
    ConstrainedDSA reduced_cdsa(reduced_dsa, num_bc, bc);
    podBasis.makeSparseBasis(meshRenumbering.reducedNodeIds(), domain->getCDSA(), &reduced_cdsa);
    {
      std::string filename = getMeshFilename(fileInfo).c_str(); 
      filename.append(".compressed.basis");
      if(domain->solInfo().newmarkBeta == 0 || domain->solInfo().useMassNormalizedBasis) filename.append(".massorthonormalized");
      else filename.append(".orthonormalized");
      filePrint(stderr," ... Writing compressed basis to file %s ...\n", filename.c_str());
      VecNodeDof6Conversion converter(reduced_cdsa);
      BasisOutputStream<6> output(filename, converter, false);

      for (int iVec = 0; iVec < podBasis.vectorCount(); ++iVec) {
        output << podBasis.compressedBasis().col(iVec);
      }
    }
#endif
  }
}

void
DistrElementSamplingDriver::addContactElems(std::vector<int> &sampleElemIds, std::map<int, double> &weights){
  // Add elements attached to Contact surfaces, needed for assigning mass to nodes in surface 
  int numMC = domain->GetnMortarConds();
  if(numMC > 0) { 
    std::set<int> activeSurfs;
    //First: loop over Mortar Constraints to identify active surfaces
    for(int mc = 0; mc<numMC; ++mc){
      if(domain->GetMortarCond(mc)->GetInteractionType() == MortarHandler::CTC){
        activeSurfs.insert(domain->GetMortarCond(mc)->GetMasterEntityId());
        activeSurfs.insert(domain->GetMortarCond(mc)->GetSlaveEntityId());
      }
    }
  
    Elemset &inputElemSet = *(geoSource->getElemSet());
    Connectivity elemToNode{ SetAccess<Elemset>{inputElemSet} };
    Connectivity nodeToElem = elemToNode.reverse();

    //Second: loop over list of surfaces to extract elements attached to surface nodes 
    for(const auto surfPtr : activeSurfs) {
      for(int surf = 0; surf < domain->getNumSurfs(); ++surf) { // loop through surfaces and check for matching ID
        if(surfPtr == domain->GetSurfaceEntity(surf)->GetId()) {
          int nNodes = domain->GetSurfaceEntity(surf)->GetnNodes();
          // loop through nodes of surface
          for(int nn = 0; nn < nNodes; ++nn) {
            int glNode = domain->GetSurfaceEntity(surf)->GetGlNodeId(nn);
            // loop through elements attached to node
            for(int j = 0; j < nodeToElem.num(glNode); ++j) {
              const int elemRank = nodeToElem[glNode][j];
              // check if element already in weighted set, if not, add to weights and reduced ElemIds; 
              if(std::find(sampleElemIds.begin(), sampleElemIds.end(), elemRank) == sampleElemIds.end()){
                sampleElemIds.push_back(elemRank);
                // now add zero weight to each local weight map
                weights.insert(std::make_pair(elemRank, 0.0));
              }
            } 
          }   
        }     
      }
    }
  }
  std::sort(sampleElemIds.begin(), sampleElemIds.end());
}

void
DistrElementSamplingDriver::computeSolution(Vector *solutions, double relativeTolerance, bool verboseFlag)
{
  solver_->relativeToleranceIs(relativeTolerance);
  solver_->verboseFlagIs(verboseFlag);
  solver_->scalingFlagIs(domain->solInfo().useScalingSpnnls);
  solver_->centerFlagIs(domain->solInfo().useCenterSpnnls);
  solver_->projectFlagIs(domain->solInfo().projectSolution);
  solver_->positivityIs(domain->solInfo().positiveElements);
  solver_->solverTypeIs(domain->solInfo().solverTypeSpnnls);
  solver_->maxSizeRatioIs(domain->solInfo().maxSizeSpnnls);
  solver_->maxIterRatioIs(domain->solInfo().maxIterSpnnls);
  solver_->npMaxIs(domain->solInfo().npMax);
  solver_->scpkMBIs(domain->solInfo().scpkMB);
  solver_->scpkNBIs(domain->solInfo().scpkNB);
  solver_->scpkMPIs(domain->solInfo().scpkMP);
  solver_->scpkNPIs(domain->solInfo().scpkNP);
  try {
    solver_->solve();
  }
  catch(std::runtime_error& e) {
    if(structCom->myID() == 0) std::cerr << " *** WARNING: " << e.what() << std::endl;
  }

  if(verboseFlag) {
/*  std::cout << "Primal solution:";
    for (int elemRank = 0; elemRank != elementCount(); ++elemRank) {
      std::cout << " " << solver_.solutionEntry(elemRank);
    }
    std::cout << "\n"; */
    StackVector trainingTarget(solver_->rhsBuffer(), solver_->equationCount());
    double glTargMagnitude = trainingTarget.norm();
    double oneNorm = 0;
    for(int i = 0; i < solver_->subdomainCount(); ++i) {
      oneNorm += std::accumulate(solver_->subdomainSolver(i)->solutionBuffer(),
                                 solver_->subdomainSolver(i)->solutionBuffer() + solver_->subdomainSolver(i)->unknownCount(), 0.0);
    }
#ifdef USE_MPI
    oneNorm = structCom->globalSum(oneNorm);
#endif
    if(!structCom || structCom->myID() == 0) {
      std::cout << "Error magnitude / Absolute tolerance = " << solver_->errorMagnitude() << " / " << solver_->relativeTolerance() * glTargMagnitude << "\n";
      std::cout << "1-norm of primal solution = " << oneNorm << "\n";
      std::cout.flush();
    }
  }

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
  // Read solution
  for(int i = 0; i < solver_->subdomainCount(); ++i) {
    solutions[i].initialize(solver_->subdomainSolver(i)->unknownCount());
    std::copy(solver_->subdomainSolver(i)->solutionBuffer(), solver_->subdomainSolver(i)->solutionBuffer() + solver_->subdomainSolver(i)->unknownCount(),
              solutions[i].data());
  }
}

void
DistrElementSamplingDriver::buildDomainCdsa()
{
  const int numdof = domain->numdof();
  SimpleBuffer<int> bc(numdof);
  SimpleBuffer<double> bcx(numdof);

  domain->make_bc(bc.array(), bcx.array());
  domain->make_constrainedDSA(bc.array());
}

void
DistrElementSamplingDriver::subMakeMass(int i, SparseMatrix **subM)
{
  AllOps<double> allOps;
  allOps.M = subM[i] = decDomain->getSubDomain(i)->constructDBSparseMatrix<double>();
  decDomain->getSubDomain(i)->makeSparseOps<double>(allOps, 0, 0, 0);
}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrElementSamplingDriverNew(Domain *domain)
{
  return new Rom::DistrElementSamplingDriver(domain, structCom);
}
