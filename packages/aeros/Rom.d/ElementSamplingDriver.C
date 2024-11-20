#include "ElementSamplingDriver.h"

#include "DistrBasisFile.h"
#include "VecBasis.h"
#include "BasisOps.h" 
#include "FileNameInfo.h"
#include "BasisFileStream.h"
#include "VecBasisFile.h"
#include "SimpleBuffer.h"
#include "RenumberingUtils.h"
#include "MeshDesc.h"

#include "VecBasisOps.h"

#include <Driver.d/GeoSource.h>
#include <Driver.d/Domain.h>
#include <Driver.d/SysState.h>
#include <Math.d/Vector.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/EiSparseMatrix.h>
#include <Math.d/DiagMatrix.h>
#include <Math.d/VectorSet.h>
#include <Timers.d/StaticTimers.h>
#include <Utils.d/Connectivity.h>
#include <Element.d/Element.h>
#include <Utils.d/Conwep.d/BlastLoading.h>
#include <Corotational.d/Corotator.h>

#include <cstddef>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <set>
#include <map>
#include <vector>
#include <memory>
#include <fstream>
#include <string>
#include <utility>

#include <cassert>
#include <iostream>

extern GeoSource *geoSource;

namespace Rom {

// Non-member functions
// ====================
template <typename Scalar>
inline
void
copy(const GenVector<Scalar> &v, Scalar *target) {
  const double *originBuf = v.data();
  std::copy(originBuf, originBuf + v.size(), target);
}

template <typename Scalar, typename OutputIterator>
inline
void
copy(const GenVector<Scalar> &v, OutputIterator target) {
  const double *originBuf = v.data();
  std::copy(originBuf, originBuf + v.size(), target);
}

std::string
getMeshFilename(const FileNameInfo &fileInfo) {
//  return fileInfo.prefix() + ".elementmesh.inc";
  std::string FName(domain->solInfo().reducedMeshFile, std::strlen(domain->solInfo().reducedMeshFile));
  return FName + ".elementmesh.inc";
}

void
outputMeshFile(const FileNameInfo &fileInfo, const MeshDesc &mesh, const int podVectorCount) {
  const std::ios_base::openmode mode = std::ios_base::out; 
  std::ofstream meshOut(getMeshFilename(fileInfo).c_str(), mode);
  filePrint(stderr," ... Writing Mesh File to %s ...\n", getMeshFilename(fileInfo).c_str());
  meshOut.precision(std::numeric_limits<double>::digits10+1);
  std::string basisfile = getMeshFilename(fileInfo).c_str();
  basisfile.append(".compressed.basis");
  meshOut << "READMODE\n";
  if(domain->solInfo().useMassNormalizedBasis || domain->solInfo().newmarkBeta == 0) {
    meshOut << "0 mnorm \"" << basisfile << ".massorthonormalized\" " << podVectorCount << "\n"; 
  }
  else {
    meshOut << "0 inorm \"" << basisfile << ".orthonormalized\" " << podVectorCount << "\n";
  }
  if (!domain->solInfo().readInDualROB.empty()) {
    std::string dualBasisfile = getMeshFilename(fileInfo).c_str();
    dualBasisfile.append(".compressed.dualbasis");
    meshOut << "1 noneg \"" << dualBasisfile << "\" " << domain->solInfo().localDualBasisSize[0] << "\n";
  }
  meshOut << "*\n";
  meshOut << mesh;
}

void
outputMeshFile(const FileNameInfo &fileInfo, const MeshDesc &mesh, const std::vector<int> &localBasisSize) {
  const std::ios_base::openmode mode = std::ios_base::out;
  std::ofstream meshOut(getMeshFilename(fileInfo).c_str(), mode);
  filePrint(stderr," ... Writing Mesh File to %s ...\n", getMeshFilename(fileInfo).c_str());
  meshOut.precision(std::numeric_limits<double>::digits10+1);
  std::string basisfile = getMeshFilename(fileInfo).c_str();
  basisfile.append(".compressed.basis");
  meshOut << "READMODE\n";
  for(int j=0; j<localBasisSize.size(); ++j) {
    if(domain->solInfo().useMassNormalizedBasis || domain->solInfo().newmarkBeta == 0) {
      meshOut << j+1 << " mnorm \"" << basisfile << (j+1) << ".massorthonormalized\" " << localBasisSize[j] << "\n";                
    }
    else {
      meshOut << j+1 << " inorm \"" << basisfile << (j+1) << ".orthonormalized\" " << localBasisSize[j] << "\n";
    }
  }
  meshOut << "*\n";
  meshOut << mesh;
}


template<typename WeightsVecType, typename ElemIdsVecType>
void
outputFullWeights(const WeightsVecType &weights, const ElemIdsVecType &elemIds, int j)
{
  assert(weights.size() == elemIds.size());

  std::string fileName = domain->solInfo().reducedMeshFile;
  if(j > -1) {
    std::ostringstream ss;
    ss << ".cluster" << j+1;
    fileName.append(ss.str());
  }
  const std::ios_base::openmode mode = std::ios_base::out;
  std::ofstream weightOut(fileName.c_str(), mode);
  weightOut.precision(std::numeric_limits<double>::digits10+1);
  bool firstTime = true;

  std::map<int, Attrib> &attrib = geoSource->getAttributes();
  int na = geoSource->getNumAttributes();
  int nMaxEle = geoSource->getElemSet()->last();

  std::map<int, Attrib>::iterator *elemAttrib = new std::map<int, Attrib>::iterator[nMaxEle];
  for(int i = 0; i < nMaxEle; ++i) elemAttrib[i] = attrib.end();
  for (std::map<int, Attrib>::iterator it = attrib.begin(); it != attrib.end(); ++it) {
    if(it->second.nele < nMaxEle)
      elemAttrib[it->second.nele] = it;
  }

  weightOut << "ATTRIBUTES";
  if(j >= 0) weightOut << " " << j+1;
  weightOut << "\n";
  for (int i = 0, iEnd = weights.size(); i != iEnd; ++i) {
    if(elemAttrib[elemIds[i]] != attrib.end() && j < 0) {
      // element has an attribute
      Attrib &a = elemAttrib[elemIds[i]]->second;
      weightOut << elemIds[i]+1 << " " << a.attr+1 << " ";
      if(a.cmp_attr >= 0) {
        if(a.cmp_frm >= 0) weightOut << a.cmp_attr+1 << " " << a.cmp_frm+1 << " ";
        else  weightOut << a.cmp_attr+1 << " THETA " << a.cmp_theta << " ";
      }
    }
    else {
      // element has no attribute, however we still need to define the weight for it
      weightOut << elemIds[i]+1 << " ";
    }
    if(domain->solInfo().reduceFollower && firstTime) {
      weightOut << "HRC " << weights[i] << " EXTFOL\n";
      firstTime = false;
    }
    else {
      weightOut << "HRC " << weights[i] << "\n";
    }
  }

  delete [] elemAttrib;
  weightOut.close();
}

int snapSize(BasisId::Type type, std::vector<int> &snapshotCounts, int j=-1, int podBasisSize = 0)
{
  // compute the number of snapshots that will be used for a given skipping strategy
  FileNameInfo fileInfo;
  snapshotCounts.clear();
  int skipFactor = std::max(domain->solInfo().skipPodRom, 1); // skipFactor must be >= 1
  const int skipOffSet = std::max(domain->solInfo().skipOffSet, 0); // skipOffSet must be >= 0
  for(int i = 0; i < FileNameInfo::size(type, BasisId::SNAPSHOTS); i++) {
    std::string fileName = BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, i);
    DistrBasisInputFile in(fileName);
    if(domain->solInfo().randomVecSampling){ // if random sampling, do this
      //const int singleSnapshotCount = std::max(domain->solInfo().randomSampleSize, podBasisSize);
      const int singleSnapshotCount = domain->solInfo().randomSampleSize;
      snapshotCounts.push_back(singleSnapshotCount);
    } else { // otherwise compute this stuff
      const double N = in.stateCount() - skipOffSet;
      const int singleSnapshotCount = (N > 0) ? 1+(N-1)/skipFactor : 0;
      snapshotCounts.push_back(singleSnapshotCount);
    }
  }
  // return the total
  return (j<0) ? std::accumulate(snapshotCounts.begin(), snapshotCounts.end(), 0) : snapshotCounts[j];
}

void readAndProjectSnapshots(BasisId::Type type, const int vectorSize, VecBasis &podBasis,
                             const VecNodeDof6Conversion &vecDofConversion,
                             std::vector<int> &snapshotCounts, std::vector<double> &timeStamps, VecBasis &config,
                             SparseMatrix *M, int j = -1)
{
  // if j == -1 then read all of the snapshot files
  // if j >= 0  the read only the (j+1)-th snapshot file
#ifdef PRINT_ESTIMERS
  double t1 = getTime();
#endif
  const int snapshotCount = snapSize(type, snapshotCounts, j, podBasis.vectorCount());
  filePrint(stderr, " ... Reading in and Projecting %d %s Snapshots ...\n", snapshotCount, toString(type).c_str());

  config.dimensionIs(snapshotCount, vectorSize);
  timeStamps.clear();
  timeStamps.reserve(snapshotCount);
  
  const int skipFactor = std::max(domain->solInfo().skipPodRom, 1); // skipFactor must be >= 1
  const int skipOffSet = std::max(domain->solInfo().skipOffSet, 0); // skipOffSet must be >= 0
  const int podVectorCount = podBasis.vectorCount();
  Vector snapshot(vectorSize), Msnapshot(vectorSize);
  Vector podComponents(podVectorCount);
  const FileNameInfo fileInfo;

  // this ain't pretty
  // allocate vector with length equal to total number of vectors 
  std::vector<int> randRead;
  {
    std::string fileName = BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, 0);
    BasisInputStream<6> in(fileName, vecDofConversion);
    int numberOfVectorsInFile = in.size();
    randRead.resize(numberOfVectorsInFile);
  }
  std::fill(randRead.begin(),randRead.end(),0);
  // if reading in N randomly selected vectors, call this portion and shuffle 
  if(domain->solInfo().randomVecSampling){
   // start at the offset, then fill in
   //int numberOfVectors = std::max(domain->solInfo().randomSampleSize,podVectorCount);// need at least as many vectors as the size of the subspace because reasons
   int numberOfVectors = domain->solInfo().randomSampleSize;
   for(int setRead = 0+skipOffSet; setRead < numberOfVectors; setRead++)
    randRead[setRead] = 1;

   std::srand(podVectorCount); // use same seed for each file to get consistent vectors
   std::random_shuffle(randRead.begin()+skipOffSet, randRead.end());
  }

  int offset = 0;
  for(int i = ((j<0)?0:j); i < ((j<0)?FileNameInfo::size(type, BasisId::SNAPSHOTS):(j+1)); i++) {
    std::string fileName = BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, i);
    filePrint(stderr, " ... Processing File: %s ...\n", fileName.c_str());
    BasisInputStream<6> in(fileName, vecDofConversion);

    double s0 = -getTime(), s1 = -51, s2 = 0;
    int count = 0;
    int skipCounter = 0;
    if(!domain->solInfo().randomVecSampling) // don't use skipCounter if using randomly selected vectors
      skipCounter = skipFactor - skipOffSet;
    std::pair<double, double *> data;
   
    long int yolo = 0; 
    while(count < snapshotCounts[i]) {
      // if doing random sampling, ignore when skipCounter equals skipFactor, only read when randRead is 1
      if((skipCounter == skipFactor) || (domain->solInfo().randomVecSampling && randRead[yolo+skipOffSet])) {
        data.second = snapshot.data();
        in >> data;
        assert(in);
        config[offset+count] = snapshot.data();
        if(domain->solInfo().useMassOrthogonalProjection) {
          M->mult(snapshot, Msnapshot);
          expand(podBasis, reduce(podBasis, Msnapshot, podComponents), config[offset+count]);
        }
        else {
          expand(podBasis, reduce(podBasis, snapshot, podComponents), config[offset+count]);
        }
        timeStamps.push_back(data.first);
        if(!domain->solInfo().randomVecSampling)
          skipCounter = 1;
        ++count;
        if((s2-s1 > 50)) { // only print to the screen every 50 milliseconds, otherwise it's too slow... (Noyse)
          s1 = s2;
          filePrint(stderr, "\r ... timeStamp = %8.2e, %3d%% done ...", data.first, (count*100)/snapshotCounts[i]);
        }
        s2 = s0+getTime();
      }
      else {
        in.file().currentStateIndexInc();
        if(!domain->solInfo().randomVecSampling)
          ++skipCounter;
      }
      ++yolo;
    }

    filePrint(stderr, "\r ... timeStamp = %8.2e, %3d%% done... \n", data.first, 100);
    offset += snapshotCounts[i];
  }
  
  assert(timeStamps.size() == snapshotCount);
#ifdef PRINT_ESTIMERS
  fprintf(stderr, "time for readAndProjectSnapshots = %f\n", (getTime()-t1)/1000.0);
#endif
}

// Member functions
// ================
template<typename MatrixBufferType, typename SizeType>
int
ElementSamplingDriver<MatrixBufferType,SizeType>::elementCount() const {
  return domain_->numElements();
}

template<typename MatrixBufferType, typename SizeType>
int
ElementSamplingDriver<MatrixBufferType,SizeType>::vectorSize() const {
  return domain_->numUncon();
}

template<typename MatrixBufferType, typename SizeType>
ElementSamplingDriver<MatrixBufferType,SizeType>::ElementSamplingDriver(Domain *d) :
  SingleDomainDynamic(d),
  domain_(d),
  corotators_(NULL),
  geomState_(NULL),
  kelArray_(NULL),
  melArray_(NULL),
  veloc_(NULL),
  accel_(NULL),
  ndscfgCoords_(NULL),
  ndscfgMassOrthogonalBases_(NULL)
{}

template<typename MatrixBufferType, typename SizeType>
ElementSamplingDriver<MatrixBufferType,SizeType>::~ElementSamplingDriver() {
  clean();
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::clean() {
  if (corotators_) {
    for (int iElem = 0; iElem != elementCount(); ++iElem) {
      const Corotator * c = corotators_[iElem];
      if (!dynamic_cast<const Element *>(c)) {
        delete corotators_[iElem];
      }
    }
    delete[] corotators_;
    corotators_ = NULL;
  }

  if(geomState_) { delete geomState_; geomState_ = NULL; }
  if(kelArray_) { delete [] kelArray_; kelArray_ = NULL; }
  if(melArray_) { delete [] melArray_; melArray_ = NULL; }

  if(veloc_) { delete veloc_; veloc_ = NULL; }
  if(accel_) { delete accel_; accel_ = NULL; }

  if(ndscfgCoords_) { delete [] ndscfgCoords_; ndscfgCoords_ = 0; }
  if(ndscfgMassOrthogonalBases_) { delete [] ndscfgMassOrthogonalBases_; ndscfgMassOrthogonalBases_ = 0; }

  timeStamps_.clear();
}

template<typename MatrixBufferType, typename SizeType>
template<typename VecBasisType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>
::assembleTrainingData(VecBasisType &podBasisInput, int podVectorCount, VecBasisType &displac,
                       VecBasisType *veloc, VecBasisType *accel, VecBasisType *ndscfgCoords,
                       VecBasisType *ndscfgMassOrthogonalBases, int j)
{
  const FileNameInfo fileInfo;
  std::vector<double>::iterator timeStampFirst = timeStamps_.begin();
  typename MatrixBufferType::iterator elemContributions = solver_.matrixBuffer();
  double *trainingTarget = solver_.rhsBuffer();
  VecBasisType *podBasis = &podBasisInput;

  // Temporary buffers shared by all iterations
  Vector elemTarget(podVectorCount);
  Vector *elemTarget2;
  SimpleBuffer<double> *elementForce2;
  if(domain_->solInfo().stackedElementSampling){
    elemTarget2   = new Vector(podVectorCount);
    elementForce2 = new SimpleBuffer<double>(domain_->maxNumDOF());
  }
  SimpleBuffer<double> elementForce(domain_->maxNumDOF());

  // load CONWEP configuration data
  BlastLoading::BlastData *conwep;
  if(domain_->solInfo().conwepConfigurations.empty()) conwep = (domain_->solInfo().ConwepOnOff) ? &BlastLoading::InputFileData : NULL;
  else if(domain_->solInfo().conwepConfigurations.size() < domain_->solInfo().statePodRomFile.size()) {
    filePrint(stderr, " *** ERROR: Must provide one Conwep configuration per state snapshot file\n");
    exit(-1);
  }
  domain_->makeElementAdjacencyLists();

  // set snapshot sampling parameters
  int kParam, kSnap, jParam = -1;
  const int skipFactor  = std::max(domain->solInfo().skipPodRom, 1); // skipFactor must be >= 1
  const int skipOffSet  = std::max(domain->solInfo().skipOffSet, 0); // skipOffSet must be >= 0
  // these flags are for training with snapshots from multiple simulations in which the mesh is deformed
  const bool poscfg     = domain->solInfo().activatePOSCFG;
  const bool sclfactor  = (domain->solInfo().xScaleFactors.size() > 0);
  const bool ndfile     = (ndscfgCoords && domain->solInfo().NodeTrainingFiles.size() > 0);
  const bool swtchbasis = (ndscfgMassOrthogonalBases && domain->solInfo().MassOrthogonalBasisFiles.size() > 0);
  double xScaleFactor, yScaleFactor, zScaleFactor;
  std::ifstream sources;

  for (int iElem = 0; iElem != elementCount(); ++iElem) { // outer loop over element set
    filePrint(stderr,"\r %4.2f%% complete", double(iElem)/double(elementCount())*100.);
    std::vector<double>::iterator timeStampIt = timeStampFirst;
    int *nodes = domain_->getElementSet()[iElem]->nodes();
    int iSnap = 0;
    for(int i = ((j<0)?0:j); i < ((j<0)?snapshotCounts_.size():(j+1)); i++) { // loop over type of snapshot
      if(!domain_->solInfo().conwepConfigurations.empty()) {
        conwep = &domain_->solInfo().conwepConfigurations[i];
      }
      if(poscfg) { // open sources file 
        std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::SNAPSHOTS, i);
        fileName.append(".sources");
        sources.open(fileName.c_str());
        // this discards everything up to the \n delimiter. called skipOffSet times to account for potential offset
        for(int k = 0; k < skipOffSet; ++k) sources.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }

      for (int jSnap = 0; jSnap != snapshotCounts_[i]; ++iSnap, ++jSnap) { // inner loop over snapshot list
        if(poscfg) { // if there is a transformation supplied, apply to mesh for current snapshot
          sources >> kParam >> kSnap;
          for(int k = 0; k < skipFactor; ++k) sources.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // account for snapshot skipping
          if(kParam != jParam) { // check to see if new transformation is required, if not continue
            if(sclfactor) { // if doing linear stretching of coordinates, do this block
              if(jParam == -1) {
                xScaleFactor = domain->solInfo().xScaleFactors[kParam-1]/domain->solInfo().xScaleFactor;
                yScaleFactor = domain->solInfo().yScaleFactors[kParam-1]/domain->solInfo().yScaleFactor;
                zScaleFactor = domain->solInfo().zScaleFactors[kParam-1]/domain->solInfo().zScaleFactor;
              }
              else {
                xScaleFactor = domain->solInfo().xScaleFactors[kParam-1]/domain->solInfo().xScaleFactors[jParam-1];
                yScaleFactor = domain->solInfo().yScaleFactors[kParam-1]/domain->solInfo().yScaleFactors[jParam-1];
                zScaleFactor = domain->solInfo().zScaleFactors[kParam-1]/domain->solInfo().zScaleFactors[jParam-1];
              }
              //std::cerr << "scaling position coordinates by " << domain->solInfo().xScaleFactor << " " << domain->solInfo().yScaleFactor
              //          << " " << domain->solInfo().zScaleFactor << std::endl;
              geomState_->transformCoords(xScaleFactor, yScaleFactor, zScaleFactor);
              if(domain_->solInfo().getNLInfo().linearelastic) {
                kelArray_[iElem].copy(domain_->getElementSet()[iElem]->stiffness(domain_->getNodes(),  kelArray_[iElem].data()));
                melArray_[iElem].copy(domain_->getElementSet()[iElem]->massMatrix(domain_->getNodes(), melArray_[iElem].data(), geoSource->getMRatio()));
              }
            }
            else if(ndfile) { // if using arbitrarily deformed mesh, do this block
              geomState_->setNewCoords(ndscfgCoords[kParam-1][0]);
            }
            jParam = kParam;
            if(swtchbasis) { // if using mass orthonormalized basis, need to switch to that basis
              int numClust = domain_->solInfo().readInROBorModes.size(); // stride through appropriate local training bases
              podBasis = &ndscfgMassOrthogonalBases[((j<0)?0:j)+(numClust*(kParam-1))];
            }
          }
        }

        geomState_->explicitUpdate(domain_->getNodes(), domain_->getElementSet()[iElem]->numNodes(),
            nodes, displac[iSnap]); // just set the state at the nodes of element iElem
        if(veloc) geomState_->setVelocity(domain_->getElementSet()[iElem]->numNodes(), nodes,
            (*veloc)[iSnap], 2); // just set the velocity at the nodes of element iElem
        if(accel) geomState_->setAcceleration(domain_->getElementSet()[iElem]->numNodes(), nodes,
            (*accel)[iSnap], 2); // just set the acceleration at the nodes of element iElem
        // Evaluate and store element contribution at training configuration
        if(corotators_[iElem] && (!domain_->solInfo().getNLInfo().linearelastic ||
          (domain_->getElementSet()[iElem]->isConstraintElement() && domain_->solInfo().getNLInfo().linearelastic == 2))) {
          domain_->getElemInternalForce(*geomState_, *timeStampIt, geomState_, *(corotators_[iElem]), elementForce.array(), kelArray_[iElem]);
        }
        else {
          Vector disp(kelArray_[iElem].dim());
          StackVector force(elementForce.array(), kelArray_[iElem].dim());
          domain->getElementDisp(iElem, *geomState_, disp);
          force.zero();
          kelArray_[iElem].multiply(disp, force, 1.0);
          if(domain_->solInfo().useConstantMassForces){
            Vector accel(melArray_[iElem].dim()); 
            StackVector *iForce;
            if(domain_->solInfo().stackedElementSampling){
              iForce = new StackVector(elementForce2->array(), melArray_[iElem].dim());; 
              iForce->zero();
            } else {
              iForce = new StackVector(elementForce.array(), melArray_[iElem].dim());
            }
            domain->getElementAccel(iElem, *geomState_, accel);
            melArray_[iElem].multiply(accel, *iForce, -1.0);
          }
        }
        if(domain_->solInfo().reduceFollower)
          domain_->getElemFollowerForce(iElem, *geomState_, elementForce.array(), elementForce.size(), (corotators_[iElem]),
                                        kelArray_[iElem], 1.0, *timeStampIt, false, conwep);
        if(domain_->getElementSet()[iElem]->hasRot() && !domain_->solInfo().getNLInfo().linearelastic) {
          domain_->transformElemStiffAndForce(*geomState_, elementForce.array(), kelArray_[iElem], iElem, false);
          domain_->getElemFictitiousForce(iElem, *geomState_, elementForce.array(), kelArray_[iElem],
              *timeStampIt, geomState_, melArray_[iElem], false);
        }

        elemTarget.zero();
        if(domain_->solInfo().stackedElementSampling)
          elemTarget2->zero();
        const int dofCount = kelArray_[iElem].dim();
        for (int iDof = 0; iDof != dofCount; ++iDof) {
          const int vecLoc = domain_->getCDSA()->getRCN((*domain_->getAllDOFs())[iElem][iDof]);
          if (vecLoc >= 0) {
            const double dofForce = elementForce[iDof];
            double dofForce2 = 0.0; 
            if(domain_->solInfo().stackedElementSampling)
              dofForce2 = (*elementForce2)[iDof];
            for (int iPod = 0; iPod != podVectorCount; ++iPod) {
              const double contrib = dofForce * (*podBasis)[iPod][vecLoc];
              double contrib2 = 0.0;
              if(domain_->solInfo().stackedElementSampling){
                contrib2 = dofForce2 * (*podBasis)[iPod][vecLoc];
                (*elemTarget2)[iPod] += contrib2;
              }
              elemTarget[iPod] += contrib;
            }
          }
        }
        for(int iPod = 0; iPod != podVectorCount; ++iPod) {
          *elemContributions = elemTarget[iPod];
          elemContributions++;
          if(domain_->solInfo().stackedElementSampling){
            *elemContributions = (*elemTarget2)[iPod];
            elemContributions++;
            trainingTarget[2*podVectorCount * iSnap + 2*iPod]   +=   elemTarget[iPod];
            trainingTarget[2*podVectorCount * iSnap + 2*iPod+1] += (*elemTarget2)[iPod];
          } else {
            trainingTarget[podVectorCount * iSnap + iPod] += elemTarget[iPod];
          }
        }
        timeStampIt++;
      }
      // reset the element internal states for the next group of snapshots, assuming for now that each group of snapshots
      // was generated by an independent simulation rather than a restart. For a restart it would be appropriate to not
      // reset the internal states. This could be controlled by a user defined switch.
      int numStates = geomState_->getNumElemStates(domain_->getElementSet()[iElem]->getGlNum());
      if(numStates > 0) {
        double *states = geomState_->getElemState(domain_->getElementSet()[iElem]->getGlNum());
        for (int j = 0; j < numStates; ++j) states[j] = 0;
        domain_->getElementSet()[iElem]->initStates(states);
      }
      if(poscfg) sources.close();
    }

    delete [] nodes;
  }
  filePrint(stderr,"\r %4.2f%% complete\n", 100.);

  if(domain_->solInfo().stackedElementSampling){
    delete elemTarget2;
    delete elementForce2;
  } 

  // finally, undo scaling to original mesh for post processing
  if(poscfg) {
    if(sclfactor) {
      xScaleFactor = domain->solInfo().xScaleFactor/domain->solInfo().xScaleFactors[kParam-1];
      yScaleFactor = domain->solInfo().yScaleFactor/domain->solInfo().yScaleFactors[kParam-1];
      zScaleFactor = domain->solInfo().zScaleFactor/domain->solInfo().zScaleFactors[kParam-1];
      geomState_->transformCoords(xScaleFactor, yScaleFactor, zScaleFactor);
    }
    else if(ndfile) {
      geomState_->setNewCoords(ndscfgCoords[domain->solInfo().NodeTrainingFiles.size()-1][0]);
    }
  }
}


template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::solve() {

  AllOps<double> allOps;
  preProcessGlobal(allOps);

  std::vector<int> sampleElemIds, packedToInput;
  std::vector<std::map<int, double> > weights(domain_->solInfo().readInROBorModes.size());
  makePackedToInput(packedToInput);

  for(int j=0; j<domain_->solInfo().readInROBorModes.size(); ++j) {
   
    preProcessLocal(allOps, j);

    Vector solution;

    if (geoSource->elementLumpingWeightSize() == 0 || geoSource->elementLumpingWeightLocalSize(0) == 0) {
        // Training target (solver_.rhsBuffer) is the sum of elementary contributions
      #ifdef PRINT_ESTIMERS
          double t2 = getTime();
      #endif
          for(int i=0; i<solver_.equationCount(); ++i) solver_.rhsBuffer()[i] = 0.0;
          assembleTrainingData(podBasis_, podBasis_.vectorCount(), displac_, veloc_, accel_, ndscfgCoords_, ndscfgMassOrthogonalBases_,
                               ((domain_->solInfo().readInROBorModes.size() == 1) ? -1 : j));
      #ifdef PRINT_ESTIMERS
          fprintf(stderr, "time for assembleTrainingData = %f\n", (getTime()-t2)/1000.0);
      #endif

          computeSolution(solution, domain_->solInfo().tolPodRom);
    }
    else {
        solution.initialize(elementCount());
        for (GeoSource::ElementWeightMap::const_iterator it = geoSource->elementLumpingWeightBegin(j),
                                                          it_end = geoSource->elementLumpingWeightEnd(j);
                it != it_end; ++it) {
                solution[domain_->glToPackElem(it->first)] = it->second;
        }
    } 
    postProcessLocal(solution, packedToInput, j, sampleElemIds, weights[j]);
  }
  addContactElems(sampleElemIds, weights);
  postProcessGlobal(sampleElemIds, weights);
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::addContactElems(std::vector<int> &sampleElemIds, std::vector<std::map<int, double> > &weights) {

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
    std::unique_ptr<Connectivity> elemToNode(new Connectivity(inputElemSet.asSet()));
    Connectivity nodeToElem = elemToNode->reverse();

    //Second: loop over list of surfaces to extract elements attached to surface nodes 
    for(std::set<int>::iterator it = activeSurfs.begin(); it != activeSurfs.end(); ++it){
      for(int surf = 0; surf < domain->getNumSurfs(); ++surf) { // loop through surfaces and check for matching ID
        if(*it == domain->GetSurfaceEntity(surf)->GetId()) {
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
                for(int lm = 0; lm < weights.size(); ++lm)
                  weights[lm].insert(std::make_pair(elemRank, 0.0)); 
              }
            }
          }
        }
      }
    }
  }
  std::sort(sampleElemIds.begin(), sampleElemIds.end());
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::computeSolution(Vector &solution, double relativeTolerance, bool verboseFlag) {

  solver_.relativeToleranceIs(relativeTolerance);
  solver_.verboseFlagIs(verboseFlag);
  solver_.scalingFlagIs(domain->solInfo().useScalingSpnnls);
  solver_.centerFlagIs(domain->solInfo().useCenterSpnnls);
  solver_.reverseFlagIs(domain->solInfo().useReverseOrder);
  solver_.projectFlagIs(domain->solInfo().projectSolution);
  solver_.positivityIs(domain->solInfo().positiveElements);
  solver_.solverTypeIs(domain->solInfo().solverTypeSpnnls);
  solver_.maxSizeRatioIs(domain->solInfo().maxSizeSpnnls);
  solver_.maxIterRatioIs(domain->solInfo().maxIterSpnnls);
  solver_.maxNumElemsIs(domain->solInfo().maxElemSpnnls);
  solver_.useHotStart(domain->solInfo().hotstartSample);
  try {
    solver_.solve();
  }
  catch(std::runtime_error& e) {
    std::cerr << " *** WARNING: " << e.what() << std::endl;
  }

  if(verboseFlag) {
    std::cout << "Primal solution:";
    for (int elemRank = 0; elemRank != elementCount(); ++elemRank) {
      std::cout << " " << solver_.solutionEntry(elemRank);
    }
    std::cout << "\n";

    StackVector trainingTarget(solver_.rhsBuffer(), solver_.equationCount());
    std::cout << "Error magnitude / Absolute tolerance = " << solver_.errorMagnitude() << " / " << solver_.relativeTolerance() * trainingTarget.norm() << "\n";
    std::cout << "1-norm of primal solution = " << std::accumulate(solver_.solutionBuffer(), solver_.solutionBuffer() + solver_.unknownCount(), 0.0) << "\n";
  }

  // Read solution
  solution.initialize(elementCount());
  std::copy(solver_.solutionBuffer(), solver_.solutionBuffer() + elementCount(), solution.data());
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::postProcess(Vector &solution, bool verboseFlag) {

  std::vector<int> sampleElemIds, packedToInput;
  std::vector<std::map<int, double> > weights(1);

  makePackedToInput(packedToInput);
  postProcessLocal(solution, packedToInput, 0, sampleElemIds, weights[0], verboseFlag);
  postProcessGlobal(sampleElemIds, weights, verboseFlag);
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::makePackedToInput(std::vector<int> &packedToInput)
{
  packedToInput.resize(elementCount());
  Elemset &inputElemSet = *(geoSource->getElemSet());
  for (int iElem = 0, iElemEnd = inputElemSet.size(); iElem != iElemEnd; ++iElem) {
    Element *elem = inputElemSet[iElem];
    if (elem) {
      const int iPackElem = domain_->glToPackElem(iElem);
      if(iPackElem >= 0) {
        assert(iPackElem < packedToInput.size());
        packedToInput[iPackElem] = iElem;
      }
    }
  }
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::postProcessLocal(Vector &solution, std::vector<int> packedToInput, int j,
                                                                   std::vector<int> &sampleElemIds, std::map<int, double> &weights,
                                                                   bool verboseFlag)
{
  std::set<int> sampleElemRanks;
  {
    for (int iElem = 0; iElem != elementCount(); ++iElem) {
      if (solution[iElem] > 0.0) {
        sampleElemRanks.insert(sampleElemRanks.end(), iElem);
      }
    }
  }

  sampleElemIds.reserve(sampleElemRanks.size());
  for (std::set<int>::const_iterator it = sampleElemRanks.begin(), it_end = sampleElemRanks.end(); it != it_end; ++it) {
    const int elemRank = packedToInput[*it];
    weights.insert(std::make_pair(elemRank, solution[*it]));
    if(domain_->solInfo().readInROBorModes.size() == 1 || std::find(sampleElemIds.begin(), sampleElemIds.end(), elemRank) == sampleElemIds.end())
      sampleElemIds.push_back(elemRank);
  }

  outputFullWeights(solution, packedToInput, ((domain_->solInfo().readInROBorModes.size() == 1) ? -1 : j));
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::postProcessGlobal(std::vector<int> &sampleElemIds, std::vector<std::map<int, double> > &weights,
                                                                    bool verboseFlag)
{
  const FileNameInfo fileInfo;

  if(domain_->solInfo().localBasisSize.size() > 1) {
    int globalBasisSize = std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.end(), 0);
    podBasis_.localBasisIs(0, globalBasisSize);
  }

  // compute the reduced forces
  Vector forceFull(SingleDomainDynamic::solVecInfo(),0.0);
  // 1) gravity
  Vector gravForceRed(podBasis_.vectorCount());
  if(domain->gravityFlag()) {
    domain->addGravityForce(forceFull);
    reduce(podBasis_, forceFull, gravForceRed);
  }
  // 2) constant force or constant part of time-dependent forces (default loadset only) TODO add support for multiple loadsets
  std::map<std::pair<int,int>,int>& loadfactor_MFTT =  domain->getLoadFactorMFTT();
  int nLoadSets = std::max((size_t)1,loadfactor_MFTT.size());
  std::vector<int> reduce_f;
  GenVectorSet<double> constForceRed(nLoadSets, podBasis_.vectorCount());
  if(loadfactor_MFTT.size() > 0) {
    for(std::map<std::pair<int,int>,int>::iterator it=loadfactor_MFTT.begin(); it!=loadfactor_MFTT.end(); it++) {
      domain_->solInfo().loadcases.clear();
      domain_->solInfo().loadcases.push_back(it->first.first);
      int j = it->first.second; // loadset_id
      if(std::find(reduce_f.begin(),reduce_f.end(),j) != reduce_f.end()) continue;
      domain->computeUnamplifiedExtForce(forceFull, j);
      if(forceFull.norm() != 0) {
        std::cerr << " ... Computing reduced forces for loadset " << j << " ...\n";
        reduce(podBasis_, forceFull, constForceRed[reduce_f.size()]);
        reduce_f.push_back(j);
      }
    }
  }
  else {
    domain->computeUnamplifiedExtForce(forceFull, 0);
    if(forceFull.norm() != 0) {
      std::cerr << " ... Computing reduced forces ..\n";
      reduce(podBasis_, forceFull, constForceRed[0]);
      reduce_f.push_back(0);
    }
  }

  // 3) set LMPC
  int numLMPC = domain->getNumLMPC(); 

  // compute the reduced initial conditions
  Vector d0Full(SingleDomainDynamic::solVecInfo()),
         v0Full(SingleDomainDynamic::solVecInfo());
  Vector tmp(SingleDomainDynamic::solVecInfo());
  SysState<Vector> inState(d0Full, v0Full, tmp, tmp);
  SingleDomainDynamic::getInitState(inState);
  Vector d0Red(podBasis_.vectorCount()),
         v0Red(podBasis_.vectorCount());
  bool reduce_idis = (d0Full.norm() != 0),
       reduce_ivel = (v0Full.norm() != 0);
#ifdef USE_EIGEN3
  GenEiSparseMatrix<double,Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> > *M;
#endif
  AllOps<double> allOps;
  if(((domain->solInfo().useMassNormalizedBasis || domain->solInfo().newmarkBeta == 0)
       && (reduce_idis || reduce_ivel || domain_->solInfo().readInROBorModes.size() > 1)) || !domain->solInfo().useMassNormalizedBasis) {
#ifdef USE_EIGEN3
    allOps.M = M = domain_->constructEiSparseMatrix<double,Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> >();
#else
    allOps.M = domain->constructDBSparseMatrix<double>();
#endif
    domain->makeSparseOps(allOps, 0, 0, 0);
  }
  if(domain->solInfo().useMassNormalizedBasis || domain->solInfo().newmarkBeta == 0) {
    if(reduce_idis) {
      allOps.M->mult(d0Full, tmp);
      reduce(podBasis_, tmp, d0Red);
    }
    if(reduce_ivel) {
      allOps.M->mult(v0Full, tmp);
      reduce(podBasis_, tmp, v0Red);
    }
  }
  else {
    if(reduce_idis) reduce(podBasis_, d0Full, d0Red);
    if(reduce_ivel) reduce(podBasis_, v0Full, v0Red);
  }

  // output the reduced mesh
  Elemset &inputElemSet = *(geoSource->getElemSet());
  std::unique_ptr<Connectivity> elemToNode(new Connectivity(inputElemSet.asSet()));
  const MeshRenumbering meshRenumbering(sampleElemIds.begin(), sampleElemIds.end(), *elemToNode, verboseFlag);
  const MeshDesc reducedMesh(domain_, geoSource, meshRenumbering, weights);
  if(domain_->solInfo().localBasisSize.size() <= 1)
    outputMeshFile(fileInfo, reducedMesh, podBasis_.vectorCount());
  else
    outputMeshFile(fileInfo, reducedMesh, domain_->solInfo().localBasisSize);

  // output the reduced forces
  std::ofstream meshOut(getMeshFilename(fileInfo).c_str(), std::ios_base::app);
  if(domain->solInfo().reduceFollower) meshOut << "*\nEXTFOL\n";
  if(domain->gravityFlag()) {
    meshOut << "*\nFORCES -1\nMODAL\n"; // note: gravity forces are put in loadset -1 so that MFTT (if present) will not be applied
    meshOut.precision(std::numeric_limits<double>::digits10+1);
    for(int i=0; i<podBasis_.vectorCount(); ++i)
      meshOut << i+1 << " " << gravForceRed[i] << std::endl;
  }

  for(std::vector<int>::iterator it = reduce_f.begin(); it != reduce_f.end(); ++it) {
    if(*it == 0) meshOut << "*\nFORCES\nMODAL\n";
    else meshOut << "*\nFORCES " << *it << "\nMODAL\n";
    meshOut.precision(std::numeric_limits<double>::digits10+1);
    int i = std::distance(reduce_f.begin(), it);
    for (int l=0; l<podBasis_.vectorCount(); ++l) 
      meshOut << l+1 << " " << constForceRed[i][l] << std::endl;
  }

#ifdef USE_EIGEN3
  if(numLMPC > 0) {

    if(!domain_->solInfo().readInDualROB.empty()) {
 
      // get dual basis sizes
      std::vector<int> locDualBasisVec;
      VecNodeDof1Conversion vecNodeDof1Conversion(numLMPC);
      for(int j = 0; j < domain->solInfo().readInDualROB.size(); ++j){
        std::string fileName = BasisFileId(fileInfo, BasisId::DUALSTATE, BasisId::POD,j);
        BasisInputStream<1> dualProjectionBasisInput(fileName, vecNodeDof1Conversion);

        locDualBasisVec.push_back(dualProjectionBasisInput.size());

        filePrint(stderr, " ... Dual Basis %d size %d ...\n", j, locDualBasisVec[j]);
      }

      meshOut << "*";
      VecBasis dualProjectionBasis_;
      // Load local dual projection basis   
      for(int j = 0 ; j < domain->solInfo().readInDualROB.size(); j++) { 
        std::string fileName = BasisFileId(fileInfo, BasisId::DUALSTATE, BasisId::POD,j);
        BasisInputStream<1> dualProjectionBasisInput(fileName, vecNodeDof1Conversion);
        const int dualProjectionSubspaceSize = locDualBasisVec[j] ?
                                               std::min(locDualBasisVec[j], dualProjectionBasisInput.size()) :
                                               dualProjectionBasisInput.size();
 
        meshOut << "\nLMPC\nMODAL " << dualProjectionSubspaceSize;

        readVectors(dualProjectionBasisInput, dualProjectionBasis_,
                  std::accumulate(locDualBasisVec.begin(), locDualBasisVec.end(), 0),
                  dualProjectionSubspaceSize,
                  std::accumulate(locDualBasisVec.begin(), locDualBasisVec.begin()+j, 0));
      } 

      meshOut << std::endl;

      std::vector<Eigen::Triplet<double> > tripletList;
      Eigen::SparseMatrix<double> C(numLMPC, domain_->getCDSA()->size());
      Eigen::Matrix<double,Eigen::Dynamic,1> g(numLMPC);
      LMPCons** lmpc(domain_->getLMPC()->data());
      
      //construct constraint matrix and right hand side
      for(int i=0; i<numLMPC; ++i) {
        for(int j=0; j<lmpc[i]->nterms; ++j) {
          int cdof = domain_->getCDSA()->locate(lmpc[i]->terms[j].nnum, 1 << lmpc[i]->terms[j].dofnum);
          if(cdof > -1) {
            tripletList.push_back(Eigen::Triplet<double>(i, cdof, double(lmpc[i]->terms[j].coef.r_value)));
          }
        }
        g[i] = lmpc[i]->rhs.r_value;
      } 
 
      C.setFromTriplets(tripletList.begin(), tripletList.end());
      const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &V = podBasis_.basis(), &W = dualProjectionBasis_.basis();
      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> reducedConstraintMatrix_ = W.transpose()*C*V;
      Eigen::Matrix<double,Eigen::Dynamic,1> reducedConstraintRhs0_ = W.transpose()*g;
      // output reduced Constraints to reduced mesh file
      for(int col=0; col<podBasis_.vectorCount(); ++col) {
         meshOut << reducedConstraintMatrix_.col(col) << std::endl;
      }
      meshOut << reducedConstraintRhs0_ << std::endl;

    } else {
      filePrint(stderr, "... No Filename given for Precomputation of reduced Constraints ...\n");
    }

  } else {

    if(!domain->solInfo().readInDualROB.empty()) {
      // first get size of dual basis files
      std::vector<int> locDualBasisVec;
      VecNodeDof1Conversion vecNodeDof1Conversion(*domain_->getCDSA());
      for(int j = 0; j < domain->solInfo().readInDualROB.size(); ++j){
        std::string fileName = BasisFileId(fileInfo, BasisId::DUALSTATE, BasisId::POD,j);
        BasisInputStream<1> dualProjectionBasisInput(fileName, vecNodeDof1Conversion);

        locDualBasisVec.push_back(dualProjectionBasisInput.size());

        filePrint(stderr, " ... Dual Basis %d size %d ...\n", j, locDualBasisVec[j]);
      }

      // next read them in
      VecBasis dualProjectionBasis_;
      for(int j = 0 ; j < domain->solInfo().readInDualROB.size(); j++) {
        std::string fileName = BasisFileId(fileInfo, BasisId::DUALSTATE, BasisId::POD,j);
        BasisInputStream<1> dualProjectionBasisInput(fileName, vecNodeDof1Conversion);
        const int dualProjectionSubspaceSize = locDualBasisVec[j] ?
                                               std::min(locDualBasisVec[j], dualProjectionBasisInput.size()) :
                                               dualProjectionBasisInput.size();

        readVectors(dualProjectionBasisInput, dualProjectionBasis_,
                  std::accumulate(locDualBasisVec.begin(), locDualBasisVec.end(), 0),
                  dualProjectionSubspaceSize,
                  std::accumulate(locDualBasisVec.begin(), locDualBasisVec.begin()+j, 0));
      }

      // then set up new compressed file to write them to
      DofSetArray reduced_dsa(reducedMesh.nodes().size(), const_cast<Elemset&>(reducedMesh.elements()));
      int num_bc = reducedMesh.dirichletBConds().size();
      BCond *bc = (num_bc > 0) ? const_cast<BCond*>(&reducedMesh.dirichletBConds()[0]) : NULL;
      ConstrainedDSA reduced_cdsa(reduced_dsa, num_bc, bc);
      dualProjectionBasis_.makeSparseBasis(meshRenumbering.reducedNodeIds(), domain_->getCDSA(), &reduced_cdsa);
      {
        VecNodeDof1Conversion converter(reduced_cdsa);
        for(int j=0; j<domain_->solInfo().readInDualROB.size(); ++j) {
          std::string filename = getMeshFilename(fileInfo).c_str();
          filename.append(".compressed.dualbasis");
          if(domain_->solInfo().readInDualROB.size() != 1) {
            std::ostringstream ss;
            ss << j+1;
            filename.append(ss.str());
          }
          filePrint(stderr," ... Writing dual compressed basis to file %s ...\n", filename.c_str());
          BasisOutputStream<1> output(filename, converter, false);
 
          int startCol = (domain_->solInfo().readInDualROB.size() == 1) ? 0 :
                        std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.begin()+j, 0);
          int blockCols = (domain_->solInfo().readInDualROB.size() == 1) ? dualProjectionBasis_.vectorCount() : domain_->solInfo().localBasisSize[j];
          for (int iVec = startCol; iVec < startCol+blockCols; ++iVec) {
            output << dualProjectionBasis_.compressedBasis().col(iVec);
          }
        }
      }

    }
  }
#endif

  // output the reduced initial conditions
  if(reduce_idis) {
    meshOut << "*\nIDISPLACEMENTS\nMODAL\n";
    meshOut.precision(std::numeric_limits<double>::digits10+1);
    for(int i=0; i<podBasis_.vectorCount(); ++i)
      meshOut << i+1 << " " << d0Red[i] << std::endl;
  }
  if(reduce_ivel) {
    meshOut << "*\nIVELOCITIES\nMODAL\n";
    meshOut.precision(std::numeric_limits<double>::digits10+1);
    for(int i=0; i<podBasis_.vectorCount(); ++i)
      meshOut << i+1 << " " << v0Red[i] << std::endl;
  }

#ifdef USE_EIGEN3
  // pre-computation required for local bases method
  if(domain_->solInfo().readInROBorModes.size() > 1 && 
     (domain->solInfo().readInLocalBasesCent.size() == domain_->solInfo().readInROBorModes.size())) {
    // read cluster centroids
    const int n = podBasis_.vectorSize();
    Eigen::MatrixXd uc(n, domain->solInfo().readInLocalBasesCent.size());
    const VecNodeDof6Conversion vecNodeDof6Conversion(*domain_->getCDSA());
    if(domain->solInfo().readInLocalBasesCent.size() == domain_->solInfo().readInROBorModes.size()) {
      for(int i=0; i<domain->solInfo().readInLocalBasesCent.size(); ++i) {
        BasisInputStream<6> in(domain->solInfo().readInLocalBasesCent[i], vecNodeDof6Conversion);
        in >> uc.col(i).data();
      }
    }
 
    // compute and output auxiliary quantities
    meshOut << "*\nLOCROB\n";
    for(int i=0; i<domain_->solInfo().readInROBorModes.size(); ++i) {
      int startColi = std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.begin()+i, 0);
      int blockColsi = domain_->solInfo().localBasisSize[i];
      for(int j=i+1; j<domain_->solInfo().readInROBorModes.size(); ++j) {
        int startColj = std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.begin()+j, 0);
        int blockColsj = domain_->solInfo().localBasisSize[j];
        Eigen::MatrixXd VtVij;
        std::string fileName = domain->solInfo().reducedMeshFile;
        std::ostringstream ss;
        ss << ".auxiliary.cluster" << i+1 << "." << j+1;
        fileName.append(ss.str());
        meshOut << "auxi " << i+1 << " " << j+1 << " \"" << fileName << "\"\n";
        filePrint(stderr," ... Writing local bases auxiliary quantities to file %s ...\n", fileName.c_str());
        std::ofstream matrixOut(fileName.c_str());
        // matrix: Vi.transpose()*Vj (or Vi.transpose()*M*Vj)
        if(domain_->solInfo().newmarkBeta == 0 || domain_->solInfo().useMassNormalizedBasis) {
          VtVij = podBasis_.basis().block(0,startColi,n,blockColsi).transpose()*
                  (M->getEigenSparse().selfadjointView<Eigen::Upper>()*podBasis_.basis().block(0,startColj,n,blockColsj)); 
        }
        else {
           VtVij = podBasis_.basis().block(0,startColi,n,blockColsi).transpose()*
                   podBasis_.basis().block(0,startColj,n,blockColsj);
        }
        matrixOut << std::setprecision(16) << VtVij << std::endl;
        // scalar: ||uci|| - ||ucj||
        double d = uc.col(i).squaredNorm() - uc.col(j).squaredNorm();
        matrixOut << d << std::endl;
        // vector: 2*V.transpose()*(ucj-uci)
        Eigen::VectorXd w = 2*podBasis_.basis().transpose()*(uc.col(j)-uc.col(i));
        matrixOut << w.transpose() << std::endl;

        matrixOut.close();
      }
    }
  }

  // build and output reduced mass matrix if mass-normalized basis is not used
  if(!(domain_->solInfo().newmarkBeta == 0 || domain_->solInfo().useMassNormalizedBasis)) {
    std::string fileName = domain->solInfo().reducedMeshFile;
    fileName.append(".reducedmass");
    meshOut << "*\nDIMASS\nMODAL\n\"" << fileName << "\"\n";
    std::ofstream matrixOut(fileName.c_str());
    filePrint(stderr," ... Writing reduced mass matrix to file %s ...\n", fileName.c_str());
    matrixOut << std::setprecision(16) 
              << podBasis_.basis().transpose()*(M->getEigenSparse().selfadjointView<Eigen::Upper>()*podBasis_.basis()) << std::endl;
    matrixOut.close();
  }

  // build and output compressed basis
  DofSetArray reduced_dsa(reducedMesh.nodes().size(), const_cast<Elemset&>(reducedMesh.elements()));
  int num_bc = reducedMesh.dirichletBConds().size();
  BCond *bc = (num_bc > 0) ? const_cast<BCond*>(&reducedMesh.dirichletBConds()[0]) : NULL;
  ConstrainedDSA reduced_cdsa(reduced_dsa, num_bc, bc);
  podBasis_.makeSparseBasis(meshRenumbering.reducedNodeIds(), domain_->getCDSA(), &reduced_cdsa);
  {
    VecNodeDof6Conversion converter(reduced_cdsa);
    for(int j=0; j<domain_->solInfo().readInROBorModes.size(); ++j) {
      std::string filename = getMeshFilename(fileInfo).c_str();
      filename.append(".compressed.basis");
      if(domain_->solInfo().readInROBorModes.size() != 1) {
        std::ostringstream ss;
        ss << j+1;
        filename.append(ss.str());
      }
      if(domain_->solInfo().newmarkBeta == 0 || domain_->solInfo().useMassNormalizedBasis) filename.append(".massorthonormalized");
      else filename.append(".orthonormalized");
      filePrint(stderr," ... Writing compressed basis to file %s ...\n", filename.c_str());
      BasisOutputStream<6> output(filename, converter, false);

      int startCol = (domain_->solInfo().readInROBorModes.size() == 1) ? 0 :
                      std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.begin()+j, 0);
      int blockCols = (domain_->solInfo().readInROBorModes.size() == 1) ? podBasis_.vectorCount() : domain_->solInfo().localBasisSize[j];
      for (int iVec = startCol; iVec < startCol+blockCols; ++iVec) {
        output << podBasis_.compressedBasis().col(iVec);
      }
    }
  }
#endif
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::preProcess()
{
  AllOps<double> allOps;
  preProcessGlobal(allOps);
  preProcessLocal(allOps, 0);
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::preProcessGlobal(AllOps<double>& allOps)
{
  domain_->preProcessing();
  buildDomainCdsa();
  domain_->makeAllDOFs();
  
  StaticTimers dummyTimes;
  GenFullSquareMatrix<double> *dummyGeomKelArray = NULL;
  const bool buildMelArray = true;
  domain_->computeGeometricPreStress(corotators_, geomState_, kelArray_, &dummyTimes, dummyGeomKelArray, melArray_, buildMelArray);
  if(domain_->nDirichlet() > 0) {
    geomState_->updatePrescribedDisplacement(domain_->getDBC(), domain_->nDirichlet(), domain_->getNodes());
  }
  if(domain_->solInfo().newmarkBeta == 0) {
    domain_->assembleNodalInertiaTensors(melArray_);
  }

  // Assemble mass matrix if necessary
  if(domain_->solInfo().useMassOrthogonalProjection) {
    if(geoSource->getMRatio() != 0) {
      allOps.M = domain_->constructDBSparseMatrix<double>();
    }
    else {
      allOps.M = new DiagMatrix(domain->getCDSA());
    }
    domain_->makeSparseOps<double>(allOps, 0.0, 1.0, 0.0);
  }
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::preProcessLocal(AllOps<double> &allOps, int j)
{
  const FileNameInfo fileInfo;
  
  // Read order reduction data
  const VecNodeDof6Conversion vecDofConversion(*domain_->getCDSA());
  assert(vectorSize() == vecDofConversion.vectorSize()); 

  // Read in basis to be used for the projection:
  // (a) if a mass-orthogonal projection is to be done, then the mass-normalized basis will be read
  // (b) if and orthogonal projection is to be done, then the identity-normalized basis will be read
  {
    std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD, j);
    if(domain_->solInfo().useMassOrthogonalProjection) fileName.append(".massorthonormalized");
    else fileName.append(".orthonormalized");
    BasisInputStream<6> in(fileName, vecDofConversion);
    if(domain_->solInfo().readInROBorModes.size() == 1) {
      filePrint(stdout," ... Reading basis from %s ...\n", fileName.c_str());
      const int podSizeMax = domain_->solInfo().maxSizePodRom;
      if(podSizeMax != 0) {
        readVectors(in, podBasis_, podSizeMax);
      } else {
        readVectors(in, podBasis_);
      }
    } else {
      int globalBasisSize = std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.end(), 0);
      if(globalBasisSize <= 0) {
        int sillyCounter = 0;
        for(std::vector<int>::iterator it = domain_->solInfo().localBasisSize.begin(); it != domain_->solInfo().localBasisSize.end(); it++) {
          std::string dummyName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD, sillyCounter); sillyCounter++;
          if(domain_->solInfo().useMassOrthogonalProjection) dummyName.append(".massorthonormalized");
          else dummyName.append(".orthonormalized");
          BasisInputStream<6> dummyIn(dummyName, vecDofConversion);
          *it = dummyIn.size();
          globalBasisSize += *it;
        }
      }
      int startCol = std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.begin()+j, 0);
      int blockCols = domain_->solInfo().localBasisSize[j];
      readVectors(in, podBasis_, globalBasisSize, blockCols, startCol);
      podBasis_.localBasisIs(startCol, blockCols);
    }
  }

  const int podVectorCount = podBasis_.vectorCount();

  // Optionally, read some displacement snapshots from one or more files and project them on to the basis
  if (!domain_->solInfo().statePodRomFile.empty())
    readAndProjectSnapshots(BasisId::STATE, vectorSize(), podBasis_, vecDofConversion,
                              snapshotCounts_, timeStamps_, displac_, allOps.M,
                              ((domain_->solInfo().readInROBorModes.size() == 1) ? -1 : j));
  else
    snapshotCounts_.push_back(0);

  // Optionally, read some velocity snapshots and project them on to the reduced order basis
  if(!domain_->solInfo().velocPodRomFile.empty()) {
    std::vector<double> velTimeStamps;
    std::vector<int> velSnapshotCounts;
    veloc_ = new VecBasis;
    readAndProjectSnapshots(BasisId::VELOCITY, vectorSize(), podBasis_, vecDofConversion,
                            velSnapshotCounts, velTimeStamps, *veloc_, allOps.M,
                            ((domain_->solInfo().readInROBorModes.size() == 1) ? -1 : j));
    if(velSnapshotCounts != snapshotCounts_) std::cerr << " *** WARNING: inconsistent velocity snapshots\n";
  }

  // Optionally, read some acceleration snapshots and project them on to the reduced order basis
  if(!domain_->solInfo().accelPodRomFile.empty()) {
    std::vector<double> accTimeStamps;
    std::vector<int> accSnapshotCounts;
    accel_ = new VecBasis;
    readAndProjectSnapshots(BasisId::ACCELERATION, vectorSize(), podBasis_, vecDofConversion,
                            accSnapshotCounts, accTimeStamps, *accel_, allOps.M,
                            ((domain_->solInfo().readInROBorModes.size() == 1) ? -1 : j));
    if(accSnapshotCounts != snapshotCounts_) std::cerr << " *** WARNING: inconsistent acceleration snapshots\n";
  }
  
  // Read in training basis if necessary (i.e. if it is different to the projection basis)
  bool useMassNormalizedBasis = (domain_->solInfo().newmarkBeta == 0 || domain_->solInfo().useMassNormalizedBasis);
  if((useMassNormalizedBasis && !domain_->solInfo().useMassOrthogonalProjection) ||
     (!useMassNormalizedBasis && domain_->solInfo().useMassOrthogonalProjection))  {
    std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD, j);
    if(useMassNormalizedBasis) fileName.append(".massorthonormalized");
    else fileName.append(".orthonormalized");
    if(domain_->solInfo().readInROBorModes.size() == 1)
      fprintf(stdout," ... reading %d basis vectors from %s ...\n",domain_->solInfo().maxSizePodRom,fileName.c_str());
    else
      fprintf(stdout," ... reading %d bases vectors from %s ...\n",domain_->solInfo().localBasisSize[j],fileName.c_str());
    BasisInputStream<6> in(fileName, vecDofConversion);
    if(domain_->solInfo().readInROBorModes.size() == 1) {
      const int podSizeMax = domain_->solInfo().maxSizePodRom;
      if(podSizeMax != 0) {
        readVectors(in, podBasis_, podSizeMax);
      } else {
        readVectors(in, podBasis_);
      }
    }
    else {
      int globalBasisSize = std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.end(), 0);
      int startCol = std::accumulate(domain_->solInfo().localBasisSize.begin(), domain_->solInfo().localBasisSize.begin()+j, 0);
      int blockCols = domain_->solInfo().localBasisSize[j];
      readVectors(in, podBasis_, globalBasisSize, blockCols, startCol);
    }
  }

  // read in ndscfg coordinates
  if(domain->solInfo().NodeTrainingFiles.size() > 0) {
    const VecNodeDofConversion<6> vecDofConversion(domain_->getCDSA()->numNodes()); // note: no constraints
    ndscfgCoords_ = new VecBasis[domain->solInfo().NodeTrainingFiles.size()];
    for(int i = 0; i < domain->solInfo().NodeTrainingFiles.size(); ++i) {
      BasisInputStream<6> in(domain->solInfo().NodeTrainingFiles[i], vecDofConversion); // XXX only need 3 dofs per node, not 6.
      readVectors(in, ndscfgCoords_[i], 1);
    }
  }

  // read in ndscfg mass-normalized bases
  if(useMassNormalizedBasis && domain->solInfo().MassOrthogonalBasisFiles.size() > 0) {
    ndscfgMassOrthogonalBases_ = new VecBasis[domain->solInfo().MassOrthogonalBasisFiles.size()];
    for(int i = 0; i < domain->solInfo().MassOrthogonalBasisFiles.size(); ++i) {
      int numClust = domain_->solInfo().readInROBorModes.size(); // stride through appropriate local training bases
      BasisInputStream<6> in(domain->solInfo().MassOrthogonalBasisFiles[i], vecDofConversion);
      readVectors(in, ndscfgMassOrthogonalBases_[i], podVectorCount); // make sure to read in number of vectors already allocated for */
    }
  }

  const int snapshotCount = (domain_->solInfo().readInROBorModes.size() == 1)
                          ? std::accumulate(snapshotCounts_.begin(), snapshotCounts_.end(), 0) : snapshotCounts_[j];
  // check whether to separate the inertial and internal forces or not
  if(domain_->solInfo().stackedElementSampling)
    solver_.problemSizeIs(2*podVectorCount*snapshotCount, elementCount());
  else 
    solver_.problemSizeIs(podVectorCount*snapshotCount, elementCount());
}

template<typename MatrixBufferType, typename SizeType>
void
ElementSamplingDriver<MatrixBufferType,SizeType>::buildDomainCdsa() {
  const int numdof = domain_->numdof();
  SimpleBuffer<int> bc(numdof);
  SimpleBuffer<double> bcx(numdof);

  domain_->make_bc(bc.array(), bcx.array());
  domain_->make_constrainedDSA(bc.array());
}

} // end namespace Rom

Rom::DriverInterface *elementSamplingDriverNew(Domain *d) {
  return new Rom::ElementSamplingDriver<std::vector<double>,size_t>(d);
}
