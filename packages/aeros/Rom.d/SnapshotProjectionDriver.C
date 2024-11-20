#include "SnapshotProjectionDriver.h"

#include "VecBasis.h"
#include "VecBasisOps.h"
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
#include <Math.d/Vector.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/DiagMatrix.h>
#include <Timers.d/StaticTimers.h>
#include <Utils.d/Connectivity.h>
#include <Element.d/Element.h>
#include <Driver.d/SysState.h>

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

// Forward declarations
// ====================
void readAndProjectSnapshots(BasisId::Type type, const int vectorSize, VecBasis &podBasis,
                             const VecNodeDof6Conversion &vecDofConversion,
                             std::vector<int> &snapshotCounts, std::vector<double> &timeStamps, VecBasis &config, SparseMatrix *M, int j=-1);

// Member functions
// ================
int
SnapshotProjectionDriver::elementCount() const {
  return domain->numElements();
}

int
SnapshotProjectionDriver::vectorSize() const {
  return domain->numUncon();
}

SnapshotProjectionDriver::SnapshotProjectionDriver(Domain *d) :
  SingleDomainDynamic(d),
  velocSnapshots(NULL),
  veloc_(NULL),
  accelSnapshots(NULL),
  accel_(NULL)
{}

SnapshotProjectionDriver::~SnapshotProjectionDriver() {
  if(veloc_) delete veloc_;
  if(accel_) delete accel_;
}

void
SnapshotProjectionDriver::postProcess() {

  const int snapshotCount = displac_.vectorCount();

  SDDynamPostProcessor *postProcessor = getPostProcessor();
  DynamMat *dMat = buildOps(1.0,0.0,0.0);
  Vector zero(solVecInfo(), 0.0);
  Vector *v = &zero, *a = &zero, *vp = &zero, *externalForce = &zero, *aeroForce = NULL;

  std::vector<double>::iterator timeStampIt = timeStamps_.begin();
  for (int iSnap = 0; iSnap != snapshotCount; ++iSnap) {
    filePrint(stderr,"\r %4.2f%% complete", double(iSnap)/double(snapshotCount)*100.);
    geomState->explicitUpdate(domain->getNodes(), displac_[iSnap]);
    if(veloc_) { 
      geomState->setVelocity((*veloc_)[iSnap], 2); 
      v = &(*veloc_)[iSnap];
    }
    if(accel_) {
      geomState->setAcceleration((*accel_)[iSnap], 2);
      a = &(*accel_)[iSnap]; 
    }
    SysState<Vector> systemState(displac_[iSnap], *v, *a, *vp);
    postProcessor->dynamOutput(iSnap, *timeStampIt, *dMat, // XXX note: iSnap is not always the correct timeStepIndex
                               *externalForce, aeroForce, systemState);
    timeStampIt++;
  }
 
  filePrint(stderr,"\r %4.2f%% complete\n", 100.);
}

void
SnapshotProjectionDriver::solve() {
  preProcess();
  postProcess();
  compProjError();
}

void
SnapshotProjectionDriver::preProcess() {

  SingleDomainDynamic::preProcess();

  const FileNameInfo fileInfo;
  
  // Read order reduction data
  const VecNodeDof6Conversion vecDofConversion(*domain->getCDSA());
  assert(vectorSize() == vecDofConversion.vectorSize());

  // Read in basis to be used for the projection:
  // (a) if a mass-orthogonal projection is to be done, then the mass-normalized basis will be read
  // (b) if and orthogonal projection is to be done, then the identity-normalized basis will be read
  {
    std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
    if(domain->solInfo().useMassOrthogonalProjection) {
      if(!domain->solInfo().useMassNormalizedBasis) {
        std::string::size_type n = fileName.rfind(".orthonormalized");
        if(n != std::string::npos) 
          fileName = fileName.substr(0,n);
      }
      fileName.append(".massorthonormalized");
    }
    BasisInputStream<6> in(fileName, vecDofConversion);
    const int podSizeMax = domain->solInfo().maxSizePodRom;
    if (podSizeMax != 0) {
      readVectors(in, podBasis_, podSizeMax);
    } else {
      readVectors(in, podBasis_);
    }
  }

  // Assemble mass matrix if necessary
  AllOps<double> allOps;
  if(domain->solInfo().useMassOrthogonalProjection) {
    if(geoSource->getMRatio() != 0) {
      allOps.M = domain->constructDBSparseMatrix<double>();
    }
    else {
      allOps.M = new DiagMatrix(domain->getCDSA());
    }
    domain->makeSparseOps<double>(allOps, 0.0, 1.0, 0.0);
  }

  const int podVectorCount = podBasis_.vectorCount();

  // Read some displacement snapshots from one or more files and project them on to the basis
  std::vector<int> snapshotCounts;
  readAndProjectSnapshots(BasisId::STATE, vectorSize(), podBasis_, vecDofConversion,
                          snapshotCounts, timeStamps_, displac_, allOps.M);

  const int snapshotCount = std::accumulate(snapshotCounts.begin(), snapshotCounts.end(), 0);

  // Optionally, read some velocity snapshots and project them on to the reduced order basis
  if(!domain->solInfo().velocPodRomFile.empty()) {
    std::vector<double> velTimeStamps;
    std::vector<int> velSnapshotCounts;
    veloc_ = new VecBasis;
    readAndProjectSnapshots(BasisId::VELOCITY, vectorSize(), podBasis_, vecDofConversion,
                            velSnapshotCounts, velTimeStamps, *veloc_, allOps.M);
    if(velSnapshotCounts != snapshotCounts) std::cerr << " *** WARNING: inconsistent velocity snapshots\n";
  }

  // Optionally, read some acceleration snapshots and project them on to the reduced order basis
  if(!domain->solInfo().accelPodRomFile.empty()) {
    std::vector<double> accTimeStamps;
    std::vector<int> accSnapshotCounts;
    accel_ = new VecBasis;
    readAndProjectSnapshots(BasisId::ACCELERATION, vectorSize(), podBasis_, vecDofConversion,
                            accSnapshotCounts, accTimeStamps, *accel_, allOps.M);
    if(accSnapshotCounts != snapshotCounts) std::cerr << " *** WARNING: inconsistent acceleration snapshots\n";
  }
}

void
SnapshotProjectionDriver::compProjError() {
#ifdef USE_EIGEN3
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> dispError(snapshots.vectorSize(),snapshots.numVec());
  Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > dispBuf(snapshots.data(),snapshots.vectorSize(),snapshots.numVec());
  Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > pdispBuf(displac_.data(),displac_.vectorSize(),displac_.numVec());

  dispError = pdispBuf - dispBuf;

  if(domain->solInfo().PODerrornorm.size() == 0) {
     std::cerr << "...No displacement file specified, exiting..." << std::endl;
     exit(-1);
    }

  FILE * dispFile = fopen(domain->solInfo().PODerrornorm[0].c_str(),"w");
  filePrint(dispFile,"Number of Training Configurations: %d\n",snapshots.numVec());

  int numNorms = domain->solInfo().PODerrornorm.size();

  filePrint(dispFile,"Displacement Projection Error:\n");
  filePrint(dispFile,"Frobenius: %1.6e\n", dispError.norm()/dispBuf.norm());
  filePrint(dispFile,"   Snapshot   |      L_inf      |      L1       |        L2   \n");
  for (int i = 0; i != snapshots.numVec(); ++i){
    filePrint(dispFile,"     %d      ",i+1);
    for(int pnorm = 0; pnorm != 3; pnorm ++){
      if((pnorm) == 0){
          double Lperror = dispError.col(i).lpNorm<Eigen::Infinity>();
          Lperror = Lperror/dispBuf.col(i).lpNorm<Eigen::Infinity>()*100.;
          filePrint(dispFile,"     %1.6e ",Lperror);
      } else if(pnorm == 1){
          double Lperror = dispError.col(i).lpNorm<1>();
          Lperror = Lperror/dispBuf.col(i).lpNorm<1>()*100.;
          filePrint(dispFile,"    %1.6e ",Lperror);
      } else if(pnorm == 2) {
          double Lperror = dispError.col(i).lpNorm<2>();
          Lperror = Lperror/dispBuf.col(i).lpNorm<2>()*100.;
          filePrint(dispFile,"    %1.6e   ",(pnorm),Lperror);
      }
    }
    filePrint(dispFile,"\n");
  }

  if(velocSnapshots) {
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> velError(velocSnapshots->vectorSize(),velocSnapshots->numVec());
    Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > velBuf(velocSnapshots->data(),velocSnapshots->vectorSize(),velocSnapshots->numVec());
    Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > pvelBuf(veloc_->data(),veloc_->vectorSize(),veloc_->numVec());

    velError = pvelBuf - velBuf;

    if(domain->solInfo().PODerrornorm.size() < 2) {
     std::cerr << "...No velocity file specified, exiting..." << std::endl;
     exit(-1);
    }

    FILE * velFile = fopen(domain->solInfo().PODerrornorm[1].c_str(),"w");
    filePrint(velFile,"Number of Training Configurations: %d\n",velocSnapshots->numVec());

    filePrint(velFile,"Velocity Projection Error:\n");
    filePrint(velFile,"Frobenius: %1.6e\n", velError.norm()/velBuf.norm());
    filePrint(velFile,"   Snapshot   |      L_inf      |      L1       |        L2   \n");
    for (int i = 0; i != velocSnapshots->numVec(); ++i){
      filePrint(velFile,"     %d      ",i+1);
      for(int pnorm = 0; pnorm != 3; pnorm ++){
        if((pnorm) == 0){
            double Lperror = velError.col(i).lpNorm<Eigen::Infinity>();
            Lperror = Lperror/velBuf.col(i).lpNorm<Eigen::Infinity>()*100.;
            filePrint(velFile,"     %1.6e ",Lperror);
        } else if(pnorm == 1) {
            double Lperror = velError.col(i).lpNorm<1>();
            Lperror = Lperror/velBuf.col(i).lpNorm<1>()*100.;
            filePrint(velFile,"    %1.6e ",Lperror);
        } else if(pnorm == 2) {
            double Lperror = velError.col(i).lpNorm<2>();
            Lperror = Lperror/velBuf.col(i).lpNorm<2>()*100.;
            filePrint(velFile,"    %1.6e   ",(pnorm),Lperror);
        }
      }
      filePrint(velFile,"\n");
    }
  }

  if(accelSnapshots) {
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> accelError(accelSnapshots->vectorSize(),accelSnapshots->numVec());
    Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > accelBuf(accelSnapshots->data(),accelSnapshots->vectorSize(),accelSnapshots->numVec());
    Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > paccelBuf(accel_->data(),accel_->vectorSize(),accel_->numVec());

    accelError = paccelBuf - accelBuf;

    if(domain->solInfo().PODerrornorm.size() < 3) {
     std::cerr << "...No acceleration file specified, exiting..." << std::endl;
     exit(-1);
    }

    FILE * accelFile = fopen(domain->solInfo().PODerrornorm[2].c_str(),"w");
    filePrint(accelFile,"Number of Training Configurations: %d\n",accelSnapshots->numVec());

    filePrint(accelFile,"Acceleration Projection Error\n");
    filePrint(accelFile,"Frobenius: %1.6e\n", accelError.norm()/accelBuf.norm());
    filePrint(accelFile,"   Snapshot   |      L_inf      |      L1       |        L2   \n");
    for (int i = 0; i != accelSnapshots->numVec(); ++i){
      filePrint(accelFile,"     %d      ",i+1);
      for(int pnorm = 0; pnorm != 3; pnorm ++){
        if((pnorm) == 0){
            double Lperror = accelError.col(i).lpNorm<Eigen::Infinity>();
            Lperror = Lperror/accelBuf.col(i).lpNorm<Eigen::Infinity>()*100.;
            filePrint(accelFile,"     %1.6e ",Lperror);
        } else if(pnorm == 1){
            double Lperror = accelError.col(i).lpNorm<1>();
            Lperror = Lperror/accelBuf.col(i).lpNorm<1>()*100.;
            filePrint(accelFile,"    %1.6e ",Lperror);
        } else if(pnorm == 2) {
            double Lperror = accelError.col(i).lpNorm<2>();
            Lperror = Lperror/accelBuf.col(i).lpNorm<2>()*100.;
            filePrint(accelFile,"    %1.6e   ",(pnorm),Lperror);
        }
      }
      filePrint(accelFile,"\n");
    }
  }
#endif
}

} // end namespace Rom

Rom::DriverInterface *snapshotProjectionDriverNew(Domain *d) {
  return new Rom::SnapshotProjectionDriver(d);
}
