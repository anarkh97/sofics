#include "DistrROMPostProcessingDriver.h"

#include "DistrDomainUtils.h"
#include <Driver.d/Domain.h>
#include <Driver.d/DecDomain.h>

#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"
#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include "DistrMasterMapping.h"
#include "DistrVecNodeDof6Conversion.h"
#include "DistrNodeDof6Buffer.h"

#include <Feti.d/DistrVector.h>
#include <Utils.d/DistHelper.h>

#include "PtrPtrIterAdapter.h"

#include <Element.d/MpcElement.d/MpcElement.h>

#include <algorithm>
#include <memory>
#include <fstream>
#include <iostream>
#include <cassert>

class Rbm;
 
namespace Rom {

DistrROMPostProcessingDriver::DistrROMPostProcessingDriver(Domain *domain_) :
MultiDomainDynam(domain_),
normalizedBasis_(),
curState(NULL), fullDispBuffer(NULL), fullVelBuffer(NULL), fullAccBuffer(NULL),
fullVel2Buffer(NULL), fullDummyBuffer(NULL),
fullConstForceBuffer(NULL), fullExtForceBuffer(NULL), fullInertForceBuffer(NULL), fullResBuffer(NULL),
dummyDynOps(NULL)
{}

DistrROMPostProcessingDriver::~DistrROMPostProcessingDriver()
{
 if(curState) delete curState;
 if(fullDispBuffer) delete fullDispBuffer;
 if(fullVelBuffer) delete fullVelBuffer;
 if(fullAccBuffer) delete fullAccBuffer;
 if(fullVel2Buffer) delete fullVel2Buffer;
 if(fullDummyBuffer) delete fullDummyBuffer;
 if(fullConstForceBuffer) delete fullConstForceBuffer;
 if(fullExtForceBuffer) delete fullExtForceBuffer;
 if(fullInertForceBuffer) delete fullInertForceBuffer;
 if(fullResBuffer) delete fullResBuffer;
 if(dummyDynOps) delete dummyDynOps;
}

void
DistrROMPostProcessingDriver::preProcess() {

  {MultiDomainDynam::preProcess();
  bufferReducedFiles();
  //initialized decDomain class for use in projection basis preprocessing

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());

  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
  DistrNodeDof6Buffer buffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());

  normalizedBasis_.dimensionIs(projectionSubspaceSize, decDomain->masterSolVecInfo());

  for(int j=0; j<domain->solInfo().readInROBorModes.size(); ++j) {
    // read in distributed POD basis
    FileNameInfo fileInfo;
    std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD, j);
    if(domain->solInfo().newmarkBeta == 0 || domain->solInfo().useMassNormalizedBasis) fileName.append(".massorthonormalized");  
    DistrBasisInputFile podBasisFile(fileName);  //read in mass-normalized basis
    if(verboseFlag) filePrint(stderr, " ... Reading basis from file %s ...\n", fileName.c_str());

    int globalBasisSize = std::accumulate(domain->solInfo().localBasisSize.begin(), domain->solInfo().localBasisSize.end(), 0);
    if(globalBasisSize <= 0 || globalBasisSize == std::numeric_limits<int>::max()) {
      int sillyCounter = 0;
      for(std::vector<int>::iterator it = domain->solInfo().localBasisSize.begin(); it != domain->solInfo().localBasisSize.end(); it++) {
        std::string dummyName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD, sillyCounter); sillyCounter++;
        if(domain->solInfo().useMassNormalizedBasis) dummyName.append(".massorthonormalized");
        DistrBasisInputFile dummyFile(dummyName);
        *it = dummyFile.stateCount();
        globalBasisSize += *it;
      }
    }

    for (DistrVecBasis::iterator it = normalizedBasis_.begin() + std::accumulate(domain->solInfo().localBasisSize.begin(), domain->solInfo().localBasisSize.begin()+j, 0),
                                 it_end = normalizedBasis_.begin() + std::accumulate(domain->solInfo().localBasisSize.begin(), domain->solInfo().localBasisSize.end(), 0);
                                 it != it_end; ++it) {
      assert(podBasisFile.validCurrentState());

      podBasisFile.currentStateBuffer(buffer);
      converter.vector(buffer, *it);

      podBasisFile.currentStateIndexInc();
    }
  }

  filePrint(stderr, " ... Proj. Subspace Dimension = %-3d ...\n", projectionSubspaceSize);
  if(domain->solInfo().readInROBorModes.size() > 1)
    filePrint(stderr, " ... Number of Local Bases = %-3d    ...\n", domain->solInfo().readInROBorModes.size());
  
  VECsize = normalizedBasis_.size();}

  {//initialize multi domain dynamic post processor
  mddPostPro = MultiDomainDynam::getPostProcessor();

  //initialize containers for full coordinates 
  fullDispBuffer  = new GenDistrVector<double>(MultiDomainDynam::solVecInfo()); fullDispBuffer->zero();
  fullVelBuffer   = new GenDistrVector<double>(MultiDomainDynam::solVecInfo()); fullVelBuffer->zero();
  fullAccBuffer   = new GenDistrVector<double>(MultiDomainDynam::solVecInfo()); fullAccBuffer->zero();
  fullVel2Buffer  = new GenDistrVector<double>(MultiDomainDynam::solVecInfo()); fullVel2Buffer->zero();
  fullDummyBuffer = new GenDistrVector<double>(MultiDomainDynam::solVecInfo()); fullDummyBuffer->zero();

  //initialize system state vector container for use in multi domain dynamic post processor
  curState = new SysState<GenDistrVector<double> >( *fullDispBuffer, *fullVelBuffer, *fullAccBuffer, *fullVel2Buffer);}
}  //end preProcessing

void 
DistrROMPostProcessingDriver::bufferReducedFiles(){

  //get output information needed for parsing reduced data
  numConversionFiles = decDomain->getDomain()->solInfo().numRODFile;
  int skipTime = decDomain->getDomain()->solInfo().skipPodRom;
  //loop over reduced coordinate files
  for(int i = 0; i < numConversionFiles; i++) {
    //have all threads parse the reduced coordinate input file
    //there should be plenty of memory per node since projectionSubspaceSize is small
    std::ifstream reducedCoordFile(decDomain->getDomain()->solInfo().RODConversionFiles[i].c_str());
    if(reducedCoordFile.is_open()) {
      if(skipTime > 1) filePrint(stderr, " ... Skipping every %3d snapshots   ...\n", skipTime);

#ifdef ANDROID
      float time, dummyVar; // XXX DEBUG ANDROID
#else
      double time, dummyVar;
#endif
      int datatype, podsize, skipCounter;
      skipCounter = skipTime; // need to include t0 
      reducedCoordFile>>datatype; reducedCoordFile>>podsize;
      DataType.push_back(std::make_pair(datatype, podsize));
          switch(DataType[i].first) {
            case 0 :   // read reduced acceleration data
              {filePrint(stderr, " ... Buffering Reduced Acceleration Data ...\n");
              std::vector<double> timestamps;
              while(reducedCoordFile>>time) {
                if(skipCounter == skipTime) {
                  skipCounter = 1;
                  timestamps.push_back(time);
                  for(int j = 0; j < podsize; j++) {
                    reducedCoordFile>>dummyVar;
                    reducedAccBuffer.push_back(dummyVar);
                  }
                  filePrint(stderr,"\r ... Timestamp = %e ... ", time);
                } else {
                  for(int j = 0; j < podsize; j++) 
                    reducedCoordFile>>dummyVar;
                  skipCounter += 1;
                }
              }
              TimeStamps.push_back(timestamps);}
              filePrint(stderr,"\n");
              break;
            case 1 :   // read reduced displacement data
              {filePrint(stderr, " ... Buffering Reduced Displacement Data ...\n");
              std::vector<double> timestamps;
              while(reducedCoordFile>>time) {
                if(skipCounter == skipTime) {
                  skipCounter = 1;
                  timestamps.push_back(time);
                  for(int j = 0; j < podsize; j++) {
                    reducedCoordFile>>dummyVar;
                    reducedDispBuffer.push_back(dummyVar);
                  }
                  filePrint(stderr,"\r ... Timestamp = %e ... ", time);
                } else {
                  for(int j = 0; j < podsize; j++) 
                    reducedCoordFile>>dummyVar;
                  skipCounter += 1;
                }
              }
              TimeStamps.push_back(timestamps);}
              filePrint(stderr,"\n");
              break;
            case 2 :   // read reduced velocity data
              {filePrint(stderr, " ... Buffering Reduced Velocity Data ...\n");
              std::vector<double> timestamps;
              while(reducedCoordFile>>time) {
                if(skipCounter == skipTime){
                  skipCounter = 1;
                  timestamps.push_back(time);
                  for(int j = 0; j < podsize; j++) {
                    reducedCoordFile>>dummyVar;
                    reducedVelBuffer.push_back(dummyVar);
                  }
                  filePrint(stderr,"\r ... Timestamp = %e ...", time);
                } else {
                  for(int j = 0; j < podsize; j++)
                    reducedCoordFile>>dummyVar;
                  skipCounter += 1;
                }
              }
              TimeStamps.push_back(timestamps);}
              filePrint(stderr,"\n");
              break;
            default :
              filePrint(stderr, "\n... ROD conversion only supports Acceleration, Displacement, and Velocity ...\n");
          }
    } else {
      filePrint(stderr,"\nFailure to open file \n");
    }

    if(i != 0) { 
      if(DataType[i].second != DataType[i-1].second) {
        filePrint(stderr,"\n *** WARNING: Incompatible Input files %f %f\n", DataType[i-1].second, DataType[i].second);
        //exit(-1);
      }
    }

  } //end loop over input files, finished reading reduced data

  if(DataType.size() > 0) projectionSubspaceSize = DataType[0].second;
}

void
DistrROMPostProcessingDriver::solve() {

   preProcess();
   std::ofstream cvout(domain->solInfo().constraintViolationFile);

   bool computeResidual = false, computeExtForce = false, computeEnergies = false; 
   double residualNorm = 0, extForceNorm = 0;
   GenAssembler<double> * assembler = decDomain->getSolVecAssembler();
   OutputInfo *oinfo = geoSource->getOutputInfo();
   for (int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++) {
     if(oinfo[iOut].type == OutputInfo::RomResidual || oinfo[iOut].type == OutputInfo::RomResidual6) {
       computeResidual = true;
       filePrint(stderr, " ... Computing ROM Residual         ...\n");
       fullConstForceBuffer = new GenDistrVector<double>(MultiDomainDynam::solVecInfo());
       fullExtForceBuffer = new GenDistrVector<double>(decDomain->masterSolVecInfo());
       fullInertForceBuffer = new GenDistrVector<double>(MultiDomainDynam::solVecInfo());
       fullResBuffer = new GenDistrVector<double>(decDomain->masterSolVecInfo());
       getConstForce(*fullConstForceBuffer);
       break;
     }
   }
   for (int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++) {
     if(oinfo[iOut].type == OutputInfo::RomExtForce || oinfo[iOut].type == OutputInfo::RomExtForce6) {
       computeExtForce = true;
       filePrint(stderr, " ... Computing ROM External Force   ...\n");
       if(!computeResidual) {
         fullConstForceBuffer = new GenDistrVector<double>(MultiDomainDynam::solVecInfo());
         fullExtForceBuffer = new GenDistrVector<double>(decDomain->masterSolVecInfo());
         getConstForce(*fullConstForceBuffer);
       }
       break;
     }
   }
   for (int iOut = 0; iOut < geoSource->getNumOutInfo(); iOut++) {
     if(oinfo[iOut].type == OutputInfo::Energies) {
       computeEnergies = true;
       if(!computeResidual && !computeExtForce) {
         fullConstForceBuffer = new GenDistrVector<double>(MultiDomainDynam::solVecInfo());
         fullExtForceBuffer = new GenDistrVector<double>(decDomain->masterSolVecInfo());
         getConstForce(*fullConstForceBuffer);
       }
       break;
     }
   }
   int counter = 0; //TODO: make this portion more general so it doesn't depend on the assumption
                    //that all files have matching timestamps
   if(TimeStamps.size() > 0)
   for(std::vector<double>::iterator it = TimeStamps[0].begin(); it != TimeStamps[0].end(); it++) {

     // load current state for output 
     for(int i = 0; i < numConversionFiles; i++) {
       std::vector<double> buffer(projectionSubspaceSize);
           switch(DataType[i].first) {
            case 0 :

              for (int j = 0; j < projectionSubspaceSize; j++) 
                buffer[j] = reducedAccBuffer[counter*projectionSubspaceSize+j];
             
              normalizedBasis_.expand(buffer, *fullAccBuffer);
              break;
            case 1 :

              for (int j = 0; j < projectionSubspaceSize; j++)
                buffer[j] = reducedDispBuffer[counter*projectionSubspaceSize+j];

              normalizedBasis_.expand(buffer, *fullDispBuffer);
              break;
            case 2 :
              if(counter != 0)
                *fullVel2Buffer = *fullVelBuffer;

              for (int j = 0; j < projectionSubspaceSize; j++) 
                buffer[j] = reducedVelBuffer[counter*projectionSubspaceSize+j];

              normalizedBasis_.expand(buffer, *fullVelBuffer);
              break;
            default :
              break; 
              //nothing
           }
     }

     geomState->explicitUpdate(decDomain, *fullDispBuffer);
     geomState->setVelocity(*fullVelBuffer, 2);
     geomState->setAcceleration(*fullAccBuffer, 2);
     execParal(decDomain->getNumSub(), this, &DistrROMPostProcessingDriver::subUpdateStates, *it);

     if(computeResidual || computeExtForce || computeEnergies) {
       if(!dummyDynOps) {
         dummyDynOps = buildOps(1,0,0);
       }
       // compute non-follower external forces
       computeExtForce2(*curState, *fullExtForceBuffer, *fullConstForceBuffer, counter, *it);
       if(computeResidual) {
         // compute linear inertial forces M*a
         dummyDynOps->M->mult(*fullAccBuffer, *fullInertForceBuffer);
         // compute internal forces and other configuration-dependent forces including external follower forces
         // and nonlinear inertial force correction
         getInternalForce(*fullDispBuffer, *fullResBuffer, *it, counter);
         // add all forces and compute norm
         *fullResBuffer -= *fullExtForceBuffer;
         *fullResBuffer += *fullInertForceBuffer;
         if(domain->solInfo().romresidType == 0) geomState->transform(*fullResBuffer, 2, true); // multiply by T^{-1}
         assembler->assemble(*fullResBuffer);
         residualNorm += fullResBuffer->sqNorm();
       }
       if(computeExtForce) {
         // compute external follower forces and add to fullExtForceBuffer
         getFollowerForce(*fullExtForceBuffer, *it, counter);
         if(domain->solInfo().romresidType == 0) geomState->transform(*fullExtForceBuffer, 2, true); // multiply by T^{-1}
         assembler->assemble(*fullExtForceBuffer);
         extForceNorm += fullExtForceBuffer->sqNorm();
       }
     } else
     if(!dummyDynOps) dummyDynOps = new MDDynamMat;

     geomState->transform(*fullVelBuffer, 0, true); // transform angular velocity from the 1st time derivative
                                                    // of the (unscaled) total rotation vector to convected
     geomState->transform(*fullAccBuffer, 4, true); // transform angular acceleration from the 2nd time derivative
                                                    // of the (unscaled) total rotation vector to convected

     mddPostPro->dynamOutput(counter, *it, *dummyDynOps, ((fullExtForceBuffer) ? *fullExtForceBuffer : *fullDummyBuffer),
                             fullDummyBuffer, *curState, fullResBuffer);

     if(geoSource->getNumConstraintElementsIeq() && decDomain->getGlobalNumSub() == 1) { // output the constraint violation
       double err = 0;
       Elemset &eset = geoSource->getPackedEsetConstraintElementIeq();
       for(int i=0; i<geoSource->getNumConstraintElementsIeq(); ++i) {
         MpcElement *ele = static_cast<MpcElement*>(eset[i]);
         if(counter==0) ele->renum(decDomain->getSubDomain(0)->getGlobalToLocalNode());
         ele->update((*geomState)[0], *((*geomState)[0]), decDomain->getSubDomain(0)->getNodes(), *it);
         err = std::max(err,ele->getError(*((*geomState)[0])));
       }
       cvout << *it << " " << err << std::endl;
     }

     filePrint(stderr,"\r ... ROM Conversion Loop: t = %9.3e, %3d%% complete ...",
                *it, int(*it/(TimeStamps[0].back())*100));

     counter += 1;
   }  //end of loop over time stamps

   if(computeResidual)
     filePrint(stderr, "\n ... ROM Residual Norm = %10.4e ...\n", sqrt(residualNorm));
   if(computeExtForce)
     filePrint(stderr, " ... Ext. Force Norm   = %10.4e ...\n", sqrt(extForceNorm));
}

void
DistrROMPostProcessingDriver::subUpdateStates(int i, double time)
{
  decDomain->getSubDomain(i)->updateStates((*geomState)[i],*((*geomState)[i]),allCorot[i],time);  
}

} //end namespace Rom

Rom::DriverInterface *distrROMPostProcessingDriverNew(Domain *D) {
  return new Rom::DistrROMPostProcessingDriver(D);
}
