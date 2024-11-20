#include "DistrExplicitPodProjectionNonLinDynamicBase.h"

#include "DistrGalerkinProjectionSolver.h"

#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include "DistrMasterMapping.h"
#include "DistrVecNodeDof6Conversion.h"
#include "DistrNodeDof6Buffer.h"

#include "PtrPtrIterAdapter.h"

#include <Threads.d/PHelper.h>

#include <Driver.d/DecDomain.h>
#include <Feti.d/DistrVector.h>
#include <Timers.d/StaticTimers.h>
#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <limits>
#include <Corotational.d/utilities.h>
#include <Corotational.d/DistrGeomState.h>
#include <Driver.d/SysState.h>

extern Communicator *structCom;
extern int verboseFlag;

namespace Rom {

MultiDomDynPodPostProcessor::MultiDomDynPodPostProcessor(DecDomain *d, StaticTimers* times, DistrGeomState *geomState = 0, Corotator ***allCorot = 0) :
    MultiDomDynPostProcessor(d, times, geomState, allCorot),
    DispSensorValues(NULL),
    AccSensorValues(NULL),
    VelSensorValues(NULL),
    all_cdsa(NULL)
{
  oinfo = geoSource->getOutputInfo();
  numOutInfo = geoSource->getNumOutInfo();

  DispSensor = false;
  AccSensor  = false;
  VelSensor  = false;

  for (int iOut = 0; iOut < numOutInfo; iOut++) {
    switch(oinfo[iOut].type) {
      case OutputInfo::Accel6 : case OutputInfo::Acceleration :
        if(oinfo[iOut].nodeNumber != -1) {
          AccSensor = true;
          if(oinfo[iOut].type == OutputInfo::Accel6 && (oinfo[iOut].angularouttype != OutputInfo::total || oinfo[iOut].rescaling)) {
            // in this case rotation vector and angular velocity are required for angular acceleration transformation
            DispSensor = true;
            VelSensor = true;
          }
        } else {
          if(!structCom || structCom->myID() == 0) {
            oinfo[iOut].filptr = fopen(oinfo[iOut].filename, "wb");
            if(!oinfo[iOut].filptr) {
              fprintf(stderr," *** ERROR: Cannot open %s, exiting...\n", oinfo[iOut].filename);
              exit(-1);
            }
          }
          filePrint(oinfo[iOut].filptr, "0\n"); 
        }
        break;
      case OutputInfo::Disp6DOF : case OutputInfo::Displacement :
        if(oinfo[iOut].nodeNumber != -1) {
          DispSensor = true;
        } else {
          if(!structCom || structCom->myID() == 0) {
            oinfo[iOut].filptr = fopen(oinfo[iOut].filename, "wb");
            if(!oinfo[iOut].filptr) {
              fprintf(stderr," *** ERROR: Cannot open %s, exiting...\n", oinfo[iOut].filename);
              exit(-1);
            }
          }
          filePrint(oinfo[iOut].filptr, "1\n");
        }
        break;
      case OutputInfo::Velocity6 : case OutputInfo::Velocity :
        if(oinfo[iOut].nodeNumber != -1) {
          VelSensor = true;
          if(oinfo[iOut].type == OutputInfo::Velocity6 && (oinfo[iOut].angularouttype != OutputInfo::total || oinfo[iOut].rescaling)) {
            // in this case rotation vector is required for angular velocity transformation
            DispSensor = true;
          }
        } else {
          if(!structCom || structCom->myID() == 0) {
            oinfo[iOut].filptr = fopen(oinfo[iOut].filename, "wb");
            if(!oinfo[iOut].filptr) {
              fprintf(stderr," *** ERROR: Cannot open %s, exiting...\n", oinfo[iOut].filename);
              exit(-1);
            }
          }
          filePrint(oinfo[iOut].filptr, "2\n");
        }
        break; 
      default:
        filePrint(stderr, " *** WARNING: Online ROM output only supports GDISPLAC, GVELOCIT, and GACCELER.\n");
        filePrint(stderr, "     output type selected is %d \n", oinfo[iOut].type);
    }
  }

  if(DispSensor || VelSensor || AccSensor) {
#ifdef DISTRIBUTED
    for(int iSub = 0; iSub < decDomain->getNumSub(); iSub++) {
      int nOutNodes = decDomain->getSubDomain(iSub)->getNumNodalOutput();
      if(nOutNodes) {
        geoSource->openOutputFiles(decDomain->getSubDomain(iSub)->getOutputNodes(),
                                   decDomain->getSubDomain(iSub)->getOutIndex(), nOutNodes);
      }
    }
#else
    geoSource->openSensorOutputFiles();
#endif
  }

  if(structCom) structCom->sync();
  nodeVector.resize(decDomain->getNumSub());
  execParal(decDomain->getNumSub(), this, &MultiDomDynPodPostProcessor::subBuildSensorNodeVector);
}

MultiDomDynPodPostProcessor::~MultiDomDynPodPostProcessor() {

  if(all_cdsa) delete [] all_cdsa;
  if(DispSensorValues) delete DispSensorValues;
  if(AccSensorValues) delete AccSensorValues;
  if(VelSensorValues) delete VelSensorValues;
}

void
MultiDomDynPodPostProcessor::printPODSize(int PODsize) {

  podSize = PODsize;

  for (int iOut = 0; iOut < numOutInfo; iOut++) {
    switch(oinfo[iOut].type) {
      case OutputInfo::Accel6 : case OutputInfo::Acceleration :
        if(oinfo[iOut].nodeNumber == -1)
          filePrint(oinfo[iOut].filptr, "%d\n", PODsize);
        break;
      case OutputInfo::Disp6DOF : case OutputInfo::Displacement :
        if(oinfo[iOut].nodeNumber == -1)
          filePrint(oinfo[iOut].filptr, "%d\n", PODsize);
        break;
      case OutputInfo::Velocity6 : case OutputInfo::Velocity :
        if(oinfo[iOut].nodeNumber == -1)
          filePrint(oinfo[iOut].filptr, "%d\n", PODsize);
        break;
      default:
        break;
    }
  }
}

void
MultiDomDynPodPostProcessor::makeSensorBasis(DistrVecBasis *fullBasis) {

  all_cdsa = new DofSetArray * [decDomain->getNumSub()];

  for(int i = 0; i < decDomain->getNumSub(); ++i) all_cdsa[i] = decDomain->getSubDomain(i)->getCDSA();

  SensorBasis = fullBasis;

  SensorBasis->makeSparseBasis2(nodeVector, all_cdsa);
  // allocate space for containers to hold sensor values
  if(DispSensor) {
    DispSensorValues = new GenDistrVector<double>(SensorBasis->vectorInfo());
    DispSensorValues->zero();
  }
  if(AccSensor) {
    AccSensorValues = new GenDistrVector<double>(SensorBasis->vectorInfo());
    AccSensorValues->zero();
  }
  if(VelSensor) {
    VelSensorValues = new GenDistrVector<double>(SensorBasis->vectorInfo());
    VelSensorValues->zero();
  }
}
  
void
MultiDomDynPodPostProcessor::subBuildSensorNodeVector(int iSub) {

  std::vector<int> &subSensorNodes = nodeVector[iSub];

  //load vector of sensor nodes converted to local numbering
  for (int iOut = 0; iOut < numOutInfo; iOut++) {
    if(oinfo[iOut].nodeNumber != -1) {
      int locNode = decDomain->getSubDomain(iSub)->globalToLocal(oinfo[iOut].nodeNumber);
      if(locNode > -1) 
        subSensorNodes.push_back(locNode);
    }
  }

  //if multiple outputs are requested for a single node, cull redundancies from node vector
  std::sort(subSensorNodes.begin(), subSensorNodes.end());
  std::vector<int>::iterator packedNodeIt = std::unique(subSensorNodes.begin(), subSensorNodes.end());
  subSensorNodes.resize(packedNodeIt-subSensorNodes.begin());
}

double
MultiDomDynPodPostProcessor::subGetPrescribedSensorValue(int iSub, int locNode, int j, int bctype)
{
  int dof = -1;
  switch(j) {
    case 0 : dof = decDomain->getSubDomain(iSub)->getDSA()->locate(locNode, DofSet::Xdisp); break;
    case 1 : dof = decDomain->getSubDomain(iSub)->getDSA()->locate(locNode, DofSet::Ydisp); break;
    case 2 : dof = decDomain->getSubDomain(iSub)->getDSA()->locate(locNode, DofSet::Zdisp); break;
    case 3 : dof = decDomain->getSubDomain(iSub)->getDSA()->locate(locNode, DofSet::Xrot); break;
    case 4 : dof = decDomain->getSubDomain(iSub)->getDSA()->locate(locNode, DofSet::Yrot); break;
    case 5 : dof = decDomain->getSubDomain(iSub)->getDSA()->locate(locNode, DofSet::Zrot); break;
  }

  double val = 0;
  if(dof != -1) {
    switch(bctype) {
      case 0 : val = decDomain->getSubDomain(iSub)->getBcx()[dof]; break;
      case 1 : val = decDomain->getSubDomain(iSub)->getVcx()[dof]; break;
      case 2 : val = decDomain->getSubDomain(iSub)->getAcx()[dof]; break;
    }
  }
  return val;
}

void
MultiDomDynPodPostProcessor::subPrintSensorValues(int iSub, GenDistrVector<double> &SensorData, OutputInfo *OINFO, double *time, int bctype) {

  // XXX prescribed values in vcx and acx are convected.
#ifdef DISTRIBUTED
  int locNode = decDomain->getSubDomain(iSub)->globalToLocal(OINFO->nodeNumber);
  if(locNode > -1) { // if node is -1, sensor node is not in this subdomain, don't print
    for(int k = 0; k < decDomain->getSubDomain(iSub)->getNumNodalOutput(); ++k) { // not the best way
      if(locNode == decDomain->getSubDomain(iSub)->getOutputNodes()[k]) {
#else
  int subI = decDomain->getGlSubToLocal()[(*decDomain->getNodeToSub())[OINFO->nodeNumber][0]];
  if(subI == iSub) {
        int locNode = decDomain->getSubDomain(iSub)->globalToLocal(OINFO->nodeNumber);
#endif
        int ndofs;
        int dofs[6];
        if(OINFO->type == OutputInfo::Disp6DOF || OINFO->type == OutputInfo::Velocity6 || OINFO->type == OutputInfo::Accel6) {
          ndofs = all_cdsa[iSub]->number(locNode, DofSet::XYZdisp | DofSet::XYZrot, dofs);
        }
        else {
          ndofs = all_cdsa[iSub]->number(locNode, DofSet::XYZdisp, dofs);
        }
        double *data = new double[ndofs];
        for(int j = 0; j < ndofs; ++j) {
          data[j] = (dofs[j] != -1) ? SensorData.subData(decDomain->getSubDomain(iSub)->localSubNum())[dofs[j]] 
                                    : subGetPrescribedSensorValue(iSub, locNode, j, bctype);
        }
        // transform rotation vector, if necessary
        if(OINFO->type == OutputInfo::Disp6DOF && (OINFO->rotvecouttype != OutputInfo::Euler || OINFO->rescaling)) {
          double rten[3][3], psi[3] = { data[3], data[4], data[5] };
          if(OINFO->rescaling) vec_to_mat(psi, rten);
          tran_rvec(rten, psi, OINFO->rescaling, OINFO->rotvecouttype, data+3);
        }
        // transform angular velocity, if necessary
        else if(OINFO->type == OutputInfo::Velocity6 && (OINFO->angularouttype != OutputInfo::total || OINFO->rescaling)) {
          double rten[3][3], psi[3], psidot[3] = { data[3], data[4], data[5] };
          for(int j = 0; j < 3; ++j) psi[j] = (dofs[3+j] != -1) ? DispSensorValues->subData(decDomain->getSubDomain(iSub)->localSubNum())[dofs[3+j]]
                                                                : subGetPrescribedSensorValue(iSub, locNode, 3+j, 0);
          if(OINFO->rescaling) vec_to_mat(psi, rten);
          tran_veloc(rten, psi, psidot, 2, OINFO->angularouttype, false, OINFO->rescaling, data+3);
        }
        // transform angular acceleration, if necessary
        else if(OINFO->type == OutputInfo::Accel6 && (OINFO->angularouttype != OutputInfo::total || OINFO->rescaling)) {
          double rten[3][3], psi[3], psidot[3], psiddot[3] = { data[3], data[4], data[5] };
          for(int j = 0; j < 3; ++j) {
            psi[j] = (dofs[3+j] != -1) ? DispSensorValues->subData(decDomain->getSubDomain(iSub)->localSubNum())[dofs[3+j]]
                                       : subGetPrescribedSensorValue(iSub, locNode, 3+j, 0);
            psidot[j] = (dofs[3+j] != -1) ? VelSensorValues->subData(decDomain->getSubDomain(iSub)->localSubNum())[dofs[3+j]]
                                          : subGetPrescribedSensorValue(iSub, locNode, 3+j, 1);
          }
          if(OINFO->rescaling) vec_to_mat(psi, rten);
          tran_accel(rten, psi, psidot, psiddot, 2, OINFO->angularouttype, false, OINFO->rescaling, data+3);
        }
        int w = OINFO->width;
        int p = OINFO->precision; 
        fprintf(OINFO->filptr, "  % *.*E  ", w, p, *time);
        for(int j = 0; j < ndofs; ++j) {
          fprintf(OINFO->filptr, " % *.*E", w, p, data[j]);
        }
        fprintf(OINFO->filptr, "\n");
        fflush(OINFO->filptr);
        delete [] data;
#ifdef DISTRIBUTED
        break;
      }
    }
#endif
  }
}

void
MultiDomDynPodPostProcessor::dynamOutput(int tIndex, double t, MDDynamMat &dynOps, DistrVector &distForce,
                                         DistrVector *distAeroF, SysState<DistrVector>& distState) {

  //all MPI processes have a full copy of reduced coordinates, only master processes needs to print
  if( structCom && distState.getDisp().size() == 0 ) return;

  int p = std::numeric_limits<double>::digits10+1;

  bool DispProjected = false, AccProjected = false, VelProjected = false;

  for(int iOut = 0; iOut < numOutInfo; iOut++) {

    if(tIndex % oinfo[iOut].interval == 0) {

      switch(oinfo[iOut].type) {
         case OutputInfo::Accel6 : case OutputInfo::Acceleration :
           {
             if(oinfo[iOut].nodeNumber == -1) {
               filePrint(oinfo[iOut].filptr, "   %.*e\n", p, t); // print timestamp
               for(int i = 0; i < podSize; i++) {
                 filePrint(oinfo[iOut].filptr, "%.*e ", p, distState.getAccel()[i]);
               }
               filePrint(oinfo[iOut].filptr, "\n");
             }
             else {
               if(!AccProjected) {
                 SensorBasis->expand2(distState.getAccel(), *AccSensorValues);
                 AccProjected = true;
               }
               if(oinfo[iOut].type == OutputInfo::Accel6 && (oinfo[iOut].angularouttype != OutputInfo::total || oinfo[iOut].rescaling)) {
                 // in this case rotation vector and angular velocity are required for angular acceleration transformation
                 if(!DispProjected) {
                   SensorBasis->expand2(distState.getDisp(), *DispSensorValues);
                   DispProjected = true;
                 }
                 if(!VelProjected) {
                   SensorBasis->expand2(distState.getVeloc(), *VelSensorValues);
                   VelProjected = true;
                 }
               }
               execParal(decDomain->getNumSub(), this, &MultiDomDynPodPostProcessor::subPrintSensorValues, *AccSensorValues, &oinfo[iOut], &t, 2);
             }
           }
           break;
         case OutputInfo::Disp6DOF : case OutputInfo::Displacement : 
           {
             if(oinfo[iOut].nodeNumber == -1) {
               filePrint(oinfo[iOut].filptr, "   %.*e\n", p, t); // print timestamp
               for(int i = 0; i < podSize; i++) {
                 filePrint(oinfo[iOut].filptr, "%.*e ", p, distState.getDisp()[i]);
               }
               filePrint(oinfo[iOut].filptr, "\n");
             }
             else {
               if(!DispProjected) {
                 SensorBasis->expand2(distState.getDisp(), *DispSensorValues);
                 DispProjected = true;
               }
               execParal(decDomain->getNumSub(), this, &MultiDomDynPodPostProcessor::subPrintSensorValues, *DispSensorValues, &oinfo[iOut], &t, 0);
             }
           }
           break;
         case OutputInfo::Velocity6 : case OutputInfo::Velocity :
           {
             if(oinfo[iOut].nodeNumber == -1) {
               filePrint(oinfo[iOut].filptr, "   %.*e\n", p, t); // print timestamp
               for(int i = 0; i < podSize; i++) {
                 filePrint(oinfo[iOut].filptr, "%.*e ", p, distState.getVeloc()[i]);
               }
               filePrint(oinfo[iOut].filptr, "\n");
             }
             else {
               if(!VelProjected) {
                 SensorBasis->expand2(distState.getVeloc(), *VelSensorValues);
                 VelProjected = true;
               }
               if(oinfo[iOut].type == OutputInfo::Velocity6 && (oinfo[iOut].angularouttype != OutputInfo::total || oinfo[iOut].rescaling)) {
                 // in this case rotation vector is required for angular velocity transformation
                 if(!DispProjected) {
                   SensorBasis->expand2(distState.getDisp(), *DispSensorValues);
                   DispProjected = true;
                 }
               }
               execParal(decDomain->getNumSub(), this, &MultiDomDynPodPostProcessor::subPrintSensorValues, *VelSensorValues, &oinfo[iOut], &t, 1);
             }
           }
           break;
         default:
           filePrint(stderr, " ... ROM output only supports Acceleration, Displacement, and Velocity ... \n");
      }
    }
  }
}

MultiDomDynPodPostProcessor *
DistrExplicitPodProjectionNonLinDynamicBase::getPostProcessor() {

  mddPostPro = new MultiDomDynPodPostProcessor(decDomain, times, geomState, allCorot);
  mddPostPro->printPODSize(normalizedBasis_.numVectors());
  mddPostPro->makeSensorBasis(&normalizedBasis_);

  return mddPostPro;
}

DistrExplicitPodProjectionNonLinDynamicBase::DistrExplicitPodProjectionNonLinDynamicBase(Domain *_domain) :
  MultiDomainDynam(_domain),
  domain(_domain),
  stableCount(0),
  localBasisId(-1)
{}

void
DistrExplicitPodProjectionNonLinDynamicBase::preProcess() {
  { 
    MultiDomainDynam::preProcess();
    //preProcessing for local/global reduced order basis
    //check how many vectors are in each file for memory allocation
    //if pod size is given as 0, then use all vectors in basis file
    for(int rob=0; rob<domain->solInfo().readInROBorModes.size(); ++rob) {
      FileNameInfo fileInfo;
      std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD,rob);
      fileName.append(".massorthonormalized");
      DistrBasisInputFile BasisFile(fileName);
      int locSize = domain->solInfo().localBasisSize[rob] ?
                    std::min(domain->solInfo().localBasisSize[rob], BasisFile.stateCount()) :
                    BasisFile.stateCount();
      locBasisVec.push_back(locSize);
      if(verboseFlag && domain->solInfo().readInROBorModes.size()>1) filePrint(stderr, " ... Local Basis %d size %-3d ...\n",rob,locBasisVec[rob]);
    }

    //initialize helper objects for reading in distributed basis vectors
    const int projectionSubspaceSize = std::accumulate(locBasisVec.begin(),locBasisVec.end(),0);
    filePrint(stderr, " ... Proj. Subspace Dimension = %-3d ...\n", projectionSubspaceSize);
    normalizedBasis_.dimensionIs(projectionSubspaceSize, decDomain->masterSolVecInfo());

    DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());

    typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
    DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                    SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
    DistrNodeDof6Buffer buffer(masterMapping.masterNodeBegin(), masterMapping.masterNodeEnd());

    // this loop reads in vectors and stores them in a single Distributed Basis structure
    DistrVecBasis::iterator it = normalizedBasis_.begin(); 
    for(int rob=0; rob<domain->solInfo().readInROBorModes.size(); ++rob){
      FileNameInfo fileInfo; 
      std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD,rob);
      fileName.append(".massorthonormalized");
      DistrBasisInputFile BasisFile(fileName); //read in mass-normalized basis
      if(verboseFlag) filePrint(stderr, " ... Reading basis from file %s ...\n", fileName.c_str());

      for (int currentVec = 0; currentVec<locBasisVec[rob]; ++currentVec) {
        assert(BasisFile.validCurrentState());

        BasisFile.currentStateBuffer(buffer);
        converter.vector(buffer, *it);
        it++;
      
        BasisFile.currentStateIndexInc();
      }
    }

    // if performing local basis analysis and centroids are provided, read in centroids for basis switching
    if(domain->solInfo().readInLocalBasesCent.size() >= 1) {
      centroids.dimensionIs(locBasisVec.size(), decDomain->masterSolVecInfo()); 
      DistrVecBasis::iterator it = centroids.begin();
      for(int rob=0; rob<locBasisVec.size(); ++rob) { 
        std::string fileName = domain->solInfo().readInLocalBasesCent[rob].c_str();
        DistrBasisInputFile CentroidFile(fileName);
        if(verboseFlag) filePrint(stderr, " ... Reading centroid from file %s ...\n", fileName.c_str());
        assert(CentroidFile.validCurrentState());
        CentroidFile.currentStateBuffer(buffer);
        converter.vector(buffer, *it);
        it++;
      }
    }
    // if auxillary bases are provided, read in data structures for fast switching 
    // every process owns a copy of VtV, w, and d
    if(domain->solInfo().readInLocalBasesAuxi.size() > 0) {
#ifdef USE_EIGEN3
      const int Nv = domain->solInfo().readInROBorModes.size();
      const int k = std::accumulate(domain->solInfo().localBasisSize.begin(), domain->solInfo().localBasisSize.end(), 0);
      VtV.resize(Nv,Nv); //grammian
      d.resize(Nv,Nv);
      w.resize(Nv,Nv);
      for(int i=0; i<Nv; ++i) {
        int ki = domain->solInfo().localBasisSize[i];
        for(int j=i+1; j<Nv; ++j) {
          int kj = domain->solInfo().localBasisSize[j];
          std::string fileName = domain->solInfo().readInLocalBasesAuxi[std::make_pair(i,j)];
          std::ifstream file(fileName.c_str());
          VtV(i,j).resize(ki,kj);
          for(int irow = 0; irow < ki; ++irow) {
            for(int jcol = 0; jcol < kj; ++jcol) {
              file >> VtV(i,j)(irow,jcol);
            }
          }
          file >> d(i,j);
          w(i,j).resize(k);
          for(int irow=0; irow<k; ++irow) {
            file >> w(i,j)[irow];
          }
          file.close();
        }
      }
#endif    
    }
  }

  //preProcessing for solution vector information
  {
   reducedInfo.domLen = new int[MultiDomainDynam::solVecInfo().numDom]; 
   reducedInfo.numDom = MultiDomainDynam::solVecInfo().numDom;
   int totLen = 0;
   for(int iSub = 0; iSub < MultiDomainDynam::solVecInfo().numDom; ++iSub) {
     reducedInfo.domLen[iSub] = (iSub==0) ? normalizedBasis_.numVec() : 0;
     totLen += reducedInfo.domLen[iSub];
   }

   reducedInfo.len = totLen;
   reducedInfo.setMasterFlag();
  }
  
  //preProcessing for dummy working variables
  {fExt      = new DistrVector(MultiDomainDynam::solVecInfo());
  fInt       = new DistrVector(MultiDomainDynam::solVecInfo());
  cnst_fBig  = new DistrVector(MultiDomainDynam::solVecInfo());
  aero_fBig  = new DistrVector(MultiDomainDynam::solVecInfo());
  d_n        = new DistrVector(MultiDomainDynam::solVecInfo());
  v_n        = new DistrVector(MultiDomainDynam::solVecInfo());
  a_n        = new DistrVector(MultiDomainDynam::solVecInfo());
  v_p        = new DistrVector(MultiDomainDynam::solVecInfo());
  tempVec    = new DistrVector(MultiDomainDynam::solVecInfo());
  ctc_f      = new DistrVector(MultiDomainDynam::solVecInfo());
  dummyState = new SysState<DistrVector>(*d_n, *v_n, *a_n, *v_p);
  times      = new StaticTimers;}

  decDomain->assembleNodalInertiaTensors(melArray);
  if(domain->solInfo().stable && !domain->solInfo().elemLumpPodRom)
    execParal(decDomain->getNumSub(), this, &DistrExplicitPodProjectionNonLinDynamicBase::subInitStiff);
}

const DistrInfo &
DistrExplicitPodProjectionNonLinDynamicBase::solVecInfo() {
  return reducedVecInfo();
}

DistrInfo &
DistrExplicitPodProjectionNonLinDynamicBase::reducedVecInfo() {
  return reducedInfo; //ROB size
}

void 
DistrExplicitPodProjectionNonLinDynamicBase::printFullNorm(DistrVector &v) {

  normalizedBasis_.expand(v,*tempVec);

  filePrint(stderr,"%1.4e\n",tempVec->norm());
}

void
DistrExplicitPodProjectionNonLinDynamicBase::getInitState(SysState<DistrVector>& curState) {

  DistrVector &dr_n = curState.getDisp();
  DistrVector &vr_n = curState.getVeloc();
  DistrVector &ar_n = curState.getAccel();
  DistrVector &vr_p = curState.getPrevVeloc();

  dr_n.zero();
  vr_n.zero();
  ar_n.zero();
  vr_p.zero();

  int numIDisModal = domain->numInitDispModal();
  if(numIDisModal) {
    filePrint(stderr, " ... Using Modal IDISPLACEMENTS     ...\n");
    BCond* iDisModal = domain->getInitDispModal();
    for(int i = 0; i < numIDisModal; ++i) {
      if(iDisModal[i].nnum < dr_n.size())
        dr_n[iDisModal[i].nnum] = iDisModal[i].val;
    }
    normalizedBasis_.expand( dr_n, *d_n);
    geomState->update(*d_n);
    initLocalBasis(dr_n);
  }

  int numIVelModal = domain->numInitVelocityModal();
  if(numIVelModal) {
    filePrint(stderr, " ... Using Modal IVELOCITIES        ...\n");
    if(!numIDisModal) initLocalBasis(dr_n);
    BCond* iVelModal = domain->getInitVelocityModal();
    int numLocalIVelModal = normalizedBasis_.numVec();   // returns blockCols_
    int startIVelModal    = normalizedBasis_.startCol(); // returns startCol_
    for(int i = startIVelModal; i < startIVelModal+numLocalIVelModal; ++i) {
      if(iVelModal[i].nnum < vr_n.size())
        vr_n[iVelModal[i].nnum] = iVelModal[i].val;
    }
    normalizedBasis_.expand( vr_n, *v_n);
    geomState->setVelocity(*v_n);
  }

  // XXX currently, if modal initial conditions are defined then any non-modal initial conditions are ignored
  if(numIDisModal == 0 && numIVelModal == 0) {

    MultiDomainDynam::getInitState(*dummyState);

    MultiDomainDynam::buildOps(1.0, 0.0, 0.0);

    if(d_n->norm() != 0) reduceDisp(*d_n, dr_n);
    initLocalBasis(dr_n);
    if(v_n->norm() != 0) reduceDisp(*v_n, vr_n);
    if(a_n->norm() != 0) reduceDisp(*a_n, ar_n);
    if(v_p->norm() != 0) reduceDisp(*v_p, vr_p);
  }
}

//start local basis functions
void
DistrExplicitPodProjectionNonLinDynamicBase::initLocalBasis(DistrVector &q){
#if defined(USE_EIGEN3) && ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11))
  if(domain->solInfo().readInLocalBasesCent.size() >= 1 || VtV.size() > 0) { // only execute this if doing local ROM
    localBasisId  = selectLocalBasis(q); 
    if(verboseFlag) filePrint(stderr," ... Initial Local Basis # %d        ...\n",localBasisId);
    int blockCols = locBasisVec[localBasisId];
    int startCol  = std::accumulate(locBasisVec.begin(),locBasisVec.begin()+localBasisId,0);
    normalizedBasis_.localBasisIs(startCol,blockCols);
    setLocalReducedMesh(localBasisId);
  }
#else
  filePrint(stderr,"*** Error: Local auxiliary bases requres Aero-S to be built with eigen library and c++11 support.\n");
  exit(-1);
#endif
}

int
DistrExplicitPodProjectionNonLinDynamicBase::selectLocalBasis(DistrVector &q){

  int cc = std::numeric_limits<int>::max(); 
  if(centroids.numVec() > 0) { // modelII: slow implementation using cluster centroids
    double minDist = std::numeric_limits<double>::max(); 
    for(int c=0; c<locBasisVec.size(); ++c){
      DistrVector centroidVec = centroids[c]; 
      DistrVector q_Big(MultiDomainDynam::solVecInfo());
      normalizedBasis_.expand(q,q_Big);
      centroidVec -= q_Big;
      double distNorm = centroidVec.norm();
      if(distNorm < minDist) { // check proximity to centroid, update if closer
        minDist = distNorm;
        cc      = c; 
      }
    }
  }
#if defined(USE_EIGEN3) && ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11))
  else if(d.size() > 0 && w.size() > 0 ) { // modelIII: fast implementation using pre-computed auxiliary quantities
    const int Nv = domain->solInfo().readInROBorModes.size(); // get number of local bases
    const int k  = q.size();                                  // get total number of vectors
    Eigen::Map<Eigen::VectorXd> qi(q.data(),k);               // get working array
    std::vector<int> s(Nv); for(int m=0; m<Nv; ++m) s[m] = m; // fill array with 1 through number of bases 
    std::sort(s.begin(), s.end(), [&](int m, int p) {
      return (p>m && (w(m,p).dot(qi) + d(m,p)) < 0) || (m>p && (w(p,m).dot(qi) + d(p,m)) > 0);
    });
    return s[0];
  }
#endif
  else {
    filePrint(stderr,"*** Error: cluster centroids or auxiliary quantities required to select local basis.\n");
    exit(-1);
  }
  return cc; 
}

void
DistrExplicitPodProjectionNonLinDynamicBase::setLocalBasis(DistrVector &q, DistrVector &qd){
#if defined(USE_EIGEN3) && ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11))
  if(domain->solInfo().readInLocalBasesCent.size() >= 1 || VtV.size() > 0) {
    int j = selectLocalBasis(q);
    if(j != localBasisId) { // if a new local basis has been selected, update the things
      if(verboseFlag) filePrint(stderr," ... Selecting Local Basis # %d     ...\n",j);
      int blockCols = locBasisVec[j]; 
      int startCol  = std::accumulate(locBasisVec.begin(),locBasisVec.begin()+j,0);

      DistrVector q_Big(MultiDomainDynam::solVecInfo());
      DistrVector qd_Big(MultiDomainDynam::solVecInfo());

      if(VtV.size() == 0){
        normalizedBasis_.expand(q,  q_Big);
        normalizedBasis_.expand(qd, qd_Big); // make sure velocities are projected

        Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > qMap(q.data(), q.size());

        //zero out generalized coordinates
        q.zero();
        qd.zero();

        // reduced coordinates have already been projected into embedding space with previous basis
        // now project tham back on to the new subspace with updated basis
        normalizedBasis_.localBasisIs(startCol,blockCols);
        reduceDisp(q_Big,  q);  // d_n has already been updated, regardless of rotational degrees of freedom
        reduceDisp(qd_Big, qd);
      } else {
        normalizedBasis_.localBasisIs(startCol,blockCols);
        projectLocalBases(localBasisId,j,q);
        projectLocalBases(localBasisId,j,qd);
      }
      setLocalReducedMesh(j);
      localBasisId = j;
    }
  }
#endif
}

void
DistrExplicitPodProjectionNonLinDynamicBase::projectLocalBases(int i, int j, DistrVector &q)
{
#ifdef USE_EIGEN3
  int ki = locBasisVec[i];
  int pi = std::accumulate(locBasisVec.begin(), locBasisVec.begin()+i, 0);
  Eigen::Map<Eigen::VectorXd> qi(q.data()+pi,ki);

  int kj = locBasisVec[j];
  int pj = std::accumulate(locBasisVec.begin(), locBasisVec.begin()+j, 0);
  Eigen::Map<Eigen::VectorXd> qj(q.data()+pj,kj);

  if(j < i) qj = VtV(j,i)*qi;
  else qj = VtV(i,j).transpose()*qi;

  qi.setZero();
#endif
}

//end local basis functions

void
DistrExplicitPodProjectionNonLinDynamicBase::reduceDisp(DistrVector &d, DistrVector &dr) const
{
  DistrVector Md(MultiDomainDynam::solVecInfo());
  dynMat->M->mult(d, Md);
  normalizedBasis_.reduce(Md, dr, false);
}

void 
DistrExplicitPodProjectionNonLinDynamicBase::updateState(double dt_n_h, DistrVector& v_n_h, DistrVector& d_n1) {
  //update geomState for Fint, but no need to update displacment vector from geometry 

  DistrVector temp1(solVecInfo());
  temp1 = dt_n_h*v_n_h; // temp1 is size of reduced coordinates

  normalizedBasis_.expand( temp1, *d_n); 

  geomState->update(*d_n, 2);

  normalizedBasis_.addLocalPart(temp1,d_n1);  //we save the increment vectors for postprocessing, add only local part

  if(haveRot) { // currently we only need to project the velocity up when there are rotation dofs
                // int the future, there may be other cases in which this is also necessary, e.g. viscoelastic materials
                // XXX this is a bit inconsistent, because the transformation to convected angular velocity
                // in GeomState::setVelocity is done with the velocity at t_{n+0.5} and the rotation vector at t_{n+1}
    normalizedBasis_.expand(v_n_h, *v_n);
    geomState->setVelocity(*v_n, 2);
  }
  setLocalBasis(d_n1,v_n_h); // this will update the local basis if multiple bases are supplied  
}

void DistrExplicitPodProjectionNonLinDynamicBase::getConstForce(DistrVector& v)
{
  int numNeumanModal = domain->nNeumannModal();
  if(numNeumanModal) {
    filePrint(stderr, " ... Using Reduced Constant Force   ...\n");
    BCond* nbcModal = domain->getNBCModal();
    v.zero();
    for(int i = 0; i < numNeumanModal; ++i) {
      if(nbcModal[i].nnum < v.size()) {
        if(!domain->getMFTT(nbcModal[i].loadsetid)) {
          double loadFactor = (nbcModal[i].loadsetid == -1) ? 1.0 : domain->getLoadFactor(nbcModal[i].loadsetid); // -1 is gravity
          v[nbcModal[i].nnum] += loadFactor*nbcModal[i].val;
        }
      }
    }
    fExt->zero();
  }
  else { // XXX currently, if modal forces are defined then any non-modal forces are ignored
    MultiDomainDynam::getConstForce(*cnst_fBig);
    normalizedBasis_.reduceAll(*cnst_fBig,v); // project using all local Bases
    cnst_fBig->zero();
  }
}

void
DistrExplicitPodProjectionNonLinDynamicBase::getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex)
{
  //Build internal force and project into reduced coordinates

  MultiDomainDynam::getInternalForce(*d_n, *fInt, t, tIndex);
  //compute residual here to prevent having to project into reduced basis twice
  *tempVec = *fInt - *fExt;
  f.zero();
  normalizedBasis_.reduce(*tempVec, f); 
}

void
DistrExplicitPodProjectionNonLinDynamicBase::computeExtForce2(SysState<DistrVector> &distState,
                        DistrVector &f, DistrVector &cnst_f, int tIndex,
                        double t, DistrVector *aero_f,
                        double gamma, double alphaf)
{
  int numNeumanModal = domain->nNeumannModal();
  if(numNeumanModal) {
    BCond* nbcModal = domain->getNBCModal();
    f.zero();
    normalizedBasis_.addLocalPart(cnst_f,f);
    if(domain->getNumMFTT() > 0) {
      for(int i = 0; i < numNeumanModal; ++i) {
        if(nbcModal[i].nnum < cnst_f.size()) {
          if(MFTTData *mftt = domain->getMFTT(nbcModal[i].loadsetid)) f[nbcModal[i].nnum] += mftt->getVal(t)*nbcModal[i].val;
        }
      }
    }
  }
  else { // currently, if modal forces are defined then any non-modal forces are ignored
    f.zero();
    normalizedBasis_.addLocalPart(cnst_f,f);
    MultiDomainDynam::computeExtForce2( *dummyState, *fExt, *cnst_fBig, tIndex, t, aero_fBig, gamma, alphaf);
  }
}

void
DistrExplicitPodProjectionNonLinDynamicBase::getContactForce(DistrVector &d, DistrVector &dinc, DistrVector &ctc_r, double t_n_p, double dt, double dt_old) {

  // first expand state increment into high dimensional container
  // could be optimized to expand only components on contact surfaces
  normalizedBasis_.expand(dinc, *tempVec);
  // next call masked function to compute contact force in embedding space
  MultiDomainDynam::getContactForce(*d_n, *tempVec, *ctc_f, t_n_p, dt, dt_old); 
  // then project contact force onto reduced subspace
  normalizedBasis_.sparseVecReduce(*ctc_f, ctc_r);

}

MDDynamMat *
DistrExplicitPodProjectionNonLinDynamicBase::buildOps(double mCoef, double cCoef, double kCoef) {
  if(!dynMat) MultiDomainDynam::buildOps(mCoef, cCoef, kCoef); // note: may be called previously in getInitState

  haveRot = geomState->getHaveRot();
  delete dynMat->dynMat;

  dynMat->dynMat = new DistrGalerkinProjectionSolver(normalizedBasis_);

  return dynMat;
}

void
DistrExplicitPodProjectionNonLinDynamicBase::computeStabilityTimeStep(double& dt, MDDynamMat& dynMat) {

  if(stableCount != 0) {
    GenMDDynamMat<double> ops;
    ops.K = dynMat.K;
    decDomain->rebuildOps(ops, 0, 0, 0, kelArray);
  }

  GenFullSquareMatrix<double> K_reduced;
  calculateReducedStiffness(*dynMat.K, normalizedBasis_, K_reduced);
  stableCount++;

  double dt_c;
  //computes stability timestep on root node then sends timestep to all other processes
  int myID = (structCom) ? structCom->myID() : 0;
  if(myID == 0) {
    dt_c = domain->computeStabilityTimeStepROM(K_reduced);
  }
  if(structCom)
    structCom->broadcast(1, &dt_c, 0);

  if(dt_c == std::numeric_limits<double>::infinity()) {
    filePrint(stderr," **************************************\n");
    filePrint(stderr," Stability max. timestep could not be  \n");
    filePrint(stderr," determined for this model.            \n");
    filePrint(stderr," Specified time step is selected\n");
    filePrint(stderr," **************************************\n");
    domain->solInfo().stable = 0;
  }
  else {
    filePrint(stderr," CONDITIONALLY STABLE NEWMARK ALGORITHM \n");
    filePrint(stderr," --------------------------------------\n");
    filePrint(stderr," Specified time step      = %10.4e\n", dt);
    filePrint(stderr," Stability max. time step = %10.4e\n", dt_c);
    filePrint(stderr," **************************************\n");
    if((domain->solInfo().stable == 1 && dt_c < dt) || domain->solInfo().stable == 2) {
      dt = dt_c;
      filePrint(stderr," Stability max. time step is selected\n");
    }
    else {
      filePrint(stderr," Specified time step is selected\n");
    }
    filePrint(stderr," **************************************\n");
  }
  
  for(int i = 0; i < decDomain->getNumSub(); ++i)
    decDomain->getSubDomain(i)->solInfo().setTimeStep(dt);
  domain->solInfo().setTimeStep(dt);
}

void
DistrExplicitPodProjectionNonLinDynamicBase::subInitStiff(int isub)
{
  SubDomain *sd = decDomain->getSubDomain(isub);

  Vector elementInternalForce(sd->maxNumDOF(), 0.0);
  Vector residual(sd->numUncon(), 0.0);
  sd->getStiffAndForce(*(*geomState)[isub], elementInternalForce, allCorot[isub], kelArray[isub], residual,
                       1.0, 0.0, (*geomState)[isub], (Vector*) NULL, ((melArray) ? melArray[isub] : NULL));
}

} // end namespace Rom
