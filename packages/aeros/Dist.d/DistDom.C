#include <cstdlib>
#include <Utils.d/dbg_alloca.h>

#include <Threads.d/PHelper.h>
#include <Feti.d/Feti.h>
#include <Timers.d/GetTime.h>
#include <Paral.d/DomainGroupTask.h>
#include <Utils.d/Memory.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/pstress.h>
#include <Driver.d/SubDomain.h>
#include <Driver.d/SysState.h>
#include <Paral.d/MDDynam.h>
#include <Driver.d/GeoSource.h>
#include <Paral.d/SubDOp.h>
#include <Mortar.d/MortarDriver.d/MortarHandler.h>
#include <Corotational.d/DistrGeomState.h>

extern FILE *debugFile;

//#define SERIALIZED_OUTPUT

template<class Scalar>
GenDistrDomain<Scalar>::GenDistrDomain(Domain *d) : GenDecDomain<Scalar>(d)
{
  initialize();
  this->myCPU = this->communicator->cpuNum();
}

template<class Scalar>
void 
GenDistrDomain<Scalar>::initialize()
{
  numRes = 0; 
  x = 0;
  masterFlag = 0; 
  numFlags = 0; 
  nodeOffsets = 0; 
  elemNodeOffsets = 0; 
  nodePat = 0;
  masterStress = 0;
  elemOffsets = 0;
  clusToCpu = 0;
}

template<class Scalar>
GenDistrDomain<Scalar>::~GenDistrDomain()
{
  if(masterFlag) { 
    for(int i=0; i<this->numSub; ++i)
      if(masterFlag[i]) { delete [] masterFlag[i]; masterFlag[i] = 0; }
    delete [] masterFlag; masterFlag = 0; 
  }
  if(numFlags) { delete [] numFlags; numFlags = 0; }
  if(nodeOffsets) { delete [] nodeOffsets; nodeOffsets = 0; }
  if(elemNodeOffsets) { delete [] elemNodeOffsets; elemNodeOffsets = 0; }
  if(numRes) { delete [] numRes; numRes = 0; }
  if(masterStress) { delete masterStress; masterStress = 0;} 
  if(nodePat) { delete nodePat; nodePat = 0; }
  if(elemOffsets) { delete [] elemOffsets; elemOffsets = 0; }
  if(clusToCpu) { delete clusToCpu; clusToCpu = 0; }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::clean()
{
  if(nodePat) { delete nodePat; nodePat = 0; }
  GenDecDomain<Scalar>::clean();
}

template<class Scalar>
void
GenDistrDomain<Scalar>::initPostPro()
{
  if(geoSource->getNumOutInfo()) {
#ifdef DISTRIBUTED
    createMasterFlag();
    createOutputOffsets();
    makeMasterInfo();
#endif
    if(!this->elemToNode && geoSource->elemOutput()) this->createElemToNode();
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::makeMasterInfo()
{
  masterInfo.domLen = new int[this->numSub];
  masterInfo.numDom = this->numSub;
  for(int iSub = 0; iSub < this->numSub; ++iSub) {
    masterInfo.domLen[iSub] = numFlags[iSub];
  }
#ifdef DISTRIBUTED
  masterInfo.computeOffsets();
#else
  masterInfo.setMasterFlag();
#endif
}

template<class Scalar>
void
GenDistrDomain<Scalar>::forceContinuity(GenDistrVector<Scalar> &u) {
  if(!masterFlag) initPostPro();

  int iSub;

  // initialize and merge displacements from subdomains into cpu array
  DistSVec<Scalar, 11> disps(*this->nodeInfo);
  DistSVec<Scalar, 11> masterDisps(masterInfo);
  disps = 0;
  for(iSub = 0; iSub < this->numSub; ++iSub) {
    Scalar (*xyz)[11] = (Scalar (*)[11]) disps.subData(iSub);
    Scalar *bcx = this->subDomain[iSub]->getBcx();
    this->subDomain[iSub]->template mergeDistributedDisp<Scalar>(xyz, u.subData(iSub), bcx);
  }
  unify(disps); // make sure master has solution for all dofs before reducing
  disps.reduce(masterDisps, masterFlag, numFlags);
  this->communicator->sync();
  for(iSub = 0; iSub < this->numSub; ++iSub) {
    Scalar (*xyz)[11] = (Scalar (*)[11]) disps.subData(iSub);
    this->subDomain[iSub]->template forceDistributedContinuity<Scalar>(u.subData(iSub), xyz);
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::postProcessing(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &f,
                                       double eigV, GenDistrVector<Scalar> *aeroF, int x,
                                       GenMDDynamMat<Scalar> *dynOps, SysState<GenDistrVector<Scalar> > *distState, int ndflag)
{
  int numOutInfo = geoSource->getNumOutInfo();
  if(numOutInfo == 0) return;

  // get output information
  OutputInfo *oinfo = geoSource->getOutputInfo();

  // check if there are any output files which need to be processed now
  int outLimit = geoSource->getOutLimit();
  if(geoSource->noOutput(x, ndflag) && !((x == domain->solInfo().initialTimeIndex) || (outLimit > 0 && x%outLimit == 0))) return;

  if(domain->outFlag && domain->nodeTable == 0) domain->makeNodeTable(domain->outFlag);
  int iOut_ffp = -1;
  int iOut_kir = -1;

  if(x == domain->solInfo().initialTimeIndex && ndflag == 0 && !(domain->solInfo().isDynam() || domain->solInfo().timeIntegration == 1))
    filePrint(stderr," ... Postprocessing                 ...\n");
  if(!masterFlag) initPostPro();

  int iSub;

  // initialize and merge displacements from subdomains into cpu array
  DistSVec<Scalar, 11> disps_glo(*this->nodeInfo);
  DistSVec<Scalar, 11> masterDisps_glo(masterInfo);
  disps_glo = 0;
  DistSVec<Scalar, 11> *disps_loc = 0, *masterDisps_loc = 0;
  if(!domain->solInfo().basicDofCoords) {
    disps_loc = new DistSVec<Scalar, 11>(*this->nodeInfo);
    masterDisps_loc = new DistSVec<Scalar, 11>(masterInfo);
  }
  for(iSub = 0; iSub < this->numSub; ++iSub) {
    Scalar (*xyz)[11] = (Scalar (*)[11]) disps_glo.subData(iSub);
    Scalar *bcx = this->subDomain[iSub]->getBcx();
    Scalar (*xyz_loc)[11] = (disps_loc) ? (Scalar (*)[11]) disps_loc->subData(iSub) : 0;
    this->subDomain[iSub]->template mergeDistributedDisp<Scalar>(xyz, u.subData(iSub), bcx, xyz_loc);
  }
  if((domain->solInfo().isCoupled && domain->solInfo().isMatching) || domain->GetnContactSurfacePairs()) {
    unify(disps_glo); // make sure master has solution for all dofs before reducing
    if(disps_loc) unify(*disps_loc);
  }
  disps_glo.reduce(masterDisps_glo, masterFlag, numFlags);
  if(disps_loc) disps_loc->reduce(*masterDisps_loc, masterFlag, numFlags);

  // initialize and merge aeroelastic forces
  DistSVec<Scalar, 6> aerof(*this->nodeInfo);
  DistSVec<Scalar, 6> masterAeroF(masterInfo);
  if(domain->solInfo().aeroFlag > -1 && aeroF) {
    GenDistrVector<Scalar> assembledAeroF(*aeroF);
    this->ba->assemble(assembledAeroF);
    aerof = 0;
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      Scalar (*mergedAeroF)[6] = (Scalar (*)[6]) aerof.subData(iSub);
      this->subDomain[iSub]->mergeDistributedForces(mergedAeroF, assembledAeroF.subData(iSub));
    }
    aerof.reduce(masterAeroF, masterFlag, numFlags);
  }

  // initialize and merge velocities & accelerations
  DistSVec<Scalar, 11> vels_glo(*this->nodeInfo), accs_glo(*this->nodeInfo);
  DistSVec<Scalar, 11> masterVels_glo(masterInfo), masterAccs_glo(masterInfo);
  DistSVec<Scalar, 11> *vels_loc = 0, *masterVels_loc = 0, *accs_loc = 0, *masterAccs_loc = 0;
  if(!domain->solInfo().basicDofCoords) {
    vels_loc = new DistSVec<Scalar, 11>(*this->nodeInfo);
    masterVels_loc = new DistSVec<Scalar, 11>(masterInfo);
    accs_loc = new DistSVec<Scalar, 11>(*this->nodeInfo);
    masterAccs_loc = new DistSVec<Scalar, 11>(masterInfo);
  }
  if(distState) {
    GenDistrVector<Scalar> *v_n = &distState->getVeloc();
    GenDistrVector<Scalar> *a_n = &distState->getAccel();
    vels_glo = 0; accs_glo = 0;
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      Scalar (*mergedVel)[11] = (Scalar (*)[11]) vels_glo.subData(iSub);
      Scalar (*mergedAcc)[11] = (Scalar (*)[11]) accs_glo.subData(iSub);
      double *vcx = this->subDomain[iSub]->getVcx();
      double *acx = this->subDomain[iSub]->getAcx();
      Scalar *vcx_scalar = new Scalar[this->subDomain[iSub]->numdof()];
      Scalar *acx_scalar = new Scalar[this->subDomain[iSub]->numdof()];
      Scalar (*mergedVel_loc)[11] = (vels_loc) ? (Scalar (*)[11]) vels_loc->subData(iSub) : 0;
      Scalar (*mergedAcc_loc)[11] = (accs_loc) ? (Scalar (*)[11]) accs_loc->subData(iSub) : 0;
      for(int i=0; i<this->subDomain[iSub]->numdof(); ++i) { vcx_scalar[i] = vcx[i]; acx_scalar[i] = acx[i]; }
      this->subDomain[iSub]->template mergeDistributedDisp<Scalar>(mergedVel, v_n->subData(iSub), vcx_scalar, mergedVel_loc);
      this->subDomain[iSub]->template mergeDistributedDisp<Scalar>(mergedAcc, a_n->subData(iSub), acx_scalar, mergedAcc_loc);
      delete [] vcx_scalar;
      delete [] acx_scalar;
    }
    vels_glo.reduce(masterVels_glo, masterFlag, numFlags);
    accs_glo.reduce(masterAccs_glo, masterFlag, numFlags);
    if(vels_loc) vels_loc->reduce(*masterVels_loc, masterFlag, numFlags);
    if(accs_loc) accs_loc->reduce(*masterAccs_loc, masterFlag, numFlags);
  }

  // compute current time (or frequency in the case of a helmholtz problem)
  double time;
  if(geoSource->isShifted() && domain->probType() != SolverInfo::Modal 
                            && domain->probType() != SolverInfo::Dynamic) {
    time = domain->getFrequencyOrWavenumber();
    if(domain->solInfo().doFreqSweep) x = this->outFreqCount++;
  } else if(domain->probType() == SolverInfo::Modal) {
    time = eigV;
    if(domain->solInfo().doEigSweep) x = this->outEigCount++;
  }
  else time = eigV;
  if (domain->solInfo().loadcases.size() > 0 && !domain->solInfo().doFreqSweep) time = domain->solInfo().loadcases.front();

// RT - serialize the OUTPUT,  PJSA - stress output doesn't work with serialized output. need to reconsider
#ifdef SERIALIZED_OUTPUT
for(int iCPU = 0; iCPU < this->communicator->size(); iCPU++) {
 this->communicator->sync();
 if(this->communicator->cpuNum() == iCPU) {
#endif

  // open binary output files
  if(x == domain->solInfo().initialTimeIndex) {
    if(!numRes) {
      numRes = new int[numOutInfo];
      for(int i=0; i<numOutInfo; ++i) numRes[i] = 0;
    }
  }

  if((x == domain->solInfo().initialTimeIndex) || (outLimit > 0 && x%outLimit == 0)) {
#ifdef DISTRIBUTED

    int *subToClus = geoSource->getSubToClus();
    int clusterId = (this->numSub > 0) ? subToClus[this->localSubToGl[0]] : 0;
    int firstCpuInCluster = (clusToCpu) ? (*clusToCpu)[clusterId][0] : 0;
    for(int iInfo = 0; iInfo < numOutInfo; iInfo++) {
      if(oinfo[iInfo].type == OutputInfo::Farfield || 
         oinfo[iInfo].type == OutputInfo::Energies ||
         oinfo[iInfo].type == OutputInfo::Kirchhoff || 
         oinfo[iInfo].type == OutputInfo::ModalDsp ||
         oinfo[iInfo].type == OutputInfo::ModalExF ||
         oinfo[iInfo].type == OutputInfo::ModalMass ||
         oinfo[iInfo].type == OutputInfo::ModalStiffness ||
         oinfo[iInfo].type == OutputInfo::ModalDamping ||
         oinfo[iInfo].type == OutputInfo::ModalDynamicMatrix ||
         oinfo[iInfo].type == OutputInfo::ModalMatrices ||
         oinfo[iInfo].type == OutputInfo::AeroForce) { 
        int oI = iInfo;
        if(this->firstOutput) { geoSource->openOutputFiles(0,&oI,1); } 
        continue;
      }
      else if(oinfo[iInfo].nodeNumber == -1 && this->firstOutput) {
        if(this->communicator->cpuNum() == firstCpuInCluster && this->numSub > 0) geoSource->createBinaryOutputFile(iInfo,this->localSubToGl[0],x);
        else geoSource->computeAndCacheHeaderLength(iInfo);
      }
    }
#ifndef SERIALIZED_OUTPUT
    this->communicator->sync();
#endif

    for(int iInfo = 0; iInfo < numOutInfo; iInfo++) {
      if(oinfo[iInfo].nodeNumber == -1 && 
         oinfo[iInfo].type != OutputInfo::Farfield && 
         oinfo[iInfo].type != OutputInfo::Energies &&
         oinfo[iInfo].type != OutputInfo::Kirchhoff && 
         oinfo[iInfo].type != OutputInfo::ModalExF &&
         oinfo[iInfo].type != OutputInfo::ModalMass &&
         oinfo[iInfo].type != OutputInfo::ModalStiffness &&
         oinfo[iInfo].type != OutputInfo::ModalDamping &&
         oinfo[iInfo].type != OutputInfo::ModalDynamicMatrix &&
         oinfo[iInfo].type != OutputInfo::ModalMatrices &&
         oinfo[iInfo].type != OutputInfo::AeroForce) {
        for(iSub = 0; iSub < this->numSub; iSub++) {
          int glSub = this->localSubToGl[iSub];
          if(oinfo[iInfo].dataType == 1) {
            geoSource->outputRange(iInfo, masterFlag[iSub],
                       this->subDomain[iSub]->numNode(), glSub, nodeOffsets[iSub], x);
          }
          else {
            geoSource->outputRange(iInfo, this->subDomain[iSub]->getGlElems(),
                       this->subDomain[iSub]->numElements(), glSub,
                       elemOffsets[iSub], x);
          }
        }
      }
    }


    if(x == domain->solInfo().initialTimeIndex) { // always put single node output in one file regardless of outLimit
      for(iSub = 0; iSub < this->numSub; iSub++)  {
        int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
        if(nOutNodes) {
          if(this->firstOutput) geoSource->openOutputFiles(this->subDomain[iSub]->getOutputNodes(), 
                   this->subDomain[iSub]->getOutIndex(), nOutNodes);                    
        }
      }
    }

#else
    if(x == domain->solInfo().initialTimeIndex && this->firstOutput) geoSource->openOutputFiles();  // opens all output files
#endif
  }

  if (!nodePat) makeNodePat();

  for(int iOut = 0; iOut < numOutInfo; iOut++) {

    if(oinfo[iOut].interval == 0 || x % oinfo[iOut].interval != 0)
      continue;

    // update number of results
    numRes[iOut]++; 

    if(oinfo[iOut].ndtype != ndflag) continue;
    if(ndflag !=0 && oinfo[iOut].type != OutputInfo::Disp6DOF && oinfo[iOut].type !=  OutputInfo::Displacement) continue;

    // set primal output states to either global or local
    DistSVec<Scalar, 11> &disps = (oinfo[iOut].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords)
                                  ? disps_glo : *disps_loc;
    DistSVec<Scalar, 11> &masterDisps = (oinfo[iOut].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords)
                                  ? masterDisps_glo : *masterDisps_loc;
    DistSVec<Scalar, 11> &vels = (oinfo[iOut].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords)
                                  ? vels_glo : *vels_loc;
    DistSVec<Scalar, 11> &masterVels = (oinfo[iOut].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords)
                                  ? masterVels_glo : *masterVels_loc;
    DistSVec<Scalar, 11> &accs = (oinfo[iOut].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords)
                                  ? accs_glo : *accs_loc;
    DistSVec<Scalar, 11> &masterAccs = (oinfo[iOut].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords)
                                  ? masterAccs_glo : *masterAccs_loc;

    switch(oinfo[iOut].type)  {

      case OutputInfo::EigenPair:
      case OutputInfo::FreqRespModes:
      case OutputInfo::Displacement:
        getPrimal(disps, masterDisps, time, x, iOut, 3, 0);
        break;
      case OutputInfo::Velocity:
        if(distState) getPrimal(vels, masterVels, time, x, iOut, 3, 0);
        break;
      case OutputInfo::Acceleration:
        if(distState) getPrimal(accs, masterAccs, time, x, iOut, 3, 0);
        break;
      case OutputInfo::EigenPair6:
      case OutputInfo::Disp6DOF:
        getPrimal(disps, masterDisps, time, x, iOut, 6, 0);
        break;
      case OutputInfo::Velocity6:
        if(distState) getPrimal(vels, masterVels, time, x, iOut, 6, 0);
        break;
      case OutputInfo::Accel6:
        if(distState) getPrimal(accs, masterAccs, time, x, iOut, 6, 0);
        break;
      case OutputInfo::Temperature:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 6);
        break;
      case OutputInfo::TemperatureFirstTimeDerivative:
        if(distState) getPrimal(vels, masterVels, time, x, iOut, 1, 6);
        break;
      case OutputInfo::AcousticPressure:
      case OutputInfo::EigenPressure:
      case OutputInfo::HelmholtzModes:
      case OutputInfo::Helmholtz:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 7);
        break;
      case OutputInfo::PressureFirstTimeDerivative:
        if(distState) getPrimal(vels, masterVels, time, x, iOut, 1, 7);
        break;
      case OutputInfo::PressureSecondTimeDerivative:
        if(distState) getPrimal(accs, masterAccs, time, x, iOut, 1, 7);
        break;
      case OutputInfo::StressXX:
        getStressStrain(u, time, x, iOut, SXX);
        break;
      case OutputInfo::StressYY:
        getStressStrain(u, time, x, iOut, SYY);
        break;
      case OutputInfo::StressZZ:
        getStressStrain(u, time, x, iOut, SZZ);
        break;
      case OutputInfo::StressXY:
        getStressStrain(u, time, x, iOut, SXY);
        break;
      case OutputInfo::StressYZ:
        getStressStrain(u, time, x, iOut, SYZ);
        break;
      case OutputInfo::StressXZ:
        getStressStrain(u, time, x, iOut, SXZ);
        break;
      case OutputInfo::StrainXX:
        getStressStrain(u, time, x, iOut, EXX);
        break;
      case OutputInfo::StrainYY:
        getStressStrain(u, time, x, iOut, EYY);
        break;
      case OutputInfo::StrainZZ:
        getStressStrain(u, time, x, iOut, EZZ);
        break;
      case OutputInfo::StrainXY:
        getStressStrain(u, time, x, iOut, EXY);
        break;
      case OutputInfo::StrainYZ:
        getStressStrain(u, time, x, iOut, EYZ);
        break;
      case OutputInfo::StrainXZ:
        getStressStrain(u, time, x, iOut, EXZ);
        break;
      case OutputInfo::StressVM:
        getStressStrain(u, time, x, iOut, VON);
        break;
      case OutputInfo::StrainVM:
        getStressStrain(u, time, x, iOut, STRAINVON);
        break;
      case OutputInfo::ContactPressure: {
        if(!domain->tdenforceFlag())
          getStressStrain(u, time, x, iOut, CONPRESS);
        else
          filePrint(stderr," *** WARNING: Output case %d not supported \n", iOut);
      } break;
      case OutputInfo::Damage:
        getStressStrain(u, time, x, iOut, DAMAGE);
        break;
      case OutputInfo::EquivalentPlasticStrain:
        getStressStrain(u, time, x, iOut, EQPLSTRN);
        break;
      case OutputInfo::StressPR1:
        getPrincipalStress(u, time, x, iOut, PSTRESS1);
        break;
      case OutputInfo::StressPR2:
        getPrincipalStress(u, time, x, iOut, PSTRESS2);
        break;
      case OutputInfo::StressPR3:
        getPrincipalStress(u, time, x, iOut, PSTRESS3);
        break;
      case OutputInfo::StrainPR1:
        getPrincipalStress(u, time, x, iOut, PSTRAIN1);
        break;
      case OutputInfo::StrainPR2:
        getPrincipalStress(u, time, x, iOut, PSTRAIN2);
        break;
      case OutputInfo::StrainPR3:
        getPrincipalStress(u, time, x, iOut, PSTRAIN3);
        break;
      case OutputInfo::InXForce:
        getElementForce(u, time, x, iOut, INX);
        break;
      case OutputInfo::InYForce:
        getElementForce(u, time, x, iOut, INY);
        break;
      case OutputInfo::InZForce:
        getElementForce(u, time, x, iOut, INZ);
        break;
      case OutputInfo::AXMoment:
        getElementForce(u, time, x, iOut, AXM);
        break;
      case OutputInfo::AYMoment:
        getElementForce(u, time, x, iOut, AYM);
        break;
      case OutputInfo::AZMoment:
        getElementForce(u, time, x, iOut, AZM);
        break;
      case OutputInfo::DispX:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 0);
        break;
      case OutputInfo::DispY:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 1);
        break;
      case OutputInfo::DispZ:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 2);
        break;
      case OutputInfo::RotX:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 3);
        break;
      case OutputInfo::RotY:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 4);
        break;
      case OutputInfo::RotZ:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 5);
        break;
      case OutputInfo::DispMod:
        if(oinfo[iOut].nodeNumber == -1) {
          for(iSub = 0; iSub < this->numSub; ++iSub) {
            int size = masterDisps.subSize(iSub);
            Scalar (*xyz)[11] = (Scalar (*)[11]) masterDisps.subData(iSub);
            Scalar *dispMod = new Scalar[size];
            for(int iNode=0; iNode<size; ++iNode) {
              dispMod[iNode] = ScalarTypes::sqrt(xyz[iNode][0]*xyz[iNode][0] +
                                                 xyz[iNode][1]*xyz[iNode][1] +
                                                 xyz[iNode][2]*xyz[iNode][2]);
            }
            geoSource->writeNodeScalarToFile(dispMod, size, this->localSubToGl[iSub], nodeOffsets[iSub],
                                             iOut, x, numRes[iOut], time, 1, masterFlag[iSub]);
            delete [] dispMod;
          }
        }
        else {
          for(iSub = 0; iSub < this->numSub; ++iSub) {
            int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
            if(nOutNodes) {
              int *outIndex = this->subDomain[iSub]->getOutIndex();
              for(int iNode = 0; iNode < nOutNodes; iNode++) {
                if(outIndex[iNode] == iOut) {
                  Scalar (*xyz)[11] = (Scalar (*)[11]) disps.subData(iSub);
                  int *outNodes = this->subDomain[iSub]->getOutputNodes();
                  Scalar dispMod = ScalarTypes::sqrt(xyz[outNodes[iNode]][0]*xyz[outNodes[iNode]][0] +
                                                     xyz[outNodes[iNode]][1]*xyz[outNodes[iNode]][1] +
                                                     xyz[outNodes[iNode]][2]*xyz[outNodes[iNode]][2]);
                  geoSource->outputNodeScalars(iOut, &dispMod, 1, time);
                }
              }
            }
          }
        }
        break;
      case OutputInfo::RotMod:
        if(oinfo[iOut].nodeNumber == -1) {
          for(iSub = 0; iSub < this->numSub; ++iSub) {
            int size = masterDisps.subSize(iSub);
            Scalar (*xyz)[11] = (Scalar (*)[11]) masterDisps.subData(iSub);
            Scalar *rotMod = new Scalar[size];
            for(int iNode=0; iNode<size; ++iNode) {
              rotMod[iNode] = ScalarTypes::sqrt(xyz[iNode][3]*xyz[iNode][3] +
                                                xyz[iNode][4]*xyz[iNode][4] +
                                                xyz[iNode][5]*xyz[iNode][5]);
            }
            geoSource->writeNodeScalarToFile(rotMod, size, this->localSubToGl[iSub], nodeOffsets[iSub],
                                             iOut, x, numRes[iOut], time, 1, masterFlag[iSub]);
            delete [] rotMod;
          }
        }
        else {
          for(iSub = 0; iSub < this->numSub; ++iSub) {
            int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
            if(nOutNodes) {
              int *outIndex = this->subDomain[iSub]->getOutIndex();
              for(int iNode = 0; iNode < nOutNodes; iNode++) {
                if(outIndex[iNode] == iOut) {
                  Scalar (*xyz)[11] = (Scalar (*)[11]) disps.subData(iSub);
                  int *outNodes = this->subDomain[iSub]->getOutputNodes();
                  Scalar rotMod = ScalarTypes::sqrt(xyz[outNodes[iNode]][3]*xyz[outNodes[iNode]][3] +
                                                    xyz[outNodes[iNode]][4]*xyz[outNodes[iNode]][4] +
                                                    xyz[outNodes[iNode]][5]*xyz[outNodes[iNode]][5]);
                  geoSource->outputNodeScalars(iOut, &rotMod, 1, time);
                }
              }
            }
          }
        }
        break;
      case OutputInfo::TotMod:
        if(oinfo[iOut].nodeNumber == -1) {
          for(iSub = 0; iSub < this->numSub; ++iSub) {
            int size = masterDisps.subSize(iSub);
            Scalar (*xyz)[11] = (Scalar (*)[11]) masterDisps.subData(iSub);
            Scalar *totMod = new Scalar[size];
            for(int iNode=0; iNode<size; ++iNode) {
              totMod[iNode] = ScalarTypes::sqrt(xyz[iNode][0]*xyz[iNode][0] +
                                                xyz[iNode][1]*xyz[iNode][1] +
                                                xyz[iNode][2]*xyz[iNode][2] +
                                                xyz[iNode][3]*xyz[iNode][3] +
                                                xyz[iNode][4]*xyz[iNode][4] +
                                                xyz[iNode][5]*xyz[iNode][5]);
            }
            geoSource->writeNodeScalarToFile(totMod, size, this->localSubToGl[iSub], nodeOffsets[iSub],
                                             iOut, x, numRes[iOut], time, 1, masterFlag[iSub]);
            delete [] totMod;
          }
        }
        else {
          for(iSub = 0; iSub < this->numSub; ++iSub) {
            int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
            if(nOutNodes) {
              int *outIndex = this->subDomain[iSub]->getOutIndex();
              for(int iNode = 0; iNode < nOutNodes; iNode++) {
                if(outIndex[iNode] == iOut) {
                  Scalar (*xyz)[11] = (Scalar (*)[11]) disps.subData(iSub);
                  int *outNodes = this->subDomain[iSub]->getOutputNodes();
                  Scalar totMod = ScalarTypes::sqrt(xyz[outNodes[iNode]][0]*xyz[outNodes[iNode]][0] +
                                                    xyz[outNodes[iNode]][1]*xyz[outNodes[iNode]][1] +
                                                    xyz[outNodes[iNode]][2]*xyz[outNodes[iNode]][2] +
                                                    xyz[outNodes[iNode]][3]*xyz[outNodes[iNode]][3] +
                                                    xyz[outNodes[iNode]][4]*xyz[outNodes[iNode]][4] +
                                                    xyz[outNodes[iNode]][5]*xyz[outNodes[iNode]][5]);
                  geoSource->outputNodeScalars(iOut, &totMod, 1, time);
                }
              }
            }
          }
        }
        break;
      case OutputInfo::Energies:
        this->getEnergies(u, f, iOut, time, distState, dynOps, aeroF);
        break;
      case OutputInfo::Farfield: 
        domain->nffp = oinfo[iOut].interval;
        iOut_ffp = iOut; // PJSA 3-1-2007 buildFFP doesn't work with serialized output
        //this->buildFFP(u,oinfo[iOut].filptr,true);
        break;
      case OutputInfo::Kirchhoff:
        domain->nffp = oinfo[iOut].interval;
        iOut_kir = iOut; // PJSA 3-1-2007 buildFFP doesn't work with serialized output
        //this->buildFFP(u,oinfo[iOut].filptr,false);
        break;
      case OutputInfo::AeroForce: break; // this is done in DistFlExchange.C
      case OutputInfo::AeroXForce:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 0);
        break;
      case OutputInfo::AeroYForce:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 1);
        break;
      case OutputInfo::AeroZForce:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 2);
        break;
      case OutputInfo::AeroXMom:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 3);
        break;
      case OutputInfo::AeroYMom:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 4);
        break;
      case OutputInfo::AeroZMom:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 5);
        break;
/* TODO
      case OutputInfo::Reactions:
        break;
      case OutputInfo::Reactions6:
        break;
*/
      case OutputInfo::EigenSlosh:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 10);
        break;
      case OutputInfo::YModulus:
        this->getElementAttr(iOut,YOUNG, time);
        break;
      case OutputInfo::MDensity:
        this->getElementAttr(iOut,MDENS, time);
        break;
      case OutputInfo::Thicknes:
        this->getElementAttr(iOut,THICK, time);
        break;
      case OutputInfo::TDEnforcement: {
        if(domain->tdenforceFlag()) {
          DistSVec<double, 1> all_data(*this->nodeInfo);
          if(oinfo[iOut].tdenforc_var == 1) all_data = 0.5;
          else all_data = 0;
          double **sub_data = new double * [this->numSub];
          for(iSub = 0; iSub < this->numSub; ++iSub) sub_data[iSub] = (double *) all_data.subData(iSub);
          for(int iMortar=0; iMortar<domain->GetnMortarConds(); iMortar++) {
            domain->GetMortarCond(iMortar)->get_plot_variable(oinfo[iOut].tdenforc_var, sub_data, this->numSub, this->subDomain);
          }
          DistSVec<double, 1> master_data(masterInfo);
          all_data.reduce(master_data, masterFlag, numFlags);
          for(iSub = 0; iSub < this->numSub; ++iSub) {
            geoSource->writeNodeScalarToFile((double *) master_data.subData(iSub), master_data.subSize(iSub), this->localSubToGl[iSub], nodeOffsets[iSub],
                                             iOut, x, numRes[iOut], time, 1, masterFlag[iSub]);
          }
          delete [] sub_data;
        }
        else filePrint(stderr," *** WARNING: Output case %d not supported \n", iOut); 
      } break;
      case OutputInfo::Statevector:
      case OutputInfo::Velocvector:
      case OutputInfo::Accelvector:
      case OutputInfo::InternalStateVar:
      case OutputInfo::DualStateVar:
      case OutputInfo::MuStateVar:
      case OutputInfo::Forcevector:
      case OutputInfo::Constraintvector:
      case OutputInfo::Constraintviolation:
      case OutputInfo::Residual:
      case OutputInfo::Jacobian:
      case OutputInfo::RobData:
      case OutputInfo::SampleMesh:
      case OutputInfo::ModalDsp:
      case OutputInfo::ModalExF:
      case OutputInfo::ModalMass:
      case OutputInfo::ModalStiffness:
      case OutputInfo::ModalDamping:
      case OutputInfo::ModalDynamicMatrix:
      case OutputInfo::ModalMatrices:
        break;
      default:
        filePrint(stderr," *** WARNING: Output case %d not implemented \n", iOut);
        break;
    }
  }
  this->firstOutput = false;
// RT - serialize the OUTPUT
#ifdef SERIALIZED_OUTPUT
  }
}
this->communicator->sync();
#endif
  if(iOut_ffp > -1) this->buildFFP(u,oinfo[iOut_ffp].filptr,true); // PJSA 3-1-2007 buildFFP doesn't work with serialized output
  if(iOut_kir > -1) this->buildFFP(u,oinfo[iOut_kir].filptr,false); // PJSA 3-1-2007 buildFFP doesn't work with serialized output

  if(disps_loc) delete disps_loc;
  if(masterDisps_loc) delete masterDisps_loc;
  if(vels_loc) delete vels_loc;
  if(masterVels_loc) delete masterVels_loc;
  if(accs_loc) delete accs_loc;
  if(masterAccs_loc) delete masterAccs_loc;
}

template<class Scalar>
template<int dim>
void
GenDistrDomain<Scalar>::getPrimal(DistSVec<Scalar, dim> &disps, DistSVec<Scalar, dim> &masterDisps, 
                                  double time, int x, int fileNumber, int ndof, int startdof)
{
  // this function outputs the primal variables: displacement, temperature, pressure
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  int iSub;

  if(oinfo.nodeNumber == -1) { // output binary data for all nodes
    for(iSub = 0; iSub < this->numSub; ++iSub)   
      geoSource->writeNodeVectorToFile(masterDisps(iSub), this->localSubToGl[iSub],
                                       nodeOffsets[iSub], fileNumber, x, numRes[fileNumber],
                                       time, ndof, startdof, masterFlag[iSub]); 
  }
  else {  
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
      if(nOutNodes) {
        int *outIndex = this->subDomain[iSub]->getOutIndex();
        for(int iNode = 0; iNode < nOutNodes; iNode++) {
          if(outIndex[iNode] == fileNumber) {
            Scalar (*nodeDisp)[dim] = (Scalar (*)[dim]) disps.subData(iSub);
            int *outNodes = this->subDomain[iSub]->getOutputNodes();
            switch(ndof) {
              case 6:
                geoSource->outputNodeVectors6(fileNumber, nodeDisp + outNodes[iNode], 1, time);
                break;
              case 3:
                geoSource->outputNodeVectors(fileNumber, nodeDisp + outNodes[iNode], 1, time);
                break;
              case 1:
                geoSource->outputNodeScalars(fileNumber, nodeDisp[outNodes[iNode]]+startdof, 1, time);
                break;
              default:
                filePrint(stderr, " *** WARNING: single node primal output not supported for ndof = %d \n", ndof);
                break;
            }
          }
        }
      }
    }
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::getAeroForceScalar(DistSVec<Scalar, 6> &aerof, DistSVec<Scalar, 6> &masterAeroF,
                                           double time, int x, int fileNumber, int dof)
{
  // this function outputs an aeroelastic force scalar  
  // dof == 0 is X Force, dof == 1 is Y Force,  etc.
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  int iSub;

  if(oinfo.nodeNumber == -1) { // output binary data for all nodes
    for(iSub = 0; iSub < this->numSub; ++iSub)
      geoSource->writeNodeVectorToFile(masterAeroF(iSub), this->localSubToGl[iSub],
                                       nodeOffsets[iSub], fileNumber, x, numRes[fileNumber],
                                       time, 1, dof, masterFlag[iSub]);
  }
  else { // output ascii data at selected node or nodes only
    for(iSub = 0; iSub < this->numSub; ++iSub)  {
      int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
      if(nOutNodes)  {
        int *outIndex = this->subDomain[iSub]->getOutIndex();
        for(int iNode = 0; iNode < nOutNodes; iNode++)
          if(outIndex[iNode] == fileNumber)  {
            Scalar (*nodeAeroF)[6] = (Scalar (*)[6]) aerof.subData(iSub);
            int *outNodes = this->subDomain[iSub]->getOutputNodes();
            geoSource->outputNodeScalars(fileNumber, nodeAeroF[outNodes[iNode]]+dof, 1, time);
          }
      }
    }
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::setsizeSfemStress(int fileNumber)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;

  if(avgnum == 1) this->sizeSfemStress = masterInfo.totLen(); 
  else if(avgnum == 0) {  // element-based output
    this->sizeSfemStress = 0;
  }
  else {
    std::cerr << "avgnum = " << avgnum << " not implemented in Domain::setsizeSfemStress()" << std::endl;
    this->sizeSfemStress = 0;
  }
}

template<class Scalar>
Scalar*
GenDistrDomain<Scalar>::getSfemStress(int fileNumber)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;

  if(avgnum == 1) return masterStress->data(); // node-based
  else if(avgnum == 0) return 0; // element-based
  else  { std::cerr << "avgnum = " << avgnum << " not implemented in Domain::getSfemStress()" << std::endl; return 0; }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::updateSfemStress(Scalar* str, int fileNumber)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;

  if(avgnum == 1)  masterStress->setNewData(str);
  else if(avgnum == 0) std::cerr << "updateSfemStress for element not yet implemented" << std::endl;
  else std::cerr << "avgnum = " << avgnum << " not implemented in Domain::updateSfemStress()" << std::endl;
}

template<class Scalar>
void 
GenDistrDomain<Scalar>::getStressStrain(GenDistrVector<Scalar> &u, double time,
                                        int x, int fileNumber, int Findex, int printFlag) 
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  if(oinfo.averageFlg == 0) {
    getElementStressStrain(u, time, x, fileNumber, Findex, printFlag);
    return;
  }

  DistVec<Scalar> stress(*this->nodeInfo);
  DistVec<Scalar> weight(*this->nodeInfo);

  stress = 0;
  weight = 0;


  int iSub;
  
  if(printFlag != 2) {
    // each subdomain computes its stress vector
    for (iSub = 0; iSub < this->numSub; ++iSub) {
      if(Findex != 16) {
        this->subDomain[iSub]->computeStressStrain(fileNumber, u.subData(iSub), 
                                                   Findex, stress.subData(iSub), weight.subData(iSub));
      }
      else {
        this->subDomain[iSub]->computeLocalContactPressure(stress.subData(iSub), weight.subData(iSub));
      }
    }

    paralApply(this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>, nodePat, &stress);
    nodePat->exchange();
    paralApply(this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, nodePat, &stress);

    paralApply(this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>, nodePat, &weight);
    nodePat->exchange();
    paralApply(this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, nodePat, &weight);

    // Divide stress by weight
    for (iSub = 0; iSub < this->numSub; ++iSub)  {
      Vec<Scalar> &locStress = stress(iSub);
      Vec<Scalar> &locWeight = weight(iSub);
      for(int i = 0; i < stress.subSize(iSub); ++i)
        if(locWeight[i] != 0.0)
          locStress[i] /= locWeight[i];
        else
          locStress[i] = 0.0;
    }

    // reduce the stress vector to just the master quantities
    if(masterStress == 0) masterStress = new DistVec<Scalar>(masterInfo);
    stress.reduce(*masterStress, masterFlag, numFlags);
  }
  

  if(printFlag != 1) {
    if(oinfo.nodeNumber == -1) { // output binary or ascii data for all nodes or node group
      for(iSub = 0; iSub < this->numSub; ++iSub) {
        geoSource->writeNodeScalarToFile(masterStress->subData(iSub), masterStress->subSize(iSub), this->localSubToGl[iSub],
                                         nodeOffsets[iSub], fileNumber, x, numRes[fileNumber], time, 1, masterFlag[iSub]);
      }
    }
    else { // output ascii data for one node
      for(iSub = 0; iSub < this->numSub; ++iSub) {
        int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
        if(nOutNodes) {
          int *outIndex = this->subDomain[iSub]->getOutIndex();
          for(int iNode = 0; iNode < nOutNodes; iNode++) {
            if(outIndex[iNode] == fileNumber) {
              Scalar *nodeStress = (Scalar *) stress.subData(iSub);
              int *outNodes = this->subDomain[iSub]->getOutputNodes();
              geoSource->outputNodeScalars(fileNumber, nodeStress+outNodes[iNode], 1, time);
            }
          }
        }
      }
    }
  }

}

template<class Scalar>
void
GenDistrDomain<Scalar>::getElementStressStrain(GenDistrVector<Scalar> &u, double time, 
                                               int iter, int fileNumber, int Findex, int printFlag) 
{
  // allocate arrays
  for(int iSub = 0; iSub < this->numSub; iSub++)  {
    int numElemNodes = this->subDomain[iSub]->countElemNodes();
    Scalar *elemStress = new Scalar[numElemNodes];
    this->subDomain[iSub]->computeStressStrain(fileNumber, u.subData(iSub), Findex, elemStress);
    geoSource->writeElemScalarToFile(elemStress, numElemNodes, this->localSubToGl[iSub], elemNodeOffsets[iSub], fileNumber, iter,
                                     numRes[fileNumber], time, this->elemToNode->numConnect(), this->subDomain[iSub]->getGlElems());
    delete [] elemStress;
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::getElementPrincipalStress(GenDistrVector<Scalar> &u, double time,
                                                  int iter, int fileNumber, int strIndex)
{
  // set stress VS. strain for element subroutines
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];
                                                                                                                                                             
  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i)
      strDir[i] = i;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i)
      strDir[i] = i+7;
  }
  else {
    fprintf(stderr," *** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // allocate arrays
  for(int iSub = 0; iSub < this->numSub; iSub++) {  // this could be parallelized
    int numElemNodes = this->subDomain[iSub]->countElemNodes();
    Scalar **elemAllStress = new Scalar * [6];
    for(i=0; i<6; ++i) elemAllStress[i] = new Scalar[numElemNodes];
                                                                                                                               
    // Compute Each Required Stress (all 6) using same routines as for
    // individual stresses
    int str_loop;
    int Findex;
    for(str_loop = 0; str_loop < 6; ++str_loop) {
      // get current stress/strain index
      Findex = strDir[str_loop];

      // compute stress vector
      this->subDomain[iSub]->computeStressStrain(fileNumber, u.subData(iSub),
                                           Findex, elemAllStress[str_loop]);
    }

    // ... CALCULATE PRINCIPALS AT EACH NODE
    Scalar svec[6], pvec[3];
    Scalar *elemPVec = new Scalar[numElemNodes];
    for(i = 0; i < numElemNodes; ++i) {
      for(j = 0; j < 6; ++j)
        svec[j] = elemAllStress[j][i];
      // Convert Engineering to Tensor Strains
      if(stressORstrain != 0) {
        svec[3] /= 2;
        svec[4] /= 2;
        svec[5] /= 2;
      }
      pstress(svec, pvec);
      elemPVec[i] = pvec[strInd];
    }

    geoSource->writeElemScalarToFile(elemPVec, numElemNodes, this->localSubToGl[iSub], elemNodeOffsets[iSub], fileNumber, iter,
                                     numRes[fileNumber], time, this->elemToNode->numConnect(), this->subDomain[iSub]->getGlElems());

    delete [] elemPVec;
    for(i=0; i<6; ++i) delete [] elemAllStress[i];
    delete [] elemAllStress;
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::getElementPrincipalStress(DistrGeomState *gs, Corotator ***allCorot, double time,
                                                  int x, int fileNumber, int strIndex, DistrGeomState *refState)
{
  // set stress VS. strain for element subroutines (Non-linear version)
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];

  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i)
      strDir[i] = i;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i)
      strDir[i] = i+7;
  }
  else {
    fprintf(stderr," *** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // allocate arrays
  for(int iSub = 0; iSub < this->numSub; iSub++) {  // this could be parallelized
    int numElemNodes = this->subDomain[iSub]->countElemNodes();
    Scalar **elemAllStress = new Scalar * [6];
    for(i=0; i<6; ++i) elemAllStress[i] = new Scalar[numElemNodes];
                                                                                                                                                           
    // Compute Each Required Stress (all 6) using same routines as for
    // individual stresses
    int str_loop;
    int Findex;
    for(str_loop = 0; str_loop < 6; ++str_loop) {
      // get current stress/strain index
      Findex = strDir[str_loop];
                                                                                                                                                           
      // compute stress vector
      GeomState *subRefState = (refState) ? (*refState)[iSub] : 0;
      this->subDomain[iSub]->computeStressStrain((*gs)[iSub], allCorot[iSub], fileNumber,
                                           Findex, elemAllStress[str_loop], (Scalar *) 0, subRefState);
    }
    // ... CALCULATE PRINCIPALS AT EACH NODE
    Scalar svec[6], pvec[3];
    Scalar *elemPVec = new Scalar[numElemNodes];
    for(i = 0; i < numElemNodes; ++i) {
      for(j = 0; j < 6; ++j)
        svec[j] = elemAllStress[j][i];
      // Convert Engineering to Tensor Strains
      if(stressORstrain != 0) {
        svec[3] /= 2;
        svec[4] /= 2;
        svec[5] /= 2;
      }
      pstress(svec, pvec);
      elemPVec[i] = pvec[strInd];
    }

    geoSource->writeElemScalarToFile(elemPVec, numElemNodes, this->localSubToGl[iSub], elemNodeOffsets[iSub], fileNumber, x,
                                     numRes[fileNumber], time, this->elemToNode->numConnect(), this->subDomain[iSub]->getGlElems());

    delete [] elemPVec;
    for(i=0; i<6; ++i) delete [] elemAllStress[i];
    delete [] elemAllStress;
  }
}

template<class Scalar>
void 
GenDistrDomain<Scalar>::getPrincipalStress(GenDistrVector<Scalar> &u, double time, int x, 
                                           int fileNumber, int strIndex)  
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  if(oinfo.averageFlg == 0) {
    getElementPrincipalStress(u, time, x, fileNumber, strIndex);
    return;
  }

  // set stress VS. strain for element subroutines
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];

  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i)
      strDir[i] = i;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i)
      strDir[i] = i+7;
  }
  else {
    fprintf(stderr," *** ERROR: Bad Principal Stress Direction\n");
    exit(-1);
  }

  // Allocate a distributed vector for stress
  DistVec<Scalar> **stress = new DistVec<Scalar>*[6];
  DistVec<Scalar> weight(*this->nodeInfo);

  int iSub;
  int str_loop;
  for(str_loop = 0; str_loop < 6; ++str_loop)
    stress[str_loop] = new DistVec<Scalar> (*this->nodeInfo);  

  // each subdomain computes its stress/strain vector

  // Compute Each Required Stress (all 6) using same routines as for
  // individual stresses
  int Findex;

  for(str_loop = 0; str_loop < 6; ++str_loop) {

    // get current stress/strain index
    Findex = strDir[str_loop];

    // Initialize distributed vector to zero
    *stress[str_loop] = 0;
    weight = 0;

    for (iSub = 0; iSub < this->numSub; ++iSub) {
      this->subDomain[iSub]->computeStressStrain(fileNumber, u.subData(iSub),
                Findex, stress[str_loop]->subData(iSub), weight.subData(iSub));
    }

    paralApply(this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>,
               this->nodePat, stress[str_loop]);
    this->nodePat->exchange();
    paralApply(this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, this->nodePat, stress[str_loop]);

    paralApply(this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>, this->nodePat, &weight);
    this->nodePat->exchange();
    paralApply(this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, this->nodePat, &weight);

    // Divide stress by weight
    for(iSub = 0; iSub < this->numSub; ++iSub)  {
      Vec<Scalar> &locStress = (*stress[str_loop])(iSub);
      Vec<Scalar> &locWeight = weight(iSub);
      for(int i = 0; i < stress[str_loop]->subSize(iSub); ++i)  {
        if (locWeight[i] != 0.0)
          locStress[i] /= locWeight[i];
        else
          locStress[i] = 0.0;
      }
    }
  }

  // Calculate Principals at each node
  Scalar svec[6], pvec[3];
  DistVec<Scalar> allPVec(*this->nodeInfo);
  for(iSub = 0; iSub < this->numSub; ++iSub)  {

    Vec<Scalar> &locPVec = allPVec(iSub);
    for(i = 0; i < this->subDomain[iSub]->numNode(); ++i) {

      for(j = 0; j < 6; ++j)
        svec[j] = stress[j]->subData(iSub)[i];
 
      // Convert Engineering to Tensor Strains
      if(stressORstrain != 0) {
        svec[3] /= 2;
        svec[4] /= 2;
        svec[5] /= 2;
      }
      pstress(svec,pvec);
      locPVec[i] = pvec[strInd];
    }
  }

  // reduce stress vector to master quantities
  DistVec<Scalar> masterPVec(masterInfo);
  allPVec.reduce(masterPVec, masterFlag, numFlags);

  if(oinfo.nodeNumber == -1) { // output binary or ascii data for all nodes or node group
    for(iSub = 0; iSub < this->numSub; ++iSub)  {
      geoSource->writeNodeScalarToFile(masterPVec.subData(iSub), masterPVec.subSize(iSub), this->localSubToGl[iSub],
                                       nodeOffsets[iSub], fileNumber, x, numRes[fileNumber], time, 1, masterFlag[iSub]); 
    }
  }
  else { // output ascii data for one node
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
      if(nOutNodes) {
        int *outIndex = this->subDomain[iSub]->getOutIndex();
        for(int iNode = 0; iNode < nOutNodes; iNode++) {
          if(outIndex[iNode] == fileNumber) {
            Scalar *nodeStress = (Scalar *) allPVec.subData(iSub);
            int *outNodes = this->subDomain[iSub]->getOutputNodes();
            geoSource->outputNodeScalars(fileNumber, nodeStress+outNodes[iNode], 1, time);
          }
        }
      }
    }
  }
}

template<class Scalar>
void 
GenDistrDomain<Scalar>::getElementForce(GenDistrVector<Scalar> &u, double time, int x, 
                                        int fileNumber, int Findex)  
{
  for(int iSub = 0; iSub < this->numSub; iSub++)  {
    int numElemNodes = this->subDomain[iSub]->countElemNodes();
    Scalar *elemForce = new Scalar[numElemNodes];
    this->subDomain[iSub]->computeElementForce(fileNumber, u.subData(iSub),
                                               Findex, elemForce);
    geoSource->writeElemScalarToFile(elemForce, numElemNodes, this->localSubToGl[iSub], elemNodeOffsets[iSub], fileNumber, x,
                                     numRes[fileNumber], time, this->elemToNode->numConnect(), this->subDomain[iSub]->getGlElems());
    delete [] elemForce;
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::getElementForce(DistrGeomState *gs, Corotator ***allCorot, double time, int x,
                                        int fileNumber, int Findex)
{
  for(int iSub = 0; iSub < this->numSub; iSub++)  {
    int numElemNodes = this->subDomain[iSub]->countElemNodes();
    Scalar *elemForce = new Scalar[numElemNodes];
    this->subDomain[iSub]->computeElementForce((*gs)[iSub], allCorot[iSub],
                                               fileNumber, Findex, elemForce);
    geoSource->writeElemScalarToFile(elemForce, numElemNodes, this->localSubToGl[iSub], elemNodeOffsets[iSub], fileNumber, x,
                                     numRes[fileNumber], time, this->elemToNode->numConnect(), this->subDomain[iSub]->getGlElems());
    delete [] elemForce;
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::createMasterFlag() 
{
  // allocate size of master flag
  masterFlag = new int *[this->numSub];
  numFlags = new int[this->numSub];

  int iSub;

  for(iSub = 0; iSub < this->numSub; iSub++)  {
    int numNodes = this->subDomain[iSub]->numNode();
    masterFlag[iSub] = new int[numNodes];
    numFlags[iSub] = numNodes;

    // initialize master flag for each subdomain
    auto &clusNodes = this->subDomain[iSub]->getGlNodes();
    for (int iNode = 0; iNode < numNodes; iNode++)
      masterFlag[iSub][iNode] = clusNodes[iNode];
  }

  int *subToClus = geoSource->getSubToClus();
  for(iSub = 0; iSub < this->numSub; iSub++) {
    // get connected subdomain information
    auto &connDoms = this->subDomain[iSub]->getSComm()->subNums;

    int thisGlSub = this->subDomain[iSub]->subNum();
    int thisCluster = subToClus ? subToClus[thisGlSub] : 0;

    // loop over connected subdomains in this cluster
    for(int jSub = 0; jSub < connDoms.size(); jSub++) {
      // set master flag if connected subdomain global number is smaller
      if(connDoms[jSub] < thisGlSub && (subToClus == 0 || subToClus[connDoms[jSub]] == thisCluster)) {

        // get shared nodes
        auto &sharedNodes = this->subDomain[iSub]->getSComm()->sharedNodes;
        for(int iNode = 0; iNode < sharedNodes->num(jSub); iNode++) {
          int shNode = (*sharedNodes)[jSub][iNode];
          if(masterFlag[iSub][shNode] >= 0) {
            masterFlag[iSub][shNode] = -1;
            numFlags[iSub]--;
          }
        }
      }
    }
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::createOutputOffsets() 
{
  const Connectivity &clusToSub = geoSource->getClusToSub();
  int *subToClus = geoSource->getSubToClus();
  nodeOffsets = new int[this->numSub];
  elemNodeOffsets = new int[this->numSub];
  elemOffsets = new int[this->numSub];

  int iSub, iCluster;
  int numClusters = geoSource->getNumClusters();
  if(numClusters == 0) return;
  // create cluster array of offsets
  int **clNodeOffsets = new int *[numClusters];
  int **clElemNodeOffsets = new int *[numClusters];
  int **clElemOffsets = new int *[numClusters]; 

  for(iCluster = 0; iCluster < numClusters; iCluster++) {
    clNodeOffsets[iCluster] = new int[clusToSub.num(iCluster)];
    clElemNodeOffsets[iCluster] = new int[clusToSub.num(iCluster)];
    clElemOffsets[iCluster] = new int[clusToSub.num(iCluster)];
    for(int j=0; j<clusToSub.num(iCluster); ++j) { // initialize to zero so global sum will work
      clNodeOffsets[iCluster][j] = 0;
      clElemNodeOffsets[iCluster][j] = 0;
      clElemOffsets[iCluster][j] = 0;
    }
  }

  // populate glOffsets w/this mpi's data
  for(iSub = 0; iSub < this->numSub; iSub++) {
    int glSub = this->subDomain[iSub]->subNum();
    int clusNum = subToClus[glSub];
    for(int jSub = 0; jSub < clusToSub.num(clusNum); jSub++)
      if(glSub == clusToSub[clusNum][jSub]) {
        clNodeOffsets[clusNum][jSub] = numFlags[iSub];
        clElemNodeOffsets[clusNum][jSub] = this->subDomain[iSub]->countElemNodes();
        clElemOffsets[clusNum][jSub] = this->subDomain[iSub]->numElements();
      }
  }

  // sum up all offsets in all mpi processes
  for(iCluster = 0; iCluster < numClusters; iCluster++) {
    this->communicator->globalSum(clusToSub.num(iCluster), clNodeOffsets[iCluster]);
    this->communicator->globalSum(clusToSub.num(iCluster), clElemNodeOffsets[iCluster]);
    this->communicator->globalSum(clusToSub.num(iCluster), clElemOffsets[iCluster]);
  }

  for(iCluster = 0; iCluster < numClusters; iCluster++) {
    int nOffset = 0;
    int enOffset = 0;
    int eOffset = 0;
    for(iSub = 0; iSub < clusToSub.num(iCluster); iSub++) {
      int tmpNOff = clNodeOffsets[iCluster][iSub];
      int tmpENOff = clElemNodeOffsets[iCluster][iSub];
      int tmpEOff = clElemOffsets[iCluster][iSub];
      clNodeOffsets[iCluster][iSub] = nOffset;
      clElemNodeOffsets[iCluster][iSub] = enOffset;
      clElemOffsets[iCluster][iSub] = eOffset;
      nOffset += tmpNOff;
      enOffset += tmpENOff;
      eOffset += tmpEOff;
    }
  }
   
  for(iSub = 0; iSub < this->numSub; iSub++) {
    int glSub = this->subDomain[iSub]->subNum();
    int clusNum = subToClus[glSub];
    for(int jSub = 0; jSub < clusToSub.num(clusNum); jSub++)
      if(glSub == clusToSub[clusNum][jSub]) {
        nodeOffsets[iSub] = clNodeOffsets[clusNum][jSub];
        elemNodeOffsets[iSub] = clElemNodeOffsets[clusNum][jSub];
        elemOffsets[iSub] = clElemOffsets[clusNum][jSub];
      }
  }

  if(clNodeOffsets) { for(int i=0; i<numClusters; i++) delete [] clNodeOffsets[i]; delete [] clNodeOffsets; }
  if(clElemNodeOffsets) { for(int i=0; i<numClusters; i++) delete [] clElemNodeOffsets[i]; delete [] clElemNodeOffsets; }
  if(clElemOffsets) { for(int i=0; i<numClusters; i++) delete [] clElemOffsets[i]; delete [] clElemOffsets; }

  // create clusToCpu connectivity, used to decide which process should initially open each output file
  Connectivity subToCpu = this->cpuToSub->reverse();
  clusToCpu = clusToSub.transcon(&subToCpu);
}

// ---------------------------------------------------------------
// Nonlinear DistrDomain functions
// ----------------------------------------------------------------

template<class Scalar>
void
GenDistrDomain<Scalar>::postProcessing(DistrGeomState *geomState, GenDistrVector<Scalar> &extF, Corotator ***allCorot, double time,
                                       SysState<GenDistrVector<Scalar> > *distState, GenDistrVector<Scalar> *aeroF, DistrGeomState *refState,
                                       GenDistrVector<Scalar> *reactions, GenMDDynamMat<Scalar> *dynOps, GenDistrVector<Scalar> *resF)
{
  int numOutInfo = geoSource->getNumOutInfo();
  if(numOutInfo == 0) return;

  // get output information
  OutputInfo *oinfo = geoSource->getOutputInfo();

  // check if there are any output files which need to be processed now
  int outLimit = geoSource->getOutLimit();
  if(geoSource->noOutput(x) && !((x == 0) || (outLimit > 0 && x%outLimit == 0))) {
    x++; return;
  }

  if(domain->outFlag && domain->nodeTable == 0) domain->makeNodeTable(domain->outFlag);

  if(numOutInfo && x == 0 && !(domain->solInfo().isDynam() || domain->solInfo().timeIntegration == 1))
    filePrint(stderr," ... Postprocessing                 ...\n");
  if(!masterFlag) initPostPro();

  int iSub;

  // initialize and merge displacements from subdomains into cpu array
  DistSVec<Scalar, 11> disps_glo(*this->nodeInfo);
  DistSVec<Scalar, 11> masterDisps_glo(masterInfo);
  disps_glo = 0;
  DistSVec<Scalar, 11> *disps_loc = 0, *masterDisps_loc = 0;
  if(!domain->solInfo().basicDofCoords) {
    disps_loc = new DistSVec<Scalar, 11>(*this->nodeInfo);
    masterDisps_loc = new DistSVec<Scalar, 11>(masterInfo);
  }
  for(iSub = 0; iSub < this->numSub; ++iSub) {
    Scalar (*xyz)[11] = (Scalar (*)[11]) disps_glo.subData(iSub);//DofSet::max_known_nonL_dof
    Scalar (*xyz_loc)[11] = (disps_loc) ? (Scalar (*)[11]) disps_loc->subData(iSub) : 0;
    this->subDomain[iSub]->mergeDistributedNLDisp(xyz, (*geomState)[iSub], xyz_loc);
  }
  if((domain->solInfo().isCoupled && domain->solInfo().isMatching) || domain->GetnContactSurfacePairs()) {
    unify(disps_glo); // make sure master has solution for all dofs before reducing
    if(disps_loc) unify(*disps_loc);
  }
  disps_glo.reduce(masterDisps_glo, masterFlag, numFlags);
  if(disps_loc) disps_loc->reduce(*masterDisps_loc, masterFlag, numFlags);
  // initialize and merge aeroelastic forces
  DistSVec<Scalar, 6> aerof(*this->nodeInfo);
  DistSVec<Scalar, 6> masterAeroF(masterInfo);
  if(domain->solInfo().aeroFlag > -1 && aeroF) {
    GenDistrVector<Scalar> assembledAeroF(*aeroF);
    this->ba->assemble(assembledAeroF);
    aerof = 0;
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      Scalar (*mergedAeroF)[6] = (Scalar (*)[6]) aerof.subData(iSub);
      this->subDomain[iSub]->mergeDistributedForces(mergedAeroF, assembledAeroF.subData(iSub));
    }
    aerof.reduce(masterAeroF, masterFlag, numFlags);
  }
  // initialize and merge velocities & accelerations
  DistSVec<Scalar, 11> vels_glo(*this->nodeInfo), accs_glo(*this->nodeInfo);
  DistSVec<Scalar, 11> masterVels_glo(masterInfo), masterAccs_glo(masterInfo);
  DistSVec<Scalar, 11> *vels_loc = 0, *masterVels_loc = 0, *accs_loc = 0, *masterAccs_loc = 0;
  if(!domain->solInfo().basicDofCoords) {
    vels_loc = new DistSVec<Scalar, 11>(*this->nodeInfo);
    masterVels_loc = new DistSVec<Scalar, 11>(masterInfo);
    accs_loc = new DistSVec<Scalar, 11>(*this->nodeInfo);
    masterAccs_loc = new DistSVec<Scalar, 11>(masterInfo);
  }
  if(distState) {
    GenDistrVector<Scalar> *v_n = &distState->getVeloc();
    GenDistrVector<Scalar> *a_n = &distState->getAccel();
    vels_glo = 0; accs_glo = 0;
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      Scalar (*mergedVel)[11] = (Scalar (*)[11]) vels_glo.subData(iSub);
      Scalar (*mergedAcc)[11] = (Scalar (*)[11]) accs_glo.subData(iSub);
      double *vcx = this->subDomain[iSub]->getVcx();
      double *acx = this->subDomain[iSub]->getAcx();
      Scalar *vcx_scalar = new Scalar[this->subDomain[iSub]->numdof()];
      Scalar *acx_scalar = new Scalar[this->subDomain[iSub]->numdof()];
      Scalar (*mergedVel_loc)[11] = (vels_loc) ? (Scalar (*)[11]) vels_loc->subData(iSub) : 0;
      Scalar (*mergedAcc_loc)[11] = (accs_loc) ? (Scalar (*)[11]) accs_loc->subData(iSub) : 0;
      for(int i=0; i<this->subDomain[iSub]->numdof(); ++i) { vcx_scalar[i] = vcx[i]; acx_scalar[i] = acx[i]; }
      this->subDomain[iSub]->template mergeDistributedDisp<Scalar>(mergedVel, v_n->subData(iSub), vcx_scalar, mergedVel_loc);
      this->subDomain[iSub]->template mergeDistributedDisp<Scalar>(mergedAcc, a_n->subData(iSub), acx_scalar, mergedAcc_loc);
      delete [] vcx_scalar;
      delete [] acx_scalar;
    }
    vels_glo.reduce(masterVels_glo, masterFlag, numFlags);
    accs_glo.reduce(masterAccs_glo, masterFlag, numFlags);
    if(vels_loc) vels_loc->reduce(*masterVels_loc, masterFlag, numFlags);
    if(accs_loc) accs_loc->reduce(*masterAccs_loc, masterFlag, numFlags);
  }
  // initialize and merge reaction forces
  DistSVec<Scalar, 11> reacts(*this->nodeInfo);
  DistSVec<Scalar, 11> masterReacts(masterInfo);
  if(reactions) {
    GenDistrVector<Scalar> assembledReactions(*reactions);
    //TODO: this->xx->assemble(assembledReactions);
    reacts = 0;
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      Scalar (*mergedReactions)[11] = (Scalar (*)[11]) reacts.subData(iSub);
      this->subDomain[iSub]->mergeDistributedReactions(mergedReactions, assembledReactions.subData(iSub));
    }
    reacts.reduce(masterReacts, masterFlag, numFlags);
  }
  // initialize and merge residual
  DistSVec<Scalar, 6> resf(*this->nodeInfo);
  DistSVec<Scalar, 6> masterResF(masterInfo);
  if(resF) {
    GenDistrVector<Scalar> assembledResF(*resF);
    this->getSolVecAssembler()->assemble(assembledResF);
    resf = 0; 
    for(iSub = 0; iSub < this->numSub; ++iSub) { 
      Scalar (*mergedResF)[6] = (Scalar (*)[6]) resf.subData(iSub);
      this->subDomain[iSub]->mergeDistributedForces(mergedResF, assembledResF.subData(iSub));
    }
    resf.reduce(masterResF, masterFlag, numFlags);
  }
  // initialize and merge external force
  DistSVec<Scalar, 6> extf(*this->nodeInfo);
  DistSVec<Scalar, 6> masterExtF(masterInfo);
  if(geoSource->romExtForceOutput()) {
    GenDistrVector<Scalar> assembledExtF(extF);
    this->getSolVecAssembler()->assemble(assembledExtF);
    extf = 0;
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      Scalar (*mergedExtF)[6] = (Scalar (*)[6]) extf.subData(iSub);
      this->subDomain[iSub]->mergeDistributedForces(mergedExtF, assembledExtF.subData(iSub));
    }
    extf.reduce(masterExtF, masterFlag, numFlags);
  }
// RT - serialize the OUTPUT, PJSA - stress output doesn't work with serialized output. need to reconsider
#ifdef SERIALIZED_OUTPUT
for(int iCPU = 0; iCPU < this->communicator->size(); iCPU++) {
 this->communicator->sync();
 if(this->communicator->cpuNum() == iCPU) {
#endif

  // open binary output files
  if(x == 0) {
    if(!numRes) numRes = new int[numOutInfo];
    for(int i=0; i<numOutInfo; ++i) numRes[i] = 0;
  }

  if((x == 0) || (outLimit > 0 && x%outLimit == 0)) {
#ifdef DISTRIBUTED

    int *subToClus = geoSource->getSubToClus();
    int clusterId = (this->numSub > 0) ? subToClus[this->localSubToGl[0]] : 0;
    int firstCpuInCluster = (clusToCpu) ? (*clusToCpu)[clusterId][0] : 0;
    for(int iInfo = 0; iInfo < numOutInfo; iInfo++) {
      if(oinfo[iInfo].type == OutputInfo::Farfield || oinfo[iInfo].type == OutputInfo::AeroForce
         || oinfo[iInfo].type == OutputInfo::Energies || oinfo[iInfo].type == OutputInfo::DissipatedEnergy
         || oinfo[iInfo].type == OutputInfo::DeletedElements) {
        int oI = iInfo;
        if(this->firstOutput) { geoSource->openOutputFiles(0,&oI,1); }
        continue;
      }
      else if(oinfo[iInfo].nodeNumber == -1 && this->firstOutput) {
        if(this->communicator->cpuNum() == firstCpuInCluster && this->numSub > 0) geoSource->createBinaryOutputFile(iInfo,this->localSubToGl[0],x);
        else geoSource->computeAndCacheHeaderLength(iInfo);
      }
    }
#ifndef SERIALIZED_OUTPUT
    this->communicator->sync();
#endif
    for(int iInfo = 0; iInfo < numOutInfo; iInfo++) {
      if(oinfo[iInfo].nodeNumber == -1 && oinfo[iInfo].type != OutputInfo::Farfield && oinfo[iInfo].type != OutputInfo::AeroForce
         && oinfo[iInfo].type != OutputInfo::Energies && oinfo[iInfo].type != OutputInfo::DissipatedEnergy
         && oinfo[iInfo].type != OutputInfo::DeletedElements) {
        numRes[iInfo] = 0;
        for(iSub = 0; iSub < this->numSub; iSub++) {
          int glSub = this->localSubToGl[iSub];
          if(oinfo[iInfo].dataType == 1) {
            geoSource->outputRange(iInfo, masterFlag[iSub],
                       this->subDomain[iSub]->numNode(), glSub, nodeOffsets[iSub], x);
          }
          else {
            geoSource->outputRange(iInfo, this->subDomain[iSub]->getGlElems(),
                       this->subDomain[iSub]->numElements(), glSub,
                       elemOffsets[iSub], x);
          }
        }
      }
    }

    if(x == 0) { // always put single node output in one file regardless of outLimit
      for(iSub = 0; iSub < this->numSub; iSub++)  {
        int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
        if(nOutNodes) {
          if(this->firstOutput) geoSource->openOutputFiles(this->subDomain[iSub]->getOutputNodes(),
                   this->subDomain[iSub]->getOutIndex(), nOutNodes);
        }
      }
    }

#else
    if(x == 0 && this->firstOutput) geoSource->openOutputFiles();  // opens all output files
#endif
  }

  if (!nodePat) makeNodePat();

  for(int iOut = 0; iOut < numOutInfo; iOut++) {

    if(oinfo[iOut].interval == 0 || x % oinfo[iOut].interval != 0)
      continue;

    // set primal output states to either global or local
    DistSVec<Scalar, 11> &disps = (oinfo[iOut].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords)
                                  ? disps_glo : *disps_loc;
    DistSVec<Scalar, 11> &masterDisps = (oinfo[iOut].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords)
                                  ? masterDisps_glo : *masterDisps_loc;
    DistSVec<Scalar, 11> &vels = (oinfo[iOut].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords)
                                  ? vels_glo : *vels_loc;
    DistSVec<Scalar, 11> &masterVels = (oinfo[iOut].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords)
                                  ? masterVels_glo : *masterVels_loc;
    DistSVec<Scalar, 11> &accs = (oinfo[iOut].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords)
                                  ? accs_glo : *accs_loc;
    DistSVec<Scalar, 11> &masterAccs = (oinfo[iOut].oframe == OutputInfo::Global || domain->solInfo().basicDofCoords)
                                  ? masterAccs_glo : *masterAccs_loc;

    // update number of results
    numRes[iOut]++;
    switch(oinfo[iOut].type)  {

      case OutputInfo::FreqRespModes:
      case OutputInfo::Displacement:
        getPrimal(disps, masterDisps, time, x, iOut, 3, 0);
        break;
      case OutputInfo::Velocity:
        if(distState) getPrimal(vels, masterVels, time, x, iOut, 3, 0);
        break;
      case OutputInfo::Acceleration:
        if(distState) getPrimal(accs, masterAccs, time, x, iOut, 3, 0);
        break;
      case OutputInfo::RomResidual:
        if(resF) getPrimal(resf, masterResF, time, x, iOut, 3, 0);
        break;
      case OutputInfo::RomExtForce:
        getPrimal(extf, masterExtF, time, x, iOut, 3, 0);
        break;
      case OutputInfo::Disp6DOF:
        if(oinfo[iOut].rotvecouttype != OutputInfo::Euler || !oinfo[iOut].rescaling) {
          filePrint(stderr," *** WARNING: Output case %d not implemented\n", iOut);
          break;
        }
        getPrimal(disps, masterDisps, time, x, iOut, 6, 0);
        break;
      case OutputInfo::Velocity6:
        if(oinfo[iOut].angularouttype != OutputInfo::convected) {
          filePrint(stderr," *** WARNING: Output case %d not implemented\n", iOut);
          break;
        }
        if(distState) getPrimal(vels, masterVels, time, x, iOut, 6, 0);
        break;
      case OutputInfo::Accel6:
        if(oinfo[iOut].angularouttype != OutputInfo::convected) {
          filePrint(stderr," *** WARNING: Output case %d not implemented\n", iOut);
          break;
        }
        if(distState) getPrimal(accs, masterAccs, time, x, iOut, 6, 0);
        break;
      case OutputInfo::RomResidual6:
        if(resF) getPrimal(resf, masterResF, time, x, iOut, 6, 0);
        break;
      case OutputInfo::RomExtForce6:
        getPrimal(extf, masterExtF, time, x, iOut, 6, 0);
        break;
      case OutputInfo::Temperature:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 0);
        break;
      case OutputInfo::TemperatureFirstTimeDerivative:
        if(distState) getPrimal(vels, masterVels, time, x, iOut, 1, 6);
        break;
      case OutputInfo::StressXX:
        getStressStrain(geomState, allCorot, time, x, iOut, SXX, refState);
        break;
      case OutputInfo::StressYY:
        getStressStrain(geomState, allCorot, time, x, iOut, SYY, refState);
        break;
      case OutputInfo::StressZZ:
        getStressStrain(geomState, allCorot, time, x, iOut, SZZ, refState);
        break;
      case OutputInfo::StressXY:
        getStressStrain(geomState, allCorot, time, x, iOut, SXY, refState);
        break;
      case OutputInfo::StressYZ:
        getStressStrain(geomState, allCorot, time, x, iOut, SYZ, refState);
        break;
      case OutputInfo::StressXZ:
        getStressStrain(geomState, allCorot, time, x, iOut, SXZ, refState);
        break;
      case OutputInfo::StrainXX:
        getStressStrain(geomState, allCorot, time, x, iOut, EXX, refState);
        break;
      case OutputInfo::StrainYY:
        getStressStrain(geomState, allCorot, time, x, iOut, EYY, refState);
        break;
      case OutputInfo::StrainZZ:
        getStressStrain(geomState, allCorot, time, x, iOut, EZZ, refState);
        break;
      case OutputInfo::StrainXY:
        getStressStrain(geomState, allCorot, time, x, iOut, EXY, refState);
        break;
      case OutputInfo::StrainYZ:
        getStressStrain(geomState, allCorot, time, x, iOut, EYZ, refState);
        break;
      case OutputInfo::StrainXZ:
        getStressStrain(geomState, allCorot, time, x, iOut, EXZ, refState);
        break;
      case OutputInfo::StressVM:
        getStressStrain(geomState, allCorot, time, x, iOut, VON, refState);
        break;
      case OutputInfo::StrainVM:
        getStressStrain(geomState, allCorot, time, x, iOut, STRAINVON, refState);
        break;
      case OutputInfo::ContactPressure: {
        if(!domain->tdenforceFlag()) 
          getStressStrain(geomState, allCorot, time, x, iOut, CONPRESS, refState);
        else
          filePrint(stderr," *** WARNING: Output case %d not supported \n", iOut);
      } break;
      case OutputInfo::Damage:
        getStressStrain(geomState, allCorot, time, x, iOut, DAMAGE, refState);
        break;
      case OutputInfo::EquivalentPlasticStrain:
        getStressStrain(geomState, allCorot, time, x, iOut, EQPLSTRN, refState);
        break;
      case OutputInfo::Energies:
        this->getEnergies_b(geomState, extF, allCorot, iOut, time, distState, dynOps, aeroF);
        break;
      case OutputInfo::DissipatedEnergy:
        this->getDissipatedEnergy(geomState, allCorot, iOut, time);
        break;
      case OutputInfo::StressPR1:
        getPrincipalStress(geomState, allCorot, time, x, iOut, PSTRESS1, refState);
        break;
      case OutputInfo::StressPR2:
        getPrincipalStress(geomState, allCorot, time, x, iOut, PSTRESS2, refState);
        break;
      case OutputInfo::StressPR3:
        getPrincipalStress(geomState, allCorot, time, x, iOut, PSTRESS3, refState);
        break;
      case OutputInfo::StrainPR1:
        getPrincipalStress(geomState, allCorot, time, x, iOut, PSTRAIN1, refState);
        break;
      case OutputInfo::StrainPR2:
        getPrincipalStress(geomState, allCorot, time, x, iOut, PSTRAIN2, refState);
        break;
      case OutputInfo::StrainPR3:
        getPrincipalStress(geomState, allCorot, time, x, iOut, PSTRAIN3, refState);
        break;
      case OutputInfo::InXForce:
        getElementForce(geomState, allCorot, time, x, iOut, INX);
        break;
      case OutputInfo::InYForce:
        getElementForce(geomState, allCorot, time, x, iOut, INY);
        break;
      case OutputInfo::InZForce:
        getElementForce(geomState, allCorot, time, x, iOut, INZ);
        break;
      case OutputInfo::AXMoment:
        getElementForce(geomState, allCorot, time, x, iOut, AXM);
        break;
      case OutputInfo::AYMoment:
        getElementForce(geomState, allCorot, time, x, iOut, AYM);
        break;
      case OutputInfo::AZMoment:
        getElementForce(geomState, allCorot, time, x, iOut, AZM);
        break;
      case OutputInfo::DispX:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 0);
        break;
      case OutputInfo::DispY:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 1);
        break;
      case OutputInfo::DispZ:
        getPrimal(disps, masterDisps, time, x, iOut, 1, 2);
        break;
      case OutputInfo::RotX:
        if(oinfo[iOut].rotvecouttype != OutputInfo::Euler || !oinfo[iOut].rescaling) {
          filePrint(stderr," *** WARNING: Output case %d not implemented\n", iOut);
          break;
        }
        getPrimal(disps, masterDisps, time, x, iOut, 1, 3);
        break;
      case OutputInfo::RotY:
        if(oinfo[iOut].rotvecouttype != OutputInfo::Euler || !oinfo[iOut].rescaling) {
          filePrint(stderr," *** WARNING: Output case %d not implemented\n", iOut);
          break;
        }
        getPrimal(disps, masterDisps, time, x, iOut, 1, 4);
        break;
      case OutputInfo::RotZ:
        if(oinfo[iOut].rotvecouttype != OutputInfo::Euler || !oinfo[iOut].rescaling) {
          filePrint(stderr," *** WARNING: Output case %d not implemented\n", iOut);
          break;
        }
        getPrimal(disps, masterDisps, time, x, iOut, 1, 5);
        break;
      case OutputInfo::DispMod:
        if(oinfo[iOut].nodeNumber == -1) {
          for(iSub = 0; iSub < this->numSub; ++iSub) {
            int size = masterDisps.subSize(iSub);
            Scalar (*xyz)[11] = (Scalar (*)[11]) masterDisps.subData(iSub);
            Scalar *dispMod = new Scalar[size];
            for(int iNode=0; iNode<size; ++iNode) {
              dispMod[iNode] = ScalarTypes::sqrt(xyz[iNode][0]*xyz[iNode][0] +
                                                 xyz[iNode][1]*xyz[iNode][1] +
                                                 xyz[iNode][2]*xyz[iNode][2]);
            }
            geoSource->writeNodeScalarToFile(dispMod, size, this->localSubToGl[iSub], nodeOffsets[iSub],
                                             iOut, x, numRes[iOut], time, 1, masterFlag[iSub]);
            delete [] dispMod;
          }
        }
        else {
          for(iSub = 0; iSub < this->numSub; ++iSub) {
            int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
            if(nOutNodes) {
              int *outIndex = this->subDomain[iSub]->getOutIndex();
              for(int iNode = 0; iNode < nOutNodes; iNode++) {
                if(outIndex[iNode] == iOut) {
                  Scalar (*xyz)[11] = (Scalar (*)[11]) disps.subData(iSub);
                  int *outNodes = this->subDomain[iSub]->getOutputNodes();
                  Scalar dispMod = ScalarTypes::sqrt(xyz[outNodes[iNode]][0]*xyz[outNodes[iNode]][0] +
                                                     xyz[outNodes[iNode]][1]*xyz[outNodes[iNode]][1] +
                                                     xyz[outNodes[iNode]][2]*xyz[outNodes[iNode]][2]);
                  geoSource->outputNodeScalars(iOut, &dispMod, 1, time);
                }
              }
            }
          }
        }
        break;
      case OutputInfo::RotMod:
        if(oinfo[iOut].rotvecouttype != OutputInfo::Euler || !oinfo[iOut].rescaling) {
          filePrint(stderr," *** WARNING: Output case %d not implemented\n", iOut);
          break;
        }
        if(oinfo[iOut].nodeNumber == -1) {
          for(iSub = 0; iSub < this->numSub; ++iSub) {
            int size = masterDisps.subSize(iSub);
            Scalar (*xyz)[11] = (Scalar (*)[11]) masterDisps.subData(iSub);
            Scalar *rotMod = new Scalar[size];
            for(int iNode=0; iNode<size; ++iNode) {
              rotMod[iNode] = ScalarTypes::sqrt(xyz[iNode][3]*xyz[iNode][3] +
                                                xyz[iNode][4]*xyz[iNode][4] +
                                                xyz[iNode][5]*xyz[iNode][5]);
            }
            geoSource->writeNodeScalarToFile(rotMod, size, this->localSubToGl[iSub], nodeOffsets[iSub],
                                             iOut, x, numRes[iOut], time, 1, masterFlag[iSub]);
            delete [] rotMod;
          }
        }
        else {
          for(iSub = 0; iSub < this->numSub; ++iSub) {
            int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
            if(nOutNodes) {
              int *outIndex = this->subDomain[iSub]->getOutIndex();
              for(int iNode = 0; iNode < nOutNodes; iNode++) {
                if(outIndex[iNode] == iOut) {
                  Scalar (*xyz)[11] = (Scalar (*)[11]) disps.subData(iSub);
                  int *outNodes = this->subDomain[iSub]->getOutputNodes();
                  Scalar rotMod = ScalarTypes::sqrt(xyz[outNodes[iNode]][3]*xyz[outNodes[iNode]][3] +
                                                    xyz[outNodes[iNode]][4]*xyz[outNodes[iNode]][4] +
                                                    xyz[outNodes[iNode]][5]*xyz[outNodes[iNode]][5]);
                  geoSource->outputNodeScalars(iOut, &rotMod, 1, time);
                }
              }
            }
          }
        }
        break;
      case OutputInfo::TotMod:
        if(oinfo[iOut].rotvecouttype != OutputInfo::Euler || !oinfo[iOut].rescaling) {
          filePrint(stderr," *** WARNING: Output case %d not implemented\n", iOut);
          break;
        }
        if(oinfo[iOut].nodeNumber == -1) {
          for(iSub = 0; iSub < this->numSub; ++iSub) {
            int size = masterDisps.subSize(iSub);
            Scalar (*xyz)[11] = (Scalar (*)[11]) masterDisps.subData(iSub);
            Scalar *totMod = new Scalar[size];
            for(int iNode=0; iNode<size; ++iNode) {
              totMod[iNode] = ScalarTypes::sqrt(xyz[iNode][0]*xyz[iNode][0] +
                                                xyz[iNode][1]*xyz[iNode][1] +
                                                xyz[iNode][2]*xyz[iNode][2] +
                                                xyz[iNode][3]*xyz[iNode][3] +
                                                xyz[iNode][4]*xyz[iNode][4] +
                                                xyz[iNode][5]*xyz[iNode][5]);
            }
            geoSource->writeNodeScalarToFile(totMod, size, this->localSubToGl[iSub], nodeOffsets[iSub],
                                             iOut, x, numRes[iOut], time, 1, masterFlag[iSub]);
            delete [] totMod;
          }
        }
        else {
          for(iSub = 0; iSub < this->numSub; ++iSub) {
            int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
            if(nOutNodes) {
              int *outIndex = this->subDomain[iSub]->getOutIndex();
              for(int iNode = 0; iNode < nOutNodes; iNode++) {
                if(outIndex[iNode] == iOut) {
                  Scalar (*xyz)[11] = (Scalar (*)[11]) disps.subData(iSub);
                  int *outNodes = this->subDomain[iSub]->getOutputNodes();
                  Scalar totMod = ScalarTypes::sqrt(xyz[outNodes[iNode]][0]*xyz[outNodes[iNode]][0] +
                                                    xyz[outNodes[iNode]][1]*xyz[outNodes[iNode]][1] +
                                                    xyz[outNodes[iNode]][2]*xyz[outNodes[iNode]][2] +
                                                    xyz[outNodes[iNode]][3]*xyz[outNodes[iNode]][3] +
                                                    xyz[outNodes[iNode]][4]*xyz[outNodes[iNode]][4] +
                                                    xyz[outNodes[iNode]][5]*xyz[outNodes[iNode]][5]);
                  geoSource->outputNodeScalars(iOut, &totMod, 1, time);
                }
              }
            }
          }
        }
        break;
      case OutputInfo::AeroForce: break; // this is done in DistFlExchange.C
      case OutputInfo::AeroXForce:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 0);
        break;
      case OutputInfo::AeroYForce:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 1);
        break;
      case OutputInfo::AeroZForce:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 2);
        break;
      case OutputInfo::AeroXMom:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 3);
        break;
      case OutputInfo::AeroYMom:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 4);
        break;
      case OutputInfo::AeroZMom:
        if(aeroF) getAeroForceScalar(aerof, masterAeroF, time, x, iOut, 5);
        break;
      case OutputInfo::Reactions:
        if(reactions) getPrimal(reacts, masterReacts, time, x, iOut, 3, 0);
        break;
      case OutputInfo::Reactions6:
        if(reactions) getPrimal(reacts, masterReacts, time, x, iOut, 6, 0);
        break;
      case OutputInfo::TDEnforcement: {
        if(domain->tdenforceFlag()) {
          DistSVec<double, 1> all_data(*this->nodeInfo);
          if(oinfo[iOut].tdenforc_var == 1) all_data = 0.5;
          else all_data = 0;
          double **sub_data = new double * [this->numSub];
          for(iSub = 0; iSub < this->numSub; ++iSub) sub_data[iSub] = (double *) all_data.subData(iSub);
          for(int iMortar=0; iMortar<domain->GetnMortarConds(); iMortar++) {
            domain->GetMortarCond(iMortar)->get_plot_variable(oinfo[iOut].tdenforc_var, sub_data, this->numSub, this->subDomain);
          }
          DistSVec<double, 1> master_data(masterInfo);
          all_data.reduce(master_data, masterFlag, numFlags);
          for(iSub = 0; iSub < this->numSub; ++iSub) {
            geoSource->writeNodeScalarToFile((double *) master_data.subData(iSub), master_data.subSize(iSub), this->localSubToGl[iSub], nodeOffsets[iSub],
                                             iOut, x, numRes[iOut], time, 1, masterFlag[iSub]);
          }
          delete [] sub_data;
        }
        else filePrint(stderr," *** WARNING: Output case %d not supported \n", iOut);
      } break;
      case OutputInfo::DeletedElements:
        getDeletedElements(iOut);
        break;
      case OutputInfo::Statevector:
      case OutputInfo::Velocvector:
      case OutputInfo::Accelvector:
      case OutputInfo::InternalStateVar:
      case OutputInfo::DualStateVar:
      case OutputInfo::MuStateVar: 
      case OutputInfo::Forcevector:
      case OutputInfo::Constraintvector:
      case OutputInfo::Constraintviolation:
      case OutputInfo::Residual:
      case OutputInfo::Jacobian:
      case OutputInfo::RobData:
      case OutputInfo::SampleMesh:
        break;
      default:
        filePrint(stderr," *** WARNING: Output case %d not implemented\n", iOut);
        break;
    }
  }
  x++;
#ifdef SERIALIZED_OUTPUT
 }
}
this->communicator->sync();
#endif
  if(disps_loc) delete disps_loc;
  if(masterDisps_loc) delete masterDisps_loc;
  if(masterVels_loc) delete masterVels_loc;
  if(accs_loc) delete accs_loc;
  if(masterAccs_loc) delete masterAccs_loc;
}

template<class Scalar>
void
GenDistrDomain<Scalar>::getStressStrain(DistrGeomState *gs, Corotator ***allCorot, double time,
                                        int x, int fileNumber, int Findex, DistrGeomState *refState)
{
  // non-linear version of getStressStrain
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];

  if(oinfo.averageFlg == 0) {
    getElementStressStrain(gs, allCorot, time, x, fileNumber, Findex, refState);
    return;
  }

  DistVec<Scalar> stress(*this->nodeInfo);
  DistVec<Scalar> weight(*this->nodeInfo);

  stress = 0;
  weight = 0;

  int iSub;

  // each subdomain computes its stress vector
  for (iSub = 0; iSub < this->numSub; ++iSub) {
    if(Findex != 16) {
      GeomState *subRefState = (refState) ? (*refState)[iSub] : 0;
      this->subDomain[iSub]->computeStressStrain((*gs)[iSub], allCorot[iSub], fileNumber,
                                                 Findex, stress.subData(iSub), weight.subData(iSub), subRefState);
    }
    else {
      this->subDomain[iSub]->computeLocalContactPressure(stress.subData(iSub), weight.subData(iSub));
    }
  }

  paralApply(this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>, this->nodePat, &stress);
  this->nodePat->exchange();
  paralApply(this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, this->nodePat, &stress);

  paralApply(this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>, this->nodePat, &weight);
  this->nodePat->exchange();
  paralApply(this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, this->nodePat, &weight);

  // Divide stress by weight
  for (iSub = 0; iSub < this->numSub; ++iSub)  {
    Vec<Scalar> &locStress = stress(iSub);
    Vec<Scalar> &locWeight = weight(iSub);
    for(int i = 0; i < stress.subSize(iSub); ++i)
      if(locWeight[i] != 0.0)
        locStress[i] /= locWeight[i];
      else
        locStress[i] = 0.0;
  }

  // reduce the stress vector to just the master quantities
  DistVec<Scalar> masterStress(masterInfo);
  stress.reduce(masterStress, masterFlag, numFlags);

  // write to file
  if(oinfo.nodeNumber == -1) { // output binary or ascii data for all nodes or node group
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      geoSource->writeNodeScalarToFile(masterStress.subData(iSub), masterStress.subSize(iSub), this->localSubToGl[iSub],
                                       nodeOffsets[iSub], fileNumber, x, numRes[fileNumber], time, 1, masterFlag[iSub]); 
    }
  }
  else { // output ascii data for one node
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
      if(nOutNodes) {
        int *outIndex = this->subDomain[iSub]->getOutIndex();
        for(int iNode = 0; iNode < nOutNodes; iNode++) {
          if(outIndex[iNode] == fileNumber) {
            Scalar *nodeStress = (Scalar *) stress.subData(iSub);
            int *outNodes = this->subDomain[iSub]->getOutputNodes();
            geoSource->outputNodeScalars(fileNumber, nodeStress+outNodes[iNode], 1, time);
          }
        }
      }
    }
  }

}
  
template<class Scalar>
void
GenDistrDomain<Scalar>::getElementStressStrain(DistrGeomState *gs, Corotator ***allCorot, double time,
                                               int iter, int fileNumber, int Findex, DistrGeomState *refState)
{
  // non-linear version of getElementStressStrain
  for(int iSub = 0; iSub < this->numSub; iSub++)  {
    int numElemNodes = this->subDomain[iSub]->countElemNodes();
    Scalar *elemStress = new Scalar[numElemNodes];
    GeomState *subRefState = (refState) ? (*refState)[iSub] : 0;
    this->subDomain[iSub]->computeStressStrain((*gs)[iSub], allCorot[iSub], fileNumber, Findex, elemStress, (Scalar *) 0, subRefState);
    geoSource->writeElemScalarToFile(elemStress, numElemNodes, this->localSubToGl[iSub], elemNodeOffsets[iSub], fileNumber, iter,
                                     numRes[fileNumber], time, this->elemToNode->numConnect(), this->subDomain[iSub]->getGlElems());
    delete [] elemStress;
  }
}

template<class Scalar>
void 
GenDistrDomain<Scalar>::getPrincipalStress(DistrGeomState *gs, Corotator ***allCorot, double time,  
                                           int x, int fileNumber, int strIndex, DistrGeomState *refState)
{
  OutputInfo &oinfo = geoSource->getOutputInfo()[fileNumber];
  if(oinfo.averageFlg == 0) {
    getElementPrincipalStress(gs, allCorot, time, x, fileNumber, strIndex, refState);
    return;
  }

  // non-linear version of getPrincipalStress
  // set stress VS. strain for element subroutines
  int i, j;
  int strInd;
  int stressORstrain;
  int strDir[6];

  if((strIndex==0) || (strIndex==1) || (strIndex==2)) {
    strInd = strIndex;
    stressORstrain = 0;
    for(i=0; i<6; ++i)
      strDir[i] = i;
  }
  else if((strIndex==3) || (strIndex==4) || (strIndex==5)) {
    strInd = strIndex-3;
    stressORstrain = 1;
    for(i=0; i<6; ++i)
      strDir[i] = i+7;
  }
  DistVec<Scalar> **stress = new DistVec<Scalar>*[6];
  DistVec<Scalar> weight(*this->nodeInfo);

  int iSub;
  int str_loop;
  for(str_loop = 0; str_loop < 6; ++str_loop)
    stress[str_loop] = new DistVec<Scalar> (*this->nodeInfo);  

  // each subdomain computes its stress/strain vector

  // Compute Each Required Stress (all 6) using same routines as for
  // individual stresses
  int Findex;

  for(str_loop = 0; str_loop < 6; ++str_loop) {

    // get current stress/strain index
    Findex = strDir[str_loop];

    // Initialize distributed vector to zero
    *stress[str_loop] = 0;
    weight = 0;

    for(iSub = 0; iSub < this->numSub; ++iSub) {
      GeomState *subRefState = (refState) ? (*refState)[iSub] : 0;
      this->subDomain[iSub]->computeStressStrain((*gs)[iSub], allCorot[iSub], fileNumber, 
                Findex, stress[str_loop]->subData(iSub), weight.subData(iSub), subRefState);
    }

    paralApply(this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>,
               this->nodePat, stress[str_loop]);
    this->nodePat->exchange();
    paralApply(this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, this->nodePat, stress[str_loop]);

    paralApply(this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<Scalar>, this->nodePat, &weight);
    this->nodePat->exchange();
    paralApply(this->subDomain, &GenSubDomain<Scalar>::template addNodalData<Scalar>, this->nodePat, &weight);

    // Divide stress by weight
    for(iSub = 0; iSub < this->numSub; ++iSub)  {
      Vec<Scalar> &locStress = (*stress[str_loop])(iSub);
      Vec<Scalar> &locWeight = weight(iSub);
      for(int i = 0; i < stress[str_loop]->subSize(iSub); ++i)  {
        if(locWeight[i] != 0.0)
          locStress[i] /= locWeight[i];
        else
          locStress[i] = 0.0;
      }
    }
  }

  // Calculate Principals at each node
  Scalar svec[6], pvec[3];
  DistVec<Scalar> allPVec(*this->nodeInfo);
  for(iSub = 0; iSub < this->numSub; ++iSub) {

    Vec<Scalar> &locPVec = allPVec(iSub);
    for(i = 0; i < this->subDomain[iSub]->numNode(); ++i) {

      for(j = 0; j < 6; ++j)
        svec[j] = stress[j]->subData(iSub)[i];
 
      // Convert Engineering to Tensor Strains
      if(stressORstrain != 0) {
        svec[3] /= 2;
        svec[4] /= 2;
        svec[5] /= 2;
      }
      pstress(svec,pvec);
      locPVec[i] = pvec[strInd];
    }
  }

  // reduce stress vector to master quantities
  DistVec<Scalar> masterPVec(masterInfo);
  allPVec.reduce(masterPVec, masterFlag, numFlags);

  // print out data
  if(oinfo.nodeNumber == -1) { // output binary or ascii data for all nodes or node group
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      geoSource->writeNodeScalarToFile(masterPVec.subData(iSub), masterPVec.subSize(iSub), this->localSubToGl[iSub],
                                       nodeOffsets[iSub], fileNumber, x, numRes[fileNumber], time, 1, masterFlag[iSub]); 
    }
  }
  else { // output ascii data for one node
    for(iSub = 0; iSub < this->numSub; ++iSub) {
      int nOutNodes = this->subDomain[iSub]->getNumNodalOutput();
      if(nOutNodes) {
        int *outIndex = this->subDomain[iSub]->getOutIndex();
        for(int iNode = 0; iNode < nOutNodes; iNode++) {
          if(outIndex[iNode] == fileNumber) {
            Scalar *nodeStress = (Scalar *) allPVec.subData(iSub);
            int *outNodes = this->subDomain[iSub]->getOutputNodes();
            geoSource->outputNodeScalars(fileNumber, nodeStress+outNodes[iNode], 1, time);
          }
        }
      }
    }
  }
}

template<class Scalar>
void
GenDistrDomain<Scalar>::makeNodePat()
{
  nodePat = new FSCommPattern<Scalar>(
      this->communicator, this->cpuToSub.get(), this->myCPU, FSCommPattern<Scalar>::CopyOnSend);
  for(int i=0; i<this->numSub; ++i)
    this->subDomain[i]->setNodeCommSize(nodePat);
  nodePat->finalize();
}

template<class Scalar>
void
GenDistrDomain<Scalar>::unify(DistSVec<Scalar, 11> &vec)
{
  // make sure that every subdomain sharing a node has the same solution for all the dofs. Actually the master sub is really the only one that needs it.
  FSCommPattern<Scalar> *pat =
      new FSCommPattern<Scalar>(this->communicator, this->cpuToSub.get(), this->myCPU, FSCommPattern<Scalar>::CopyOnSend);
  for(int i=0; i<this->numSub; ++i)
    this->subDomain[i]->setNodeCommSize(pat, 11);
  pat->finalize();
  for(int i=0; i<this->numSub; ++i)
    this->subDomain[i]->sendNode((Scalar (*)[11]) vec.subData(i), pat);
  pat->exchange();
  for(int i=0; i<this->numSub; ++i)
    this->subDomain[i]->collectNode((Scalar (*)[11]) vec.subData(i), pat);
  delete pat;
}

template<class Scalar>
void GenDistrDomain<Scalar>::getElementAttr(int fileNumber,int iAttr, double time)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
  int avgnum = oinfo[fileNumber].averageFlg;
  if(avgnum == 0) 
    {
      assert(0); // not implemented
      return;
    }

  DistVec<double> props(*this->nodeInfo);
  DistVec<double> weight(*this->nodeInfo);

  // Initialize distributed vector to zero
  props  = 0;
  weight = 0;
  for(int iSub = 0; iSub < this->numSub; ++iSub) 
    { this->getSubDomain(iSub)->mergeElemProps(props.subData(iSub), weight.subData(iSub), iAttr); }
  
  paralApply(this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<double>,
	     this->nodePat, &props);
  this->nodePat->exchange();
  paralApply(this->subDomain, &GenSubDomain<Scalar>::template addNodalData<double>, this->nodePat, &props);
  
  paralApply(this->subDomain, &GenSubDomain<Scalar>::template dispatchNodalData<double>, this->nodePat, &weight);
  this->nodePat->exchange();
  paralApply(this->subDomain, &GenSubDomain<Scalar>::template addNodalData<double>, this->nodePat, &weight);
  
  // average
  for(int iSub = 0; iSub < this->numSub; ++iSub)  
    {
      Vec<double> &locProps =  props (iSub);
      Vec<double> &locWeight = weight(iSub);
      for(int i = 0; i < props.subSize(iSub); ++i)  
	{
	  if (locWeight[i] != 0.0)
	    { locProps[i] /= locWeight[i]; }
	  else
	    { locProps[i] = 0.0; }
	}
    }
  // reduce props vector to master quantities
  DistVec<double> masterVec(masterInfo);
  props.reduce(masterVec, masterFlag, numFlags);
  
  // print out data
  for(int iSub = 0; iSub < this->numSub; ++iSub) {
    geoSource->writeNodeScalarToFile(masterVec.subData(iSub),
                                     masterVec.subSize(iSub), this->localSubToGl[iSub],
                                     nodeOffsets[iSub], fileNumber, x, numRes[fileNumber], time, 1, masterFlag[iSub]);
  }
  return;
}

template<class Scalar>
void
GenDistrDomain<Scalar>::getDeletedElements(int iOut)
{
  OutputInfo *oinfo = geoSource->getOutputInfo();
#ifdef SERIALIZED_OUTPUT
  for(int iSub = 0; iSub < this->numSub; ++iSub) {
    std::vector<std::pair<double,int> > &deletedElements = this->subDomain[iSub]->getDeletedElements();
    for(std::vector<std::pair<double,int> >::iterator it = deletedElements.begin(); it != deletedElements.end(); ++it) {
      filePrint(oinfo[i].filptr, " %12.6e  %9d          Undetermined\n", it->first, it->second+1);
      fflush(oinfo[i].filptr);
    }
    deletedElements.clear();
  }
#else
  std::vector<std::pair<double,int> > localDeletedElements;
  for(int iSub = 0; iSub < this->numSub; ++iSub) {
    std::vector<std::pair<double,int> > &deletedElements = this->subDomain[iSub]->getDeletedElements();
    localDeletedElements.insert(localDeletedElements.end(), deletedElements.begin(), deletedElements.end());
    deletedElements.clear();
  }
  int localCount = localDeletedElements.size();
  int *recvbuf = new int[this->communicator->size()];
  this->communicator->gather(&localCount, 1, recvbuf, 1);
  int globalCount;
  if(this->communicator->cpuNum() == 0) {
    globalCount = 0;
    for(int i=0; i<this->communicator->size(); ++i) globalCount += recvbuf[i];
  }
  this->communicator->broadcast(1, &globalCount);
  if(globalCount > 0) {
    int *sendbuf2 = new int[localCount];
    double *sendbuf3 = new double[localCount];
    if(localCount > 0) {
      int i=0;
      for(std::vector<std::pair<double,int> >::iterator it = localDeletedElements.begin(); it != localDeletedElements.end(); ++it, ++i) {
        sendbuf2[i] = it->second;
        sendbuf3[i] = it->first;
      }
    }
    int *recvbuf2, *displs;
    double *recvbuf3;
    if(this->communicator->cpuNum() == 0) {
      recvbuf2 = new int[globalCount];
      recvbuf3 = new double[globalCount];
      displs = new int[this->communicator->size()];
      displs[0] = 0;
      for(int i=1; i<this->communicator->size(); ++i) {
        displs[i] = displs[i-1] + recvbuf[i-1];
      }
    }
    this->communicator->gatherv(sendbuf2, localCount, recvbuf2, recvbuf, displs);
    this->communicator->gatherv(sendbuf3, localCount, recvbuf3, recvbuf, displs);
    delete [] sendbuf2;
    delete [] sendbuf3;
    if(this->communicator->cpuNum() == 0) {
      for(int i=0; i<globalCount; ++i) {
        filePrint(oinfo[iOut].filptr, " %12.6e  %9d          Undetermined\n", recvbuf3[i], recvbuf2[i]+1);
        fflush(oinfo[iOut].filptr);
      }
      delete [] recvbuf2;
      delete [] recvbuf3;
      delete [] displs;
    }
  }
  delete [] recvbuf;
#endif
}

