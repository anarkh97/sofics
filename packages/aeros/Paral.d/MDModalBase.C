#include <Paral.d/MDModalBase.h>
#include <Paral.d/MDOp.h>
#include <Problems.d/ModalBase.h>
#include <Driver.d/Domain.h>
#include <Utils.d/ModeData.h>
#include <Driver.d/DecDomain.h>
#ifdef DISTRIBUTED
#include <Dist.d/DistDom.h>
#endif
#include <Driver.d/SysState.h>
#include <Utils.d/DistHelper.h>
#include <Solvers.d/MultiDomainRbm.h>

extern ModeData modeData;

//------------------------------------------------------------------------------
//******************************************************************************
//------------------------------------------------------------------------------

void MDModalBase::preProcessBase(){
/*PRE: none
 POST: initialize and populate various data members
*/
  numModes = numRBM = numFlex = 0;
#ifdef DISTRIBUTED
  decDomain = new GenDistrDomain<double>(domain);
#else
  decDomain = new GenDecDomain<double>(domain);
#endif
  decDomain->preProcess();
  MultiDomainOp mdop(&MultiDomainOp::makeAllDOFs, decDomain->getAllSubDomains());
  execParal(decDomain->getNumSub(), &mdop, &MultiDomainOp::runFor);

  numModes = modeData.numModes;  // this is assigment is temporary; it is
                                 //   overrulled in populateFlexModes
  fullTmpF = new DistrVector(decDomain->solVecInfo());
  fullTmpGrav = new DistrVector(decDomain->solVecInfo());
  fullAeroF = new DistrVector(decDomain->solVecInfo());
  fullDsp = new DistrVector(decDomain->solVecInfo());
  fullVel = new DistrVector(decDomain->solVecInfo());
  fullAcc = new DistrVector(decDomain->solVecInfo());
  fullPrevVel = new DistrVector(decDomain->solVecInfo());

  fullTmpF->zero();
  fullTmpGrav->zero();
  fullAeroF->zero();
  fullDsp->zero();
  fullVel->zero();
  fullAcc->zero();
  fullPrevVel->zero();

  prevFrc = new DistrVector(decDomain->solVecInfo());
  prevFrcBackup = new DistrVector(decDomain->solVecInfo());
}

//------------------------------------------------------------------------------

void MDModalBase::populateRBModes(){

  // compute the rigid body modes
  MultiDomainRbm<double> *rbm = decDomain->constructRbm();
  numRBM   = rbm->numRBM();
  modesRB  = new GenDistrVectorSet<double>(numRBM, decDomain->solVecInfo());

  rbm->getRBMs(*modesRB);
}

//------------------------------------------------------------------------------

void MDModalBase::populateFlexModes(double scale, bool readAll){
/*PRE: preProcessBase has been called
       scale has default value of 1.0; readAll has default value of 0
       the modes in modeData are : mode^T.M.mode = identity
 POST: populate modesFl with data from modeData multiplied by scale
       populate freqs with the circular frequency of each flexible mode
 NOTE: if readAll, then also include in modesFL those modes with zero frequency
*/
  // count the number of flexible modes
  numFlex = 0;
  int iMode;
  for(iMode = 0; iMode < numModes; ++iMode){
    if( readAll || (modeData.frequencies[iMode] != 0.0) )
      ++numFlex;
  }
  modesFl = new GenDistrVectorSet<double>(numFlex, decDomain->solVecInfo());
  freqs   = new double[numFlex];

  int iModeFl = 0, iNode, dof;
  for(iMode = 0; iMode < numModes; ++iMode){
    if( readAll || (modeData.frequencies[iMode] != 0.0) ){
      (*modesFl)[iModeFl].zero();
      freqs[iModeFl] = modeData.frequencies[iMode] * 8 * atan(1.);
      for(int iSub = 0; iSub < decDomain->getNumSub(); ++iSub) {

        SubDomain *sd = decDomain->getSubDomain(iSub);
        for(iNode = 0; iNode < sd->numNodes(); ++iNode){

          dof = sd->getCDSA()->locate(iNode, DofSet::Xdisp);
          if(dof >= 0) (*modesFl)[iModeFl].subData(iSub)[dof] = scale*modeData.modes[iMode][sd->localToGlobal(iNode)][0];

          dof = sd->getCDSA()->locate(iNode, DofSet::Ydisp);
          if(dof >= 0) (*modesFl)[iModeFl].subData(iSub)[dof] = scale*modeData.modes[iMode][sd->localToGlobal(iNode)][1];

          dof = sd->getCDSA()->locate(iNode, DofSet::Zdisp);
          if(dof >= 0) (*modesFl)[iModeFl].subData(iSub)[dof] = scale*modeData.modes[iMode][sd->localToGlobal(iNode)][2];

          dof = sd->getCDSA()->locate(iNode, DofSet::Xrot);
          if(dof >= 0) (*modesFl)[iModeFl].subData(iSub)[dof] = scale*modeData.modes[iMode][sd->localToGlobal(iNode)][3];

          dof = sd->getCDSA()->locate(iNode, DofSet::Yrot);
          if(dof >= 0) (*modesFl)[iModeFl].subData(iSub)[dof] = scale*modeData.modes[iMode][sd->localToGlobal(iNode)][4];

          dof = sd->getCDSA()->locate(iNode, DofSet::Zrot);
          if(dof >= 0) (*modesFl)[iModeFl].subData(iSub)[dof] = scale*modeData.modes[iMode][sd->localToGlobal(iNode)][5];
        }
      }
      ++iModeFl;
    }
  }
}

//------------------------------------------------------------------------------

void MDModalBase::initStateBase(Vector& dsp, Vector& vel,
  Vector& acc, Vector& vel_p, int idxOffset){

  dsp.zero();
  vel.zero();
  acc.zero();
  vel_p.zero();

  ControlInfo *cinfo = geoSource->getCheckFileInfo();
  SolverInfo &sinfo = domain->solInfo();

  if (cinfo->lastRestartFile) {
     filePrint(stderr, " ... Restarting From a Previous Run ...\n");
     int fn = open(cinfo->lastRestartFile,O_RDONLY );
     if(fn >= 0) {
       int vsize, restartTIndex;
       double restartT;
       int readSize = read(fn, &vsize, sizeof(int));
       if(vsize != dsp.size() || readSize != sizeof(int))
         filePrint(stderr," *** ERROR: Inconsistent restart file\n");
       readSize = read(fn, &restartTIndex, sizeof(int));
       if(readSize != sizeof(int))
         filePrint(stderr," *** ERROR: Inconsistent restart file\n");
       if(strcmp(cinfo->FlagRST,"new") == 0)
         sinfo.initialTimeIndex = 0;
       else
         sinfo.initialTimeIndex = restartTIndex;
       readSize = read(fn, &restartT, sizeof(double));
       if(readSize != sizeof(double))
         filePrint(stderr," *** ERROR: Inconsistent restart file\n");
       if(strcmp(cinfo->FlagRST,"new") == 0)
         sinfo.initialTime = 0.0;
       else
         sinfo.initialTime = restartT;
       readSize = read(fn, dsp.data(), vsize*sizeof(double));
       if(int(readSize) != int(vsize*sizeof(double)))
         filePrint(stderr," *** ERROR: Inconsistent restart file\n");

       readSize = read(fn, vel.data(), vsize*sizeof(double));
       if(int(readSize) != int(vsize*sizeof(double)))
         filePrint(stderr," *** ERROR: Inconsistent restart file\n");

       readSize = read(fn, vel_p.data(), vsize*sizeof(double));
       if(int(readSize) != int(vsize*sizeof(double))) {
         filePrint(stderr," *** WARNING: Older version of restart file"
                          " -- Missing velocity field is set to zero\n");
         vel_p.zero();
       }

       close(fn);
     } 
     else {
       filePrint(stderr, " *** ERROR: Restart file could not be opened\n");
     }
  }
  else  {

    for(int j = 0; j <  domain->numInitVelocityModal(); ++j) {
      vel[domain->getInitVelocityModal()[j].nnum + idxOffset] += domain->getInitVelocityModal()[j].val;
    }

    for(int j = 0; j <  domain->numInitDispModal(); ++j) {
      dsp[domain->getInitDispModal()[j].nnum + idxOffset] += domain->getInitDispModal()[j].val;
    }

    // superimpose the non-modal initial velocity and/or displacement
    if(domain->numInitVelocity() > 0 || ((domain->numInitDisp() > 0 || domain->numInitDisp6() > 0) && sinfo.zeroInitialDisp == 0)) {
      // check if basis is mass-normalized or not
      ModalParams &modalParams = domain->solInfo().readInModes[domain->solInfo().modal_id.front()];
      bool massNormalizedBasis = (modalParams.type == ModalParams::Eigen || modalParams.type == ModalParams::Mnorm ||
                                  modalParams.type == ModalParams::Undefined);

      GenDistrVectorSet<double> tPhiM(numFlex+numRBM, decDomain->solVecInfo());

      if(massNormalizedBasis) {
        // construct and assemble full mass matrix
        SparseMatrix **subM = new SparseMatrix * [decDomain->getNumSub()];
        for(int iSub = 0; iSub < decDomain->getNumSub(); ++iSub) {
          AllOps<double> allOps;
          allOps.M = subM[iSub] = decDomain->getSubDomain(iSub)->constructDBSparseMatrix<double>();
          decDomain->getSubDomain(iSub)->makeSparseOps<double>(allOps, 0, 0, 0);
        }
        SubDOp M(decDomain->getNumSub(), subM);

        // taking advantage of symmetry of M and computing M*Phi_i instead of transpose(Phi_i)*M
        for(int i = 0; i<numRBM; ++i)
          M.mult((*modesRB)[i], tPhiM[i]);
        for(int i = 0; i<numFlex; ++i)
          M.mult((*modesFl)[i], tPhiM[numRBM+i]);
      }
      else {
        for(int i = 0; i<numRBM; ++i) tPhiM[i] = (*modesRB)[i];
        for(int i = 0; i<numFlex; ++i) tPhiM[numRBM+i] = (*modesFl)[i];
      }

      MultiDomainOp mdop(&MultiDomainOp::getInitState, decDomain->getAllSubDomains(),
                         fullDsp, fullVel, fullAcc, fullPrevVel);
      threadManager->execParal(decDomain->getNumSub(), &mdop);

      if(domain->numInitVelocity() > 0) {
        filePrint(stderr, " ... Compute initial velocity in generalized coordinate system ... \n");
        for(int j = 0; j < vel.size(); ++j)
          vel[j] += tPhiM[j]*(*fullVel);
      }
      if((domain->numInitDisp() > 0 || domain->numInitDisp6() > 0) && sinfo.zeroInitialDisp == 0) {
        filePrint(stderr, " ... Compute initial displacement in generalized coordinate system ... \n");
        for(int j = 0; j < dsp.size(); ++j)
          dsp[j] += tPhiM[j]*(*fullDsp);
      }
    }

    // TODO consider the case when both idisp and idisp6 are present
  }
}

//------------------------------------------------------------------------------

void MDModalBase::outputModal(SysState<Vector>& state, Vector& extF, int tIndex, ModalOps &ops){
/*PRE:
 POST:
*/
  if(decDomain->getCommunicator()->cpuNum() != 0) return;
  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();
  SolverInfo &sinfo = domain->solInfo();
  int i, w, p, iMode;
  double time = tIndex * domain->solInfo().getTimeStep();
  Vector &mDsp = state.getDisp();
  Vector &mVel = state.getVeloc();
  Vector &mVel_p = state.getPrevVeloc();

  if (sinfo.nRestart > 0) {
    domain->writeRestartFile(time, tIndex, mDsp, mVel, mVel_p);
  }

  for(i = 0; i < numOutInfo; ++i){
    if((oinfo[i].interval != 0) && (tIndex % oinfo[i].interval == 0)){

      w = oinfo[i].width;
      p = oinfo[i].precision;
      ModalParams &modalParams = domain->solInfo().readInModes[domain->solInfo().modal_id.front()];

      switch(oinfo[i].type){

        case OutputInfo::ModalDsp:
//          filePrint(oinfo[i].filptr, "# modal displacements\n");
          filePrint(oinfo[i].filptr, " % *.*E ", w, p, time);
          for(iMode = 0; iMode < mDsp.size(); ++iMode){
            filePrint(oinfo[i].filptr, "  % *.*E", w, p, mDsp[iMode]);
          }
          filePrint(oinfo[i].filptr, "\n");
          fflush(oinfo[i].filptr);
          break;

        case OutputInfo::ModalExF:
//          filePrint(oinfo[i].filptr, "# modal external forces\n");
          filePrint(oinfo[i].filptr, " % *.*E ", w, p, time);
          for(iMode = 0; iMode < extF.size(); ++iMode){
            filePrint(oinfo[i].filptr, "  % *.*E", w, p, extF[iMode]);
          }
          filePrint(oinfo[i].filptr, "\n");
          fflush(oinfo[i].filptr);
          break;

        case OutputInfo::ModalMatrices:
        case OutputInfo::ModalMass:
          if(time == 0) {
            switch(modalParams.type) {
              case ModalParams::Undefined :
              case ModalParams::Eigen :
              case ModalParams::Mnorm : {
                filePrint(oinfo[i].filptr, "* Mr: Diagonal reduced-order mass matrix\n");
                filePrint(oinfo[i].filptr, "%s\n", (modalParams.type == ModalParams::Mnorm) ? "mnorm" : "eigen");
                filePrint(oinfo[i].filptr, "%d\n", numFlex);
                for(iMode = 0; iMode < numFlex; ++iMode) {
                  filePrint(oinfo[i].filptr, "% *.*E", w, p, ops.M->diag(iMode));
                }
                filePrint(oinfo[i].filptr, "\n");
              } break;
              case ModalParams::Inorm : {
                filePrint(oinfo[i].filptr, "* Mr: Full reduced-order mass matrix\n");
                filePrint(oinfo[i].filptr, "%s\n", "inorm");
                filePrint(oinfo[i].filptr, "%d %d\n", numFlex, numFlex);
                for(iMode = 0; iMode < numFlex; ++iMode) {
                  for(int jMode = 0; jMode < numFlex; ++jMode) {
                    filePrint(oinfo[i].filptr, "% *.*E", w, p, (*ops.M)[iMode+numModes*jMode]);
                  }
                  filePrint(oinfo[i].filptr, "\n");
                }
              } break;
            }
          }
          if(oinfo[i].type != OutputInfo::ModalMatrices) break;

        case OutputInfo::ModalStiffness:
          if(time == 0) {
            switch(modalParams.type) {
              case ModalParams::Undefined :
              case ModalParams::Eigen : {
                filePrint(oinfo[i].filptr, "* Mr: Diagonal reduced-order stiffness matrix\n");
                filePrint(oinfo[i].filptr, "%s\n", "eigen");
                filePrint(oinfo[i].filptr, "%d\n", numFlex);
                for(iMode = 0; iMode < numFlex; ++iMode) {
                  filePrint(oinfo[i].filptr, "% *.*E", w, p, ops.K->diag(iMode));
                }
                filePrint(oinfo[i].filptr, "\n");
              } break;
              case ModalParams::Inorm : {
                filePrint(oinfo[i].filptr, "* Mr: Full reduced-order stiffness matrix\n");
                filePrint(oinfo[i].filptr, "%s\n", (modalParams.type == ModalParams::Mnorm) ? "mnorm" : "inorm");
                filePrint(oinfo[i].filptr, "%d %d\n", numFlex, numFlex);
                for(iMode = 0; iMode < numFlex; ++iMode) {
                  for(int jMode = 0; jMode < numFlex; ++jMode) {
                    filePrint(oinfo[i].filptr, "% *.*E", w, p, (*ops.K)[iMode+numModes*jMode]);
                  }
                  filePrint(oinfo[i].filptr, "\n");
                }
              } break;
            }
          }
          if(oinfo[i].type != OutputInfo::ModalMatrices) break;

        case OutputInfo::ModalDamping:
          if(time == 0) {
            switch(modalParams.type) {
              case ModalParams::Undefined :
              case ModalParams::Eigen : {
                filePrint(oinfo[i].filptr, "* Mr: Diagonal reduced-order damping matrix\n");
                filePrint(oinfo[i].filptr, "%s\n", "eigen");
                filePrint(oinfo[i].filptr, "%d\n", numFlex);
                for(iMode = 0; iMode < numFlex; ++iMode) {
                  filePrint(oinfo[i].filptr, "% *.*E", w, p, (ops.C) ? ops.C->diag(iMode) : 0);
                }
                filePrint(oinfo[i].filptr, "\n");
              } break;
              case ModalParams::Inorm : {
                filePrint(oinfo[i].filptr, "* Mr: Full reduced-order damping matrix\n");
                filePrint(oinfo[i].filptr, "%s\n", (modalParams.type == ModalParams::Mnorm) ? "mnorm" : "inorm");
                filePrint(oinfo[i].filptr, "%d %d\n", numFlex, numFlex);
                for(iMode = 0; iMode < numFlex; ++iMode) {
                  for(int jMode = 0; jMode < numFlex; ++jMode) {
                    filePrint(oinfo[i].filptr, "% *.*E", w, p, (ops.C) ? (*ops.C)[iMode+numModes*jMode] : 0);
                  }
                  filePrint(oinfo[i].filptr, "\n");
                }
              } break;
            }
          }
          break;

          case OutputInfo::ModalDynamicMatrix:
          if(time == 0) {
            filePrint(oinfo[i].filptr, "%d %d\n", numFlex, numFlex);
            for(iMode = 0; iMode < numFlex; ++iMode){
              for(int jMode = 0; jMode < numFlex; ++jMode){
                switch(modalParams.type) {
                  case ModalParams::Undefined :
                  case ModalParams::Eigen : {
                    filePrint(oinfo[i].filptr, "% *.*E", w, p, (iMode==jMode) ? ops.dynMat->diag(iMode) : 0.);
                  } break;
                  default: {
                    filePrint(oinfo[i].filptr, "% *.*E", w, p, (*ops.dynMat)[iMode+numModes*jMode]);
                  } break;
                }
              }
              filePrint(oinfo[i].filptr, "\n");
            }
          }
          break;

        default:
          break;
      }
    }
  }
}
