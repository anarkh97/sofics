#include <Problems.d/ModalBase.h>
#include <Driver.d/Domain.h>
#include <Solvers.d/Rbm.h>
#include <Utils.d/ModeData.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/SysState.h>
#include <Utils.d/DistHelper.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/EiSparseMatrix.h>

extern ModeData modeData;

void DiagonalMatrix::setDiag(double val){
/*PRE: *d is allocated to double[neq]
 POST: set all entires of d to val
*/
  for(int i = 0; i < neq; ++i)
    d[i] = val;
}

//------------------------------------------------------------------------------

void DiagonalMatrix::mult(Vector &v, Vector &Av){
/*PRE: the data in Av have been allocated
 POST: return in Av, the matrix-vector product (*this)*v
*/
  for(int i = 0; i < neq; ++i)
    Av[i] = d[i] * v[i];
}

//------------------------------------------------------------------------------

void DiagonalMatrix::invertDiag(){
/*PRE: d is populated, preferably with non-zero entries
 POST: inversion of the diagonal done in place, ie entries in d are over-written
*/
  for(int i = 0; i < neq; ++i)
    d[i] = 1.0 / d[i];
}

//------------------------------------------------------------------------------

void DiagonalMatrix::reSolve(Vector &rhs) {
/*PRE: d is populated with the inverse of A
 POST: returns solution to Ax = rhs; rhs is over-written with x
*/
  for(int i = 0; i < neq; ++i)
    rhs[i] *= d[i];
}

#ifdef USE_EIGEN3
void DenseMatrix::mult(Vector &v, Vector &Av) {

  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > vMap(v.data(), v.size());
  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > AvMap(Av.data(), Av.size());

  AvMap = denseMat*vMap; 

}

void DenseMatrix::invertDiag() {
  // reusing invertDiag to avoid proliferation of virtual functions
  // compute cholesky factorization for use in reSolve function
  llt_.compute(denseMat);
}

void DenseMatrix::reSolve(Vector &rhs) {

  Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, 1> > rhsMap(rhs.data(), rhs.size());
  llt_.solveInPlace(rhsMap);
}
#endif

//------------------------------------------------------------------------------
//******************************************************************************
//------------------------------------------------------------------------------

void ModalBase::preProcessBase(){
/*PRE: none
 POST: initialize and populate various data members
*/
  numModes = numRBM = numFlex = 0;
  domain->preProcessing();

  const int numDofs = domain->numdof();

  int *bc = (int *) dbg_alloca(sizeof(int)*numDofs);
  bcx = new double[numDofs];
  vcx = new double[numDofs];

  int i;
  for(i = 0; i < numDofs; ++i) vcx[i] = 0.0;

  domain->make_bc(bc, bcx);
  domain->make_constrainedDSA(1);
  domain->makeAllDOFs();

  numModes = modeData.numModes;  // this is assigment is temporary; it is
                                 //   overrulled in populateFlexModes
  fullTmpF.setData(new double[numDofs], numDofs);
  fullTmpGrav.setData(new double[numDofs], numDofs);
  fullAeroF.setData(new double[numDofs], numDofs);
  fullDsp.setData(new double[numDofs], numDofs);
  fullVel.setData(new double[numDofs], numDofs);
  fullAcc.setData(new double[numDofs], numDofs);
  fullPrevVel.setData(new double[numDofs], numDofs);

  fullTmpF.zeroAll();
  fullTmpGrav.zeroAll();
  fullAeroF.zeroAll();
  fullDsp.zeroAll();
  fullVel.zeroAll();
  fullAcc.zeroAll();
  fullPrevVel.zeroAll();

  prevFrc = new PrevFrc(numDofs);
  prevFrcBackup = new PrevFrc(numDofs);

  numConstr = domain->nDirichlet();
  cDofIdx = new int[numConstr];
  for(i = 0; i < numConstr; ++i){
    cDofIdx[i] = domain->getCDSA()->locate(domain->getDBC()[i].nnum,
      1 << domain->getDBC()[i].dofnum);
  }
}

//------------------------------------------------------------------------------

void ModalBase::populateRBModes(){

  // compute the rigid body modes
  Rbm *rbm = domain->constructRbm();
  numRBM   = rbm->numRBM();
  modesRB  = new Vector[numRBM];

  rbm->getRBMs(modesRB);
  rbm->getxyzRot(0, cg);      // cg temporarily stores the point about
                              //   which rotation RBMs were calculated
}

//------------------------------------------------------------------------------

void ModalBase::populateFlexModes(double scale, bool readAll){
/*PRE: preProcessBase has been called
       scale has default value of 1.0; readAll has default value of 0
       the modes in modeData are : mode^T.M.mode = identity if type
       is eigen or Mnorm, dense otherwise
 POST: populate modesFl with data from modeData multiplied by scale
       populate freqs with the circular frequency of each flexible mode
 NOTE: if readAll, then also include in modesFL those modes with zero frequency
*/
  const int numNodes   = modeData.numNodes;
  const int numDofs  = domain->numdof();

  // count the number of flexible modes
  numFlex = 0;
  int iMode;
  for(iMode = 0; iMode < numModes; ++iMode){
    if( readAll || (modeData.frequencies[iMode] != 0.0) )
      ++numFlex;
  }

  // allocate space for modal basis and eigenfrequencies
  modesFl = new Vector[numFlex];
  freqs   = new double[numFlex];

  // load scaled modal data
  int iModeFl = 0, iNode, dof;
  for(iMode = 0; iMode < numModes; ++iMode){
    if( readAll || (modeData.frequencies[iMode] != 0.0) ){
      modesFl[iModeFl].setData(new double[numDofs], numDofs);
      modesFl[iModeFl].zeroAll();
      freqs[iModeFl] = modeData.frequencies[iMode] * 8 * atan(1.);

      for(iNode = 0; iNode < numNodes; ++iNode){

        dof = domain->getCDSA()->locate(iNode, DofSet::Xdisp);
        if(dof >= 0) modesFl[iModeFl][dof] = scale*modeData.modes[iMode][iNode][0];

        dof = domain->getCDSA()->locate(iNode, DofSet::Ydisp);
        if(dof >= 0) modesFl[iModeFl][dof] = scale*modeData.modes[iMode][iNode][1];

        dof = domain->getCDSA()->locate(iNode, DofSet::Zdisp);
        if(dof >= 0) modesFl[iModeFl][dof] = scale*modeData.modes[iMode][iNode][2];

        dof = domain->getCDSA()->locate(iNode, DofSet::Xrot);
        if(dof >= 0) modesFl[iModeFl][dof] = scale*modeData.modes[iMode][iNode][3];

        dof = domain->getCDSA()->locate(iNode, DofSet::Yrot);
        if(dof >= 0) modesFl[iModeFl][dof] = scale*modeData.modes[iMode][iNode][4];

        dof = domain->getCDSA()->locate(iNode, DofSet::Zrot);
        if(dof >= 0) modesFl[iModeFl][dof] = scale*modeData.modes[iMode][iNode][5];
      }
      ++iModeFl;
    }
  }
  
}

//------------------------------------------------------------------------------

void ModalBase::initStateBase(Vector& dsp, Vector& vel,
  Vector& acc, Vector& vel_p, int idxOffset){

  dsp.zero();
  vel.zero();
  acc.zero();
  vel_p.zero();

  ControlInfo *cinfo = geoSource->getCheckFileInfo();
  SolverInfo &sinfo = domain->solInfo();

  if (cinfo->lastRestartFile) {
     fprintf(stderr, " ... Restarting From a Previous Run ...\n");
     int fn = open(cinfo->lastRestartFile,O_RDONLY );
     if(fn >= 0) {
       int vsize, restartTIndex;
       double restartT;
       int readSize = read(fn, &vsize, sizeof(int));
       if(vsize != dsp.size() || readSize != sizeof(int))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
       readSize = read(fn, &restartTIndex, sizeof(int));
       if(readSize != sizeof(int))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
       if(strcmp(cinfo->FlagRST,"new") == 0)
         sinfo.initialTimeIndex = 0;
       else
         sinfo.initialTimeIndex = restartTIndex;
       readSize = read(fn, &restartT, sizeof(double));
       if(readSize != sizeof(double))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");
       if(strcmp(cinfo->FlagRST,"new") == 0)
         sinfo.initialTime = 0.0;
       else
         sinfo.initialTime = restartT;
       readSize = read(fn, dsp.data(), vsize*sizeof(double));
       if(int(readSize) != int(vsize*sizeof(double)))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");

       readSize = read(fn, vel.data(), vsize*sizeof(double));
       if(int(readSize) != int(vsize*sizeof(double)))
         fprintf(stderr," *** ERROR: Inconsistent restart file\n");

       readSize = read(fn, vel_p.data(), vsize*sizeof(double));
       if(int(readSize) != int(vsize*sizeof(double))) {
         fprintf(stderr," *** WARNING: Older version of restart file"
                        " -- Missing velocity field is set to zero\n");
         vel_p.zero();
       }

       close(fn);
     } 
     else {
       fprintf(stderr, " *** ERROR: Restart file could not be opened\n");
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

      double **tPhiM = new double*[numFlex+numRBM];
      if(massNormalizedBasis) {
        for(int i = 0; i < numFlex+numRBM; ++i)
          tPhiM[i] = new double[domain->numdof()];

        // construct and assemble full mass matrix
        AllOps<double> allOps;
        allOps.M = domain->constructDBSparseMatrix<double>();
        domain->makeSparseOps(allOps, 0, 0, 0, (SparseMatrix *) NULL, (FullSquareMatrix *) NULL, (FullSquareMatrix *) NULL);

        // taking advantage of symmetry of M and computing M*Phi_i instead of transpose(Phi_i)*M
        for(int i = 0 ; i<numRBM; ++i)
          allOps.M->mult(modesRB[i].data(), tPhiM[i]);
        for(int i = 0 ; i<numFlex; ++i)
          allOps.M->mult(modesFl[i].data(), tPhiM[numRBM+i]);
        delete allOps.M;
      }
      else {
        for(int i = 0 ; i<numRBM; ++i)
          tPhiM[i] = modesRB[i].data();
        for(int i = 0 ; i<numFlex; ++i)
          tPhiM[numRBM+i] = modesFl[i].data();
      }

      if(domain->numInitVelocity() > 0) {
        filePrint(stderr, " ... Compute initial velocity in generalized coordinate system ... \n");
        Vector fullVel(domain->numdof(), 0.0);
        for(int j = 0; j < domain->numInitVelocity(); ++j) {
          int k = domain->getCDSA()->locate(domain->getInitVelocity()[j].nnum, 1 << domain->getInitVelocity()[j].dofnum);
          if(k > -1) fullVel[k] = domain->getInitVelocity()[j].val;
        }
        for(int i = 0; i < numConstr; ++i){
          fullVel[cDofIdx[i]] = 0; // just in case an initial velocity is given for a constrained node
        }
        for(int j = 0; j < vel.size(); ++j)
          for(int k = 0; k < fullVel.size(); ++k)
            vel[j] += tPhiM[j][k]*fullVel[k];
      }
      if(sinfo.zeroInitialDisp == 0) {
        if(domain->numInitDisp() > 0 && (domain->numInitDisp6() == 0 || sinfo.gepsFlg == 1)) {
          filePrint(stderr, " ... Compute initial displacement in generalized coordinate system ... \n");
          Vector fullDsp(domain->numdof(), 0.0);
          for(int j = 0; j <  domain->numInitDisp(); ++j) {
            int k = domain->getCDSA()->locate(domain->getInitDisp()[j].nnum, 1 << domain->getInitDisp()[j].dofnum);
            if(k > -1) fullDsp[k] = domain->getInitDisp()[j].val;
          }
          for(int j = 0; j < dsp.size(); ++j)
            for(int k = 0; k < fullDsp.size(); ++k)
              dsp[j] += tPhiM[j][k]*fullDsp[k];
        }
        if(domain->numInitDisp6() > 0 && sinfo.gepsFlg == 0) {
          filePrint(stderr, " ... Compute initial displacement in generalized coordinate system ... \n");
          Vector fullDsp(domain->numdof(), 0.0);
          for(int j = 0; j <  domain->numInitDisp6(); ++j) {
            int k = domain->getCDSA()->locate(domain->getInitDisp6()[j].nnum, 1 << domain->getInitDisp6()[j].dofnum);
            if(k > -1) fullDsp[k] = domain->getInitDisp6()[j].val;
          }
          for(int j = 0; j < dsp.size(); ++j)
            for(int k = 0; k < fullDsp.size(); ++k)
              dsp[j] += tPhiM[j][k]*fullDsp[k];
        }
      }

      if(massNormalizedBasis) {
        for(int i = 0; i < numFlex+numRBM; ++i)
          delete [] tPhiM[i];
      }
      delete [] tPhiM;
    }

    // TODO consider the case when both idisp and idisp6 are present
  }
}

//------------------------------------------------------------------------------

void ModalBase::outputModal(SysState<Vector>& state, Vector& extF, int tIndex, ModalOps &ops){
/*PRE:
 POST:
*/
  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();
  SolverInfo &sinfo = domain->solInfo();
  int i, w, p, iMode;
  double time = tIndex * domain->solInfo().getTimeStep();
  Vector &mDsp = state.getDisp();
  Vector &mVel = state.getVeloc();
  Vector &mVel_p = state.getPrevVeloc();

  if (sinfo.nRestart > 0)
    domain->writeRestartFile(time, tIndex, mDsp, mVel, mVel_p);

  for(i = 0; i < numOutInfo; ++i){
    if((oinfo[i].interval != 0) && (tIndex % oinfo[i].interval == 0)){

      w = oinfo[i].width;
      p = oinfo[i].precision;
      ModalParams &modalParams = domain->solInfo().readInModes[domain->solInfo().modal_id.front()];

      switch(oinfo[i].type){

        case OutputInfo::ModalDsp:
//          fprintf(oinfo[i].filptr, "# modal displacements\n");
          fprintf(oinfo[i].filptr, " % *.*E ", w, p, time);
          for(iMode = 0; iMode < mDsp.size(); ++iMode){
            fprintf(oinfo[i].filptr, "  % *.*E", w, p, mDsp[iMode]);
          }
          fprintf(oinfo[i].filptr, "\n");
          fflush(oinfo[i].filptr);
          break;

        case OutputInfo::ModalExF:
//          fprintf(oinfo[i].filptr, "# modal external forces\n");
          fprintf(oinfo[i].filptr, " % *.*E ", w, p, time);
          for(iMode = 0; iMode < extF.size(); ++iMode){
            fprintf(oinfo[i].filptr, "  % *.*E", w, p, extF[iMode]);
          }
          fprintf(oinfo[i].filptr, "\n");
          fflush(oinfo[i].filptr);
          break;

        case OutputInfo::ModalMatrices:
        case OutputInfo::ModalMass:
          if(time == 0) {
            switch(modalParams.type) {
              case ModalParams::Undefined :
              case ModalParams::Eigen :
              case ModalParams::Mnorm : {
                fprintf(oinfo[i].filptr, "* Mr: Diagonal reduced-order mass matrix\n");
                fprintf(oinfo[i].filptr, "%s\n", (modalParams.type == ModalParams::Mnorm) ? "mnorm" : "eigen");
                fprintf(oinfo[i].filptr, "%d\n", numFlex);
                for(iMode = 0; iMode < numFlex; ++iMode) {
                  fprintf(oinfo[i].filptr, "% *.*E", w, p, ops.M->diag(iMode));
                }
                fprintf(oinfo[i].filptr, "\n");
              } break;
              case ModalParams::Inorm : {
                fprintf(oinfo[i].filptr, "* Mr: Full reduced-order mass matrix\n");
                fprintf(oinfo[i].filptr, "%s\n", "inorm");
                fprintf(oinfo[i].filptr, "%d %d\n", numFlex, numFlex);
                for(iMode = 0; iMode < numFlex; ++iMode) {
                  for(int jMode = 0; jMode < numFlex; ++jMode) {
                    fprintf(oinfo[i].filptr, "% *.*E", w, p, (*ops.M)[iMode+numModes*jMode]);
                  }
                  fprintf(oinfo[i].filptr, "\n");
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
                fprintf(oinfo[i].filptr, "* Mr: Diagonal reduced-order stiffness matrix\n");
                fprintf(oinfo[i].filptr, "%s\n", "eigen");
                fprintf(oinfo[i].filptr, "%d\n", numFlex);
                for(iMode = 0; iMode < numFlex; ++iMode) {
                  fprintf(oinfo[i].filptr, "% *.*E", w, p, ops.K->diag(iMode));
                }
                fprintf(oinfo[i].filptr, "\n");
              } break;
              case ModalParams::Inorm : {
                fprintf(oinfo[i].filptr, "* Mr: Full reduced-order stiffness matrix\n");
                fprintf(oinfo[i].filptr, "%s\n", (modalParams.type == ModalParams::Mnorm) ? "mnorm" : "inorm");
                fprintf(oinfo[i].filptr, "%d %d\n", numFlex, numFlex);
                for(iMode = 0; iMode < numFlex; ++iMode) {
                  for(int jMode = 0; jMode < numFlex; ++jMode) {
                    fprintf(oinfo[i].filptr, "% *.*E", w, p, (*ops.K)[iMode+numModes*jMode]);
                  }
                  fprintf(oinfo[i].filptr, "\n");
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
                fprintf(oinfo[i].filptr, "* Mr: Diagonal reduced-order damping matrix\n");
                fprintf(oinfo[i].filptr, "%s\n", "eigen");
                fprintf(oinfo[i].filptr, "%d\n", numFlex);
                for(iMode = 0; iMode < numFlex; ++iMode) { 
                  fprintf(oinfo[i].filptr, "% *.*E", w, p, (ops.C) ? ops.C->diag(iMode) : 0);
                }
                fprintf(oinfo[i].filptr, "\n");
              } break;
              case ModalParams::Inorm : {
                fprintf(oinfo[i].filptr, "* Mr: Full reduced-order damping matrix\n");
                fprintf(oinfo[i].filptr, "%s\n", (modalParams.type == ModalParams::Mnorm) ? "mnorm" : "inorm");
                fprintf(oinfo[i].filptr, "%d %d\n", numFlex, numFlex);
                for(iMode = 0; iMode < numFlex; ++iMode) {
                  for(int jMode = 0; jMode < numFlex; ++jMode) {
                    fprintf(oinfo[i].filptr, "% *.*E", w, p, (ops.C) ? (*ops.C)[iMode+numModes*jMode] : 0);
                  }
                  fprintf(oinfo[i].filptr, "\n");
                }
              } break;
            }
          }
          break;

          case OutputInfo::ModalDynamicMatrix:
          if(time == 0) {
            fprintf(oinfo[i].filptr, "%d %d\n", numFlex, numFlex);
            for(iMode = 0; iMode < numFlex; ++iMode){
              for(int jMode = 0; jMode < numFlex; ++jMode){
                switch(modalParams.type) {
                  case ModalParams::Undefined :
                  case ModalParams::Eigen : {
                    fprintf(oinfo[i].filptr, "% *.*E", w, p, (iMode==jMode) ? ops.dynMat->diag(iMode) : 0.);
                  } break;
                  default: {
                    fprintf(oinfo[i].filptr, "% *.*E", w, p, (*ops.dynMat)[iMode+numModes*jMode]);
                  } break;
                }
              }
              fprintf(oinfo[i].filptr, "\n");
            }
          }
          break;          

          case OutputInfo::ErrorIndicator:
          if(time == 0) {
#ifdef USE_EIGEN3
            std::ifstream fin("EIGENMODES");
            if(fin.is_open()) {
              // 1. read ROM eigenvalues and eigenvectors from file
              int Nev;
              fin >> Nev;
              if(Nev != numFlex) { fprintf(stderr, " *** Error: Bad data in EIGENMODES file %d %d ***\n", Nev, numFlex); exit(-1); }
              fprintf(stderr," ... Reading %d modes from EIGENMODES file ...\n", Nev);
              Eigen::VectorXd lambda(Nev);
              Eigen::MatrixXd y(Nev, Nev);
              for(int j = 0; j < Nev; ++j) {
                fin >> lambda[j];
                for(int k = 0; k < Nev; ++k) {
                  fin >> y.col(j)[k];
                }
              }

              // 2. construct and assemble full mass and stiffness matrices
              AllOps<double> allOps;
              allOps.M = domain->constructEiSparseMatrix<double, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper>>();
              allOps.K = domain->constructEiSparseMatrix<double, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper>>();
              domain->makeSparseOps(allOps, 0, 0, 0, (SparseMatrix *) NULL, (FullSquareMatrix *) NULL, (FullSquareMatrix *) NULL);
              Eigen::MappedSparseMatrix<double, Eigen::ColMajor, int>& M =
                static_cast<GenEiSparseMatrix<double, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper>>&>(*allOps.M).getEigenSparse();
              Eigen::MappedSparseMatrix<double, Eigen::ColMajor, int>& K =
                static_cast<GenEiSparseMatrix<double, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper>>&>(*allOps.K).getEigenSparse();

              // 3. copy flexural modes into basis V
              Eigen::MatrixXd V(M.rows(), Nev);
              for(int j = 0; j < Nev; ++j)
                V.col(j) = Eigen::Map<Eigen::VectorXd>(modesFl[j].data(), V.rows());

              // 4. precompute the matrix 1-norms and matrix products
              double Mnorm = (M.cwiseAbs() * Eigen::VectorXd::Ones(M.cols())).maxCoeff();
              double Knorm = (K.cwiseAbs() * Eigen::VectorXd::Ones(K.cols())).maxCoeff();
              Eigen::MatrixXd Vy = V*y;
              Eigen::MatrixXd KVy = K*Vy;
              Eigen::MatrixXd MVy = M*Vy;

              // 5. compute the error indicator and write to file
              double num = 0, den = 0;
              for(int j = 0; j < Nev; ++j) {
                num += (KVy.col(j) - lambda[j]*MVy.col(j)).norm();
                den += (Knorm + lambda[j]*Mnorm)*(Vy.col(j)).norm();
              }
              fprintf(oinfo[i].filptr, "% *.*E", w, p, num/den);
            }
            else {
              fprintf(stderr, " *** Error: Failed to open EIGENMODES file ***\n");
              exit(-1);
            }
#endif
          }
          break;

        default:
          break;
      }
    }
  }
}
