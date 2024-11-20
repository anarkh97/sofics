#include <Driver.d/DecDomain.h>
#include <algorithm>
#include <Math.d/BLAS.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Utils.d/DistHelper.h>

extern int verboseFlag;

template<class Scalar>
MultiDomainRbm<Scalar>::MultiDomainRbm(GenDecDomain<Scalar> *_decDomain, double _tolgrb)
  : decDomain(_decDomain), tolgrb(_tolgrb), numGtGsing(0)
{
  computeRbms();
}

template<class Scalar>
void
MultiDomainRbm<Scalar>::computeRbms()
{
  int i, nGroups, nGroups1;
  int *groups = 0, *ngrbmGr = 0;
  const Connectivity *groupToSub;
  Connectivity subToGroup;
  auto bodyToSub = decDomain->getGroupToSub().get();
  if(!bodyToSub) {
    int nsubGl = decDomain->getGlobalNumSub();
    int *pointer = new int[nsubGl+1];
    int *target = new int[nsubGl];
    for(int i=0; i<nsubGl; ++i) { pointer[i] = i; target[i] = i; }
    pointer[nsubGl] = nsubGl;
    bodyToSub = new Connectivity(nsubGl, pointer, target);
  }
  Connectivity subToBody = bodyToSub->reverse();
  int nsub = decDomain->getNumSub();
  GenSubDomain<Scalar> **sd = decDomain->getAllSubDomains();
  FSCommunicator *com = decDomain->getCommunicator();
  int numCPUs = com->size();
  int myCPU = com->cpuNum();
  Connectivity *mpcToSub_primal = decDomain->getMpcToSub_primal();
  int glNumMpc_primal = (mpcToSub_primal) ? mpcToSub_primal->csize() : 0;

  // PJSA: start new code **************************************************

  int nBodies = (bodyToSub) ? bodyToSub->csize() : 0;
  bool mbgflag = false;
  Connectivity *groupToBody = 0;

  if(glNumMpc_primal > 0) {

    // Find out which bodies are connected together by mpcs 
    // a collection of inter-connected bodies is referred to as a group
    Connectivity subToMpc = mpcToSub_primal->reverse();
    Connectivity bodyToMpc = bodyToSub->transcon(subToMpc);
    Connectivity mpcToBody = bodyToMpc.reverse();
    Connectivity bodyToBodyTmp = bodyToMpc.transcon(mpcToBody);
    Connectivity bodyToBody = bodyToBodyTmp.modify();
    compStruct renumber = bodyToBody.renumByComponent(1);  // 1 = sloan renumbering
    nGroups = renumber.numComp;
    if(nGroups < nBodies) { // at least one multi-body group exists
      mbgflag = true;
      // make groupToBody connectivity
      renumber.order = new int[nBodies];
      for(i = 0; i < nBodies; ++i)
        renumber.order[renumber.renum[i]] = i;
      int *pointer = new int[nGroups + 1];
      pointer[0] = 0;
      for(i = 0; i < nGroups; ++i) {
        int nbod = renumber.xcomp[i + 1] - renumber.xcomp[i];
        pointer[i + 1] = pointer[i] + nbod;
      }
      groupToBody = new Connectivity(nGroups, pointer, renumber.order);
      groupToSub = groupToBody->transcon(bodyToSub);
      subToGroup = groupToSub->reverse();
    }
    if(renumber.xcomp) delete [] renumber.xcomp;
    if(renumber.renum) delete [] renumber.renum;
  }
  if(!mbgflag) { // one body per group
    groupToSub = bodyToSub;
    subToGroup = subToBody;
    nGroups = nBodies;
  }

  // tell each subDomain what group it is in and find groups on each processor
  paralApply(nsub, sd, &BaseSub::setGroup, &subToGroup);
  groups = new int[nGroups];  // groups represented on this processor
#ifdef  DISTRIBUTED
  if(sd && nsub > 0) {
    groups[0] = subToGroup[sd[0]->subNum()][0];
    int n = 1;
    for(int i = 1; i < nsub; ++i) {
      int group = subToGroup[sd[i]->subNum()][0];
      int j;
      for(j = 0; j < n; ++j) if(group == groups[j]) j = n+1;
      if(j == n) groups[n++] = group;
    }
    nGroups1 = n;  // number of groups represented on this processor
  }
  else nGroups1 = 0;
#else
  for(int i = 0; i < nGroups; ++i) groups[i] = i; 
  nGroups1 = nGroups;  
#endif
  ngrbmGr = new int[nGroups];
  for(int i = 0; i < nGroups; ++i) ngrbmGr[i] = 0;
  if(verboseFlag) filePrint(stderr, " ... Number of bodies = %3d         ...\n", nGroups);
  
  int ngrbm = 0;

/* XXX
  if((this->fetiInfo->corners == FetiInfo::noCorners) && (glNumMpc_primal > 0)) {
    // subdomain ZEMs need to be projected or eliminated
    // I think adding the dofs of primal mpcs to the "c" dofs would resolve this issue
    filePrint(stderr, " *** ERROR: mpc_type 2 is not supported with no corners *** \n"); 
    exit(-1);
  }
*/
  if(verboseFlag) filePrint(stderr, " ... Computing Multi-domain GRBMs   ...\n");
  // calculate the centroid of each group
  double *centroid = new double[nGroups*3];   // pseudo centroid of each group
  double *nNodes = new double[nGroups];       // number of nodes in each group;
  for(int i = 0; i < nGroups; ++i) {
    nNodes[i] = 0.0;
    for(int j = 0; j < 3; ++j) centroid[3*i+j] = 0.0;
  }
  for(int i = 0; i < nsub; ++i) sd[i]->addNodeXYZ(centroid, nNodes);  // groups could be done in parallel
#ifdef DISTRIBUTED
  com->globalSum(nGroups, nNodes);
  com->globalSum(nGroups*3, centroid);
#endif
  for(int i = 0; i < nGroups; ++i) {
    if(nNodes[i] > 0) 
      for(int j = 0; j < 3; ++j) 
        centroid[3*i+j] = centroid[3*i+j]/nNodes[i];
  }
  // note: this centroid calculation assumes nodes are equally spaced, and also
  // interface nodes from a group split over more than one interface are used more than once.
  // however, the objective is only to find a point inside the group to be used as 
  // a reference for the geometric rbm calculation.  it is not necessary to use the exact
  // geometric centroid.
  delete [] nNodes;
 
  // make Zstar and R matrices for each subdomain
  paralApply(nsub, sd, &GenSubDomain<Scalar>::makeZstarAndR, centroid);
  delete [] centroid;
 
  Connectivity *groupToMpc = 0;
  if(glNumMpc_primal > 0) {
    Connectivity *subToMpc = mpcToSub_primal->alloc_reverse();
    groupToMpc = groupToSub->transcon(subToMpc);
    paralApply(nsub, sd, &BaseSub::makeLocalToGroupMPC, groupToMpc);
    delete subToMpc;
  }

  // assemble global Zstar matrix for each body
  std::vector<FullM> globalZstar;
  globalZstar.reserve(nGroups);
  int *zRow = new int[nGroups];
  int *zRowDim = new int[nGroups];
  int *zColDim = new int[nGroups];
  int *zColOffset = new int[nBodies];
  int zColDim1 = (sd && nsub > 0) ? sd[0]->zColDim() : 0;  // (6 for 3D, 3 for 2D)
#ifdef DISTRIBUTED
  zColDim1 = com->globalMax(zColDim1);  // enforce it to be the same
#endif
  for(int i = 0; i < nGroups; ++i) {
    zRowDim[i] = 0; 
    if(mbgflag) {
      int cOff = 0;
      for(int j = 0; j < groupToBody->num(i); ++j) {
        int body = (*groupToBody)[i][j];
        zColOffset[body] = cOff;
        cOff += zColDim1;
      }
      zColDim[i] = cOff;
    }
    else {
      zColOffset[i] = 0;
      zColDim[i] = zColDim1;
    }
  }
  for(int iSub = 0; iSub < nsub; ++iSub) {
    int subGroup = subToGroup[sd[iSub]->subNum()][0];
    zRowDim[subGroup] += sd[iSub]->zRowDim();
  }
  if(glNumMpc_primal > 0) {
    for(int i = 0; i < nGroups; ++i) zRowDim[i] += groupToMpc->num(i);
  }
  if(groupToBody) delete groupToBody;

#ifdef DISTRIBUTED
  int *zRowOffset = new int[numCPUs*nGroups];
  for(int i = 0; i < numCPUs*nGroups; ++i) zRowOffset[i] = 0;
  for(int i = 0; i < nGroups1; ++i) {
    int iGroup = groups[i];
    for(int j = myCPU+1; j < numCPUs; ++j) zRowOffset[iGroup*numCPUs +j] = zRowDim[iGroup];
  }
  com->globalSum(nGroups, zRowDim);
  com->globalSum(numCPUs*nGroups, zRowOffset);
  for(int i = 0; i < nGroups; ++i) zRow[i] = zRowOffset[i*numCPUs + myCPU];
  delete [] zRowOffset;
#else
  for(int i = 0; i < nGroups; ++i) zRow[i] = 0;
#endif
  for(int i = 0; i < nGroups; ++i) {
    globalZstar.emplace_back(zRowDim[i], zColDim[i]);
    globalZstar[i].zero();
  }
  // could do this in parallel (by groups)
  for(int iSub = 0; iSub < nsub; ++iSub) {
    int subBody = subToBody[sd[iSub]->subNum()][0];
    int subGroup = subToGroup[sd[iSub]->subNum()][0];
    if(sd[iSub]->zRowDim() > 0) 
      sd[iSub]->addSPCsToGlobalZstar(globalZstar[subGroup], zRow[subGroup], zColOffset[subBody]);
    if(sd[iSub]->numMPCs_primal() > 0) {
      int startRow = zRowDim[subGroup] - groupToMpc->num(subGroup);
      sd[iSub]->addMPCsToGlobalZstar(globalZstar[subGroup], startRow, zColOffset[subBody], zColDim1);
    }
  }
  if(glNumMpc_primal > 0) execParal(nsub, this, &MultiDomainRbm<Scalar>::setBodyRBMoffset, zColOffset, &subToBody);
  delete [] zColOffset;
  if(groupToMpc) delete groupToMpc;
 
  int *groupProc = new int[nGroups];
#ifdef DISTRIBUTED
  for(int i = 0; i < nGroups; ++i) {
    com->globalSum(zRowDim[i]*zColDim[i], globalZstar[i].data());
    groupProc[i] = -1;
  }
  for(int i = 0; i < nGroups1; ++i) groupProc[groups[i]] = myCPU;
  for(int i = 0; i < nGroups; ++i) groupProc[i] = com->globalMax(groupProc[i]);
#else
  for(int i = 0; i < nGroups; ++i) groupProc[i] = myCPU;
#endif
 
  // now do svd on globalZstar for each group to get globalQ for each group
  ngrbm = 0;  // total of all groups
  FullM  **Qtranspose;
  Qtranspose = new FullM * [nGroups];
  for(int i = 0; i < nGroups1; ++i) {
    int iGroup = groups[i];
    int ncol = zColDim[iGroup];  
    int nrow = zRowDim[iGroup];
    FullM U(ncol,ncol); U.zero(); 
    int rank = 0;
    singularValueDecomposition(globalZstar[iGroup], U, ncol, nrow, rank, tolgrb);
    int ngrbmGrTmp = ncol - rank;
    globalZstar[iGroup].clean_up();
    if(groupProc[iGroup] == myCPU) {
      ngrbmGr[iGroup] = ngrbmGrTmp;
      ngrbm += ngrbmGr[iGroup];
      if(verboseFlag) fprintf(stderr, " ... Number of GRBMs for Body %-2d= %-2d...\n", iGroup, ngrbmGrTmp);
    }
    Qtranspose[iGroup] = new FullM(U, ngrbmGrTmp, rank, ncol, 0);
  }
#ifdef DISTRIBUTED
  int ngrbms = com->globalSum(ngrbm);  // total number of rigid body modes for all processes
#else
  int ngrbms = ngrbm;
#endif
  if(verboseFlag) filePrint(stderr, " ... Total Number of GRBMs = %-4d   ...\n", ngrbms);
 
  delete [] groupProc;
  delete [] zRow;
  delete [] zRowDim;
  delete [] zColDim;
 
  // make local Rstar
  paralApply(nsub, sd, &GenSubDomain<Scalar>::makeLocalRstar, Qtranspose);
  for(int i = 0; i < nGroups1; ++i) delete Qtranspose[groups[i]];
  delete [] Qtranspose;
 
  if(ngrbms) {
#ifdef DISTRIBUTED
    com->globalSum(nGroups, ngrbmGr);
#endif
    paralApply(nsub, sd, &BaseSub::setNumGroupRBM, ngrbmGr);
    paralApply(nsub, sd, &FetiSub<Scalar>::deleteLocalRBMs);
  }

  // end new code ********************************************************

  if(ngrbms) {

    // 1. make local G = B * R
    paralApply(nsub, sd, &GenSubDomain<Scalar>::makeG);

    // 2. exchange G's between neighboring subdomains
    Connectivity *cpuToSub = decDomain->getCpuToSub().get();
    FSCommPattern<int> *sPat = new FSCommPattern<int>(com, cpuToSub, myCPU, FSCommPattern<int>::CopyOnSend);
    for(int i = 0; i < nsub; ++i) sd[i]->setMpcNeighbCommSize(sPat, 2);
    sPat->finalize();
    paralApply(nsub, sd, &BaseSub::sendNeighbGrbmInfo, sPat);
    sPat->exchange();
    paralApply(nsub, sd, &BaseSub::receiveNeighbGrbmInfo, sPat);
    delete sPat;
    FSCommPattern<Scalar> *gPat = new FSCommPattern<Scalar>(com, cpuToSub, myCPU,
                                       FSCommPattern<Scalar>::CopyOnSend, FSCommPattern<Scalar>::NonSym);
    for(int i = 0; i < nsub; ++i) sd[i]->setGCommSize(gPat);
    gPat->finalize();
    paralApply(nsub, sd, &GenSubDomain<Scalar>::sendG, gPat);
    gPat->exchange();
    paralApply(nsub, sd, &GenSubDomain<Scalar>::receiveG, gPat);
    delete gPat;

    // 3. build coarse connectivity and equation numberer
    Connectivity coarseConnectGtG;
    const Connectivity &mpcToSub = decDomain->getMpcToSub();
    if(mpcToSub.csize() > 0) {
      Connectivity mpcToBody = mpcToSub.transcon(subToGroup);
      Connectivity bodyToMpc = mpcToBody.reverse();
      Connectivity bodyToBody_mpc = bodyToMpc.transcon(mpcToBody);
      coarseConnectGtG = bodyToBody_mpc.modify();
    }
    else {
      coarseConnectGtG = bodyToSub->transcon(subToBody);
    }

    SimpleNumberer *eqNumsGtG = new SimpleNumberer(nGroups);
    for(int i = 0; i < nGroups; ++i) eqNumsGtG->setWeight(i, ngrbmGr[i]);
    eqNumsGtG->makeOffset();

    // 4. create, assemble and factorize GtG
    GenSkyMatrix<Scalar> *GtGtilda = new GenSkyMatrix<Scalar>(&coarseConnectGtG, eqNumsGtG, tolgrb);
    execParal(nGroups1, this, &MultiDomainRbm<Scalar>::assembleGtG, groups, groupToSub, 
                static_cast<GenSparseMatrix<Scalar> *>(GtGtilda));
#ifdef DISTRIBUTED
    GtGtilda->unify(com);
#endif
    GtGtilda->setPrintNullity(false);
    GtGtilda->parallelFactor();

    // 5. check for singularities in GtG (representing global RBMs)
    numGtGsing = GtGtilda->numRBM();
    if(numGtGsing > 0) {
      // get null space of GtGtilda
      Scalar *zem = new Scalar[numGtGsing*ngrbms];
      GtGtilda->getNullSpace(zem);
      GenFullM<Scalar> X(zem, ngrbms, numGtGsing, 1);
      // build global RBMs
      paralApply(nsub, sd, &GenSubDomain<Scalar>::buildGlobalRBMs, X, (Connectivity*) NULL);
      delete [] zem;
    }

    delete eqNumsGtG;
    delete GtGtilda;
  }
  if(!decDomain->getGroupToSub()) delete bodyToSub;
  delete [] groups;
  delete [] ngrbmGr;
}

template<class Scalar>
void
MultiDomainRbm<Scalar>::setBodyRBMoffset(int iSub, int *zColOffset, Connectivity *subToBody)
{
  GenSubDomain<Scalar> *sd = decDomain->getSubDomain(iSub);
  int subBody = (*subToBody)[sd->subNum()][0];
  sd->setBodyRBMoffset(zColOffset[subBody]);
}

template<class Scalar>
void
MultiDomainRbm<Scalar>::assembleGtG(int iGroup, const int *groups, const Connectivity *groupToSub,
                                    GenSparseMatrix<Scalar> *GtGsparse)
{
  // assembles groups in parallel, subdomains with same group sequentially 
  // threadsafe implementation - avoids simultaneous writing to same memory
  // note: the distributed version will work for shared memory too, but the 
  // alternative code is a bit more efficient
  int nsub = decDomain->getNumSub();
  GenSubDomain<Scalar> **sd = decDomain->getAllSubDomains();
#ifdef DISTRIBUTED
  for(int i = 0; i < nsub; ++i) {
    if(sd[i]->group == groups[iGroup])
      sd[i]->assembleGtGsolver(GtGsparse);
  }
#else
  auto grsubs = (*groupToSub)[iGroup];
  for(int i = 0; i < groupToSub->num(iGroup); ++i) {
    int iSub = grsubs[i];
    sd[iSub]->assembleGtGsolver(GtGsparse);
  }
#endif
}

template<class Scalar>
int
MultiDomainRbm<Scalar>::numRBM()
{
  return numGtGsing;
}

template<class Scalar>
void
MultiDomainRbm<Scalar>::getRBMs(GenDistrVectorSet<Scalar>& rigidBodyModes)
{
  int nsub = decDomain->getNumSub();
  for(int i = 0; i < numRBM(); ++i) {
    execParal(nsub, this, &MultiDomainRbm<Scalar>::getGlobalRBM, i, rigidBodyModes[i]);
  }
}

template<class Scalar>
void
MultiDomainRbm<Scalar>::getRBMs(GenDistrVectorSet<Scalar>& rigidBodyModes, std::set<int> &rbmFilters)
{
  int i = 0, nsub = decDomain->getNumSub();
  for(std::set<int>::iterator it = rbmFilters.begin(); it != rbmFilters.end(); ++it) {
    int iMode = *it;
    if(iMode < numRBM()) {
      execParal(nsub, this, &MultiDomainRbm<Scalar>::getGlobalRBM, iMode, rigidBodyModes[i]);
      i++;
    }
  }
}

template<class Scalar>
void
MultiDomainRbm<Scalar>::getGlobalRBM(int iSub, int &iRBM, GenDistrVector<Scalar> &R)
{
  GenSubDomain<Scalar> *sd = decDomain->getSubDomain(iSub);
  Scalar *localRvec = R.subData(sd->localSubNum());
  sd->getGlobalRBM(iRBM, localRvec);
}




template<class Scalar>
void
MultiDomainRbm<Scalar>::singularValueDecomposition(FullM &A, FullM &U, int ncol, int nrow, int &rank, double tol, FullM *V)
{
  int info = 0;
  int mindim = std::min(nrow,ncol);
  int maxdim = std::max(nrow,ncol);
  double max_value = A.maxAbs();
#ifdef FILERING
  for(int i=0; i<A.numCol()*A.numRow(); i++) //HB
    if(fabs((A.data())[i])<1.E-16*max_values+1.E-20) (A.data())[i] = 0.0;
#endif
  double *w    = (double*) dbg_alloca(sizeof(double)* maxdim);
  double *e    = (double*) dbg_alloca(sizeof(double)* maxdim);
  double *work = (double*) dbg_alloca(sizeof(double)* maxdim);

  for(int i=0; i<maxdim; i++) { w[i] = 0.0; e[i] = 0.0; work[i] = 0.0; }

  if(V == 0)
    Tsvdc(A.data(), ncol, ncol, nrow, w, e, U.data(), ncol, U.data(), ncol, work, 10, info);
  else
    Tsvdc(A.data(), ncol, ncol, nrow, w, e, U.data(), ncol, V->data(), nrow, work, 11, info);

  //  ... CHECK RETURN STATUS OF DSVDC:
  if(info <  0)
    filePrint(stderr," *** ERROR: Illegal value in argument #%d in FetiDPSolver::singularValueDecomposition()\n",info);
  if(info >  0)
    filePrint(stderr," *** WARNING: %d diagonals did not converge in FetiDPSolver::singularValueDecomposition(). Check your result.\n",info);

  //  ... DETERMINE RANK
  rank = 0;
  double tolerance = max_value*tol;
  for(int i=0; i<mindim; ++i) {
    if(fabs(w[i]) > tolerance) rank += 1;
  }
}

template<class Scalar>
DistrInfo &
MultiDomainRbm<Scalar>::solVecInfo()
{
  return decDomain->solVecInfo();
}
