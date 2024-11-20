template<class Scalar>
SubBlockCCtSolver<Scalar>::SubBlockCCtSolver(const Connectivity *mpcToMpc, const Connectivity *_mpcToSub,
                                             int _numSubsWithMpcs, std::vector<FetiSub<Scalar> *> subsWithMpcs,
                                             FSCommunicator *_fetiCom, const Connectivity *_cpuToSub) :
		CCtSolver<Scalar>(std::move(subsWithMpcs))
{
  mpcToSub = _mpcToSub;
  this->numSubsWithMpcs = _numSubsWithMpcs;
  this->fetiCom = _fetiCom;
  cpuToSub = _cpuToSub;
  myCPU = _fetiCom->cpuNum();
  
  paralApply(this->subsWithMpcs, &FetiSub<Scalar>::makeLocalMpcToGlobalMpc, mpcToMpc);
  mpcvPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, myCPU, FSCommPattern<Scalar>::CopyOnSend);
  for(int i=0; i<this->numSubsWithMpcs; ++i) this->subsWithMpcs[i]->setMpcCommSize(mpcvPat);  // this is used for combineMpcInterfaceVec()
  mpcvPat->finalize();
  paralApplyToAll(this->numSubsWithMpcs, this->subsWithMpcs, &FetiSub<Scalar>::constructLocalCCtsolver);
  cctPat = 0;
}

template<class Scalar>
SubBlockCCtSolver<Scalar>::~SubBlockCCtSolver()
{
  paralApplyToAll(this->subsWithMpcs, &FetiSub<Scalar>::deleteLocalCCtsolver);
  delete mpcvPat;
  delete cctPat;
}

template<class Scalar>
void
SubBlockCCtSolver<Scalar>::assemble()
{
  paralApplyToAll(this->subsWithMpcs, &FetiSub<Scalar>::assembleLocalCCtsolver);
  if(!cctPat) {
    cctPat = new FSCommPattern<Scalar>(this->fetiCom, this->cpuToSub, this->myCPU, FSCommPattern<Scalar>::CopyOnSend,
                                       FSCommPattern<Scalar>::NonSym);
    for(int i=0; i<this->numSubsWithMpcs; ++i) this->subsWithMpcs[i]->setCCtCommSize(cctPat);
    cctPat->finalize();
  }
  paralApply(this->subsWithMpcs, &FetiSub<Scalar>::sendNeighbCCtsolver, cctPat, mpcToSub);
  cctPat->exchange();
  paralApply(this->subsWithMpcs, &FetiSub<Scalar>::recNeighbCCtsolver, cctPat, mpcToSub);
}

template<class Scalar>
void
SubBlockCCtSolver<Scalar>::factor()
{
  paralApplyToAll(this->numSubsWithMpcs, this->subsWithMpcs, &FetiSub<Scalar>::factorLocalCCtsolver);
}

template<class Scalar>
void
SubBlockCCtSolver<Scalar>::zeroAll()
{
  paralApplyToAll(this->numSubsWithMpcs, this->subsWithMpcs, &FetiSub<Scalar>::zeroLocalCCtsolver);
}

template<class Scalar>
void
SubBlockCCtSolver<Scalar>::reSolve(GenDistrVector<Scalar> &v)
{
  execParal(this->numSubsWithMpcs, this, &SubBlockCCtSolver<Scalar>::solveLocalCCt, v);
  execParal(this->numSubsWithMpcs, this, &SubBlockCCtSolver<Scalar>::sendMpcInterfaceVec, v);
  mpcvPat->exchange();
  execParal(this->numSubsWithMpcs, this, &SubBlockCCtSolver<Scalar>::combineMpcInterfaceVec, v);
}

template<class Scalar>
void
SubBlockCCtSolver<Scalar>::solveLocalCCt(int iSub, GenDistrVector<Scalar> &v)
{
  Scalar *subv = v.subData(this->subsWithMpcs[iSub]->localSubNum());
  this->subsWithMpcs[iSub]->solveLocalCCt(subv);
}

template<class Scalar>
void
SubBlockCCtSolver<Scalar>::sendMpcInterfaceVec(int iSub, GenDistrVector<Scalar> &v)
{
  Scalar *subv = v.subData(this->subsWithMpcs[iSub]->localSubNum());
  this->subsWithMpcs[iSub]->sendMpcInterfaceVec(mpcvPat, subv);
}
template<class Scalar>
void
SubBlockCCtSolver<Scalar>::combineMpcInterfaceVec(int iSub, GenDistrVector<Scalar> &v)
{
  Scalar *subv = v.subData(this->subsWithMpcs[iSub]->localSubNum());
  this->subsWithMpcs[iSub]->combineMpcInterfaceVec(mpcvPat, subv);
}
