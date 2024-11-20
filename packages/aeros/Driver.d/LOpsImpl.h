template<class Scalar>
void
GenSubDomain<Scalar>::makeLocalToGlobalDofMap()
{
  l2g = new int[c_dsa->size()];
  int localdofs[DofSet::max_known_dof];
  int globaldofs[DofSet::max_known_dof];
  for(int i = 0; i < c_dsa->numNodes(); ++i) {
    c_dsa->number(i, (*c_dsa)[i], localdofs);
    domain->getCDSA()->number(glNums[i], (*c_dsa)[i], globaldofs);
    for(int j = 0; j < (*c_dsa)[i].count(); ++j) l2g[localdofs[j]] = globaldofs[j];
  }
}

template<class Scalar>
void
GenSubDomain<Scalar>::multAddLT(const Scalar *localvec, Scalar *globalvec)
{
  // globalvec += L^T * localvec
  if(!l2g) makeLocalToGlobalDofMap();
  for(int i = 0; i < c_dsa->size(); ++i)
    globalvec[l2g[i]] += localvec[i];
}

template<class Scalar>
void
GenSubDomain<Scalar>::multAddLinv(const Scalar *localvec, Scalar *globalvec)
{
  // globalvec += L^{-1} * localvec
  if(!l2g) makeLocalToGlobalDofMap();
  for(int i = 0; i < c_dsa->size(); ++i)
    globalvec[l2g[i]] += localvec[i]/double(dofWeight(i));
}

template<class Scalar>
void
GenSubDomain<Scalar>::multLTinv(const Scalar *globalvec, Scalar *localvec)
{
  // localvec = L^{-T} * globalvec
  if(!l2g) makeLocalToGlobalDofMap();
  for(int i = 0; i < c_dsa->size(); ++i)
    localvec[i] = globalvec[l2g[i]]/double(dofWeight(i));
}

template<class Scalar>
void
GenSubDomain<Scalar>::multL(const Scalar *globalvec, Scalar *localvec)
{
  // localvec = L * globalvec
  if(!l2g) makeLocalToGlobalDofMap();
  for(int i = 0; i < c_dsa->size(); ++i)
    localvec[i] = globalvec[l2g[i]];
}
