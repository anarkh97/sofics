#include <Math.d/matrix.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include <Math.d/Vector.h>
#include <Math.d/FullSquareMatrix.h>
#include <Driver.d/Mpc.h>

template<class Scalar>
GenMpcSparse<Scalar>::GenMpcSparse(int numMPC, const std::vector<std::unique_ptr<SubLMPCons<Scalar> > > &_mpc,
                                   const DofSetArray *_dsa) :
	mpc(_mpc),
	NumCol(numMPC),
	NumRow(_dsa->size()),
	dsa(_dsa)
{
}

template<class Scalar>
GenMpcSparse<Scalar>::GenMpcSparse(int numMPC, const std::vector<std::unique_ptr<SubLMPCons<Scalar> > > &_mpc, const DofSetArray *_dsa,
                                   const DofSetArray *_DSA, int *_wetInterfaceMap, int _mpcOffset) :
	mpc(_mpc),
	NumCol(numMPC),
	NumRow(_dsa->size()),
	dsa(_dsa)
{
  // PJSA 10-19-04 for wet interface / mpc interaction
  DSA    = _DSA;
  wetInterfaceMap = _wetInterfaceMap;
  mpcOffset = _mpcOffset;
}

template<class Scalar>
GenMpcSparse<Scalar>::~GenMpcSparse()
{
}

template<class Scalar>
void
GenMpcSparse<Scalar>::zeroAll()
{
}

template<class Scalar>
void
GenMpcSparse<Scalar>::negate()
{
}

template<class Scalar>
void
GenMpcSparse<Scalar>::print(FILE *file,const char*)
{
  int i,iMPC;
  for(iMPC=0; iMPC<NumCol; ++iMPC) {
    fprintf(file," MPC %d\n",iMPC+1);
    for(i=0; i<mpc[iMPC]->nterms; ++i)
      fprintf(file," %d  node %d  dof %d coef %g    rhs %g \n",
                   i+1,mpc[iMPC]->terms[i].nnum+1,
                   mpc[iMPC]->terms[i].dofnum,
                   mpc[iMPC]->terms[i].coef,
                   mpc[iMPC]->rhs);
  }
}

template<class Scalar>
double
GenMpcSparse<Scalar>::getMemoryUsed() const
{
  return 0.0;
}

template<class Scalar>
void
GenMpcSparse<Scalar>::multSubtract(const GenVector<Scalar> &, GenVector<Scalar> &) const
{
  fprintf(stderr,"GenMpcSparse<Scalar> does not support multSubtract(const GenVector<Scalar> &, GenVector<Scalar> &) \n");
}

template<class Scalar>
void
GenMpcSparse<Scalar>::multSubtract(const Scalar *, Scalar *) const
{
  fprintf(stderr,"GenMpcSparse<Scalar> does not support multSubtract(const Scalar *, Scalar *) \n");
}

template<class Scalar>
void
GenMpcSparse<Scalar>::mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const
{
  mult(rhs.data(),result.data() );
}

template<class Scalar>
void
GenMpcSparse<Scalar>::mult(const Scalar *rhs, Scalar *result) const
{
  int i,iMPC;
  for(iMPC=0; iMPC<NumCol; ++iMPC) {
    result[iMPC]=0.0;
    for(i=0; i<mpc[iMPC]->nterms; ++i) {
      int dof = dsa->locate(mpc[iMPC]->terms[i].nnum, (1 << mpc[iMPC]->terms[i].dofnum));
      if(dof >= 0) 
        result[iMPC] += mpc[iMPC]->terms[i].coef * rhs[dof];
    }
  }
}

template<class Scalar>
void
GenMpcSparse<Scalar>::multIdentity(int iMPC, Scalar *result) const
{
  int i;
  for(i=0; i<mpc[iMPC]->nterms; ++i) {
    int dof = dsa->locate(mpc[iMPC]->terms[i].nnum, (1 << mpc[iMPC]->terms[i].dofnum));
    if(dof >= 0)
      result[dof] += mpc[iMPC]->terms[i].coef;
  }
}

template<class Scalar>
void
GenMpcSparse<Scalar>::multIdentity(Scalar **result) const
{
  int i,iMPC;
  for(iMPC=0; iMPC<NumCol; ++iMPC) {
    for(i=0; i<mpc[iMPC]->nterms; ++i) {
      int dof = dsa->locate(mpc[iMPC]->terms[i].nnum, (1 << mpc[iMPC]->terms[i].dofnum));
      if(dof >= 0)
        result[iMPC][dof] += mpc[iMPC]->terms[i].coef;
    }
  }
}

template<class Scalar>
void
GenMpcSparse<Scalar>::multIdentity(Scalar *) const
{
  fprintf(stderr,"GenMpcSparse<Scalar> does not support multIdentity(double*)\n");
}

// vc = Kcr Krr^-1 Krc
//
// vr = Krr^-1 Krc (one column of Krc at a time)
// vc = Kcr vr
//
template<class Scalar>
void
GenMpcSparse<Scalar>::multSub(const Scalar *rhs, Scalar *result) const
{
  int i,iMPC;
  for(iMPC=0; iMPC<NumCol; ++iMPC) {
    for(i=0; i<mpc[iMPC]->nterms; ++i) {
      int dof = dsa->locate(mpc[iMPC]->terms[i].nnum, (1 << mpc[iMPC]->terms[i].dofnum));
      if(dof >= 0)
        result[iMPC] -= mpc[iMPC]->terms[i].coef * rhs[dof];
    }
  }
}

template<class Scalar>
void
GenMpcSparse<Scalar>::multSubWI(const Scalar *wi_rhs, Scalar *result) const
{
  // PJSA: for wet interface / mpc interaction
  int i,iMPC;
  for(iMPC=0; iMPC<NumCol; ++iMPC) {
    for(i=0; i<mpc[iMPC]->nterms; ++i) {
      int dof = DSA->locate(mpc[iMPC]->terms[i].nnum, (1 << mpc[iMPC]->terms[i].dofnum));
      if((dof >= 0) && (wetInterfaceMap[dof] >= 0)) 
        result[iMPC+mpcOffset] -= mpc[iMPC]->terms[i].coef * wi_rhs[wetInterfaceMap[dof]];
    }
  }
}

template<class Scalar>
void
GenMpcSparse<Scalar>::transposeMultSubtract(const Scalar *rhs, Scalar *result) const
{
  int i,iMPC;
  for(iMPC=0; iMPC<NumCol; ++iMPC) {
    for(i=0; i<mpc[iMPC]->nterms; ++i) {
      int dof = dsa->locate(mpc[iMPC]->terms[i].nnum, (1 << mpc[iMPC]->terms[i].dofnum));
      if(dof >= 0)
        result[dof] -= mpc[iMPC]->terms[i].coef * rhs[iMPC];
    }
  }
}

template<class Scalar>
void
GenMpcSparse<Scalar>::transposeMultSubtractWI(const Scalar *rhs, Scalar *wi_result) const
{
  // PJSA: for wet interface / mpc interaction
  int i,iMPC;
  for(iMPC=0; iMPC<NumCol; ++iMPC) {
    for(i=0; i<mpc[iMPC]->nterms; ++i) {
      int dof = DSA->locate(mpc[iMPC]->terms[i].nnum, (1 << mpc[iMPC]->terms[i].dofnum));
      if((dof >= 0) && (wetInterfaceMap[dof] >= 0)) 
        wi_result[wetInterfaceMap[dof]] -= mpc[iMPC]->terms[i].coef * rhs[iMPC+mpcOffset];
    }
  }
}

template<class Scalar>
void
GenMpcSparse<Scalar>::transposeMult(const Scalar *, Scalar *) const
{
  fprintf(stderr,"GenMpcSparse<Scalar> does not support transposeMult(const Scalar *, Scalar *)\n");
}

template<class Scalar>
void
GenMpcSparse<Scalar>::multSub(int, Scalar **, Scalar **) const
{
  fprintf(stderr,"GenMpcSparse<Scalar> does not support multSub(int nRHS, Scalar **, Scalar **)\n");
}

template<class Scalar>
void
GenMpcSparse<Scalar>::multAdd(const Scalar *rhs, Scalar *result) const
{
  int i,iMPC;
  for(iMPC=0; iMPC<NumCol; ++iMPC) {
    for(i=0; i<mpc[iMPC]->nterms; ++i) {
      int dof = dsa->locate(mpc[iMPC]->terms[i].nnum, 
                            (1 << mpc[iMPC]->terms[i].dofnum));
      if(dof >= 0)
        result[iMPC] += mpc[iMPC]->terms[i].coef * rhs[dof];
    }
  }
}

template<class Scalar>
void
GenMpcSparse<Scalar>::transposeMultAdd(const Scalar *rhs, Scalar *result) const
{
   int i,iMPC;
   for(iMPC=0; iMPC<NumCol; ++iMPC) {
     for(i=0; i<mpc[iMPC]->nterms; ++i) {
       int dof = dsa->locate(mpc[iMPC]->terms[i].nnum, (1 << mpc[iMPC]->terms[i].dofnum));
       if(dof >= 0)
         result[dof] += mpc[iMPC]->terms[i].coef * rhs[iMPC];
     }
  }
}

template<class Scalar>
void
GenMpcSparse<Scalar>::addBoeing(int, const int *, const int *, const double *, const int *, Scalar)
{
  fprintf(stderr,"GenMpcSparse<Scalar> does not support addBoeing(...) \n");
}

template<class Scalar>
void
GenMpcSparse<Scalar>::add(const FullSquareMatrix&, const int*)
{
  fprintf(stderr,"GenMpcSparse<Scalar> does not support add(FullSquareMatrix&, int *)\n");
}

template<class Scalar>
void
GenMpcSparse<Scalar>::add(const double *const *, const int *, int)
{
  fprintf(stderr,"GenMpcSparse<Scalar> does not support add(const double *const *, int *, int)\n");
}

template<class Scalar>
void
GenMpcSparse<Scalar>::add(const FullM&, int, int)
{
  fprintf(stderr,"GenMpcSparse<Scalar> does not support add(FullM&, int, int)\n");
}

template<class Scalar>
void
GenMpcSparse<Scalar>::add(const GenAssembledFullM<Scalar>&, const int*)
{
  fprintf(stderr,"GenMpcSparse<Scalar> does not support add(GenAssembledFullM<Scalar>&, int*)\n");
}
