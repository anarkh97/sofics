#include <Sfem.d/cijk.h>

template<class Scalar>
DistrBlockVector<Scalar>::~DistrBlockVector()
{
 for(int i=0; i<inf.nblocks; ++i) delete v[i];
 delete [] v;
}


template<class Scalar>
int DistrBlockVector<Scalar>::size() const
{
 int curlen=0;
/* DistrInfo dinfo;                  
 for(int i=0; i<inf.nblocks; ++i) {
   dinfo = v[i]->info();
   curlen = curlen + dinfo.totLen();
 }*/
 for(int i=0; i<inf.nblocks; ++i) { curlen = curlen + v[i]->size();}  
 return curlen;
}


template<class Scalar>
void DistrBlockVector<Scalar>::zero()
{
 for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) v[i]->zero();
}


template<class Scalar>
double DistrBlockVector<Scalar>::norm()
{
  double tempnorm = 0;
  for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) tempnorm = tempnorm + pow(v[i]->norm(),2);
  tempnorm = sqrt(tempnorm);
  return tempnorm; 
}


template<class Scalar>
double DistrBlockVector<Scalar>::sqNorm()
{
  double tempnorm = 0;
  for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) tempnorm = tempnorm + pow(v[i]->norm(),2);
  return tempnorm;
}


template<class Scalar>
Scalar DistrBlockVector<Scalar>::operator * (DistrBlockVector<Scalar> &x)
{
 if(this->size() != x.size()) std::cerr << "Dimensions don't match in * operation " << std::endl;
 Scalar temp=0;
 for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) temp=temp+(*(v[i]))*(x.getBlock(i));
 return temp;
}


template<class Scalar>
DistrBlockVector<Scalar>&  DistrBlockVector<Scalar>::operator=(const DistrBlockVector<Scalar> &x)
{
 if(this->size() != x.size()) std::cerr << "Dimensions don't match in = operation " << std::endl;
 for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) *(v[i])= *(x.getv()[i]);
 return *this;
}


template<class Scalar>
DistrBlockVector<Scalar>&  DistrBlockVector<Scalar>::operator=(Scalar c)
{
 for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) *(v[i]) = c;
 return *this;
}


template<class Scalar>
DistrBlockVector<Scalar>&  DistrBlockVector<Scalar>::operator*=(Scalar c)
{
 for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) *(v[i]) *= c;
 return *this;
}


template<class Scalar>
DistrBlockVector<Scalar>&  DistrBlockVector<Scalar>::operator+=(DistrBlockVector<Scalar> &x)
{
 if(size() != x.size()) std::cerr << "Dimensions don't match in +=" << std::endl;
 for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) *(v[i]) += x.getBlock(i);
 return *this;
}


template<class Scalar>
DistrBlockVector<Scalar>&  DistrBlockVector<Scalar>::operator-=(DistrBlockVector<Scalar> &x)
{
 if(size() != x.size()) std::cerr << "Dimensions don't match in -=" << std::endl;
 for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) *(v[i]) -= x.getBlock(i);
 return *this;
}


template<class Scalar>
DistrBlockVector<Scalar>&  DistrBlockVector<Scalar>::linAdd(DistrBlockVector<Scalar> &x)
{
 if(size() != x.size()) std::cerr << "Dimensions don't match in linAdd" << std::endl;
 for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) v[i]->linAdd(x.getBlock(i));
 return *this;
}


template<class Scalar>
DistrBlockVector<Scalar>&  DistrBlockVector<Scalar>::linAdd(Scalar c, DistrBlockVector<Scalar> &x)
{
 if(size() != x.size()) std::cerr << "Dimensions don't match in linAdd" << std::endl;
 for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) v[i]->linAdd(c,x.getBlock(i));
 return *this;
}


template<class Scalar>
DistrBlockVector<Scalar>&  DistrBlockVector<Scalar>::linAdd(Scalar c1, DistrBlockVector<Scalar> &x, Scalar c2, DistrBlockVector<Scalar> &y)
{
 if(x.size() != y.size()) std::cerr << "Dimensions don't match in linAdd" << std::endl;
 for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) v[i]->linAdd(c1,x.getBlock(i),c2,y.getBlock(i)); 
 return *this;
}


template<class Scalar>
DistrBlockVector<Scalar>&  DistrBlockVector<Scalar>::linC(const DistrBlockVector<Scalar> &x, Scalar c)
{
 if(size() != x.size()) std::cerr << "Dimensions don't match in linC" << std::endl;
 for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) v[i]->linC(c,x.getBlock(i));
 return *this;
}


template<class Scalar>
DistrBlockVector<Scalar>&  DistrBlockVector<Scalar>::linC(const DistrBlockVector<Scalar> &x, Scalar c, const DistrBlockVector<Scalar> &y)
{
 if(x.size() != y.size()) std::cerr << "Dimensions don't match in linC" << std::endl;
 for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) v[i]->linC(x.getBlock(i),c,y.getBlock(i)); 
 return *this;
}


template<class Scalar>
DistrBlockVector<Scalar>&  DistrBlockVector<Scalar>::linC(Scalar c1, const DistrBlockVector<Scalar> &x, Scalar c2, const DistrBlockVector<Scalar> &y)
{
 if(x.size() != y.size()) std::cerr << "Dimensions don't match in linC" << std::endl;
 for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) v[i]->linC(c1,x.getBlock(i),c2,y.getBlock(i)); 
 return *this;
}


template<class Scalar>
DistrBlockVector<Scalar>&  DistrBlockVector<Scalar>::swap(DistrBlockVector<Scalar> &x)
{
 if(size() != x.size()) std::cerr << "Dimensions don't match in swap" << std::endl;
 for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) v[i]->swap(x.getBlock(i));
 return *this;
}


template<class Scalar>
Scalar DistrBlockVector<Scalar>::sum()
{
 Scalar temp=0;
 for(int i=0; i<inf.nblocks; ++i) if(inf.nnzblkindex[i]==1) temp=temp+v[i]->sum();
 return temp;
}

template<class Scalar>
void DistrBlockVector<Scalar>::updateBlock(int ii, Scalar c, DistrBlockVector<Scalar> &x)
{
 v[ii]->linAdd(c,x.getBlock(0));
}


template<class Scalar>
void DistrBlockVector<Scalar>::copyBlock(DistrBlockVector<Scalar> &x, int ii)
{
 v[0]->zero();
 v[0]->linAdd(1,x.getBlock(ii));
}


template<class Scalar>
void DistrBlockVector<Scalar>::addBlockSqr(int ii, Scalar c, DistrBlockVector<Scalar> &x)
{ 
 v[0]->addBlockSqr(ii, c, x.getBlock(ii)); // v = v + pow(x[ii],2)*c
}


template<class Scalar>
void DistrBlockVector<Scalar>::computeSqrt()
{
  v[0]->computeSqrt();
}


template<class Scalar>
void DistrBlockVector<Scalar>::computeRealz(int ii, Scalar c, DistrBlockVector<Scalar> &x)
{
 v[0]->linAdd(c,x.getBlock(ii));
}



template<class Scalar>
void DistrBlockVector<Scalar>::print()
{
  for(int i=0; i<inf.nblocks; ++i) { std::cerr << "Block number = " << i << std::endl; v[i]->print();}
}


template<class Scalar>
void DistrBlockVector<Scalar>::printNonZeroTerms()
{
  std::cerr << "DistrBlockVector<Scalar>::printNonZeroTerms() not implemented" << std::endl;
}


template<class Scalar>
void DistrBlockVector<Scalar>::printAll()
{
  std::cerr << "DistrBlockVector<Scalar>::printAll() not implemented" << std::endl;
}


template<class Scalar>
void DistrBlockVector<Scalar>::computeBlockNorms()
{
 for(int i=0; i<inf.nblocks; ++i) blocknorms[i]= pow(v[i]->norm(),2); // YYY DG actually the square of the norm
}


template<class Scalar>
void DistrBlockVector<Scalar>::printBlockNorms()
{
 std::cerr << "DistrBlockVector::printBlockNorms() called" << std::endl;
 std::cerr << "blocknorms are = ";
 for(int i=0; i<inf.nblocks; ++i) std::cerr << "  " << blocknorms[i];
 std::cerr << std::endl;
}

template<class Scalar>
void DistrBlockVector<Scalar>::printBlockDetails()   
{
 std::cerr << "inf.nblocks = "; 
 for (int i=0; i<inf.nblocks; ++i) std::cerr << inf.nnzblkindex[i] << " "; 
 std::cerr << std::endl;
} 


// ########

template<class Scalar>
DistrSfemBlockMatrix<Scalar>::DistrSfemBlockMatrix(int __L, int __P, int _ndim, int output_order, DistrBlockInfo & info) : inf(info)
{
  int i,j,k;
  
  P=__P;
  n=inf.blockinfo[0].totLen();
  L=__L;
  ndim=_ndim;  

  inpsi= new double**[P];
  for (i=0; i<P; i++) {
    inpsi[i]=new double*[P];
    for (j=0; j<P; j++) {
      inpsi[i][j]= new double[P];
    }
  }
  Cijk cijks(ndim, output_order, P);
  for (i=0; i<P; i++) {
    for (j=0; j<P; j++) {
      for (k=0; k<P; k++) {
        inpsi[i][j][k] = cijks.expectation(i,j,k);
      }
    }
  }

  int *ptr = new int[P+1]; ptr[0] = 0;
  int *tgt = new int[P*P];
  int count = 0;
  for (k=0; k<P; k++) { // loop over rows
    for (j=0; j<P; j++) { // loop over columns
      for (i=0; i<L; i++) { // look for non-zero term in expansion
        if(inpsi[i][j][k] != 0.0) {
          tgt[count++] = j;
          break;
        }
      }
    }
    ptr[k+1] = count;
  } 


  blockToBlock = new Connectivity(P, ptr, tgt);
//  blockToBlock->print(); // YYY
//  allK = new SubDOp *[L];
  allK = new SubDOp *[ndim+1];
  diags = 0;
  firstdof=0;
  //diagblocks=0;
}


template<class Scalar>
void DistrSfemBlockMatrix<Scalar>::setKi(SubDOp* Ki, int i) 
{
  allK[i] = Ki;
}

template<class Scalar>
void DistrSfemBlockMatrix<Scalar>::matvec_sfem_block(DistrBlockVector<Scalar> &u)  
{
  int i,j;
  // For a given set {u_j}_{j=0}^{P-1}, mat-vec K_i*u_j
//  kiuj = new GenDistrVector<Scalar> **[L]; // defined as 3-d array for convenience in accessibility
  kiuj = new GenDistrVector<Scalar> **[ndim+1]; // defined as 3-d array for convenience in accessibility
//  for (i=0; i<L; i++) {
  for (i=0; i<ndim+1; i++) {
    kiuj[i]= new GenDistrVector<Scalar> *[P];
    for (j=0; j<P; j++) {
     if(inf.nnzblkindex[j]==1) kiuj[i][j]= new GenDistrVector<Scalar>(inf.blockinfo[j]);
    }
  }

//  for (i=0; i<L; i++) {
  for (i=0; i<ndim+1; i++) {
    for (j=0; j<P; j++) {
      if(inf.nnzblkindex[j]==1) allK[i]->multNoAssemble(u.getBlock(j), *(kiuj[i][j])); // mult of SubDOp class 
      // Use the "mult" function defined in every storage/solve class
    }
  }
// Note : kiuj deleted in the mult fn
}

template<class Scalar>
void DistrSfemBlockMatrix<Scalar>::mult(DistrBlockVector<Scalar> &u, DistrBlockVector<Scalar> &ku) 
{
  // Final big mat-vec   ku = K_big*u  (u is an nP dimensional vector, 
  // given by the iterative linear solver code)
  matvec_sfem_block(u);

  ku.zero();
  int i,j,k,J;
  int* nonzindex = sfem->getnonzindex();
  GenDistrVector<Scalar> tempvec(inf.blockinfo[0]); 
  for(k=0; k<blockToBlock->csize(); k++) {
    tempvec.zero();
    for(J=0; J<blockToBlock->num(k); J++) {
      j = (*blockToBlock)[k][J];
	for (i=0; i<L; i++) {
          if(inf.nnzblkindex[j]==1) if(nonzindex[i]!=-1) tempvec.linAdd(inpsi[i][j][k],*(kiuj[nonzindex[i]][j])); 
	}
    }
    allK[0]->assemble(tempvec); 
    ku.setBlock(k, tempvec);
  }


  //delete []tempvec;

//  for (i=0; i<L; i++) {
  for (i=0; i<ndim+1; i++) {
    for (j=0; j<P; j++) {
      if(inf.nnzblkindex[j]==1) delete kiuj[i][j];
    }
    delete [] kiuj[i];
  }
  delete [] kiuj;
}

template<class Scalar>
DistrBlockInfo & DistrSfemBlockMatrix<Scalar>::dim()
{  
  return inf;
}

template<class Scalar>
int* DistrSfemBlockMatrix<Scalar>::getFirstDof()
{
  if(!firstdof) {
    firstdof = new int[P];
    for(int i=0; i<P; i++) firstdof[i]=n*i;
  }
  return firstdof;
}


template<class Scalar>
int DistrSfemBlockMatrix<Scalar>::numNodes() const
{
 return P;
}

template<class Scalar>
GenFullM<Scalar>* DistrSfemBlockMatrix<Scalar>::getDiagMatrix(int i)
{
 std::cerr << " DistrSfemBlockMatrix<Scalar>::getDiagMatrix(int i) is not implemented \n";
}

template<class Scalar>
Scalar* DistrSfemBlockMatrix<Scalar>::getBlockScalarMultipliers()
{
 scalarfactors=new Scalar[P];
 for(int j=0; j<P; ++j)  scalarfactors[j] = inpsi[0][j][j]; 
 return scalarfactors;
}


template<class Scalar>
void DistrSfemBlockMatrix<Scalar>::setMeanSolver(GenFetiSolver<Scalar> *prc)
{
 meansolver= prc;
}

template<class Scalar>
GenFetiSolver<Scalar>* DistrSfemBlockMatrix<Scalar>::getMeanSolver()
{
 return meansolver;
}




