#include <Sfem.d/cijk.h>

extern Sfem *sfem;

template<class Scalar>
SfemBlockMatrix<Scalar>::SfemBlockMatrix(int in_L, int in_n, int in_P, int in_ndim, int output_order)
{
  int i,j,k;
  
  P=in_P;
  n=in_n;
  L=in_L;
  ndim=in_ndim;

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
//  blockToBlock->print(); // YYY DG
//  allK = new GenSparseMatrix<Scalar> *[L];
  allK = new GenSparseMatrix<Scalar> *[ndim+1];
  diags = 0;
  firstdof=0;
  diagblocks=0;
}

template<class Scalar>
void SfemBlockMatrix<Scalar>::setKi(GenSparseMatrix<Scalar>* Ki, int i) 
{
  allK[i] = Ki;
}

template<class Scalar>
void SfemBlockMatrix<Scalar>::matvec_sfem_block(GenVector<Scalar> &u)  
{
  int i,j;
  // For a given set {u_j}_{j=0}^{P-1}, mat-vec K_i*u_j
//  kiuj = new Scalar **[L]; // defined as 3-d array for convenience in accessibility
  kiuj = new Scalar **[ndim+1]; // defined as 3-d array for convenience in accessibility
//  for (i=0; i<L; i++) {
  for (i=0; i<ndim+1; i++) {
    kiuj[i]= new Scalar *[P];
    for (j=0; j<P; j++) {
      kiuj[i][j]= new Scalar[n];
    }
  }

//  for (i=0; i<L; i++) {
  for (i=0; i<ndim+1; i++) {
    for (j=0; j<P; j++) {
//      allK[i]->mult(u.getBlock(j, n), kiuj[i][j]); // mult of dbsparse class
      allK[i]->mult(u.getBlock(j), kiuj[i][j]); // mult of dbsparse class
    }
  }
// Note : kiuj deleted in the mult fn
}

template<class Scalar>
void SfemBlockMatrix<Scalar>::mult(GenVector<Scalar> &u, GenVector<Scalar> &ku) 
{
  // Final big mat-vec   ku = K_big*u  (u is an nP dimensional vector, 
  // given by the iterative linear solver code)
  matvec_sfem_block(u);
  ku.zero();
//  cerr << "within SfemBlockMatrix<Scalar>::mult, n = " << n << endl;
//  ku.setBlockSize(n); 
  int i,j,k,J;
  int* nonzindex = sfem->getnonzindex();
  Scalar *tempvec = new Scalar[n];
  for(k=0; k<blockToBlock->csize(); k++) {
    for (i=0;i<n;i++) {tempvec[i]=0.0;}
    for(J=0; J<blockToBlock->num(k); J++) {
      j = (*blockToBlock)[k][J];
      for (i=0; i<L; i++) {
        if(nonzindex[i]!=-1) {
          for (int ij=0; ij<n; ij++) {
//	      tempvec[ij] += kiuj[i][j][ij]*inpsi[i][j][k];
            tempvec[ij] += kiuj[nonzindex[i]][j][ij]*inpsi[i][j][k];
          }
        }
      }
    }

    for(i=0;i<n;i++) { 
//       ku.setBlockValue(k, i, n, tempvec[i]);
       ku.setBlockValue(k, i, tempvec[i]);
    }
  }

  delete []tempvec;

//  for (i=0; i<L; i++) {
  for (i=0; i<ndim+1; i++) {
    for (j=0; j<P; j++) {
      delete [] kiuj[i][j];
    }
    delete [] kiuj[i];
  }
  delete [] kiuj;
}

template<class Scalar>
Scalar SfemBlockMatrix<Scalar>::diag(int i)
{
  if(!diags) {
    diags = new Scalar[n*P];
    for(int i=0; i<n*P; ++i) diags[i] = 0.0;
    // build the diags vector 
    for(int i=0; i<P; ++i) 
       for(int k=0; k<allK[0]->dim(); ++k) 
         diags[n*i+k] += allK[0]->diag(k)*inpsi[0][i][i]; 
  }
  return diags[i];
}

template<class Scalar>
BlockInfo & SfemBlockMatrix<Scalar>::dim()
{  
  inf.blocklen = n;
  inf.numblocks = P;
  return inf;
}

template<class Scalar>
int* SfemBlockMatrix<Scalar>::getFirstDof()
{
  if(!firstdof) {
    firstdof = new int[P];
    for(int i=0; i<P; i++) firstdof[i]=n*i;
  }
  return firstdof;
}


template<class Scalar>
int SfemBlockMatrix<Scalar>::numNodes() const
{
 return P;
}

template<class Scalar>
GenFullM<Scalar>* SfemBlockMatrix<Scalar>::getDiagMatrix(int i) // Preconditioner with each block-diagonal in GenFullM format 
{
  if(!diagblocks) {
    diagblocks = new GenFullM<Scalar> * [P];

//    cerr << "Block-Jacobi of SFEM called" << endl;

    diagblocks[0] = allK[0]->getFullMatrix(); 
    for(int j=1; j<P; ++j) diagblocks[j] = new GenFullM<Scalar>((*diagblocks[0]),inpsi[0][j][j]);
  }
  diagblocks[i]->symmetrize_from_uptriag();
  return diagblocks[i];

}


template<class Scalar>
Scalar* SfemBlockMatrix<Scalar>::getBlockScalarMultipliers()
{
 scalarfactors=new Scalar[P];
 for(int j=0; j<P; ++j)  scalarfactors[j] = inpsi[0][j][j]; 
 return scalarfactors;
}


template<class Scalar>
void SfemBlockMatrix<Scalar>::setMeanSolver(GenSolver<Scalar> *prc)
{
 meansolver= prc;
}

template<class Scalar>
GenSolver<Scalar>* SfemBlockMatrix<Scalar>::getMeanSolver()
{
 return meansolver;
}




