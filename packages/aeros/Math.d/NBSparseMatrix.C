#include <Utils.d/dbg_alloca.h>
#include <Utils.d/Connectivity.h>
#include <Math.d/Vector.h>
#include <Math.d/NBSparseMatrix.h>

template<class Scalar>
GenFullM<Scalar> *
GenNBSparseMatrix<Scalar>::getDiagMatrix(int i)
{
 int offset = con->offset(i,i);
 return (offset < 0) ? 0 : allM+offset;
}

template<class Scalar>
void
GenNBSparseMatrix<Scalar>::zeroAll()
{
 int inode,j,k,l;
 int im = 0;
 for(inode = 0; inode < numnodes; ++inode) {
   int numRows = dsa->weight(inode);
   for(j = 0; j < con->num(inode); ++j) {
       int jnode  = (*con)[inode][j];
       int numCol = dsa->weight(jnode);
       for(k=0; k < numRows; ++k)
         for(l=0; l< numCol; ++l)
           allM[im][k][l] = 0.0;
       im++;
     }
  }
}

template<class Scalar>
void
GenNBSparseMatrix<Scalar>::add(const FullSquareMatrix &kel, const int *dofs)
{
 int i, j, dof;

 // For each dof figure out to which node it belongs and
 // its unconstrained number

 int kndof = kel.dim();

 // memory for dof to node tables
 int *dton = (int *) dbg_alloca(sizeof(int)*kndof);
 int *dn   = (int *) dbg_alloca(sizeof(int)*kndof);

 for(i = 0; i < kndof; ++i) {
   dof = dsa->getRCN(dofs[i]);
   if(dof == -1)
     dn[i] = -1;
   else {
     dton[i] = dofToNode[ dof ];
     dn[i]   = dof - dsa->firstdof( dton[i] );
   }
 }

// Now add the contribution to the matrix
  for(i = 0; i < kndof; ++i) {
    if(dn[i] == -1) continue;
    for(j=0; j < kndof; ++j) {
      if(dn[j] == -1) continue;
      int matid = con->offset(dton[i],dton[j]);
      allM[matid][dn[i]][dn[j]] += kel[i][j];
    }
  }
}

template<class Scalar>
GenNBSparseMatrix<Scalar>::GenNBSparseMatrix(const Connectivity *_cn, const ConstrainedDSA *c_dsa)
{
  con = _cn;
  dsa = c_dsa;

  int inode, j;

  numnodes = dsa->numNodes();
  if(con->csize() < numnodes) numnodes = con->csize();

  dofToNode = new int[dsa->size()];

  for(inode=0; inode < numnodes; ++inode) {
   int fDof = dsa->firstdof(inode);
   int lastDof = fDof + dsa->weight(inode);
   for(j = fDof; j < lastDof; ++j)
     dofToNode[j] = inode;
  }

  allM     = new GenFullM<Scalar>[con->numConnect()];
  firstDof = new int[numnodes];

  int im=0;
  for(inode = 0; inode < numnodes; ++inode) {
    firstDof[inode] = dsa->firstdof(inode);
    int numRows = dsa->weight(inode);
    for(j = 0; j < con->num(inode); ++j) {
      int jnode  = (*con)[inode][j];
      int numCol = dsa->weight(jnode);
      allM[im].setNewSize(numRows,numCol);
      im++;
    }
  }
}

template<class Scalar>
GenNBSparseMatrix<Scalar>::~GenNBSparseMatrix()
{
  if(dofToNode) { delete [] dofToNode; dofToNode = 0; }
  if(firstDof)  { delete [] firstDof;  firstDof  = 0; }
  if(allM) { delete [] allM; allM = 0; }
}

template<class Scalar>
void
GenNBSparseMatrix<Scalar>::mult(const Scalar *rhs, Scalar *result) const
{
   int im = 0;
   int i,j,k,l;

   // zero the results vector
   for(i=0; i<dsa->size(); ++i)
     result[i] = 0.0;

   for(i = 0; i < numnodes; ++i) {
     int fi = dsa->firstdof(i);
     for(j = 0; j < con->num(i); ++j) {
       int fj = dsa->firstdof( (*con)[i][j] );
       int nk = allM[im].numRow();
       int nl = allM[im].numCol();
       for(k = 0; k < nk; ++k)
         for(l = 0; l < nl; ++l)
           result[fi+k] += rhs[fj+l]*allM[im][k][l];
       ++im;
     }
   }
}

template<class Scalar>
void
GenNBSparseMatrix<Scalar>::multAdd(const Scalar *rhs, Scalar *result) const
{
   int im = 0;
   int i,j,k,l;

   for(i = 0; i < numnodes; ++i) {
     int fi = dsa->firstdof(i);
     for(j = 0; j < con->num(i); ++j) {
       int fj = dsa->firstdof( (*con)[i][j] );
       int nk = allM[im].numRow();
       int nl = allM[im].numCol();
       for(k = 0; k < nk; ++k)
         for(l = 0; l < nl; ++l)
           result[fi+k] += rhs[fj+l]*allM[im][k][l];
       ++im;
     }
   }
}

template<class Scalar>
void
GenNBSparseMatrix<Scalar>::mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const
{
   int im = 0;
   int i,j,k,l;

   // zero the results vector
   result.zero();

   for(i = 0; i < numnodes; ++i) {
     int fi = dsa->firstdof(i);
     for(j = 0; j < con->num(i); ++j) {
       int fj = dsa->firstdof( (*con)[i][j] );
       int nk = allM[im].numRow();
       int nl = allM[im].numCol();
       for(k = 0; k < nk; ++k)
         for(l = 0; l < nl; ++l)
           result[fi+k] += rhs[fj+l]*allM[im][k][l];
       ++im;
     }
   }
}

template<class Scalar>
Scalar
GenNBSparseMatrix<Scalar>::diag(int dof) const
{
// For this dof figure out to which node it belongs and
// its number as unconstrained

  int dton = dofToNode[dof];
  int dn = dof - dsa->firstdof(dton);
  int matid = con->offset(dton, dton);

// Check to see if the diagonal is 0 and if the bc is fixed
// Then return 1 or the diagonal entry of the matrix

// if( allM[matid][dn][dn] == 0.0 && bc[dof] == BCFIXED )

  if(allM[matid][dn][dn] == 0.0)
    return (1);
  else
    return allM[matid][dn][dn];
}

template<class Scalar>
Scalar &
GenNBSparseMatrix<Scalar>::diag(int dof)
{
// For this dof figure out to which node it belongs and
// its number as unconstrained

  int dton = dofToNode[dof];
  int dn = dof - dsa->firstdof(dton);
  int matid = con->offset(dton, dton);

    return allM[matid][dn][dn];
}

