#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <Math.d/BLKSparseMatrix.h>
#include <Utils.d/Memory.h>
#include <Utils.d/linkfc.h>
#include <Element.d/Element.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/matrix.h>
#include <Solvers.d/Rbm.h>
#include <Timers.d/GetTime.h>
#include <Solvers.d/SolverCntl.h>
#include <Driver.d/Communicator.h>

#define MIN_MEMORY

#ifdef USE_METIS
//extern "C" void METIS_NodeND(int *n, int * xadj, int *adj, int *numflag, int *options, int *perm, int* iperm);
#include <metis.h>

int callMetis(idx_t *nvtxs, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
              idx_t *options, idx_t *perm, idx_t *iperm)
{
    return METIS_NodeND(nvtxs, xadj, adjncy, vwgt, options, perm, iperm);
}

template <typename T>
int callMetis(T *n, T *xadj, T *adj, T *weights, idx_t *options, T *perm,
    T *iperm)
{
    idx_t m = *n;
    std::vector<idx_t> mxadj(xadj, xadj + *n + 1);
    std::vector<idx_t> madj (adj, adj+xadj[*n + 1]);
    std::vector<idx_t> mweights = weights != nullptr
                                  ? std::vector<idx_t>(weights, weights + *n)
                                  : std::vector<idx_t>() ;
    std::vector<idx_t> mperm(*n);
    std::vector<idx_t> miperm(*n);
    auto status = METIS_NodeND(&m, mxadj.data(), madj.data(), weights ? mweights.data() : nullptr,
        options, mperm.data(), miperm.data());
    for(T i = 0; i < *n; ++i) {
        perm[i] = mperm[i];
        iperm[i] = miperm[i];
    }
    return status;
}

template <typename T>
int callMetis(T *n, T *xadj, T *adj, std::nullptr_t weights, idx_t *options, T *perm,
              T *iperm)
{
    idx_t m = *n;
    std::vector<idx_t> mxadj(xadj, xadj + *n +1);
    std::vector<idx_t> madj (adj, adj+xadj[*n]);
    std::vector<idx_t> mperm(*n);
    std::vector<idx_t> miperm(*n);
    auto status = METIS_NodeND(&m, mxadj.data(), madj.data(), nullptr,
                               options, mperm.data(), miperm.data());
    for(T i = 0; i < *n; ++i) {
        perm[i] = mperm[i];
        iperm[i] = miperm[i];
    }
    return status;
}
#endif

extern "C" {

void _FORTRAN(ordmmd2)(int& n, int* xadj, int* adj, int* invp, int* perm,
                      int& iwsize, int* iwork, int& nofsub, int& iflag);

void _FORTRAN(sfinit)(int& n,   int& nnza,  int* xadj,  int* adj,   int* perm,
                      int* invp,int& maxsup,int& defblk,int* colcnt,int& nnzl,
                      int& nsub,int& nsuper,int* xsuper,int* snode, int& iwsize,
                      int* iwork, int& iflag);

void _FORTRAN(symfct)(int& n, int& nnza, int* xadj, int* adj, int* perm,
                      int* invp, int* colcnt, int& nsuper, int* xsuper,
                      int* snode, int& nofsub, int* xlindx, int* lindx,
                      int* xlnz, int& iwsize, int* iwork, int& iflag);

void _FORTRAN(bfinit)(int& nsuper, int* xsuper, int* snode, int* xlindx,
                      int* lindx, int& tmpsiz, int& rwsize);

void _FORTRAN(blkldl)(int& nsuper, int* xsuper, int *pnode, int* xlindx,
                      int* lindx,  int* xlnz, double *lnz, int& defblk,
                      int &asdef,  int& numZEM, int& lbdef, int *def,
                      double& tol,
                      int *iprow, int* ipcol, int& tmpsiz, double *temp,
                      int& iwsize, int* iwork, int& rwsize, double *rwork,
                      int& iflag);

void _FORTRAN(zblkldl)(int& nsuper, int* xsuper, int *pnode, int* xlindx,
                      int* lindx,  int* xlnz, complex<double> *lnz, int& defblk,
                      int &asdef,  int& numZEM, int& lbdef, int *def,
                      double& tol,
                      int *iprow, int* ipcol, int& tmpsiz, complex<double> *temp,
                      int& iwsize, int* iwork, int& rwsize, complex<double> *rwork,
                      int& iflag);

void _FORTRAN(blkslv)(const int &nsuper, const int* xsuper, const int* xlindx, const int *lindx,
                      const int* xlnz,
                      const double *lnx, const int& defblk, const int& numZEM, const int& lbdef,
                      const int* def, const int* iprow, const int* ipcol, const int* perm, const int* invp,
                      const double *rhs, double *sol, double *temp);

void _FORTRAN(zblkslv)(const int &nsuper, const int* xsuper, const int* xlindx, const int *lindx,
                      const int* xlnz,
                      const complex<double> *lnx, const int& defblk, const int& numZEM, const int& lbdef,
                      const int* def, const int* iprow, const int* ipcol, const int* perm, const int* invp,
                      const complex<double> *rhs, complex<double> *sol, complex<double> *temp);

void _FORTRAN(zblkslv2)(const int &nsuper, const int* xsuper, const int* xlindx, const int *lindx,
                      const int* xlnz,
                      const DComplex *lnx, const int& defblk, const int& numZEM, const int& lbdef,
                      const int* def, const int* iprow, const int* ipcol, const int* perm, const int* invp,
                      const DComplex *r1, const DComplex *r2,
                      DComplex *s1, DComplex *s2,
                      DComplex *t1, DComplex *t2);

void _FORTRAN(zblkslv3)(const int &nsuper, const int* xsuper, const int* xlindx, const int *lindx,
                      const int* xlnz,
                      const DComplex *lnx, const int& defblk, const int& numZEM, const int& lbdef,
                      const int* def, const int* iprow, const int* ipcol, const int* perm, const int* invp,
                      const DComplex *r1, const DComplex *r2, const DComplex *r3,
                      DComplex *s1, DComplex *s2, DComplex *s3,
                      DComplex *t1, DComplex *t2, DComplex *t3);

void _FORTRAN(zblkslv4)(const int &nsuper, const int* xsuper, const int* xlindx, const int *lindx,
                      const int* xlnz,
                      const DComplex *lnx, const int& defblk, const int& numZEM, const int& lbdef,
                      const int* def, const int* iprow, const int* ipcol, const int* perm, const int* invp,
                      const DComplex *r1, const DComplex *r2, const DComplex *r3, const DComplex *r4,
                      DComplex *s1, DComplex *s2, DComplex *s3, DComplex *s4,
                      DComplex *t1, DComplex *t2, DComplex *t3, DComplex *t4);

void _FORTRAN(blkslv2)(const int &nsuper, const int* xsuper, const int* xlindx, const int *lindx,
                      const int* xlnz,
                      const double *lnx, const int& defblk, const int& numZEM, const int& lbdef,
                      const int* def, const int* iprow, const int* ipcol, const int* perm, const int* invp,
                      const double *r1, const double *r2,
                      double *s1, double *s2,
                      double *t1, double *t2);

void _FORTRAN(blkslv3)(const int &nsuper, const int* xsuper, const int* xlindx, const int *lindx,
                      const int* xlnz,
                      const double *lnx, const int& defblk, const int& numZEM, const int& lbdef,
                      const int* def, const int* iprow, const int* ipcol, const int* perm, const int* invp,
                      const double *r1, const double *r2, const double *r3,
                      double *s1, double *s2, double *s3,
                      double *t1, double *t2, double *t3);

void _FORTRAN(blkslv4)(const int &nsuper, const int* xsuper, const int* xlindx, const int *lindx,
                      const int* xlnz,
                      const double *lnx, const int& defblk, const int& numZEM, const int& lbdef,
                      const int* def, const int* iprow, const int* ipcol, const int* perm, const int* invp,
                      const double *r1, const double *r2, const double *r3, const double *r4,
                      double *s1, double *s2, double *s3, double *s4,
                      double *t1, double *t2, double *t3, double *t4);

void _FORTRAN(zeromat)(int &numUncon, int &nsuper, int *xsuper,
                       int *xlnz, double *lnz, int *invsuper);

void _FORTRAN(zzeromat)(int &numUncon, int &nsuper, int *xsuper,
                       int *xlnz, complex<double> *lnz, int *invsuper);

void _FORTRAN(blkns)(int &nsuper, int *xsuper, int *xlindx, int *lindx,
                     int *xlnz, double *lnz, int &defblk, int &nrbm,
                     int &lbdef, int *def, int *ipcol, int *invp, double *ns,
                     int &numUncon, double *temp);

void _FORTRAN(zblkns)(int &nsuper, int *xsuper, int *xlindx, int *lindx,
                     int *xlnz, complex<double> *lnz, int &defblk, int &nrbm,
                     int &lbdef, int *def, int *ipcol, int *invp, complex<double> *ns,
                     int &numUncon, complex<double> *temp);

/*
// Not Used
void _FORTRAN(addmat)(int &numUncon, int *rowidx, double *values, int *perm,
                      int *invp, int &nsuper, int *xsuper, int *xlindx,
                      int *lindx , int *xlnz, double *lnz, int *offset,
                      int *invsuper);

void _FORTRAN(addone)(int &rowidx, int &colidx, double &value,
                      int *invp,   int *xsuper,   int *xlindx,
                      int *lindx,  int *xlnz,   double *lnz, int *invsuper);
*/

void _FORTRAN(blkslvp)(const int &nsuper, const int* xsuper, const int* xlindx, const int *lindx, const int* xlnz,
                      const double *lnx, const int& defblk, const int& numZEM, const int& lbdef,
                      const int* def, const int* iprow, const int* ipcol, const int* perm, const int* invp,
                      const int &nrhs, const double *rhs, const int &ldr, double *sol, double *temp);

void _FORTRAN(zblkslvp)(const int &nsuper, const int* xsuper, const int* xlindx, const int *lindx, const int* xlnz,
                      const DComplex *lnx, const int& defblk, const int& numZEM, const int& lbdef,
                      const int* def, const int* iprow, const int* ipcol, const int* perm, const int* invp,
                      const int &nrhs, const DComplex *rhs, const int &ldr, DComplex *sol, DComplex *temp);
}

inline void Tblkslvp(const int &a, const int* b, const int* c, const int *d, const int* e,
                      const DComplex *f, const int& g, const int& h, const int& i,
                      const int* j, const int* k, const int* l, const int* m, const int* n,
                      const int &o, const DComplex *p, const int &q, DComplex *r, DComplex *s)
{ _FORTRAN(zblkslvp)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s); }

inline void Tblkslvp(const int &a, const int* b, const int* c, const int *d, const int* e,
                     const double *f, const int& g, const int& h, const int& i,
                     const int* j, const int* k, const int* l, const int* m, const int* n,
                     const int &o, const double *p, const int &q, double *r, double *s)
{ _FORTRAN(blkslvp)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s); }

inline void Tblkldl(int& a, int* b, int *c, int* d,
                    int* e,  int* f, double *g, int& h,
                    int &i,  int& j, int& k, int *l, double& m,
                    int *n, int* o, int& p, double *q,
                    int& r, int* s, int& t, double *u, int& v)
{
_FORTRAN(blkldl)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v);
}

inline void Tblkldl(int& a, int* b, int *c, int* d,
                    int* e,  int* f, complex<double> *g, int& h,
                    int &i,  int& j, int& k, int *l, double& m,
                    int *n, int* o, int& p, complex<double> *q,
                    int& r, int* s, int& t, complex<double> *u,
                    int& v)
{
_FORTRAN(zblkldl)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v);
}

inline void Tblkslv(const int &a, const int* b, const int* c, const int *d, const int* e,
                    const double *f, const int& g, const int& h, const int& i,
                    const int* j, const int* k, const int* l, const int* m, const int* n,
                    const double *o, double *p, double *q)
{
_FORTRAN(blkslv)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q);
}

inline void Tblkslv(const int &a, const int* b, const int* c, const int *d, const int* e,
                    const complex<double> *f, const int& g, const int& h, const int& i,
                    const int* j, const int* k, const int* l, const int* m, const int* n,
                    const complex<double> *o, complex<double> *p, complex<double> *q)
{
_FORTRAN(zblkslv)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q);
}


inline void Tblkslv2(const int &a, const int* b, const int* c, const int *d, const int* e,
                    const complex<double> *f, const int& g, const int& h, const int& i,
                    const int* j, const int* k, const int* l, const int* m, const int* n,
                    const DComplex *o1, const DComplex *o2,
                    DComplex *p1, DComplex *p2,
                    DComplex *q1, DComplex *q2)
{
_FORTRAN(zblkslv2)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o1,o2,p1,p2,q1,q2);
}

inline void Tblkslv3(const int &a, const int* b, const int* c, const int *d, const int* e,
                    const complex<double> *f, const int& g, const int& h, const int& i,
                    const int* j, const int* k, const int* l, const int* m, const int* n,
                    const DComplex *o1, const DComplex *o2, const DComplex *o3,
                    DComplex *p1, DComplex *p2, DComplex *p3,
                    DComplex *q1, DComplex *q2, DComplex *q3)
{
_FORTRAN(zblkslv3)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o1,o2,o3,p1,p2,p3,q1,q2,q3);
}

inline void Tblkslv4(const int &a, const int* b, const int* c, const int *d, const int* e,
                    const complex<double> *f, const int& g, const int& h, const int& i,
                    const int* j, const int* k, const int* l, const int* m, const int* n,
                    const DComplex *o1, const DComplex *o2, const DComplex *o3, const DComplex *o4,
                    DComplex *p1, DComplex *p2, DComplex *p3, DComplex *p4,
                    DComplex *q1, DComplex *q2, DComplex *q3, DComplex *q4)
{
_FORTRAN(zblkslv4)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o1,o2,o3,o4,p1,p2,p3,p4,q1,q2,q3,q4);
}

inline void Tblkslv2(const int &a, const int* b, const int* c, const int *d, const int* e,
                    const double *f, const int& g, const int& h, const int& i,
                    const int* j, const int* k, const int* l, const int* m, const int* n,
                    const double *o1, const double *o2,
                    double *p1, double *p2,
                    double *q1, double *q2)
{
_FORTRAN(blkslv2)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o1,o2,p1,p2,q1,q2);
}

inline void Tblkslv3(const int &a, const int* b, const int* c, const int *d, const int* e,
                    const double *f, const int& g, const int& h, const int& i,
                    const int* j, const int* k, const int* l, const int* m, const int* n,
                    const double *o1, const double *o2, const double *o3,
                    double *p1, double *p2, double *p3,
                    double *q1, double *q2, double *q3)
{
_FORTRAN(blkslv3)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o1,o2,o3,p1,p2,p3,q1,q2,q3);
}

inline void Tblkslv4(const int &a, const int* b, const int* c, const int *d, const int* e,
                    const double *f, const int& g, const int& h, const int& i,
                    const int* j, const int* k, const int* l, const int* m, const int* n,
                    const double *o1, const double *o2, const double *o3, const double *o4,
                    double *p1, double *p2, double *p3, double *p4,
                    double *q1, double *q2, double *q3, double *q4)
{
_FORTRAN(blkslv4)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o1,o2,o3,o4,p1,p2,p3,p4,q1,q2,q3,q4);
}



inline void Tzeromat(int &a, int &b, int *c,
                     int *d, double *e, int *f)
{
_FORTRAN(zeromat)(a,b,c,d,e,f);
}

inline void Tzeromat(int &a, int &b, int *c,
                     int *d, complex<double> *e, int *f)
{
_FORTRAN(zzeromat)(a,b,c,d,e,f);
}

inline void Tblkns(int &a, int *b, int *c, int *d,
                   int *e, double *f, int &g, int &h,
                   int &i, int *j, int *k, int *l, double *m,
                   int &n, double *o)
{
_FORTRAN(blkns)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o);
}

inline void Tblkns(int &a, int *b, int *c, int *d,
                   int *e, complex<double> *f, int &g, int &h,
                   int &i, int *j, int *k, int *l, complex<double> *m,
                   int &n, complex<double> *o)
{
_FORTRAN(zblkns)(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o);
}

template<class Scalar>
GenBLKSparseMatrix<Scalar>::~GenBLKSparseMatrix()
{
 if(def) { delete [] def; def=0; }
 if(lindx) { delete [] lindx; lindx    = 0; }
 if(xlindx) { delete [] xlindx; xlindx   = 0; }
 if(snode) { delete [] snode; snode    = 0; }
 if(xsuper) { delete [] xsuper; xsuper   = 0; }
 if(xlnz) { delete [] xlnz; xlnz     = 0; }
 if(perm) { delete [] perm; perm     = 0; }
 if(invp) { delete [] invp; invp     = 0; }
 if(invsuper) { delete [] invsuper; invsuper=0; }
 if(lnz) { delete [] lnz; lnz=0; }
 if(iwork) { delete [] iwork; iwork = 0; }
 if(myRbm && rbm) { delete rbm; rbm = 0; }
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::init()
{
  tol     = 1.0E-6;
  rbm     = 0;
  numrbm  = 0;
  ngrbm   = 0;
  def     = 0;
  lindx   = 0;
  xlindx  = 0;
  snode   = 0;
  xsuper  = 0;
  xlnz    = 0;
  perm    = 0;
  invp    = 0;
  invsuper= 0;
  lnz     = 0;
  iwork   = 0;
  lbdef   = 0;
  myRbm   = false;
}

template<class Scalar>
GenBLKSparseMatrix<Scalar>::GenBLKSparseMatrix(const Connectivity *cn, const DofSetArray *_dsa,
                                               const DofSetArray *c_dsa, double _tol,
                                               const SolverCntl &_scntl, Rbm *_rbm) :
 SparseData(_dsa,c_dsa,cn,1), scntl(_scntl)
{
  init();
  rbm   = _rbm;
  tol   = _tol;

  // Geometric RBM
  if(rbm) ngrbm  = rbm->numRBM();

  allocateMemory();
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::factor()
{
  if(numUncon == 0) return;

  // reallocate this memory if we are doing multiple factor/solves like
  // in nonlinear statics or nonlinear dynamics.

  if(iwork==0) iwork = new int[7*numUncon+3];

  int iflag; // Error flag

//       ***************************************************
//       Numerical input into data structure for sparse LDL'
//       factorization.
//       ***************************************************
//
//       --------------------------------------------------------
//       INPNV ...   input numerical values into data structures.
//
//       Input:      N, COLPTR, ROWIDX, VALUES, PERM, INVP,
//                   NSUPER, XSUPER, XLINDX, LINDX, XLNZ,
//       Output:     LNZ
//       Work:       IWORK(N)
//       --------------------------------------------------------

//       ************************
//       Numerical factorization.
//       ************************
//       ---------------------------------------------------
//       BFINIT ...  initialization for block factorization.
//
//       Input:      NSUPER, XSUPER, SNODE, XLINDX, LINDX
//       Output:     TMPSIZ, RWSIZE
//       ---------------------------------------------------

  int rwsize;

  _FORTRAN(bfinit)(nsuper, xsuper, snode, xlindx,
                   lindx,    tmpsiz, rwsize);

  // TESTING:
  tmpsiz = tmpsiz * 2;

  Scalar *tmpvec = new Scalar[tmpsiz];

  Scalar *rwork  = new Scalar[rwsize];

  iwsiz  =  3 * numUncon + 2 * nsuper;

//       -------------------------------------------------------
//       BLKLDL ...  numerical factorization.
//
//       Input:      NSUPER, XSUPER, SNODE, XLINDX, LINDX, XLNZ,
//                   LNZ, DEFBLK, TOL, TMPSIZ, IWSIZE, RWSIZE
//       Output:     LNZ, NDEF, LBDEF, DEF, IPROW, IPCOL, IFLAG
//       Work:       TMPVEC(TMPSIZ), IWORK(2*N+2*NSUPER),
//                   RWORK(RWSIZE)
//       -------------------------------------------------------


  // if ngrbm > 0, lbdef = number of rbm in last block

  // perform full pivoting on last defblk x defblk block of sparse matrix
  // def = integer array, dimension numUncon, it identifies
  //       the columns of the matrix that are linearly dependent.
  int* deftemp = new int[numUncon];

  iprow = 0;
  ipcol = 0;
  int asdef = 0;
  if(defblk > 0) {
    iprow = new int[defblk]; // contains row pivoting sequence for last block.
    ipcol = new int[defblk]; // contains col pivoting sequence for last block.
    asdef = ngrbm; // assumed deficiency of the last block.
  }
  Tblkldl(nsuper, xsuper, snode, xlindx, lindx,
          xlnz, lnz, defblk, asdef, numrbm, lbdef,
          deftemp, tol, iprow, ipcol, tmpsiz,
          tmpvec, iwsiz, iwork, rwsize, rwork, iflag);

  // IFLAG =  0: successful factorization.
  // IFLAG = 31: insufficient work space in tmpvec.
  // IFLAG = 32: insufficient work space in iwork.
  // IFLAG = 33: insufficient work space in rwork.
  if(iflag != 0)
    throw std::runtime_error("Error during sparse factor");

  if(def) delete [] def; def = 0;
  if(numrbm > 0) {
    def = new int[numrbm];
    for(int i=0; i<numrbm; ++i) def[i] = deftemp[i];
  }

  if(this->print_nullity && numrbm > 0) {
     std::cerr << " ... Matrix is singular: size = " << numUncon << ", rank = " << numUncon-numrbm << ", nullity = " << numrbm;
     if(rbm) std::cerr << " (" << ngrbm << " grbm/hzem + " << numrbm-ngrbm << " other)";
     std::cerr << " ...\n";
  }

  delete [] rwork;
  delete [] tmpvec;
  delete [] deftemp;

  computeRBMs();
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::computeRBMs()
{
  if((ngrbm != numrbm) && (numrbm > 0)) {
    //filePrint(stderr," ... Computing %d Sparse RBM(s), tolerance = %e\n",numrbm,tol);

    // compute rigid body modes
    Scalar *ns  = new Scalar[numrbm*numUncon];
    Scalar *tempvec = new Scalar[numUncon];
    Tblkns(nsuper, xsuper, xlindx,  lindx,
           xlnz,    lnz, defblk, numrbm,
           lbdef,    def,  ipcol,   invp,
           ns,   numUncon, tempvec);

    GenVector<Scalar> *zem = new GenVector<Scalar>[numrbm];
    GenVector<Scalar> v(numUncon,0.0);
    for(int m = 0; m < numrbm; ++m) {
      zem[m] = v;
      for(int i=0; i<numUncon; ++i)
        zem[m][i] = ns[i+m*numUncon];
    }
    if(rbm && myRbm) delete rbm;
    rbm = new Rbm(zem,numrbm,numUncon,1);
    myRbm = true;
    delete [] ns; delete [] tempvec;
  }
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::getNullSpace(Scalar *ns)
{
  if(numrbm > 0) {
    Scalar *tempvec = new Scalar[numUncon];
    Tblkns(nsuper, xsuper, xlindx,  lindx,
           xlnz,    lnz, defblk, numrbm,
           lbdef,    def,  ipcol,   invp,
           ns,   numUncon, tempvec);
    delete [] tempvec;
  }
}

template<class Scalar>
double
GenBLKSparseMatrix<Scalar>::getMemoryUsed() const
{
 // Figure this out later!
 return 0.0;
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::solve(const Scalar *rhs, Scalar *solution)
{
 this->solveTime -= getTime();

 Scalar *temp = new Scalar[numUncon];

//       ********************
//       Triangular solution.
//       ********************
//
//       --------------------------------------------------------------
//       BLKSLV ...  numerical forward/backward solution.
//
//       Input:      nsuper, xsuper, xlindx, lindx, xlnz, lnz, defblk,
//                   numrbm, lbdef, def, iprow, ipcol, perm, invp, solution
//       Output:     solution
//       Work:       temp[numUncon]
//       --------------------------------------------------------------

 Tblkslv(nsuper, xsuper, xlindx, lindx, xlnz,
         lnz, defblk, numrbm, lbdef, def,
         iprow,  ipcol, perm, invp, rhs,
         solution, temp);

 delete [] temp;

 this->solveTime += getTime();
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::solve(const GenVector<Scalar> &rhs, GenVector<Scalar> &solution)
{
 solve((const Scalar *)rhs.data(), solution.data());
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::reSolve(Scalar *rhs)
{
 if(numUncon==0) return;

 this->solveTime -= getTime();

 Scalar *temp     = new Scalar[numUncon+1];
 Scalar *solution = new Scalar[numUncon];

 Tblkslv(nsuper, xsuper, xlindx, lindx, xlnz, lnz, defblk, numrbm, lbdef, def,
         iprow, ipcol, perm, invp, rhs, solution, temp);

 for(int i=0; i < numUncon; i++) rhs[i] = solution[i];

 delete [] temp;
 delete [] solution;

 this->solveTime += getTime();
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::reSolve(GenVector<Scalar> &rhs)
{
 reSolve(rhs.data());
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::reSolve(int nRHS, Scalar **RHS)
{
#ifdef MIN_MEMORY
 for(int n=0; n<nRHS; ++n) {
   reSolve(RHS[n]);
}
#else
 if(numUncon==0) return;
 this->solveTime -= getTime();
 Scalar *s = new Scalar[numUncon*nRHS];
 Scalar *t = new Scalar[numUncon*nRHS];
 Scalar *a = new Scalar[numUncon*nRHS];

 for(int j=0; j<nRHS; ++j)
   for(int i=0; i<numUncon; ++i) a[i+j*numUncon] = RHS[j][i];

 Tblkslvp(nsuper, xsuper, xlindx, lindx, xlnz, lnz, defblk, numrbm, lbdef,  def,
          iprow, ipcol, perm, invp, nRHS, a, numUncon, s, t);

 for(int j=0; j<nRHS; ++j)
   for(int i=0; i<numUncon; ++i) RHS[j][i] = s[i+j*numUncon];

 delete [] s; delete [] t; delete [] a;
 this->solveTime += getTime();
#endif
}

/*
template<>
void
GenBLKSparseMatrix<double>::reSolve(int nRHS, double **RHS)
{
 if(numUncon==0) return;
 this->solveTime -= getTime();
 for(int n=0; n<nRHS; ++n) reSolve(RHS[n]);
 this->solveTime += getTime();
}


template<>
void
GenBLKSparseMatrix<DComplex>::reSolve(int nRHS, DComplex **RHS) {

 this->solveTime -= getTime();

 int i = 0;
 int multiple = 4;
 int j;

 if(nRHS <= 4)
   multiple = nRHS;
// RT: new would be preferred, but causes crash on delete from leak somewhere
// else; if I ever have time, I will try to track it down
// DComplex *t1 = new DComplex[8*numUncon];
 DComplex *t1 = (DComplex*)alloca(8*numUncon*sizeof(DComplex));
 DComplex *t2 = t1 + numUncon;
 DComplex *t3 = t1 + 2*numUncon;
 DComplex *t4 = t1 + 3*numUncon;

 DComplex *s1 = t1 + 4*numUncon;
 DComplex *s2 = t1 + 5*numUncon;
 DComplex *s3 = t1 + 6*numUncon;
 DComplex *s4 = t1 + 7*numUncon;

 switch (multiple) {
   default:
   case 4:
     {
       // zblkslv4 = forward/backward routine for 4 RHS vectors
       for( ; i < nRHS-3; i += 4) {
         _FORTRAN(zblkslv4)(nsuper, xsuper, xlindx, lindx, xlnz,
                   lnz, defblk, numrbm, lbdef,  def,
                   iprow,  ipcol,   perm,  invp,  RHS[i], RHS[i+1],
                   RHS[i+2], RHS[i+3], s1, s2, s3, s4, t1, t2, t3, t4);
         for (j=0; j<numUncon; j++) {
            RHS[i][j] = s1[j];
            RHS[i+1][j] = s2[j];
            RHS[i+2][j] = s3[j];
            RHS[i+3][j] = s4[j];
         }
       }
     }
   case 3:
     {
       // zblkslv3 = forward/backward routine for 3 RHS vectors
       for( ; i < nRHS-2; i += 3) {
         _FORTRAN(zblkslv3)(nsuper, xsuper, xlindx, lindx, xlnz,
                   lnz, defblk, numrbm, lbdef,  def,
                   iprow,  ipcol,   perm,  invp,  RHS[i], RHS[i+1],
                   RHS[i+2], s1, s2, s3, t1, t2, t3);
         for (j=0; j<numUncon; j++) {
            RHS[i][j] = s1[j];
            RHS[i+1][j] = s2[j];
            RHS[i+2][j] = s3[j];
         }
       }
     }
   case 2:
     {
       // zblkslv2 = forward/backward routine for 2 RHS vectors
       for( ; i < nRHS-1; i += 2) {
         _FORTRAN(zblkslv2)(nsuper, xsuper, xlindx, lindx, xlnz,
                   lnz, defblk, numrbm, lbdef,  def,
                   iprow,  ipcol,   perm,  invp,  RHS[i], RHS[i+1],
                   s1, s2, t1, t2);
         for (j=0; j<numUncon; j++) {
            RHS[i][j] = s1[j];
            RHS[i+1][j] = s2[j];
         }
       }
     }
   case 1:
     {
       // zblkslv = forward/backward routine for 1 RHS vectors
       for ( ; i < nRHS; ++i) {
         _FORTRAN(zblkslv)(nsuper, xsuper, xlindx, lindx, xlnz,
                   lnz, defblk, numrbm, lbdef,  def,
                   iprow,  ipcol,   perm,  invp,  RHS[i],
                   s1,   t1);
         for (j=0; j<numUncon; j++)
            RHS[i][j] = s1[j];
       }
     }
     break;
 }

// delete [] t1;

 this->solveTime += getTime();

}
*/

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::reSolve(int nRHS, GenVector<Scalar> *RHS)
{
#ifdef MIN_MEMORY
 for(int n=0; n<nRHS; ++n) reSolve(RHS[n].data());
#else
 if(numUncon==0) return;
 this->solveTime -= getTime();
 Scalar *s = new Scalar[numUncon*nRHS];
 Scalar *t = new Scalar[numUncon*nRHS];
 Scalar *a = new Scalar[numUncon*nRHS];

 for(int j=0; j<nRHS; ++j)
   for(int i=0; i<numUncon; ++i) a[i+j*numUncon] = RHS[j][i];

 Tblkslvp(nsuper, xsuper, xlindx, lindx, xlnz, lnz, defblk, numrbm, lbdef,  def,
          iprow, ipcol, perm, invp, nRHS, a, numUncon, s, t);

 for(int j=0; j<nRHS; ++j)
   for(int i=0; i<numUncon; ++i) RHS[j][i] = s[i+j*numUncon];

 delete [] s; delete [] t; delete [] a;
 this->solveTime += getTime();
#endif
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::zeroAll()
{
  numrbm=0;
  lbdef=0;
  if(numUncon > 0) {
    for(int i=0; i < xlnz[numUncon]; ++i)
      ScalarTypes::initScalar(lnz[i], 0.0, 0.0);
  }
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::clean_up()
{
 if(lnz) {
   delete [] lnz;
   lnz = 0;
 }
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::unify(FSCommunicator *communicator)
{
#ifdef DISTRIBUTED
 communicator->globalSum(xlnz[numUncon], lnz);
#endif
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::add(const FullSquareMatrix &kel, const int *dofs)
{
 if(numUncon == 0) return;

 int i, j, k, rowi, colj, offset, p1, position, csuper, fstcol, lxbeg, lxend;

 int kndof = kel.dim();
 for(i = 0; i < kndof; ++i) {
   if((rowi = unconstrNum[dofs[i]]) == -1) continue;
   p1     = invp[rowi] - 1;
   position = xlnz[p1+1];
   csuper = invsuper[p1] - 1;
   fstcol = xsuper[csuper] - 1;
   lxbeg  = xlindx[csuper]-1;
   lxend  = xlindx[csuper+1]-1;
   for(j = 0; j < kndof; ++j) {
     if((colj = unconstrNum[dofs[j]]) == -1) continue;
     int irow   = invp[colj] - 1;
     if(irow >= fstcol) {
       offset = lxend - lxbeg;
       for(k=lxbeg; k<lxend; ++k) {
         offset -= 1;
         if(lindx[k]-1 == irow) {
           lnz[position - 2 - offset] += kel[i][j];
           break;
         }
       }
     }
   }
 }
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::add(const GenAssembledFullM<Scalar> &kel, const int *dofs)
{
 if(numUncon == 0) return;

 int i, j, k, rowi, colj, offset, p1, position, csuper, fstcol, lxbeg, lxend;

 int kndof = kel.dim();
 for( i = 0; i < kndof; ++i ) {
   if( (rowi = dofs[i]) == -1 ) continue;
   p1     = invp[rowi] - 1;
   position = xlnz[p1+1];
   csuper = invsuper[p1] - 1;
   fstcol = xsuper[csuper] - 1;
   lxbeg  = xlindx[csuper]-1;
   lxend  = xlindx[csuper+1]-1;
   for( j = 0; j < kndof; ++j ) {
     if( (colj = dofs[j]) == -1 ) continue;
     int irow   = invp[colj] - 1;
     if ( irow >= fstcol ) {
       offset = lxend - lxbeg;
       for(k=lxbeg; k<lxend; ++k) {
         offset -= 1;
         if(lindx[k]-1 == irow) {
           lnz[ position - 2 - offset] += kel[i][j];
           break;
         }
       }
     }
   }
 }
}

// This routine has to be converted to input to the LNZ Scalar
// precision array instead of the temporary unonz Scalar precision
// array.

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::add(const FullM &knd, int fRow, int fCol)
{
  int offset, k, p1, position, csuper,fstcol, lxbeg, lxend;
  int iCol, iRow;
  int nrow = knd.numRow(); // number of rows
  int ncol = knd.numCol(); // number of columns

  for(iRow = 0; iRow < nrow; ++iRow) {
    int rowi = fRow + iRow;
    p1     = invp[rowi] - 1;
    position = xlnz[p1+1];
    csuper = invsuper[p1] - 1;
    fstcol = xsuper[csuper] - 1;
    lxbeg  = xlindx[csuper]-1;
    lxend  = xlindx[csuper+1]-1;
    for(iCol = 0; iCol < ncol; ++iCol) {
      int colj = fCol + iCol ;
      int irow   = invp[colj] - 1;
      if(irow >= fstcol) {
        offset = lxend - lxbeg;
        for(k=lxbeg; k<lxend; ++k) {
          offset -= 1;
          if(lindx[k]-1 == irow) {
            lnz[ position - 2 - offset] += knd[iRow][iCol];
            break;
          }
        }
      }
    }
  }
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::add(Scalar *_lnz)
{
  int i;
  for(i=0; i < xlnz[numUncon]; ++i)
    lnz[i] += _lnz[i];
}

// Assembly from a Boeing Sparse data structure in fortran
template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::addBoeing(int nl, const int *Kai, const int *Kaj,
                                      const double *nz, const int *map, Scalar multiplier)
{
 int i, j, k, offset, p1,position,csuper,fstcol,lxbeg,lxend;

 // FIRST pass to fill upper half of sparse matrix.
 for(i = 0; i < nl; ++i) {
   if(map[i] >= neq)
      { fprintf(stderr, "Out of bounds, %d %d\n",map[i],numUncon); return ;}
   if(map[i] == -1) continue;
   int rowi = unconstrNum[map[i]];
   if(rowi < 0) continue;
   p1       = invp[rowi] - 1;
   position = xlnz[p1+1];
   csuper   = invsuper[p1] - 1;
   fstcol   = xsuper[csuper] - 1;
   lxbeg    = xlindx[csuper]-1;
   lxend    = xlindx[csuper+1]-1;
   for(j = Kai[i]; j < Kai[i+1]; ++j) {
     if(map[Kaj[j-1]-1] == -1) continue;
     int rowj = (unconstrNum.size() != 0) ? unconstrNum[map[Kaj[j-1]-1]] : map[Kaj[j-1]-1];
     if(rowj < 0) continue;
     int irow   = invp[rowj] - 1;
     if(irow >= fstcol) {
       offset = lxend - lxbeg;
       for(k=lxbeg; k<lxend; ++k) {
         offset -= 1;
         if(lindx[k]-1 == irow) {
           lnz[ position - 2 - offset] += (nz[j-1]*multiplier);
           break;
         }
       }
     }
   }
 }

 // SECOND pass to fill lower half of sparse matrix.
 for(i = 0; i < nl; ++i) {
   if(map[i] >= neq)
      { fprintf(stderr, "Out of bounds, %d %d\n",map[i],numUncon); return ;}
   if(map[i] == -1) continue;
   int rowj = unconstrNum[map[i]];
   if(rowj < 0) continue;
   for(j = Kai[i]; j < Kai[i+1]; ++j) {
     if(map[Kaj[j-1]-1] == -1) continue;
     int rowi = unconstrNum[map[Kaj[j-1]-1]];
     if(rowi < 0) continue;
     if(rowi == rowj) continue; // skip the second pass on the diagonal entries
     p1       = invp[rowi] - 1;
     position = xlnz[p1+1];
     csuper   = invsuper[p1] - 1;
     fstcol   = xsuper[csuper] - 1;
     lxbeg    = xlindx[csuper]-1;
     lxend    = xlindx[csuper+1]-1;
     int irow   = invp[rowj] - 1;
     if(irow >= fstcol) {
       offset = lxend - lxbeg;
       for(k=lxbeg; k<lxend; ++k) {
         offset -= 1;
         if(lindx[k]-1 == irow) {
           lnz[ position - 2 - offset] += (nz[j-1]*multiplier);
           break;
         }
       }
     }
   }
 }
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::addone(Scalar d, int dofi, int dofj)
{
 // WARNING: this adds both [i,j] and [j,i]
 if((dofi < 0) || (dofj < 0)) return;
 int k, irow, rowi, colj, offset, p1,position,csuper,fstcol,lxbeg,lxend;
 if ( unconstrNum.size() != 0 ) {
   if((rowi = unconstrNum[dofi]) == -1 || (colj = unconstrNum[dofj]) == -1) return;
 }
 else { rowi = dofi; colj = dofj; }

 p1      = invp[rowi] - 1;
 position = xlnz[p1+1];
 csuper = invsuper[p1] - 1;
 fstcol = xsuper[csuper] - 1;
 lxbeg  = xlindx[csuper]-1;
 lxend  = xlindx[csuper+1]-1;
 irow   = invp[colj] - 1;
 if ( irow >= fstcol ) {
   offset = lxend - lxbeg;
   for(k=lxbeg; k<lxend; ++k) {
     offset -= 1;
     if(lindx[k]-1 == irow) {
        lnz[ position - 2 - offset] += d;
        break;
     }
   }
 }

 if (dofi == dofj) return;
 p1      = invp[colj] - 1;
 position = xlnz[p1+1];
 csuper = invsuper[p1] - 1;
 fstcol = xsuper[csuper] - 1;
 lxbeg  = xlindx[csuper]-1;
 lxend  = xlindx[csuper+1]-1;
 irow   = invp[rowi] - 1;
 if(irow >= fstcol) {
   offset = lxend - lxbeg;
   for(k=lxbeg; k<lxend; ++k) {
     offset -= 1;
     if(lindx[k]-1 == irow) {
       lnz[ position - 2 - offset] += d;
       break;
     }
   }
 }
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::mult(const Scalar *rhs, Scalar *result) const
{
  std::cerr << " *** WARNING: GenBLKSparseMatrix<Scalar>::mult(const Scalar *rhs, Scalar *result) is not implemented \n";
}

template<class Scalar>
Scalar
GenBLKSparseMatrix<Scalar>::getone(int dofi, int dofj)
{
 Scalar d; // return value

 int k, irow, rowi, colj, offset, p1,position,csuper,fstcol,lxbeg,lxend;
 if(unconstrNum.size() != 0) {
   if((rowi = unconstrNum[dofi]) == -1 || (colj = unconstrNum[dofj]) == -1) return 0.0;
 }
 else { rowi = dofi; colj = dofj; }

 p1      = invp[rowi] - 1;
 position = xlnz[p1+1];
 csuper = invsuper[p1] - 1;
 fstcol = xsuper[csuper] - 1;
 lxbeg  = xlindx[csuper]-1;
 lxend  = xlindx[csuper+1]-1;
 irow   = invp[colj] - 1;
 if(irow >= fstcol) {
   offset = lxend - lxbeg;
   for(k=lxbeg; k<lxend; ++k) {
     offset -= 1;
     if(lindx[k]-1 == irow) {
        d = lnz[position - 2 - offset];
        return d;
     }
   }
 }

 p1      = invp[colj] - 1;
 position = xlnz[p1+1];
 csuper = invsuper[p1] - 1;
 fstcol = xsuper[csuper] - 1;
 lxbeg  = xlindx[csuper]-1;
 lxend  = xlindx[csuper+1]-1;
 irow   = invp[rowi] - 1;
 if(irow >= fstcol) {
   offset = lxend - lxbeg;
   for(k=lxbeg; k<lxend; ++k) {
     offset -= 1;
     if(lindx[k]-1 == irow) {
       d = lnz[position - 2 - offset];
       return d;
     }
   }
 }
 return 0.0;
}


template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::print()
{
 fprintf(stdout,"--- BLK MATRIX --- numUncon %d xlnz %d\n",
         numUncon, xlnz[numUncon]);
 int i;
 for(i=0; i<xlnz[numUncon]-1; ++i)
   std::cerr << lnz[i] << "  ";
 std::cerr << std::endl;

 fprintf(stdout,"============\n"); fflush(stdout);
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::printAll()
{
  fprintf(stderr, " numUncon       = %d\n",numUncon);
  fprintf(stderr, " neq            = %d\n",neq);
  fprintf(stderr, " numConstrained = %d\n",numConstrained);
}

template<class Scalar>
Scalar
GenBLKSparseMatrix<Scalar>::diag(int dof) const
{
 int k,colj, offset, p1,position,csuper,fstcol,lxbeg,lxend;
 p1 = invp[dof]-1;
 position = xlnz[p1+1];
 csuper = invsuper[p1] - 1;
 fstcol = xsuper[csuper] - 1;
 lxbeg  = xlindx[csuper]-1;
 lxend  = xlindx[csuper+1]-1;
 colj = dof;
 int irow   = invp[colj] - 1;
 if ( irow >= fstcol ) {
   offset = lxend - lxbeg;
   for(k=lxbeg; k<lxend; ++k) {
      offset -= 1;
      if(lindx[k]-1 == irow)
            return lnz[ position - 2 - offset];
   }
 }
 return 0.0;
}

template<class Scalar>
Scalar &
GenBLKSparseMatrix<Scalar>::diag(int dof)
{
 int k,colj, offset, p1,position,csuper,fstcol,lxbeg,lxend;
 p1 = invp[dof]-1;
 position = xlnz[p1+1];
 csuper = invsuper[p1] - 1;
 fstcol = xsuper[csuper] - 1;
 lxbeg  = xlindx[csuper]-1;
 lxend  = xlindx[csuper+1]-1;
 colj = dof;
 int irow   = invp[colj] - 1;
 if ( irow >= fstcol ) {
   offset = lxend - lxbeg;
   for(k=lxbeg; k<lxend; ++k) {
      offset -= 1;
      if(lindx[k]-1 == irow)
            return lnz[ position - 2 - offset];
   }
 }
 throw "Diagonal does not exist";
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::addDiscreteMass(int dof, Scalar mass)
{
 if(dof < 0) return;
 int cdof;
 if(unconstrNum.size() != 0)
 	cdof = unconstrNum[dof]; // dof is now in unconstrained numbering
 else
 	cdof = dof;
 if(cdof < 0) return;

 int k,colj, offset, p1,position,csuper,fstcol,lxbeg,lxend;
 p1 = invp[cdof]-1;
 position = xlnz[p1+1];
 csuper = invsuper[p1] - 1;
 fstcol = xsuper[csuper] - 1;
 lxbeg  = xlindx[csuper]-1;
 lxend  = xlindx[csuper+1]-1;
 colj = cdof;
 int irow   = invp[colj] - 1;
 if ( irow >= fstcol ) {
   offset = lxend - lxbeg;
   for(k=lxbeg; k<lxend; ++k) {
      offset -= 1;
      if(lindx[k]-1 == irow) {
        lnz[ position - 2 - offset] += mass;
      }
   }
 }
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::getRBMs(double *rigidBodyModes)
{
  if(rbm) rbm->getRBMs(rigidBodyModes);
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::getRBMs(Vector *rigidBodyModes)
{
  if(rbm) rbm->getRBMs(rigidBodyModes);
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::getRBMs(VectorSet &rigidBodyModes)
{
 if(rbm) rbm->getRBMs(rigidBodyModes);
}

// Constructor Kii, FETI Subdomain preconditioner solver

template<class Scalar>
GenBLKSparseMatrix<Scalar>::GenBLKSparseMatrix(const Connectivity *cn, const DofSetArray *_dsa,
                                               int *glInternalMap, double _tol, const SolverCntl &_scntl) :
  SparseData(_dsa, glInternalMap, cn, 1), scntl(_scntl)
{
  init();
  tol    = _tol;
  allocateMemory();
}

// Constructor for GtG Solver (First level coarse problem solver in FETI)

template<class Scalar>
GenBLKSparseMatrix<Scalar>::GenBLKSparseMatrix(const Connectivity *cn, const EqNumberer *_dsa, double _tol,
                                               const SolverCntl &_scntl, int _ngrbm)
  : SparseData(cn,_dsa,_tol), scntl(_scntl)
{
  init();
  tol    = _tol;
  ngrbm = _ngrbm;
  allocateMemory();
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::allocateMemory()
{
  using std::min;
  this->solveTime = 0.0;
  if(numUncon == 0) return;

  double t1 = -getTime();

  snode  = new int[numUncon];
  xsuper = new int[numUncon+1];
  xlindx = new int[numUncon+1];
  xlnz   = new int[numUncon+1];
  iwork  = new int[7*numUncon+3];
  perm   = new int[numUncon];
  invp   = new int[numUncon];
  colcnt = new int[numUncon];
  lindx  = new int[xunonz[numUncon]];

  int lxsize = xunonz[numUncon];
  iwsiz = 7 * numUncon + 3;

  int nnza = xunonz[numUncon] - 1 - numUncon;

  int i,j,k,iflag;

  // maxsup = maximum number of columns in each supernode (parameter)
  int maxsup = scntl.sparse_maxsup;  // default is 100, but may need larger maxsup for big subdomains
                                     // set using sparse_maxsup parameter in fem input file

  adj = new int[nnza];

  for (k=0, j=0; k<numUncon; k++)
     for (i=xunonz[k]; i<=xunonz[k+1]-1; i++)
        if ( rowu[i-1] == k+1 );
        else  {
          adj[j++] = rowu[i-1];
        }

  xadj = new int[numUncon + 1];

  for (i=0; i<numUncon+1; i++)
     xadj[i] = xunonz[i] - i ;

	xunonz.clear(); xunonz.shrink_to_fit();
	rowu.clear();   rowu.shrink_to_fit();

//-------------------------------------------
// Copy matrix structure from (XADJ,ADJ) to
// (XLINDX,LINDX) (because matrix structure is
// destroyed by the minimum degree ordering
// subroutine).
// -------------------------------------------

  for (k=0; k<numUncon+1; k++)
     xlindx[k] = xadj[k];

  for (k=0; k<nnza; k++)
     lindx[k] = adj[k];

// -------------------------------------------------------
// STEP 1:   ORDMMD ...  multiple minimum degree ordering.
//
//           Input:      NUMUNCON, XLINDX, LINDX, IWSIZE
//           Output:     PERM, INVP, NSUB, IFLAG
//           Work:       IWORK(4*N)
// -------------------------------------------------------

  // iwsiz = 4 * numUncon;  // iwsiz should be the dimension of iwork (not changed)
#ifdef USE_METIS
    if(scntl.sparse_renum == 1) {
        std::vector<idx_t> metis_options(METIS_NOPTIONS, 0);
	    METIS_SetDefaultOptions(metis_options.data());
        metis_options[METIS_OPTION_NUMBERING] = 1;
        metis_options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
        metis_options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_EDGE;
        metis_options[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP2SIDED;
        metis_options[METIS_OPTION_COMPRESS] = 1;
        metis_options[METIS_OPTION_CCORDER ] = 1;
        metis_options[METIS_OPTION_PFACTOR] = 400;
        metis_options[METIS_OPTION_NSEPS] = 4;

        int res = callMetis(&numUncon, xlindx, lindx, nullptr, metis_options.data(), perm, invp);
        if(res != METIS_OK) {
            std::cerr << "Metis failed result: " << res << std::endl;
            exit(1);
        }
    }
  else {
#endif
    _FORTRAN(ordmmd2)(numUncon, xlindx, lindx, invp, perm,
                      iwsiz, iwork, nsub, iflag);

//               IFLAG =  0: successful ordering.
//               IFLAG = 11: insufficient work space in IWORK.

    if(iflag == 11)
      fprintf(stderr,"*** ERROR: not enough space for iwork temperary array\n");
#ifdef USE_METIS
  }
#endif

  // iwsiz = 7 * numUncon + 3;

  // Rank deficiency information
  // defblk = size of last block to perform full pivoting on
  if(ngrbm > 0) {
    defblk = min(numUncon-1,scntl.sparse_defblk); // note: sparse_defblk default is 10
    // restore original numbering in defblk
    // note: perm maps from new index to original index and invp maps from original index to new index
    for(int i = numUncon-defblk; i < numUncon; ++i) {
      int p = i;
      int q = invp[i]-1;
      int j = perm[p]-1;
      // swap positions of equations i and j
      invp[i] = p+1;
      invp[j] = q+1;
      perm[p] = i+1;
      perm[q] = j+1;
    }
  }
  else
    defblk = 0;

//      ***********************
//      Symbolic factorization.
//      ***********************
//      -------------------------------------------------------------
//      STEP 2:
//      SFINIT ...   symbolic factorization initialization, which
//                   computes supernode partition and storage
//                   requirements for symbolic factorization;
//                   new ordering is a postordering of the nodal
//                   elimination tree.
//
//       Input:      N, NADJ, XADJ, ADJ, PERM, INVP, MAXSUP, DEFBLK,
//                   IWSIZE
//       Output:     PERM, INVP, COLCNT, NNZL, NSUB, NSUPER, XSUPER,
//                   SNODE, IFLAG
//       Work:       IWORK(7*N+3) ... the max required any subroutine.
//       -------------------------------------------------------------

  _FORTRAN(sfinit)(numUncon, nnza,    xadj,   adj,    perm,
                   invp,     maxsup,  defblk, colcnt, nnzl,
                   nsub,     nsuper,  xsuper, snode,  iwsiz,
                   iwork,    iflag);
  // warning: if defblk outputted by sfinit is less that the input value
  //          then the grbm method can fail.
  // IFLAG = 0: successful symbolic factorization initialization
  // IFLAG = 1: insufficent work space in IWORK.
  if(iflag == 1) {
    fprintf(stderr," *** ERROR: Insufficient working storage for SFINIT\n");
    exit(1);
  }

//       ------------------------------------------------------
//       STEP 3:
//       SYMFCT ...  perform supernodal symbolic factorization.
//
//       Input:      N, NADJ, XADJ, ADJ, PERM, INVP, COLCNT,
//                   NSUPER, XSUPER, SNODE , NSUB, IWSIZE
//       Output:     XLINDX, LINDX, XLNZ, IFLAG
//       WORK:       IWORK(NSUPER+2*N+1)
//
//       No longer needed: ADJ, XADJ, COLCNT
//       ------------------------------------------------------

  // iwsiz = nsuper + 2 * numUncon + 1;  // iwsize is the dimension of iwork
  if(nsub > lxsize) { delete [] lindx; lindx = new int[nsub]; }

  _FORTRAN(symfct)(numUncon, nnza,   xadj,   adj,    perm,
                   invp,     colcnt, nsuper, xsuper, snode,
                   nsub,     xlindx, lindx, xlnz,    iwsiz,
                   iwork,    iflag);

  // IFLAG =  0: no error.
  // IFLAG = 22: insufficient work space in IWORK.
  // IFLAG = 23: inconsistancy in the input.

  // ERROR Checking
  if (iflag == 22) {
     fprintf(stderr," *** ERROR: (SYMFCT) Memory for iwork is too small\n");
     exit(1);
  } else if (iflag == 23) {
     fprintf(stderr," *** ERROR: (SYMFCT) Inconsistency in input\n");
     exit(1);
  }

  // Delete memory that is not used any longer
  delete [] xadj;   xadj   = 0;
  delete [] adj;    adj    = 0;
  delete [] colcnt; colcnt = 0;

  lnz = new Scalar[xlnz[numUncon]];
  if(lnz==0) {
    fprintf(stderr,"ERROR: Cannot allocate enough memory for Sparse Matrix\n");
    fprintf(stderr,"Trying to allocate %f Mb\n",
                    8.0*xlnz[numUncon]/(1024.0*1024.0));
    exit(-1);
  }

  // zero lnz and create invsuper array
  invsuper = new int[numUncon];

  Tzeromat(numUncon,nsuper,xsuper,xlnz,lnz,invsuper);
  t1 += getTime();
}

template<class Scalar>
int
GenBLKSparseMatrix<Scalar>::numRBM() const
{
  return numrbm; // return the total number of zems (both geometric and otherwise) same as SkyMatrix
}

template<class Scalar>
void
GenBLKSparseMatrix<Scalar>::reSolve(GenFullM<Scalar> *mat)
{
 this->solveTime -= getTime();

 int nRHS = mat->numCol();
 int length = mat->numRow();

 Scalar *s = new Scalar[length*nRHS];
 Scalar *t = new Scalar[length*nRHS];
 Scalar *a = new Scalar[length*nRHS];

 for(int j=0; j<nRHS; ++j)
   for(int i=0; i<length; ++i) a[i+j*length] = (*mat)[i][j];

 Tblkslvp(nsuper, xsuper, xlindx, lindx, xlnz, lnz, defblk, numrbm, lbdef,  def,
          iprow,  ipcol,   perm,  invp, nRHS, a, length, s, t);

 for(int j=0; j<nRHS; ++j)
   for(int i=0; i<length; ++i) (*mat)[i][j] = s[i+j*length];

 delete [] s; delete [] t; delete [] a;

 this->solveTime += getTime();
}

