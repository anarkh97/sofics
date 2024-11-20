#ifndef _COARSESET_H_
#define _COARSESET_H_

#include <Utils.d/dofset.h>
#include <Feti.d/FetiSub.h>

template <class Scalar> class GenFullM;
template <class Scalar> class GenSymFullMatrix;
class IntFullM;
template <class Scalar> class GenSubDomain;
template <class Scalar> class GenSparseMatrix;
template <class Scalar> class GenSolver;
template <class Scalar> class FSCommPattern;

template<class Scalar>
class GenCoarseSet {
 public:
   FetiSub<Scalar> *sd;
   int myNum;
   int numGs;
   Scalar *locGs;
   Scalar *locQGs;	// for Feti Static, locQGs = locGs
   int numNeighb;
   const int *neighbs;
   int *neighbNumRBMs;
   Scalar **neighbQGs;
   int *subOffset; // indicates the size of interface to each subdomains
   int gSize;	   // interfSize, leading dimension of locGs, locQGs ...
   Scalar *locFGs;     // for Feti Dynamics, locFGs == locQGs;
   int numBCs;           // # of column of locBCs. i.e. # of lambda
   IntFullM *locBCs;     // Only non-zero values are stored.
   Scalar   *locFBCs;    // For Feti2, Q is Fi
   Scalar   **neighbFBCs2;
   int *neighbNumCRNs;

   // To store neighbor's Gs and leading dimensions of those Gs
   Scalar **neighbGs;
   int *leadingDimGs;
   int isDynamic; // 0 = static 1 = Dynamic

   Scalar *locInterfCRNs; // working only.
   int    *neighbCSubIndex; // Pointer into neighbCs for the neighbor's Cs
   int    (*neighbCs)[3];   // vector of all the indices of the neighbor;s Cs
                    // first column is the corner's dof index in the interface.
                    // second column is the master sub number
                    // third column is the corner number in the master numbering

   // GtC
   GenFullM<Scalar> * coarseGtCMult(int i, Scalar *remG, int ldG);
   GenFullM<Scalar> * coarseGtCMult();

   // Nonlinear variables to reduce reallocations
   Scalar *pool1;
   Scalar *pool2;

   void clean_up();
   void getNeighbGs(GenCoarseSet<Scalar> *csets, GenSubDomain<Scalar> *sd, FSCommPattern<Scalar> *rbmPat);
   void getNeighbQGs(GenCoarseSet<Scalar> *csets, GenSubDomain<Scalar> *sd, int kc, 
                     FSCommPattern<Scalar> *rbmPat, GenSolver<Scalar> * = 0);
   // void getNeighbQGs(GenSubDomain<Scalar> *sd);
   void getNeighbFBCs(GenCoarseSet<Scalar> *csets, GenSubDomain<Scalar> *sd, GenSolver<Scalar> * = 0);
   void computeNeighbFBCs(GenSubDomain<Scalar> *sd, GenSolver<Scalar> *solver);
   void buildBNeighbCs(GenCoarseSet<Scalar> *csets);

   int getNumCs(int sub);
   void getCOffset(int sub, int (*offsets)[3]);

   void addQContrib(GenSparseMatrix<Scalar> *GtQG, int iSub,EqNumberer *, int hasCorners=0);
   void addGContrib(GenSparseMatrix<Scalar> *GtQG, int iSub,EqNumberer *, int gOffset);
   void addCContrib(GenSparseMatrix<Scalar> *coarseMat, int iSub,
                    EqNumberer* eqNumber, int cOffset, int gOffset);
   void addCGContrib(GenSparseMatrix<Scalar> *coarseMat,
                     EqNumberer* eqNumber, int cOffset, int gOffset);
   void addFContrib(GenSymFullMatrix<Scalar> *GtFG, int iSub, EqNumberer *);
   void addCtFCContrib(GenFullM<Scalar> *CtFC, int iSub, EqNumberer *);
   void addGtFCContrib(GenFullM<Scalar> *GtFC,int iSub, EqNumberer *, EqNumberer *);
   void addCtGContrib(GenFullM<Scalar> *CtG, EqNumberer *,EqNumberer *);
   void addMPCContrib(int iMPC, GenSparseMatrix<Scalar> *coarseMat,
                      EqNumberer* eqNumber, int cOffset, int gOffset, int mOffset,
                      double *locR, GenSolver<Scalar> *s);

   void getGtMult(const Scalar *vec, Scalar *alpha) const;
   void getGtQMult(const Scalar *vec, Scalar *alpha) const;
   void getGtQMult(int iSub, const Scalar *vec, Scalar *alpha) const;
   void getGtFMult(const Scalar *vec, Scalar *alpha) const;
   void getGtFMult(int iSub, const Scalar *vec, Scalar *alpha) const;
   void addLocAlphaG(Scalar *vec, Scalar *alpha);
   void subLocAlphaG(Scalar *vec, Scalar *alpha);
   void subLocAlphaQG(Scalar *vec, Scalar *alpha);
   void subLocAlphaFG(Scalar *vec, Scalar *alpha);
   void subAlphaGtQ(int, const Scalar *vec, Scalar *alpha) const;
   void subAlphaNeighbQG(int, Scalar *vec, Scalar *alpha);
   void subAlphaNeighbFG(int, Scalar *vec, Scalar *alpha);
   void subNuNeighbFC(Scalar *vec, Scalar *alpha, int *);
   void subAlphaNeighbG(int, Scalar *vec, Scalar *neighbG,
                        int lda, Scalar *alpha);
   void addAlphaNeighbG(int, Scalar *vec, Scalar *neighbG,
                        int lda, Scalar *alpha);

   void addNuC(Scalar *vec, Scalar *sv, EqNumberer *eqNums, int cOffset);
   void subLocNuC(Scalar *vec, Scalar *nu);
   void subLocNuFC(Scalar *vec, Scalar *nu);
   void subNuNeighbC(Scalar *vec, Scalar *nu, int *offsets);

   void getCtMult(const Scalar *vec, Scalar *beta) const;
   void getCtFMult(const Scalar *vec, Scalar *beta) const;
   void getCtFMult(int iSub, const Scalar *vec, Scalar *beta) const;

   int offset(int);
   Scalar *subG(int iSub)  { return locGs  + offset(iSub); }
   Scalar *subQG(int iSub) { return locQGs + offset(iSub); }
   Scalar *subBC(int iSub);
   void makeBCs();

   int subSize(int i)  { return subOffset[i+1]-subOffset[i]; }
   int edgeSize(int i);

   // Nonlinear
   void reGetNeighbGs(GenCoarseSet<Scalar> *csets, GenSubDomain<Scalar> *sd);
   void reGetNeighbQGs(GenCoarseSet<Scalar> *csets, GenSubDomain<Scalar> *sd, int kc, GenSolver<Scalar> * = 0);
   void reComputeNeighbFBCs(GenSubDomain<Scalar> *sd, GenSolver<Scalar> *solver);
};

typedef GenCoarseSet<double> CoarseSet;

#endif


