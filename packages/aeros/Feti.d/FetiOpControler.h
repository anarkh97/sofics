#ifndef _FETIOPCONTROLER_H_
#define _FETIOPCONTROLER_H_

#include <Utils.d/dofset.h>

template <class Scalar> class GenFullM;
template <class Scalar> class GenSymFullMatrix;
template <class Scalar> class GenSparseMatrix;
template <class Scalar> class GenDistrVector;
template <class Scalar> class GenCoarseSet;

template<class Scalar>
class GenFetiOpControler {
  public:

    GenCoarseSet<Scalar> *cset;
    GenDistrVector<Scalar> *dv1, *dv2, *dv3, *dv4;
    Scalar *vec1, *vec2, *vec3, *vec4;
    Scalar r;
    int nonLocalQ; // Whether or not Q has local effect
    int nQ;

    GenFullM<Scalar> *coarseGtCMult(int, int);
    GenFullM<Scalar> *coarseGtCMult(int);

    GenSparseMatrix<Scalar> *sparseGtG;
    SimpleNumberer *eqNums;
    SimpleNumberer *PFcNums;

    GenFullM<Scalar> *GtFC;
    GenFullM<Scalar> *CtFCs;
    GenSymFullMatrix<Scalar> *GtQGs;

    Scalar *Rgc;    // Rgc is stored in columnwise!!!!
    int    crns;    // number of column of Rgc
    int    rbms;    // number of row of Rgc

    void **globalSumBuffer;

    int numNeighbSubs(int sub); // How many neighbors sub has
    int neighbSubId(int sub, int isub); // number of the isub'th neighbor to sub
    int index(int sub, int subJ);  // identify the order of subJ in sub's list

    GenFetiOpControler(int nsub, int loc) {
        cset = new GenCoarseSet<Scalar>[nsub];
        nonLocalQ =!loc;

        globalSumBuffer = new void *[nsub];

        CtFCs = (GenFullM<Scalar> *)0;
        GtQGs = (GenSymFullMatrix<Scalar> *) 0;
        GtFC  = (GenFullM<Scalar> *) 0;
    }

    void (GenFetiOp<Scalar>::*operation)();
};

typedef GenFetiOpControler<double> FetiOpControler;

#ifdef _TEMPLATE_FIX_
  #include <Feti.d/FetiOpControler.C>
#endif

#endif


