#ifndef _CGORTHOSET_H_
#define _CGORTHOSET_H_

#include <Feti.d/OrthoSet.h>

template <class Scalar> class CGOrthoOp;

template<class Scalar>
class GenCGOrthoSet : public GenOrthoSet<Scalar>
{
    Scalar *allPFiP;

    int numOrthoSet;
    int *kindex;
    int op4;

    int kSize;
    int lastindex;

    void (CGOrthoOp<Scalar>::*operation)();

 public:
    GenCGOrthoSet(int _len, int maxSize, FSCommunicator *fetiCom);
    ~GenCGOrthoSet();

    void orthoAdd(Scalar *, Scalar *, Scalar);
    void orthoAddTimed(DistTimer &, Scalar *, Scalar *, Scalar);
    void orthogonalize(Scalar *, Scalar *);
    void orthogonalizeTimed(DistTimer &, Scalar *, Scalar *, bool hermitian = false);
    void precondition(Scalar *, Scalar *);
    void predict(Scalar *, Scalar *);
    int numDir()    { return kindex[numOrthoSet]-kindex[numOrthoSet-1]; }
    int lastIndex() { return kindex[numOrthoSet-1]; }
    int numOrthoSets() { return numOrthoSet; }
    void reset();
    void newOrthoSet();
    void shiftOrthoSets();
    void clean_up();

    friend class CGOrthoOp<Scalar>;
};

template<class Scalar>
class CGOrthoOp : public GenOrthoOp<Scalar>
{
    GenCGOrthoSet<Scalar> *os;
    Scalar *locAllD, *locAllFiP;

  public:
    CGOrthoOp(GenCGOrthoSet<Scalar> *, int, int);
    ~CGOrthoOp();

    void addVec() override;
    void dot() override;
    void mult() override;
    void multAdd() override;
    void run() override;
    void runFor(int) override { throw "Illegal operation called on CGOrthoOp"; }

    void addVecFric();
    void FdotSingleDir();
    void multAddSingleD();
    void multAddSingleP();
    void RdotSingleD();
    void multSubtractSingleP();
    void multSubtractSingleFp();
    void wtr();
    void wtKpr();
    void shift();
    void addWy();
    void Fdot();
    void FdotH();
    void multSubFp();
    void reset();
    void clean_up();

    long int length() { return this->loclen; }
    int offset() { return this->idx;  }
    Scalar* getAllP()   { return this->locAllP; }
    Scalar* getAllFiP() { return locAllFiP; }
    Scalar* getAllD()   { return locAllD; }
};

#endif
