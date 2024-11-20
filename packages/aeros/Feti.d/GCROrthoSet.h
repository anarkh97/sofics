#ifndef _GCRORTHOSET_H_
#define _GCRORTHOSET_H_

#include <Feti.d/OrthoSet.h>

template <class Scalar> class GCROp;

template<class Scalar>
class GenGCROrthoSet : public GenOrthoSet<Scalar>
{
    Scalar *allFPFiP;

    void (GCROp<Scalar>::*operation)();

 public:
    GenGCROrthoSet(int _len, int maxsize, FSCommunicator *fetiCom);
    ~GenGCROrthoSet();

    void orthoAdd(Scalar *, Scalar *, Scalar);
    void orthoAddTimed(DistTimer &, Scalar *, Scalar *, Scalar);
    void orthogonalize(Scalar *, Scalar *, Scalar *, Scalar *);
    void orthogonalizeTimed(DistTimer &, Scalar *, Scalar *, Scalar *, Scalar *);
    void predict(Scalar *, Scalar *);
    void reset();

    friend class GCROp<Scalar>;
};

template<class Scalar>
class GCROp : public GenOrthoOp<Scalar>
{
    GenGCROrthoSet<Scalar> *os;
    Scalar *locAllFiP;
#ifdef HB_USE_MGS
    Scalar *currentvector;
    int currentindex;
#endif

  public:
    GCROp(GenGCROrthoSet<Scalar> *, int, int);
    ~GCROp();

    void addVec() override;
    void dot() override;
    void mult() override;
    void multAdd() override;
    void run() override;
    void runFor(int) override { throw "Illegal operation called on GCROp"; }

    void Fdot();
    void Fdot_mod();
    void multFAdd();
    void reset() { this->numP = 0; }
};

#ifdef _TEMPLATE_FIX_
 #include <Feti.d/GCROrthoSet.C>
#endif

#endif
