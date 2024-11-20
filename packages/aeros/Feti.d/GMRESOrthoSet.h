#ifndef _GMRESORTHOSET_H_
#define _GMRESORTHOSET_H_

#include <Feti.d/OrthoSet.h>

class DistTimer;
template <class Scalar> class GMRESOp;

template<class Scalar>
class GenGMRESOrthoSet : public GenOrthoSet<Scalar>
{
    Scalar *givensC;
    Scalar *givensS;
    Scalar *g;
    Scalar *y;
    Scalar *matrixH;
 
    void (GMRESOp<Scalar>::*operation)();

 public:
    GenGMRESOrthoSet(int _len, int maxsize, FSCommunicator *fetiCom);
    ~GenGMRESOrthoSet();

    double orthoAdd(Scalar *, Scalar *);
    double orthoAddTimed(DistTimer &, Scalar *, Scalar *);
    void generateRotation(Scalar a, Scalar b, Scalar &cs, Scalar &ss);
    void applyRotation(Scalar &a, Scalar &b, Scalar cs, Scalar ss);
    void reInit();
    void init(Scalar *v0, double beta);
    double ModorthoAddTimed(DistTimer &timer, Scalar *Fv, Scalar *v);
    void solution(Scalar *u);

    friend class GMRESOp<Scalar>;
};

template<class Scalar>
class GMRESOp : public GenOrthoOp<Scalar>
{
    GenGMRESOrthoSet<Scalar> *os;
    Scalar *currentvector;
    int currentindex;

  public:
    GMRESOp(GenGMRESOrthoSet<Scalar> *, int, int);
    ~GMRESOp();

    void addVec() override;
    void dot() override;
    void mult() override;
    void multAdd() override;
    void run() override;
    void runFor(int) override { throw "Illegal operation called on GMRESOp"; }

    void reInit();
    void addVecAndNorm();
    void Moddot();  
    void ModmultAdd();  
    void takecurrent();
    void reInitcurrent();
};

#ifdef _TEMPLATE_FIX_
 #include <Feti.d/GMRESOrthoSet.C>
#endif

#endif
