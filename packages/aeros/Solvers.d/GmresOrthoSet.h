#ifndef __GMRESORTHOSET_H__
#define __GMRESORTHOSET_H__

#include <Threads.d/Paral.h>

template <class Scalar> class GmresOp;
class FSCommunicator;

template<class Scalar>
class GmresOrthoSet
{
    Scalar *givensC;
    Scalar *givensS;
    Scalar *g;
    Scalar *y;
    Scalar *matrixH;

    int len;
    int numP;
    int maxP;

    Scalar *op1;
    Scalar *op2;
    Scalar *op3;

    int numTasks;
    TaskDescr **oos;
    std::mutex lock;

    FSCommunicator *com;
 
    void (GmresOp<Scalar>::*operation)();

 public:
    GmresOrthoSet(int _len, int maxsize, FSCommunicator *);
    ~GmresOrthoSet();

    int numDir() { return numP; }
    double orthoAdd(Scalar *, Scalar *);
    void generateRotation(Scalar a, Scalar b, Scalar &cs, Scalar &ss);
    void applyRotation(Scalar &a, Scalar &b, Scalar cs, Scalar ss);
    void reset();
    void init(Scalar *v0, double beta);
    double init(Scalar *p);
    double ModorthoAdd(Scalar *Fv, Scalar *v);
    void solution(Scalar *u);

    friend class GmresOp<Scalar>;
};

template<class Scalar>
class GmresOp : public TaskDescr
{
    GmresOrthoSet<Scalar> *os;
    Scalar *currentvector;
    int currentindex;
    Scalar *locAllP;
    int loclen;
    int numP;
    int idx;

  public:
    GmresOp(GmresOrthoSet<Scalar> *, int, int);
    ~GmresOp();

    void addVec();
    void dot();
    void mult();
    void multAdd();
    void run() override;
    void runFor(int) override { throw "Illegal operation called on GmresOp"; }

    void reset();
    void addVecAndNorm();
    void Moddot();  
    void ModmultAdd();  
    void takecurrent();
    void reInitcurrent();
};

#ifdef _TEMPLATE_FIX_
 #include <Solvers.d/GmresOrthoSet.C>
#endif

#endif
