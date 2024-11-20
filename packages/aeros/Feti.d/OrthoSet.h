#ifndef _ORTHO_SET_H_
#define _ORTHO_SET_H_

#include <Threads.d/Paral.h>
#include <Utils.d/MyComplex.h>

template<class Scalar> class GenOrthoOp;

// perhaps these classes can be used to generalize some of the ortho set functions
// but this isn't yet implemented

template<class Scalar>
class GenOrthoSet 
{
  protected:
    int len;
    int numP;
    int maxP;

    Scalar *op1;
    Scalar *op2;
    Scalar *op3;

    int numTasks;
    TaskDescr **oos;
    std::mutex lock;

    FSCommunicator *fetiCom;

 public:
    GenOrthoSet() { };
    virtual ~GenOrthoSet() { };
    virtual int numDir() { return numP; }

    friend class GenOrthoOp<Scalar>;

};

template<class Scalar>
class GenOrthoOp : public TaskDescr 
{
  protected:
    Scalar *locAllP;
    long int loclen;
    int numP;
    int idx;

  public:
    GenOrthoOp() { };
    virtual ~GenOrthoOp() { };

    virtual void addVec() = 0;
    virtual void dot() = 0;
    virtual void mult() = 0;
    virtual void multAdd() = 0;
    virtual void run() = 0;

};


#endif
