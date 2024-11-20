#ifndef _ASSEMBLER_H_
#define _ASSEMBLER_H_

#include <Threads.d/Paral.h>

template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
class SComm;
template <class Scalar> class GenSubDomain;
template <class Scalar> class FSCommPattern;

template<class Scalar>
class GenAssembler 
{
  public:
    virtual void assemble(GenDistrVector<Scalar> &v, int flag = 0) = 0;
    virtual void split(GenDistrVector<Scalar> &v) = 0;
    virtual ~GenAssembler() { /* Nothing to do */ };
};

template<class Scalar>
class GenBasicAssembler : public GenAssembler<Scalar> 
{
     int numSub;
     GenSubDomain<Scalar> **sd;
     FSCommPattern<Scalar> *pat;
  public:
     GenBasicAssembler(int _numSub, GenSubDomain<Scalar> **_sd, FSCommPattern<Scalar> *_pat);
     virtual ~GenBasicAssembler() { delete pat; }
     void assemble(GenDistrVector<Scalar> &v, int flag = 0);
     void split(GenDistrVector<Scalar> &v);
  private:
     void sendSubVec(int iSub, GenDistrVector<Scalar> &v);
     void assembleSubVec(int iSub, GenDistrVector<Scalar> &v, int flag);
     void splitSubVec(int iSub, GenDistrVector<Scalar> &v);
};

typedef GenAssembler<double> Assembler;
typedef GenBasicAssembler<double> BasicAssembler;

#ifdef _TEMPLATE_FIX_
  #include <Paral.d/Assembler.C>
#endif

#endif
