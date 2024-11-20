#ifndef _MULTI_DOMAIN_BASE_H_
#define _MULTI_DOMAIN_BASE_H_

#include <Utils.d/SolverInfo.h>
#include <complex>

template <class Scalar> class GenDistrVectorSet;
typedef GenDistrVectorSet<double> DistrVectorSet;
template <class Sclaar> class MultiDomainRbm;
template <class Scalar> class GenSubDOp;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
typedef GenDistrVector<std::complex<double> > ComplexDistrVector;
struct DistrInfo;

class MultiDomainBase
{
 protected:
    SolverInfo &sinfo;
    DistrVectorSet *X, *R;
    int numR;

 public:
    MultiDomainBase(SolverInfo &);
    virtual ~MultiDomainBase();

    void projector_prep(MultiDomainRbm<double> *rbms, GenSubDOp<double> *M);
    void projector_prep(MultiDomainRbm<std::complex<double> > *rbms, GenSubDOp<std::complex<double> > *M);
    void eigmode_projector_prep(DistrInfo &);

    void trProject(DistrVector &f);
    void trProject(ComplexDistrVector &f);
    void project(DistrVector &v);
    void project(ComplexDistrVector &v);
};

#endif
