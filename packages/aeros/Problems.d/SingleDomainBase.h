#ifndef _SINGLE_DOMAIN_BASE_H_
#define _SINGLE_DOMAIN_BASE_H_

#include <Utils.d/SolverInfo.h>
#include <complex>

class Rbm;
template <class Scalar> class GenFSFullMatrix;
typedef GenFSFullMatrix<double> FSFullMatrix;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
typedef GenSparseMatrix<std::complex<double> > ComplexSparseMatrix;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
typedef GenVector<std::complex<double> > ComplexVector;

class SingleDomainBase
{
 protected:
    SolverInfo &sinfo;
    FSFullMatrix *X;     // pre-calculated projector
    double *Rmem;        // global rigid body modes (numdof X numR)
    int numR;            // number of rigid body modes

 public:
    SingleDomainBase(SolverInfo &);
    virtual ~SingleDomainBase();

    void projector_prep(Rbm *rbms, SparseMatrix *M);
    void projector_prep(Rbm *rbms, ComplexSparseMatrix *M);
    void eigmode_projector_prep();

    void trProject(Vector &f);
    void trProject(ComplexVector &f);
    void project(Vector &v);
    void project(ComplexVector &v);
};

#endif
