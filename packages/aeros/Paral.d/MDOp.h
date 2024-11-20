#ifndef _MD_OP_H_
#define _MD_OP_H_

//#include <Driver.d/SubDomain.h>
//#include <Feti.d/DistrVector.h>

template <class Scalar> class GenSubDomain;
typedef GenSubDomain<double> SubDomain;
template <class Scalar> class GenCuCSparse;
typedef GenCuCSparse<double> CuCSparse;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
class DistrGeomState;
template <class Scalar> class GenSubDOp;
typedef GenSubDOp<double> SubDOp;
class ControlInterface;

class MultiDomainOp : public TaskDescr {
    SubDomain **sd;
    DistrVector *v1, *v2, *v3, *v4;
    double c1, c2;
    ControlInterface *userSupFunc;
    double **temprcvd;
    SubDOp *Kuc, *Cuc, *Muc;

    void  (MultiDomainOp::*f)(int);
 public:
    MultiDomainOp(void (MultiDomainOp::*_f)(int),  SubDomain **,
                  DistrVector*, DistrVector*, double, SubDOp*,
                  ControlInterface*, SubDOp*, double, SubDOp*);
    MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **,
                  DistrVector*, DistrVector*, DistrVector*, DistrVector*);
    MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **,
                  DistrVector*, DistrVector*, DistrVector*);
    MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **,
                  DistrVector*, SubDOp*);
    MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **);
    MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **,
                  DistrVector*, DistrVector*);

    void computeExtForce(int);
    void getConstForce(int);
    void getInitState(int);
    void runFor(int);
    void makeAllDOFs(int);

};

#endif
