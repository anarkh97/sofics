#ifndef _CONDDESCR_H_
#define _CONDDESCR_H_

class Domain;
template <class Scalar> class GenDynamMat;
typedef GenDynamMat<double> DynamMat;

// Single Domain Condition Number Post Processor

class SDCondPostProcessor {
    Domain *domain;
  public:
    // Constructor
    SDCondPostProcessor(Domain *_domain) { domain = _domain; }

    void conditionOutput(double eigmin, double eigmax);
};


// Single Domain Condition Number Problem Description

class SingleDomainCond {
    Domain *domain;
 public:
    // Constructor
    SingleDomainCond(Domain *d) { domain = d; }

    int  solVecInfo();
    void preProcess();

    SDCondPostProcessor *getPostProcessor();

    DynamMat buildCondOps();
    void getConditionInfo(double &tolerance, int &maxit, int &ndof);
};

#endif
