#ifndef _MULTI_DOMAIN_RBM_H_
#define _MULTI_DOMAIN_RBM_H_

#include <set>

template <class Scalar> class GenDecDomain;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
template <class Scalar> class GenDistrVectorSet;
template <class Scalar> class GenSubDomain;
class Connectivity;
template <class Scalar> class GenSparseMatrix;
struct DistrInfo;

template<class Scalar>
class MultiDomainRbm
{
    GenDecDomain<Scalar> *decDomain;
    double tolgrb;
    int numGtGsing;

  public:
    MultiDomainRbm(GenDecDomain<Scalar> *decDomain, double tolgrb);
    int numRBM();
    void getRBMs(GenDistrVectorSet<Scalar>& rigidBodyModes);
    void getRBMs(GenDistrVectorSet<Scalar>& rigidBodyModes, std::set<int> &rbmFilters);
    DistrInfo &solVecInfo();

  private:
    void computeRbms();
    void setBodyRBMoffset(int iSub, int *zColOffset, Connectivity *subToBody);
    void assembleGtG(int iGroup, const int *groups, const Connectivity *groupToSub, GenSparseMatrix<Scalar> *GtGsparse);
    void getGlobalRBM(int iSub, int &iRBM, GenDistrVector<Scalar> &R);
    void singularValueDecomposition(FullM &A, FullM &U, int ncol, int nrow, int &rank, double tol, FullM *V = 0);
};

#ifdef _TEMPLATE_FIX_
  #include <Solvers.d/MultiDomainRbm.C>
#endif

#endif
