#ifndef _CCTSOLVER_H_
#define _CCTSOLVER_H_

#include <vector>
template <typename Scalar>
class GenDistrVector;
template <typename Scalar>
class FetiSub;

class FSCommunicator;
class Connectivity;

template<class Scalar>
class CCtSolver
{
  public:
	CCtSolver(std::vector<FetiSub<Scalar>*> subsWithMpcs) : subsWithMpcs(std::move(subsWithMpcs)) {}
    virtual void reSolve(GenDistrVector<Scalar> &v) = 0;
    virtual void zeroAll() = 0;
    virtual void assemble() = 0;
    virtual void factor() = 0;
    virtual ~CCtSolver() { };
  protected:
    FSCommunicator *fetiCom;
    const Connectivity *mpcToCpu;
    int numSubsWithMpcs;
    std::vector<FetiSub<Scalar>*> subsWithMpcs;
    int glNumMpc;
};

#endif
