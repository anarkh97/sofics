//
// Created by Michel Lesoinne on 2018-11-01.
//

#include <complex>
#include "Threads.d/PHelper.h"
#include "SetOfSubs.h"
#include "Math.d/CuCSparse.h"
#include "Math.d/SparseSet.h"

template<typename Scalar>
void SetOfSubs<Scalar>::getSharedDOFs()
{
	auto myCPU = communicator->cpuNum();
	auto numSub = communicator->size();
	FSCommPattern<int> nodeIntPat(communicator, cpuToSub.get(), myCPU, FSCommPattern<int>::CopyOnSend);
	for (auto &sub: subDomain)
		sub->setNodeCommSize(&nodeIntPat);
	nodeIntPat.finalize();

	paralApply(subDomain, &FetiBaseSub::sendDOFList, &nodeIntPat);
	nodeIntPat.exchange();
	paralApply(subDomain, &FetiSub<Scalar>::gatherDOFList, &nodeIntPat);
	paralApply(subDomain, &FetiBaseSub::gatherDOFListPlus, &nodeIntPat);

	paralApply(subDomain, &FetiBaseSub::mergeInterfaces);
}

template <typename Ftor>
void
makeBasicDistrInfo(DistrInfo &info, int numSub, Ftor countFunc )
{
	info.domLen = new int[numSub];
	info.numDom = numSub;
	int totLen = 0;
	for(int iSub = 0; iSub < numSub; ++iSub) {
		info.domLen[iSub] = countFunc(iSub);
		totLen += info.domLen[iSub];
	}
	info.len = totLen;
}

template<typename Scalar>
void SetOfSubs<Scalar>::makeInternalInfo() {
	// Create internal Distributed information, only for unconstrained dofs
	if(!internalInfo) {
		internalInfo = std::make_unique<DistrInfo>();
		makeBasicDistrInfo(*internalInfo, subDomain.size(),
		                   [this](int iSub) {
			                   auto dsa = subDomain[iSub]->getDsa();
			                   auto c_dsa = subDomain[iSub]->get_c_dsa();
			                   return c_dsa ? c_dsa->size() : dsa ? dsa->size() : 0;
		                   });

//		if(domain->solInfo().inpc || domain->solInfo().timeIntegration == SolverInfo::Qstatic) {
//			setNonTrivialMasterFlag(*internalInfo);
//		} else {
			internalInfo->setMasterFlag();
//		}
	}

	// Create internal Distributed information for all dofs, both constrained and unconstrained
	if(!internalInfo2) {
		internalInfo2 = std::make_unique<DistrInfo>();
		makeBasicDistrInfo(*internalInfo2, subDomain.size(),
		                   [this](int iSub) {
			                   auto dsa = subDomain[iSub]->getDsa();
			                   return dsa ? dsa->size() : -1;
		                   } );
		internalInfo2->setMasterFlag();
	}
}

template<typename Scalar>
SetOfSubs<Scalar>::SetOfSubs(FSCommunicator *communicator,
                             std::vector<std::unique_ptr<FetiSub<Scalar> >>  subDomain,
                             const std::shared_ptr<Connectivity> &cpuToSub) :
	communicator(communicator),
	subDomain(std::move(subDomain)), cpuToSub(cpuToSub)
{
	getSharedDOFs();
}

template
class SetOfSubs<double>;
template
class SetOfSubs<std::complex<double>>;