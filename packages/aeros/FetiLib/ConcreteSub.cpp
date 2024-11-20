//
// Created by Michel Lesoinne on 2/23/18.
//

#include <Math.d/CuCSparse.h>
#include <Math.d/SparseSet.h>
#include "ConcreteSub.h"
#include <FetiLib/Types.h>
#include <FetiLib/Subdomain.h>
#include <FetiLib/Details/SubImpl.h>

namespace FetiLib {

namespace {

/** \brief A pseudo element representing a whole subdomain for the purpose of corner selection..
 *
 * \details FetiLib requires that a subdomain be mechanism-less and without any fake elements.
 */
class SafeElement : public Element {
public:
	int numNodes() const override {
		return static_cast<int>(n.size());
	}

	int *nodes(int *pInt) const override {
		int *res = new int[n.size()];
		for(int i = 0; i < n.size(); ++i)
			res[i] = n[i];
		return res;
	}

	bool isSafe() const override {
		return true;
	}

	// TODO add the possibility of changing the dimensionality of the problem if necessary.
	int  dim() const override { return 3; }

	SafeElement(std::vector<local_node_index> n) : n(std::move(n)) {}

private:
	std::vector<local_node_index> n;
};

}


//const FetiInfo &ConcreteBaseSub::getFetiInfo() const {
//	return <#initializer#>;
//}
//
//void test() {
//	ConcreteSub<double> cs;
//}

const FetiInfo &ConcreteBaseSub::getFetiInfo() const {
	return fetiInfo;
}

ConcreteBaseSub::ConcreteBaseSub(global_subdomain_index globSubIndex, local_subdomain_index locSubIndex,
                                 const SubImpl &subdomain, SComm *sComm)
{
	subNumber = globSubIndex;
	localSubNumber = locSubIndex;
	glNums = subdomain.glNodes;
	scomm = sComm;
//	glNumNodes = dom.numNode(); ??? TODO Do we need this?
	glToLocalNode.initialize(subdomain.glNodes.size(),subdomain.glNodes.data());
}

void ConcreteBaseSub::computeWaveNumbers() {

}

void ConcreteBaseSub::averageMatProps() {

}

const CoordSet &ConcreteBaseSub::getNodeSet() const {
	return coordinates;
}

double ConcreteBaseSub::getShiftVal() const {
	return 0.0;
}

template <typename S>
SubImpl &getter(const Subdomain<S> &s) { return  *s.subImpl; }

template<typename Scalar>
ConcreteSub<Scalar>::ConcreteSub(global_subdomain_index globSubIndex, local_subdomain_index locSubIndex,
                                 const Subdomain<Scalar> &subdomain)
	: ConcreteBaseSub(globSubIndex, locSubIndex, getter(subdomain), &getter(subdomain).scomm)
{


}

}