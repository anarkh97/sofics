//
// Created by Michel Lesoinne on 6/25/18.
//
#include <sstream>
#include <algorithm>
#include <memory>
#include <tuple>

#include <gsl/span>

// Let's get rid of these....
#include <Math.d/CuCSparse.h>
#include <Math.d/SparseSet.h>

#include <Paral.d/MDDynam.h>
#include <Element.d/MatrixElement.d/MatrixElement.h>
#include "Feti.h"
#include <Driver.d/DecDomain.h>
#include <FetiLib/Subdomain.h>
#include <Types.h>
#include "SetOfSubs.h"
#include "CornerSelector.h"

extern double t5;

using std::vector;
using std::unique_ptr;

namespace MultiLevel {

template <typename Scalar>
class MLSub : virtual public FetiSub<Scalar> {
public:
	MLSub(gl_sub_idx glSub,
	      FetiInfo fetiInfo,
	      CoordSet coordinates,
	      Connectivity nodeToNode,
	      DofSetArray dsa,
	      GenDBSparseMatrix<Scalar> sysMat)
		: fetiInfo(std::move(fetiInfo)),
		  coordinates(std::move(coordinates)),
		  dsa(std::move(dsa)),
		  cdsa(this->dsa, 0, (const BCond *)nullptr),
		  sysMat(std::move(sysMat)),
		  nodeToNode(nodeToNode)
		  {
		this->subNumber = glSub;
	}

	const FetiInfo &getFetiInfo() const override { return fetiInfo; }

	int getNumUncon() const override { return dsa.size(); }
	const CoordSet& getNodeSet() const override { return coordinates; }
	double getShiftVal() const override {
		std::cout << "Warning, calling getShiftVal" << std::endl;
		return 0.0;
	}

	const Connectivity *getNodeToNode() const override { return &nodeToNode; }

	const DofSetArray * getDsa() const override { return &dsa; }; //!< Why do we need this?
	const ConstrainedDSA * get_c_dsa() const override { return &cdsa; /*return &dsa;*/ };
	void computeWaveNumbers() override {};
	void averageMatProps() override {};

private:
	FetiInfo fetiInfo;
	CoordSet coordinates;
	DofSetArray dsa;
	ConstrainedDSA cdsa;
	GenDBSparseMatrix<Scalar> sysMat;
	Connectivity nodeToNode;
};

}

namespace { // Anonymous workspace

/** \brief Representation of the basic data of a subdomain as a super-element.
 *
 * @tparam Scalar The scalar type for the element matrix (double or complex).
 */
template<typename Scalar>
struct SuperElement {
	// We start with global numbers but it would be good to have
	// local numbers fairly soon.
	std::vector<gl_num_t> nodes;
	std::vector<DofSet> dofs;
	GenAssembledFullM<Scalar> *Kel;
};

} // End anonymous workspace.

/** \brief Helper class for building connectivities  of SuperElements to local nodes (local to the MPI process).
 *
 * @tparam S The scalar type for the element matrix (double or complex).
 */
template <typename S>
class SetAccess<SuperElement<S>> {
public:
	SetAccess(gsl::span<const gl_sub_idx> elemIdx,
	          const lc_sub_idx *glSubToLoc,
	          const vector<SuperElement<S>> &elems,
	          const map<gl_node_idx, lc_node_idx> &glToLoc) :
		elemIdx(elemIdx),
		glSubToLoc(glSubToLoc),
		elems(elems),
		glToLoc(glToLoc) {}

	auto size() const { return elemIdx.size(); }
	auto numNodes(lc_sub_idx i) const { return ele(i).nodes.size(); }
	void nodes(lc_sub_idx i, int *nd) const {
		for(int j = 0; j < numNodes(i); ++j) {
			const auto &el = ele(i);
			auto n = el.nodes[j];
			nd[j] = glToLoc.at(n);
		}
	}
private:
	gsl::span<const gl_sub_idx> elemIdx;
	const lc_sub_idx *glSubToLoc;
	const std::vector<SuperElement<S>> &elems;
	const std::map<gl_node_idx, lc_node_idx> &glToLoc;

	const SuperElement<S> &ele(lc_sub_idx i) const { return elems[glSubToLoc[elemIdx[i]]]; }
};

namespace {

template <typename S>
std::tuple<Connectivity, DofSetArray, GenDBSparseMatrix<S> >
formMatrix(const std::vector<SuperElement<S>> &allElements, const std::map<gl_node_idx, lc_node_idx> &glToLocNodes,
           gsl::span<const gl_sub_idx> elems, const int *glSubToLoc)
{
	DofSetArray dsa(glToLocNodes.size());

	for (auto elNum : elems) {
		lc_sub_idx subIdx = glSubToLoc[elNum];
		if(subIdx < 0 || subIdx >= allElements.size())
			std::cout << "Exceeded the size: " <<  subIdx << std::endl;
		const auto &element = allElements[subIdx];
		for (int i = 0; i < element.nodes.size(); ++i) {
			auto glN = element.nodes[i];
			if(glToLocNodes.find(glN) == glToLocNodes.end()) {
				std::cerr << "Could not find " << glN << " had";
				for(auto nd: glToLocNodes)
					std::cerr << " ( " << nd.first << " : " << nd.second << ")";
				std::cerr << std::endl;
				std::cerr << "I was working on global " << elNum << " which is local " << subIdx << std::endl;
				for(auto gN : element.nodes)
					std:: cerr << " " << gN;
				std::cerr << std::endl;
				throw std::logic_error("Could not find node in formMatrix");
			}
			dsa.mark(glToLocNodes.at(element.nodes[i]), element.dofs[i]);
		}
	}
	dsa.finish();
	Connectivity eToN( SetAccess<SuperElement<S>> { elems, glSubToLoc, allElements, glToLocNodes } );

	auto nToE = eToN.reverse();
	auto nToN = nToE.transcon(eToN);

	GenDBSparseMatrix<S> K(&nToN, &dsa);
	// Assemble the matrix
	for(auto elNum : elems) {
		lc_sub_idx subIdx = glSubToLoc[elNum];
		auto &subDomElem = allElements[subIdx];
		int nDofs = 0;
		for(auto &dofs : subDomElem.dofs)
			nDofs += dofs.count();
		std::vector<int> elemDOFs(nDofs,-1);
		int *ptr = elemDOFs.data();
		for(int i = 0; i < subDomElem.nodes.size(); ++i)
			ptr += dsa.number(glToLocNodes.at(subDomElem.nodes[i]), subDomElem.dofs[i], ptr);
		K.add(*subDomElem.Kel, {elemDOFs.data(), static_cast<gsl::span<const int>::index_type> (elemDOFs.size())});
	}
#if defined(__GNUC__) && (__GNUC__ <= 5)
        return std::tuple<Connectivity, DofSetArray, GenDBSparseMatrix<S>>(std::move(nToN), std::move(dsa), std::move(K));
#else
	return { std::move(nToN), std::move(dsa), std::move(K) } ;
#endif
}

/// \brief Make a vector out of a node.
Eigen::Map<const Eigen::Vector3d> v(const Node &node) {
	return Eigen::Map<const Eigen::Vector3d>(&node.x, 3, 1);
}

//struct SubRenumAndNeighb {
//	std::vector<gl_sub_idx> neighbors;
//	std::map<gl_node_idx, lc_node_idx> glToLocNode;
//	Connectivity sharedNodes;
//};

using SubRenumAndNeighb = std::tuple<
	std::vector<gl_sub_idx>,
	std::map<gl_node_idx, lc_node_idx>,
	Connectivity
>;

/** \brief Form the list of neighbor subdomains and the set of nodes shared with each one.
 *
 * @param subIdx The subdomain for which the result is sought.
 * @param nodes The global node numbers of this subdomain.
 * @param nodeToSub The connectivity from global node to global subdomain.
 * @return A triplet with the array of global neighboring subdomains a global to local map of the nodes in the subdomain,
 * and a connectivity with the shared nodes.
 */
SubRenumAndNeighb
subSharedNodes(gl_sub_idx subIdx, gsl::span<const gl_node_idx> nodes, const Connectivity &nodeToSub)
{
	// Build the list of shared nodes with neighbors with global numbering
	std::vector<std::pair<gl_sub_idx, gl_node_idx>> sharedNodes;
	std::map<gl_node_idx, lc_node_idx> glToLocNode;
	for(auto node : nodes) {
		if (nodeToSub.num(node) > 1) // If the node is in at least one other subdomain.
			for (auto sub : nodeToSub[node])
				if (sub != subIdx) // Only account for other subs.
					sharedNodes.emplace_back(sub, node);
		glToLocNode.insert({ node, static_cast<lc_node_idx>(glToLocNode.size()) });
	}
	// Sort the list to have the nodes per neighbor and in the same order the neighbor will
	// compute them.
	std::sort(sharedNodes.begin(), sharedNodes.end());

	// Map global node numbers to local numbers.
	std::vector<std::pair<lc_sub_idx, lc_node_idx>> locShared;
	locShared.reserve(sharedNodes.size());

	std::vector<gl_sub_idx> globalNeighbors;
	for(const auto &subAndNode : sharedNodes) {
		// Check if this is a new global neighboring subdomain number.
		if(globalNeighbors.size() == 0 || subAndNode.first != globalNeighbors.back())
			globalNeighbors.push_back(subAndNode.first);
		locShared.emplace_back(globalNeighbors.size()-1, glToLocNode.at(subAndNode.second));
	}

	// Sanity check. A subdomain should not have itself amongst its neighbors. TODO REMOVE
	if(std::find(globalNeighbors.begin(), globalNeighbors.end(), subIdx) != globalNeighbors.end())
		throw "Bad apple";
#if defined(__GNUC__) && (__GNUC__ <= 5)
	return std::tuple<std::vector<gl_sub_idx>, std::map<gl_node_idx, lc_node_idx>, Connectivity>(std::move(globalNeighbors),
		 std::move(glToLocNode), Connectivity::fromLinkRange(locShared));
#else
	return { std::move(globalNeighbors), std::move(glToLocNode), Connectivity::fromLinkRange(locShared) };
#endif
};

/** \brief Build the meta subdomains.
 *
 * @tparam Scalar
 * @param cpuMetasubIndices Array of the meta-subdomains in this CPU.
 * @param metaDec Decomposition of the original subdomains into meta-subdomains.
 * @param superElems Super elements representing the original subdomains.
 * @param subToNode Connectivity from original subdomain to coarse nodes, including on neighboring CPUs.
 * @return
 */
template <typename Scalar>
vector<unique_ptr<FetiSub<Scalar>>>
formSubdomains(const gsl::span<gl_sub_idx> &cpuMetasubIndices, const Connectivity &metaDec,
               std::vector<SuperElement<Scalar>> &superElems, const Connectivity &subToNode,
               const std::map<gl_num_t, Eigen::Vector3d> &nodes, const int *glSubToLoc) {
	// Create the matrices for the super-subdomains with nodal information and
	// nodal DOFs

	vector<unique_ptr<FetiSub<Scalar>>> superSubs;
	superSubs.reserve(cpuMetasubIndices.size());

	Connectivity metaSubToNode = metaDec.transcon(subToNode);
	auto nodeToMetaSub = metaSubToNode.reverse();

	// Build the meta subs assigned to this subdomain.
	for(gl_sub_idx glSub : cpuMetasubIndices) {
		std::vector<gl_sub_idx> neighbors;
		std::map<gl_node_idx , lc_node_idx > glToLocNodes;
		Connectivity subToLocNodes;

		std::tie(neighbors, glToLocNodes, subToLocNodes) =
			subSharedNodes(glSub, metaSubToNode[glSub], nodeToMetaSub);
		int nn = static_cast<int>(neighbors.size());
		auto sComm = std::make_unique<SComm>(nn, std::move(neighbors), std::vector<lc_sub_idx>{},
		                         std::make_unique<Connectivity>(std::move(subToLocNodes)));
		// Form the meta subdomain's matrix
		auto dsaAndK = formMatrix(superElems, glToLocNodes, metaDec[glSub], glSubToLoc);

		// Form the vector of nodes.
		CoordSet coordSet;
		Eigen::Matrix<double,3, Eigen::Dynamic> subNodes(3, glToLocNodes.size());
		for(auto globLoc : glToLocNodes) {
			if(nodes.find(globLoc.first) == nodes.end()) {
				std::cout << "Could not find " << globLoc.first << std::endl;
				throw std::logic_error("Could not find node");
			}

			subNodes.col(globLoc.second) = nodes.at(globLoc.first);
			coordSet.nodeadd(globLoc.second, nodes.at(globLoc.first).data());
		}
		FetiInfo fetiInfo; // TODO get something set here.
		superSubs.push_back(
			std::make_unique<MultiLevel::MLSub<Scalar>>(glSub,
			                                            fetiInfo,
			                                            std::move(coordSet),
			                                            std::move(std::get<0>(dsaAndK)),
			                                            std::move(std::get<1>(dsaAndK)),
			                                            std::move(std::get<2>(dsaAndK))
			)
		);
		superSubs.back()->setSComm(sComm.release());
	}

	return superSubs;
}

} // End of anonymous workspace.


/*
 * Solver building does:
 buildSharedNodeComm(nodeToSub, subToNode);

 makeCorners();// Corners for FETI-DP

 getSharedDOFs();

 preProcessMPCs();//Multi-Point Constraint   N/A right now

 getSharedFSIs();  // N/A right now

 getSharedMPCs(); // N/A right now

  paralApply(subDomain, &BaseSub::mergeInterfaces);
  // Splitting needs to be done on the force, not on applied BC....!!!!
 paralApply(subDomain, &GenSubDomain<Scalar>::applySplitting);

 //paralApply(subDomain, &GenSubDomain<Scalar>::initSrc);
 makeInternalInfo();

 makeNodeInfo();

 */

/** \brief Get a vector with the 'nodes' of the subdomain created super-element
 *
 * @param sub The subdomain for which the nodes are sought.
 * @param corners Corner nodes of the subdomain, not including edge generalized nodes..
 * @param augOffset
 * @param subdomainEdgeIndices
 * @param withEdgeAugmentation
 * @return
 */
std::vector<gl_num_t> subSuperNodes(const FetiBaseSub &sub,
                                    gsl::span<const gl_num_t> corners,
                                    gl_num_t augOffset,
                                    gsl::span<const gl_num_t> subdomainEdgeIndices,
                                    bool withEdgeAugmentation) {
	std::vector<gl_num_t> nodes {corners.begin(), corners.end()};

	if(withEdgeAugmentation) {
		int iEdgeN = 0;
		for (int iNeighb = 0; iNeighb < sub.numNeighbors(); ++iNeighb) {
			if (sub.isEdgeNeighbor(iNeighb)) {
				if (sub.edgeDofs[iNeighb].count() != 0)
					nodes.push_back(augOffset + subdomainEdgeIndices[iEdgeN]);
				iEdgeN++;
			}
		}
		// MLX
		if(subdomainEdgeIndices.size() != iEdgeN)
			throw "Inconsistent number of eddges.";
	}
	return nodes;
}

template <typename T>
void
selectCorners(std::vector<std::unique_ptr<FetiSub<T>>> &subs, FSCommunicator *communicator, const Connectivity *cpuToSub,
              gl_sub_idx glNumSub) {
	vector<unique_ptr<FetiSubCornerHandler>> cornerHandlers;
	cornerHandlers.reserve(subs.size());
	vector<FetiSubCornerHandler *> fsh;
	fsh.reserve(subs.size());
	for(auto &sub: subs) {
		cornerHandlers.emplace_back(
			new FetiSubCornerHandler(sub.get())
			);
		fsh.emplace_back(cornerHandlers.back().get());
	}
	std::cout << "I am " << communicator->cpuNum() << "My cpuToSub is " << cpuToSub->csize() << " to " << cpuToSub->getNumTarget() << std::endl;
	FSCommPattern<int> cpat(communicator, cpuToSub, communicator->cpuNum(), FSCommPattern<int>::CopyOnSend);
	for(int i=0; i<subs.size(); ++i)
		subs[i]->setNodeCommSize(&cpat);
	cpat.finalize();
	CornerSelector cornerSelector(glNumSub, subs.size(), std::move(fsh), &cpat, communicator);
	int res = cornerSelector.makeCorners();
	std::cout << "Result of corner selection is " << res << std::endl;
	for(int i = 0; i < subs.size(); ++i)
		subs[i]->setCorners( cornerHandlers[i]->getCorners() );
}

/** \brief Build an upper level FetiDP Solver.
 *
 * \details The danger here is that subToCorner has to be consistent with the rest of data. Not a naturally
 * occuring fact.
 *
 * @tparam Scalar
 * @param subToCorner (global-)subdomain to global corner nodes. Not including edge-representing nodes.
 */
template <typename Scalar>
void GenFetiDPSolver<Scalar>::makeMultiLevelDPNew(const Connectivity &subToCorner) {
	int numSub = subToCorner.csize();

	// Get the decomposition of subdomains into super-subdomains.
	const Connectivity *decCoarse = this->cpuToSub;
	std::stringstream fn;
	fn << "decomposition." << numSub;
	FILE *f;
	auto filename = fn.str();
	std::cout << "Filename is " << filename << std::endl;
	if ((f = fopen(filename.c_str(),"r")) != NULL) {
		if(verboseFlag)
			filePrint(stderr, " ... Reading Decomposition from file %s !!! ...\n", filename.c_str());
		decCoarse = new Connectivity(f,numSub);
		fclose(f);
	}

	Connectivity subToCoarse = decCoarse->reverse();
	Connectivity cpuToCoarse = this->cpuToSub->transcon(subToCoarse);

	// Build the supersubdomains.
//	std::vector<FetiLib::Subdomain<Scalar>> coarseSubdomains;

	std::vector<SuperElement<Scalar>> superElements;
	superElements.reserve(this->subdomains.size());
    // Create a super-element for each subdomain and put it in the coarseSubdomains arrray.
	for(int iSub = 0; iSub < this->nsub; ++iSub) {
		// Get the node numbers.
		auto &sub = *(this->subdomains[iSub]);
		auto globalSubIndex = sub.subNum();
		auto edges = (*(this->subToEdge))[globalSubIndex];

		std::vector<gl_num_t> coarsenodes = subSuperNodes(sub, subToCorner[globalSubIndex], augOffset, edges,
			fetiInfo->augmentimpl == FetiInfo::Primal);
		// Get the DOFs
		std::vector<DofSet> dofs;
		dofs.reserve(coarsenodes.size());
		for(auto &ds : this->subdomains[iSub]->cornerDofs)
			dofs.push_back(ds);
		if (fetiInfo->augmentimpl == FetiInfo::Primal) {
			for(int iNeighb = 0; iNeighb < this->subdomains[iSub]->numNeighbors(); ++iNeighb)
				if (this->subdomains[iSub]->isEdgeNeighbor(iNeighb) &&
				    this->subdomains[iSub]->edgeDofs[iNeighb].count())
					dofs.push_back( this->subdomains[iSub]->edgeDofs[iNeighb] );
		}
		// Get the matrix
		auto Kel = this->subdomains[iSub]->Kcc.get();
		if(dofs.size() != coarsenodes.size())
			std::cerr << "BAD SIZE!\n";
		// Build the super-element
		superElements.push_back({std::move(coarsenodes), std::move(dofs), Kel});
	}

	// Set of connections from subdomain to all coarse nodes. Note however that we only know about
	// subdomains that are in this CPU or touch a subdomain on this CPU.
	std::vector<std::pair<gl_sub_idx, gl_node_idx>> glSubToAllCoarse;
	// Now form the set of nodes.
	std::map<gl_num_t, Eigen::Vector3d> nodes;

	for(int i = 0; i < this->nsub; i++) {
		Eigen::Vector3d xyz;
		gl_sub_idx glSub = this->subdomains[i]->subNum();
		int numCorner = this->subdomains[i]->numCorners();
		const auto &localCornerNodes = this->subdomains[i]->getLocalCornerNodes();
		for(int iCorner = 0; iCorner < numCorner; ++iCorner) {
			auto lcn = localCornerNodes[iCorner];
			const Node *node = this->subdomains[i]->getNodeSet()[localCornerNodes[iCorner]];
			int cornerNum = subToCorner[glSub][iCorner];
			nodes.insert({cornerNum, v(*node)});
		}
		if (fetiInfo->augmentimpl == FetiInfo::Primal) {
			const CoordSet &subnodes = this->subdomains[i]->getNodeSet();
			Connectivity &sharedNodes = *(this->subdomains[i]->getSComm()->sharedNodes);
			int iEdgeN = 0;
			for(int iNeighb = 0; iNeighb < this->subdomains[i]->numNeighbors(); ++iNeighb) {
				if(this->subdomains[i]->isEdgeNeighbor(iNeighb)) {
					if (this->subdomains[i]->edgeDofs[iNeighb].count() != 0) {
						int edgeNum = augOffset + (*(this->subToEdge))[glSub][iEdgeN];
						if (nodes.find(edgeNum) == nodes.end()) {
							xyz.setZero();
							for (auto node : sharedNodes[iNeighb])
								xyz += v(*subnodes[node]);
							xyz /= sharedNodes.num(iNeighb);
							nodes.insert({edgeNum, xyz});
						}
						glSubToAllCoarse.emplace_back(glSub, edgeNum);
					}
					++iEdgeN;
				}
			}
		}
	}

	for(gl_sub_idx iSub = 0; iSub < subToCorner.csize(); ++iSub) {
		auto nlist = subToCorner[iSub];
		for (gl_node_idx cn : nlist)
			glSubToAllCoarse.emplace_back(iSub, cn);
	}

	Connectivity subToAllCoarseNodes = Connectivity::fromLinkRange(glSubToAllCoarse);

	std::cout << "cpuToCoarse has " << cpuToCoarse.csize() << " sources and "
	<< cpuToCoarse.numConnect() << " targets." << std::endl;
	auto subs = formSubdomains(cpuToCoarse[this->myCPU], *decCoarse,
	                           superElements, subToAllCoarseNodes, nodes, this->glSubToLoc);

	selectCorners(subs, this->fetiCom, &cpuToCoarse, decCoarse->csize());

	SetOfSubs<Scalar> setOfSubs{ this->fetiCom, std::move(subs), std::make_shared<Connectivity>(cpuToCoarse)};

//	std::vector<FetiBaseSub *> baseSubs(decCoarseDomain->getAllSubDomains(),
//	                                    decCoarseDomain->getAllSubDomains()+decCoarseDomain->getNumSub());
//	paralApply(this->nsub, this->subdomains.data(), &FetiSub<Scalar>::makeKccDofsExp2,
//	           decCoarseDomain->getNumSub(),
//	           baseSubs.data(),
//	           augOffset, this->subToEdge);
}

template <typename Scalar>
void GenFetiDPSolver<Scalar>::makeMultiLevelDP(unique_ptr<const Connectivity> subToCorner) {
	makeMultiLevelDPNew(*subToCorner);

	Domain *coarseDomain = new Domain();
	coarseDomain->solInfo().solvercntl = fetiInfo->coarse_cntl;
	int numSub = subToCorner->csize();
	coarseDomain->setNumElements(subToCorner->csize());

	// Get the decomposition of subdomains into super-subdomains.
	std::unique_ptr<Connectivity> decCoarse{std::make_unique<Connectivity>(*this->cpuToSub)};
	char fn[65];
	FILE *f;
	sprintf(fn,"decomposition.%d",numSub);
	if ((f = fopen(fn,"r")) != NULL) {
		if(verboseFlag) filePrint(stderr, " ... Reading Decomposition from file %s ...\n", fn);
		decCoarse = std::make_unique<Connectivity>(f,numSub);
		fclose(f);
	}

	Connectivity elemToSubCoarse = decCoarse->reverse();
	auto CPUMapCoarse = std::make_unique<Connectivity>(this->cpuToSub->transcon(elemToSubCoarse));

	std::vector<size_t> pointer(this->glNumSub+1, 0);

	// Create a super-element for each subdomain and put it in the coarseDomain.

	// Phase one only works to get node numbers.
	for(int i = 0; i < this->nsub; ++i) {
		int globalSubIndex = this->subdomains[i]->subNum();
		auto &sub = *(this->subdomains[i]);
		auto edges = (*(this->subToEdge))[globalSubIndex];
		std::vector<gl_num_t> coarsenodes = subSuperNodes(sub, (*subToCorner)[globalSubIndex], augOffset, edges,
		                                                  fetiInfo->augmentimpl == FetiInfo::Primal);
		int n = static_cast<int> ( coarsenodes.size() );
		// TODO Fix this leak!!!!
		int *elem = new int[n];
		for(int j = 0; j < n; ++j)
			elem[j] = coarsenodes[j];
		coarseDomain->addElem(globalSubIndex,0,n,elem);//.data());
		pointer[globalSubIndex] = n;
	}

	// Set the dofs and stiffness of each super-element.
	Elemset& elems = coarseDomain->getElementSet();
	for(int i = 0; i < this->nsub; ++i) {
		int s = this->subdomains[i]->subNum();
		int n = pointer[s];
		int nc = subToCorner->num(s);
		DofSet *coarseDofs = new DofSet[n];
		for(n = 0; n < nc; n++)
			coarseDofs[n] = this->subdomains[i]->cornerDofs[n];
		if (fetiInfo->augmentimpl == FetiInfo::Primal) {
			for(int iNeighb = 0; iNeighb < this->subdomains[i]->numNeighbors(); ++iNeighb)
				if(this->subdomains[i]->isEdgeNeighbor(iNeighb) &&
				   this->subdomains[i]->edgeDofs[iNeighb].count())
					coarseDofs[n++] = this->subdomains[i]->edgeDofs[iNeighb];
		}
		((MatrixElement*)elems[s])->setDofs(coarseDofs);
		((MatrixElement*)elems[s])->setStiffness(this->subdomains[i]->Kcc.get());
	}

	// Create the nodeset.
	CoordSet& nodes = coarseDomain->getNodes();
	// Add the corner and edge nodes.
	for(int i = 0; i < this->nsub; i++) {
		Eigen::Vector3d xyz;
		int s = this->subdomains[i]->subNum();
		int numCorner = this->subdomains[i]->numCorners();
		const auto &localCornerNodes = this->subdomains[i]->getLocalCornerNodes();
		for(int iCorner = 0; iCorner < numCorner; ++iCorner) {
			auto lcn = localCornerNodes[iCorner];
//			auto gcn = this->subdomains[i]->
			const Node *node = this->subdomains[i]->getNodeSet()[localCornerNodes[iCorner]];
			int cornerNum = (*subToCorner)[s][iCorner];
			if (!nodes[cornerNum]) {
				xyz[0] = node->x; xyz[1] = node->y; xyz[2] = node->z;
				nodes.nodeadd(cornerNum, xyz.data());
			}
		}
		if (fetiInfo->augmentimpl == FetiInfo::Primal) { // 020314 JAT
			const CoordSet &subnodes = this->subdomains[i]->getNodeSet();
			Connectivity &sharedNodes = *(this->subdomains[i]->getSComm()->sharedNodes);
			int iEdgeN = 0, edgeNum;
			for(int iNeighb = 0; iNeighb < this->subdomains[i]->numNeighbors(); ++iNeighb) {
				if(this->subdomains[i]->isEdgeNeighbor(iNeighb)) {
					edgeNum = augOffset + (*(this->subToEdge))[s][iEdgeN++];
					if (!nodes[edgeNum]) {
						xyz.setZero();
						for(auto node : sharedNodes[iNeighb])
							xyz += v(*subnodes[node]);
						xyz /= sharedNodes.num(iNeighb);
						nodes.nodeadd(edgeNum, xyz.data());
					}
				}
			}
		}
	}


	coarseDomain->setNumNodes(cornerToSub->csize());

#ifdef USE_MPI
	Communicator *structCom = new Communicator(CommunicatorHandle{this->fetiCom->getComm()});
#else
	Communicator *structCom = NULL;
#endif
	GenDecDomain<Scalar> *decCoarseDomain = new GenDecDomain<Scalar>(coarseDomain, structCom, false, true);

	std::unique_ptr<const Connectivity> elemToNode; // JAT 220216
	if (fetiInfo->augmentimpl == FetiInfo::Primal) { // JAT 041114
#ifdef DISTRIBUTED
		this->fetiCom->globalSum(this->glNumSub, pointer.data());
#endif
		int total = 0;
		for(int iSub = 0; iSub < this->glNumSub; ++iSub) {
			int tmp = pointer[iSub];
			pointer[iSub] = total;
			total += tmp;
		}
		pointer[this->glNumSub] = total;

		std::vector<int> target(total, 0);
		for(int i = 0; i < this->nsub; i++) {
			int s = this->subdomains[i]->subNum();
			int nc = subToCorner->num(s);
			int n, k = pointer[s];
			for(n = 0; n < nc; n++)
				target[k+n] = (*subToCorner)[s][n];
			int iEdgeN = 0;
			for(int iNeighb = 0; iNeighb < this->subdomains[i]->numNeighbors(); ++iNeighb) {
				if(this->subdomains[i]->isEdgeNeighbor(iNeighb)) {
					if(this->subdomains[i]->edgeDofs[iNeighb].count())
						target[k+n++] = augOffset + (*(this->subToEdge))[s][iEdgeN];
					iEdgeN++;
				}
			}
		}
#ifdef DISTRIBUTED
		this->fetiCom->globalSum(total, target.data());
#endif
		elemToNode = std::make_unique<Connectivity>(this->glNumSub, std::move(pointer), std::move(target));
	} else
		elemToNode = std::move(subToCorner);

	decCoarseDomain->setElemToNode(std::move(elemToNode));
	decCoarseDomain->setSubToElem(std::move(decCoarse));
	decCoarseDomain->setCPUMap(std::move(CPUMapCoarse));
	decCoarseDomain->preProcess();
	for(int i=0; i<decCoarseDomain->getNumSub(); ++i) decCoarseDomain->getSubDomain(i)->makeAllDOFs();
	GenMDDynamMat<Scalar> ops;
	if(verboseFlag) filePrint(stderr, " ... Factor Kcc solver              ...\n");
	decCoarseDomain->buildOps(ops, 0.0, 0.0, 1.0);
	coarseInfo = &(decCoarseDomain->solVecInfo());
	KccParallelSolver = ops.dynMat;
	std::vector<FetiBaseSub *> baseSubs(decCoarseDomain->getAllSubDomains(),
	                                    decCoarseDomain->getAllSubDomains()+decCoarseDomain->getNumSub());
	paralApply(this->subdomains, &FetiSub<Scalar>::makeKccDofsMultiLevel,
	           baseSubs,
	           augOffset, this->subToEdge); // JAT 101816
}

template
void GenFetiDPSolver<double>::makeMultiLevelDP(unique_ptr<const Connectivity> subToCorner);

template
void GenFetiDPSolver<std::complex<double>>::makeMultiLevelDP(unique_ptr<const Connectivity> subToCorner);
