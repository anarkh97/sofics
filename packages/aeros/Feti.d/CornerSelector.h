//
// Created by Michel Lesoinne on 4/12/18.
//

#ifndef FEM_CORNERSELECTOR_H
#define FEM_CORNERSELECTOR_H


#include <vector>
#include <gsl/span>
#include "FetiSub.h"

class FetiSubCornerHandler;
class Connectivity;
class FSCommunicator;
class CoordSet;
class DofSetArray;
class ConstrainedDSA;
class BaseSub;
template <class Type> class FSCommPattern;

class CornerSelector
{
	int glNumSub;
	int nSub;
	std::unique_ptr<const Connectivity> grToSub;
	int *glSubGroup;
	std::vector<FetiSubCornerHandler*> cornerHandler;
	FSCommPattern<int> *cpat;
	FSCommunicator *communicator;
	int dim;
	int dims[4];

	void chooseCorners(char *glCornerList, double (*xyz)[3],
	                   Connectivity &cNConnect, Connectivity &subToRotCrn,
	                   int *glCrnGroup);
public:
	CornerSelector(int nGlobSub, int nLocSub,
	               std::vector<FetiSubCornerHandler *> handlers,
	               FSCommPattern<int> *commPattern, FSCommunicator *communicator);
	~CornerSelector();
	/** \brief Pick the corners for FETI-DP.
	 *
	 * @return The total number of corners picked in the complete problem.
	 */
	int makeCorners();
	std::unique_ptr<const Connectivity> yieldGrToSub() { return std::move(grToSub); }
	auto &handlers() { return cornerHandler; }
};

// This class is the algorithmic class to select corner nodes
class FetiSubCornerHandler
{

public:
	/**
	 *
	 * @param sub Global number of this subdomain.
	 * @param nn
	 * @param n Coordinates of the subdomain's nodes.
	 * @param nTn Node to node connectivity.
	 * @param d Unconstrained set of DOFs for the subdomain.
	 * @param sh Shared nodes
	 * @param nsb Neighboring subdomains
	 * @param c_dsa
	 * @param _subPre
	 */
	FetiSubCornerHandler(gl_sub_idx sub, int nn, const CoordSet &n, const Connectivity &nTn, const DofSetArray &d,
	                     const Connectivity &sh,
	                     const std::vector<gl_sub_idx> &nsb,
	                     const ConstrainedDSA *c_dsa, FetiBaseSub *_subPre);
	FetiSubCornerHandler(FetiBaseSub *sub);
	~FetiSubCornerHandler() = default;
	void markMultiDegNodes();
	void dispatchSafeNodes(FSCommPattern<int> *);
	void markSafeNodes(FSCommPattern<int> *);
	void dispatchRotCorners(FSCommPattern<int> *);
	void markRotCorners(FSCommPattern<int> *);
	void pickAnyCorners();
	void countAndMarkCornerCand(int *mync, int *totnc);
	void getCornerXYZ(int *, double (*)[3], char *essential, int *cTsP, int *cTsT);
	void dispatchNumbering(FSCommPattern<int> *pat, char *crnMrk,
	                       int *allOrigFC, int *allNewFC, int, int *cntOff);
	void dispatchInitialNumbering(FSCommPattern<int> *pat, int *firstC);
	void recNumbering(FSCommPattern<int> *, int *fM);
	void recInitialNumbering(FSCommPattern<int> *pat, int *numRotCrn);
	void listRotCorners(int *fN, int *crnNum);
	void countContact(int *, char *crnMrk);
	void markDims(int *_dims);

	void resendNumbers(FSCommPattern<int> *pat);
	void checkNumbers(FSCommPattern<int> *pat);

	gsl::span<const lc_node_idx> getCorners() const { return crnList; }
	int getNumCorners() { return totNC; }

protected:
	int glSubNum;
	int nnodes;
	std::vector<int> deg;
	std::vector<lc_node_idx> crnList;
	int totNC;
	std::vector<int> weight;
	std::vector<bool> isCorner;
	std::vector<bool> isSafe;
	std::vector<bool> glSafe;
	int nNeighb;
	const Connectivity &sharedNodes;
	const Connectivity &nToN;
	const std::vector<int> &neighbSubs;
	const CoordSet &nodes;
	const DofSetArray &dsa;
	std::vector<bool> isRotMidSideNode;
	int dim;
	int dims[4];
	bool allSafe;
	int nTC; // Total number of corner candidates for this sub
	bool checkForColinearCrossPoints(int numCornerPoints,
	                                 int *localCornerPoints);
	using bit_iterator = std::vector<bool>::const_iterator;
	bool addPotCornerPoints(gsl::span<const int> allNodes, std::vector<bool>::const_iterator isSafe);
	bool mixed; // true if subdomain has active fluid and structure dofs
	FetiBaseSub *subPre;
};

#endif //FEM_CORNERSELECTOR_H
