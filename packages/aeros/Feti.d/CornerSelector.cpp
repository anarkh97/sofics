//
// Created by Michel Lesoinne on 4/12/18.
//

#include "CornerSelector.h"

#include <Feti.d/FetiSub.h>
#include <cfloat>
#include <algorithm>

//#define DEBUG_CORNER

extern SolverInfo &solInfo;

FILE *myOut;
void crossprod(double [3], double [3], double [3]);

double magnitude(double[3]);

double magnitude(double v[3]) {
	return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// variable used to select corner nodes for the 2D-case
DofSet XYDofs = DofSet::Xdisp | DofSet::Ydisp;

static double checkArea(int np, gsl::span<const int> pt, double pxyz[][3], int dim =3)
{
	double area = 0.0;
	if(dim < 3) {
		for(int i = 1; i < np; ++i) {
			double dx1 = pxyz[pt[i-1]][0] - pxyz[pt[0]][0];
			double dy1 = pxyz[pt[i-1]][1] - pxyz[pt[0]][1];
			double dz1 = pxyz[pt[i-1]][2] - pxyz[pt[0]][2];
			area += fabs(sqrt(dx1*dx1+dy1*dy1+dz1*dz1));
		}
		return area;
	}
	for(int i = 2; i < np; ++i)
	{
		double dx1 = pxyz[pt[i-1]][0] - pxyz[pt[0]][0];
		double dy1 = pxyz[pt[i-1]][1] - pxyz[pt[0]][1];
		double dz1 = pxyz[pt[i-1]][2] - pxyz[pt[0]][2];
		double dx2 = pxyz[pt[i]][0] - pxyz[pt[0]][0];
		double dy2 = pxyz[pt[i]][1] - pxyz[pt[0]][1];
		double dz2 = pxyz[pt[i]][2] - pxyz[pt[0]][2];
		area += fabs(sqrt(
				(dy1*dz2-dy2*dz1)*(dy1*dz2-dy2*dz1) +
				(dz1*dx2-dz2*dx1)*(dz1*dx2-dz2*dx1) +
				(dx1*dy2-dx2*dy1)*(dx1*dy2-dx2*dy1)));
	}
	return area;
}


FetiSubCornerHandler::FetiSubCornerHandler(gl_sub_idx _glSubNum, int _nnodes, const CoordSet &_nodes,
                                           const Connectivity &_nToN,
                                           const DofSetArray &_dsa,
                                           const Connectivity &_sharedNodes, const std::vector<gl_sub_idx> &_neighbSubs,
                                           const ConstrainedDSA *c_dsa, FetiBaseSub *_subPre)
		:
		nnodes(_nnodes),
		isCorner(_nnodes, false),
		glSafe(_nnodes, true),
		sharedNodes(_sharedNodes), nToN(_nToN), neighbSubs(_neighbSubs),
		nodes(_nodes), dsa(_dsa),  isRotMidSideNode(_nnodes, false), dim(3)
{
	glSubNum = _glSubNum;
	nnodes = _nnodes;
	subPre = _subPre;
	deg.resize(nnodes);
	weight.resize(nnodes);
	for(int iNode = 0; iNode < nnodes; ++iNode) isCorner[iNode] = false;
	nNeighb = sharedNodes.csize();

	// dimension
	DofSet dofs = DofSet();
	for(int i=0; i<dsa.numNodes(); ++i) dofs.mark(dsa[i].list());
	dims[0] = (dofs.contains(DofSet::Xdisp)) ? 1 : 0;
	dims[1] = (dofs.contains(DofSet::Ydisp)) ? 1 : 0;
	dims[2] = (dofs.contains(DofSet::Zdisp)) ? 1 : 0;
	dims[3] = (dofs.contains(DofSet::Helm) || dofs.contains(DofSet::Temp)) ? 1 : 0;
	int sdim = dims[0]+dims[1]+dims[2]; // structure
	int fdim = dims[3]; // fluid
	mixed = (fdim && sdim);

	// active dimensions
	int cdims[4];
	DofSet cdofs = DofSet();
	for(int i=0; i<c_dsa->numNodes(); ++i) cdofs.mark((*c_dsa)[i].list());
	cdims[0] = (cdofs.contains(DofSet::Xdisp)) ? 1 : 0;
	cdims[1] = (cdofs.contains(DofSet::Ydisp)) ? 1 : 0;
	cdims[2] = (cdofs.contains(DofSet::Zdisp)) ? 1 : 0;
	cdims[3] = (cdofs.contains(DofSet::Helm) || cdofs.contains(DofSet::Temp)) ? 1 : 0;
	int cdim = cdims[0]+cdims[1]+cdims[2]+cdims[3];
	allSafe = ((cdim == 1) || !solInfo.getFetiInfo().pick_unsafe_corners);  // if there is only one active dimension then all elements are safe
	// can ignore unsafe corners by setting pick_unsafe_corners to false (maybe ok if using pivoting local solver like spooles or mumps)

#ifdef DEBUG_CORNER
	cerr << "glSubNum = " << glSubNum << ", dim = " << dim << ", dims = " << dims[0] << " " << dims[1] << " "
       << dims[2] << " " << dims[3] << ", mixed = " << mixed << ", allSafe = " << allSafe << ", nnodes = " << nnodes << endl;
#endif
}


FetiSubCornerHandler::FetiSubCornerHandler(FetiBaseSub *sub) :
	FetiSubCornerHandler(sub->subNum(),
	                     sub->getNodeToNode()->csize(),
	                     sub->getNodeSet(),
	                     *sub->getNodeToNode(),
	                     *sub->getDsa(),
	                     *sub->getSComm()->sharedNodes,
	                     sub->getSComm()->subNums,
	                     sub->get_c_dsa(),
	                     sub)
{

}

void
FetiSubCornerHandler::dispatchSafeNodes(FSCommPattern<int> *cpat)
{
	// Dispatch the nodes that are not safe
	for(int iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		FSSubRecInfo<int> rInfo = cpat->getSendBuffer(glSubNum, neighbSubs[iNeighb]);
		for(int iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode) {
			lc_node_idx node = sharedNodes[iNeighb][iNode];
			rInfo.data[iNode] = int(glSafe[node]);
		}
	}
}

void
FetiSubCornerHandler::markSafeNodes(FSCommPattern<int> *cpat)
{
	// if node is unsafe in at least one subdomain, mark it as unsafe in all subdomains
	for(int iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		FSSubRecInfo<int> ri = cpat->recData(neighbSubs[iNeighb], glSubNum);
		for(int iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode) {
			if(ri.data[iNode] == 0)
				glSafe[sharedNodes[iNeighb][iNode]] = false; // node is unsafe in neighboring subdomain
		}
	}

#ifdef DEBUG_CORNER
	int cnt = 0;
  for(int iNode = 0; iNode < nnodes; ++iNode) if(!glSafe[iNode]) cnt++;
  cerr << "sub " << glSubNum << " has " << cnt << " unsafe nodes and " << nnodes-cnt << " safe nodes\n";
#endif

	// fill isSafe array, indexed by sharedNodes
	isSafe.resize(sharedNodes.numConnect());
	for(int iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		auto lSafe = isSafe.begin() + sharedNodes.offset(iNeighb);
		for(int iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode)
			lSafe[iNode] = glSafe[sharedNodes[iNeighb][iNode]];
	}
}

void
FetiSubCornerHandler::markMultiDegNodes()
{
	for(int iNode = 0; iNode < nnodes; ++iNode) deg[iNode] = weight[iNode] = 0;

	// count the number of neighbors sharing each safe node (not including this subdomain) and store in deg
	// To be perfectly correct, we should also look at the sharedDofs really
	// this may be done if we reorganize the set up of the shared dofs XML
	for(int iNeighb = 0; iNeighb < nNeighb; ++iNeighb)
		for(int iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode)
			deg[sharedNodes[iNeighb][iNode]]++;

	// Now find faces that could contribute to troubles
	// set flag = true if face has one or more nodes shared by only one other neighbor (ie, deg == 1)
	// then count again number of neighbors sharing each safe node on multiDeg edges
	for(int iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		bool flag = false;
		for(int iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode)
			if(deg[sharedNodes[iNeighb][iNode]] == 1) { flag = true; break; }
		if(!flag) { // face is entirely made up of nodes with more than one neighbor
			for(int iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode)
				weight[sharedNodes[iNeighb][iNode]]++;
		}
	}

	if(!solInfo.isCoupled || solInfo.getFetiInfo().fsi_corner == 0) { // JLchange
		for(int iNode = 0; iNode < nnodes; ++iNode) {
			if(weight[iNode] >= 2) isCorner[iNode] = true;
		}
	}
	else {
		for(int iNode = 0; iNode < nnodes; ++iNode)
			if((weight[iNode] >= 2) && !subPre->onWetInterface(iNode))
				isCorner[iNode] = true;
	}

#ifdef DEBUG_CORNER
	int cnt = 0;
  for(int iNode = 0; iNode < nnodes; ++iNode) if(isCorner[iNode]) cnt++;
  cerr << "sub " << glSubNum << " has " << cnt << " multi-degree corners\n";
#endif
}

void
FetiSubCornerHandler::dispatchRotCorners(FSCommPattern<int> *cpat)
{
	for(int iNode = 0; iNode < nnodes; ++iNode) { weight[iNode] = 0; deg[iNode] = -1; }

	for(int iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		for(int iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode)
			deg[sharedNodes[iNeighb][iNode]] = iNeighb;
		for(int iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode) {
			int count = 0;
			int node = sharedNodes[iNeighb][iNode];
			if(isRotMidSideNode[node]) continue;  // PJSA: fix for 6 node tri shell
			//  include end points of all subdomain edges as corner nodes for the 2D-case XXXX not general enough
			if(dsa[iNode].contains(DofSet::XYZrot) || (dsa[iNode] == XYDofs) || (dsa[iNode] == DofSet::Helm) || (dsa[iNode] == DofSet::Temp)) {
				for(int i = 0; i < nToN.num(node); ++i)
					if((deg[nToN[node][i]] == iNeighb) && !isRotMidSideNode[nToN[node][i]]) count++; // fix for 6 node tri shell
				if(count < 3) {
					weight[node] = 2;
					if(!solInfo.isCoupled || solInfo.getFetiInfo().fsi_corner == 0) // JLchange
						isCorner[node] = true;
					else { if ( !(subPre->onWetInterface(node)) ) isCorner[node] = true; }
				}
			}
		}
	}

	for(int iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		FSSubRecInfo<int> rInfo = cpat->getSendBuffer(glSubNum, neighbSubs[iNeighb]);
		for(int iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode)
			rInfo.data[iNode] = weight[sharedNodes[iNeighb][iNode]];
	}
}

void
FetiSubCornerHandler::markRotCorners(FSCommPattern<int> *cpat)
{
	for(int iNode = 0; iNode < nnodes; ++iNode) deg[iNode] = -1;

	for(int iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		FSSubRecInfo<int> ri = cpat->recData(neighbSubs[iNeighb], glSubNum);
		for(int iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode)
			if(ri.data[iNode] > 1) {
				if(!solInfo.isCoupled || solInfo.getFetiInfo().fsi_corner == 0) // JLchange
					isCorner[sharedNodes[iNeighb][iNode]] = true;
				else {
					if(!subPre->onWetInterface(sharedNodes[iNeighb][iNode]))
						isCorner[sharedNodes[iNeighb][iNode]] = true;
				}
			}
	}

#ifdef DEBUG_CORNER
	int cnt = 0;
  for(int iNode = 0; iNode < nnodes; ++iNode) if(isCorner[iNode]) cnt++;
  cerr << "sub " << glSubNum << " has " << cnt << " rotational and multi-degree corners\n";
#endif
}

void
FetiSubCornerHandler::pickAnyCorners()
{
#ifdef DEBUG_CORNER
	cerr << "glSubNum = " << glSubNum << " in pickAnyCorner\n";
#endif
	if(nnodes < dim) {
		for(int i = 0; i < nnodes; ++i) {
			if(!solInfo.isCoupled || solInfo.getFetiInfo().fsi_corner == 0) // JLchange
				isCorner[i] = true;
			else { if ( !(subPre->onWetInterface(i)) ) isCorner[i] = true; }
			weight[i] = 3;
		}
		fprintf(stderr, "Really not enough nodes: nnodes = %d dim = %d\n", nnodes, dim);
		return;
	}
	int n1 = 0;
	while(n1 <nnodes && glSafe[ n1 ] == false )
		n1++;
	// fprintf(stderr, "N1 is %d out of %d in %d dim\n", n1, nnodes, dim);
	if(n1+dim > nnodes) return;
	double maxDist = 0.0, dist;
	Node &nd1 = *nodes[n1];
	int n2 = n1;
	int n;
	for(n = n1+1; n < nodes.size(); ++n) { // XXXX
		Node &nd2 = *nodes[n];
		dist = nd2.distance(nd1);
		if(glSafe[n] && dist > maxDist) { maxDist = dist; n2 = n; }
	}
	if(!solInfo.isCoupled || solInfo.getFetiInfo().fsi_corner == 0) { // JLchange
		isCorner[n1] = true;
		isCorner[n2] = true;
	}
	else {
		if(!subPre->onWetInterface(n1)) isCorner[n1] = true;
		if(!subPre->onWetInterface(n2)) isCorner[n2] = true;
	}
	weight[n1] = 3;
	weight[n2] = 3;
	if((dist > 0 && dim <= 2) || dist == 0) return;

	Node &nd2 = *nodes[n2];
	double dx= nd2.x - nd1.x;
	double dy= nd2.y - nd1.y;
	double dz= nd2.z - nd1.z;

	// finally, find the node that maximizes the area with nodes 1 & 2
	double maxArea=0.0, area;
	int n3 = n1;
	for(n=0; n<nodes.size(); ++n) { // XXXX
		Node &nd3 = *nodes[n];
		double dx2= nd3.x - nd1.x;
		double dy2= nd3.y - nd1.y;
		double dz2= nd3.z - nd1.z;
		double cross[3] = {dy*dz2 - dz*dy2, dx2*dz-dz2*dx, dx*dy2 - dx2*dy};
		area = (cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
		if(glSafe[n] && area > maxArea) { maxArea = area; n3 = n; }
	}
	if(maxArea > 0) {
		if(!solInfo.isCoupled || solInfo.getFetiInfo().fsi_corner == 0) // JLchange
			isCorner[n3] = true;
		else { if (!subPre->onWetInterface(n3)) isCorner[n3] = true; }
		weight[n3] = 3;
	}
	// fprintf(stderr, "Picked %d %d %d\n", n1,n2,n3);
	{
		//Node &nd1 = *nodes[n1];
		//Node &nd2 = *nodes[n2];
		//Node &nd3 = *nodes[n3];
		// fprintf(stderr, "N1 %e %e %e\n", nd1.x,nd1.y,nd1.z);
		// fprintf(stderr, "N2 %e %e %e\n", nd2.x,nd2.y,nd2.z);
		// fprintf(stderr, "N3 %e %e %e\n", nd3.x,nd3.y,nd3.z);
	}
}

// This routine internaly marks corner candidates
// and returns the numbers of corners in 2 numbers. The ones for which
// this subdomain is master and the number of corners as well as the number
// of corners already chosen and the total number of subdomains connected
// to all the corners (including candidate corners) of this subdomain.
// only safe corners can be marked
void
FetiSubCornerHandler::countAndMarkCornerCand(int *mync, int *totnc)
{
	bool gotEnough = false;
	for(int iNode = 0; iNode < nnodes; ++iNode) deg[iNode] = weight[iNode] = 0;

	for(int iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		for(int iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode)
			if(glSafe[sharedNodes[iNeighb][iNode]]) weight[sharedNodes[iNeighb][iNode]]++;
	}

	for(int iNode=0; iNode<nnodes; iNode++) {
		if(!glSafe[iNode]) continue;
		if(isCorner[iNode]) {
			gotEnough = true; // rotational corners remove singularity XXXX dangerous: also could be multi degree corner
		}
		else if(weight[iNode] > 2) {
			if(!solInfo.isCoupled || solInfo.getFetiInfo().fsi_corner == 0) // JLchange
				deg[iNode] = 1;
			else if(!subPre->onWetInterface(iNode))
				deg[iNode] = 1;
		}
	}

/* Colorado changes
// this array needs to be sized to allow for each face to add
 // 3 cross points
 int *currentC = (int *) dbg_alloca((numC+3*nNeighb)*sizeof(int));
 for(iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
    int numSharedNodes = sharedNodes.num(iNeighb);
    numC = 0;*/
	for(int iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		int numC = 0;
		int *currentC = (int *) alloca(sharedNodes.num(iNeighb)*sizeof(int));
		auto lSafe = isSafe.cbegin() + sharedNodes.offset(iNeighb);
		for(int iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode) {
			if(lSafe[iNode] && (deg[sharedNodes[iNeighb][iNode]] || isCorner[sharedNodes[iNeighb][iNode]])) {
				if(!solInfo.isCoupled || solInfo.getFetiInfo().fsi_corner == 0) { // JLchange
					currentC[numC] = sharedNodes[iNeighb][iNode];
					numC++;
				}
				else if(!subPre->onWetInterface(sharedNodes[iNeighb][iNode])) {
					currentC[numC] = sharedNodes[iNeighb][iNode];
					numC++;
				}
			}
		}
		bool needsMoreCrossPoints = checkForColinearCrossPoints(numC, currentC);
		if(needsMoreCrossPoints) {
			if(addPotCornerPoints(sharedNodes[iNeighb], lSafe))
				gotEnough = true;
		}
		else gotEnough = true;
	}

	mync[glSubNum] = 0; totnc[glSubNum] = 0;
	// use the weight array to see if this subdomain is master of each shared node (ie, have larger glSubNum than all neighbors)
	for(int iNode = 0; iNode < nnodes; ++iNode) weight[iNode] = 1;
	for(int iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		if(neighbSubs[iNeighb] < glSubNum)
			for(int iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode) {
				if(sharedNodes[iNeighb][iNode] >= nnodes) {
					std::cerr << "Shared node is out of bounds: " << sharedNodes[iNeighb][iNode] <<
					" vs " << nnodes << std::endl;
					throw  "Bad index";
				}

				weight[sharedNodes[iNeighb][iNode]] = 0;
			}
	}

	if(gotEnough == false && solInfo.getFetiInfo().pickAnyCorner) { // pick dim corners
		pickAnyCorners();
		// XXXX pickAnyCorners adjusts weight array but not deg ... is this ok?
		// what if 2 subdomains pick the same corner here ???
	}

	nTC = 0;
	for(int iNode = 0; iNode < nnodes; ++iNode) {
		if(!glSafe[iNode]) { weight[iNode] = -1; continue; } // PJSA don't mark unsafe nodes
		if(weight[iNode]) { // this subdomain is the master of this node
			if(isCorner[iNode] || deg[iNode]) weight[iNode] = nTC++;
			else weight[iNode] = -1;
		}
		else weight[iNode] = -1;
	}
	mync[glSubNum] = nTC;

	// Now count how many subs are connected to all the corners for which this subdomain is a master
	totnc[glSubNum] = nTC;
	for(int iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		for(int iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode) {
			if(weight[sharedNodes[iNeighb][iNode]] >= 0) totnc[glSubNum]++;
		}
	}
#ifdef DEBUG_CORNER
	cerr << "sub " << glSubNum << " is master of " << nTC << " candidate corners \n";
#endif
}

void
FetiSubCornerHandler::getCornerXYZ(int *allFCNum, double (*xyz)[3],
                               char *essential, int *cTsP, int *cTsT)
{
	int iNeighb, iNode;
	int firstCNum = allFCNum[glSubNum];
	double (*myXYZ)[3] = xyz + firstCNum;
	int *myCTsP = cTsP + firstCNum;
	int firstP = myCTsP[0];
	char *myEssent = essential + firstCNum;

	// Count how many subs each corner touches. Start by counting ourself
	for(iNode = 0; iNode < nTC; ++iNode) myCTsP[iNode] = 1;
	// By now, weight[iNode] is a local corner number for the corners
	// controlled by this subdomain
	for(iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		int numSharedNodes = sharedNodes.num(iNeighb);
		auto lSafe = isSafe.begin() + sharedNodes.offset(iNeighb);
		for(iNode = 0; iNode < numSharedNodes; ++iNode) {
			int nd = sharedNodes[iNeighb][iNode];
			if(lSafe[iNode] && weight[ nd ] >= 0) {
				myCTsP[weight[ nd ]]++;
				if(weight[ nd ] >= nTC) fprintf(stderr, "Inconsistency found in get Corner XYZ\n"); // XML
				myEssent[ weight[ nd ] ] = (isCorner[nd]) ? 1 : 0;
			}
		}
	}
	for(iNode = 0; iNode < nnodes; ++iNode) // some corners may not be in sharedNodes (eg pickAnyCorners)
		if(glSafe[iNode] && (weight[iNode] >= 0)) {
			myEssent[weight[iNode]] = (isCorner[iNode]) ? 1 : 0;
		}

	int ptr = firstP;
	for(iNode = 0; iNode < nTC; ++iNode) {
		ptr += myCTsP[iNode];
		myCTsP[iNode] = ptr;
	}

	// Now fill the targets
	for(iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		int numSharedNodes = sharedNodes.num(iNeighb);
		auto lSafe = isSafe.begin() + sharedNodes.offset(iNeighb);
		for(iNode = 0; iNode < numSharedNodes; ++iNode) {
			int nd = sharedNodes[iNeighb][iNode];
			if(lSafe[iNode] && weight[ nd ] >= 0) {
				cTsT[ --myCTsP[weight[ nd ]] ] = neighbSubs[iNeighb];
			}
		}
	}

	for(iNode = 0; iNode < nTC; ++iNode)
		cTsT[ --myCTsP[iNode] ] = glSubNum;
	for(iNode = 0; iNode < nnodes; ++iNode) {
		if((weight[iNode] >= 0) && glSafe[iNode]) {
			int nd = weight[iNode];
			myXYZ[nd][0] = nodes[iNode]->x;
			myXYZ[nd][1] = nodes[iNode]->y;
			myXYZ[nd][2] = nodes[iNode]->z;
			//if(isCorner[iNode] && dsa[iNode].contains(DofSet::XYZrot))
			//numRotCrn[glSubNum]++;
		}
	}
}

bool
FetiSubCornerHandler::addPotCornerPoints(gsl::span<const int> allNodes, std::vector<bool>::const_iterator safe)
{
	int numShared = allNodes.size();
	if(allNodes.size() < dim) {
		for(int i = 0; i < numShared; ++i)
			if(safe[i]) deg[allNodes[i]] = 3;
		return false;
	}
	int first = 0;
	while(first < numShared && safe[first] == false) first++;
	if(first+dim > numShared) {
		return false;
	}

	int n1 = allNodes[first];
	double maxDist = 0.0, dist;
	Node &nd1 = *nodes[n1];
	int n2 = n1;
	int n;
	// find the most distant node from node 1 (second node)
	double roundoff_tol = 1.0e-6;  // PJSA: this helps to make corners same on different platforms
	for(n=first+1; n<numShared; ++n)
		if(allNodes[n] < nnodes) {
			Node &nd2 = *nodes[allNodes[n]];
			dist = nd2.distance(nd1);
			if(safe[n] && dist > maxDist*(1.0+roundoff_tol)) { maxDist = dist; n2 = allNodes[n]; } // PJSA
		}
	if(maxDist > 0 && dim <= 2) {  // PJSA
		deg[n1] = 3;
		deg[n2] = 3;
		return true;
	}
	Node &nd2 = *nodes[n2];
	double dx= nd2.x - nd1.x;
	double dy= nd2.y - nd1.y;
	double dz= nd2.z - nd1.z;

	// finally, find the node that maximizes the area with nodes 1 & 2
	double maxArea=0.0, area;
	int n3 = n1;
	for(n=0; n<numShared; ++n)
		if(allNodes[n] < nnodes) {
			Node &nd3 = *nodes[allNodes[n]];
			double dx2= nd3.x - nd1.x;
			double dy2= nd3.y - nd1.y;
			double dz2= nd3.z - nd1.z;
			double cross[3] = {dy*dz2 - dz*dy2, dx2*dz-dz2*dx, dx*dy2 - dx2*dy};
			area = (cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
			if(safe[n] && area > maxArea*(1.0+roundoff_tol)) { maxArea = area; n3 = allNodes[n]; }  // PJSA
		}
	deg[n1] = 3;
	deg[n2] = 3;
	//cerr << (dsa[n1] == DofSet::Helm) << " " << (dsa[n2] == DofSet::Helm) << " " << (dsa[n3] == DofSet::Helm) << endl;
	if(maxArea > 0) {
		deg[n3] = 3;
		return true;
	} else {
		return false;
	}
	// XXXX for mixed subdomain add three structure potential corners and one fluid
}

bool
FetiSubCornerHandler::checkForColinearCrossPoints(int numCornerPoints, int *localCornerPoints)
{
	// If we are 3D we need 3 points, 2D, two, 1d, 1... hence:
	if(numCornerPoints < dim) {
		return true;
	}
	if(dim == 1) {
		return false;
	}

	//int *allNodes = sharedNodes[0];

	// find the most distant cross point from the first cross point.
	int n1 = localCornerPoints[0];
	Node &nd1 = *nodes[n1];
	int n, n2;
	double maxDist=0.0, dist;
	for(n=1; n<numCornerPoints; ++n) {
		Node &nd2 = *nodes[localCornerPoints[n]];
		dist = nd2.distance(nd1);
		if(dist > maxDist) { maxDist = dist; n2 = localCornerPoints[n]; }
	}
	Node &nd2 = *nodes[n2];
	if(maxDist==0.0) {
		fprintf(stderr,"Coincident Cross Points! Very Very Bad!\n");
		for(n=0; n<numCornerPoints; ++n) {
			Node *nd = nodes[ localCornerPoints[n] ];
			fprintf(stderr, "%d  %e %e %e\n",
			        localCornerPoints[n], nd->x, nd->y, nd->z);
		}
	}
	if(dim < 3) {
		return false;
	}

	// loop over all cross points checking for colinearity between 3 points
	double v1[3];
	double v2[3];
	double v3[3];
	v1[0] = nd2.x - nd1.x;
	v1[1] = nd2.y - nd1.y;
	v1[2] = nd2.z - nd1.z;
	double v1Norm = (v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);

	bool answer = true;
	for(n=1; n<numCornerPoints; ++n) {
		if(n==n2) continue;
		//int n3 = localCornerPoints[n];
		Node &nd3 = *nodes[localCornerPoints[n]];
		v2[0] = nd3.x - nd1.x;
		v2[1] = nd3.y - nd1.y;
		v2[2] = nd3.z - nd1.z;
		crossprod(v1,v2,v3);
		if(magnitude(v3) >= 0.01*v1Norm) {
			//cerr << (dsa[n1] == DofSet::Helm) << " " << (dsa[n2] == DofSet::Helm) << " " << (dsa[n3] == DofSet::Helm) << endl;
			answer = false;
			break;
		}
	}
	// XXXX for mixed subdomain, check if subdomain has 3 structure nodes which are not colinear and one fluid node
	return answer;
}

void
FetiSubCornerHandler::countContact(int * cntNum, char *crnMrk)
{
	int iNode, iNeighb;
	int cnt = 0;
	bool *isContact = (bool *) dbg_alloca(sizeof(bool)*nnodes);
	for(iNode = 0; iNode < nnodes; ++iNode)
		isContact[iNode] = false;

	for(iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		auto lSafe = isSafe.begin() + sharedNodes.offset(iNeighb);
		for(iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode)
			if(lSafe[iNode] == false)
				isContact[sharedNodes[iNeighb][iNode]] = true;
	}

	for(iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		// only count the contacts for which I am master
		// eliminate the possible contact that are accounted for by a neighbor
		if(glSubNum > neighbSubs[iNeighb])
			for(iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode) {
				isContact[sharedNodes[iNeighb][iNode]] = false;
			}
		else
			for(iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode) {
				int node = deg[ sharedNodes[iNeighb][iNode] ];
				if(node >= 0)
					if(crnMrk[node] == 1) {
						isContact[sharedNodes[iNeighb][iNode]] = false;
					}
			}
	}

	for(iNode = 0; iNode < nnodes; ++iNode)
		if(isContact[iNode]) cnt++;
	cntNum[glSubNum] = cnt;
}

void
FetiSubCornerHandler::dispatchNumbering(FSCommPattern<int> *pat, char *crnMrk,
                                    int *allOrigFC, int *allNewFC, int, int *cntOff)
{
	// Dispatching corner numbering and add the contact numbers
	int origFirstC = allOrigFC[glSubNum];
	int newFirstC = allNewFC[glSubNum];
	// fprintf(stderr, "%d numberings %d and %d with cnt %d\n", glSubNum, origFirstC, newFirstC, cntOff[glSubNum]);
	// fprintf(myOut, "My first Cs are %d and %d\n", origFirstC, newFirstC);
	int crnNum = newFirstC;
	int iNode, iNeighb;
	int cnt = 0;
	//int cnt2 = 0;

	// It is important that this array be built before 'deg' is touched!!!
	bool *isContact = (bool *) dbg_alloca(sizeof(bool)*nnodes);
	for(iNode = 0; iNode < nnodes; ++iNode)
		isContact[iNode] = false;

	for(iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		auto lSafe = isSafe.begin() + sharedNodes.offset(iNeighb);
		for(iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode)
			if(lSafe[iNode] == false)
				isContact[sharedNodes[iNeighb][iNode]] = true;
	}

	for(iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		// only count the contacts for which I am master
		// eliminate the possible contact that are accounted for by a neighbor
		if(glSubNum > neighbSubs[iNeighb])
			for(iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode) {
				isContact[sharedNodes[iNeighb][iNode]] = false;
			}
		else
			for(iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode) {
				int node = deg[ sharedNodes[iNeighb][iNode] ];
				if(node >= 0)
					if(crnMrk[node] == 1) {
						isContact[sharedNodes[iNeighb][iNode]] = false;
					}
			}
	}

	for(iNode = 0; iNode < nnodes; ++iNode) {
		int node = deg[iNode];
		if(node >= origFirstC && crnMrk[node] == 1) {
			deg[iNode] = crnNum++;
		} else
			deg[iNode] = -1;
	}

	/* Why does this not give the same result???
	int cnt2 = 0;
	for(iNode = 0; iNode < nnodes; ++iNode)
	  if( isContact[iNode]) {
		int node = deg[ iNode ];
		if(node >= 0 && crnMrk[ node ] == 1){
		  isContact[iNode] = false;
		  cnt2++;
		}
	  }
	  */

	for(iNode = 0; iNode < nnodes; ++iNode)
		if(isContact[iNode]) {
			deg[ iNode ] = cntOff[glSubNum] + cnt;
			cnt++;
		}

	for(iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		FSSubRecInfo<int> rInfo = pat->getSendBuffer(glSubNum, neighbSubs[iNeighb]);
		for(iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode)
			rInfo.data[iNode] = deg[ sharedNodes[iNeighb][iNode] ];
	}
}

void
FetiSubCornerHandler::dispatchInitialNumbering(FSCommPattern<int> *pat, int *firstC)
{
	int myFirstC = firstC[glSubNum];
	int iNode, iNeighb;
	for(iNode = 0; iNode < nnodes; ++iNode)
		deg[iNode] = (weight[iNode] >= 0) ? weight[iNode]+myFirstC : -1;
	for(iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		FSSubRecInfo<int> rInfo = pat->getSendBuffer(glSubNum, neighbSubs[iNeighb]);
		for(iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode)
			rInfo.data[iNode] = deg[ sharedNodes[iNeighb][iNode] ];
	}
}

void
FetiSubCornerHandler::recInitialNumbering(FSCommPattern<int> *pat, int *numRotCrn)
{
	int iNeighb, iNode;
	numRotCrn[glSubNum] = 0;
	for(iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		FSSubRecInfo<int> ri = pat->recData(neighbSubs[iNeighb], glSubNum);
		for(iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode)
			if(ri.data[iNode] >= 0) {
				if(deg[ sharedNodes[iNeighb][iNode] ] >= 0)
					fprintf(stderr, "(a) Found a duplicate corner\n");
				deg[ sharedNodes[iNeighb][iNode] ] = ri.data[iNode];
			}
	}
	for(iNode = 0; iNode < nnodes; ++iNode)
		//  GR: include end points of all subdomain edges as corner nodes for the 2D-case
		if(deg[iNode] >= 0 && isCorner[iNode] &&
		   (dsa[iNode].contains(DofSet::XYZrot) || (dsa[iNode] == XYDofs) || (dsa[iNode] == DofSet::Helm) || (dsa[iNode] == DofSet::Temp)))
			numRotCrn[glSubNum]++;
}

void
FetiSubCornerHandler::recNumbering(FSCommPattern<int> *pat, int *fMaster)
{
	int iNeighb, iNode;
	for(iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		FSSubRecInfo<int> ri = pat->recData(neighbSubs[iNeighb], glSubNum);
		for(iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode) {
			if(ri.data[iNode] >= 0) {
				if(deg[ sharedNodes[iNeighb][iNode] ] >= 0)
					fprintf(stderr, "(b) Found a duplicate corner\n");
				deg[ sharedNodes[iNeighb][iNode] ] = ri.data[iNode];
			}
		}
	}
	if(!solInfo.isCoupled || solInfo.getFetiInfo().fsi_corner == 0) { // JLchange
		totNC = 0;
		for(iNode = 0; iNode < nnodes; ++iNode)
			if(deg[iNode] >= 0)
				totNC++;
		crnList.reserve(totNC);
		for(iNode = 0; iNode < nnodes; ++iNode) {
			if(deg[iNode] >= 0)
				crnList.push_back(iNode);
		}
	}
	else {
		totNC = 0;
		for(iNode = 0; iNode < nnodes; ++iNode) {
			if ( ( (deg[iNode] >= 0) && !(subPre->onWetInterface(iNode)) )
			    || subPre->isWetInterfaceCorner(iNode))
				totNC++;
		}
		crnList.reserve(totNC);
		for(iNode = 0; iNode < nnodes; ++iNode) {
			if( ( (deg[iNode] >= 0) && !(subPre->onWetInterface(iNode)) )
			    || subPre->isWetInterfaceCorner(iNode))
				crnList.push_back(iNode);

		}
	}
}

void
FetiSubCornerHandler::listRotCorners(int *fN, int *crnNum)
{
	int *myNum = crnNum+fN[glSubNum];
	int iNode, iCrn = 0;
	for(iNode = 0; iNode < nnodes; ++iNode)
		// GR: include end points of all subdomain edges as corner nodes for the 2D-case
		if(deg[iNode] >= 0 && isCorner[iNode] &&
		   (dsa[iNode].contains(DofSet::XYZrot) || (dsa[iNode] == XYDofs) || (dsa[iNode] == DofSet::Helm) || (dsa[iNode] == DofSet::Temp))) {
			myNum[iCrn++] = deg[iNode];
		}
}

void
FetiSubCornerHandler::resendNumbers(FSCommPattern<int> *pat)
{
	int iNeighb, iNode;
	for(iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		FSSubRecInfo<int> rInfo = pat->getSendBuffer(glSubNum, neighbSubs[iNeighb]);
		for(iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode)
			rInfo.data[iNode] = deg[ sharedNodes[iNeighb][iNode] ];
	}
}

void
FetiSubCornerHandler::checkNumbers(FSCommPattern<int> *pat)
{
	int iNeighb, iNode;
	for(iNeighb = 0; iNeighb < nNeighb; ++iNeighb) {
		FSSubRecInfo<int> rInfo = pat->recData(neighbSubs[iNeighb], glSubNum);
		for(iNode = 0; iNode < sharedNodes.num(iNeighb); ++iNode)
			if(rInfo.data[iNode] != deg[ sharedNodes[iNeighb][iNode] ])
				fprintf(stderr, "We disagree: %d %d\n", rInfo.data[iNode],
				        deg[ sharedNodes[iNeighb][iNode] ]);
	}
}

void
FetiSubCornerHandler::markDims(int *_dims)
{
	for(int i=0; i<4; ++i)
		if(dims[i]) _dims[i] = 1;
}

CornerSelector::CornerSelector(int _glNumSub, int _nSub, std::vector<FetiSubCornerHandler *> _cornerHandler,
                               FSCommPattern<int> *_cpat, FSCommunicator *_communicator) :
		cornerHandler(std::move(_cornerHandler))
{
	glNumSub = _glNumSub;
	nSub = _nSub;
	communicator = _communicator;
	cpat = _cpat;

	for(int i=0; i<4; ++i) dims[i] = 0;
	paralApply(nSub, cornerHandler.data(), &FetiSubCornerHandler::markDims, dims);
#ifdef DISTRIBUTED
	communicator->globalMax(4, dims);
#endif
	int sdim = dims[0]+dims[1]+dims[2]; // structure
	int fdim = dims[3]; // fluid
	bool mixed = (fdim && sdim);
	dim = (mixed || (fdim==0)) ? sdim : fdim; // use the structure dimension for a mixed subdomain XXXX
#ifdef DEBUG_CORNER
	cerr << "global dim = " << dim << endl;
#endif
}

CornerSelector::~CornerSelector()
{
	delete [] glSubGroup;
	// don't delete grToSub (GenDecDomain)
}

int
CornerSelector::makeCorners()
{
	int i, iSub;

	// First locate ``unsafe'' nodes. These nodes must be corner
	// nodes (to eliminate subdomain ZEMs) but do not work in tying two subdomains together
	paralApply(nSub, cornerHandler.data(), &FetiSubCornerHandler::dispatchSafeNodes, cpat);
	cpat->exchange();
	paralApply(nSub, cornerHandler.data(), &FetiSubCornerHandler::markSafeNodes, cpat);

	// Take care of the multi degreed nodes that could upset the coarse problem
	paralApply(nSub, cornerHandler.data(), &FetiSubCornerHandler::markMultiDegNodes);

	// Take care of the corners for 4th order problems
	paralApply(nSub, cornerHandler.data(), &FetiSubCornerHandler::dispatchRotCorners, cpat);
	cpat->exchange();
	paralApply(nSub, cornerHandler.data(), &FetiSubCornerHandler::markRotCorners, cpat);

	// Build the list of potential corner and mark those already certainly chosen
	// Phase 1, count how many corners each sub has and to how many subdomains these corners connect
	int *cPerSub = new int[glNumSub+1];
	int *nTot = new int[glNumSub];
	for(i = 0; i < glNumSub; ++i) cPerSub[i] = nTot[i] = 0;

	// countAndMarkCornerCand internaly marks corner candidates and returns the numbers of corners
	// for which this subdomain is master and the total number of corners connectivities
	// weight[inode] contains the local numbering of corner candidates for which this sub is a master
	paralApply(nSub, cornerHandler.data(), &FetiSubCornerHandler::countAndMarkCornerCand, cPerSub, nTot);
	communicator->globalSum(glNumSub, cPerSub);
	communicator->globalSum(glNumSub, nTot);

	int totNC, tmpNC;
	totNC = 0;
	// start creating the corner candidate to sub (master) connectivity
	for(i = 0; i < glNumSub; ++i) {
		tmpNC = cPerSub[i];
		cPerSub[i] = totNC;
		totNC += tmpNC;
	}
	cPerSub[glNumSub] = totNC;

#ifdef DEBUG_CORNER
	cerr << "Number of corner candidates = " << totNC << endl;
#endif

	int *cPtr = new int[totNC+1];
	int tot = 0;
	for(i = 0; i < glNumSub; ++i) {
		for(int j = cPerSub[i]; j < cPerSub[i+1]; ++j) cPtr[j] = tot;
		tot += nTot[i];
	}
	cPtr[cPerSub[glNumSub]] = tot;

	int *cTg = new int[tot];
	double (*xyz)[3] = new double[tot][3];
	int *essentPtr = new int[glNumSub+1];
	for(i = 0; i < tot; ++i) {
		cTg[i] = 0;
		xyz[i][0] = xyz[i][1] = xyz[i][2] = 0.0;
	}
	for(i = 0; i < glNumSub; ++i) essentPtr[i] = 0;

	char *essential = new char[totNC];
	for(i = 0; i < totNC; ++i) essential[i] = 0;

	// getCornerXYZ collects the coordinates of the corner candidates and in addition lets the
	// caller know which ones are already essential corners and finally fills the corner to sub
	// connectivity data and count the number of corner with rot dofs
	paralApply(nSub, cornerHandler.data(), &FetiSubCornerHandler::getCornerXYZ,
	           cPerSub, xyz, essential, cPtr, cTg);
	communicator->globalMax(totNC, cPtr);
	communicator->globalSum(3*tot, (double *) xyz);
	communicator->globalSum(tot, cTg);
#if defined(LAM_MPI) || defined(SUN10) || defined(OPEN_MPI)
	// temporary fix, LAM 7.1.1 MPI_Allreduce isn't working with character so converting to int
  int *essential_tmp = new int[totNC];
  for(i = 0; i < totNC; ++i) essential_tmp[i] = int(essential[i]);
  communicator->globalSum(totNC, essential_tmp);
  for(i = 0; i < totNC; ++i) essential[i] = char(essential_tmp[i]);
  delete [] essential_tmp;
#else
	communicator->globalSum(totNC, essential);
#endif
	// now we have the "corner candidates" to subdomain connectivity
	Connectivity *cToSub = new Connectivity(totNC, cPtr, cTg);

	// exchange initial global numbering of the corner candidates (because we need
	// it to build the "rotational corner" connectivity
	paralApply(nSub, cornerHandler.data(), &FetiSubCornerHandler::dispatchInitialNumbering, cpat, cPerSub);
	cpat->exchange();
	paralApply(nSub, cornerHandler.data(), &FetiSubCornerHandler::recInitialNumbering, cpat, essentPtr);
	communicator->globalSum(glNumSub, essentPtr);

	int cNum = 0;
	for(iSub = 0; iSub < glNumSub; ++iSub) {
		int tmp = cNum;
		cNum += essentPtr[iSub];
		essentPtr[iSub] = tmp;
	}
	essentPtr[glNumSub] = cNum;

	int *essentTg = new int[cNum];
	for(i = 0; i < cNum; ++i) essentTg[i] = 0;

	paralApply(nSub, cornerHandler.data(),&FetiSubCornerHandler::listRotCorners, essentPtr, essentTg);
	communicator->globalSum(cNum, essentTg);

	Connectivity subToRotCrn(glNumSub, essentPtr, essentTg);

	// initialize glCrnGroup and glSubGroup
	int *glCrnGroup = new int[totNC];
	for(i = 0; i < totNC; ++i) glCrnGroup[i] = -1;
	glSubGroup = new int[glNumSub];
	for(i = 0; i < glNumSub; ++i) glSubGroup[i] = -1;

	// Now chose the corners
	chooseCorners(essential, xyz, *cToSub, subToRotCrn, glCrnGroup);

#ifdef DEBUG_CORNER
	int xx = 0;
  for(i = 0; i < totNC; ++i) if(essential[i]) xx++;
  fprintf(stderr, "Found %d essential safe corners\n", xx);
#endif

	// Count contact corners
	int *numCnt = new int[glNumSub];
	for(iSub = 0; iSub < glNumSub; ++iSub)
		numCnt[iSub] = 0;

	paralApply(nSub, cornerHandler.data(),&FetiSubCornerHandler::countContact,
	           numCnt, essential);
	communicator->globalSum(glNumSub, numCnt);

	int *newCPerSub = new int[glNumSub];
	int cCount = 0;
	for(iSub = 0; iSub < glNumSub; ++iSub) {
		newCPerSub[iSub] = cCount;
		for(int jj = cPerSub[iSub]; jj < cPerSub[iSub+1]; ++jj)
			if(essential[jj]) cCount++;
	}

	int cntCount = cCount;
	for(iSub = 0; iSub < glNumSub; ++iSub) {
		int tmpCntCount = numCnt[iSub];
		numCnt[iSub] = cntCount;
		cntCount += tmpCntCount;
	}
#ifdef DEBUG_CORNER
	fprintf(stderr, "Found %d essential safe+unsafe corners\n", cntCount);
#endif
	paralApply(nSub, cornerHandler.data(), &FetiSubCornerHandler::dispatchNumbering,
	           cpat, essential, cPerSub, newCPerSub, cCount, numCnt);
	cpat->exchange();

	paralApply(nSub, cornerHandler.data(), &FetiSubCornerHandler::recNumbering,
	           cpat, newCPerSub);

	paralApply(nSub, cornerHandler.data(), &FetiSubCornerHandler::resendNumbers, cpat);
	cpat->exchange();
	paralApply(nSub, cornerHandler.data(), &FetiSubCornerHandler::checkNumbers, cpat);

	delete [] essential;
	delete [] numCnt;
	delete [] newCPerSub;
	delete [] nTot;
	delete [] cPerSub;
	delete [] glCrnGroup;
	delete cToSub;

	// return the total number of corners including the ones for contact
	return cntCount;
}

void
CornerSelector::chooseCorners(char *glCornerList, double (*crnXYZ)[3],
                           Connectivity &cNConnect, Connectivity &subToRotCrn,
                           int *glCrnGroup)
{
	Connectivity subToCrn = cNConnect.reverse();
	Connectivity rotCrnToSub = subToRotCrn.reverse();
	// fprintf(stderr, "We have %d rot corners\n", rotCrnToSub->csize());
	Connectivity rotSubToSub = subToRotCrn.transcon(rotCrnToSub);
	grToSub = std::make_unique<Connectivity> ( rotSubToSub.connexGroups() );

	double aamin = DBL_MAX, aamax = -DBL_MAX;

	// guarantee we will loop at least once
	int lastGrSize = grToSub->csize()+1;
	while(grToSub->csize() > 1 && grToSub->csize() != lastGrSize) {
		lastGrSize = grToSub->csize();
		int iGr;
		// fprintf(stderr, "Group has size %d\n", grToSub->csize());
		Connectivity grToCrn = grToSub->transcon(subToCrn);
		Connectivity crnToGr = grToCrn.reverse();
		// Eliminate the corners that are not corners anymore
		//grToCrn = grToCrn.trim(crnToGr);
		grToCrn.sortTargets();
		Connectivity grToGr = grToCrn.transcon(crnToGr);
		// Now examine each group and establish a set of connection
		// priorities
		int (*fav)[2] = new int[grToSub->csize()][2];
		double (*bamax)[2] = new double[grToSub->csize()][2];
		int (*choices)[2][3] = new int[grToSub->csize()][2][3];
		int *finalFav = new int[grToSub->csize()];
		bool *tied = new bool[grToSub->csize()];

		for(iGr = 0; iGr < grToCrn.csize(); ++iGr) {
			tied[iGr] = false;
			int jGr;
			finalFav[iGr] = fav[iGr][0] = fav[iGr][1]  = -1;
			double tamax = 0, grAmax = 0.0;
			int ntmax = -1;

			//double amax = 0, amax2 = 0;
			bamax[iGr][0] = bamax[iGr][1] = 0.0;
			for(jGr = 0; jGr < grToGr.num(iGr); ++jGr) {
				int grJ = (grToGr)[iGr][jGr];
				if(grJ == iGr) continue;
				// Find the comon nodes
				int maxCm = std::min( grToCrn.num(iGr), grToCrn.num(grJ) );
				// XML this needs to be changed
				std::vector<int>cnode(maxCm, -1);
				int nc = 0;
				int n1 = 0, n2 = 0;
				while(n1 < grToCrn.num(iGr) && n2 < grToCrn.num(grJ)) {
					if( (grToCrn)[iGr][n1] == grToCrn[grJ][n2] ) {
						cnode[nc++] = grToCrn[iGr][n1];
						n1++; n2++;
					} else {
						if( grToCrn[iGr][n1] < grToCrn[grJ][n2] )
							n1++;
						else
							n2++;
					}
				}
				if(nc < dim) continue;
				int iC;
				// pick the candidate with most number of groups connected to it
				// Or a candidate already chosen
				int bc = cnode[0];
				for(iC = 1; iC < nc; ++iC)
					if(glCornerList[cnode[iC]]) {
						if(glCornerList[bc] == 0 ||
						   crnToGr.num(cnode[iC]) > crnToGr.num(bc))
							bc = cnode[iC];
					} else if(glCornerList[bc] == 0
					          && crnToGr.num(cnode[iC]) > crnToGr.num(bc))
						bc = cnode[iC];

				// Now identify the furthest nodes from one another
				double mxDist2 = 0.0, mxDist2b = 0.0;
				int c2 = -1;
				int c2b = -1;
				for(iC = 0; iC < nc; ++iC) {
					int ndI = cnode[iC];
					int ndB = bc;
					double dst2 =
							(crnXYZ[ndB][0]-crnXYZ[ndI][0])*
							(crnXYZ[ndB][0]-crnXYZ[ndI][0]) +
							(crnXYZ[ndB][1]-crnXYZ[ndI][1])*
							(crnXYZ[ndB][1]-crnXYZ[ndI][1]) +
							(crnXYZ[ndB][2]-crnXYZ[ndI][2])*
							(crnXYZ[ndB][2]-crnXYZ[ndI][2]);
					if(dst2 > mxDist2) {
						mxDist2 = dst2;
						c2 = cnode[iC];
					}
					if(glCornerList[ndI] && dst2 > mxDist2b) {
						mxDist2b = dst2;
						c2b = ndI;
					}
				}
				// XML put it back
				if(mxDist2b >= 1.6 * mxDist2)
					c2 = c2b;
				if(c2 < 0) {
					if(nc > 1)
						fprintf(stderr, "no line %d %e %e, %e %e\n", nc,
						        crnXYZ[cnode[0]][0],
						        crnXYZ[cnode[1]][0], crnXYZ[cnode[0]][1],
						        crnXYZ[cnode[1]][1]);
					else
						fprintf(stderr,"no line %d\n",nc);
					continue; }
				// finally, find the node that maximizes the area with nodes 1 & 2
				int c3 = -1;
				if(dim >= 3) {
					double maxArea=0.0, area;
					double dx = crnXYZ[bc][0] - crnXYZ[c2][0];
					double dy = crnXYZ[bc][1] - crnXYZ[c2][1];
					double dz = crnXYZ[bc][2] - crnXYZ[c2][2];
					for(iC = 0; iC < nc; ++iC) {
						int ndI = cnode[iC];
						int ndC = c2;
						double dx2= crnXYZ[ndC][0] - crnXYZ[ndI][0];
						double dy2= crnXYZ[ndC][1] - crnXYZ[ndI][1];
						double dz2= crnXYZ[ndC][2] - crnXYZ[ndI][2];
						double cross[3] = {dy*dz2 - dz*dy2, dx2*dz-dz2*dx, dx*dy2 - dx2*dy};
						area = (cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
						if(area > maxArea) { maxArea = area; c3 = cnode[iC]; }
						//if(area > maxAreab) { maxAreab = area; c3b = cnode[iC]; }
					}

					if(maxArea > bamax[iGr][0]) {
						fav[iGr][1] = fav[iGr][0];
						bamax[iGr][1] = bamax[iGr][0];
						choices[iGr][1][0] = choices[iGr][0][0];
						choices[iGr][1][1] = choices[iGr][0][1];
						choices[iGr][1][2] = choices[iGr][0][2];
						fav[iGr][0] = grJ;
						choices[iGr][0][0] = bc;
						choices[iGr][0][1] = c2;
						choices[iGr][0][2] = c3;
						bamax[iGr][0] = maxArea;
					} else if(maxArea > bamax[iGr][1]) {
						fav[iGr][1] = grJ;
						choices[iGr][1][0] = bc;
						choices[iGr][1][1] = c2;
						choices[iGr][1][2] = c3;
						bamax[iGr][1] = maxArea;
					}
				} else {
					if(mxDist2 > bamax[iGr][0]) {
						fav[iGr][1] = fav[iGr][0];
						bamax[iGr][1] = bamax[iGr][0];
						choices[iGr][1][0] = choices[iGr][0][0];
						choices[iGr][1][1] = choices[iGr][0][1];
						choices[iGr][1][2] = choices[iGr][0][2];
						bamax[iGr][0] = mxDist2;
						fav[iGr][0] = grJ;
						choices[iGr][0][0] = bc;
						choices[iGr][0][1] = c2;
						choices[iGr][0][2] = -1;
					}
				}

				// Now double check if we already had the domain tied
				int ncc = 0;
				for(iC = 0; iC < nc; ++iC)
					if(glCornerList[cnode[iC]])
						cnode[ncc++] = cnode[iC];
				if(ncc >= 3) {
					double area = checkArea(ncc, cnode, crnXYZ);
					area *= area;
					//fprintf(stderr, "Area pretied is %e vs %e %e\n", area, amax, maxArea);
					if(area > tamax) {
						ntmax = grJ;
						tamax = area;
					}
					// Keeping track of the largest area for this group
				}
				if(bamax[iGr][0] > grAmax)
					grAmax = bamax[iGr][0];
			}
			if(tamax > 0.4*grAmax) {
				//fprintf(stderr, "Doing a tie\n");
				finalFav[iGr] = ntmax;
				tied[iGr] = true;
				tied[ntmax] = true;
			} //else
			//fprintf(stderr, "No tie %e %e\n", tamax, grAmax);
			if(grAmax > aamax)
				aamax = grAmax;
			if(0.6*grAmax < aamin)
				aamin = 0.6*grAmax;
		}
		// Make pairs of groups
		for(iGr = 0; iGr < grToCrn.csize(); ++iGr) {
			if(tied[iGr]) continue;
			if(fav[iGr][0] < 0) {
				//fprintf(stderr, "Warning: Group %d cannot be tied\n", iGr);
				int in;
				for(in = 0; in < grToSub->num(iGr); ++in) {
//           fprintf(stderr, "Sub %d\n", (*grToSub)[iGr][in]);
					glSubGroup[(*grToSub)[iGr][in]] = iGr;
				}
				for(in = 0; in < grToCrn.num(iGr); ++in) {
//           fprintf(stderr, "Corner %d\n", grToCrn[iGr][in]);
					glCrnGroup[grToCrn[iGr][in]] = iGr;
				}
			} else {
				// check if one of our favorites is not tied yet
				int pref = -1;
				if(tied[ fav[iGr][0] ] == false)
					pref = 0;
				else if(fav[iGr][1] >= 0 && tied[ fav[iGr][1] ] == false) {
					if(bamax[iGr][1] > 0.25*bamax[iGr][0])
						pref = 1;
					else
						pref = 0;
				}
				if(pref >= 0) {
					glCornerList[choices[iGr][pref][0]] = 1;
					if(choices[iGr][pref][1] > -1) glCornerList[choices[iGr][pref][1]] = 1; // PJSA
					if(choices[iGr][pref][2] > -1) glCornerList[choices[iGr][pref][2]] = 1; // PJSA 1-16-07
					tied[ fav[iGr][pref] ] = true;
					tied[ iGr ] = true;
					finalFav[iGr] = fav[iGr][pref];
				}
			}
		}
		// now attach the groups that had not been attached yet
		for(iGr = 0; iGr < grToCrn.csize(); ++iGr) {
			if(tied[iGr]) continue;
			if(fav[iGr][0] >= 0) {
				glCornerList[choices[iGr][0][0]] = 1;
				if(choices[iGr][0][1] > -1) glCornerList[choices[iGr][0][1]] = 1;
				if(choices[iGr][0][2] > -1) glCornerList[choices[iGr][0][2]] = 1; // PJSA 1-16-07
				finalFav[iGr] = fav[iGr][0];
			}
		}
		// Now create a new group to group connectivity
		int *count = new int[grToCrn.csize()+1];
		for(iGr = 0; iGr < grToCrn.csize(); ++iGr)
			count[iGr] = 1;
		count[grToCrn.csize()] = 0;
		for(iGr = 0; iGr < grToCrn.csize(); ++iGr)
			if(finalFav[iGr] >= 0)
				count[finalFav[iGr]]++;
		for(iGr = 0; iGr < grToCrn.csize(); ++iGr)
			count[iGr+1] += count[iGr];
		int *target = new int[count[grToCrn.csize()] ];
		for(iGr = 0; iGr < grToCrn.csize(); ++iGr) {
			target[--count[iGr]] = finalFav[iGr];
			if(finalFav[iGr] >= 0)
				target[--count[ finalFav[iGr] ] ] = iGr;
		}
		Connectivity *gr2Gr = new Connectivity(grToCrn.csize(), count, target);

		// Find the disconnected parts as separate groups.
		Connectivity grToPrevGr = gr2Gr->connexGroups();
		// get the group To Sub graph
		grToSub = std::make_unique<Connectivity>(grToPrevGr.transcon(*grToSub));
		delete [] finalFav;
		delete [] tied;
		delete [] choices;
		delete [] fav;
		delete [] bamax;
		delete gr2Gr;
		
	}
	//communicator->sync();
	//fprintf(stderr, "Possible extremes: %e %e\n",aamax, aamin);
	delete [] crnXYZ;


}

