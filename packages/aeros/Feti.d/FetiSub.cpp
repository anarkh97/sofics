//
// Created by Michel Lesoinne on 12/6/17.
//
#include <Solvers.d/Rbm.h>
#include <Driver.d/SubDomain.h>
#include <Solvers.d/SolverFactory.h>
#include <Element.d/Element.h>
#include <complex>
#include <Eigen/Sparse>
#include <Utils.d/dofset.h>
#include "FetiSub.h"
#include <Driver.d/Mpc.h>
#include <Utils.d/SolverInfo.h>
#include <Utils.d/dbg_alloca.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/SparseSet.h>
#include <Math.d/CuCSparse.h>
#include <Solvers.d/Solver.h>
#include <Math.d/BLAS.h>
#include <Driver.d/HData.h>


const ConstrainedDSA * FetiBaseSub::getCCDSA() const {
	return cc_dsa ? cc_dsa.get() : get_c_dsa();
}


void
FetiBaseSub::markCornerDofs(gsl::span<int> glCornerDofs) const
{
	for(int i=0; i<numCRN; ++i)
		glCornerDofs[glCornerNodes[i]] |= cornerDofs[i].list();
}

void
FetiBaseSub::makeKccDofs(DofSetArray *cornerEqs, int augOffset,
                                  Connectivity *subToEdge, int mpcOffset) {
	int numC = numCoarseDofs();
	cornerEqNums.resize(numC);

	// numbers the corner equations
	int offset = 0;
	for (int i = 0; i < numCRN; ++i)
		offset += cornerEqs->number(glCornerNodes[i], cornerDofs[i].list(), cornerEqNums.data() + offset);

	// number the mpc equations
	for (int i = 0; i < numMPC_primal; ++i) {
		int fDof = cornerEqs->firstdof(mpcOffset + localToGlobalMPC_primal[i]);
		cornerEqNums[offset++] = fDof;
	}

	// number the augmentation equations
	if (getFetiInfo().augment == FetiInfo::Gs) {
		int fDof = cornerEqs->firstdof(augOffset + subNum());
		for (int i = 0; i < nGrbm; ++i)
			cornerEqNums[offset++] = fDof + i;
		for (int iNeighb = 0; iNeighb < scomm->numNeighb; ++iNeighb) {
			int fDof = cornerEqs->firstdof(augOffset + scomm->subNums[iNeighb]);
			for (int i = 0; i < neighbNumGRBMs[iNeighb]; ++i)
				cornerEqNums[offset++] = fDof + i;
		}
	} else if (getFetiInfo().isEdgeAugmentationOn()) {
		int iEdgeN = 0;
		for (int iNeighb = 0; iNeighb < scomm->numNeighb; ++iNeighb) {
			if (scomm->isEdgeNeighb[iNeighb]) {
				int fDof = cornerEqs->firstdof(augOffset + (*subToEdge)[subNum()][iEdgeN]);
				if (isMixedSub) {
					for (int i = 0; i < edgeDofSize[iNeighb] - edgeDofSizeTmp[iNeighb]; ++i)  // fluid
						cornerEqNums[offset++] = fDof + i;
				} else {
					for (int i = 0; i < edgeDofSize[iNeighb]; ++i)
						cornerEqNums[offset++] = fDof + i;
				}
				iEdgeN++;
			}
		}
		iEdgeN = 0;
		if (isMixedSub) {
			for (int iNeighb = 0; iNeighb < scomm->numNeighb; ++iNeighb) {
				if (scomm->isEdgeNeighb[iNeighb]) {
					int fDof = cornerEqs->firstdof(augOffset + (*subToEdge)[subNum()][iEdgeN]) +
					           edgeDofSize[iNeighb] - edgeDofSizeTmp[iNeighb];
					for (int i = 0; i < edgeDofSizeTmp[iNeighb]; ++i)  // structure
						cornerEqNums[offset++] = fDof + i;
					iEdgeN++;
				}
			}
		}
	}
}

void
FetiBaseSub::setSComm(SComm *sc)
{
	scomm = sc;
}


void FetiBaseSub::setLocalMpcToBlock(const Connectivity *mpcToBlock, const Connectivity *blockToMpc)
{
	if(numMPC > 0) {
		int i, j;

		std::vector<size_t> ptr(numMPC+1); for(i=0; i<=numMPC; ++i) ptr[i] = i;
		std::vector<int> tgt(numMPC); for(i=0; i<numMPC; ++i) tgt[i] = localToGlobalMPC[i];
		localMpcToBlock = Connectivity(numMPC, std::move(ptr), std::move(tgt)).transcon(mpcToBlock);

		std::vector<int> target;
		target.reserve(localMpcToBlock->numConnect());
		std::vector<size_t> pointer(numMPC+1);
		for(i=0; i<numMPC; ++i) {
			pointer[i] = target.size();
			int gi = localToGlobalMPC[i];
			for(j=0; j<localMpcToBlock->num(i); ++j) {
				int jBlock = (*localMpcToBlock)[i][j];
				target.push_back(blockToMpc->cOffset(jBlock, gi));
			}
		}
		pointer[numMPC] = target.size();
		localMpcToBlockMpc = new Connectivity(numMPC, std::move(pointer), std::move(target));

		// HB & PJSA: for assembleBlockCCtsolver & extractBlockMpcResidual
		blockToLocalMpc = localMpcToBlock->alloc_reverse();

		if(localMpcToBlock->isDiagonal()) { // <=> a local lmpc belong to only ONE block
			blockToBlockMpc = std::make_unique<Connectivity>( blockToLocalMpc->transcon(*localMpcToBlockMpc) );
		} else {
			// HB & PJSA: to deal with possible overlapping (a local lmpc can belong to several blocks)
			blockToBlockMpc = std::make_unique<Connectivity>(*blockToLocalMpc);
			// over-write the target array
			auto &array = blockToBlockMpc->tgt();
			size_t count = 0;
			for(int iblk=0; iblk<blockToLocalMpc->csize(); ++iblk) {
				for(i=0; i<blockToLocalMpc->num(iblk); ++i) {
					int j = (*blockToLocalMpc)[iblk][i];
					int jb= localMpcToBlock->cOffset(j,iblk);
					array[count++] = (*localMpcToBlockMpc)[j][jb];
				}
			}
		}

		int size = totalInterfSize;
		std::vector<size_t> point(size + 1);
		std::vector<int> targ;
		targ.reserve(numMPC);
		for(i=0; i<size; ++i) {
			point[i] = targ.size();
			if(boundDofFlag[i] == 2) {
				int li = -1 - scomm->boundDofT(SComm::all,i);
				targ.push_back(li);
			}
		}
		point[size] = targ.size();
		Connectivity boundDofToMpc(size, std::move(point), std::move(targ));
		mpcToBoundDof = boundDofToMpc.alloc_reverse();
	}
	else {
		localMpcToBlock    = nullptr;
		localMpcToBlockMpc = nullptr;
		mpcToBoundDof      = nullptr;
		blockToLocalMpc    = nullptr;
		blockToBlockMpc.reset();
	}
}

int FetiBaseSub::getLocalMPCIndex(int globalMpcIndex) const {
	return globalToLocalMPC[globalMpcIndex];
}

int FetiBaseSub::getGlobalMPCIndex(int localMpcIndex) const {
	return localToGlobalMPC[localMpcIndex];
}

void FetiBaseSub::makeLocalMpcToGlobalMpc(const Connectivity *mpcToMpc)
{
	// PJSA: make a different localMpcToGlobalMpc that includes connections inside neighbors
	int i, j;
	int size = numMPC;
	// step 1: find size of target
	int numtarget = 0;
	for(i=0; i<size; i++) {
		int gi = localToGlobalMPC[i];
		for(j=0; j<mpcToMpc->num(gi); ++j) {
			int gj = (*mpcToMpc)[gi][j];
			if(globalToLocalMPC[gj] > -1) numtarget++;
		}
	}
	// step 2: fill target
	int *pointer = new int[size+1];
	int *target  = new int[numtarget];
	int count = 0;
	for(i=0; i<size; i++) {
		pointer[i] = count;
		int gi = localToGlobalMPC[i];
		for(j=0; j<mpcToMpc->num(gi); ++j) {
			int gj = (*mpcToMpc)[gi][j];
			int lj = globalToLocalMPC[gj];
			if(lj > -1) {
				target[count] = lj;
				count++;
			}
		}
	}
	pointer[i] = numtarget;

	// step 3: construct localMpcToGlobalMpc connectivity
	localMpcToGlobalMpc = new Connectivity(size, pointer, target);
}

void
FetiBaseSub::GramSchmidt(double *Q, bool *isUsed, DofSet desired, int nQPerNeighb, bool isPrimalAugmentation)
{
	double rtol = getFetiInfo().orthotol;
	double atol = getFetiInfo().orthotol2;
	if(rtol == 0.0) return;

	int i, j, k, l, m, iSub;
	Connectivity &sharedNodes = *(scomm->sharedNodes);
	int numNeighb = sharedNodes.csize();
	int numInterfNodes = sharedNodes.numConnect();

	int numDofPerNode;
	if((desired.contains(DofSet::Helm)) || (desired.contains(DofSet::Temp)))
		numDofPerNode = 1;
	else if((desired.contains(DofSet::XYZdisp)))
		numDofPerNode = (desired.contains(DofSet::XYZrot))?6:3;
	else {
		std::cerr << " *** WARNING: numDofPerNode = 0 in BaseSub::GramSchmidt(...) sub " << subNum() << std::endl;
		return;
	}

	int r_numdofs = cc_dsa->size();
	int *count = (int *) dbg_alloca(sizeof(int)*r_numdofs);
	for(i = 0; i < r_numdofs; ++i) count[i] = 0;
	for(iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
		for(i = 0; i < scomm->lenT(SComm::all,iSub); ++i) {
			int dof = scomm->boundDofT(SComm::all,iSub,i);
			if(dof > -1)
				if(count[dof] < 2) count[dof]++;
		}
	}

	// Apply modified Gram Schmidt twice to orthogonalize the column vectors of Q matrix
	int *dxyz = (int *) dbg_alloca(sizeof(int)*numDofPerNode);
	for(i = 0; i < 2; ++i) {
		for(iSub = 0; iSub < numNeighb; ++iSub) {
			int vlen = numDofPerNode*sharedNodes.num(iSub);
			double *tmpV = (double *) dbg_alloca(sizeof(double)*vlen);
			int sOffset = sharedNodes.offset(iSub);
			for(j = 0; j < nQPerNeighb; ++j) {
				if(isUsed[j+iSub*nQPerNeighb]) {
					StackVector vj(vlen, Q + numDofPerNode*(numInterfNodes*j + sOffset));
					StackVector vjt(vlen, tmpV);
					for(l = 0; l < sharedNodes.num(iSub); ++l) {
						cc_dsa->number(sharedNodes[iSub][l], desired, dxyz);
						for(m = 0; m < numDofPerNode; ++m)
							if((dxyz[m] >= 0) && (count[dxyz[m]] < 2))
								vjt[numDofPerNode*l+m] = vj[numDofPerNode*l+m];
							else {
								// Edge modification for primal augmentation 082213 JAT
								if(isPrimalAugmentation)
									vj[numDofPerNode*l+m] = 0.0;
								vjt[numDofPerNode*l+m] = 0.0;
							}
					}

					double initNorm = vjt.norm();
					for(k = 0; k < j; ++k) {
						if(!isUsed[k+iSub*nQPerNeighb])
							continue;
						StackVector vk(vlen, Q + numDofPerNode*(numInterfNodes*k + sOffset));
						StackVector vkt(vlen, tmpV);
						for(l = 0; l < sharedNodes.num(iSub); ++l) {
							cc_dsa->number(sharedNodes[iSub][l], desired, dxyz);
							for(m = 0; m < numDofPerNode; ++m)
								if((dxyz[m] >= 0) && (count[dxyz[m]] < 2))
									vkt[numDofPerNode*l+m] = vk[numDofPerNode*l+m];
								else
									vkt[numDofPerNode*l+m] = 0.0;
						}
						vj.linAdd(-(vkt*vj), vk);
					}

					for(l = 0; l < sharedNodes.num(iSub); ++l) {
						cc_dsa->number(sharedNodes[iSub][l], desired, dxyz);
						for(int m = 0; m < numDofPerNode; ++m)
							if((dxyz[m] >= 0) && (count[dxyz[m]] < 2))
								vjt[numDofPerNode*l+m] = vj[numDofPerNode*l+m];
							else
								vjt[numDofPerNode*l+m] = 0.0;
					}
					double newNorm = vjt.norm();
					// if((newNorm <= rtol*initNorm) || ((i == 1) && (newNorm <= 0.4*initNorm)) { 
					if((newNorm <= rtol*initNorm) || ((i == 1) && (newNorm <= 0.4*initNorm)) || (newNorm <= atol)) { // PJSA
						isUsed[j+iSub*nQPerNeighb] = false;
					}
					else
						vj *= 1.0/newNorm;
				}
			}
		}
	}
}

void
FetiBaseSub::addMPCsToGlobalZstar(FullM &globalZstar, int startRow, int startCol, int numCol) const
{
	FullM *Zmpc = rigidBodyModesG->Zmpc;
	for(int i=0; i<numMPC_primal; ++i)
		for(int j=0; j<numCol; ++j)
			globalZstar[startRow+localToGroupMPC[i]][startCol+j] += (*Zmpc)[i][j];
}

void
FetiBaseSub::addSPCsToGlobalZstar(FullM &globalZstar, int &zRow, int zColOffset) const
{
	FullM *Zstar = rigidBodyModesG->Zstar;
	globalZstar.add(*Zstar, zRow, zColOffset);
	zRow += Zstar->numRow();
}


void
FetiBaseSub::setWIoneCommSize(FSCommStructure *pat) const
{
	for(int i = 0; i < scomm->numT(SComm::fsi); ++i)
		if(subNum() != scomm->neighbT(SComm::fsi,i))
			pat->setLen(subNum(), scomm->neighbT(SComm::fsi,i), 1);
}

void
FetiBaseSub::setWICommSize(FSCommStructure *pat) {
	for (int i = 0; i < scomm->numT(SComm::fsi); ++i)
		pat->setLen(subNum(), scomm->neighbT(SComm::fsi, i), numNeighbWIdof[i]);
}

void
FetiBaseSub::setWImapCommSize(FSCommPattern<int> *pat)
{
	for(int i = 0; i < scomm->numT(SComm::fsi); ++i)
		if(subNum() != scomm->neighbT(SComm::fsi,i))
			glToLocalWImap.setCommSize(pat, subNum(), scomm->neighbT(SComm::fsi,i));
}


void FetiBaseSub::setNumGroupRBM(int *ngrbmGr)
{
	groupRBMoffset = 0;
	for(int i=0; i<group; ++i) groupRBMoffset += ngrbmGr[i];
	numGroupRBM = ngrbmGr[group];
}

void FetiBaseSub::getNumGroupRBM(int *ngrbmGr)
{
	//cerr << "in getNumGroupRBM " << subNum() << " " << group << " " << numGroupRBM << std::endl;
	ngrbmGr[group] = numGroupRBM;
}

void
FetiBaseSub::sendNumWIdof(FSCommPattern<int> *pat) const
{
	for(int i = 0; i < scomm->numT(SComm::fsi); ++i) {
		if(subNum() != scomm->neighbT(SComm::fsi,i)) {
			FSSubRecInfo<int> sInfo = pat->getSendBuffer(subNum(), scomm->neighbT(SComm::fsi,i));
			sInfo.data[0] = numWIdof;
		}
	}
}

void
FetiBaseSub::recvNumWIdof(FSCommPattern<int> *pat)
{
	numNeighbWIdof.resize(scomm->numT(SComm::fsi));
	for(int i = 0; i < scomm->numT(SComm::fsi); ++i) {
		if(subNum() != scomm->neighbT(SComm::fsi,i)) {
			FSSubRecInfo<int> rInfo = pat->recData(scomm->neighbT(SComm::fsi,i), subNum());
			numNeighbWIdof[i] = rInfo.data[0];
		}
		else numNeighbWIdof[i] = 0;
	}
}

void
FetiBaseSub::sendWImap(FSCommPattern<int> *pat)
{
	for(int i=0; i< scomm->numT(SComm::fsi); ++i)
		if(subNum() != scomm->neighbT(SComm::fsi,i))
			glToLocalWImap.pack(pat,subNum(), scomm->neighbT(SComm::fsi,i));
}

void
FetiBaseSub::recvWImap(FSCommPattern<int> *pat)
{
	if(neighbGlToLocalWImap) delete [] neighbGlToLocalWImap;
	neighbGlToLocalWImap = new GlobalToLocalMap[scomm->numT(SComm::fsi)];
	for(int i=0; i<scomm->numT(SComm::fsi); ++i)
		if(subNum() != scomm->neighbT(SComm::fsi,i))
			neighbGlToLocalWImap[i].unpack(pat, scomm->neighbT(SComm::fsi,i), subNum());
}


void
FetiBaseSub::sendNeighbGrbmInfo(FSCommPattern<int> *pat)
{
	// send number of group GRBMs and the group GRBM offset to each potential contact neighbor
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc, i);
		FSSubRecInfo<int> sInfo = pat->getSendBuffer(subNum(), neighb);
		sInfo.data[0] = numGroupRBM;
		sInfo.data[1] = groupRBMoffset;
	}
}

void
FetiBaseSub::receiveNeighbGrbmInfo(FSCommPattern<int> *pat)
{
	if(neighbNumGroupGrbm) delete [] neighbNumGroupGrbm;
	neighbNumGroupGrbm = new int[scomm->numT(SComm::mpc)];
	if(neighbGroupGrbmOffset) delete [] neighbGroupGrbmOffset;
	neighbGroupGrbmOffset = new int[scomm->numT(SComm::mpc)];
	// get number of group GRBMs and the group GRBM offset for each potential contact neighbor
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc, i);
		FSSubRecInfo<int> rInfo = pat->recData(neighb, subNum());
		neighbNumGroupGrbm[i] = rInfo.data[0];
		neighbGroupGrbmOffset[i] = rInfo.data[1];
	}
}


void
FetiBaseSub::sendNumNeighbGrbm(FSCommPattern<int> *pat)
{
	// send Number of RBMs for each neighbor, used for augmentation
	for(int i = 0; i < scomm->numT(SComm::std); ++i) {
		FSSubRecInfo<int> sInfo = pat->getSendBuffer(subNum(), scomm->neighbT(SComm::std,i));
		sInfo.data[0] = nGrbm;
	}
}

void
FetiBaseSub::recvNumNeighbGrbm(FSCommPattern<int> *pat)
{
	neighbNumGRBMs = new int[scomm->numT(SComm::std)];
	// get Number of RBMs for each neighbor, used for augmentation
	for(int i = 0; i < scomm->numT(SComm::std); ++i) {
		FSSubRecInfo<int> rInfo = pat->recData(scomm->neighbT(SComm::std,i), subNum());
		neighbNumGRBMs[i] = rInfo.data[0];
	}
}

void
FetiBaseSub::setNodeCommSize(FSCommStructure *pt, int d) const
{
	for(int iSub = 0; iSub < scomm->numNeighb; ++iSub)
		pt->setLen(subNumber, scomm->subNums[iSub], scomm->sharedNodes->num(iSub)*d);
}

void FetiBaseSub::makeLocalToGroupMPC(const Connectivity &groupToMPC)
{
	// PJSA: new version for multi-body mpc compatability
	int i;
	if(numMPC_primal > 0) {
		localToGroupMPC = new int[numMPC_primal];
		int groupOffset = groupToMPC.offset(group);
		for(i=0; i<groupToMPC.num(group); ++i) {
			int glMpcID = groupToMPC.getTargetValue(groupOffset+i);
			int localMpcID = globalToLocalMPC_primal[glMpcID];
			if(localMpcID > -1) localToGroupMPC[localMpcID] = i;
		}
	}
}

template <typename Scalar>
void FetiSub<Scalar>::makeBs() {
	std::vector<Eigen::Triplet<double>> b;
	std::vector<Eigen::Triplet<double>> bw;
	std::vector<Eigen::Triplet<Scalar>> bm;
	std::vector<Eigen::Triplet<Scalar>> bc;
	std::vector<bool> mpcFlag(numMPC, true);
	for (int iDof = 0; iDof < totalInterfSize; ++iDof) {
		switch (boundDofFlag[iDof]) {
			case 0: // note B is used with a - sign.
				b.emplace_back(allBoundDofs[iDof], iDof, 1.0);
				break;
			case 1:  // wet interface Note Bw is used with a + sign
				bw.emplace_back(-1 - allBoundDofs[iDof], iDof, 1.0);
				break;
			case 2: { // dual mpc or contact. Note Bm is used with a - sign.
				int locMpcNb = -1 - allBoundDofs[iDof];
				if (mpcFlag[locMpcNb]) {
					const auto &m = mpc[locMpcNb];
					for (int k = 0; k < m->nterms; k++) {
						int ccdof = (m->terms)[k].ccdof;
						bm.emplace_back(ccdof, iDof, (m->terms)[k].coef);
					}
					mpcFlag[locMpcNb] = false;
				}
			}
				break;
		}
	}
	B.resize(localLen(), totalInterfSize);
	B.setZero();
	B.setFromTriplets(b.begin(), b.end()); 
	Bw.resize(numWIdof, totalInterfSize);
	Bw.setZero();
	Bw.setFromTriplets(bw.begin(), bw.end());
	Bm.resize(localLen(), totalInterfSize);
	Bm.setZero();
	Bm.setFromTriplets(bm.begin(), bm.end());

	mpcFlag.assign(numMPC, true);
	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		if (mpcFlag[locMpcNb]) {
			const auto &m = mpc[locMpcNb];
			for (int k = 0; k < m->nterms; k++) {
				int dof = (m->terms)[k].dof;
				if ((dof >= 0) && (cornerMap[dof] >= 0))
					bc.emplace_back(cornerMap[dof], scomm->mapT(SComm::mpc, i), (m->terms)[k].coef);
			}
			mpcFlag[locMpcNb] = false;
		}
	}
	Bc.resize(Src->numCol(), totalInterfSize);
	Bc.setFromTriplets(bc.begin(), bc.end());
}

template<typename Scalar>
double FetiSub<Scalar>::getMpcError() const {
	double ret = 0;
	for(int i = 0; i < numMPC; ++i) {
		if(mpcMaster[i]) {
			if(mpc[i]->type == 0) {
				ret += ScalarTypes::sqNorm(mpc[i]->rhs);
			}
			else if(mpc[i]->type == 1 && ScalarTypes::lessThan(mpc[i]->rhs, 0.)) ret += ScalarTypes::sqNorm(mpc[i]->rhs);
		}
	}
	return ret;
}

template<class Scalar>
void
FetiSub<Scalar>::makeLocalMpcToDof() {
	auto &mpc = this->mpc;
	if (mpcToDof) delete mpcToDof;
	mpcToDof = 0;

	// step 1: make mpcToDof
	int size = numMPC;
	// step 1.1: find size of target: total number of coefficients involving a different dof
	int numtarget = 0;
	int i, j, jj;
	for (i = 0; i < size; i++) {
		for (j = 0; j < mpc[i]->nterms; j++) {
			int dofj = mpc[i]->terms[j].dof;
			for (jj = 0; jj < j; jj++) {
				int dofjj = mpc[i]->terms[jj].dof;
				if (dofj == dofjj) break;
			}
			if ((jj == j) && (dofj >= 0)) numtarget++;
		}
	}
	// step 1.2: fill target with coefficient dofs
	int *pointer = new int[size + 1];
	int *target = new int[numtarget];
	int count = 0;
	for (i = 0; i < size; i++) {
		pointer[i] = count;
		for (j = 0; j < mpc[i]->nterms; j++) {
			int dofj = mpc[i]->terms[j].dof;
			for (jj = 0; jj < j; jj++) {
				int dofjj = mpc[i]->terms[jj].dof;
				if (dofj == dofjj) break;
			}
			if ((jj == j) && (dofj >= 0)) {
				target[count] = dofj;
				count++;
			}
		}
	}
	pointer[i] = numtarget;
	// step 1.3: construct mpcToDof connectivity
	mpcToDof = new Connectivity(size, pointer, target);
}

template<class Scalar>
void
FetiSub<Scalar>::makeLocalMpcToMpc() {
	// step 1: make mpcToDof
	makeLocalMpcToDof();

	// step 2: make localMpcToMpc connectivity
	Connectivity *dofToMpc = mpcToDof->alloc_reverse();
	localMpcToMpc = mpcToDof->transcon(dofToMpc);
	delete dofToMpc;
}


template<class Scalar>
void
FetiSub<Scalar>::applyMpcSplitting()
{
	// adjust discrete masses, forces and mpcs using subdomain multiplicity
	// num = number of subdomains touching a dof
	int cdof, num;

	// mpcs (NOTE: optional kscaling is done later, hhs is not split)
	if(getFetiInfo().mpc_scaling == FetiInfo::tscaling) {
		for(int iMPC = 0; iMPC < numMPC; ++iMPC) { // dual mpcs
			if(mpc[iMPC]->type == 2) continue; // bmpc
			for(int i = 0; i < mpc[iMPC]->nterms; ++i) {
				if((cdof = mpc[iMPC]->terms[i].cdof) > -1 && (num = weightPlus[cdof]) > 1)
					mpc[iMPC]->terms[i].coef /= double(num);
			}
		}
	}
	// XXXX kscaling currently not supported for primal mpcs
	for(int iMPC = 0; iMPC < numMPC_primal; ++iMPC) { // primal mpcs
		for(int i = 0; i < mpc_primal[iMPC]->nterms; ++i) {
			if((cdof = mpc_primal[iMPC]->terms[i].cdof) > -1 && (num = weightPlus[cdof]) > 1)
				mpc_primal[iMPC]->terms[i].coef /= double(num);
		}
	}

}


template<class Scalar>
void
FetiSub<Scalar>::subtractMpcRhs(Scalar *interfvec)
{
	for(int i=0; i < scomm->lenT(SComm::mpc); ++i) {
		interfvec[scomm->mapT(SComm::mpc,i)] -= mpc[scomm->mpcNb(i)]->rhs;
	}
}

template<class Scalar>
void
FetiSub<Scalar>::assembleGtGsolver(GenSparseMatrix<Scalar> *GtGsolver)
{
	if(numGroupRBM == 0) return;
	bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool)*numMPC);
	for(int i = 0; i < numMPC; ++i) mpcFlag[i] = true;
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int numGroupRBM2 = neighbNumGroupGrbm[i];
		int groupRBMoffset2 = neighbGroupGrbmOffset[i];
		GenVector<Scalar> d(scomm->lenT(SComm::mpc,i));
		for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
			int locMpcNb = scomm->mpcNb(i,j);
			d[j] = (mpc[locMpcNb]->active) ? 0.0 : 1.0;
		}
		if((numGroupRBM2 > 0) && (subNum() != scomm->neighbT(SComm::mpc,i))) {
			GenFullM<Scalar> tmp2(numGroupRBM, numGroupRBM2);  // coupling term
			G[i]->transposeMultD(*neighbG[i], d, tmp2); // tmp2 = G^T * D * neighbG
			GtGsolver->add(tmp2, groupRBMoffset, groupRBMoffset2);
		}
		for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
			int locMpcNb = scomm->mpcNb(i,j);
			if(!mpcFlag[locMpcNb]) d[j] = 0.0; // prevents duplication for mpc shared between more than 2 subdomains
			else mpcFlag[locMpcNb] = false;
		}
		GenFullM<Scalar> tmp(numGroupRBM, numGroupRBM);
		G[i]->transposeMultD(*G[i], d, tmp); // tmp = G^T * D * G
		GtGsolver->add(tmp, groupRBMoffset, groupRBMoffset);
	}
}


template<class Scalar>
void
FetiSub<Scalar>::getLocalMpcForces(double *mpcLambda, DofSetArray *cornerEqs,
                                        int mpcOffset, GenVector<Scalar> &uc) {
// XXXX needs some work to map both dual and primal into single mpcLambda array
	if (numMPC > 0 && numMPC_primal > 0) std::cerr << "unsupported feature in FetiSub::getLocalMpcForces \n";
	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) { // dual mpcs
		int locMpcNb = scomm->mpcNb(i);
		if (localLambda) mpcLambda[locMpcNb] = localLambda[scomm->mapT(SComm::mpc, i)];
		else mpcLambda[locMpcNb] = 0;
	}
	for (int i = 0; i < numMPC_primal; ++i) {
		int glMpcNb = localToGlobalMPC_primal[i];
		int dof = cornerEqs->firstdof(mpcOffset + glMpcNb);
		mpcLambda[i] = -ScalarTypes::Real(uc[dof]);
	}
//	if (salinasFlag)
//		for (int i = 0; i < numMPC + numMPC_primal; ++i)
//			mpcLambda[i] = -mpcLambda[i];  // different sign convention
}


template<class Scalar>
void
FetiSub<Scalar>::useKrrNullspace()
{
	// EXPERMENTAL... use alaebraic null space of Krr with no corners
	int neq = this->Krr->neqs();
	int nzem = this->Krr->numRBM();
	Rstar.setNewSize(neq, nzem);
	if(nzem > 0) {
		std::vector<Scalar> rbmv(neq*nzem);
		this->Krr->getNullSpace(rbmv.data());

		// Copy rigid body modes (rbm) into Rstar
		for(int m=0; m<nzem; ++m) {
			for(int i=0; i<neq; ++i)
				Rstar[i][m] = rbmv[i+m*neq];
		}
	}
}

template<>
void
FetiSub<double>::makeLocalRstar(FullM **Qtranspose)
{
	FullM &R = rigidBodyModesG->R;
	FullM *Rc = rigidBodyModesG->Rc;
	if(numMPC_primal > 0) {
		FullM Qbody(Qtranspose[group]->transpose(), R.numCol(), bodyRBMoffset, Qtranspose[group]->numRow(), 0);
		Rstar = R * Qbody;
	}
	else {
		Rstar = R % *(Qtranspose[group]);
	}
}

template<>
void
FetiSub<DComplex>::makeLocalRstar(FullM **Qtranspose)
{
	std::cerr << "FetiSub<DComplex>::makeLocalRstar(FullM **Qtranspose) is not implemented\n";
}

template<class Scalar>
void
FetiSub<Scalar>::addRalpha(Scalar *u, GenVector<Scalar> &alpha) const
{
	int i, j;
	for(i=0; i<Rstar.numRow(); ++i)
		for(j=0; j<Rstar.numCol(); ++j)
			u[i] += Rstar[i][j] * alpha[groupRBMoffset + j];
}

template<class Scalar>
void
FetiSub<Scalar>::assembleE(GenVector<Scalar> &e, Scalar *f) const
{
	if(numGroupRBM > 0) {
		GenVector<Scalar> local_e(numGroupRBM, 0.0);
		GenVector<Scalar> fvec(f, Rstar.numRow());
		local_e = Rstar ^ fvec; // = Rtranspose * fvec
		e.add(local_e, groupRBMoffset);
	}
}

template<class Scalar>
void
FetiSub<Scalar>::makeG()
{
	// make G for each potential contact/mpc neighbour
	G.resize(scomm->numT(SComm::mpc));
	neighbG.resize(scomm->numT(SComm::mpc));
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {  // loop through all potential contact/mpc neighbours
		neighbG[i] = 0;
		G[i] = std::make_unique<GenFullM<Scalar>>(scomm->lenT(SComm::mpc,i), numGroupRBM);
		if(numGroupRBM > 0) {
			G[i]->zero();
			for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {  // loop through potential contact/mpc nodes
				int locMpcNb = scomm->mpcNb(i,j);
				const auto &m = mpc[locMpcNb];
				for(int k = 0; k < m->nterms; k++) {
					int cDof = (m->terms)[k].cdof;
					if(cDof > -1) {
						for(int iRbm = 0; iRbm < numGroupRBM; ++iRbm)
							(*G[i])[j][iRbm] += Rstar[cDof][iRbm]*(m->terms)[k].coef;
					}
				}
			}
		}
	}
}


template<class Scalar>
void
FetiSub<Scalar>::makeTrbmG(Scalar *rbms, int nrbm, int size)
{
	auto &Src = this->Src;
	auto &mpc = this->mpc;
	// rbms is the null space of the global Kcc^* matrix
	// nrbm is the nullity of the global Kcc^* matrix
	// size is the number of rows and columns of the global Kcc^* matrix
	// TODO what about augmentation
	int numCDofs = (Src) ? Src->numCol() : 0;
	Scalar *localc = (Scalar *) dbg_alloca(sizeof(Scalar)*numCDofs);

	int nrbms_local = 0; int first = 0;
	std::map<int, int> localToGlobalRBM;
	std::map<int, int> globalToLocalRBM;
	for(int iRbm = 0; iRbm < nrbm; ++iRbm) {
		Scalar *rbm = rbms + size*iRbm;
		Scalar dot = 0.0;
		for(int i=0; i<numCDofs; ++i) {
			if(cornerEqNums[i] > -1) { dot += rbm[cornerEqNums[i]]*rbm[cornerEqNums[i]]; }
		}
		if(dot != 0.0) {
			localToGlobalRBM[nrbms_local] = iRbm;
			globalToLocalRBM[iRbm] = nrbms_local;
			if(nrbms_local == 0) first = iRbm;
			nrbms_local++;
		}
	}

	numGroupRBM = nrbms_local; // TODO this isn't general since global trbms may not be grouped like grbms
	groupRBMoffset = first;    // TODO this isn't general

	G.resize(scomm->numT(SComm::mpc));
	neighbG.resize(scomm->numT(SComm::mpc));
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		neighbG[i] = 0;
		G[i] = std::make_unique<GenFullM<Scalar>>(scomm->lenT(SComm::mpc,i), nrbms_local);
		G[i]->zero();
	}

	if(nrbms_local == 0 || numMPC == 0) return;

	Scalar *localr = new Scalar[localLen()];

	for(int iRbm = 0; iRbm < nrbms_local; ++iRbm) {
		int glRbmId = localToGlobalRBM[iRbm];
		Scalar *rbm = rbms + size*glRbmId;
		for(int i=0; i<numCDofs; ++i) {
			if(cornerEqNums[i] > -1) { localc[i] = rbm[cornerEqNums[i]]; }
			else localc[i] = 0.0;
		}

		// G = (-Br^(s)Krr^{-1}Krc + Bc)Lcc Nc
		for(int i=0; i<localLen(); ++i) localr[i] = 0.0;
		if(Src) Src->transposeMultSubtract(localc, localr);
		if(this->Krr) this->Krr->reSolve(localr);
		for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
			for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
				int locMpcNb = scomm->mpcNb(i,j);
				const auto &m = mpc[locMpcNb];
				for(int k = 0; k < m->nterms; k++) {
					int cc_dof = (m->terms)[k].ccdof;
					if(cc_dof >= 0) (*G[i])[j][iRbm] += localr[cc_dof]*(m->terms)[k].coef;
					else {
						int dof = (m->terms)[k].dof;
						if((dof >= 0) && (cornerMap[dof] >= 0))
							(*G[i])[j][iRbm] += localc[cornerMap[dof]] * (m->terms)[k].coef;
					}
				}
			}
		}
	}
	delete [] localr;
}

// ****************************************************************************************************

template<class Scalar>
void
FetiSub<Scalar>::setGCommSize(FSCommStructure *pat) const
{
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int nRow = G[i]->numRow();
		int nCol = numGroupRBM;
		pat->setLen(subNum(), scomm->neighbT(SComm::mpc, i), nRow*nCol);
	}
}

template<class Scalar>
void
FetiSub<Scalar>::sendG(FSCommPattern<Scalar> *rbmPat)
{
	if(numGroupRBM == 0) return;
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i)
		rbmPat->sendData(subNum(), scomm->neighbT(SComm::mpc, i),  G[i]->data());
}

template<class Scalar>
void
FetiSub<Scalar>::receiveG(FSCommPattern<Scalar> *rbmPat)
{
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		FSSubRecInfo<Scalar> rInfo = rbmPat->recData(scomm->neighbT(SComm::mpc, i), subNum());
		int nRow = G[i]->numRow();  // number of potential contact dofs on interface with neighb
		int nCol = neighbNumGroupGrbm[i];  //number of rbms for neighb's group
		neighbG[i] = std::make_unique<GenFullM<Scalar>>(rInfo.data.data(), nRow, nCol);
	}
}

template<class Scalar>
void FetiSub<Scalar>::zeroG()
{
	if(G.size() != 0) {
		for(int i=0; i<scomm->numT(SComm::mpc); ++i)
			if(G[i]) G[i]->zero();
	}
	if(neighbG.size() != 0) {
		for(int i=0; i<scomm->numT(SComm::mpc); ++i)
			if(neighbG[i]) neighbG[i]->zero();
	}
}

template<class Scalar>
void FetiSub<Scalar>::deleteG()
{
	G.clear();
	neighbG.clear();
	if(neighbNumGroupGrbm) { delete [] neighbNumGroupGrbm; neighbNumGroupGrbm = 0; }
	if(neighbGroupGrbmOffset) { delete [] neighbGroupGrbmOffset; neighbGroupGrbmOffset = 0; }
}

template<class Scalar>
void
FetiSub<Scalar>::multG(const GenVector<Scalar> &x, Scalar *y, Scalar alpha) const
{
	// y += alpha * G * x
	Scalar *mpcvec = new Scalar[numMPC];
	for(int i = 0; i < numMPC; ++i) mpcvec[i] = 0.0;
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc, i);
		for(int j = 0; j < scomm->lenT(SComm::mpc, i); ++j) {
			int locMpcNb = scomm->mpcNb(i,j);
			if(mpcvec[locMpcNb] == 0.0)
				for(int k = 0; k < numGroupRBM; ++k)
					mpcvec[locMpcNb] += (*G[i])[j][k] * x[k + groupRBMoffset];
			if(subNum() != neighb)
				for(int k = 0; k < neighbNumGroupGrbm[i]; ++k)
					mpcvec[locMpcNb] += (*neighbG[i])[j][k] * x[k + neighbGroupGrbmOffset[i]];
		}
	}
	for(int i = 0; i < scomm->lenT(SComm::mpc); ++i)
		y[scomm->mapT(SComm::mpc,i)] += alpha*mpcvec[scomm->boundDofT(SComm::mpc,i)];
	delete [] mpcvec;
}

template<class Scalar>
void
FetiSub<Scalar>::trMultG(const Scalar *x, GenVector<Scalar> &y, Scalar alpha) const
{
	// compute y += alpha * G^t * x
	bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool)*numMPC);
	for(int i = 0; i < numMPC; ++i) mpcFlag[i] = true;
	for(int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		for(int j = 0; j < scomm->lenT(SComm::mpc,i); ++j) {
			int locMpcNb = scomm->mpcNb(i,j);
			int iDof = scomm->mapT(SComm::mpc,i,j);
			if(mpcFlag[locMpcNb]) {
				for(int k = 0; k < numGroupRBM; ++k)
					y[k + groupRBMoffset] += alpha * (*G[i])[j][k] * x[iDof];
				mpcFlag[locMpcNb] = false;
			}
		}
	}
}


template<class Scalar>
void
FetiSub<Scalar>::multMFi(GenSolver<Scalar> *s, Scalar *u, Scalar *Fiu, int nRHS) const {
	// multMFi is never called in DP. Otherwise it will crash in contact
	int iRHS, iDof;
	int ndof = localLen();

	Scalar **iDisp = new Scalar *[nRHS];
	Scalar *firstpointer = new Scalar[nRHS * ndof];

	for (iRHS = 0; iRHS < nRHS; iRHS++) {
		iDisp[iRHS] = firstpointer + iRHS * ndof;
		for (iDof = 0; iDof < ndof; ++iDof)
			iDisp[iRHS][iDof] = 0.0;
	}

	// Add the interfvec contribution to localvec
	for (iRHS = 0; iRHS < nRHS; iRHS++) {
		for (iDof = 0; iDof < totalInterfSize; ++iDof)
			iDisp[iRHS][allBoundDofs[iDof]] += u[iDof + iRHS * totalInterfSize];
	}

	// solve for localvec
	if (s) s->reSolve(nRHS, iDisp);

	// redistribute the solution to the interface
	for (iRHS = 0; iRHS < nRHS; iRHS++) {
		for (iDof = 0; iDof < totalInterfSize; ++iDof)
			Fiu[iDof + iRHS * totalInterfSize] = iDisp[iRHS][allBoundDofs[iDof]];
	}

	delete[] firstpointer;
	delete[] iDisp;
}

template<class Scalar>
void
FetiSub<Scalar>::getQtKQ(GenSolver<Scalar> *s) {
	const auto &c_dsa = get_c_dsa();
	if (numMPC == 0) return;

	int numDOFs = localLen();
	std::vector<Scalar> locKpQ(numMPC *numDOFs,
	0.0);

	// loop over mpc structure and fill coefficients
	for (int iMPC = 0; iMPC < numMPC; ++iMPC) {
		for (int i = 0; i < mpc[iMPC]->nterms; ++i) {
			int dof = c_dsa->locate(mpc[iMPC]->terms[i].nnum,
			                        (1 << mpc[iMPC]->terms[i].dofnum));
			if (dof >= 0) {
				locKpQ[dof + iMPC * numDOFs] = mpc[iMPC]->terms[i].coef;
			}
		}
	}

	for (int iMPC = 0; iMPC < numMPC; ++iMPC)
		if (s) s->reSolve(locKpQ.data() + iMPC * numDOFs);

	QtKpBt.resize(numMPC * totalInterfSize);

	for (int i = 0; i < numMPC * totalInterfSize; ++i)
		QtKpBt[i] = 0.0;

	for (int i = 0; i < numMPC; ++i)
		for (int j = 0; j < totalInterfSize; ++j)
			QtKpBt[j + i * totalInterfSize] = locKpQ[allBoundDofs[j] + i * numDOFs];

	qtkq = std::make_unique<GenFullM<Scalar>>(numMPC, numMPC);
	qtkq->zero();

	for (int iMPC = 0; iMPC < numMPC; ++iMPC)
		for (int jMPC = 0; jMPC < numMPC; ++jMPC)
			for (int i = 0; i < mpc[iMPC]->nterms; ++i) {
				int dof = c_dsa->locate(mpc[iMPC]->terms[i].nnum,
				                        (1 << mpc[iMPC]->terms[i].dofnum));
				if (dof < 0) continue;

				(*qtkq)[iMPC][jMPC] += mpc[iMPC]->terms[i].coef * locKpQ[dof + jMPC * numDOFs];
			}

}

template<class Scalar>
void
FetiSub<Scalar>::getQtKQ(int glMPCnum, Scalar *QtKQ) {
	int thisMPC = globalToLocalMPC[glMPCnum];
	int iMPC;
	for (iMPC = 0; iMPC < numMPC; ++iMPC)
		QtKQ[iMPC] = (*qtkq)[thisMPC][iMPC];
}

template<class Scalar>
void
FetiSub<Scalar>::multQtKBt(int glMPCnum, const Scalar *G, Scalar *QtKBtG,
                                Scalar alpha, Scalar beta) const {
	int iMPC = globalToLocalMPC[glMPCnum];
	// QtKBtG = QtKBt*G
	Tgemv('T', totalInterfSize, 1, alpha, QtKpBt.data() + iMPC * totalInterfSize,
	      totalInterfSize, G, 1, beta, QtKBtG, 1);
}

template<class Scalar>
void
FetiSub<Scalar>::multQt(int glMPCnum, const Scalar *V, int numV, Scalar *QtV) const {
	auto *c_dsa = get_c_dsa();
	int numDofs = localLen();
	int iMPC = globalToLocalMPC[glMPCnum];
	for (int n = 0; n < numV; ++n) {
		for (int i = 0; i < mpc[iMPC]->nterms; ++i) {
			int dof = c_dsa->locate(mpc[iMPC]->terms[i].nnum,
			                        (1 << mpc[iMPC]->terms[i].dofnum));
			if (dof < 0) continue;
			const Scalar *beta = V + n * numDofs;
			QtV[n] += mpc[iMPC]->terms[i].coef * beta[dof];
		}
	}
}
// ************************************************************************************************


template<class Scalar>
void
FetiSub<Scalar>::split(const Scalar *v, Scalar *v_f, Scalar *v_c) const {
	// split v into free (v_f) and chopped (v_c) components
	for (int i = 0; i < totalInterfSize; ++i) {
		v_f[i] = v[i];
		v_c[i] = 0.0;
	}
	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		if (mpc[locMpcNb]->type == 1 && mpc[locMpcNb]->active) {
			int iDof = scomm->mapT(SComm::mpc, i);
			if (ScalarTypes::greaterThan(v[iDof], 0.0)) v_c[iDof] = v[iDof];
			v_f[iDof] = 0.0;
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::sendInterf(const Scalar *interfvec, FSCommPattern<Scalar> *vPat) const {
	int iDof = 0;
	for (int i = 0; i < scomm->numT(SComm::all); ++i) {
		FSSubRecInfo<Scalar> sInfo = vPat->getSendBuffer(subNum(), scomm->neighbT(SComm::all, i));
		for (int j = 0; j < scomm->lenT(SComm::all, i); ++j) {
			sInfo.data[j] = interfvec[iDof++];
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::buildGlobalRBMs(GenFullM<Scalar> &Xmatrix, const Connectivity *cornerToSub)
{
	int i,j,k;
	if(numGroupRBM == 0) numGlobalRBMs = 0;
	else {
		numGlobalRBMs = Xmatrix.numCol();
		GenFullM<Scalar> groupX(Xmatrix, numGroupRBM, groupRBMoffset, numGlobalRBMs, 0);
		Rstar_g = Rstar * groupX;

		if(!cornerToSub) return;
		std::vector<double> sharedUse(Rstar.numRow(), 1.0);
		// if i is a shared dof set sharedUse[i] = 0.0
		// for all but one of the subdomains sharing it in this body

		for(i=0; i<scomm->numT(SComm::std); ++i) {  // check non-corner dofs
			if(subNum() > scomm->neighbT(SComm::std,i))
				for(j=0; j<scomm->lenT(SComm::std,i); ++j)
					sharedUse[ccToC[scomm->boundDofT(SComm::std,i,j)]] = 0.0;
		}
		for(i=0; i<numCRN; ++i) { // check corner dofs
			if(subNum() != (*cornerToSub)[glCornerNodes[i]][0]) {
				int lDof[6];
				get_c_dsa()->number(cornerNodes[i], cornerDofs[i].list(), lDof);
				for(j=0; j<cornerDofs[i].count(); ++j)
					if(lDof[j] >= 0) sharedUse[lDof[j]] = 0.0;
			}
		}

		// if i is a shared dof sharedRstar_g[i] is set to zero 
		// for all but one of the subdomains sharing it in this body
		// (used to prevent duplication in construction of RtR)
		sharedRstar_g = std::make_unique<GenFullM<Scalar>>(Rstar_g.numRow(), Rstar_g.numCol());
		for(i=0; i<Rstar_g.numRow(); ++i)
			for(j=0; j<Rstar_g.numCol(); ++j) (*sharedRstar_g)[i][j] = Rstar_g[i][j] * sharedUse[i];

		// if i is a shared dof tmpRstar_g[i] is set to Rstar_g[i]/n 
		// where n is the number of subdomains (is this body) sharing this dof
		// used to compute a distributed vector, since distvec[i]*n = actual value at dof i  
		tmpRstar_g = std::make_unique<GenFullM<Scalar>>(Rstar_g);
		for(i=0; i<scomm->numT(SComm::all); ++i) {
			for(j=0; j<scomm->lenT(SComm::all,i); ++j) {
				int bdof = scomm->boundDofT(SComm::all,i,j);
				if(bdof >= 0)  // not a contact dof 
					for(k=0; k<numGlobalRBMs; ++k)
						(*tmpRstar_g)[ccToC[bdof]][k] /= getWeights()[bdof];
			}
		}
		for(i=0; i<numCRN; ++i) {
			int lDof[6];
			get_c_dsa()->number(cornerNodes[i], cornerDofs[i].list(), lDof);
			for(j=0; j<cornerDofs[i].count(); ++j)
				if(lDof[j] >= 0)
					for(k=0; k<numGlobalRBMs; ++k)
						(*tmpRstar_g)[lDof[j]][k] /= (double) (cornerToSub->num(glCornerNodes[i]));
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::getGlobalRBM(int iRBM, Scalar *Rvec) const
{
	if(numGlobalRBMs > 0)
		for(int iRow=0; iRow<Rstar_g.numRow(); ++iRow)
			Rvec[iRow] = Rstar_g[iRow][iRBM];
}

template<class Scalar>
void
FetiSub<Scalar>::subtractRstar_g(Scalar *u, GenVector<Scalar> &beta) const
{
	int i;
	if(numGlobalRBMs > 0) {
		// compute u = u - Rstar_g * beta  (second part of displacement projection)
		GenVector<Scalar> tmpu(Rstar_g.numRow());
		tmpu = Rstar_g * beta;
		for(i=0; i<Rstar_g.numRow(); ++i)
			u[i] = u[i] - tmpu[i];
	}
}

template<class Scalar>
void
FetiSub<Scalar>::addRstar_gT(Scalar *u, GenVector<Scalar> &beta) const
{
	if(numGlobalRBMs > 0) {
		// compute beta += Rstar_g^t * u  (first part of displacement projection)
		GenVector<Scalar> uvec(u, Rstar_g.numRow());
		beta += *tmpRstar_g ^ uvec;
	}
}

template<class Scalar>
void
FetiSub<Scalar>::assembleRtR(GenFullM<Scalar> &RtR)
{
	// builds RtR for u projection
	if(numGlobalRBMs > 0) {
		GenFullM<Scalar> tmp(numGlobalRBMs, numGlobalRBMs);
		sharedRstar_g->transposeMult((*sharedRstar_g), tmp);
		RtR.add(tmp, 0, 0);
	}
}


template<class Scalar>
void
FetiSub<Scalar>::multAddCT(const Scalar *interfvec, Scalar *localvec) const {
	auto &mpc = this->mpc;
	// localvec += C^T * interfvec
	bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool) * numMPC);
	for (int i = 0; i < numMPC; ++i) mpcFlag[i] = true;

	for (int i = 0; i < scomm->lenT(SComm::mpc); i++) {
		int locMpcNb = scomm->mpcNb(i);
		if (mpcFlag[locMpcNb]) {
			const auto &m = mpc[locMpcNb];
			int iDof = scomm->mapT(SComm::mpc, i);
			for (int k = 0; k < m->nterms; k++) {
				int c_dof = (m->terms)[k].cdof;
				if (c_dof >= 0) localvec[c_dof] += interfvec[iDof] * (m->terms)[k].coef;
			}
			mpcFlag[locMpcNb] = false;
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::multC(const Scalar *localvec, Scalar *interfvec) const {
	auto &mpc = this->mpc;
	// interfvec = C * localvec
	for (int i = 0; i < scomm->lenT(SComm::mpc); i++) {
		int locMpcNb = scomm->mpcNb(i);
		const auto &m = mpc[locMpcNb];
		int iDof = scomm->mapT(SComm::mpc, i);
		interfvec[iDof] = 0;
		for (int k = 0; k < m->nterms; k++) {
			int c_dof = (m->terms)[k].cdof;
			if (c_dof >= 0) interfvec[iDof] += localvec[c_dof] * (m->terms)[k].coef;
		}
	}
}


template<class Scalar>
void
FetiSub<Scalar>::addTrbmRalpha(Scalar *rbms, int nrbms, int glNumCDofs, Scalar *alpha, Scalar *ur) const
{
	auto &Src = this->Src;
	int numCDofs = (Src) ? Src->numCol() : 0;
	Scalar *localc = (Scalar *) dbg_alloca(sizeof(Scalar)*numCDofs);
	Scalar *localr = new Scalar[localLen()];

	for(int iRbm = 0; iRbm < nrbms; ++iRbm) {
		Scalar *rbm = rbms + glNumCDofs*iRbm;
		for(int i=0; i<numCDofs; ++i) {
			if(cornerEqNums[i] > -1) { localc[i] = rbm[cornerEqNums[i]]*alpha[iRbm]; }
			else localc[i] = 0.0;
		}

		for(int i = 0; i < localLen(); ++i) localr[i] = 0.0;
		if(Src) Src->transposeMultAdd(localc, localr);
		if(this->Krr) this->Krr->reSolve(localr);

		for(int i = 0; i < localLen(); ++i) ur[i] -= localr[i];
	}

	delete [] localr;
}

template<class Scalar>
void
FetiSub<Scalar>::assembleTrbmE(Scalar *rbms, int nrbms, int size, Scalar *e, Scalar *fr) const
{
	auto &Src = this->Src;
	int numCDofs = (Src) ? Src->numCol() : 0;
	Scalar *localc = (Scalar *) dbg_alloca(sizeof(Scalar)*numCDofs);
	for(int i=0; i<numCDofs; ++i) localc[i] = 0.0;

	Scalar *localr = new Scalar[localLen()];
	for(int i = 0; i < localLen(); ++i) localr[i] = -fr[i];
	if(this->Krr) this->Krr->reSolve(localr);
	if(Src) Src->multAdd(localr, localc); // localc = - (Krr^-1 Krc)^T fr
	delete [] localr;

	for(int iRbm = 0; iRbm < nrbms; ++iRbm) {
		Scalar *rbm = rbms + size*iRbm;
		for(int i=0; i<numCDofs; ++i) {
			if(cornerEqNums[i] > -1) { e[iRbm] += rbm[cornerEqNums[i]]*localc[i]; } // e += -N^T (Krr^-1 Krc)^T fr
		}
	}
}

// ****************************************************************************************************

template<class Scalar>
Scalar
FetiSub<Scalar>::getMpcRhs(int glMPCnum) const {
	return mpc[globalToLocalMPC[glMPCnum]]->rhs;
}

template<class Scalar>
Scalar
FetiSub<Scalar>::getMpcRhs_primal(int glMPCnum) const {
	return mpc_primal[globalToLocalMPC_primal[glMPCnum]]->rhs;
}

/**
* @param[out] fr Force on remainder DOFs.
* @param[in] uc  Corner displacements.
*/
template<class Scalar>
void
FetiSub<Scalar>::multKrc(Scalar *fr, const Scalar *uc) const {
	const auto &Src = this->Src;
	int numCDofs = Src->numCol();
	Scalar *ucLocal = (Scalar *) dbg_alloca(sizeof(Scalar) * numCDofs);

	for (int i = 0; i < numCDofs; ++i) {
		if (cornerEqNums[i] > -1) ucLocal[i] = uc[cornerEqNums[i]];
		else ucLocal[i] = 0.0;
	}

	// Perform fr = fr - Krc uc
	if (Src) Src->transposeMultSubtract(ucLocal, fr);
}


template<class Scalar>
void
FetiSub<Scalar>::multKbbMpc(const Scalar *u, Scalar *Pu, Scalar *deltaU, Scalar *deltaF, bool errorFlag) {
	// KHP and DJR: 3-26-98
	// multKbb has been modified to compute subdomain primal residual using
	// deltaF (in addition to deltaU which is the displacement correction
	// and of course the lumped or dirichlet preconditioning)

	// If we are computing lumped preconditioner, we compute deltaFi using Kib
	// but deltaUi is equal to zero. For the dirichlet preconditioner, deltaUi
	// is computed but deltaFi is set equal to zero.

	// deltaFi = internal primal residual
	// deltaUi = internal displacement correction
	// boundMap = from boundary number to subdomain number
	// internalMap = from internal number to subdomain number
	// invBoundMap = from all subdomain dof to a unique boundary number
	// invInternalMap =  from all subdomain dof to a unique internal number
	// allBoundDofs = indices of B, from lambda numbering to numbering of entire subdomain
	// dualToBoundary = from lambda numbering directly to boundary numbering

	// Kii works only on the numbering of the internal dofs
	// Kbb works only on the numbering of the boundary dofs
	// Kib operates on the boundary numbering and returns with internal numbering

	auto &localw = this->localw;
	auto &Kib = this->Kib;
	auto &KiiSolver = this->KiiSolver;

	Scalar *v = (Scalar *) dbg_alloca(sizeof(Scalar) * boundLen);
	Scalar *res = (Scalar *) dbg_alloca(sizeof(Scalar) * boundLen);
	this->deltaFwi.resize(numWIdof); // coupled_dph

	int i, iDof;
	for (iDof = 0; iDof < boundLen; ++iDof) v[iDof] = res[iDof] = 0.0;
	for (iDof = 0; iDof < localLen(); ++iDof) {
		if (deltaU) deltaU[iDof] = 0.0;
		if (deltaF) deltaF[iDof] = 0.0;
	}
	for (i = 0; i < numWIdof; ++i) localw[i] = 0;

	this->applyBtransposeAndScaling(u, v, deltaU, localw.data());

	//Scalar norm = 0; for(i=0; i<boundLen; ++i) norm += v[i]*v[i]; std::cerr << "1. norm = " << sqrt(norm) << std::endl;
	if ((getFetiInfo().precno == FetiInfo::lumped) ||
	    (getFetiInfo().precno == FetiInfo::dirichlet) || errorFlag)
		this->Kbb->mult(v, res);  // res = this->Kbb * v
	//norm = 0; for(i=0; i<boundLen; ++i) norm += res[i]*res[i]; std::cerr << "2. norm = " << sqrt(norm) << std::endl;

	if ((getFetiInfo().precno == FetiInfo::dirichlet) || errorFlag) {
		Scalar *iDisp = new Scalar[internalLen];
		for (iDof = 0; iDof < internalLen; ++iDof) iDisp[iDof] = 0.0;
		for (i = 0; i < numWIdof; ++i) iDisp[wiInternalMap[i]] = localw[i]; // coupled_dph
		if (Kib) Kib->transposeMultAdd(v, iDisp); // iDisp += Kib^T * v

		if (getFetiInfo().precno == FetiInfo::dirichlet) {
			//norm = 0; for(i=0; i<internalLen; ++i) norm += iDisp[i]*iDisp[i]; std::cerr << "3. norm = " << sqrt(norm) << std::endl;
			if (KiiSolver) KiiSolver->reSolve(iDisp);
			//norm = 0; for(i=0; i<internalLen; ++i) norm += iDisp[i]*iDisp[i]; std::cerr << "4. norm = " << sqrt(norm) << std::endl;
			for (i = 0; i < numWIdof; ++i) localw[i] = iDisp[wiInternalMap[i]]; // coupled_dph
			if (Kib) Kib->multSub(iDisp, res); // res -= Kib*iDisp
		} else if (deltaF) { // improves estimate of error
			for (iDof = 0; iDof < internalLen; ++iDof) {
				if (internalMap[iDof] > -1) {
					int ccDof = cToCC[internalMap[iDof]];
					if (ccDof > -1) deltaF[ccDof] = iDisp[iDof];
				}
			}
		}
		delete[] iDisp;

		if (deltaF) {
			for (iDof = 0; iDof < totalInterfSize; ++iDof)
				if (allBoundDofs[iDof] >= 0) { // deltaF of ctc and mpc nodes computed in applyScalingAndB below
					deltaF[allBoundDofs[iDof]] = res[dualToBoundary[iDof]];
				}
		}
	}

	if (getFetiInfo().precno == FetiInfo::identity) {
		for (iDof = 0; iDof < boundLen; ++iDof) res[iDof] = v[iDof];
	}

	// Return preconditioned u
	this->applyScalingAndB(res, Pu, localw.data());
}

template<class Scalar>
void
FetiSub<Scalar>::getFr(const Scalar *f, Scalar *fr) const {
	Scalar v[Ave.cols()];
	VectorView<Scalar> t(v, Ave.cols(), 1);
	VectorView<Scalar> l(fr, Ave.rows(), 1);


	for(int dof = 0; dof < ccToC.size(); ++dof)
		fr[dof] = f[ccToC[dof]];

	if (Ave.cols() > 0) { // Averages to zero 072513 JAT
		t = Ave.transpose() * l;
		l.noalias() -= Ave * t;
	}
}

template<class Scalar>
void
FetiSub<Scalar>::getFc(const Scalar *f, Scalar *Fc) const {
	int i, j;
	int dNum[DofSet::max_known_dof];
	int iOff = 0;
	const auto &c_dsa = get_c_dsa();
	for (i = 0; i < numCRN; ++i) {
		int nd = c_dsa->number(cornerNodes[i], cornerDofs[i], dNum);
		for (j = 0; j < nd; ++j) Fc[iOff + j] = f[dNum[j]];
		iOff += nd;
	}

	if (Ave.cols() > 0) { // Average corners 072513 JAT
		int numEquations = this->Krr->neqs();
		Scalar fr[numEquations];
		for(int dof = 0; dof < ccToC.size(); ++dof)
			fr[dof] = f[ccToC[dof]];

		int nAve, nCor, k;
		Scalar s;
		nCor = this->Krc ? this->Krc->numCol() : 0;
		nAve = Src->numCol() - nCor;
		for (int i = 0; i < nAve; ++i) {
			s = 0.0;
			for (int k = 0; k < numEquations; ++k)
				s += Ave[i][k] * fr[k];
			Fc[nCor + i] = s;
		}
	}
}


template<class Scalar>
void
FetiSub<Scalar>::projectActiveIneq(Scalar *v) const {
	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		if (mpc[locMpcNb]->type == 1 && mpc[locMpcNb]->active)
			v[scomm->mapT(SComm::mpc, i)] = 0.0;
	}
}


// ref: Dostal, Horak and Stefanica (IMACS 2005) "A scalable FETI-DP algorithm for coercive variational inequalities"
// every row of C matrix should have a unit norm to improve the condition number of the Feti operator
template<class Scalar>
void
FetiSub<Scalar>::normalizeCstep1(Scalar *cnorm) {
	auto &mpc = this->mpc;
	for (int i = 0; i < numMPC; ++i)
		for (int j = 0; j < mpc[i]->nterms; ++j)
			cnorm[localToGlobalMPC[i]] += mpc[i]->terms[j].coef * mpc[i]->terms[j].coef;
}

template<class Scalar>
void
FetiSub<Scalar>::normalizeCstep2(Scalar *cnorm) {
	auto &mpc = this->mpc;
	for (int i = 0; i < numMPC; ++i)
		for (int j = 0; j < mpc[i]->nterms; ++j)
			mpc[i]->terms[j].coef /= cnorm[localToGlobalMPC[i]];
}

template<class Scalar>
void
FetiSub<Scalar>::getFw(const Scalar *f, Scalar *fw) const {
	auto *dsa = getDsa();
	auto *c_dsa = get_c_dsa();
	if (numWIdof) {
		int i, j;
		int dofs[DofSet::max_known_dof];
		int cdofs[DofSet::max_known_dof];
		for (i = 0; i < numWInodes; ++i) {
			DofSet thisDofSet = wetInterfaceDofs[i];
			int nd = thisDofSet.count();
			dsa->number(wetInterfaceNodes[i], thisDofSet, dofs);
			c_dsa->number(wetInterfaceNodes[i], thisDofSet, cdofs);
			for (j = 0; j < nd; ++j)
				fw[wetInterfaceMap[dofs[j]]] = f[cdofs[j]];
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::recvMpcStatus(FSCommPattern<int> *mpcPat, int flag, bool &statusChange) {
	auto &mpc = this->mpc;
	// this function is to make sure that the status of an mpc is the same in all subdomains which share it
	// needed to due to possible roundoff error
	// if flag == 1 then make dual constraint not active in all subdomains if not active in at least one (use in proportioning step)
	// if flag == 0 then make dual constraint active in all subdomains if it is active in at least one (use in expansion step)
	// if flag == -1 use mpc master status in all subdomains
	// note: could use SComm::ieq list
	int i, j;
	bool *tmpStatus = (bool *) alloca(sizeof(bool) * numMPC);
	for (int i = 0; i < numMPC; ++i) tmpStatus[i] = !mpc[i]->active;
	for (i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc, i);
		if (subNum() != neighb) {
			FSSubRecInfo<int> rInfo = mpcPat->recData(neighb, subNum());
			for (j = 0; j < scomm->lenT(SComm::mpc, i); ++j) {
				int locMpcNb = scomm->mpcNb(i, j);
				if (flag == -1) tmpStatus[locMpcNb] = (rInfo.data[j] > -1) ? bool(rInfo.data[j]) : tmpStatus[locMpcNb];
				else // XXXX
					tmpStatus[locMpcNb] = (flag == 1) ? (tmpStatus[locMpcNb] || bool(rInfo.data[j])) : (
							tmpStatus[locMpcNb] && bool(rInfo.data[j]));
			}
		}
	}

	bool print_debug = false;
	statusChange = false;
	for (i = 0; i < numMPC; ++i) {
		if (getFetiInfo().contactPrintFlag && mpcMaster[i]) {
			if (!mpc[i]->active && !tmpStatus[i]) {
				std::cerr << "-";
				if (print_debug)
					std::cerr << " recvMpcStatus: sub = " << subNum() << ", mpc = " << localToGlobalMPC[i]
					          << std::endl;
			}
			else if (mpc[i]->active && tmpStatus[i]) {
				std::cerr << "+";
				if (print_debug)
					std::cerr << " recvMpcStatus: sub = " << subNum() << ", mpc = " << localToGlobalMPC[i]
					          << std::endl;
			}
		}
		mpc[i]->active = !tmpStatus[i];
		if (mpcStatus2[i] == mpc[i]->active) statusChange = true;
	}
}

template<class Scalar>
void
FetiSub<Scalar>::updateActiveSet(Scalar *v, double tol, int flag, bool &statusChange) {
	// flag = 0 : dual planing
	// flag = 1 : primal planing
	int *chgstatus = (int *) alloca(numMPC * sizeof(int));
	for (int i = 0; i < numMPC; ++i) chgstatus[i] = -1;  // set to 0 to remove, 1 to add

	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		if (mpc[locMpcNb]->type == 1) { // inequality constraint requiring planing

			if (flag == 0) { // dual planing
				if (mpcStatus1[locMpcNb] ==
				    1) { // active set expansion only: if constraint was initially active then it will not change status
					// if active and lambda < 0 then remove from active set
					if (mpc[locMpcNb]->active && ScalarTypes::lessThan(v[scomm->mapT(SComm::mpc, i)], tol))
						chgstatus[locMpcNb] = 1;
					// if not active and lambda >= 0 then add to the active set
					if (!mpc[locMpcNb]->active &&
					    ScalarTypes::greaterThanEq(v[scomm->mapT(SComm::mpc, i)], tol))
						chgstatus[locMpcNb] = 0;
				}
			} else { // primal planing
				if (mpcStatus1[locMpcNb] ==
				    0) { // active set contraction only: if constraint was initially inactive then it will not change status
					// if not active and w <= 0 then add to active set
					if (!mpc[locMpcNb]->active && ScalarTypes::lessThanEq(v[scomm->mapT(SComm::mpc, i)], tol))
						chgstatus[locMpcNb] = 0;
					// if active and w > 0 then remove from the active set
					if (mpc[locMpcNb]->active && ScalarTypes::greaterThan(v[scomm->mapT(SComm::mpc, i)], tol))
						chgstatus[locMpcNb] = 1;
				}
			}

		}
	}

	statusChange = false;
	for (int i = 0; i < numMPC; ++i) {
		if (chgstatus[i] > -1) {
			statusChange = true;
			mpc[i]->active = !bool(chgstatus[i]);
			if (getFetiInfo().contactPrintFlag && mpcMaster[i]) {
				if (chgstatus[i] == 0)
					std::cerr << "-";
				else std::cerr << "+";
			}
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::mergeUr(Scalar *ur, Scalar *uc, Scalar *u, Scalar *lambda) {
	int i, iNode;
	int rDofs[DofSet::max_known_dof];
	int oDofs[DofSet::max_known_dof];
	for(int dof = 0; dof < ccToC.size(); ++dof)
		u[ccToC[dof]] = ur[dof];

	int j;
	int iOff = 0;
	for (i = 0; i < numCRN; ++i) {
		int nd = get_c_dsa()->number(cornerNodes[i], cornerDofs[i], oDofs);
		for (j = 0; j < nd; ++j) {
			if (cornerEqNums[iOff + j] > -1)
				u[oDofs[j]] = uc[cornerEqNums[iOff + j]];
			else
				u[oDofs[j]] = 0.0;
		}
		iOff += nd;
	}

	// Primal augmentation 030314 JAT
	if(Ave.cols() > 0) {
		int nCor = this->Krc?this->Krc->numCol() : 0;
		int nAve = Ave.cols();
		for(int i = 0; i < ccToC.size(); ++i)
			for(int j = 0; j < nAve; ++j)
				u[ccToC[i]] += Ave[j][i]*uc[cornerEqNums[nCor+j]];
	}

	// extract uw
	Scalar *uw = (Scalar *) dbg_alloca(numWIdof * sizeof(Scalar));
	for (i = 0; i < scomm->lenT(SComm::wet); ++i)
		uw[scomm->wetDofNb(i)] = lambda[scomm->mapT(SComm::wet, i)];

//	for (i = 0; i < numWInodes; ++i) {
//		DofSet thisDofSet = wetInterfaceDofs[i]; // (*c_dsa)[wetInterfaceNodes[i]];
//		int nd = thisDofSet.count();
//		dsa->number(wetInterfaceNodes[i], thisDofSet, rDofs);
//		c_dsa->number(wetInterfaceNodes[i], thisDofSet, oDofs);
//		for (j = 0; j < nd; ++j) {
//			u[oDofs[j]] = uw[wetInterfaceMap[rDofs[j]]];
//		}
//	}
	if(numWInodes != 0)
		throw "Wet interface is not supported anymore. Work on the above lines if you want it.";

	// keep a local copy of the lagrange multipliers
	setLocalLambda(lambda);
}

template<class Scalar>
void
FetiSub<Scalar>::factorKii() {
	if (KiiSolver) KiiSolver->factor();
}

template<class Scalar>
void
FetiSub<Scalar>::setLocalLambda(Scalar *_localLambda) {
	if (localLambda) delete[] localLambda;
	localLambda = new double[totalInterfSize];
	for (int i = 0; i < totalInterfSize; ++i) localLambda[i] = ScalarTypes::Real(_localLambda[i]);
}

template<class Scalar>
void
FetiSub<Scalar>::multfc(const VectorView<const Scalar> &fr, /*Scalar *fc,*/ const VectorView<const Scalar> &lambda) const {
	Scalar v[Ave.cols()];
	VectorView<Scalar> t(v, Ave.cols(), 1);
	Eigen::Matrix<Scalar, Eigen::Dynamic, 1> force(localLen());

	force = -fr;
	force -= B*lambda;
	force -= Bm*lambda;


	if (numWIdof) Krw->multAddNew(localw.data(), force.data());  // coupled_dph: force += Krw * uw

	// Primal augmentation 072513 JAT
	if (Ave.cols() > 0) {
		t = this->Eve.transpose() * force;
		force.noalias() -= Ave * t;
	}

	if (this->Krr) this->Krr->solveInPlace(force);

	// Extra orthogonaliztion for stability  072216 JAT
	if (Ave.cols() > 0) {
		t = Ave.transpose() * force;
		force.noalias() -= Ave * t;
	}

	// re-initialization required for mpc/contact
	fcstar.assign(Src->numCol(), 0.0);

	// fcstar = - (Krr^-1 Krc)^T fr
	//        = - Krc^T (Krr^-1 fr)
	//        = Src force
	if (Src) Src->multAdd(force.data(), fcstar.data());

	// for coupled_dph add fcstar -= Kcw Bw uw
	if (numWIdof) {
		if (Kcw) Kcw->mult(localw.data(), fcstar.data(), -1.0, 1.0);
		if (Kcw_mpc) Kcw_mpc->multSubWI(localw.data(), fcstar.data());
	}

	VectorView<Scalar> fcs(fcstar.data(), fcstar.size());
	// add Bc^(s)^T lambda
	fcs += Bc*lambda;
}

template<class Scalar>
void
FetiSub<Scalar>::multAddBrT(const Scalar *interfvec, Scalar *localvec, Scalar *uw) const {
	VectorView<Scalar> locF(localvec, localLen());
	VectorView<const Scalar> lambda(interfvec, totalInterfSize);
	VectorView<Scalar> w(uw, numWIdof);
	locF += B*lambda;
	w -= Bw*lambda;
	locF += Bm*lambda;

	// coupled_dph: localvec -= Krw * uw
	if (Krw) Krw->multAddNew(uw, localvec);
}

template<class Scalar>
void
FetiSub<Scalar>::multBr(const Scalar *localvec, Scalar *interfvec, const Scalar *_uc, const Scalar *uw) const {
	VectorView<const Scalar> loc_u(localvec, localLen());
	VectorView<Scalar> u_interf(interfvec, totalInterfSize);
	VectorView<const Scalar> w(uw, numWIdof);
	VectorView<const Scalar> uc(_uc, Src->numCol());

	u_interf = B.transpose()*loc_u;
	u_interf -= Bw.transpose()*w;
	u_interf += Bc.transpose()*uc;
}

template<class Scalar>
void
FetiSub<Scalar>::fetiBaseOp(GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec) const {
	// Add the interface vector (interfvec) contribution to localvec
	int iDof;
	for (iDof = 0; iDof < totalInterfSize; ++iDof)
		localvec[allBoundDofs[iDof]] += interfvec[iDof];

	// solve for localvec
	if (s) s->reSolve(localvec);

	// redistribute the solution to the interface
	for (iDof = 0; iDof < totalInterfSize; ++iDof)
		interfvec[iDof] = localvec[allBoundDofs[iDof]];

}

template<class Scalar>
void
FetiSub<Scalar>::fetiBaseOp(GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec,
                                 Scalar *beta) const {
	auto &mpc = this->mpc;
	// multiply Q*beta to get the MPC force contribution
	if (numMPC > 0) {
		int i, iMPC;
		for (iMPC = 0; iMPC < numMPC; ++iMPC)
			for (i = 0; i < mpc[iMPC]->nterms; ++i) {
				int cdof = mpc[iMPC]->terms[i].cdof;
				if (cdof < 0) continue;
				localvec[cdof] += mpc[iMPC]->terms[i].coef * beta[localToGlobalMPC[iMPC]];
			}
	}

	// Add the interface vector (interfvec) contribution to localvec
	int iDof;
	for (iDof = 0; iDof < totalInterfSize; ++iDof)
		localvec[allBoundDofs[iDof]] += interfvec[iDof];

	// solve for localvec
	if (s) s->reSolve(localvec);

	// redistribute the solution to the interface
	for (iDof = 0; iDof < totalInterfSize; ++iDof)
		interfvec[iDof] = localvec[allBoundDofs[iDof]];
}

template<class Scalar>
void
FetiSub<Scalar>::fetiBaseOp(Scalar *uc, GenSolver<Scalar> *s, Scalar *localvec, Scalar *interfvec) const {
	Scalar v[this->Ave.cols()];
	VectorView<Scalar> t(v, this->Ave.cols(), 1);
	VectorView<Scalar> l(localvec, this->Ave.rows(), 1);

	// localvec += Br^T * interfvec
	multAddBrT(interfvec, localvec);

	// Primal augmentation 072513 JAT
	if (this->Ave.cols() > 0) {
		t = this->Eve.transpose() * l;
		l.noalias() -= this->Ave * t;
	}

	// localvec = Krr^-1 * localvec
	if (s) s->reSolve(localvec);

	// Extra orthogonaliztion for stability  072216 JAT
	if (this->Ave.cols() > 0) {
		t = this->Ave.transpose() * l;
		l.noalias() -= this->Ave * t;
	}

	// interfvec = Br * localvec
	multBr(localvec, interfvec, uc);
}

template<class Scalar>
void
FetiSub<Scalar>::fetiBaseOpCoupled1(GenSolver<Scalar> *s, Scalar *localvec, const Scalar *interfvec,
                                         FSCommPattern<Scalar> *wiPat) const {
	auto &localw = this->localw;
	// localvec += Br^T * interfvec
	multAddBrT(interfvec, localvec, localw.data());

	// solve for localvec
	if (s) s->reSolve(localvec);

	if (numWIdof) {
		int i, j;
		// compute Kww uw for this subdomain
		for (i = 0; i < numWIdof; ++i) localw_copy[i] = 0.0;
		this->Kww->mult(localw.data(), localw_copy.data());  // localw_copy = - Kww * uw

		// compute Kww uw to send to neighbors
//		if (getFetiInfo().fsi_corner == 0)
		// TODO Bring this back with wweight
		if(false)
			for (i = 0; i < scomm->numT(SComm::fsi); ++i) {
				if (subNum() != scomm->neighbT(SComm::fsi, i)) {
					FSSubRecInfo<Scalar> sInfo = wiPat->getSendBuffer(subNum(), scomm->neighbT(SComm::fsi, i));
					for (j = 0; j < numNeighbWIdof[i]; ++j) sInfo.data[j] = 0.0;
					neighbKww->multAdd(localw.data(), sInfo.data.data(), glToLocalWImap, neighbGlToLocalWImap[i]);
				} else {
					neighbKww->multAdd(localw.data(), localw_copy.data(), glToLocalWImap);
				}
			}
	}
}


template<class Scalar>
void
FetiSub<Scalar>::fetiBaseOpCoupled2(const Scalar *uc, const Scalar *localvec, Scalar *interfvec,
                                         FSCommPattern<Scalar> *wiPat, const Scalar *fw) const {
	// coupled_dph
	if (numWIdof) {
		int i, j;
		auto &Krw = this->Krw;
		auto &Kcw = this->Kcw;
		auto &Kcw_mpc = this->Kcw_mpc;


//		// TODO Bring this back with wweight
//		if (getFetiInfo().fsi_corner == 0)
//			for (i = 0; i < scomm->numT(SComm::fsi); ++i) {
//				if (subNum() != scomm->neighbT(SComm::fsi, i)) {
//					FSSubRecInfo<Scalar> rInfo = wiPat->recData(scomm->neighbT(SComm::fsi, i), subNum());
//					for (j = 0; j < numWIdof; ++j) localw_copy[j] += rInfo.data[j] / wweight[j];
//				}
//			}

		if (Krw) Krw->transposeMultSubNew(localvec, localw_copy.data()); // localw_copy -= Krw^T * localvec
		int numCDofs = Src->numCol();
		Scalar *ucLocal = (Scalar *) dbg_alloca(sizeof(Scalar) * numCDofs);
		for (i = 0; i < numCDofs; ++i) {
			if (cornerEqNums[i] > -1) ucLocal[i] = uc[cornerEqNums[i]];
			else ucLocal[i] = 0.0;
		}
		if (Kcw) Kcw->trMult(ucLocal, localw_copy.data(), -1.0, 1.0);  // localw_copy -= Kcw^T Bc uc
		if (Kcw_mpc) Kcw_mpc->transposeMultSubtractWI(ucLocal, localw_copy.data());
		if (fw) for (i = 0; i < numWIdof; ++i) localw_copy[i] += fw[i];  // localw_copy += fw
	}

	// interfvec = Br * localvec
	multBr(localvec, interfvec, uc, localw_copy.data());
}

template<class Scalar>
void
FetiSub<Scalar>::sendDiag(GenSparseMatrix<Scalar> *s, FSCommPattern<Scalar> *vPat) {
	int iDof = 0;
	for (int i = 0; i < scomm->numT(SComm::all); ++i) {
		FSSubRecInfo<Scalar> sInfo = vPat->getSendBuffer(subNum(), scomm->neighbT(SComm::all, i));
		for (int j = 0; j < scomm->lenT(SComm::all, i); ++j) {
			switch (boundDofFlag[iDof]) {
				case 0:
					scaling[iDof] = sInfo.data[j] = (s) ? s->diag(scomm->boundDofT(SComm::all, i, j)) : 1.0;
					break;
				case 1: // wet interface
					scaling[iDof] = sInfo.data[j] = (this->Kww) ? this->Kww->diag(-1 - scomm->boundDofT(SComm::all, i, j)) : 1.0;
					break;
				case 2:  // dual mpc
					scaling[iDof] = sInfo.data[j] = 1.0;
					break;
			}
			iDof++;
		}
	}

	int ndof = localLen();
	// use the kweight array also for LMPCs stiffness scaling/splitting for the primal method
	if ((getFetiInfo().augment == FetiInfo::WeightedEdges) ||
	    ((getFetiInfo().mpc_scaling == FetiInfo::kscaling) && (numMPC_primal > 0))) {
		kweight = new Scalar[ndof];  // used for WeightedEdges augmentation
		for (iDof = 0; iDof < ndof; ++iDof)
			kweight[iDof] = (s) ? s->diag(iDof) : 1.0;
	}
}

template<class Scalar>
void
FetiSub<Scalar>::multFcB(Scalar *p) {
	if (Src->numCol() == 0) return;
	if ((totalInterfSize == 0) || (localLen() == 0)) {
		for (int i = 0; i < Src->numCol(); ++i) fcstar[i] = 0.0;
		return;
	}

	// fcstar = - (Krr^-1 Krc)^T p
	//        = - Krc^T Krr^-1 p
	//        = - Acr p
	// TODO Change to fully use Eigen.
	GenStackFullM<Scalar> Acr(Src->numCol(), totalInterfSize, BKrrKrc.data());
	fcstar.resize(Src->numCol());
	Acr.mult(p, fcstar.data(), -1.0, 0.0);

	// for coupled_dph add fcstar += Kcw Bw uw
	if (numWIdof) {
		for (int i = 0; i < scomm->lenT(SComm::wet); ++i)
			localw[scomm->wetDofNb(i)] = p[scomm->mapT(SComm::wet, i)];
		if (Kcw) Kcw->mult(localw.data(), fcstar.data(), -1.0, 1.0);
		if (Kcw_mpc) Kcw_mpc->multSubWI(localw.data(), fcstar.data());
	}

	// fcstar += Bc^(s)^T p
	bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool) * numMPC);
	for (int i = 0; i < numMPC; ++i) mpcFlag[i] = true;
	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		if (mpcFlag[locMpcNb]) {
			const auto &m = mpc[locMpcNb];
			for (int k = 0; k < m->nterms; k++) {
				int dof = (m->terms)[k].dof;
				if ((dof >= 0) && (cornerMap[dof] >= 0))
					fcstar[cornerMap[dof]] += p[scomm->mapT(SComm::mpc, i)] * (m->terms)[k].coef;
			}
			mpcFlag[locMpcNb] = false;
		}
	}
}

template <typename S>
using DynMat = Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>;
template <typename S>
using DynVec = Eigen::Matrix<S, Eigen::Dynamic, 1>;

template<class Scalar>
void
FetiSub<Scalar>::formKccStar() {
	auto &BKrrKrc = this->BKrrKrc;
	auto &Krw = this->Krw;
	auto &localw = this->localw;
	// resize Kcc if necessary to add space for "primal" mpcs and augmentation
	if (this->Src->numCol() != this->Kcc->dim()) {
		std::unique_ptr<GenAssembledFullM<Scalar>> Kcc_copy = std::move(this->Kcc);
		this->Kcc = std::make_unique<GenAssembledFullM<Scalar>>(this->Src->numCol(), cornerMap);
		for (int i = 0; i < numCRNdof; ++i) for (int j = 0; j < numCRNdof; ++j) (*this->Kcc)[i][j] = (*Kcc_copy)[i][j];
	}


	// add in MPC coefficient contributions for Kcc^(s).
	if (numMPC_primal > 0)
		assembleMpcIntoKcc();

	if (this->Krr == 0) return;

	// Kcc* -> Kcc - Krc^T Krr^-1 Krc

	// first, perform Krc I = iDisp
	// which extracts the correct rhs vectors for the forward/backwards

	int nRHS = this->Src->numCol();
	Scalar **iDisp = new Scalar *[nRHS];
	Scalar *firstpointer = new Scalar[nRHS * nRHS];

	int numEquations = this->Krr->neqs();
	BKrrKrc.resize(totalInterfSize, nRHS);
	Scalar *thirdpointer = new Scalar[nRHS * numEquations];
	Scalar **KrrKrc = (Scalar **) dbg_alloca(nRHS * sizeof(Scalar *));
	//if(nRHS*numEquations == 0)
	//  fprintf(stderr, "We have a zero size %d %d %d\n",numEquations,totalInterfSize,nRHS);

	int iRHS, iDof;
	for (iRHS = 0; iRHS < nRHS; ++iRHS) {
		iDisp[iRHS] = firstpointer + iRHS * nRHS;
		KrrKrc[iRHS] = thirdpointer + iRHS * numEquations;
		for (iDof = 0; iDof < numEquations; ++iDof)
			KrrKrc[iRHS][iDof] = 0.0;
	}
	if (this->Src) this->Src->multIdentity(KrrKrc);

	// 070213 JAT
	if (this->Src && (getFetiInfo().augmentimpl == FetiInfo::Primal)) {
		int nAve, nCor;
		nCor = this->Krc ? this->Krc->numCol() : 0;
		nAve = this->Src->numCol() - nCor;
		if (nAve) {
			_AVMatrix<Scalar> Kave;
			Kave.resize(numEquations, nAve);
			DynVec<Scalar> pv(nAve);

			DynMat<Scalar> AKA(nAve, nAve);
			Scalar s;
			this->Ave.resize(numEquations, nAve);
			for (int i = 0; i < nAve; ++i) {
				for (int j = 0; j < numEquations; ++j)
					this->Ave[i][j] = KrrKrc[nCor + i][j];
				this->Ave.col(i).normalize();
			}
			this->Eve = this->Ave;
			this->Krr->reSolve(nAve, this->Eve[0]);
			AKA.noalias() = this->Eve.transpose() * this->Ave;

			// This does in-place solve.
			this->Eve.transpose() = AKA.ldlt().solve(this->Eve.transpose());

			for (int i = 0; i < nAve; ++i)
				for (int j = 0; j < nCor; ++j) {
					s = 0.0;
					for (int k = 0; k < numEquations; ++k)
						s += KrrKrc[j][k] * this->Ave[i][k];
					(*this->Kcc)[nCor + i][j] = s;
					(*this->Kcc)[j][nCor + i] = s;
				}
			for (int i = 0; i < nAve; ++i) {
				this->KrrSparse->mult(this->Ave[i], Kave.col(i).data());
				for (int j = 0; j < nAve; ++j) {
					s = 0.0;
					for (int k = 0; k < numEquations; ++k)
						s += Kave[i][k] * this->Ave[j][k];
					(*this->Kcc)[nCor + i][nCor + j] = s;
				}
			}

			int nz = 0;
			for (int i = 0; i < nAve; ++i)
				for (int k = 0; k < numEquations; ++k)
					if (std::abs(Kave[i][k]) > 0.0) nz++;

			std::vector<int> KACount(nAve);
			std::vector<int> KAList(nz);
			std::vector<Scalar> KACoefs(nz);
			nz = 0;
			for (int i = 0; i < nAve; ++i) {
				KACount[i] = 0;
				for (int k = 0; k < numEquations; ++k)
					if (std::abs(Kave[i][k]) > 0.0) {
						KACount[i]++;
						KAList[nz] = k;
						KACoefs[nz] = Kave[i][k];
						nz++;
					}
			}

			this->Grc = std::make_unique<GenCuCSparse<Scalar>>(nAve, numEquations, KACount, KAList, std::move(KACoefs));

			if (this->Src->num() == 2)
				this->Src->setSparseMatrix(1, this->Grc);
			else if (this->Src->num() == 1 && nCor == 0)
				this->Src->setSparseMatrix(0, this->Grc);
			else {
				fprintf(stderr, "unsupported number of blocks in Src\n");
				exit(1);
			}

			for (int i = 0; i < nRHS; ++i)
				for (int j = 0; j < numEquations; ++j)
					KrrKrc[i][j] = 0.0;

			this->Src->multIdentity(KrrKrc);

			Scalar vt[nAve];
			VectorView<Scalar> v{vt, nAve};
			for (int j = 0; j < nRHS; j++) {
				v = this->Eve.transpose() * VectorView<Scalar>{KrrKrc[j], numEquations};
				VectorView<Scalar>{KrrKrc[j], numEquations} -= this->Ave * v;
			}
		}
	}

	// KrrKrc <- Krr^-1 Krc
	if (this->Krr) this->Krr->reSolve(nRHS, KrrKrc); // this can be expensive when nRHS is large eg for coupled

	// -Krc^T KrrKrc
	for (iRHS = 0; iRHS < nRHS; ++iRHS)
		for (iDof = 0; iDof < nRHS; ++iDof)
			iDisp[iRHS][iDof] = 0.0;

	// Multiple RHS version of multSub: iDisp <- -Krc^T KrrKrc
	if (this->Src) this->Src->multSub(nRHS, KrrKrc, iDisp);

	if (this->Kcc) this->Kcc->add(iDisp);

	delete[] iDisp;
	delete[] firstpointer;
	auto &mpc = this->mpc;

	int k;
	for (iRHS = 0; iRHS < nRHS; ++iRHS) {
		bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool) * numMPC);
		for (int i = 0; i < numMPC; ++i) mpcFlag[i] = true;
		bool *wiFlag = (bool *) dbg_alloca(sizeof(bool) * numWIdof);
		for (int i = 0; i < numWIdof; ++i) wiFlag[i] = true;

		if (Krw) Krw->transposeMultNew(KrrKrc[iRHS], localw.data());

		for (iDof = 0; iDof < totalInterfSize; iDof++) {
			switch (boundDofFlag[iDof]) {
				case 0:
					BKrrKrc(iDof, iRHS) = KrrKrc[iRHS][allBoundDofs[iDof]];
					break;
				case 1: { // wet interface
					int windex = -1 - allBoundDofs[iDof];
					if (wiFlag[windex]) {
						BKrrKrc(iDof, iRHS) = -localw[-1 - allBoundDofs[iDof]];
						wiFlag[windex] = false;
					} else BKrrKrc(iDof, iRHS) = 0.0;
				}
					break;
				case 2: { // dual mpc
					int locMpcNb = -1 - allBoundDofs[iDof];
					const auto &m = mpc[locMpcNb];
					BKrrKrc(iDof, iRHS) = 0.0;
					if (mpcFlag[locMpcNb]) {
						for (k = 0; k < m->nterms; k++) {
							int cc_dof = (m->terms)[k].ccdof;
							if (cc_dof >= 0) BKrrKrc(iDof, iRHS) += KrrKrc[iRHS][cc_dof] * (m->terms)[k].coef;
						}
						mpcFlag[locMpcNb] = false;
					}
				}
					break;
			}
		}
	}
	delete[] thirdpointer;
}

template<class Scalar>
void
FetiSub<Scalar>::multDiagKbb(const Scalar *u, Scalar *Pu) const {
	// KHP: 02-02-99
	// boundMap    = from boundary number to subdomain number
	// internalMap = from internal number to subdomain number
	// invBoundMap =  from all subdomain dof to a unique boundary number
	// allBoundDofs = indices of B, from lambda numbering to numbering of entire
	//                subdomain
	// invBoundMap[allBoundDofs[iDof]] = from lambda numbering directly to
	//                                   boundary numbering
	// Kbb = works only on the numbering of the boundary dofs

	Scalar *v = (Scalar *) dbg_alloca(sizeof(Scalar) * boundLen);
	Scalar *res = (Scalar *) dbg_alloca(sizeof(Scalar) * boundLen);

	// XML Karim noted that his does not work for contact...

	int iDof;
	for (iDof = 0; iDof < boundLen; ++iDof)
		v[iDof] = res[iDof] = 0.0;

	for (iDof = 0; iDof < totalInterfSize; ++iDof)
		v[dualToBoundary[iDof]] += u[iDof] * scaling[iDof];

	// Perform diagonal multiplication
	Kbb->multDiag(v, res);

	for (iDof = 0; iDof < totalInterfSize; ++iDof)
		Pu[iDof] = res[dualToBoundary[iDof]] * scaling[iDof];
}

template<class Scalar>
void
FetiSub<Scalar>::multKbb(const Scalar *u, Scalar *Pu, Scalar *deltaU, Scalar *deltaF, bool errorFlag) {
	// KHP and DJR: 3-26-98
	// multKbb has been modified to compute subdomain primal residual using
	// deltaF (in addition to deltaU which is the displacement correction
	// and of course the lumped or dirichlet preconditioning)

	// If we are computing lumped preconditioner, we compute deltaFi using Kib
	// but deltaUi is equal to zero. For the dirichlet preconditioner, deltaUi
	// is computed but deltaFi is set equal to zero.

	// deltaFi = internal primal residual
	// deltaUi = internal displacement correction

	// boundMap    = from boundary number to subdomain number
	// internalMap = from internal number to subdomain number
	// invBoundMap  = from all subdomain dof to a unique boundary number
	// invInternalMap  =  from all subdomain dof to a unique internal number
	// allBoundDofs = indices of B, from lambda numbering to numbering of entire
	//                subdomain

	// dualToBoundary = from lambda numbering directly to boundary numbering

	// Kii = works only on the numbering of the internal dofs
	// Kbb = works only on the numbering of the boundary dofs
	// Kib = operates on the boundary numbering and returns with internal numbering

	Scalar *v = (Scalar *) dbg_alloca(sizeof(Scalar) * boundLen);
	Scalar *res = (Scalar *) dbg_alloca(sizeof(Scalar) * boundLen);

	int iDof;
	for (iDof = 0; iDof < boundLen; ++iDof) {
		v[iDof] = res[iDof] = 0.0;
	}
	// Karim added the following lines. I Don't think they are necessary XML
	if (deltaF && deltaU)
		for (iDof = 0; iDof < localLen(); ++iDof)
			deltaU[iDof] = deltaF[iDof] = 0;

	for (iDof = 0; iDof < totalInterfSize; ++iDof) {
		v[dualToBoundary[iDof]] += u[iDof] * scaling[iDof];
		if (deltaU)
			deltaU[allBoundDofs[iDof]] = -v[dualToBoundary[iDof]];
	}

	this->Kbb->mult(v, res);

	Scalar *iDisp = 0; // only allocate if necessary and don't use alloca (too big)

	if ((getFetiInfo().precno == FetiInfo::dirichlet) || deltaF) {
		iDisp = new Scalar[internalLen];
		for (iDof = 0; iDof < internalLen; ++iDof) iDisp[iDof] = 0.0;
		if (Kib) Kib->transposeMultAdd(v, iDisp);
	}

	if (getFetiInfo().precno == FetiInfo::dirichlet) {
		if (KiiSolver) KiiSolver->reSolve(iDisp);
		if (Kib) Kib->multSub(iDisp, res);
		if (deltaU)
			for (iDof = 0; iDof < internalLen; ++iDof)
				deltaU[internalMap[iDof]] = iDisp[iDof];
	} else {
		if (deltaF)
			for (iDof = 0; iDof < internalLen; ++iDof) {
				deltaF[internalMap[iDof]] = iDisp[iDof];
			}
	}
	if (iDisp) delete[] iDisp;

	if (deltaF) {
		for (iDof = 0; iDof < boundLen; ++iDof)
			deltaF[boundMap[iDof]] = res[iDof];
	}

	// Return preconditioned u
	for (iDof = 0; iDof < totalInterfSize; ++iDof)
		Pu[iDof] = res[dualToBoundary[iDof]] * scaling[iDof];
}

template<class Scalar>
void
FetiSub<Scalar>::multKbbCoupled(const Scalar *u, Scalar *Pu, Scalar *deltaF, bool errorFlag) {
	// modified version of multKbb for coupled_dph with primal wet interface dofs
	// included in Kii. Currently preconditioner is un-coupled, ie fluid-structure interaction
	// is ignored in the Kii solve.

	auto &localw = this->localw;
	Scalar *v = (Scalar *) dbg_alloca(sizeof(Scalar) * boundLen);
	Scalar *res = (Scalar *) dbg_alloca(sizeof(Scalar) * boundLen);
	deltaFwi.resize(numWIdof);

	int i, iDof;
	for (iDof = 0; iDof < boundLen; ++iDof) { v[iDof] = res[iDof] = 0.0; }
	for (iDof = 0; iDof < localLen(); ++iDof) deltaF[iDof] = 0;
	for (i = 0; i < numWIdof; ++i) localw[i] = 0;

	for (iDof = 0; iDof < totalInterfSize; ++iDof) {
		if (allBoundDofs[iDof] >= 0)
			v[dualToBoundary[iDof]] += u[iDof] * scaling[iDof];
		else if (boundDofFlag[iDof] == 1) { // wet interface
			int windex = -1 - allBoundDofs[iDof];
#ifdef HB_COUPLED_PRECOND
			localw[windex] = (getFetiInfo().splitLocalFsi) ? -u[iDof]*scaling[iDof] : -u[iDof];
#else
			localw[windex] = -u[iDof] * scaling[iDof];
#endif
			if (getFetiInfo().precno == FetiInfo::noPrec)
				localw[windex] = u[iDof];
			deltaFwi[windex] = -u[iDof];
		}
	}

	this->Kbb->mult(v, res);

	Scalar *iDisp = (Scalar *) dbg_alloca(sizeof(Scalar) * internalLen);
	for (iDof = 0; iDof < internalLen; ++iDof) iDisp[iDof] = 0.0;
	for (i = 0; i < numWIdof; ++i) iDisp[wiInternalMap[i]] = localw[i]; // coupled_dph
	if (Kib) Kib->transposeMultAdd(v, iDisp);

	if (getFetiInfo().precno == FetiInfo::dirichlet) {
		if (KiiSolver) KiiSolver->reSolve(iDisp);
		for (i = 0; i < numWIdof; ++i) localw[i] = -iDisp[wiInternalMap[i]];
		if (Kib) Kib->multSub(iDisp, res);
	} else {
		for (iDof = 0; iDof < internalLen; ++iDof)
			if (internalMap[iDof] > 0) deltaF[internalMap[iDof]] = iDisp[iDof]; // not including wet interface in deltaF
	}

	for (iDof = 0; iDof < boundLen; ++iDof) { deltaF[boundMap[iDof]] = res[iDof]; }

	// Return preconditioned u
	for (iDof = 0; iDof < totalInterfSize; ++iDof) {
		if (allBoundDofs[iDof] >= 0)
			Pu[iDof] = res[dualToBoundary[iDof]] * scaling[iDof];
		else if (boundDofFlag[iDof] == 1)
			Pu[iDof] = localw[-1 - allBoundDofs[iDof]] * scaling[iDof];
	}
}

template<class Scalar>
void
FetiSub<Scalar>::assembleMpcIntoKcc() {
	// compute mpcOffset for my Subdomain
	int mpcOffset = numCRNdof;

	int i, iMPC;
	for (iMPC = 0; iMPC < numMPC_primal; ++iMPC) {
		for (i = 0; i < mpc_primal[iMPC]->nterms; ++i) {
			int d = mpc_primal[iMPC]->terms[i].dof;
			int dof = mpc_primal[iMPC]->terms[i].ccdof;
			if ((dof < 0) && (d >= 0) && !isWetInterfaceDof(d)) {
				int column = cornerMap[d];
				int row = mpcOffset + iMPC;
				if (row > this->Kcc->dim())
					std::cout << " *** ERROR: Dimension Error Row = " << row << " > " << this->Kcc->dim() << std::endl;
				if (column > this->Kcc->dim())
					std::cout << " *** ERROR: Dimension Error Col = " << column << " > " << this->Kcc->dim()
					          << std::endl;
				// i.e. an mpc touches a corner node that also has DBCs
				if (column >= 0) {
					(*this->Kcc)[row][column] += mpc_primal[iMPC]->terms[i].coef;
					(*this->Kcc)[column][row] += mpc_primal[iMPC]->terms[i].coef;
				}
			}
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::constructKcc() {
	Kcc = std::make_unique<GenAssembledFullM<Scalar>>(numCRNdof, cornerMap);
//	memK += numCRNdof * numCRNdof; // TODO Move/duplicate memory use variable ???
}

template<class Scalar>
void
FetiSub<Scalar>::constructKcw() {
	if (numWIdof && numCoarseDofs()) {
		Kcw = std::make_unique<GenCuCSparse<Scalar>>(getNodeToNode(), getDsa(), wetInterfaceMap.data(), cornerMap);
		Kcw->zeroAll();
//		memK += (Kcw) ? Kcw->size() : 0;
	}
}

template<class Scalar>
void
FetiSub<Scalar>::setMpcDiagCommSize(FSCommStructure *mpcDiagPat) const {
	for (int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc, i);
		int len = 0;
		if (subNum() != neighb) {
			for (int j = 0; j < scomm->lenT(SComm::mpc, i); ++j)
				len += this->mpc[scomm->mpcNb(i, j)]->gsize;
		}
		mpcDiagPat->setLen(subNum(), neighb, len);
	}
}

template<class Scalar>
void
FetiSub<Scalar>::sendMpcDiag(FSCommPattern<Scalar> *mpcDiagPat) {
	auto &mpc = this->mpc;
	for (int i = 0; i < numMPC; ++i) mpc[i]->initKsum();

	// Get the trace of the subdomain interfaces
	int iNeighb, iDof, j;
	for (iNeighb = 0; iNeighb < scomm->numT(SComm::mpc); ++iNeighb) {
		int neighb = scomm->neighbT(SComm::mpc, iNeighb);
		FSSubRecInfo<Scalar> sInfo = mpcDiagPat->getSendBuffer(subNum(), neighb);
		int nOff = 0;
		for (iDof = 0; iDof < scomm->lenT(SComm::mpc, iNeighb); ++iDof) {
			int locMpcNb = scomm->mpcNb(iNeighb, iDof);
			if (subNum() != neighb)
				for (j = 0; j < mpc[locMpcNb]->gsize; ++j) sInfo.data[nOff + j] = 0.0;
			for (j = 0; j < mpc[locMpcNb]->nterms; ++j) {
				int c_dof = mpc[locMpcNb]->terms[j].cdof;
				if (c_dof > -1) {
					int b_dof = invBoundMap[c_dof];
					mpc[locMpcNb]->k[j] = (Kbb) ? Kbb->diag(b_dof) : 1.0;
					if (ScalarTypes::norm(mpc[locMpcNb]->k[j]) < 1.0e-12)
						std::cerr << " *** WARNING: Kbb diagonal < 1.0e-12 \n";
					if (subNum() != neighb)
						sInfo.data[nOff + mpc[locMpcNb]->gi[j]] = mpc[locMpcNb]->k[j];
					mpc[locMpcNb]->ksum[j] = mpc[locMpcNb]->k[j];
				} else {
					if (subNum() != neighb)
						sInfo.data[nOff + mpc[locMpcNb]->gi[j]] = 0.0;
					mpc[locMpcNb]->k[j] = 0.0;
					mpc[locMpcNb]->ksum[j] = 0.0;
				}
			}
			nOff += mpc[locMpcNb]->gsize;
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::collectScaling(FSCommPattern<Scalar> *vPat) {
	int offset = 0;
	int locLen = localLen();
	//Scalar *kSum = (Scalar *) dbg_alloca(sizeof(Scalar)*locLen);
	Scalar *kSum = new Scalar[locLen];
	for (int i = 0; i < locLen; ++i) kSum[i] = 0.0;
#ifdef HB_COUPLED_PRECOND
	kSumWI = new Scalar[numWIdof]; //stores the sum of the neighbourg stiffess otherwise 0.0 if not shared
#else
	Scalar *kSumWI = (Scalar *) dbg_alloca(numWIdof * sizeof(Scalar));
#endif
	bool *wflag = (bool *) dbg_alloca(numWIdof * sizeof(bool)); //HB: to avoid setting wweight multiple times
	for (int i = 0; i < numWIdof; ++i) {
		kSumWI[i] = 0.0;
		wflag[i] = false;
	}

	int iSub, iDof;
	for (iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
		if (subNum() != scomm->neighbT(SComm::all, iSub)) {
			FSSubRecInfo<Scalar> rInfo = vPat->recData(scomm->neighbT(SComm::all, iSub), subNum());
			for (iDof = 0; iDof < scomm->lenT(SComm::all, iSub); ++iDof) {
				int bdof = scomm->boundDofT(SComm::all, iSub, iDof);
				if (bdof >= 0)
					kSum[bdof] += rInfo.data[iDof];
				else if (boundDofFlag[offset + iDof] == 1)
					kSumWI[-1 - bdof] += rInfo.data[iDof];
			}
		}
		offset += scomm->lenT(SComm::all, iSub);
	}

	offset = 0;
	for (iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
		FSSubRecInfo<Scalar> rInfo = vPat->recData(scomm->neighbT(SComm::all, iSub), subNum());
		for (iDof = 0; iDof < scomm->lenT(SComm::all, iSub); ++iDof) {
			int bdof = scomm->boundDofT(SComm::all, iSub, iDof);
			if (bdof >= 0)
				scaling[offset + iDof] = rInfo.data[iDof] / (scaling[offset + iDof] + kSum[bdof]);
			else if (boundDofFlag[offset + iDof] == 1) {
				scaling[offset + iDof] = scaling[offset + iDof] / (scaling[offset + iDof] + kSumWI[-1 - bdof]);
				if ((getFetiInfo().fsi_scaling == FetiInfo::kscaling) & !wflag[-1 - bdof]) {
					wweight[-1 - bdof] = (1. / scaling[offset + iDof]) / wweight[-1 - bdof];
					wflag[-1 - bdof] = true;
				}
			}
		}
		offset += scomm->lenT(SComm::all, iSub);
	}

	if (getFetiInfo().augment == FetiInfo::WeightedEdges) {
		for (iDof = 0; iDof < locLen; ++iDof)
			kweight[iDof] += kSum[iDof];
	}
	delete[] kSum;
}

template<class Scalar>
void
FetiSub<Scalar>::fScale(Scalar *locF, FSCommPattern<Scalar> *vPat, Scalar *locFw) {
	bool *isShared = (bool *) dbg_alloca(sizeof(bool) * totalInterfSize);

	int iDof, iSub;
	int offset = 0;
	for (iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
		FSSubRecInfo<Scalar> rInfo = vPat->recData(scomm->neighbT(SComm::all, iSub), subNum());
		for (iDof = 0; iDof < scomm->lenT(SComm::all, iSub); ++iDof) {
			int bdof = scomm->boundDofT(SComm::all, iSub, iDof);
			if (bdof >= 0)
				locF[bdof] += rInfo.data[iDof];
			else if (boundDofFlag[offset + iDof] == 1) {
				if (scomm->neighbT(SComm::all, iSub) != subNum()) { // wet interface
					locFw[-1 - bdof] += rInfo.data[iDof];
					isShared[offset + iDof] = true;
				} else isShared[offset + iDof] = false;
			}
		}
		offset += scomm->lenT(SComm::all, iSub);
	}

	Scalar *interfF = (Scalar *) dbg_alloca(sizeof(Scalar) * totalInterfSize);

	// Get the trace of the subdomain interfaces
	for (iDof = 0; iDof < totalInterfSize; ++iDof)
		if (allBoundDofs[iDof] >= 0)
			interfF[iDof] = locF[allBoundDofs[iDof]] * scaling[iDof];
		else if ((boundDofFlag[iDof] == 1) && isShared[iDof]) // wet interface
			interfF[iDof] = locFw[-1 - allBoundDofs[iDof]] * scaling[iDof];

	for (iDof = 0; iDof < totalInterfSize; ++iDof)
		if (allBoundDofs[iDof] >= 0)
			locF[allBoundDofs[iDof]] -= interfF[iDof];
		else if ((boundDofFlag[iDof] == 1) && isShared[iDof])  // wet interface
			locFw[-1 - allBoundDofs[iDof]] -= interfF[iDof];
}

template<class Scalar>
void
FetiSub<Scalar>::sendDeltaF(const Scalar *deltaF, FSCommPattern<Scalar> *vPat) {
	auto &deltaFmpc = this->deltaFmpc;
	int iDof = 0;
	for (int i = 0; i < scomm->numT(SComm::all); ++i) {
		FSSubRecInfo<Scalar> sInfo = vPat->getSendBuffer(subNum(), scomm->neighbT(SComm::all, i));
		for (int j = 0; j < scomm->lenT(SComm::all, i); ++j) {
			int bdof = scomm->boundDofT(SComm::all, i, j);
			switch (boundDofFlag[iDof]) {
				case 0: {
					if (deltaF) sInfo.data[j] = deltaF[bdof];
					else sInfo.data[j] = 0.0;
				}
					break;
				case 1: {  // wet interface
					int windex = -1 - bdof;
					sInfo.data[j] = this->deltaFwi[windex];
				}
					break;
				case 2: {  // dual mpc or contact
					int locMpcNb = -1 - bdof;
					sInfo.data[j] = (masterFlag[iDof]) ? deltaFmpc[locMpcNb] : -deltaFmpc[locMpcNb];
				}
					break;
			}
			iDof++;
		}
	}
}

template<class Scalar>
double
FetiSub<Scalar>::collectAndDotDeltaF(Scalar *deltaF, FSCommPattern<Scalar> *vPat) {
	// if there are more than 2 subdomains sharing a mpc define the norm
	// as (f1 - f2 - f3)^2 = f1^2 + f2^2 + f3^2 - 2f1f2 - 2f1f3 + 2f2f3 --> currently implemented
	auto &deltaFmpc = this->deltaFmpc;

	Scalar dot = 0;
	int i, iSub, jDof;

	if (deltaF) {
		for (i = 0; i < localLen(); ++i) {
			double dPrScal = 1.0 / this->densProjCoefficient(i);
			dot += dPrScal * dPrScal * deltaF[i] * ScalarTypes::conj(deltaF[i]);
		}
	}

	for (i = 0; i < numMPC; ++i)
		dot += deltaFmpc[i] * ScalarTypes::conj(deltaFmpc[i]);

	for (i = 0; i < numWIdof; ++i) //HB ... to be checked ...
		dot += this->deltaFwi[i] * ScalarTypes::conj(this->deltaFwi[i]) / this->wweight[i];

	int nbdofs = 0;
	for (iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
		FSSubRecInfo<Scalar> rInfo = vPat->recData(scomm->neighbT(SComm::all, iSub), subNum());
		for (jDof = 0; jDof < scomm->lenT(SComm::all, iSub); ++jDof) {
			int bdof = scomm->boundDofT(SComm::all, iSub, jDof);
			switch (boundDofFlag[nbdofs]) {
				case 0:
					if (deltaF) {
						double dPrScal = 1.0 / this->densProjCoefficient(bdof);
						dot += dPrScal * dPrScal * deltaF[bdof] * ScalarTypes::conj(rInfo.data[jDof]);
					}
					break;
					// do nothing for case 1 (wet interface)
				case 2: { // dual mpc
					if (subNum() != scomm->neighbT(SComm::all, iSub)) {
						int locMpcNb = -1 - bdof;
						dot += deltaFmpc[locMpcNb] * ScalarTypes::conj(rInfo.data[jDof]);
					}
				}
					break;
			}
			nbdofs++;
		}
	}
	return ScalarTypes::Real(dot);
}

template<class Scalar>
void
FetiSub<Scalar>::splitInterf(Scalar *subvec) const {
	Scalar *interfF = (Scalar *) alloca(scomm->sharedDOFsPlus->numConnect() * sizeof(Scalar));
	for (int iDof = 0; iDof < scomm->sharedDOFsPlus->numConnect(); ++iDof)
		interfF[iDof] =
				subvec[(*scomm->sharedDOFsPlus)[0][iDof]] / double(weightPlus[(*scomm->sharedDOFsPlus)[0][iDof]]);

	for (int iDof = 0; iDof < scomm->sharedDOFsPlus->numConnect(); ++iDof)
		subvec[(*scomm->sharedDOFsPlus)[0][iDof]] -= interfF[iDof];
}

template<class Scalar>
void FetiSub<Scalar>::initMpcScaling() {
	// sets scaling = 1.0 for all mpc virtual dofs, for use with generalized preconditioner
	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i)
		scaling[scomm->mapT(SComm::mpc, i)] = 1.0;
}

template<class Scalar>
void
FetiSub<Scalar>::initScaling() {
	int i;

	scaling.resize(totalInterfSize);
	for (i = 0; i < totalInterfSize; ++i) scaling[i] = Scalar(1.0);
	if (getFetiInfo().scaling == FetiInfo::tscaling)
		for (i = 0; i < scomm->lenT(SComm::std); ++i)
			scaling[scomm->mapT(SComm::std, i)] = 1.0 / getWeights()[scomm->boundDofT(SComm::std, i)];

	if (getFetiInfo().fsi_scaling == FetiInfo::tscaling)
		for (i = 0; i < scomm->lenT(SComm::wet); ++i)
			scaling[scomm->mapT(SComm::wet, i)] = 1.0 / wweight[scomm->wetDofNb(i)];

	if (getFetiInfo().mpc_precno == FetiInfo::diagCCt)
		for (i = 0; i < scomm->lenT(SComm::mpc); ++i)
			scaling[scomm->mapT(SComm::mpc, i)] = 1.0;
}

template<class Scalar>
void
FetiSub<Scalar>::collectMpcDiag(FSCommPattern<Scalar> *mpcDiagPat) {
	auto &mpc = this->mpc;
	int iNeighb, iDof, j;
	// 1) gets & sum the stiffness contribution from neighbourg subds
	for (iNeighb = 0; iNeighb < scomm->numT(SComm::mpc); ++iNeighb) {
		// ksum already contains its own contribution (i.e. initialized in sendMpcDiag)
		// -> loop only on neighboring subds
		int neighb = scomm->neighbT(SComm::mpc, iNeighb);
		if (subNum() != neighb) {
			FSSubRecInfo<Scalar> rInfo = mpcDiagPat->recData(neighb, subNum());
			int nOff = 0;
			for (iDof = 0; iDof < scomm->lenT(SComm::mpc, iNeighb); ++iDof) {
				int locMpcNb = scomm->mpcNb(iNeighb, iDof);
				for (j = 0; j < mpc[locMpcNb]->nterms; ++j) {
					mpc[locMpcNb]->ksum[j] += rInfo.data[nOff + mpc[locMpcNb]->gi[j]];
				}
				nOff += mpc[locMpcNb]->gsize;
			}
		}
	}

	// 2) adjust dual mpcs using subdomain multiplicity
	//    c(i) -> c(i).k(i,i)/sum[k(j,j)]
	if (getFetiInfo().mpc_scaling == FetiInfo::kscaling) {
		int iMPC, i;
#ifdef DEBUG_MPC
		std::cerr << "before k scaling: \n";
   for(iMPC = 0; iMPC < numMPC; ++iMPC) mpc[iMPC]->print();
#endif
		for (iMPC = 0; iMPC < numMPC; ++iMPC) {
			if (mpc[iMPC]->type == 2) continue; // bmpc
			for (i = 0; i < mpc[iMPC]->nterms; ++i) {
				if (ScalarTypes::norm(mpc[iMPC]->ksum[i]) < 1.0e-12) {
					//std::cerr << " *** WARNING: ksum = " << mpc[iMPC]->ksum[i] << ", cdof = " << mpc[iMPC]->terms[i].cdof << ", coef = " << mpc[iMPC]->terms[i].coef << std::endl;
					mpc[iMPC]->ksum[i] = 1.0;
					mpc[iMPC]->k[i] = 1.0;
				}
				mpc[iMPC]->terms[i].coef *= (mpc[iMPC]->k[i] / mpc[iMPC]->ksum[i]);
			}
		}
#ifdef DEBUG_MPC
		std::cerr << "after k scaling: \n";
   for(iMPC = 0; iMPC < numMPC; ++iMPC) mpc[iMPC]->print();
#endif
	}
}

// TODO Reuse the B matrices.
template<class Scalar>
void
FetiSub<Scalar>::applyBtransposeAndScaling(const Scalar *u, Scalar *v, Scalar *deltaU, Scalar *localw) const {
	int i, iDof, k;
	bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool) * numMPC);
	for (i = 0; i < numMPC; ++i) mpcFlag[i] = true;

	for (iDof = 0; iDof < totalInterfSize; ++iDof) {
		switch (boundDofFlag[iDof]) {
			case 0:
				v[dualToBoundary[iDof]] += u[iDof] * scaling[iDof];
				if (deltaU) deltaU[allBoundDofs[iDof]] = -v[dualToBoundary[iDof]];
				break;
			case 1: { // wet interface
				int windex = -1 - allBoundDofs[iDof];
				localw[windex] = u[iDof] * scaling[iDof];
				deltaFwi[windex] = u[iDof];
			}
				break;
			case 2: { // dual mpc
				int locMpcNb = -1 - allBoundDofs[iDof];
				if (mpcFlag[locMpcNb]) {
					const auto &m = mpc[locMpcNb];
					if (!mpc[locMpcNb]->active) {
						for (k = 0; k < m->nterms; k++) {
							int cdof = (m->terms)[k].cdof;
							if (cdof >= 0) { // mpc dof that exists
								Scalar coef =
										(m->terms)[k].coef / m->k[k]; // 1/m->k[k] = A, see generalized preconditioner
								if (invBoundMap[cdof] < 0)
									std::cerr << "error here in FetiSub<Scalar>::applyBtransposeAndScaling\n";
								v[invBoundMap[cdof]] += u[iDof] * coef * scaling[iDof];
							}
						}
					}
					mpcFlag[locMpcNb] = false;
				}
			}
				break;
		}
	}
}

// TODO Make use of the B matrices
template<class Scalar>
void
FetiSub<Scalar>::applyScalingAndB(const Scalar *res, Scalar *Pu, Scalar *localw) const {

	deltaFmpc.resize(numMPC); // only need to allocate 1st time (unless numMPC changes)
	bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool) * numMPC);
	for (int i = 0; i < numMPC; ++i) mpcFlag[i] = true;

	// Return preconditioned u
	for (int iDof = 0; iDof < totalInterfSize; ++iDof) {
		switch (boundDofFlag[iDof]) {
			case 0:
				Pu[iDof] = res[dualToBoundary[iDof]] * scaling[iDof];
				break;
			case 1: // wet interface
				Pu[iDof] = localw[-1 - allBoundDofs[iDof]] * scaling[iDof];
				break;
			case 2: { // dual mpc or contact
				int locMpcNb = -1 - allBoundDofs[iDof];
				const auto &m = mpc[locMpcNb];
				Pu[iDof] = 0.0;
				if (mpcFlag[locMpcNb]) deltaFmpc[locMpcNb] = 0.0;
				if (!mpc[locMpcNb]->active) {
					for (int k = 0; k < m->nterms; k++) {
						int cdof = (m->terms)[k].cdof;
						if (cdof > -1) { // mpc dof that exists
							Scalar coef = (m->terms)[k].coef;
							Pu[iDof] += res[invBoundMap[cdof]] * coef * scaling[iDof] /
							            m->k[k]; // 1/m->k[k] = A, see generalized preconditioner
							if (mpcFlag[locMpcNb]) deltaFmpc[locMpcNb] += res[invBoundMap[cdof]] * coef;
						}
					}
				}
				mpcFlag[locMpcNb] = false;
			}
				break;
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::sendMpcScaling(FSCommPattern<Scalar> *mpcPat) {
	int i, j;
	diagCCt.resize(numMPC); // compute weights for mpcs also
	for (i = 0; i < numMPC; ++i) {
		diagCCt[i] = 0.0;
		for (j = 0; j < mpc[i]->nterms; ++j)
			diagCCt[i] += mpc[i]->terms[j].coef * mpc[i]->terms[j].coef / mpc[i]->k[j]; //HB: for mpc kscaling;
	}

	int neighb;
	for (i = 0; i < scomm->numT(SComm::mpc); ++i) {
		if (subNum() != (neighb = scomm->neighbT(SComm::mpc, i))) {
			FSSubRecInfo<Scalar> sInfo = mpcPat->getSendBuffer(subNum(), neighb);
			for (j = 0; j < scomm->lenT(SComm::mpc, i); ++j)
				sInfo.data[j] = diagCCt[scomm->mpcNb(i, j)];
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::collectMpcScaling(FSCommPattern<Scalar> *mpcPat) {
	int i, j;
	Scalar *mpcSum = (Scalar *) dbg_alloca(sizeof(Scalar) * numMPC);
	for (i = 0; i < numMPC; ++i) mpcSum[i] = 0.0;
	int neighb;
	for (i = 0; i < scomm->numT(SComm::mpc); ++i) {
		if (subNum() != (neighb = scomm->neighbT(SComm::mpc, i))) {
			FSSubRecInfo<Scalar> rInfo = mpcPat->recData(neighb, subNum());
			for (j = 0; j < scomm->lenT(SComm::mpc, i); ++j)
				mpcSum[scomm->mpcNb(i, j)] += rInfo.data[j];
		}
	}
	for (i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		if (diagCCt[locMpcNb] + mpcSum[locMpcNb] != 0.0)
			scaling[scomm->mapT(SComm::mpc, i)] /= (diagCCt[locMpcNb] + mpcSum[locMpcNb]);
	}
	diagCCt.resize(0);
}

template<class Scalar>
void
FetiSub<Scalar>::initMpcStatus() {
	for (int i = 0; i < numMPC; ++i) {
		mpc[i]->active = false;
	}
}

template<class Scalar>
void
FetiSub<Scalar>::saveMpcStatus() {
	// this saves the status before first update iteration so it can be reset if nonmonotic
	if (!mpcStatus) mpcStatus = new int[numMPC];
	for (int i = 0; i < numMPC; ++i) {
		mpcStatus[i] = int(!mpc[i]->active);
	}
}

template<class Scalar>
void
FetiSub<Scalar>::restoreMpcStatus() {
	for (int i = 0; i < numMPC; ++i) {
		if (getFetiInfo().contactPrintFlag && mpcMaster[i]) {
			if (!mpc[i]->active && !mpcStatus[i]) std::cerr << "-";
			else if (mpc[i]->active && mpcStatus[i]) std::cerr << "+";
		}
		mpc[i]->active = bool(!mpcStatus[i]);
	}
}

template<class Scalar>
void
FetiSub<Scalar>::saveMpcStatus1() const {
	auto &mpc = this->mpc;
	if (!mpcStatus1) mpcStatus1 = new bool[numMPC];
	for (int i = 0; i < numMPC; ++i) mpcStatus1[i] = !mpc[i]->active;
}

template<class Scalar>
void
FetiSub<Scalar>::saveMpcStatus2() {
	auto &mpc = this->mpc;
	if (!mpcStatus2) mpcStatus2 = new bool[numMPC];
	for (int i = 0; i < numMPC; ++i) mpcStatus2[i] = !mpc[i]->active;
}

template<class Scalar>
void
FetiSub<Scalar>::cleanMpcData() {
	if (mpcStatus) {
		delete[] mpcStatus;
		mpcStatus = 0;
	}
	if (mpcStatus1) {
		delete[] mpcStatus1;
		mpcStatus1 = 0;
	}
	if (mpcStatus2) {
		delete[] mpcStatus2;
		mpcStatus2 = 0;
	}
}


template<class Scalar>
void
FetiSub<Scalar>::sendMpcStatus(FSCommPattern<int> *mpcPat, int flag) {
	auto &mpc = this->mpc;
	// if flag = -1 only send status of mpcMaster otherwise send -1
	// note: could use SComm::ieq list
	for (int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc, i);
		if (subNum() != neighb) {
			FSSubRecInfo<int> sInfo = mpcPat->getSendBuffer(subNum(), neighb);
			for (int j = 0; j < scomm->lenT(SComm::mpc, i); ++j) {
				int locMpcNb = scomm->mpcNb(i, j);
				if (flag == -1) sInfo.data[j] = (this->mpcMaster[locMpcNb]) ? int(!mpc[locMpcNb]->active) : -1;
				else
					sInfo.data[j] = int(!mpc[locMpcNb]->active);
			}
		}
	}
}


/** \details The degrees of freedom are relative to*/
void
FetiBaseSub::makeKccDofsMultiLevel(gsl::span<FetiBaseSub *> sd,
                                 int augOffset, Connectivity *subToEdge) {
	int numC = numCoarseDofs();
	cornerEqNums.resize(numC);

	// numbers the corner equations
	int offset = 0;
	for (int i = 0; i < numCRN; ++i) {
		int offset2 = 0;
		for (int j = 0; j < sd.size(); ++j) {
			auto &nodeMap = sd[j]->getGlobalToLocalNode();
			auto cornerEqs = sd[j]->get_c_dsa();
			if (nodeMap[glCornerNodes[i]] > -1) {
				int count = cornerEqs->number(nodeMap[glCornerNodes[i]], cornerDofs[i].list(),
				                              cornerEqNums.data() + offset);
				for (int k = 0; k < count; ++k)
					cornerEqNums[offset + k] += offset2;
				offset += count;
				break;
			}
			offset2 += cornerEqs->size();
		}
	}

	if (getFetiInfo().augmentimpl == FetiInfo::Primal) {
		int iEdgeN = 0;
		for (int iNeighb = 0; iNeighb < scomm->numNeighb; ++iNeighb) {
			if (scomm->isEdgeNeighb[iNeighb]) {
				int edgeIdx = augOffset + (*subToEdge)[subNum()][iEdgeN];
				int offset2 = 0;
				for (auto & iSub : sd) {
					GlobalToLocalMap &nodeMap = iSub->getGlobalToLocalNode();
					auto cornerEqs = iSub->get_c_dsa();
					if (nodeMap[edgeIdx] > -1) {
						int fDof = cornerEqs->firstdof(nodeMap[edgeIdx]);
						int count = edgeDofSize[iNeighb];
						for (int k = 0; k < count; ++k)
							cornerEqNums[offset + k] = fDof + k + offset2;
						offset += count;
						break;
					}
					offset2 += cornerEqs->size();
				}
				iEdgeN++;
			}
		}
	}
}

void
getOneDirection(double d, int i,int j, int k, int &nnum, int numWaves,
                std::vector<double> &wDir_x, std::vector<double> &wDir_y, std::vector<double> &wDir_z)
{
	double x = -0.5 + d*i;
	double y = -0.5 + d*j;
	double z = -0.5 + d*k;

	double l = sqrt(x*x+y*y+z*z);
	if (l == 0.0) l = 1.0;
	wDir_x[nnum] = x/l;
	wDir_y[nnum] = y/l;
	wDir_z[nnum] = z/l;
	nnum++;
	if(numWaves == 3) {
		double lyz = sqrt(y*y+z*z);
		if(lyz == 0) {
			wDir_x[nnum] = 0.0;
			wDir_y[nnum] = x/l;
			wDir_z[nnum] = 0.0;
			nnum++;
			wDir_x[nnum] = 0.0;
			wDir_y[nnum] = 0.0;
			wDir_z[nnum] = fabs(x/l);
			nnum++;
		}
		else {
			wDir_x[nnum] = 0.0;
			wDir_y[nnum] = z/lyz;
			wDir_z[nnum] = -y/lyz;
			nnum++;
			wDir_x[nnum] = -(y*y+z*z)/lyz/l;
			wDir_y[nnum] = x*y/lyz/l;
			wDir_z[nnum] = x*z/lyz/l;
			nnum++;
		}
	}
}

void
getDirections(int numDirec, int numWaves,
	std::vector<double> &wDir_x, std::vector<double> &wDir_y, std::vector<double> &wDir_z)
{
	if(numDirec == 0) return;

/* PJSA: this is JingLi's old method for 13 or less directions and numWaves = 3
  if(numDirec <= 13) {
    wDir_x = new double[3*numDirec];
    wDir_y = new double[3*numDirec];
    wDir_z = new double[3*numDirec];
    getDirections13(numDirec, wDir_x, wDir_y, wDir_z);
    return;
  }
*/

/*  PJSA: I think this is optimum if numDirec == 4, directions are verticies of a regular tetrahedron
  if(numDirec == 4) {
    wDir_x = new double[3*numDirec];
    wDir_y = new double[3*numDirec];
    wDir_z = new double[3*numDirec];
    int nnum = 0;
    getOneDirection(1.0,0.0,0.0,nnum,numWaves,wDir_x,wDir_y,wDir_z);
    getOneDirection(-0.5,sqrt(3.0)/2.0,0.0,nnum,numWaves,wDir_x,wDir_y,wDir_z);
    getOneDirection(-0.5,-sqrt(3.0)/2.0,0.0,nnum,numWaves,wDir_x,wDir_y,wDir_z);
    getOneDirection(0.0,0.0,1.0,nnum,numWaves,wDir_x,wDir_y,wDir_z);
    return;
  }
*/
	int n = 0;
	int fullNumDirec = 0;
	int tmp = 0;  // number of layers excluding middle and top
	while(numDirec > fullNumDirec) {
		n++;
		tmp = (n+1)/2 - 1;
		fullNumDirec = (n+1)*(n+1) + 4*n*tmp + 2*n*((n+1)%2);
	}

	double d = 1.0/n;
	wDir_x.resize(numWaves*fullNumDirec);
	wDir_y.resize(numWaves*fullNumDirec);
	wDir_z.resize(numWaves*fullNumDirec);

	int i,j,k;
	int nnum = 0;

	// middle layer (half only)
	if(n%2 == 0) {
		k=n/2;
		i = n;
		for(j=0;j<n;j++)
			getOneDirection(d,i,j,k,nnum,numWaves,wDir_x,wDir_y,wDir_z);
		j = n;
		for(i=n;i>0;i--)
			getOneDirection(d,i,j,k,nnum,numWaves,wDir_x,wDir_y,wDir_z);
	}
	// top layer
	k=n;
	for(i=0;i<=n;i++)
		for(j=0;j<=n;j++)
			getOneDirection(d,i,j,k,nnum,numWaves,wDir_x,wDir_y,wDir_z);
	// other layers
	for(int l=1; l<=tmp ; ++l) {
		k = n-l;
		i = n;
		for(j=0;j<n;j++)
			getOneDirection(d,i,j,k,nnum,numWaves,wDir_x,wDir_y,wDir_z);
		j = n;
		for(i=n;i>0;i--)
			getOneDirection(d,i,j,k,nnum,numWaves,wDir_x,wDir_y,wDir_z);
		i = 0;
		for(j=n;j>0;j--)
			getOneDirection(d,i,j,k,nnum,numWaves,wDir_x,wDir_y,wDir_z);
		j = 0;
		for(i=0;i<n;i++)
			getOneDirection(d,i,j,k,nnum,numWaves,wDir_x,wDir_y,wDir_z);
	}
}

template<class Scalar>
void
FetiSub<Scalar>::makeQ() {
	switch (getFetiInfo().augment) {
		default:
			break;
		case FetiInfo::Gs: {
			nGrbm = std::min(this->rigidBodyModes->numRBM(), getFetiInfo().nGs);
			int iLen = scomm->lenT(SComm::std);
			interfaceRBMs.resize(nGrbm * iLen);
			rbms.resize(nGrbm * iLen);
			int *boundToCDSA = (int *) dbg_alloca(sizeof(int) * iLen);
			int i;
			for (i = 0; i < iLen; ++i)
				if (scomm->boundDofT(SComm::std, i) >= 0)
					boundToCDSA[i] = ccToC[scomm->boundDofT(SComm::std, i)];
				else
					boundToCDSA[i] = -1;
			int offset = 0;
			if ((this->rigidBodyModes->numRBM() == 6) && (getFetiInfo().rbmType == FetiInfo::rotation))
				offset = 3;
			this->rigidBodyModes->getRBMs(interfaceRBMs.data(), iLen, boundToCDSA, nGrbm, offset);
		}
			break;
		case FetiInfo::WeightedEdges:
		case FetiInfo::Edges: {
			switch (getFetiInfo().rbmType) {
				default:
				case FetiInfo::translation:
				case FetiInfo::rotation:
				case FetiInfo::all:
				case FetiInfo::None:
					if (isMixedSub) {
						this->makeEdgeVectorsPlus(true); // build augmentation for fluid dofs
						this->makeEdgeVectorsPlus(false);  // build augmentation for structure dofs
					} else {
						this->makeEdgeVectorsPlus(isFluidSub,
						                          isThermalSub,
						                          isUndefinedSub);
					}
					break;
				case FetiInfo::averageTran:
				case FetiInfo::averageRot:
				case FetiInfo::averageAll:
					makeAverageEdgeVectors();
					break;
			}
		}
			break;
	}
}

// I've changed this routine to compute the following Q vectors:
//
// averageTran = [Qx + Qy]  i.e. [ 1 1 0 ]^T at a node
// averageRot  = [Qy + Qz]  i.e. [ 0 1 1 ]^T at a node
// averageAll  = computes both of these Q vectors
//
// These 2 vectors are better than [Qx Qy Qz] together, which seems too
// "weak" of a constraint to help the convergence very much. Maybe there
// are other "stronger" vectors that will help more.
//
template<class Scalar>
void
FetiSub<Scalar>::makeAverageEdgeVectors() {
	int i;
	int numR = 1;
	if (getFetiInfo().rbmType == FetiInfo::averageAll)
		numR = 2;
	int totalLengthGrc = 0;

	Connectivity &sharedNodes = *(scomm->sharedNodes);
	edgeDofSize.resize(scomm->numNeighb);

	// 1. first count number of edge dofs
	int iSub, iNode, nE = 0;
	for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		edgeDofSize[iSub] = 0;
		if(std::find_if(boundaryDOFs[iSub].begin(), boundaryDOFs[iSub].end(),
		                [](auto dofs) {
			                return dofs.contains(DofSet::XYZdisp | DofSet::XYZrot);
		                }) != boundaryDOFs[iSub].end())
			edgeDofSize[iSub] = numR;
		if (edgeDofSize[iSub] > 0) nE++;
	}

	std::vector<int> xyzCount(nE * numR, 0);

	nE = 0;
	int nDof = (numR > 1) ? 1 : 0;

	for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		if (edgeDofSize[iSub] == 0) continue;
		int xCount = 0, yCount = 0, zCount = 0, xrCount = 0, yrCount = 0, zrCount = 0;
		for (iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
			if (boundaryDOFs[iSub][iNode].contains(DofSet::Xdisp)) xCount++;
			if (boundaryDOFs[iSub][iNode].contains(DofSet::Ydisp)) yCount++;
			if (boundaryDOFs[iSub][iNode].contains(DofSet::Zdisp)) zCount++;
			if (getFetiInfo().rbmType == FetiInfo::averageRot || (numR > 1)) {
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Xrot)) xrCount++;
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Yrot)) yrCount++;
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Zrot)) zrCount++;
			}
		}
		int edgeLength = 0;
		if (getFetiInfo().rbmType == FetiInfo::averageTran || (numR > 1)) {
			xyzCount[numR * nE + 0] = xCount + yCount + zCount;
			edgeLength += xCount + yCount + zCount;
		}
		if (getFetiInfo().rbmType == FetiInfo::averageRot || (numR > 1)) {
			// KHP
			//xyzCount[numR*nE+nDof]  = xrCount + yrCount + zrCount;
			//edgeLength             += xrCount + yrCount + zrCount;
			xyzCount[numR * nE + nDof] = xCount + yCount + zCount;
			edgeLength += xCount + yCount + zCount;
		}
		totalLengthGrc += edgeLength;

		nE++;
	}

// --------------------------------------------------------
	std::vector<int> xyzList(totalLengthGrc);
	std::vector<Scalar> xyzCoefs(totalLengthGrc);
	int xOffset = 0;
	int xrOffset = 0;
	nE = 0;
	for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		if (edgeDofSize[iSub] == 0) continue;
		if (numR > 1) {
			xrOffset = xOffset + xyzCount[numR * nE + 0];
		}
		Scalar sign = (scomm->subNums[iSub] < subNum()) ? 1.0 : -1.0;
		int used = 0;
		int middleNode = sharedNodes.num(iSub) / 2;
		for (iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
			if (boundaryDOFs[iSub][iNode].contains(DofSet::Xdisp)) {
				int xDof = cc_dsa->locate(sharedNodes[iSub][iNode], DofSet::Xdisp);
				xyzList[xOffset] = xDof;
				xyzCoefs[xOffset++] = (middleNode == iNode) ? sign * 1.0 : 0.0;
				used = 1;
			}
			if (boundaryDOFs[iSub][iNode].contains(DofSet::Ydisp)) {
				int yDof = cc_dsa->locate(sharedNodes[iSub][iNode], DofSet::Ydisp);
				xyzList[xOffset] = yDof;
				xyzCoefs[xOffset++] = (middleNode == iNode) ? sign * 0.0 : 0.0;
				used = 1;
			}
			if (boundaryDOFs[iSub][iNode].contains(DofSet::Zdisp)) {
				int zDof = cc_dsa->locate(sharedNodes[iSub][iNode], DofSet::Zdisp);
				xyzList[xOffset] = zDof;
				xyzCoefs[xOffset++] = (middleNode == iNode) ? sign * 1.0 : 0.0;
				used = 1;
			}

			if (getFetiInfo().rbmType == FetiInfo::averageRot || (numR > 1)) {
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Xdisp)) {
					int xDof = cc_dsa->locate(sharedNodes[iSub][iNode], DofSet::Xdisp);
					xyzList[xrOffset] = xDof;
					xyzCoefs[xrOffset++] = (middleNode == iNode) ? sign * 1.0 : 0.0;
					used = 1;
				}
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Ydisp)) {
					int yDof = cc_dsa->locate(sharedNodes[iSub][iNode], DofSet::Ydisp);
					xyzList[xrOffset] = yDof;
					xyzCoefs[xrOffset++] = (middleNode == iNode) ? sign * 1.0 : 0.0;
					used = 1;
				}
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Zdisp)) {
					int zDof = cc_dsa->locate(sharedNodes[iSub][iNode], DofSet::Zdisp);
					xyzList[xrOffset] = zDof;
					xyzCoefs[xrOffset++] = (middleNode == iNode) ? sign * 0.0 : 0.0;
					used = 1;
				}
			}

		}
		int off = 0;
		if (numR > 1)
			off += xyzCount[numR * nE + 1];

		xOffset += off;

		if (used) nE++;
	}
	this->Grc = std::make_unique<GenCuCSparse<Scalar>>(numR * nE, cc_dsa->size(),
	                                                   xyzCount, xyzList, std::move(xyzCoefs));
	// Src->setSparseMatrices(1, Grc);
	this->Src->addSparseMatrix(this->Grc);
}

template<class Scalar>
void
FetiSub<Scalar>::makeEdgeVectorsPlus(bool isFluidSub, bool isThermalSub,
                                          bool isUndefinedSub) {
	int i, iSub, iNode;
	//bool isCoupled = isCoupled;
	int spaceDim = getFetiInfo().spaceDimension;
	// Number of directions in the coarse problem, choose from 0, 1, ..., 13
	int numDirec = getFetiInfo().numdir;
	int numWaves = spaceDim;      // Number of long and trans waves
	Connectivity &sharedNodes = *(scomm->sharedNodes);

	if (isUndefinedSub) {
		int nhelm = 0, ntemp = 0, ndisp = 0;
		for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
			if (scomm->subNums[iSub] == subNum()) continue;
			for (iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Helm)) nhelm++;
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Temp)) ntemp++;
				if (boundaryDOFs[iSub][iNode].contains(DofSet::XYZdisp)) ndisp++;
			}
		}
		isFluidSub = nhelm > 0;
		isThermalSub = ntemp > 0;
		if (nhelm * ntemp > 0 || nhelm * ndisp || ntemp * ndisp)
			std::cerr << " *** WARNING: multiple types of dofs with undefined type\n";
	}

	if (isFluidSub || isThermalSub)
		numWaves = 1;
	int numCS = 2;  // Choose cos or/and sin mode

	// If no direction added, put numCS = 0, and numWaves = 0.
	if (numDirec == 0) {
		numWaves = 0;
		numCS = 0;
	}

	// To store the Directions of the long and transverse waves
	std::vector<double> d_x(numDirec);
	std::vector<double> d_y(numDirec);
	std::vector<double> d_z(numDirec);
	std::vector<double> t_x(numDirec);
	std::vector<double> t_y(numDirec);
	std::vector<double> t_z(numDirec);

	std::vector<double> wDir_x;
	std::vector<double> wDir_y;
	std::vector<double> wDir_z;

	double pi = 3.141592653589793;
	if (spaceDim == 2) {
		for (int i = 0; i < numDirec; i++) {
			d_x[i] = cos((pi * i) / (numDirec));
			d_y[i] = sin((pi * i) / (numDirec));
			d_z[i] = 0.0;
			t_x[i] = -d_y[i];
			t_y[i] = d_x[i];
			t_z[i] = 0.0;
		}
	} else {
		getDirections(numDirec, numWaves, wDir_x, wDir_y, wDir_z);
		for (i = 0; i < numDirec; i++) {
			d_x[i] = wDir_x[numWaves * i + 0];
			d_y[i] = wDir_y[numWaves * i + 0];
			d_z[i] = wDir_z[numWaves * i + 0];
		}
	}

	int numR = getFetiInfo().nGs;
	int totalLengthGrc = 0;
	DofSet desired;
	if (isFluidSub) {
		if (numR > 0) numR = 1;
		desired = DofSet::Helm;
	} else if (isThermalSub) {
		if (numR > 0) numR = 1;
		desired = DofSet::Temp;
	} else {
		if (getFetiInfo().rbmType == FetiInfo::rotation || (numR > 3)) {
			numR = 6;
			desired = DofSet::XYZdisp | DofSet::XYZrot;
		} else desired = DofSet::XYZdisp;
	}

	// edgeDofSize: number of augmentation degree of freedom per edge
	if (edgeDofSize.size() == 0) {
		edgeDofSize.resize(scomm->numNeighb);
		for (i = 0; i < scomm->numNeighb; ++i) edgeDofSize[i] = 0;
	}
	edgeDofSizeTmp.resize(scomm->numNeighb);
	// total: total number of augmentaions per sub
	int total = 0;

	int numdofperNode = 0;
	if (isFluidSub || isThermalSub)
		numdofperNode = 1;
	else if (getFetiInfo().waveType == FetiInfo::solid)
		numdofperNode = 3;
	else numdofperNode = 6;

	int numInterfNodes = sharedNodes.numConnect();
	int nQPerNeighb = numR + numDirec * numWaves * numCS;
	int lQ = nQPerNeighb * numdofperNode * numInterfNodes;
	std::vector<double> Q(lQ, 0.0);
	bool *isUsed = (bool *) dbg_alloca(sizeof(bool) * (nQPerNeighb * scomm->numNeighb));
	for (i = 0; i < nQPerNeighb * scomm->numNeighb; ++i)
		isUsed[i] = false;

	// 1. first count number of edge dofs
	edgeDofs.resize(scomm->numNeighb);
	for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		edgeDofSizeTmp[iSub] = 0;
		if (scomm->subNums[iSub] == subNum()) continue;
		int nhelm = 0, ntemp = 0;
		int nx = 0, ny = 0, nz = 0;
		int nxr = 0, nyr = 0, nzr = 0;
		int label = 0;
		for ( auto nodeDOFs : boundaryDOFs[iSub]) {
			if ((isFluidSub && !nodeDOFs.contains(DofSet::Helm)) ||
			    (!isFluidSub && nodeDOFs.contains(DofSet::Helm)))
				continue;
			if ((isThermalSub && !nodeDOFs.contains(DofSet::Temp)) ||
			    (!isThermalSub && nodeDOFs.contains(DofSet::Temp)))
				continue;

			if (nodeDOFs.contains(DofSet::Helm))
				nhelm++;
			if (nodeDOFs.contains(DofSet::Temp))
				ntemp++;
			if (nodeDOFs.contains(DofSet::Xdisp))
				nx++;
			if (nodeDOFs.contains(DofSet::Ydisp))
				ny++;
			if (nodeDOFs.contains(DofSet::Zdisp))
				nz++;
			if (nodeDOFs.contains(DofSet::Xrot))
				nxr++;
			if (nodeDOFs.contains(DofSet::Yrot))
				nyr++;
			if (nodeDOFs.contains(DofSet::Zrot))
				nzr++;

			if (nodeDOFs.containsAllDisp(spaceDim)) label++;

			edgeDofs[iSub] |= nodeDOFs & desired;
		}

		// Check if we should add rotation for problems that do not
		// define rotational degrees of freedom
		if (desired.contains(DofSet::XYZrot)) {
			// if ny +nz is bigger than 2 (at least 3) there is no danger in putting a Xrot
			if (ny + nz > 2) edgeDofs[iSub] |= DofSet::Xrot;
			// if nx+nz > 2 and there are some y AND some z, we add the Yrot
			if (nx + nz > 2 && (ny * nx) != 0) edgeDofs[iSub] |= DofSet::Yrot;
			if (nx + ny > 2 && (nx * ny * nz) != 0) edgeDofs[iSub] |= DofSet::Zrot;
		}

		if (numR > 0) {
			if (nhelm > 0)
				isUsed[iSub * nQPerNeighb] = true;
			if (ntemp > 0)
				isUsed[iSub * nQPerNeighb] = true;
			if (nx > 0)
				isUsed[0 + iSub * nQPerNeighb] = true;
			if (ny > 0)
				isUsed[1 + iSub * nQPerNeighb] = true;
			if (nz > 0)
				isUsed[2 + iSub * nQPerNeighb] = true;
			if (desired.contains(DofSet::XYZrot)) {
				if ((nxr > 0) || (ny + nz > 2))
					isUsed[3 + iSub * nQPerNeighb] = true;
				if ((nyr > 0) || (nx + nz > 2))
					isUsed[4 + iSub * nQPerNeighb] = true;
				if ((nzr > 0) || (nx + ny > 2))
					isUsed[5 + iSub * nQPerNeighb] = true;
			}
		}

		if ((label > 0) || (nhelm > 0))
			for (i = 0; i < numDirec * numWaves * numCS; i++)
				isUsed[numR + i + iSub * nQPerNeighb] = true;
		if ((label == 0) && (nhelm == 0))
			edgeDofSizeTmp[iSub] = edgeDofs[iSub].count();
		else {
			if (numR > 0)
				edgeDofSizeTmp[iSub] = edgeDofs[iSub].count() + numDirec * numWaves * numCS;
			else
				edgeDofSizeTmp[iSub] = numDirec * numWaves * numCS;
		}
		total += edgeDofSizeTmp[iSub];
	}

	// Form Q matrix
	for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		if (scomm->subNums[iSub] == subNum()) continue;
		int sOffset = sharedNodes.offset(iSub); // Offset of the neighbor in the boundary nodes.

		// compute the wave numbers for EdgeWs augmentation
		double k_pSolid = 0.0, k_sSolid = 0.0, k_pShell = 0.0, k_sShell = 0.0, k_pFluid = 0.0;
		//if(!isFluidSub && (numDirec > 0)) {
		if (numDirec > 0) { // support for multiple fluids
			if (getFetiInfo().waveMethod == FetiInfo::uniform) {
				k_pSolid = k_pShell = k_p;
				k_sSolid = k_sShell = k_s;
				k_pFluid = k_f;
			}
			if (getFetiInfo().waveMethod == FetiInfo::averageK) {
				k_pSolid = k_pShell = neighbK_p[iSub];
				k_sSolid = k_sShell = neighbK_s[iSub];
				k_pFluid = neighbK_f[iSub];
			} else if (getFetiInfo().waveMethod == FetiInfo::averageMat) {
				double omega2 = getShiftVal();
				double lambda = (neighbPrat[iSub] * neighbYmod[iSub]) /
				                (1.0 + neighbPrat[iSub]) / (1.0 - 2.0 * neighbPrat[iSub]);
				double mu = neighbYmod[iSub] / 2.0 / (1.0 + neighbPrat[iSub]);
				k_pSolid = sqrt(omega2 * neighbDens[iSub]) / sqrt(lambda + 2 * mu);
				k_sSolid = sqrt(omega2 * neighbDens[iSub]) / sqrt(mu);
				k_pFluid = sqrt(omega2) / neighbSspe[iSub]; // should use bulk modulus and density?
				if (neighbThih[iSub] > 0.0) {
					double di = neighbYmod[iSub] * neighbThih[iSub] * neighbThih[iSub] * neighbThih[iSub] /
					            (12.0 * (1.0 - neighbPrat[iSub] * neighbPrat[iSub]));
					double beta4 = omega2 * neighbDens[iSub] * neighbThih[iSub] / di;
					k_pShell = k_sShell = sqrt(sqrt(beta4));
				} else {
					k_pShell = k_sShell = 0.0;
				}
			}
		}

		auto &nodes = getNodeSet();
		// Find the center of the edge/face for EdgeGs augmentation
		double xc = 0, yc = 0, zc = 0;
		if (numR > 0) {
			for ( auto nd : sharedNodes[iSub]) {
				xc += nodes[ nd ]->x;
				yc += nodes[ nd ]->y;
				zc += nodes[ nd ]->z;
			}
			xc /= sharedNodes.num(iSub);
			yc /= sharedNodes.num(iSub);
			zc /= sharedNodes.num(iSub);
		}

		// loop over all the nodes on the interface and fill Q
		for (iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
			if ((isFluidSub && !boundaryDOFs[iSub][iNode].contains(DofSet::Helm)) ||
			    (!isFluidSub && boundaryDOFs[iSub][iNode].contains(DofSet::Helm)))
				continue;
			if ((isThermalSub && !boundaryDOFs[iSub][iNode].contains(DofSet::Temp)) ||
			    (!isThermalSub && boundaryDOFs[iSub][iNode].contains(DofSet::Temp)))
				continue;

			double x = nodes[sharedNodes[iSub][iNode]]->x;
			double y = nodes[sharedNodes[iSub][iNode]]->y;
			double z = nodes[sharedNodes[iSub][iNode]]->z;

			if (numR > 0) {
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Helm))
					Q[numdofperNode * (iNode + sOffset) + 0] = 1.0;
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Temp))
					Q[numdofperNode * (iNode + sOffset) + 0] = 1.0;
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Xdisp)) {
					Q[numdofperNode * (iNode + sOffset) + 0] = 1.0;
					if (edgeDofs[iSub].contains(DofSet::Yrot)) {
						Q[numdofperNode * (iNode + sOffset + 4 * numInterfNodes) + 0] = (z - zc);
					}
					if (edgeDofs[iSub].contains(DofSet::Zrot)) {
						Q[numdofperNode * (iNode + sOffset + 5 * numInterfNodes) + 0] = -(y - yc);
					}
				}
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Ydisp)) {
					Q[numdofperNode * (iNode + sOffset + numInterfNodes) + 1] = 1.0;
					if (edgeDofs[iSub].contains(DofSet::Xrot)) {
						Q[numdofperNode * (iNode + sOffset + 3 * numInterfNodes) + 1] = -(z - zc);
					}
					if (edgeDofs[iSub].contains(DofSet::Zrot)) {
						Q[numdofperNode * (iNode + sOffset + 5 * numInterfNodes) + 1] = (x - xc);
					}
				}
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Zdisp)) {
					Q[numdofperNode * (iNode + sOffset + 2 * numInterfNodes) + 2] = 1.0;
					if (edgeDofs[iSub].contains(DofSet::Xrot)) {
						Q[numdofperNode * (iNode + sOffset + 3 * numInterfNodes) + 2] = (y - yc);
					}
					if (edgeDofs[iSub].contains(DofSet::Yrot)) {
						Q[numdofperNode * (iNode + sOffset + 4 * numInterfNodes) + 2] = -(x - xc);
					}
				}
				if (numdofperNode == 6) {
					if ((boundaryDOFs[iSub][iNode].contains(DofSet::Xrot)) && edgeDofs[iSub].contains(DofSet::Xrot))
						Q[numdofperNode * (iNode + sOffset + 3 * numInterfNodes) + 3] = 1.0;
					if ((boundaryDOFs[iSub][iNode].contains(DofSet::Yrot)) && edgeDofs[iSub].contains(DofSet::Yrot))
						Q[numdofperNode * (iNode + sOffset + 4 * numInterfNodes) + 4] = 1.0;
					if ((boundaryDOFs[iSub][iNode].contains(DofSet::Zrot)) && edgeDofs[iSub].contains(DofSet::Zrot))
						Q[numdofperNode * (iNode + sOffset + 5 * numInterfNodes) + 5] = 1.0;
				}
			}

			if (boundaryDOFs[iSub][iNode].containsAllDisp(spaceDim) ||
			    boundaryDOFs[iSub][iNode].contains(DofSet::Helm) ||
			    boundaryDOFs[iSub][iNode].contains(DofSet::Temp))
				for (int iDir = 0; iDir < numDirec; iDir++) {
					double k_p, k_s;
					if (boundaryDOFs[iSub][iNode].containsAnyRot()) {
						k_p = k_pShell;
						k_s = k_sShell;
					} else {
						k_p = k_pSolid;
						k_s = k_sSolid;
					}
					double cosValP, cosValS, sinValP, sinValS, cosValH, sinValH;
					double ddotx = d_x[iDir] * x + d_y[iDir] * y + d_z[iDir] * z;
					cosValP = cos(k_p * ddotx);
					cosValS = cos(k_s * ddotx);
					sinValP = sin(k_p * ddotx);
					sinValS = sin(k_s * ddotx);
					sinValH = sin(k_pFluid * ddotx);
					cosValH = cos(k_pFluid * ddotx);
					double cosVal;
					double sinVal;
					for (int iW = 0; iW < numWaves; iW++) {
						if (iW == 0) {
							cosVal = cosValP;
							sinVal = sinValP;
						} else {
							cosVal = cosValS;
							sinVal = sinValS;
						}
						for (int iCS = 0; iCS < numCS; iCS++) {
							int waveOffset =
									(numR + iDir * numWaves * numCS + numCS * iW + iCS) * numInterfNodes + sOffset;
							if (iCS == 0) {
								if (numdofperNode == 1) Q[numdofperNode * (waveOffset + iNode) + 0] = cosValH;
								else {
									if (spaceDim == 3) {
										Q[numdofperNode * (waveOffset + iNode) + 0] =
												cosVal * wDir_x[iDir * numWaves + iW];
										Q[numdofperNode * (waveOffset + iNode) + 1] =
												cosVal * wDir_y[iDir * numWaves + iW];
										Q[numdofperNode * (waveOffset + iNode) + 2] =
												cosVal * wDir_z[iDir * numWaves + iW];
									} else {
										if (iW == 0) {
											Q[numdofperNode * (waveOffset + iNode) + 0] = cosVal * d_x[iDir];
											Q[numdofperNode * (waveOffset + iNode) + 1] = cosVal * d_y[iDir];
										} else {
											Q[numdofperNode * (waveOffset + iNode) + 0] = cosVal * t_x[iDir];
											Q[numdofperNode * (waveOffset + iNode) + 1] = cosVal * t_y[iDir];
										}
									}
								}
							}
							if (iCS == 1) {
								if (numdofperNode == 1) Q[numdofperNode * (waveOffset + iNode) + 0] = sinValH;
								else {
									if (spaceDim == 3) {
										Q[numdofperNode * (waveOffset + iNode) + 0] =
												sinVal * wDir_x[iDir * numWaves + iW];
										Q[numdofperNode * (waveOffset + iNode) + 1] =
												sinVal * wDir_y[iDir * numWaves + iW];
										Q[numdofperNode * (waveOffset + iNode) + 2] =
												sinVal * wDir_z[iDir * numWaves + iW];
									} else {
										if (iW == 0) {
											Q[numdofperNode * (waveOffset + iNode) + 0] = sinVal * d_x[iDir];
											Q[numdofperNode * (waveOffset + iNode) + 1] = sinVal * d_y[iDir];
										} else {
											Q[numdofperNode * (waveOffset + iNode) + 0] = sinVal * t_x[iDir];
											Q[numdofperNode * (waveOffset + iNode) + 1] = sinVal * t_y[iDir];
										}
									}
								}
							}
						}
					}
				}
		}
	}

	if (numdofperNode == 6)
		GramSchmidt(Q.data(), isUsed, DofSet::XYZdisp | DofSet::XYZrot, nQPerNeighb,
		            getFetiInfo().augmentimpl == FetiInfo::Primal);
	else
		GramSchmidt(Q.data(), isUsed, desired, nQPerNeighb, getFetiInfo().augmentimpl == FetiInfo::Primal);

	int *nQAdd = (int *) dbg_alloca(sizeof(int) * (scomm->numNeighb));
	int *nQAddWaves = (int *) dbg_alloca(sizeof(int) * (scomm->numNeighb));

	for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		nQAdd[iSub] = 0;
		nQAddWaves[iSub] = 0;
	}

	// Count the number of Qs that each neighbor will have, and reset total, edgeDofSizeTmp[]
	total = 0;
	for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		if (scomm->subNums[iSub] == subNum()) continue;
		for (int iQ = 0; iQ < nQPerNeighb; ++iQ)
			if (isUsed[iQ + iSub * nQPerNeighb]) {
				nQAdd[iSub]++;
				if (iQ >= numR)
					nQAddWaves[iSub]++;
			}
		edgeDofSizeTmp[iSub] = nQAdd[iSub];
		total += edgeDofSizeTmp[iSub];
	}

	std::vector<int> HelmCount(total,0);
	std::vector<int> TempCount(total,0);
	std::vector<int> xyzCount(total,0);

	int oldTot = total;

	total = 0;
	for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		if (edgeDofSizeTmp[iSub] == 0) continue;
		int hCount = 0, tCount = 0, xCount = 0, yCount = 0, zCount = 0, xrCount = 0, yrCount = 0, zrCount = 0;
		int waveCount = 0;
		for (iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
			if ((isFluidSub && !boundaryDOFs[iSub][iNode].contains(DofSet::Helm)) ||
			    (!isFluidSub && boundaryDOFs[iSub][iNode].contains(DofSet::Helm)))
				continue;
			if ((isThermalSub && !boundaryDOFs[iSub][iNode].contains(DofSet::Temp)) ||
			    (!isThermalSub && boundaryDOFs[iSub][iNode].contains(DofSet::Temp)))
				continue;

			if (boundaryDOFs[iSub][iNode].contains(DofSet::Helm)) hCount++;
			if (boundaryDOFs[iSub][iNode].contains(DofSet::Temp)) tCount++;
			if (numR > 0) {
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Xdisp)) xCount++;
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Ydisp)) yCount++;
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Zdisp)) zCount++;
				if (numdofperNode == 6) {
					if (boundaryDOFs[iSub][iNode].contains(DofSet::Xrot)) xrCount++;
					if (boundaryDOFs[iSub][iNode].contains(DofSet::Yrot)) yrCount++;
					if (boundaryDOFs[iSub][iNode].contains(DofSet::Zrot)) zrCount++;
				}
			}
			if (boundaryDOFs[iSub][iNode].containsAllDisp(spaceDim)) waveCount++;
		}

		int edgeLength = 0;
		if (numR > 0) {
			if (desired.contains(DofSet::Helm))
				if (edgeDofs[iSub].contains(DofSet::Helm) && isUsed[0 + iSub * nQPerNeighb])
					edgeLength += HelmCount[total++] = hCount;
			if (desired.contains(DofSet::Temp))
				if (edgeDofs[iSub].contains(DofSet::Temp) && isUsed[0 + iSub * nQPerNeighb])
					edgeLength += TempCount[total++] = tCount;
			if (desired.contains(DofSet::XYZdisp)) {
				if (edgeDofs[iSub].contains(DofSet::Xdisp) && isUsed[0 + iSub * nQPerNeighb])
					edgeLength += xyzCount[total++] = xCount;
				if (edgeDofs[iSub].contains(DofSet::Ydisp) && isUsed[1 + iSub * nQPerNeighb])
					edgeLength += xyzCount[total++] = yCount;
				if (edgeDofs[iSub].contains(DofSet::Zdisp) && isUsed[2 + iSub * nQPerNeighb])
					edgeLength += xyzCount[total++] = zCount;
			}
			if (desired.contains(DofSet::XYZrot)) {
				if (edgeDofs[iSub].contains(DofSet::Xrot) && isUsed[3 + iSub * nQPerNeighb])
					edgeLength += xyzCount[total++] = xrCount + yCount + zCount;
				if (edgeDofs[iSub].contains(DofSet::Yrot) && isUsed[4 + iSub * nQPerNeighb])
					edgeLength += xyzCount[total++] = yrCount + xCount + zCount;
				if (edgeDofs[iSub].contains(DofSet::Zrot) && isUsed[5 + iSub * nQPerNeighb])
					edgeLength += xyzCount[total++] = zrCount + xCount + yCount;
			}
		}
		if (desired.contains(DofSet::XYZdisp)) {
			if (waveCount > 0) {
				if (edgeDofs[iSub].containsAllDisp(spaceDim)) {
					for (i = 0; i < numDirec * numWaves * numCS; i++)
						if (isUsed[numR + i + iSub * nQPerNeighb])
							edgeLength += xyzCount[total++] = numWaves * waveCount;
				}
			}
		}
		if (desired.contains(DofSet::Helm))
			if ((hCount > 0) && (edgeDofs[iSub].contains(DofSet::Helm)))
				for (i = 0; i < numDirec * numWaves * numCS; i++)
					if (isUsed[numR + i + iSub * nQPerNeighb])
						edgeLength += HelmCount[total++] = 1 * hCount;
		if (desired.contains(DofSet::Temp))
			if ((tCount > 0) && (edgeDofs[iSub].contains(DofSet::Temp)))
				for (i = 0; i < numDirec * numWaves * numCS; i++)
					if (isUsed[numR + i + iSub * nQPerNeighb])
						edgeLength += TempCount[total++] = 1 * tCount;

		totalLengthGrc += edgeLength;
	}
	if (total != oldTot) fprintf(stderr, "Non match %d %d\n", total, oldTot);

	std::vector<int> HelmList(totalLengthGrc, 0);
	std::vector<Scalar> HelmCoefs(totalLengthGrc, 0.0);
	std::vector<int> TempList(totalLengthGrc, 0);
	std::vector<Scalar> TempCoefs(totalLengthGrc, 0.0);
	std::vector<int> xyzList(totalLengthGrc, 0);
	std::vector<Scalar> xyzCoefs(totalLengthGrc, 0.0);

	int hOffset = 0, tOffset = 0, waveOffset = 0;
	int xOffset = 0, yOffset = 0, zOffset = 0, xrOffset = 0, yrOffset = 0, zrOffset = 0;
	total = 0;
	for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		if (edgeDofSizeTmp[iSub] == 0) continue;

		int index = 0;
		int off;
		if ((desired.contains(DofSet::XYZdisp)) || (desired.contains(DofSet::Helm)) ||
		    (desired.contains(DofSet::Temp))) {
			if (numR > 0) {
				yOffset = xOffset +
				          ((edgeDofs[iSub].contains(DofSet::Xdisp) && isUsed[0 + iSub * nQPerNeighb])
				           ? xyzCount[total++] : 0);
				zOffset = yOffset +
				          ((edgeDofs[iSub].contains(DofSet::Ydisp) && isUsed[1 + iSub * nQPerNeighb])
				           ? xyzCount[total++] : 0);
				xrOffset = zOffset +
				           ((edgeDofs[iSub].contains(DofSet::Zdisp) && isUsed[2 + iSub * nQPerNeighb])
				            ? xyzCount[total++] : 0);
				yrOffset = xrOffset +
				           ((edgeDofs[iSub].contains(DofSet::Xrot) && isUsed[3 + iSub * nQPerNeighb])
				            ? xyzCount[total++] : 0);
				zrOffset = yrOffset +
				           ((edgeDofs[iSub].contains(DofSet::Yrot) && isUsed[4 + iSub * nQPerNeighb])
				            ? xyzCount[total++] : 0);
				if (numdofperNode == 1) {
					if (isFluidSub)
						waveOffset = hOffset +
						             ((edgeDofs[iSub].contains(DofSet::Helm) && isUsed[0 + iSub * nQPerNeighb])
						              ? HelmCount[total++] : 0);
					else
						waveOffset = tOffset +
						             ((edgeDofs[iSub].contains(DofSet::Temp) && isUsed[0 + iSub * nQPerNeighb])
						              ? TempCount[total++] : 0);
				} else {
					waveOffset = zrOffset +
					             ((edgeDofs[iSub].contains(DofSet::Zrot) && isUsed[5 + iSub * nQPerNeighb])
					              ? xyzCount[total++] : 0);
				}
			}
			if (nQAddWaves[iSub] > 0) {
				if (numdofperNode == 1)
					if (isFluidSub)
						index = HelmCount[(total += nQAddWaves[iSub]) - 1];
					else
						index = TempCount[(total += nQAddWaves[iSub]) - 1];
				else
					index = xyzCount[(total += nQAddWaves[iSub]) - 1];
			} else index = 0;
			off = waveOffset + ((nQAddWaves[iSub] > 0) ? nQAddWaves[iSub] * index : 0);
		}

		double sign = 1.0;
		if (getFetiInfo().augmentimpl != FetiInfo::Primal) // 012314 JAT
			sign = (scomm->subNums[iSub] < subNum()) ? 1.0 : -1.0;
		int sOffset = sharedNodes.offset(iSub);
		for (iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
			if ((isFluidSub && !boundaryDOFs[iSub][iNode].contains(DofSet::Helm)) ||
			    (!isFluidSub && boundaryDOFs[iSub][iNode].contains(DofSet::Helm)))
				continue;
			if ((isThermalSub && !boundaryDOFs[iSub][iNode].contains(DofSet::Temp)) ||
			    (!isThermalSub && boundaryDOFs[iSub][iNode].contains(DofSet::Temp)))
				continue;

			int qOff = sOffset + iNode;
			int hDof = -1, tDof = -1, xDof = -1, yDof = -1, zDof = -1;
			if (boundaryDOFs[iSub][iNode].contains(DofSet::Helm)) {
				hDof = cc_dsa->locate(sharedNodes[iSub][iNode], DofSet::Helm);
				if (numR > 0) {
					if (edgeDofs[iSub].contains(DofSet::Helm) && isUsed[0 + iSub * nQPerNeighb]) {
						HelmList[hOffset] = hDof;
						HelmCoefs[hOffset++] = sign * Q[numdofperNode * qOff + 0];
					}
				}
			}
			if (boundaryDOFs[iSub][iNode].contains(DofSet::Temp)) {
				tDof = cc_dsa->locate(sharedNodes[iSub][iNode], DofSet::Temp);
				if (numR > 0) {
					if (edgeDofs[iSub].contains(DofSet::Temp) && isUsed[0 + iSub * nQPerNeighb]) {
						TempList[tOffset] = tDof;
						TempCoefs[tOffset++] = sign * Q[numdofperNode * qOff + 0];
					}
				}
			}
			if (boundaryDOFs[iSub][iNode].contains(DofSet::Xdisp)) {
				xDof = cc_dsa->locate(sharedNodes[iSub][iNode], DofSet::Xdisp);
				if (numR > 0) {
					if (edgeDofs[iSub].contains(DofSet::Xdisp) && isUsed[0 + iSub * nQPerNeighb]) {
						xyzList[xOffset] = xDof;
						xyzCoefs[xOffset++] = sign * Q[numdofperNode * qOff];
					}
					if (edgeDofs[iSub].contains(DofSet::Yrot) && isUsed[4 + iSub * nQPerNeighb]) {
						xyzList[yrOffset] = xDof;
						xyzCoefs[yrOffset++] = sign * Q[numdofperNode * (qOff + 4 * numInterfNodes)];
					}
					if (edgeDofs[iSub].contains(DofSet::Zrot) && isUsed[5 + iSub * nQPerNeighb]) {
						xyzList[zrOffset] = xDof;
						xyzCoefs[zrOffset++] = sign * Q[numdofperNode * (qOff + 5 * numInterfNodes)];
					}
				}
			}
			if (boundaryDOFs[iSub][iNode].contains(DofSet::Ydisp)) {
				yDof = cc_dsa->locate(sharedNodes[iSub][iNode], DofSet::Ydisp);
				if (numR > 0) {
					if (edgeDofs[iSub].contains(DofSet::Ydisp) && isUsed[1 + iSub * nQPerNeighb]) {
						xyzList[yOffset] = yDof;
						xyzCoefs[yOffset++] = sign * Q[numdofperNode * (qOff + numInterfNodes) + 1];
					}
					if (edgeDofs[iSub].contains(DofSet::Xrot) && isUsed[3 + iSub * nQPerNeighb]) {
						xyzList[xrOffset] = yDof;
						xyzCoefs[xrOffset++] = sign * Q[numdofperNode * (qOff + 3 * numInterfNodes) + 1];
					}
					if (edgeDofs[iSub].contains(DofSet::Zrot) && isUsed[5 + iSub * nQPerNeighb]) {
						xyzList[zrOffset] = yDof;
						xyzCoefs[zrOffset++] = sign * Q[numdofperNode * (qOff + 5 * numInterfNodes) + 1];
					}
				}
			}
			if (boundaryDOFs[iSub][iNode].contains(DofSet::Zdisp)) {
				zDof = cc_dsa->locate(sharedNodes[iSub][iNode], DofSet::Zdisp);
				if (numR > 0) {
					if (edgeDofs[iSub].contains(DofSet::Zdisp) && isUsed[2 + iSub * nQPerNeighb]) {
						xyzList[zOffset] = zDof;
						xyzCoefs[zOffset++] = sign * Q[numdofperNode * (qOff + 2 * numInterfNodes) + 2];
					}
					if (edgeDofs[iSub].contains(DofSet::Xrot) && isUsed[3 + iSub * nQPerNeighb]) {
						xyzList[xrOffset] = zDof;
						xyzCoefs[xrOffset++] = sign * Q[numdofperNode * (qOff + 3 * numInterfNodes) + 2];
					}
					if (edgeDofs[iSub].contains(DofSet::Yrot) && isUsed[4 + iSub * nQPerNeighb]) {
						xyzList[yrOffset] = zDof;
						xyzCoefs[yrOffset++] = sign * Q[numdofperNode * (qOff + 4 * numInterfNodes) + 2];
					}
				}
			}
			if ((numR > 3) && (numdofperNode == 6)) {
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Xrot)) {
					int xrDof = cc_dsa->locate(sharedNodes[iSub][iNode], DofSet::Xrot);
					if (edgeDofs[iSub].contains(DofSet::Xrot) && isUsed[3 + iSub * nQPerNeighb]) {
						xyzList[xrOffset] = xrDof;
						xyzCoefs[xrOffset++] = sign * Q[numdofperNode * (qOff + 3 * numInterfNodes) + 3];
					}
				}
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Yrot)) {
					int yrDof = cc_dsa->locate(sharedNodes[iSub][iNode], DofSet::Yrot);
					if (edgeDofs[iSub].contains(DofSet::Yrot) && isUsed[4 + iSub * nQPerNeighb]) {
						xyzList[yrOffset] = yrDof;
						xyzCoefs[yrOffset++] = sign * Q[numdofperNode * (qOff + 4 * numInterfNodes) + 4];
					}
				}
				if (boundaryDOFs[iSub][iNode].contains(DofSet::Zrot)) {
					int zrDof = cc_dsa->locate(sharedNodes[iSub][iNode], DofSet::Zrot);
					if (edgeDofs[iSub].contains(DofSet::Zrot) && isUsed[5 + iSub * nQPerNeighb]) {
						xyzList[zrOffset] = zrDof;
						xyzCoefs[zrOffset++] = sign * Q[numdofperNode * (qOff + 5 * numInterfNodes) + 5];
					}
				}
			}
			if (nQAddWaves[iSub] > 0) {
				if (boundaryDOFs[iSub][iNode].containsAllDisp(spaceDim) ||
				    boundaryDOFs[iSub][iNode].contains(DofSet::Helm) ||
				    boundaryDOFs[iSub][iNode].contains(DofSet::Temp)) {
					int UsedWaves = 0;
					for (int iDir = 0; iDir < numDirec; iDir++)
						for (int iW = 0; iW < numWaves; iW++)
							for (int iCS = 0; iCS < numCS; iCS++) {
								if (isUsed[(numR + iDir * numWaves * numCS + iW * numCS + iCS) + iSub * nQPerNeighb] ==
								    false)
									continue;
								qOff = (numR + iDir * numWaves * numCS + iW * numCS + iCS) * numInterfNodes + sOffset +
								       iNode;
								if (UsedWaves == 0) {
									if (numdofperNode == 1) {
										if (hDof > -1) {
											HelmList[waveOffset] = hDof;
											HelmCoefs[waveOffset++] = sign * Q[numdofperNode * qOff + 0];
										}
										if (tDof > -1) {
											TempList[waveOffset] = tDof;
											TempCoefs[waveOffset++] = sign * Q[numdofperNode * qOff + 0];
										}
									} else {
										if (xDof > -1) {
											xyzList[waveOffset] = xDof;
											xyzCoefs[waveOffset++] = sign * Q[numdofperNode * qOff + 0];
										}
										if (yDof > -1) {
											xyzList[waveOffset] = yDof;
											xyzCoefs[waveOffset++] = sign * Q[numdofperNode * qOff + 1];
										}
										if (zDof > -1) {
											xyzList[waveOffset] = zDof;
											xyzCoefs[waveOffset++] = sign * Q[numdofperNode * qOff + 2];
										}
									}
								} else {
									if (numdofperNode == 1) {
										if (hDof > -1) {
											HelmList[UsedWaves * index + waveOffset - 1] = hDof;
											HelmCoefs[UsedWaves * index + waveOffset - 1] =
													sign * Q[numdofperNode * qOff + 0];
										}
										if (tDof > -1) {
											TempList[UsedWaves * index + waveOffset - 1] = tDof;
											TempCoefs[UsedWaves * index + waveOffset - 1] =
													sign * Q[numdofperNode * qOff + 0];
										}
									} else {
										if (xDof > -1) {
											xyzList[UsedWaves * index + waveOffset - 3] = xDof;
											xyzCoefs[UsedWaves * index + waveOffset - 3] =
													sign * Q[numdofperNode * qOff + 0];
										}
										if (yDof > -1) {
											xyzList[UsedWaves * index + waveOffset - 2] = yDof;
											xyzCoefs[UsedWaves * index + waveOffset - 2] =
													sign * Q[numdofperNode * qOff + 1];
										}
										if (zDof > -1) {
											xyzList[UsedWaves * index + waveOffset - 1] = zDof;
											xyzCoefs[UsedWaves * index + waveOffset - 1] =
													sign * Q[numdofperNode * qOff + 2];
										}
									}
								}
								UsedWaves++;
							}
					if (UsedWaves != nQAddWaves[iSub])
						std::cerr << " Something is wrong for the number of used waves " << std::endl;
				}
			}
		}
		if (numR > 0)
			if (numdofperNode == 1)
				if (isFluidSub)
					hOffset = off;
				else
					tOffset = off;
			else
				xOffset = off;
		else
			waveOffset = off;
	}
	if (oldTot != total)
		fprintf(stderr, " *** ERROR: total is incorrect %d %d\n", oldTot, total);

	if (numdofperNode == 1) {
		if (isFluidSub)
			this->Grc = std::make_unique<GenCuCSparse<Scalar>>(total, cc_dsa->size(), HelmCount,
				HelmList, std::move(HelmCoefs));
		else
			this->Grc = std::make_unique<GenCuCSparse<Scalar>>(total, cc_dsa->size(), TempCount,
				TempList, std::move(TempCoefs));
	} else {
		this->Grc = std::make_shared<GenCuCSparse<Scalar>>(total, cc_dsa->size(), xyzCount,
			xyzList, std::move(xyzCoefs));
	}

	int ii = (isFluidSub || isThermalSub) ? 1 : 0;
	if (edgeQindex[ii] == -1)
		edgeQindex[ii] = this->Src->addSparseMatrix(
				this->Grc);  // store index for possible rebuild (multiple LHS freq sweep)
	else
		this->Src->setSparseMatrix(edgeQindex[ii], this->Grc);

	for (i = 0; i < scomm->numNeighb; ++i) edgeDofSize[i] += edgeDofSizeTmp[i];

	// Clear the unused dofs JAT 112113
	for (iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		if (scomm->subNums[iSub] == subNum()) continue;
		if (numR > 0) {
			if (desired.contains(DofSet::Temp))
				if (edgeDofs[iSub].contains(DofSet::Temp) && !isUsed[0 + iSub * nQPerNeighb])
					edgeDofs[iSub].unmark(DofSet::Temp);
			if (desired.contains(DofSet::Helm))
				if (edgeDofs[iSub].contains(DofSet::Helm) && !isUsed[0 + iSub * nQPerNeighb])
					edgeDofs[iSub].unmark(DofSet::Helm);
			if (desired.contains(DofSet::XYZdisp)) {
				if (edgeDofs[iSub].contains(DofSet::Xdisp) && !isUsed[0 + iSub * nQPerNeighb])
					edgeDofs[iSub].unmark(DofSet::Xdisp);
				if (edgeDofs[iSub].contains(DofSet::Ydisp) && !isUsed[1 + iSub * nQPerNeighb])
					edgeDofs[iSub].unmark(DofSet::Ydisp);
				if (edgeDofs[iSub].contains(DofSet::Zdisp) && !isUsed[2 + iSub * nQPerNeighb])
					edgeDofs[iSub].unmark(DofSet::Zdisp);
			}
			if (desired.contains(DofSet::XYZrot)) {
				if (edgeDofs[iSub].contains(DofSet::Xrot) && !isUsed[3 + iSub * nQPerNeighb])
					edgeDofs[iSub].unmark(DofSet::Xrot);
				if (edgeDofs[iSub].contains(DofSet::Yrot) && !isUsed[4 + iSub * nQPerNeighb])
					edgeDofs[iSub].unmark(DofSet::Yrot);
				if (edgeDofs[iSub].contains(DofSet::Zrot) && !isUsed[5 + iSub * nQPerNeighb])
					edgeDofs[iSub].unmark(DofSet::Zrot);
			}
		}
	}

}

template<class Scalar>
void
FetiSub<Scalar>::extractInterfRBMs(int numRBM, Scalar *locRBMs, Scalar *locInterfRBMs) {
	int iDof, iRBM;
	int locLen = localLen();

	// locInterfRBMs are ordered by rbm
	for (iRBM = 0; iRBM < numRBM; ++iRBM) {
		int off = iRBM * scomm->numT(SComm::std);
		int locOff = iRBM * locLen;
		for (iDof = 0; iDof < scomm->numT(SComm::std); ++iDof)
			locInterfRBMs[off + iDof] = locRBMs[locOff + scomm->boundDofT(SComm::std, iDof)];
	}
}

template<class Scalar>
void
FetiSub<Scalar>::sendInterfRBMs(int numRBM, Scalar *locInterfRBMs, FSCommPattern<Scalar> *rbmPat) {
	// locInterfRBMs are ordered by-RBM, need to convert to by-neighbor ordering
	for (int iSub = 0; iSub < scomm->numT(SComm::std); ++iSub) {
		FSSubRecInfo<Scalar> sInfo = rbmPat->getSendBuffer(subNum(), scomm->neighbT(SComm::std, iSub));
		int off = 0;
		for (int iRBM = 0; iRBM < nGrbm; ++iRBM) {
			int locOff = iRBM * totalInterfSize + scomm->offsetT(SComm::std, iSub);
			for (int iDof = 0; iDof < scomm->lenT(SComm::std, iSub); ++iDof) {
				sInfo.data[off + iDof] = interfaceRBMs[locOff + iDof];
			}
			off += scomm->lenT(SComm::std, iSub);
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::recvInterfRBMs(int iNeighb, int numNeighbRBM, Scalar *neighbInterfRBMs,
                                     FSCommPattern<Scalar> *rbmPat) {
	// this is for just one neighbor
	FSSubRecInfo<Scalar> rInfo = rbmPat->recData(scomm->neighbT(SComm::std, iNeighb), subNum());
	int len = scomm->lenT(SComm::std, iNeighb) * numNeighbRBM;
	for (int j = 0; j < len; ++j) neighbInterfRBMs[j] = rInfo.data[j];
}

template<class Scalar>
void
FetiSub<Scalar>::sendInterfaceGrbm(FSCommPattern<Scalar> *rbmPat) {
	// Sub-Domain based augmented preconditioner for DP
	this->Grc = std::make_unique<GenCuCSparse<Scalar>>(scomm->lenT(SComm::std), scomm->boundDofsT(SComm::std), nGrbm,
	                                                   interfaceRBMs.data());
	this->Src->addSparseMatrix(this->Grc);

	sendInterfRBMs(nGrbm, interfaceRBMs.data(), rbmPat);
}

template<class Scalar>
void
FetiSub<Scalar>::receiveInterfaceGrbm(FSCommPattern<Scalar> *rbmPat) {
	// Sub-Domain based augmented preconditioner for DP
	for (int iNeighb = 0; iNeighb < scomm->numT(SComm::std); ++iNeighb) {
		int leadingDimGs = getSComm()->lenT(SComm::std, iNeighb);
		Scalar *neighbGs = new Scalar[leadingDimGs * neighbNumGRBMs[iNeighb]];
		recvInterfRBMs(iNeighb, neighbNumGRBMs[iNeighb], neighbGs, rbmPat);
		// TODO FIX THIS LEAK! Src does not own the pointers. In other cases they are owned by FetiSub.
		auto Grc = std::make_shared<GenCuCSparse<Scalar>>(leadingDimGs, scomm->boundDofsT(SComm::std, iNeighb),
		                                                  neighbNumGRBMs[iNeighb], neighbGs, leadingDimGs);
		Grc->negate();
		this->Src->addSparseMatrix(Grc);
	}
}

template<class Scalar>
void
FetiSub<Scalar>::weightEdgeGs() {
	// for WeightedEdge augmentation
	if (getFetiInfo().scaling == FetiInfo::tscaling)
		this->Grc->doWeighting(weight.data());
	else if (getFetiInfo().scaling == FetiInfo::kscaling)
		this->Grc->doWeighting(this->kweight.data());
}


template<class Scalar>
void
FetiSub<Scalar>::rebuildKbb() {
	this->Kbb.reset();
	this->KiiSparse.reset();
	this->Kib.reset();
	if (numMPC) makeKbbMpc(); else makeKbb(getCCDSA());
}

template<class Scalar>
void
FetiSub<Scalar>::scaleAndSplitKww() {
	if (neighbKww != 0)
		neighbKww->scale(HData::cscale_factor);

	if (neighbKww != 0)
		neighbKww->split(glToLocalWImap, this->wweight.data());

#ifdef HB_COUPLED_PRECOND
	if(solInfo().isCoupled & isMixedSub & this->neighbKww!=0 & KiiSparse!=0) {
   fprintf(stderr," ... Assemble localFsi into Kii in sub %2d\n",subNum());
   if(solInfo().getFetiInfo().splitLocalFsi) {
     neighbKww->splitLocalFsi(glToLocalWImap, this->wweight);
     neighbKww->addLocalFsiToMatrix(KiiSparse.get(), dsa, glToLocalNode);
   } else {
     fprintf(stderr," ... No local Fsi spliting in sub %2d\n",subNum());
     neighbKww->addLocalFsiToMatrix(KiiSparse.get(), dsa, glToLocalNode, kSumWI);
   }
 }
#endif
	prev_cscale_factor = HData::cscale_factor;
}

template<class Scalar>
void
FetiSub<Scalar>::reScaleAndReSplitKww() {
	double rescale_factor = HData::cscale_factor / prev_cscale_factor;

	if (this->neighbKww != 0)
		this->neighbKww->scale(rescale_factor);

	if (this->neighbKww != 0)
		if (getFetiInfo().fsi_scaling == FetiInfo::kscaling)
			this->neighbKww->split(glToLocalWImap, this->wweight.data());

#ifdef HB_COUPLED_PRECOND
	if(solInfo().isCoupled & isMixedSub & this->neighbKww!=0) {
   fprintf(stderr," ... Assemble localFsi into Kii in sub %2d\n",subNum());
   if(KiiSparse) this->neighbKww->addLocalFsiToMatrix(KiiSparse.get(), dsa, glToLocalNode);
 }
#endif
	prev_cscale_factor = HData::cscale_factor;
}

template<class Scalar>
void
FetiSub<Scalar>::makeZstarAndR(double *centroid)
{
	rigidBodyModesG = std::make_unique<Rbm>(getDsa(), get_c_dsa(), getNodeSet(), tolsvd,
	                                        centroid+3*group, cornerNodes, numCRN, numCRNdof, cornerDofs,
	                                        numMPC_primal, this->mpc_primal);
}

extern void
getJacobi(double *kappa, double * mu, FullSquareMatrix &xx,
          double *eigVal,int nsmax, int subSpaceSize, double tolJac);

template<>
void
FetiSub<DComplex>::precondGrbm()
{
	fprintf(stderr, " *** WARNING: FetiSub<DComplex>::precondGrbm() not implemented \n");
}

template<>
void
FetiSub<double>::precondGrbm()
{
	// WARNING: this routine has to be modified to test any new FETI-DP
	//          coarse grid ideas, if this return is uncommented, then
	//          the augmented coarse grid used the geometric rbms to enforce
	//          G^t Bu=0
	// return;
	// Just checking

	// 1. Augmenting Kcc by the Geometric Gs
	if(nGrbm == 0)
		return;

	auto boundDofs = scomm->boundDofsT(SComm::std);
	int iRBM, iDof;
	int iLen = scomm->lenT(SComm::std);
	int locLen = localRLen();

	//double *v = (double *) dbg_alloca(iLen*sizeof(double));
	double *y = (double *) dbg_alloca(nGrbm*locLen*sizeof(double));

// 2. Preconditioned Gs
/*
 for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
   for(iDof = 0; iDof < iLen; ++iDof)
   v[iDof]=interfaceRBMs[iRBM*iLen+iDof];
   multKbb(v,interfaceRBMs+iRBM*iLen);
 }
 return;
*/

/*
// 3. Try multipling by Kbb's diagonal values?
 for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
   for(iDof = 0; iDof < iLen; ++iDof)
      v[iDof]=interfaceRBMs[iRBM*iLen+iDof];
   multDiagKbb(v,interfaceRBMs+iRBM*iLen);
 }
 return;
*/

	int j;

// 4. just the local solver eigen vectors
/*
 for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
   for(iDof=0; iDof<locLen; ++iDof)
     y[iDof+iRBM*locLen]=0.0;
   for(iDof=0; iDof<iLen; ++iDof)
     y[boundDofs[iDof]+iRBM*locLen] += interfaceRBMs[iRBM*iLen+iDof];
 }
 for(j=0; j<5; ++j) {
   for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
     if(Krr) Krr->reSolve(y+iRBM*locLen);
   }
   // orthonormalize interfaceRBMs
   ortho(y, nGrbm, locLen);
 }

 for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
   for(iDof=0; iDof<iLen; ++iDof)
     interfaceRBMs[iRBM*iLen+iDof] = y[boundDofs[iDof]+iRBM*locLen];
 }
*/

	// Michel's new idea
	// 5. generalized eigen value problem
	// orthogonalize the RBMS
	// ortho(interfaceRBMs,nGrbm,iLen);
	// copy RBMS into some vectors
	//double *origR = (double *) dbg_alloca(nGrbm*locLen*sizeof(double));
	//for(iRBM = 0; iRBM < nGrbm; ++iRBM)
	// for(iDof=0; iDof<iLen; ++iDof)
	//   origR[iDof+iRBM*iLen] = interfaceRBMs[iRBM*iLen+iDof];

	// within the iteration, orthogonalize wrt original RBMS

	// Order 1: Krr^-1, Preconditioner
	double *mu    = (double *) dbg_alloca(sizeof(double)*nGrbm*(nGrbm+1)/2);
	double *kappa = (double *) dbg_alloca(sizeof(double)*nGrbm*(nGrbm+1)/2);
	double *x = (double *) dbg_alloca(sizeof(double)*iLen*nGrbm);

/*
 for(iRBM=nGrbm/2; iRBM<nGrbm; ++iRBM)
 for(iDof=0; iDof<iLen; ++iDof) {
   interfaceRBMs[iRBM*iLen+iDof]=(iDof+iRBM)%nGrbm;
 }
*/

	int numIterations=6;
	int i;
	GenFullSquareMatrix<double> xx(nGrbm);
	double *eval = (double *) dbg_alloca(sizeof(double)*nGrbm);

	for(int jj=0; jj<numIterations; ++jj) {
		int count=0;
		for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
			for(iDof=0; iDof<locLen; ++iDof)
				y[iDof+iRBM*locLen]=0.0;
			for(iDof=0; iDof<iLen; ++iDof) {
				y[boundDofs[iDof]+iRBM*locLen] += interfaceRBMs[iRBM*iLen+iDof];
//     y[boundDofs[iDof]+iRBM*locLen] +=
//              interfaceRBMs[iRBM*iLen+iDof]/scaling[iDof];
				x[iRBM*iLen+iDof] = interfaceRBMs[iRBM*iLen+iDof];
			}
		}

		if(Krr) Krr->reSolve(nGrbm, y);

		for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
			for(iDof=0; iDof<iLen; ++iDof)
				interfaceRBMs[iRBM*iLen+iDof] =
						y[boundDofs[iDof]+iRBM*locLen];
////                       y[boundDofs[iDof]+iRBM*locLen]/scaling[iDof];

			for(j=0; j<=iRBM; ++j)
			{ GenStackVector<double> Z(x+iRBM*iLen,iLen);
				GenStackVector<double> Q(interfaceRBMs.data()+j*iLen,iLen);
				kappa[count++] = Z*Q; }
		}
		count = 0;
		for(iRBM = 0; iRBM < nGrbm; ++iRBM) {

			for(iDof = 0; iDof < iLen; ++iDof)
				x[iRBM*iLen+iDof] = interfaceRBMs[iRBM*iLen+iDof];
			multKbb(x+iRBM*iLen,interfaceRBMs.data()+iRBM*iLen);

			for(j=0; j<=iRBM; ++j)
			{ GenStackVector<double> Z(x+iRBM*iLen,iLen);
				GenStackVector<double> Q(interfaceRBMs.data()+j*iLen,iLen);
				mu[count++] = Z*Q;
			}

		}
		for(iRBM = 0; iRBM < nGrbm; ++iRBM)
			for(iDof = 0; iDof < iLen; ++iDof)
				x[iRBM*iLen+iDof] = interfaceRBMs[iRBM*iLen+iDof];
		//orthoWithOrigR(origR,interfaceRBMs, nGrbm, iLen);
		//ortho(interfaceRBMs,nGrbm,iLen);
		getJacobi(kappa, mu, xx, eval, 20, nGrbm, 1e-4);

		for(j = 0; j < nGrbm; ++j) {
			GenStackVector<double> Q(interfaceRBMs.data()+j*iLen,iLen);
			Q.zero();
			for(i = 0; i < nGrbm; ++i) {
				GenStackVector<double> Z(x+i*iLen,iLen);
				Q.linAdd(xx[j][i], Z);
			}
		}
	}

/*
 // remultiply by scaling and then Kbb
   for(iRBM = 0; iRBM < nGrbm; ++iRBM) {

     for(iDof = 0; iDof < iLen; ++iDof)
       x[iRBM*iLen+iDof] = interfaceRBMs[iRBM*iLen+iDof];
       //x[iRBM*iLen+iDof] = interfaceRBMs[iRBM*iLen+iDof]*scaling[iDof];
     multKbb(x+iRBM*iLen,interfaceRBMs+iRBM*iLen);
   }
*/
// add Krr^-1 here
	for(iRBM = 0; iRBM < nGrbm; ++iRBM) {
		for(iDof=0; iDof<locLen; ++iDof)
			y[iDof+iRBM*locLen]=0.0;
		for(iDof=0; iDof<iLen; ++iDof) {
			y[boundDofs[iDof]+iRBM*locLen] += interfaceRBMs[iRBM*iLen+iDof];
		}

		Krr->reSolve(y+iRBM*locLen);

		for(iDof=0; iDof<iLen; ++iDof)
			interfaceRBMs[iRBM*iLen+iDof] =
					y[boundDofs[iDof]+iRBM*locLen];
	}
}

template<class Scalar>
void
FetiSub<Scalar>::assembleGlobalCCtsolver(GenSolver<Scalar> *CCtsolver, SimpleNumberer *mpcEqNums) {
	auto &mpc = this->mpc;

	int i, j, k, l;
	for (i = 0; i < numMPC; ++i) {
		Scalar dotii = 0.0;
		int gi = localToGlobalMPC[i];
		int renum_gi = mpcEqNums->firstdof(gi);
		if (mpc[i]->active) {
			CCtsolver->addone(1.0, renum_gi, renum_gi);
			continue;
		} // trick to prevent singularities in CCt when rebuit for contact
		for (k = 0; k < mpc[i]->nterms; ++k) {
			int dof = (mpc[i]->terms)[k].cdof;
			if (dof >= 0)
				dotii += mpc[i]->terms[k].coef * mpc[i]->terms[k].coef / mpc[i]->k[k]; // for mpc kscaling
		}
		CCtsolver->addone(dotii, renum_gi, renum_gi);
		for (j = 0; j < localMpcToMpc->num(i); ++j) {
			Scalar dotij = 0.0;
			int lj = (*localMpcToMpc)[i][j];
			if (mpc[lj]->active) continue;
			int gj = localToGlobalMPC[lj];
			int renum_gj = mpcEqNums->firstdof(gj);
			if (renum_gj > renum_gi) {  // work with upper symmetric half
				// now find matching dof/s
				for (k = 0; k < mpc[i]->nterms; ++k) {
					int dofk = (mpc[i]->terms)[k].cdof;
					if (dofk >= 0) {
						for (l = 0; l < mpc[lj]->nterms; ++l) {
							int dofl = (mpc[lj]->terms)[l].cdof;
							if (dofk == dofl) {
								dotij +=
										mpc[i]->terms[k].coef * mpc[lj]->terms[l].coef / mpc[lj]->k[l]; // for mpc kscaling
							}
						}
					}
				}
				CCtsolver->addone(dotij, renum_gi, renum_gj);
			}
		}
	}
}

// HB: this method add/assemble the locally stored contributions into the global CCt matrix.
//     MUST be called SEQUENTIALLY to avoid writting concurrently at the same memory location.
template<class Scalar>
void
FetiSub<Scalar>::assembleGlobalCCtsolver(GenSolver<Scalar> *CCtsolver) {
	for (int i = 0; i < lengthCCtData; i++)
		CCtsolver->addone(CCtval[i], CCtrow[i], CCtcol[i]);

	CCtrow.clear();
	CCtcol.clear();
	CCtval.clear();
	lengthCCtData = 0;
}


template<class Scalar>
void
FetiSub<Scalar>::constructLocalCCtsolver() {
	// Step 1. initialize solver object
	SimpleNumberer *mpcEqNums = new SimpleNumberer(numMPC);
	for (int i = 0; i < numMPC; ++i) mpcEqNums->setWeight(i, 1);
	mpcEqNums->makeOffset();
	localCCtsolver.reset(GenSolverFactory<Scalar>::getFactory()->createSolver(localMpcToGlobalMpc, mpcEqNums,
	                                                                          *getFetiInfo().cct_cntl,
	                                                                          localCCtsparse));
}


//HB: attempt to improve parallel efficiency of assembling the global CCt matrix
//    This method ONLY compute the subdomain contributions to the global CCt matrix,
//    and store them (values and their respective row & col position in the global CCt matrix)
//    LOCALLY. Thus this method can be executed CONCURRENTLY.
//    The assembling in the global CCt matrix is done SEQUENTIALLY (see method below this one).
template<class Scalar>
void
FetiSub<Scalar>::computeSubContributionToGlobalCCt(SimpleNumberer *mpcEqNums) {

	int i, j, k, l;
	// Step 1. Determine the size of the array & allocate array
	lengthCCtData = 0;
	// this is an upper estimate of the required array size
	// -> nearly 2x the required size

	// this is the exact required array size
	for (i = 0; i < numMPC; ++i) {
		if (mpc[i]->active) continue;
		int gi = localToGlobalMPC[i];
		int renum_gi = mpcEqNums->firstdof(gi);
		lengthCCtData++; //diagonal term
		for (j = 0; j < localMpcToMpc->num(i); ++j) {
			int lj = (*localMpcToMpc)[i][j];
			if (mpc[lj]->active) continue;
			int gj = localToGlobalMPC[lj];
			int renum_gj = mpcEqNums->firstdof(gj);
			if (renum_gj > renum_gi)   // work with upper symmetric half
				lengthCCtData++;
		}
	}
	CCtrow.resize(lengthCCtData);
	CCtcol.resize(lengthCCtData);
	CCtval.resize(lengthCCtData);

	// Step 2. Fill the array
	lengthCCtData = 0; // use it as counter (at the end, it should be the exact number of contributions)
	for (i = 0; i < numMPC; ++i) {
		if (mpc[i]->active) continue;
		Scalar dotii = 0.0;
		int gi = localToGlobalMPC[i];
		int renum_gi = mpcEqNums->firstdof(gi);
		for (k = 0; k < mpc[i]->nterms; ++k) {
			int dof = (mpc[i]->terms)[k].cdof;
			if (dof >= 0)
				dotii += mpc[i]->terms[k].coef * mpc[i]->terms[k].coef / mpc[i]->k[k]; // for mpc kscaling
		}
		CCtrow[lengthCCtData] = renum_gi;
		CCtcol[lengthCCtData] = renum_gi;
		CCtval[lengthCCtData] = dotii;
		lengthCCtData++;
		for (j = 0; j < localMpcToMpc->num(i); ++j) {
			Scalar dotij = 0.0;
			int lj = (*localMpcToMpc)[i][j];
			if (mpc[lj]->active) continue;
			int gj = localToGlobalMPC[lj];
			int renum_gj = mpcEqNums->firstdof(gj);
			if (renum_gj > renum_gi) {  // work with upper symmetric half
				// now find matching dof/s
				for (k = 0; k < mpc[i]->nterms; ++k) {
					int dofk = (mpc[i]->terms)[k].cdof;
					if (dofk >= 0) {
						for (l = 0; l < mpc[lj]->nterms; ++l) {
							int dofl = (mpc[lj]->terms)[l].cdof;
							if (dofk == dofl) {
								dotij +=
										mpc[i]->terms[k].coef * mpc[lj]->terms[l].coef / mpc[lj]->k[l]; // for mpc kscaling
							}
						}
					}
				}
				CCtrow[lengthCCtData] = renum_gi;
				CCtcol[lengthCCtData] = renum_gj;
				CCtval[lengthCCtData] = dotij;
				lengthCCtData++;
			}
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::solveLocalCCt(Scalar *subv) {
	// used for block diagonal CCt preconditioner for MPCs
	if (numMPC > 0) {
		int i;

		// Step 1: extract the mpc residual components from subv and inserts them in a local vector (mpcv)
		GenVector<Scalar> mpcv(numMPC, 0.0);
		for (i = 0; i < scomm->lenT(SComm::mpc); ++i)
			mpcv[scomm->mpcNb(i)] = subv[scomm->mapT(SComm::mpc, i)];

		// Step 2: solve CCt^-1 * mpcv
		localCCtsolver->reSolve(mpcv);

		// Step 3: redistribute mpcv to the interface vector subv
		for (i = 0; i < scomm->lenT(SComm::mpc); ++i)
			subv[scomm->mapT(SComm::mpc, i)] = mpcv[scomm->mpcNb(i)];
	}
}

template<class Scalar>
void
FetiSub<Scalar>::assembleLocalCCtsolver() {
	auto &mpc = this->mpc;
	// Step 2. add local mpc CC^t terms to solver
	int i, j, k, l;
	for (i = 0; i < numMPC; ++i) {
		if (mpc[i]->active) continue;
		Scalar dotii = 0.0;
		for (k = 0; k < mpc[i]->nterms; ++k) {
			int dof = (mpc[i]->terms)[k].cdof;
			if (dof >= 0)
				dotii += mpc[i]->terms[k].coef * mpc[i]->terms[k].coef / mpc[i]->k[k]; // for mpc kscaling;
		}
		localCCtsolver->addone(dotii, i, i);
		for (j = 0; j < localMpcToMpc->num(i); ++j) {
			Scalar dotij = 0.0;
			int lj = (*localMpcToMpc)[i][j];
			if (mpc[lj]->active) continue;
			if (lj > i) { // work with upper symmetric half only
				// now find matching dof/s
				for (k = 0; k < mpc[i]->nterms; ++k) {
					int dofk = (mpc[i]->terms)[k].cdof;
					if (dofk >= 0) {
						for (l = 0; l < mpc[lj]->nterms; ++l) {
							int dofl = (mpc[lj]->terms)[l].cdof;
							if (dofk == dofl) {
								dotij +=
										mpc[i]->terms[k].coef * mpc[lj]->terms[l].coef / mpc[lj]->k[l]; // for mpc kscaling
								//HB: can probably be optimized by assuming that no dupplicate term (i.e. dof) exist in a lmpc
								//    so that we can stop the l loop (break) when dofk == dofl ?? The non dupplicate assumption
								//    is or can be enforced in a preprocess step
							}
						}
					}
				}
				localCCtsolver->addone(dotij, i, lj);
			}
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::setCCtCommSize(FSCommPattern<Scalar> *cctPat) {
	// note: this is an upper bound, the required length will be less
	for (int i = 0; i < scomm->numT(SComm::mpc); ++i)
		cctPat->setLen(subNum(), scomm->neighbT(SComm::mpc, i), localMpcToGlobalMpc->numConnect());
}

template<class Scalar>
void
FetiSub<Scalar>::sendNeighbCCtsolver(FSCommPattern<Scalar> *cctPat, const Connectivity *mpcToSub) {
	int i, j, k;
	for (i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc, i);
		if (subNum() != neighb) {
			int count = 0;
			FSSubRecInfo<Scalar> sInfo = cctPat->getSendBuffer(subNum(), neighb);
			for (j = 0; j < numMPC; ++j) {
				int gj = localToGlobalMPC[j];
				if (mpcToSub->offset(gj, neighb) != ((size_t) -1) ) {
					sInfo.data[count++] = localCCtsolver->getone(j, j);
					for (k = j + 1; k < numMPC; ++k) {
						int gk = localToGlobalMPC[k];
						if ((mpcToSub->offset(gk, neighb) > -1) && (localMpcToGlobalMpc->offset(j, k) > -1)) {
							sInfo.data[count] = localCCtsolver->getone(j, k);
							count++;
						}
					}
				}
			}
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::recNeighbCCtsolver(FSCommPattern<Scalar> *cctPat, const Connectivity *mpcToSub) {
	int i, j, k;
	for (i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc, i);
		if (subNum() != neighb) {
			int count = 0;
			FSSubRecInfo<Scalar> rInfo = cctPat->recData(neighb, subNum());
			for (j = 0; j < numMPC; ++j) {
				int gj = localToGlobalMPC[j];
				if (mpcToSub->offset(gj, neighb) > -1) {
					localCCtsolver->addone(rInfo.data[count++], j, j);
					for (k = j + 1; k < numMPC; ++k) {
						int gk = localToGlobalMPC[k];
						if ((mpcToSub->offset(gk, neighb) > -1) && (localMpcToGlobalMpc->offset(j, k) > -1)) {
							localCCtsolver->addone(rInfo.data[count], j, k);
							count++;
						}
					}
				}
			}
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::factorLocalCCtsolver() {
	localCCtsolver->factor();
	int numCCtSing = localCCtsolver->numRBM();
	if (numCCtSing > 0)
		std::cerr << "sub = " << subNum() << ", Number of singularities in CCt = "
		          << numCCtSing << std::endl;
}

template<class Scalar>
void
FetiSub<Scalar>::zeroLocalCCtsolver() {
	localCCtsolver->zeroAll();
}

template<class Scalar>
void
FetiSub<Scalar>::assembleBlockCCtsolver(int iBlock, GenSolver<Scalar> *CCtsolver, SimpleNumberer *blockMpcEqNums) {
	auto &mpc = this->mpc;
	// optimize by looping ONLY over the lmpc eqs contributed to this iBlock
	// Make sure that the input block (global) Id iBlock is in the range of blocks this subdomain
	// contributes to
	int i, j, k, l, p;
	for (p = 0; p < blockToLocalMpc->num(iBlock); ++p) {
		i = (*blockToLocalMpc)[iBlock][p];
		if (mpc[i]->active) continue;
		int bi = (*blockToBlockMpc)[iBlock][p];
		int renum_bi = blockMpcEqNums->firstdof(bi);
		Scalar dotii = 0.0;
		for (k = 0; k < mpc[i]->nterms; ++k) {
			int dof = (mpc[i]->terms)[k].cdof;
			if (dof >= 0)
				dotii += mpc[i]->terms[k].coef * mpc[i]->terms[k].coef / mpc[i]->k[k]; // for mpc kscaling;
		}
		CCtsolver->addone(dotii, renum_bi, renum_bi);
		for (j = 0; j < localMpcToMpc->num(i); ++j) {
			int lj = (*localMpcToMpc)[i][j];
			if (mpc[lj]->active) continue;
			int jb = localMpcToBlock->cOffset(lj, iBlock);
			if (jb > -1) {
				int bj = (*localMpcToBlockMpc)[lj][jb];
				int renum_bj = blockMpcEqNums->firstdof(bj);
				if (renum_bj > renum_bi) { // work with upper symmetric part only
					Scalar dotij = 0.0;
					// now find matching dof/s
					for (k = 0; k < mpc[i]->nterms; ++k) {
						int dofk = (mpc[i]->terms)[k].cdof;
						if (dofk >= 0) {
							for (l = 0; l < mpc[lj]->nterms; ++l) {
								int dofl = (mpc[lj]->terms)[l].cdof;
								if (dofk == dofl)
									dotij += mpc[i]->terms[k].coef * mpc[lj]->terms[l].coef /
									         mpc[lj]->k[l]; // for mpc kscaling;
							}
						}
					}
					CCtsolver->addone(dotij, renum_bi, renum_bj);
				}
			}
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::extractMpcResidual(Scalar *subv, GenVector<Scalar> &mpcv, SimpleNumberer *mpcEqNums) {
	// extracts the mpc residual components from subv and inserts them in global vector (mpcv)
	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		mpcv[mpcEqNums->firstdof(localToGlobalMPC[locMpcNb])] = subv[scomm->mapT(SComm::mpc, i)];
	}
}

template<class Scalar>
void
FetiSub<Scalar>::insertMpcResidual(Scalar *subv, GenVector<Scalar> &mpcv, SimpleNumberer *mpcEqNums) {
	// extracts the mpc residual components from mpcv and inserts them in the interface vector subv
	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		subv[scomm->mapT(SComm::mpc, i)] = mpcv[mpcEqNums->firstdof(localToGlobalMPC[locMpcNb])];
	}
}

template<class Scalar>
void
FetiSub<Scalar>::extractBlockMpcResidual(int block, Scalar *subv, GenVector<Scalar> *mpcv,
                                              SimpleNumberer *blockMpcEqNums) {
	// extracts the mpc residual components from subv and inserts them in global vector (mpcv)
	// modified to loop only over the subdomain lmpcs belonging to block
	// Make sure that the input block (global) Id is in the range of blocks Id that have contribution
	// from this subdomain
	int i, j;
	for (j = 0; j < blockToLocalMpc->num(block); ++j) {
		i = (*blockToLocalMpc)[block][j];
		int iDof = (*mpcToBoundDof)[i][0];
		int bij = (*blockToBlockMpc)[block][j];
		(*mpcv)[blockMpcEqNums->firstdof(bij)] = subv[iDof];
	}
}

template<class Scalar>
void
FetiSub<Scalar>::insertBlockMpcResidual(Scalar *subv, GenVector<Scalar> **mpcv,
                                        const Connectivity *mpcToBlock,
                                        SimpleNumberer **blockMpcEqNums)
{
	// extracts the mpc residual components from mpcv and inserts them in the interface vector subv
	int i, j, k;
	for (i = 0; i < numMPC; ++i) {
		int iDof = (*mpcToBoundDof)[i][0];
		int gi = localToGlobalMPC[i];
		double w = double(mpcToBlock->num(gi));
		subv[iDof] = 0.0;
		for (j = 0; j < localMpcToBlock->num(i); ++j) {
			int jBlock = (*localMpcToBlock)[i][j];
			int bij = (*localMpcToBlockMpc)[i][j];
			subv[iDof] += (*mpcv[jBlock])[blockMpcEqNums[jBlock]->firstdof(bij)] / w;
		}
		for (k = 1; k < mpcToBoundDof->num(i); ++k) {
			int kDof = (*mpcToBoundDof)[i][k];
			subv[kDof] = subv[(*mpcToBoundDof)[i][0]];
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::sendMpcInterfaceVec(FSCommPattern<Scalar> *mpcPat, Scalar *interfvec) {
	for (int i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc, i);
		if (subNum() != neighb) {
			FSSubRecInfo<Scalar> sInfo = mpcPat->getSendBuffer(subNum(), neighb);
			for (int j = 0; j < scomm->lenT(SComm::mpc, i); ++j)
				sInfo.data[j] = interfvec[scomm->mapT(SComm::mpc, i, j)];
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::combineMpcInterfaceVec(FSCommPattern<Scalar> *mpcPat, Scalar *interfvec) {
	// this function combines the interfvec values corresponding to mpcs
	// that are shared between subdomains. Used with approximated blockDiag and diag CCt preconditioners
	int i, j;
	Scalar *mpcCombo = (Scalar *) dbg_alloca(sizeof(Scalar) * numMPC);
	int *mpcCount = (int *) dbg_alloca(sizeof(int) * numMPC);
	for (i = 0; i < numMPC; ++i) mpcCount[i] = 1;
	for (i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc, i);
		if (subNum() != neighb) {
			FSSubRecInfo<Scalar> rInfo = mpcPat->recData(neighb, subNum());
			for (j = 0; j < scomm->lenT(SComm::mpc, i); ++j) {
				int locMpcNb = scomm->mpcNb(i, j);
				if (mpcCount[locMpcNb] == 1)
					mpcCombo[locMpcNb] = interfvec[scomm->mapT(SComm::mpc, i, j)];
				mpcCount[locMpcNb]++;
				mpcCombo[locMpcNb] += rInfo.data[j];
			}
		}
	}
	// now add mpcCombo to interfvec
	for (i = 0; i < scomm->numT(SComm::mpc); ++i) {
		int neighb = scomm->neighbT(SComm::mpc, i);
		if (subNum() != neighb) {
			for (j = 0; j < scomm->lenT(SComm::mpc, i); ++j) {
				int locMpcNb = scomm->mpcNb(i, j);
				interfvec[scomm->mapT(SComm::mpc, i, j)] = mpcCombo[locMpcNb] / double(mpcCount[locMpcNb]);
			}
		}
	}
}

template<class Scalar>
void
FetiSub<Scalar>::getHalfInterf(const Scalar *s, Scalar *t) const {
	int iTg = 0;
	int i;
	for (i = 0; i < totalInterfSize; ++i)
		if (masterFlag[i]) t[iTg++] = s[i];
}

template<class Scalar>
void
FetiSub<Scalar>::getHalfInterf(const Scalar *s, Scalar *t, const Scalar *ss, Scalar *tt) const {
	int iTg = 0;
	int i;
	for (i = 0; i < totalInterfSize; ++i)
		if (masterFlag[i]) {
			t[iTg] = s[i];
			tt[iTg++] = ss[i];
		}
}

template<class Scalar>
void
FetiSub<Scalar>::scatterHalfInterf(const Scalar *s, Scalar *buffer) const {
	// note: need to send s[iTg] of mpc master to all neighbors
	int iTg = 0;
	int i;

	Scalar *mpcbuff = (Scalar *) dbg_alloca(numMPC * sizeof(Scalar));
	Scalar *wibuff = (Scalar *) dbg_alloca(numWIdof * sizeof(Scalar));

	// note: if numMPC changes (salinas) need to delete [] mpcMaster
	for (i = 0; i < totalInterfSize; ++i) {
		if (masterFlag[i]) {
			switch (boundDofFlag[i]) {
				case 0:
					iTg++;
					break;
				case 1: { // wet interface
					int windex = -1 - allBoundDofs[i];
					wibuff[windex] = s[iTg++];
				}
					break;
				case 2: { // dual mpc
					int locMpcNb = -1 - allBoundDofs[i];
					mpcbuff[locMpcNb] = s[iTg++];
				}
					break;
			}
		}
	}

	iTg = 0;
	for (i = 0; i < totalInterfSize; ++i) {
		if (masterFlag[i]) buffer[i] = s[iTg++];
		else {
			switch (boundDofFlag[i]) {
				case 1: { // wet interface
					int windex = -1 - allBoundDofs[i];
					buffer[i] = (wiMaster[windex]) ? wibuff[windex] : 0.0;
				}
					break;
				case 2: { // dual mpc
					int locMpcNb = -1 - allBoundDofs[i];
					buffer[i] = (mpcMaster[locMpcNb]) ? mpcbuff[locMpcNb] : 0.0;
				}
					break;
			}
		}
	}
}
template<class Scalar>
void
FetiSub<Scalar>::rebuildInterf(Scalar *v, FSCommPattern<Scalar> *vPat) const {
	int iSub, i;
	int iOff = 0;
	Scalar *mpcv = (Scalar *) dbg_alloca(sizeof(Scalar) * numMPC);
	for (i = 0; i < numMPC; ++i) mpcv[i] = 0.0;
	Scalar *wiv = (Scalar *) dbg_alloca(numWIdof * sizeof(Scalar));
	for (i = 0; i < numWIdof; ++i) wiv[i] = 0.0;

	for (iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
		FSSubRecInfo<Scalar> rInfo = vPat->recData(scomm->neighbT(SComm::all, iSub), subNum());
		for (i = 0; i < scomm->lenT(SComm::all, iSub); ++i) {
			int bdof = scomm->boundDofT(SComm::all, iSub, i);
			switch (boundDofFlag[iOff + i]) {
				case 0:
					if (!masterFlag[iOff + i]) v[iOff + i] = -rInfo.data[i];
					break;
				case 1: { // wet interface
					int windex = -1 - bdof;
					if (!masterFlag[iOff + i]) {
						if (rInfo.data[i] != 0.0) wiv[windex] = rInfo.data[i];
					} else wiv[windex] = v[i + iOff];
				}
					break;
				case 2: { // dual mpc
					int locMpcNb = -1 - bdof;
					if (!masterFlag[iOff + i]) {
						if (rInfo.data[i] != 0.0) mpcv[locMpcNb] = rInfo.data[i];
					} else mpcv[locMpcNb] = v[i + iOff];
				}
					break;
			}
		}
		iOff += scomm->lenT(SComm::all, iSub);
	}

	for (i = 0; i < scomm->lenT(SComm::mpc); ++i)
		v[scomm->mapT(SComm::mpc, i)] = mpcv[scomm->mpcNb(i)];

	for (i = 0; i < scomm->lenT(SComm::wet); ++i)
		v[scomm->mapT(SComm::wet, i)] = wiv[scomm->wetDofNb(i)];
}

template<class Scalar>
void
FetiSub<Scalar>::interfaceJump(Scalar *interfvec, FSCommPattern<Scalar> *vPat) const {
	int i, iSub, iDof;
	int offset = 0;

	Scalar *mpcJump = (numMPC) ? new Scalar[numMPC] : 0;
	bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool) * numMPC);
	for (i = 0; i < numMPC; ++i) mpcFlag[i] = true;

	Scalar *wiJump = (Scalar *) dbg_alloca(sizeof(Scalar) * numWIdof);
	bool *wiFlag = (bool *) dbg_alloca(sizeof(bool) * numWIdof);
	for (i = 0; i < numWIdof; ++i) wiFlag[i] = true;

	for (iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
		FSSubRecInfo<Scalar> rInfo = vPat->recData(scomm->neighbT(SComm::all, iSub), subNum());
		for (iDof = 0; iDof < scomm->lenT(SComm::all, iSub); ++iDof) {
			int bdof = scomm->boundDofT(SComm::all, iSub, iDof);
			switch (boundDofFlag[offset + iDof]) {
				case 0:
					interfvec[offset + iDof] -= rInfo.data[iDof];
					break;
				case 1: { // wet interface
					int windex = -1 - bdof;
					if (wiFlag[windex]) {
						wiJump[windex] = interfvec[offset + iDof];
						wiFlag[windex] = false;
					}
					if (subNum() != scomm->neighbT(SComm::all, iSub)) {
						wiJump[windex] += rInfo.data[iDof];
					}
				}
					break;
				case 2: {  // dual mpc
					int locMpcNb = -1 - bdof;
					if (mpcFlag[locMpcNb]) {  // do this 1st time only
						mpcJump[locMpcNb] = interfvec[offset + iDof];
						mpcFlag[locMpcNb] = false;
					}
					if (subNum() != scomm->neighbT(SComm::all, iSub))
						mpcJump[locMpcNb] += rInfo.data[iDof];
				}
					break;
			}
		}
		offset += scomm->lenT(SComm::all, iSub);
	}

	// add mpcJump to interfvec
	for (i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		interfvec[scomm->mapT(SComm::mpc, i)] = mpcJump[scomm->mpcNb(i)];
	}
	// add wiJump to interfvec
	for (i = 0; i < scomm->lenT(SComm::wet); ++i) {
		interfvec[scomm->mapT(SComm::wet, i)] = wiJump[scomm->wetDofNb(i)];
	}
	if (mpcJump) delete[] mpcJump;
}

template<class Scalar>
void
FetiSub<Scalar>::fSend(const Scalar *locF, FSCommPattern<Scalar> *vPat, Scalar *locFw) const {
	// Get the trace of the subdomain interfaces
	int iDof = 0;
	for (int i = 0; i < scomm->numT(SComm::all); ++i) {
		FSSubRecInfo<Scalar> sInfo = vPat->getSendBuffer(subNum(), scomm->neighbT(SComm::all, i));
		for (int j = 0; j < scomm->lenT(SComm::all, i); ++j) {
			int bdof = scomm->boundDofT(SComm::all, i, j);
			switch (boundDofFlag[iDof++]) {
				case 0:
					sInfo.data[j] = locF[bdof];
					break;
				case 1: // wet interface
					sInfo.data[j] = locFw[-1 - bdof];
					break;
			}
		}
	}
}

template<typename Scalar>
void FetiSub<Scalar>::makeKbb(const DofSetArray *dofsetarray) {
	auto nodeToNode = getNodeToNode();
	const auto &dsa = getDsa();
	glBoundMap = makeBMaps(dofsetarray);
	glInternalMap = makeIMaps(dofsetarray);

	this->Kbb = std::make_unique<GenDBSparseMatrix<Scalar>>(nodeToNode, dsa, glBoundMap);
	memPrec += this->Kbb->size();
	// KHP: 3-26-98 Modified this code to always construct Kib

	if (internalLen > 0) {
#ifdef HB_COUPLED_PRECOND
		if(solInfo().isCoupled & isMixedSub & this->neighbKww!=0)
       Kib = new GenCuCSparse<Scalar>(precNodeToNode, dsa, glBoundMap, glInternalMap);
    else
#endif
		Kib = std::make_unique<GenCuCSparse<Scalar>>(nodeToNode, dsa, glBoundMap, glInternalMap.data());
		memPrec += Kib->size();
		Kib->zeroAll();
	} else Kib = 0;

	if ((getFetiInfo().precno == FetiInfo::dirichlet) && (internalLen > 0)) {
#ifdef HB_COUPLED_PRECOND
		if(solInfo().isCoupled & isMixedSub & neighbKww!=0)
      KiiSolver = GenSolverFactory<Scalar>::getFactory()->createSolver(precNodeToNode, dsa, glInternalMap, *sinfo.solvercntl->fetiInfo.kii_cntl, KiiSparse.get());
    else
#endif
		auto solverAndMat = GenSolverFactory<Scalar>::getFactory()->createSolver(nodeToNode, dsa,
		                                                                         glInternalMap.data(),
		                                                                         *getFetiInfo().kii_cntl);

		KiiSolver = std::move(solverAndMat.solver);
		KiiSparse = std::move(solverAndMat.sparseMatrix);

		KiiSolver->setPrintNullity(false);
		// Find approximate preconditioner size
		memPrec += KiiSolver->size();
	} else {
		KiiSparse.reset(nullptr);
		KiiSolver = 0;
	}
}

template<typename Scalar>
void FetiSub<Scalar>::makeKbbMpc() {
	const auto &c_dsa = get_c_dsa();
	// make the mpc dofs be boundary dofs
	std::vector<int> weight_mpc(c_dsa->size(), 0);
	for (int i = 0; i < cc_dsa->size(); ++i) weight_mpc[ccToC[i]] = dofWeight(i);

	for (int i = 0; i < numMPC; i++) {
		const auto &m = mpc[i];
		//XXXX if(mpc[i]->active) continue;
		for (int k = 0; k < m->nterms; k++) {
			int cdof = (m->terms)[k].cdof;
			if (cdof >= 0) weight_mpc[cdof] = 2; // > 1 hence will be in boundary (see makeBmaps and makeImaps)
		}
	}

	// construct the mapping and the matrices
	std::swap(getWeights(), weight_mpc);
	makeKbb(c_dsa);
	// recover
	std::swap(getWeights(), weight_mpc);
}

// TODO Move this to FetiBaseSub by separating the wweight part.
template<class Scalar>
void
FetiSub<Scalar>::gatherDOFList(FSCommPattern<int> *pat) {
	int iSub, iNode, i;
	Connectivity &sharedNodes = *(FetiBaseSub::scomm->sharedNodes); // contains both dry boundary nodes and wet interface nodes
	FetiBaseSub::boundaryDOFs.resize(FetiBaseSub::scomm->numNeighb);
	int nbdofs = 0, nwdofs = 0;
	int nbneighb = 0, nwneighb = 0;
	bool isbneighb, iswneighb;
	bool isCoupled = FetiBaseSub::isCoupled;
	for (iSub = 0; iSub < FetiBaseSub::scomm->numNeighb; ++iSub) {
		FSSubRecInfo<int> rInfo = pat->recData(FetiBaseSub::scomm->subNums[iSub], FetiBaseSub::subNumber);
		FetiBaseSub::boundaryDOFs[iSub].resize(sharedNodes.num(iSub));
		isbneighb = false;
		iswneighb = false;
		for (iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
			FetiBaseSub::boundaryDOFs[iSub][iNode] = DofSet(rInfo.data[iNode]); // temporarily store neighb c_dsa in boundaryDOFs
			// or w dofs if neighb is self
			if (isCoupled && (FetiBaseSub::isWetInterfaceNode(sharedNodes[iSub][iNode]))) { // wet interface node
				DofSet shared_wdofs =
					FetiBaseSub::boundaryDOFs[iSub][iNode] & FetiBaseSub::wetInterfaceDofs[FetiBaseSub::wetInterfaceNodeMap[sharedNodes[iSub][iNode]]];
				int wcount = shared_wdofs.count();
				if (wcount > 0) {
					nwdofs += wcount;
					iswneighb = true;
				}
			}
			DofSet shared_rdofs = FetiBaseSub::boundaryDOFs[iSub][iNode] & (*FetiBaseSub::getCCDSA())[sharedNodes[iSub][iNode]];
			int bcount = shared_rdofs.count();
			nbdofs += bcount;
			isbneighb = true;  // make every "sharedNode" neighbor a "std sharedDOF neighb" also (simplifies augmentation)
		}
		if (iswneighb) nwneighb++;
		if (isbneighb) nbneighb++;
	}

	std::vector<int> boundDofs(nbdofs);
	std::vector<size_t> boundDofPointer(nbneighb + 1);
	boundDofPointer[0] = 0;
	std::vector<int> boundNeighbs(nbneighb);
	std::vector<int> wetDofs(nwdofs);
	std::vector<size_t> wetDofPointer(nwneighb + 1);
	wetDofPointer[0] = 0;
	std::vector<int> wetNeighbs(nwneighb);
	nbdofs = 0;
	nwdofs = 0;
	nbneighb = 0;
	nwneighb = 0;
	const auto &dsa = getDsa();
	for (iSub = 0; iSub < FetiBaseSub::scomm->numNeighb; ++iSub) {
		isbneighb = false;
		iswneighb = false;
		for (iNode = 0; iNode < sharedNodes.num(iSub); ++iNode) {
			// XXXXX THere is uninitialized data involved here!!! TODO FIX
			if (isCoupled && (FetiBaseSub::isWetInterfaceNode(sharedNodes[iSub][iNode]))) { // wet interface node
				DofSet shared_wdofs =
					FetiBaseSub::boundaryDOFs[iSub][iNode] & FetiBaseSub::wetInterfaceDofs[FetiBaseSub::wetInterfaceNodeMap[sharedNodes[iSub][iNode]]];
				int dofs[7]; //, cdofs[7];
				dsa->number(sharedNodes[iSub][iNode], shared_wdofs, dofs);
				int wcount = shared_wdofs.count();
				if (wcount > 0) {
					for (int i = 0; i < wcount; ++i)
						wetDofs[nwdofs++] = FetiBaseSub::wetInterfaceMap[dofs[i]]; // WHAT ABOUT CONSTRAINTS ???
					iswneighb = true;
				}
			}
			FetiBaseSub::boundaryDOFs[iSub][iNode] &= (*FetiBaseSub::getCCDSA())[sharedNodes[iSub][iNode]]; // ~> shared_rdofs
			int bcount = FetiBaseSub::getCCDSA()->number(sharedNodes[iSub][iNode],
			                                             FetiBaseSub::boundaryDOFs[iSub][iNode], boundDofs.data() + nbdofs);
			nbdofs += bcount;
			isbneighb = true; // make every "sharedNode" neighbor a "std sharedDOF neighb" also (simplifies augmentation)
		}
		if (iswneighb) {
			wetNeighbs[nwneighb++] = FetiBaseSub::scomm->subNums[iSub];
			wetDofPointer[nwneighb] = nwdofs;
		}
		if (isbneighb) {
			boundNeighbs[nbneighb++] = FetiBaseSub::scomm->subNums[iSub];
			boundDofPointer[nbneighb] = nbdofs;
		}
	}

	auto stdSharedDOFs = std::make_unique<Connectivity>(nbneighb, std::move(boundDofPointer), std::move(boundDofs));
	int ndof = FetiBaseSub::localLen();
	auto &weight = getWeights();
	weight.assign(ndof, 1);
	for (auto n : stdSharedDOFs->allTargets())
		weight[n] += 1;

	if (!getFetiInfo().bmpc)
		FetiBaseSub::scomm->setTypeSpecificList(SComm::std, std::move(boundNeighbs), std::move(stdSharedDOFs) );
	auto wetSharedDOFs = std::make_unique<Connectivity>(nwneighb, std::move(wetDofPointer), std::move(wetDofs));
	FetiBaseSub::scomm->setTypeSpecificList(SComm::wet, wetNeighbs, std::move(wetSharedDOFs) );


	if (FetiBaseSub::numWIdof) {
		this->wweight.resize(FetiBaseSub::numWIdof);
		for (i = 0; i < FetiBaseSub::numWIdof; ++i) this->wweight[i] = 0.0;
		for (i = 0; i < FetiBaseSub::scomm->lenT(SComm::wet); ++i)
			this->wweight[FetiBaseSub::scomm->boundDofT(SComm::wet, i)] += 1.0;
	}
}

template class FetiSub<double>;
template class FetiSub<std::complex<double>>;
