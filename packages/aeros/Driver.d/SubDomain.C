#include <typeinfo>
#include <cstdio>
#include <algorithm>

#ifdef SUN10
#include <typeinfo.h>
#endif

#include <cmath>
#include <Utils.d/dbg_alloca.h>
#include <set>
#include <climits> //--- UH

#include <Element.d/Element.h>
#include <Utils.d/dofset.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/CuCSparse.h>
#include <Math.d/VectorSet.h>
#include <Math.d/matrix.h>
#include <Math.d/SparseSet.h>
#include <Utils.d/Connectivity.h>
#include <Solvers.d/Rbm.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>
#include <Math.d/MpcSparse.h>
#include <Corotational.d/TemperatureState.h>
#include <Solvers.d/SolverFactory.h>
#include <Driver.d/SubDomain.h>

//#define DEBUG_MPC

extern int verboseFlag;
extern int isFeti3;
extern Domain *domain;
extern int salinasFlag;

template<class Scalar>
GenSubDomain<Scalar>::GenSubDomain(int sn, int lsn) :
	Domain(), // virtual base first
	BaseSub() {
	initialize();
	subNumber = sn;
	localSubNumber = lsn;
}

template<class Scalar>
GenSubDomain<Scalar>::GenSubDomain(Domain &dom, int sn, Connectivity &con, Connectivity &nds, int gn) :
	Domain(dom, con.num(gn), con[gn].data(), nds.num(gn), nds[gn].data()),  // virtual base first
	BaseSub(dom, sn, con, nds, gn) {
	initialize();
}

template<class Scalar>
GenSubDomain<Scalar>::GenSubDomain(Domain &dom, int sn, int nNodes, int *nds, int nElems, int *elems, int gn) :
	Domain(dom, nElems, elems, nNodes, nds), //virtual base first
	BaseSub(dom, sn, nNodes, nds, nElems, elems, gn) {
	initialize();
}

template<class Scalar>
GenSubDomain<Scalar>::GenSubDomain(Domain &dom, int sn, CoordSet *_nodes, Elemset *_elems, int *glNodeNums,
                                   int *glElemNums, int gn) :
	Domain(dom, _elems, _nodes), //virtual base first
	BaseSub(dom, sn, _nodes, _elems, glNodeNums, glElemNums, gn) {
	initialize();
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeAllDisp(Scalar (*xyz)[11], Scalar *u, Scalar (*xyz_loc)[11]) {
	// this is for either Helmholtz or other solver types
	// note: coupledScaling always has a default value of 1.0
	// xyz should be initialized to zero before being passed into this function
	// note: u is already scaled
	int inode, nodeI;
	for (inode = 0; inode < numnodes; ++inode) {
		nodeI = (domain->outFlag) ? domain->nodeTable[glNums[inode]] - 1 : glNums[inode];

		int xLoc = c_dsa->locate(inode, DofSet::Xdisp);
		int xLoc1 = dsa->locate(inode, DofSet::Xdisp);

		if (xLoc >= 0)
			xyz[nodeI][0] = u[xLoc];           // free
		else if (xLoc1 >= 0)
			xyz[nodeI][0] = Bcx(xLoc1);    // constrained

		int yLoc = c_dsa->locate(inode, DofSet::Ydisp);
		int yLoc1 = dsa->locate(inode, DofSet::Ydisp);

		if (yLoc >= 0)
			xyz[nodeI][1] = u[yLoc];
		else if (yLoc1 >= 0)
			xyz[nodeI][1] = Bcx(yLoc1);

		int zLoc = c_dsa->locate(inode, DofSet::Zdisp);
		int zLoc1 = dsa->locate(inode, DofSet::Zdisp);

		if (zLoc >= 0)
			xyz[nodeI][2] = u[zLoc];
		else if (zLoc1 >= 0)
			xyz[nodeI][2] = Bcx(zLoc1);

		int xRot = c_dsa->locate(inode, DofSet::Xrot);
		int xRot1 = dsa->locate(inode, DofSet::Xrot);

		if (xRot >= 0)
			xyz[nodeI][3] = u[xRot];
		else if (xRot1 >= 0)
			xyz[nodeI][3] = Bcx(xRot1);

		int yRot = c_dsa->locate(inode, DofSet::Yrot);
		int yRot1 = dsa->locate(inode, DofSet::Yrot);

		if (yRot >= 0)
			xyz[nodeI][4] = u[yRot];
		else if (yRot1 >= 0)
			xyz[nodeI][4] = Bcx(yRot1);

		int zRot = c_dsa->locate(inode, DofSet::Zrot);
		int zRot1 = dsa->locate(inode, DofSet::Zrot);

		if (zRot >= 0)
			xyz[nodeI][5] = u[zRot];
		else if (zRot1 >= 0)
			xyz[nodeI][5] = Bcx(zRot1);

		int xTemp = c_dsa->locate(inode, DofSet::Temp);
		int xTemp1 = dsa->locate(inode, DofSet::Temp);

		if (xTemp >= 0)
			xyz[nodeI][6] = u[xTemp];
		else if (xTemp1 >= 0)
			xyz[nodeI][6] = Bcx(xTemp1);

		int xHelm = c_dsa->locate(inode, DofSet::Helm);
		int xHelm1 = dsa->locate(inode, DofSet::Helm);

		if (xHelm >= 0)
			xyz[nodeI][7] = u[xHelm];
		else if (xHelm1 >= 0)
			xyz[nodeI][7] = Bcx(xHelm1);

		// transform displacements and rotations (if present) from DOF_FRM to basic coordinates
		// and keep a copy
		if (!domain->solInfo().basicDofCoords && c_dsa->locate(inode, DofSet::LagrangeE) < 0
		    && c_dsa->locate(inode, DofSet::LagrangeI) < 0) {
			if (xyz_loc) for (int j = 0; j < 11; ++j) xyz_loc[nodeI][j] = xyz[nodeI][j];
			bool hasRot = (xRot >= 0 || xRot1 >= 0 || yRot >= 0 || yRot1 >= 0 || zRot >= 0 || zRot1 >= 0);
			transformVectorInv(&(xyz[nodeI][0]), inode, hasRot);
		}
	}
}


template<class Scalar>
void
GenSubDomain<Scalar>::forceContinuity(Scalar *u, Scalar (*xyz)[11])//DofSet::max_known_nonL_dof
{
	int inode;
	for (inode = 0; inode < numnodes; ++inode) {
		int xLoc = c_dsa->locate(inode, DofSet::Xdisp);

		if (xLoc >= 0)
			u[xLoc] = xyz[glNums[inode]][0];

		int yLoc = c_dsa->locate(inode, DofSet::Ydisp);

		if (yLoc >= 0)
			u[yLoc] = xyz[glNums[inode]][1];

		int zLoc = c_dsa->locate(inode, DofSet::Zdisp);

		if (zLoc >= 0)
			u[zLoc] = xyz[glNums[inode]][2];

		int xRot = c_dsa->locate(inode, DofSet::Xrot);

		if (xRot >= 0)
			u[xRot] = xyz[glNums[inode]][3];

		int yRot = c_dsa->locate(inode, DofSet::Yrot);

		if (yRot >= 0)
			u[yRot] = xyz[glNums[inode]][4];

		int zRot = c_dsa->locate(inode, DofSet::Zrot);

		if (zRot >= 0)
			u[zRot] = xyz[glNums[inode]][5];

		int xTemp = c_dsa->locate(inode, DofSet::Temp);

		if (xTemp >= 0)
			u[xTemp] = xyz[glNums[inode]][6];

		int xHelm = c_dsa->locate(inode, DofSet::Helm);

		if (xHelm >= 0)
			u[xHelm] = xyz[glNums[inode]][7];
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeAllVeloc(Scalar (*xyz)[11], Scalar *v, Scalar (*xyz_loc)[11]) {
	// xyz should be initialized to zero before being passed into this function
	int inode, nodeI;
	for (inode = 0; inode < nodes.size(); ++inode) {
		if (nodes[inode] == NULL) continue;
		nodeI = (domain->outFlag) ? domain->nodeTable[glNums[inode]] - 1 : glNums[inode];

		int xLoc = c_dsa->locate(inode, DofSet::Xdisp);
		int xLoc1 = dsa->locate(inode, DofSet::Xdisp);

		if (xLoc >= 0)
			xyz[nodeI][0] = v[xLoc];           // free
		else if (xLoc1 >= 0)
			xyz[nodeI][0] = vcx[xLoc1];        // constrained

		int yLoc = c_dsa->locate(inode, DofSet::Ydisp);
		int yLoc1 = dsa->locate(inode, DofSet::Ydisp);

		if (yLoc >= 0)
			xyz[nodeI][1] = v[yLoc];
		else if (yLoc1 >= 0)
			xyz[nodeI][1] = vcx[yLoc1];

		int zLoc = c_dsa->locate(inode, DofSet::Zdisp);
		int zLoc1 = dsa->locate(inode, DofSet::Zdisp);

		if (zLoc >= 0)
			xyz[nodeI][2] = v[zLoc];
		else if (zLoc1 >= 0)
			xyz[nodeI][2] = vcx[zLoc1];

		int xRot = c_dsa->locate(inode, DofSet::Xrot);
		int xRot1 = dsa->locate(inode, DofSet::Xrot);

		if (xRot >= 0)
			xyz[nodeI][3] = v[xRot];
		else if (xRot1 >= 0)
			xyz[nodeI][3] = vcx[xRot1];

		int yRot = c_dsa->locate(inode, DofSet::Yrot);
		int yRot1 = dsa->locate(inode, DofSet::Yrot);

		if (yRot >= 0)
			xyz[nodeI][4] = v[yRot];
		else if (yRot1 >= 0)
			xyz[nodeI][4] = vcx[yRot1];

		int zRot = c_dsa->locate(inode, DofSet::Zrot);
		int zRot1 = dsa->locate(inode, DofSet::Zrot);

		if (zRot >= 0)
			xyz[nodeI][5] = v[zRot];
		else if (zRot1 >= 0)
			xyz[nodeI][5] = vcx[zRot1];

		int xTemp = c_dsa->locate(inode, DofSet::Temp);
		int xTemp1 = dsa->locate(inode, DofSet::Temp);

		if (xTemp >= 0)
			xyz[nodeI][6] = v[xTemp];
		else if (xTemp1 >= 0)
			xyz[nodeI][6] = vcx[xTemp1];

		int xHelm = c_dsa->locate(inode, DofSet::Helm);
		int xHelm1 = dsa->locate(inode, DofSet::Helm);

		if (xHelm >= 0)
			xyz[nodeI][7] = v[xHelm];
		else if (xHelm1 >= 0)
			xyz[nodeI][7] = vcx[xHelm1];

		// transform velocities and angular velocities (if present) from DOF_FRM to basic coordinates
		if (!domain->solInfo().basicDofCoords && c_dsa->locate(inode, DofSet::LagrangeE) < 0
		    && c_dsa->locate(inode, DofSet::LagrangeI) < 0) {
			if (xyz_loc) for (int j = 0; j < 11; ++j) xyz_loc[nodeI][j] = xyz[nodeI][j];
			bool hasRot = (xRot >= 0 || xRot1 >= 0 || yRot >= 0 || yRot1 >= 0 || zRot >= 0 || zRot1 >= 0);
			transformVectorInv(&(xyz[nodeI][0]), inode, hasRot);
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeAllAccel(Scalar (*xyz)[11], Scalar *a, Scalar (*xyz_loc)[11]) {
	// xyz should be initialized to zero before being passed into this function
	int inode, nodeI;
	for (inode = 0; inode < nodes.size(); ++inode) {
		if (nodes[inode] == NULL) continue;
		nodeI = (domain->outFlag) ? domain->nodeTable[glNums[inode]] - 1 : glNums[inode];

		int xLoc = c_dsa->locate(inode, DofSet::Xdisp);
		int xLoc1 = dsa->locate(inode, DofSet::Xdisp);

		if (xLoc >= 0)
			xyz[nodeI][0] = a[xLoc];           // free
		else if (xLoc1 >= 0)
			xyz[nodeI][0] = acx[xLoc1];        // constrained

		int yLoc = c_dsa->locate(inode, DofSet::Ydisp);
		int yLoc1 = dsa->locate(inode, DofSet::Ydisp);

		if (yLoc >= 0)
			xyz[nodeI][1] = a[yLoc];
		else if (yLoc1 >= 0)
			xyz[nodeI][1] = acx[yLoc1];

		int zLoc = c_dsa->locate(inode, DofSet::Zdisp);
		int zLoc1 = dsa->locate(inode, DofSet::Zdisp);

		if (zLoc >= 0)
			xyz[nodeI][2] = a[zLoc];
		else if (zLoc1 >= 0)
			xyz[nodeI][2] = acx[zLoc1];

		int xRot = c_dsa->locate(inode, DofSet::Xrot);
		int xRot1 = dsa->locate(inode, DofSet::Xrot);

		if (xRot >= 0)
			xyz[nodeI][3] = a[xRot];
		else if (xRot1 >= 0)
			xyz[nodeI][3] = acx[xRot1];

		int yRot = c_dsa->locate(inode, DofSet::Yrot);
		int yRot1 = dsa->locate(inode, DofSet::Yrot);

		if (yRot >= 0)
			xyz[nodeI][4] = a[yRot];
		else if (yRot1 >= 0)
			xyz[nodeI][4] = acx[yRot1];

		int zRot = c_dsa->locate(inode, DofSet::Zrot);
		int zRot1 = dsa->locate(inode, DofSet::Zrot);

		if (zRot >= 0)
			xyz[nodeI][5] = a[zRot];
		else if (zRot1 >= 0)
			xyz[nodeI][5] = acx[zRot1];

		int xTemp = c_dsa->locate(inode, DofSet::Temp);
		int xTemp1 = dsa->locate(inode, DofSet::Temp);

		if (xTemp >= 0)
			xyz[nodeI][6] = a[xTemp];
		else if (xTemp1 >= 0)
			xyz[nodeI][6] = acx[xTemp1];

		int xHelm = c_dsa->locate(inode, DofSet::Helm);
		int xHelm1 = dsa->locate(inode, DofSet::Helm);

		if (xHelm >= 0)
			xyz[nodeI][7] = a[xHelm];
		else if (xHelm1 >= 0)
			xyz[nodeI][7] = acx[xHelm1];

		// transform accelerations and angular accelerations (if present) from DOF_FRM to basic coordinates
		if (!domain->solInfo().basicDofCoords && c_dsa->locate(inode, DofSet::LagrangeE) < 0
		    && c_dsa->locate(inode, DofSet::LagrangeI) < 0) {
			if (xyz_loc) for (int j = 0; j < 11; ++j) xyz_loc[nodeI][j] = xyz[nodeI][j];
			bool hasRot = (xRot >= 0 || xRot1 >= 0 || yRot >= 0 || yRot1 >= 0 || zRot >= 0 || zRot1 >= 0);
			transformVectorInv(&(xyz[nodeI][0]), inode, hasRot);
		}
	}
}

template<class Scalar>
void GenSubDomain<Scalar>::addUserForce(Scalar *extForce, Scalar *usrDefForce) {
	int i;
	for (i = 0; i < claw->numUserForce; ++i) {
		int dof = c_dsa->locate(claw->userForce[i].nnum, 1 << claw->userForce[i].dofnum);
		if (dof > -1)
			extForce[dof] += usrDefForce[locToGlUserForceMap[i]] / static_cast<double>(dofWeight(dof));
	}
}

template<class Scalar>
void GenSubDomain<Scalar>::addCtrl(Scalar *force, Scalar *ctrfrc) {
	int i;
	for (i = 0; i < claw->numActuator; ++i) {
		int dof = c_dsa->locate(claw->actuator[i].nnum, 1 << claw->actuator[i].dofnum);
		if (dof > -1)
			force[dof] += ctrfrc[locToGlActuatorMap[i]] / static_cast<double>(dofWeight(dof));
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeStress(Scalar *locStress, Scalar *locWeight,
                                  Scalar *globStress, Scalar *globWeight, int glNumNodes) {
	int inode, nodeI;
	for (inode = 0; inode < numnodes; ++inode) {
		nodeI = (domain->outFlag) ? domain->nodeTable[glNums[inode]] - 1 : glNums[inode];
		if (nodeI >= glNumNodes) continue;
		globWeight[nodeI] += locWeight[inode];
		globStress[nodeI] += locStress[inode];
	}
}

template<class Scalar>
void GenSubDomain<Scalar>::mergeElemStress(Scalar *locStress, Scalar *globStress, const Connectivity *glElemToNode) {
	const auto &glOffset = glElemToNode->ptr();
	const auto &locOffset = elemToNode->ptr();
	for (int iElem = 0; iElem < numele; iElem++) {
		for (int iNode = 0; iNode < packedEset[iElem]->numNodes(); iNode++) {
			int glOff = glOffset[glElems[iElem]] + iNode;
			int locOff = locOffset[iElem] + iNode;
			globStress[glOff] = locStress[locOff];
		}
	}
}

inline double square(double x) { return x * x; }

inline double square(std::complex<double> x) { return x.real() * x.real() + x.imag() * x.imag(); }

template<class Scalar>
void
GenSubDomain<Scalar>::mergePrimalError(Scalar *error, Scalar *primal) {
	for (int inode = 0; inode < numnodes; ++inode) {
		double nd = 0.0;
		double totP = 0.0;
		int xLoc = c_dsa->locate(inode, DofSet::Xdisp);
		if (xLoc >= 0) {
			totP += square(primal[xLoc] / static_cast<double>(dofWeight(xLoc)));
			nd += 1.0;
		}
		int yLoc = c_dsa->locate(inode, DofSet::Ydisp);
		if (yLoc >= 0) {
			totP += square(primal[yLoc] / static_cast<double>(dofWeight(yLoc)));
			nd += 1.0;
		}
		int zLoc = c_dsa->locate(inode, DofSet::Zdisp);
		if (zLoc >= 0) {
			totP += square(primal[zLoc] / static_cast<double>(dofWeight(zLoc)));
			nd += 1.0;
		}
		if (nd != 0) totP = sqrt(totP / nd);
		error[glNums[inode]] += totP;
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::addDMass(int glNum, int dof, double v) {
	Domain::addDMass(glNum, dof, v);
}

template<class Scalar>
void
GenSubDomain<Scalar>::applySplitting() {
	// adjust discrete masses, forces and mpcs using subdomain multiplicity
	applyDmassSplitting();
	applyForceSplitting();
	applyMpcSplitting();
}

template<class Scalar>
void
GenSubDomain<Scalar>::applyDmassSplitting() {
	// adjust discrete masses using subdomain multiplicity
	// num = number of subdomains touching a dof
	int cdof, num;

	// discrete masses
	DMassData *cmass = firstDiMass;
	while (cmass != 0) {
		if ((cdof = c_dsa->locate(cmass->node, (1 << cmass->dof))) > -1 && (num = weightPlus[cdof]) > 1)
			cmass->diMass /= double(num);
		cmass = cmass->next;
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::applyForceSplitting() {
	// adjust forces using subdomain multiplicity
	// num = number of subdomains touching a dof
	int cdof, num;

	// forces
	for (int i = 0; i < numNeuman; ++i) {
		if ((cdof = c_dsa->locate(nbc[i].nnum, (1 << nbc[i].dofnum))) > -1 && (num = weightPlus[cdof]) > 1)
			nbc[i].val /= double(num);
	}
	for (int i = 0; i < numComplexNeuman; ++i) {
		if ((cdof = c_dsa->locate(cnbc[i].nnum, (1 << cnbc[i].dofnum))) > -1 && (num = weightPlus[cdof]) > 1) {
			cnbc[i].reval /= double(num);
			cnbc[i].imval /= double(num);
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::applyMpcSplitting() {
	// adjust discrete masses, forces and mpcs using subdomain multiplicity
	// num = number of subdomains touching a dof
	int cdof, num;

	auto &mpc_primal = this->mpc_primal;

	// mpcs (NOTE: optional kscaling is done later, hhs is not split)
	if (solInfo().getFetiInfo().mpc_scaling == FetiInfo::tscaling) {
		for (int iMPC = 0; iMPC < numMPC; ++iMPC) { // dual mpcs
			if (this->mpc[iMPC]->type == 2) continue; // bmpc
			for (int i = 0; i < this->mpc[iMPC]->nterms; ++i) {
				if ((cdof = this->mpc[iMPC]->terms[i].cdof) > -1 && (num = weightPlus[cdof]) > 1)
					this->mpc[iMPC]->terms[i].coef /= double(num);
			}
		}
	}
	// XXXX kscaling currently not supported for primal mpcs
	for (int iMPC = 0; iMPC < numMPC_primal; ++iMPC) { // primal mpcs
		for (int i = 0; i < mpc_primal[iMPC]->nterms; ++i) {
			if ((cdof = mpc_primal[iMPC]->terms[i].cdof) > -1 && (num = weightPlus[cdof]) > 1)
				mpc_primal[iMPC]->terms[i].coef /= double(num);
		}
	}

}

template<class Scalar>
void
GenSubDomain<Scalar>::sendNode(Scalar (*subvec)[11], FSCommPattern<Scalar> *pat) {
	for (int i = 0; i < scomm->numNeighb; ++i) {
		FSSubRecInfo<Scalar> sInfo = pat->getSendBuffer(subNumber, scomm->subNums[i]);
		for (int j = 0; j < scomm->sharedNodes->num(i); ++j) {
			int n = (*(scomm->sharedNodes))[i][j];
			for (int k = 0; k < 11; ++k) sInfo.data[11 * j + k] = subvec[n][k];
		}
	}
}

#ifndef MAXABS
#define MAXABS(X, Y) ((ScalarTypes::norm(X) > ScalarTypes::norm(Y)) ? X : Y)
#endif

template<class Scalar>
void
GenSubDomain<Scalar>::collectNode(Scalar (*subvec)[11], FSCommPattern<Scalar> *pat) {
	// use the one with the largest absolute value
	for (int i = 0; i < scomm->numNeighb; ++i) {
		FSSubRecInfo<Scalar> rInfo = pat->recData(scomm->subNums[i], subNumber);
		for (int j = 0; j < scomm->sharedNodes->num(i); ++j) {
			int n = (*(scomm->sharedNodes))[i][j];
			for (int k = 0; k < 11; ++k) subvec[n][k] = MAXABS(subvec[n][k], rInfo.data[11 * j + k]);
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::extractAndSendInterf(const Scalar *subvec, FSCommPattern<Scalar> *pat) const {
	for (int i = 0; i < scomm->numNeighb; ++i) {
		FSSubRecInfo<Scalar> sInfo = pat->getSendBuffer(subNumber, scomm->subNums[i]);
		for (int j = 0; j < scomm->sharedDOFsPlus->num(i); ++j) {
			sInfo.data[j] = subvec[(*scomm->sharedDOFsPlus)[i][j]];
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::assembleInterf(Scalar *subvec, FSCommPattern<Scalar> *pat) const {
	for (int i = 0; i < scomm->numNeighb; ++i) {
		FSSubRecInfo<Scalar> rInfo = pat->recData(scomm->subNums[i], subNumber);
		for (int j = 0; j < scomm->sharedDOFsPlus->num(i); ++j) {
			subvec[(*scomm->sharedDOFsPlus)[i][j]] += rInfo.data[j];
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::assembleInterfInvert(Scalar *subvec, FSCommPattern<Scalar> *pat) const {
	for (int i = 0; i < numUncon(); ++i) subvec[i] = 1.0 / subvec[i];
	for (int i = 0; i < scomm->numNeighb; ++i) {
		FSSubRecInfo<Scalar> rInfo = pat->recData(scomm->subNums[i], subNumber);
		for (int j = 0; j < scomm->sharedDOFsPlus->num(i); ++j) {
			subvec[(*scomm->sharedDOFsPlus)[i][j]] += 1.0 / rInfo.data[j];
		}
	}
	for (int i = 0; i < numUncon(); ++i) subvec[i] = 1.0 / subvec[i];
}

template<class Scalar>
void
GenSubDomain<Scalar>::sendDeltaF(const Scalar *deltaF, FSCommPattern<Scalar> *vPat) {
	auto &deltaFmpc = this->deltaFmpc;
	int iDof = 0;
	for (int i = 0; i < scomm->numT(SComm::all); ++i) {
		FSSubRecInfo<Scalar> sInfo = vPat->getSendBuffer(subNumber, scomm->neighbT(SComm::all, i));
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
GenSubDomain<Scalar>::collectAndDotDeltaF(Scalar *deltaF, FSCommPattern<Scalar> *vPat) {
	// if there are more than 2 subdomains sharing a mpc define the norm
	// as (f1 - f2 - f3)^2 = f1^2 + f2^2 + f3^2 - 2f1f2 - 2f1f3 + 2f2f3 --> currently implemented
	auto &deltaFmpc = this->deltaFmpc;

	Scalar dot = 0;
	int i, iSub, jDof;

	if (deltaF) {
		for (i = 0; i < localLen(); ++i) {
			double dPrScal = 1.0 / this->densProjCoeff(i);
			dot += dPrScal * dPrScal * deltaF[i] * ScalarTypes::conj(deltaF[i]);
		}
	}

	for (i = 0; i < numMPC; ++i)
		dot += deltaFmpc[i] * ScalarTypes::conj(deltaFmpc[i]);

	for (i = 0; i < numWIdof; ++i) //HB ... to be checked ...
		dot += this->deltaFwi[i] * ScalarTypes::conj(this->deltaFwi[i]) / this->wweight[i];

	int nbdofs = 0;
	for (iSub = 0; iSub < scomm->numT(SComm::all); ++iSub) {
		FSSubRecInfo<Scalar> rInfo = vPat->recData(scomm->neighbT(SComm::all, iSub), subNumber);
		for (jDof = 0; jDof < scomm->lenT(SComm::all, iSub); ++jDof) {
			int bdof = scomm->boundDofT(SComm::all, iSub, jDof);
			switch (boundDofFlag[nbdofs]) {
				case 0:
					if (deltaF) {
						double dPrScal = 1.0 / this->densProjCoeff(bdof);
						dot += dPrScal * dPrScal * deltaF[bdof] * ScalarTypes::conj(rInfo.data[jDof]);
					}
					break;
					// do nothing for case 1 (wet interface)
				case 2: { // dual mpc
					if (subNumber != scomm->neighbT(SComm::all, iSub)) {
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
void GenSubDomain<Scalar>::extractControlData(Scalar *disp, Scalar *vel,
                                              Scalar *acc, Scalar *ctrdsp,
                                              Scalar *ctrvel, Scalar *ctracc) {
	for (int i = 0; i < claw->numSensor; ++i) {
		int dof = c_dsa->locate(claw->sensor[i].nnum, 1 << claw->sensor[i].dofnum);
		int gi = locToGlSensorMap[i]; // local index
		if (dof >= 0) { // free
			ctrdsp[gi] = disp[dof];
			ctrvel[gi] = vel[dof];
			ctracc[gi] = acc[dof];
		} else { // either constrained or non-existant
			int dof2 = dsa->locate(claw->sensor[i].nnum, 1 << claw->sensor[i].dofnum);
			if (dof2 >= 0) { // constrained
				ctrdsp[gi] = bcx[dof2];
				ctrvel[gi] = vcx[dof2];
				ctracc[gi] = acx[dof2];
			}
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::constructKrc() {
	this->Src = std::make_unique<GenSparseSet<Scalar>>();
	if (numCRNdof) {
		this->Krc = std::make_unique<GenCuCSparse<Scalar>>(nodeToNode.get(), dsa, cornerMap, cc_dsa->getUnconstrNum());
		this->Src->addSparseMatrix(this->Krc);
	}
	setMpcSparseMatrix();
}

template<>
void
GenSubDomain<DComplex>::getSRMult(const DComplex *lvec, const DComplex *interfvec, int nRBM,
                                  const double *locRBMs, DComplex *alpha) const;

template<>
void
GenSubDomain<double>::getSRMult(const double *lvec, const double *interfvec, int nRBM,
                                const double *locRBMs, double *alpha) const;

template<class Scalar>
void
GenSubDomain<Scalar>::multFi(GenSolver<Scalar> *s, Scalar *u, Scalar *Fiu) {
	// multFi is never called in DP. Otherwise it will crash in contact

	int iDof;
	int ndof = localLen();
	Scalar *iDisp = (Scalar *) dbg_alloca(sizeof(Scalar) * ndof);
	for (iDof = 0; iDof < ndof; ++iDof)
		iDisp[iDof] = 0.0;

	// Add the interface vector contribution to local vector
	for (iDof = 0; iDof < totalInterfSize; ++iDof)
		iDisp[allBoundDofs[iDof]] += u[iDof];

	// solve for local vector
	if (s) s->reSolve(iDisp);

	// redistribute solution to the interface
	for (iDof = 0; iDof < totalInterfSize; ++iDof)
		Fiu[iDof] = iDisp[allBoundDofs[iDof]];

}

template<class Scalar>
void
GenSubDomain<Scalar>::renumberElements() {
	for (int i = 0; i < numele; ++i)
		packedEset[i]->renum(glToLocalNode);
}

template<class Scalar>
void
GenSubDomain<Scalar>::renumberElementsGlobal() {
	for (int i = 0; i < numele; ++i)
		packedEset[i]->renum(glNums.data());

	for (int i = 0; i < numNeum; ++i) {
		neum[i]->renum(glNums.data());
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::renumberSharedNodes() {
	auto &allC = scomm->sharedNodes->tgt();
	for (int i = 0; i < scomm->sharedNodes->numConnect(); ++i) {
		int gi = allC[i];
		//if(globalToLocal(gi) < 0) std::cerr << "error in GenSubDomain::renumberSharedNodes() \n";
		allC[i] = globalToLocal(gi);
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::renumberDirichlet() {
	for (int i = 0; i < numDirichlet; ++i)
		dbc[i].nnum = glToLocalNode[dbc[i].nnum];
}

template<class Scalar>
void
GenSubDomain<Scalar>::renumberBCsEtc() {
	int i;

	for (i = 0; i < numDirichlet; ++i)
		dbc[i].nnum = glToLocalNode[dbc[i].nnum];

	for (i = 0; i < numNeuman; ++i)
		nbc[i].nnum = glToLocalNode[nbc[i].nnum];

	for (i = 0; i < numIDis; ++i)
		iDis[i].nnum = glToLocalNode[iDis[i].nnum];

	for (i = 0; i < numIDis6; ++i)
		iDis6[i].nnum = glToLocalNode[iDis6[i].nnum];

	for (i = 0; i < numIVel; ++i)
		iVel[i].nnum = glToLocalNode[iVel[i].nnum];

	for (i = 0; i < numComplexDirichlet; ++i)
		cdbc[i].nnum = glToLocalNode[cdbc[i].nnum];

	for (i = 0; i < numComplexNeuman; ++i)
		cnbc[i].nnum = glToLocalNode[cnbc[i].nnum];

	for (i = 0; i < numScatter; ++i)
		scatter[i]->renum(glToLocalNode);

	for (i = 0; i < numWet; ++i)
		wet[i]->renum(glToLocalNode);

	if (sinfo.ATDARBFlag == -2.0)
		for (i = 0; i < numSommer; ++i)
			sommer[i]->renum(glToLocalNode);

	for (i = 0; i < numNeum; ++i)
		neum[i]->renum(glToLocalNode);

	for (i = 0; i < numSBoundNodes; i++) {
		int tmp = glToLocalNode[sBoundNodes[i]];
		sBoundNodes[i] = tmp;
	}

	DMassData *cmass = firstDiMass;
	while (cmass != 0) {
		cmass->node = glToLocalNode[cmass->node];
		cmass = cmass->next;
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::renumberControlLaw() {
	int i;
	if (claw) {
		for (i = 0; i < claw->numSensor; i++)
			claw->sensor[i].nnum = glToLocalNode[claw->sensor[i].nnum];

		for (i = 0; i < claw->numActuator; i++)
			claw->actuator[i].nnum = glToLocalNode[claw->actuator[i].nnum];

		for (i = 0; i < claw->numUserDisp; i++)
			claw->userDisp[i].nnum = glToLocalNode[claw->userDisp[i].nnum];

		for (i = 0; i < claw->numUserForce; i++)
			claw->userForce[i].nnum = glToLocalNode[claw->userForce[i].nnum];
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::renumberMPCs() {
	auto &mpc = this->mpc;
	auto &mpc_primal = this->mpc_primal;

	for (int iMPC = 0; iMPC < numMPC; ++iMPC) {
		for (int i = 0; i < mpc[iMPC]->nterms; ++i)
			mpc[iMPC]->terms[i].nnum = globalToLocal(mpc[iMPC]->terms[i].nnum);
	}
	for (int iMPC = 0; iMPC < numMPC_primal; ++iMPC) {
		for (int i = 0; i < mpc_primal[iMPC]->nterms; ++i)
			mpc_primal[iMPC]->terms[i].nnum = globalToLocal(mpc_primal[iMPC]->terms[i].nnum);
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::computeElementForce(int fileNumber, Scalar *u, int forceIndex, Scalar *elemForce) {
	if (elemToNode == 0) elemToNode = new Connectivity(packedEset.asSet());
	if (elDisp == 0) elDisp = new Vector(maxNumDOFs, 0.0);
	double *nodalTemperatures = getNodalTemperatures();
	int *nodeNumbers = new int[maxNumNodes];
	Vector elemNodeTemps(maxNumNodes, 0.0);
	Vector elForce(maxNumNodes, 0.0);

	int iele, k;
	for (iele = 0; iele < numele; ++iele) {

		packedEset[iele]->nodes(nodeNumbers);
		int NodesPerElement = packedEset[iele]->numNodes();

		for (k = 0; k < allDOFs->num(iele); ++k) {
			int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
			if (cn >= 0)
				(*elDisp)[k] = ScalarTypes::Real(u[cn]);
			else
				(*elDisp)[k] = 0.0;
		}

		if (packedEset[iele]->getProperty()) {
			for (int iNode = 0; iNode < NodesPerElement; ++iNode) {
				if (!nodalTemperatures || nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
					elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
				else
					elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
			}
		}

		// transform displacements from DOF_FRM to basic coordinates
		transformVectorInv(*elDisp, iele);

		if (geoSource->getOutputInfo()[fileNumber].oframe == OutputInfo::Local) {
			Vector fx(NodesPerElement, 0.0), fy(NodesPerElement, 0.0), fz(NodesPerElement, 0.0);
			if (forceIndex == INX || forceIndex == INY || forceIndex == INZ) {
				packedEset[iele]->getIntrnForce(fx, nodes, elDisp->data(), INX, elemNodeTemps.data());
				packedEset[iele]->getIntrnForce(fy, nodes, elDisp->data(), INY, elemNodeTemps.data());
				packedEset[iele]->getIntrnForce(fz, nodes, elDisp->data(), INZ, elemNodeTemps.data());
			} else {
				packedEset[iele]->getIntrnForce(fx, nodes, elDisp->data(), AXM, elemNodeTemps.data());
				packedEset[iele]->getIntrnForce(fy, nodes, elDisp->data(), AYM, elemNodeTemps.data());
				packedEset[iele]->getIntrnForce(fz, nodes, elDisp->data(), AZM, elemNodeTemps.data());
			}
			Vector f(3 * NodesPerElement);
			for (int iNode = 0; iNode < NodesPerElement; iNode++) {
				double data[3] = {fx[iNode], fy[iNode], fz[iNode]};
				transformVector(data, nodeNumbers[iNode], false);
				switch (forceIndex) {
					case INX:
					case AXM:
						elstress[iNode] = data[0];
						break;
					case INY:
					case AYM:
						elstress[iNode] = data[1];
						break;
					case INZ:
					case AZM:
						elstress[iNode] = data[2];
						break;
				}
			}
		} else {
			packedEset[iele]->getIntrnForce(elForce, nodes, elDisp->data(), forceIndex, elemNodeTemps.data());
		}

		for (k = 0; k < packedEset[iele]->numNodes(); ++k)
			elemForce[elemToNode->offset(iele) + k] = elForce[k];
	}
	delete[] nodeNumbers;
}

template<class Scalar>
void
GenSubDomain<Scalar>::computeStressStrain(int fileNumber,
                                          Scalar *u, int stressIndex,
                                          Scalar *glStress, Scalar *glWeight) {
	OutputInfo *oinfo = geoSource->getOutputInfo();
	if (elemToNode == 0) elemToNode = new Connectivity(packedEset.asSet());

	// avgnum = 2 --> do not include stress/strain of bar/beam element in averaging
	// avgnum = 3 --> only include elements whose nodes are in the ouputGroup in averaging
	int avgnum = oinfo[fileNumber].averageFlg;

	// ylayer and zlayer are needed when calculating the axial stress/strain in a beam element
	double ylayer = oinfo[fileNumber].ylayer;
	double zlayer = oinfo[fileNumber].zlayer;
	OutputInfo::FrameType oframe = oinfo[fileNumber].oframe;

	int *nodeNumbers = new int[maxNumNodes];
	Vector elemNodeTemps(maxNumNodes);
	elemNodeTemps.zero();
	double *nodalTemperatures = getNodalTemperatures();

	int k;
	int surface = oinfo[fileNumber].surface;

	int iele, iNode;
	GenVector<Scalar> *elstress = new GenVector<Scalar>(maxNumNodes);
	GenVector<double> *elweight = new GenVector<double>(maxNumNodes);
	GenVector<Scalar> *elDisp = new GenVector<Scalar>(maxNumDOFs);
	GenFullM<Scalar> *p_elstress = 0;
	if (oframe == OutputInfo::Local) p_elstress = new GenFullM<Scalar>(maxNumNodes, 9);

	for (iele = 0; iele < numele; ++iele) {

		// Don't do anything if element is a phantom
		if (packedEset[iele]->isPhantomElement()) continue;

		// Don't include beams or bars in the averaging if nodalpartial (avgnum = 2) is requested
		if ((avgnum == 2 && packedEset[iele]->getElementType() == 6) ||
		    (avgnum == 2 && packedEset[iele]->getElementType() == 7) ||
		    (avgnum == 2 && packedEset[iele]->getElementType() == 1))
			continue;

		int NodesPerElement = packedEset[iele]->numNodes();
		packedEset[iele]->nodes(nodeNumbers);

		// Don't include elements with one or more nodes not in the group if nodalpartialgroup (avgnum = 3) is requested
		if (avgnum == 3) {
			int groupId = oinfo[fileNumber].groupNumber;
			if (groupId > 0) {
				std::set<int> &groupNodes = geoSource->getNodeGroup(groupId);
				std::set<int>::iterator it;
				for (iNode = 0; iNode < NodesPerElement; ++iNode)
					if ((it = groupNodes.find(nodeNumbers[iNode])) == groupNodes.end()) break;
				if (it == groupNodes.end()) continue;
			}
		}

		elstress->zero();
		elweight->zero();
		elDisp->zero();

		if (!isFluid(iele)) {

			for (k = 0; k < allDOFs->num(iele); ++k) {
				int cn = c_dsa->getRCN((*allDOFs)[iele][k]);
				if (cn >= 0)
					(*elDisp)[k] = u[cn];
				else
					(*elDisp)[k] = Bcx((*allDOFs)[iele][k]);
			}

			if (packedEset[iele]->getProperty()) {
				for (iNode = 0; iNode < NodesPerElement; ++iNode) {
					if (!nodalTemperatures || nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
						elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
					else
						elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
				}

				// transform displacements from DOF_FRM to basic coordinates
				transformVectorInv(*elDisp, iele);

				// transform non-invariant stresses/strains from basic frame to DOF_FRM
				if (oframe == OutputInfo::Local &&
				    ((stressIndex >= 0 && stressIndex <= 5) || (stressIndex >= 7 && stressIndex <= 12))) {

					// FIRST, CALCULATE STRESS/STRAIN TENSOR FOR EACH NODE OF THE ELEMENT
					p_elstress->zero();
					int strInd = (stressIndex >= 0 && stressIndex <= 5) ? 0 : 1;
					packedEset[iele]->getAllStress(*p_elstress, *elweight, nodes,
					                               *elDisp, strInd, surface,
					                               elemNodeTemps.data());

					// second, transform stress/strain tensor to nodal frame coordinates
					transformStressStrain(*p_elstress, iele);

					// third, extract the requested stress/strain value from the stress/strain tensor
					for (iNode = 0; iNode < NodesPerElement; ++iNode) {
						if (strInd == 0)
							(*elstress)[iNode] = (*p_elstress)[iNode][stressIndex];
						else
							(*elstress)[iNode] = (*p_elstress)[iNode][stressIndex - 7];
					}

				} else {
					packedEset[iele]->getVonMises(*elstress, *elweight, nodes,
					                              *elDisp, stressIndex, surface,
					                              elemNodeTemps.data(), ylayer, zlayer, avgnum);
				}
			}
		}

		if (glWeight) {
			for (k = 0; k < NodesPerElement; ++k) {
#ifdef DISTRIBUTED
				glStress[nodeNumbers[k]] += (*elstress)[k];
				glWeight[nodeNumbers[k]] += (*elweight)[k];
#else
																																		glStress[(*elemToNode)[iele][k]] += (*elstress)[k];
        glWeight[(*elemToNode)[iele][k]] += (*elweight)[k];
#endif
			}
		} else  // non-averaged stresses
			for (k = 0; k < packedEset[iele]->numNodes(); ++k)
				glStress[elemToNode->offset(iele) + k] = (*elstress)[k];

	}

	delete[] nodeNumbers;
	delete elstress;
	delete elweight;
	delete elDisp;
	if (p_elstress) delete p_elstress;
}

// -----------------------------
//
// Nonlinear Subdomain functions
//
// -----------------------------

template<class Scalar>
void
GenSubDomain<Scalar>::computeElementForce(GeomState *gs, Corotator **allCorot,
                                          int fileNumber, int forceIndex, Scalar *elemForce)
{
  if(elemToNode==0) elemToNode = new Connectivity(packedEset.asSet());
  if(elDisp == 0) elDisp = new Vector(maxNumDOFs, 0.0);
  double *nodalTemperatures = getNodalTemperatures();
  int *nodeNumbers = new int[maxNumNodes];
  Vector elemNodeTemps(maxNumNodes, 0.0);
  Vector elForce(maxNumNodes, 0.0);

  int iele,k,flag;
  for(iele=0; iele<numele; ++iele) {

    packedEset[iele]->nodes(nodeNumbers);
    int NodesPerElement = packedEset[iele]->numNodes();

    allCorot[iele]->extractDeformations(*gs, nodes, elDisp->data(), flag);

    if(packedEset[iele]->getProperty()) {
      for(int iNode=0; iNode<NodesPerElement; ++iNode) {
        if(!nodalTemperatures || nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
          elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
        else
          elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
      }
    }

    // transform displacements from DOF_FRM to basic coordinates
    transformVectorInv(*elDisp, iele);

    if(geoSource->getOutputInfo()[fileNumber].oframe == OutputInfo::Local) {
      Vector fx(NodesPerElement,0.0), fy(NodesPerElement,0.0), fz(NodesPerElement,0.0);
      if(forceIndex == INX || forceIndex == INY || forceIndex == INZ) {
        packedEset[iele]->getIntrnForce(fx,nodes,elDisp->data(),INX,elemNodeTemps.data());
        packedEset[iele]->getIntrnForce(fy,nodes,elDisp->data(),INY,elemNodeTemps.data());
        packedEset[iele]->getIntrnForce(fz,nodes,elDisp->data(),INZ,elemNodeTemps.data());
      }
      else {
        packedEset[iele]->getIntrnForce(fx,nodes,elDisp->data(),AXM,elemNodeTemps.data());
        packedEset[iele]->getIntrnForce(fy,nodes,elDisp->data(),AYM,elemNodeTemps.data());
        packedEset[iele]->getIntrnForce(fz,nodes,elDisp->data(),AZM,elemNodeTemps.data());
      }
      Vector f(3*NodesPerElement);
      for(int iNode=0; iNode<NodesPerElement; iNode++) {
        double data[3] = { fx[iNode], fy[iNode], fz[iNode] };
        transformVector(data, nodeNumbers[iNode], false);
        switch(forceIndex) {
          case INX: case AXM: elstress[iNode] = data[0]; break;
          case INY: case AYM: elstress[iNode] = data[1]; break;
          case INZ: case AZM: elstress[iNode] = data[2]; break;
        }
      }
    }
    else {
      packedEset[iele]->getIntrnForce(elForce, nodes, elDisp->data(), forceIndex, elemNodeTemps.data());
    }

    for(k = 0; k < packedEset[iele]->numNodes(); ++k)
      elemForce[elemToNode->offset(iele) + k] = elForce[k];
  }
  delete [] nodeNumbers;
}

template<class Scalar>
void
GenSubDomain<Scalar>::computeStressStrain(GeomState *gs, Corotator **allCorot,
                                          int fileNumber, int stressIndex, Scalar *glStress, Scalar *glWeight,
                                          GeomState *refState) {
	OutputInfo *oinfo = geoSource->getOutputInfo();

	if (elemToNode == 0) elemToNode = new Connectivity(packedEset.asSet());

	// avgnum = 2 --> do not include stress/strain of bar/beam element in averaging
	// avgnum = 3 --> only include elements whose nodes are in the ouputGroup in averaging
	int avgnum = oinfo[fileNumber].averageFlg;

	// ylayer and zlayer are needed when calculating the axial stress/strain in a beam element
	double ylayer = oinfo[fileNumber].ylayer;
	double zlayer = oinfo[fileNumber].zlayer;
	OutputInfo::FrameType oframe = oinfo[fileNumber].oframe;

	int *nodeNumbers = new int[maxNumNodes];
	Vector elemNodeTemps(maxNumNodes);
	elemNodeTemps.zero();
	double *nodalTemperatures = getNodalTemperatures();

	int k;
	int surface = oinfo[fileNumber].surface;

	int iele, iNode;
	GenVector<Scalar> *elstress = new GenVector<Scalar>(maxNumNodes);
	GenVector<double> *elweight = new GenVector<double>(maxNumNodes);
	GenVector<Scalar> *elDisp = new GenVector<Scalar>(maxNumDOFs);
	GenFullM<Scalar> *p_elstress = 0;
	if (oframe == OutputInfo::Local) p_elstress = new GenFullM<Scalar>(maxNumNodes, 9);

	int flag;

	for (iele = 0; iele < numele; ++iele) {

		// Don't do anything if element is a phantom
		if (packedEset[iele]->isPhantomElement() || packedEset[iele]->isMpcElement()) continue;

		// Don't include beams or bars in the averaging if nodalpartial (avgnum = 2) is requested
		if ((avgnum == 2 && packedEset[iele]->getElementType() == 6) ||
		    (avgnum == 2 && packedEset[iele]->getElementType() == 7) ||
		    (avgnum == 2 && packedEset[iele]->getElementType() == 1))
			continue;

		int NodesPerElement = elemToNode->num(iele);
		packedEset[iele]->nodes(nodeNumbers);

		// Don't include elements with one or more nodes not in the group if nodalpartialgroup (avgnum = 3) is requested
		if (avgnum == 3) {
			int groupId = oinfo[fileNumber].groupNumber;
			if (groupId > 0) {
				std::set<int> &groupNodes = geoSource->getNodeGroup(groupId);
				std::set<int>::iterator it;
				for (iNode = 0; iNode < NodesPerElement; ++iNode)
					if ((it = groupNodes.find(nodeNumbers[iNode])) == groupNodes.end()) break;
				if (it == groupNodes.end()) continue;
			}
		}

		elstress->zero();
		elweight->zero();
		elDisp->zero();

		allCorot[iele]->extractDeformations(*gs, nodes, elDisp->data(), flag);

		if (flag == 1 && packedEset[iele]->getProperty()) {
			for (iNode = 0; iNode < NodesPerElement; ++iNode) {
				if (!nodalTemperatures || nodalTemperatures[nodeNumbers[iNode]] == defaultTemp)
					elemNodeTemps[iNode] = packedEset[iele]->getProperty()->Ta;
				else
					elemNodeTemps[iNode] = nodalTemperatures[nodeNumbers[iNode]];
			}
		}

		if (oframe == OutputInfo::Local &&
		    ((stressIndex >= 0 && stressIndex <= 5) || (stressIndex >= 7 && stressIndex <= 12))
		    && (flag == 1 || flag == 2)) { // transform non-invariant stresses/strains from basic frame to DOF_FRM

			// First, calculate stress/strain tensor for each node of the element
			p_elstress->zero();
			int strInd = (stressIndex >= 0 && stressIndex <= 5) ? 0 : 1;
			if (flag == 1) {
				// USE LINEAR STRESS ROUTINE
				packedEset[iele]->getAllStress(*p_elstress, *elweight, nodes,
				                               *elDisp, strInd, surface,
				                               elemNodeTemps.data());
			} else {
				// USE NON-LINEAR STRESS ROUTINE
				// note: in this case the element nodal temperatures are extracted from geomState inside the function
				allCorot[iele]->getNLAllStress(*p_elstress, *elweight, *gs,
				                               refState, nodes, strInd, surface);
			}

			// Second, transform stress/strain tensor to nodal frame coordinates
			transformStressStrain(*p_elstress, iele);

			// Third, extract the requested stress/strain value from the stress/strain tensor
			for (iNode = 0; iNode < NodesPerElement; ++iNode) {
				if (strInd == 0)
					(*elstress)[iNode] = (*p_elstress)[iNode][stressIndex];
				else
					(*elstress)[iNode] = (*p_elstress)[iNode][stressIndex - 7];
			}
		} else {
			if (flag == 1) {
				// USE LINEAR STRESS ROUTINE
				packedEset[iele]->getVonMises(*elstress, *elweight, nodes,
				                              *elDisp, stressIndex, surface,
				                              elemNodeTemps.data(), ylayer, zlayer, avgnum);
			} else if (flag == 2) {
				// USE NON-LINEAR STRESS ROUTINE
				// note: in this case the element nodal temperatures are extracted from geomState inside the function
				allCorot[iele]->getNLVonMises(*elstress, *elweight, *gs, refState,
				                              nodes, stressIndex, surface, ylayer, zlayer, avgnum);
			} else {
				// NO STRESS RECOVERY
				elstress->zero();
				elweight->zero();
			}
		}

		if (glWeight) {
			for (k = 0; k < NodesPerElement; ++k) {
#ifdef DISTRIBUTED
				glStress[nodeNumbers[k]] += (*elstress)[k];
				glWeight[nodeNumbers[k]] += (*elweight)[k];
#else
																																		glStress[(*elemToNode)[iele][k]] += (*elstress)[k];
        glWeight[(*elemToNode)[iele][k]] += (*elweight)[k];
#endif
			}
		} else  // non-averaged stresses
			for (k = 0; k < packedEset[iele]->numNodes(); ++k)
				glStress[elemToNode->offset(iele) + k] = (*elstress)[k];

	}
	delete elstress;
	delete elweight;
	delete elDisp;
	delete[] nodeNumbers;
	if (p_elstress) delete p_elstress;
}

template<class Scalar>
void
GenSubDomain<Scalar>::reBuildKbb(FullSquareMatrix *kel) {
	// Zero sparse matrices
	if (this->Kcc) this->Kcc->zero();
	if (this->Krc) this->Krc->zeroAll();
	if (this->Kbb) this->Kbb->zeroAll();
	if (this->Kib) this->Kib->zeroAll();
	if (this->KiiSparse) this->KiiSparse->zeroAll();

	// Assemble new subdomain sparse matrices
	int iele;
	for (iele = 0; iele < numele; ++iele) {
		if (this->Kcc) this->Kcc->add(kel[iele], (*allDOFs)[iele]);
		if (this->Krc) this->Krc->add(kel[iele], (*allDOFs)[iele].data());
		if (this->Kbb) this->Kbb->add(kel[iele], (*allDOFs)[iele].data());
		if (this->Kib) this->Kib->add(kel[iele], (*allDOFs)[iele].data());
		if (this->KiiSparse) this->KiiSparse->add(kel[iele], (*allDOFs)[iele].data());
	}

	// Factor Kii if necessary (used for dirichlet preconditioner)
	if (this->KiiSolver) this->KiiSolver->factor();
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeDisp(Scalar (*xyz)[11], GeomState *u, Scalar (*xyz_loc)[11]) //DOfSet::max_known_nonL_dof
{
	// Loop over all nodes, Subtract the initial configuration from the
	// Current configuration to get the displacement
	bool IsTemperatureState = (dynamic_cast<TemperatureState *>(u) != NULL);
	int i, nodeI;
	for (i = 0; i < numnodes; ++i) {
		if (!nodes[i] || i >= u->numNodes()) continue;
		nodeI = (domain->outFlag) ? domain->nodeTable[glNums[i]] - 1 : glNums[i];

		if (IsTemperatureState) {
			xyz[nodeI][0] = (*u)[i].x;
			continue;
		}
		xyz[nodeI][0] = ((*u)[i].x - nodes[i]->x);
		xyz[nodeI][1] = ((*u)[i].y - nodes[i]->y);
		xyz[nodeI][2] = ((*u)[i].z - nodes[i]->z);
		double rot[3];
		mat_to_vec((*u)[i].R, rot);
		xyz[nodeI][3] = rot[0];
		xyz[nodeI][4] = rot[1];
		xyz[nodeI][5] = rot[2];

		if (xyz_loc) {
			for (int j = 0; j < 6; ++j) xyz_loc[nodeI][j] = xyz[nodeI][j];
			transformVector(xyz_loc[nodeI], i, true);
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeDistributedNLDisp(Scalar (*xyz)[11], GeomState *u,
                                             Scalar (*xyz_loc)[11]) //DofSet::max_known_nonL_dof
{
	// Loop over all nodes, Subtract the initial configuration from the
	// Current configuration to get the displacement
	bool IsTemperatureState = (dynamic_cast<TemperatureState *>(u) != NULL);
	int i;
	for (i = 0; i < numnodes; ++i) {
		if (!nodes[i]) continue;
		if (IsTemperatureState) {
			xyz[i][0] = (*u)[i].x;
			continue;
		}
		xyz[i][0] = ((*u)[i].x - nodes[i]->x);
		xyz[i][1] = ((*u)[i].y - nodes[i]->y);
		xyz[i][2] = ((*u)[i].z - nodes[i]->z);
		double rot[3];
		mat_to_vec((*u)[i].R, rot);
		xyz[i][3] = rot[0];
		xyz[i][4] = rot[1];
		xyz[i][5] = rot[2];

		if (xyz_loc) {
			for (int j = 0; j < 6; ++j) xyz_loc[i][j] = xyz[i][j];
			transformVector(xyz_loc[i], i, true);
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeForces(Scalar (*mergedF)[6], Scalar *subF) {
	for (int inode = 0; inode < numnodes; ++inode) {
		int nodeI = (domain->outFlag) ? domain->nodeTable[glNums[inode]] - 1 : glNums[inode];
		for (int jdof = 0; jdof < 6; ++jdof) {
			int cdof = c_dsa->locate(inode, 1 << jdof);
			if (cdof >= 0)
				mergedF[nodeI][jdof] += subF[cdof]; // free: assemble into global array
			else
				mergedF[nodeI][jdof] = 0.0;         // constrained or doesn't exist
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeReactions(Scalar (*mergedF)[11], Scalar *subF) {
	// TODO: just loop over dbc, cdbc
	for (int inode = 0; inode < numnodes; ++inode) {
		int nodeI = (domain->outFlag) ? domain->nodeTable[glNums[inode]] - 1 : glNums[inode];
		for (int jdof = 0; jdof < 11; ++jdof) {
			int dof = dsa->locate(inode, 1 << jdof);
			if (dof >= 0) {
				int cdof = c_dsa->invRCN(dof);
				if (cdof >= 0) mergedF[nodeI][jdof] += subF[cdof]; // constrained
			}
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeDistributedForces(Scalar (*mergedF)[6], Scalar *subF) {
	for (int inode = 0; inode < numnodes; ++inode) {
		for (int jdof = 0; jdof < 6; ++jdof) {
			int cdof = c_dsa->locate(inode, 1 << jdof);
			if (cdof >= 0)
				mergedF[inode][jdof] = subF[cdof];  // free
			else
				mergedF[inode][jdof] = 0.0;         // constrained or doesn't exist
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::mergeDistributedReactions(Scalar (*mergedF)[11], Scalar *subF) {
	// TODO: just loop over dbc, cdbc
	for (int inode = 0; inode < numnodes; ++inode) {
		for (int jdof = 0; jdof < 11; ++jdof) {
			int dof = dsa->locate(inode, 1 << jdof);
			if (dof >= 0) {
				int cdof = c_dsa->invRCN(dof);
				if (cdof >= 0) mergedF[inode][jdof] += subF[cdof]; // constrained
			}
		}
	}
}


template<>
void
GenSubDomain<double>::updatePrescribedDisp(GeomState *geomState, double deltaLambda) {
	if (numDirichlet > 0)
		geomState->updatePrescribedDisplacement(dbc, numDirichlet, deltaLambda);
}

template<class Scalar>
void
GenSubDomain<Scalar>::updatePrescribedDisp(GeomState *geomState) {
	if (domain->solInfo().initialTime == 0.0) {
		// note 1: "if both IDISP and IDISP6 are present in the input file, FEM selects IDISP6 to initialize the displacement field"
		if ((domain->numInitDisp() > 0) && (domain->numInitDisp6() == 0)) {
			if (this->numInitDisp() > 0) {
				geomState->updatePrescribedDisplacement(this->getInitDisp(), this->numInitDisp(), nodes);
			}
		}

		if (domain->numInitDisp6() > 0) {
			if (this->numInitDisp6() > 0)
				geomState->updatePrescribedDisplacement(this->getInitDisp6(), this->numInitDisp6(), nodes);
		}

		if (numDirichlet > 0)
			geomState->updatePrescribedDisplacement(dbc, numDirichlet, nodes);
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::setMpcSparseMatrix() {
	const auto &mpc_primal = this->mpc_primal;

	// should always be the second term added to Src
	if (numMPC_primal > 0) {
		MPCsparse = std::make_shared<GenMpcSparse<Scalar>>(numMPC_primal, mpc_primal, cc_dsa.get());
		this->Src->addSparseMatrix(MPCsparse);

		// to take into account wet interface / mpc interaction
		if (numWIdof) {
			int mpcOffset = numCRNdof;
			this->Kcw_mpc = std::make_unique<GenMpcSparse<Scalar>>(numMPC_primal, mpc_primal,
			                                                       cc_dsa.get(), dsa, wetInterfaceMap.data(), mpcOffset);
		}
	}
}

template<class Scalar>
void GenSubDomain<Scalar>::setUserDefBC(double *usrDefDisp, double *usrDefVel, double *usrDefAcc, bool nlflag) {
	auto &mpc = this->mpc;
	// modify boundary condition values for output purposes
	int i;
	for (i = 0; i < claw->numUserDisp; ++i) {
		int dof = dsa->locate(claw->userDisp[i].nnum, 1 << claw->userDisp[i].dofnum);
		if (dof >= 0) {
			bcx[dof] = usrDefDisp[locToGlUserDispMap[i]];
			if (bcx_scalar) bcx_scalar[dof] = bcx[dof];
			vcx[dof] = usrDefVel[locToGlUserDispMap[i]];
			acx[dof] = usrDefAcc[locToGlUserDispMap[i]];
		}
	}
	updateUsddInDbc(usrDefDisp, locToGlUserDispMap); // CHECK

	if (nlflag) return; // don't need to adjust rhs for nonlinear dynamics

	for (int i = 0; i < numMPC; ++i) {
		mpc[i]->rhs = mpc[i]->original_rhs;
		for (int j = 0; j < mpc[i]->nterms; ++j) {
			int mpc_node = mpc[i]->terms[j].nnum;
			int mpc_dof = mpc[i]->terms[j].dofnum;
			for (int k = 0; k < numDirichlet; ++k) {
				int dbc_node = dbc[k].nnum;
				int dbc_dof = dbc[k].dofnum;
				if ((dbc_node == mpc_node) && (dbc_dof == mpc_dof)) {
					mpc[i]->rhs -= mpc[i]->terms[j].coef * dbc[k].val;
				}
			}
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::deleteKcc() {
	this->Kcc.reset(nullptr);
}

template<class Scalar>
void
GenSubDomain<Scalar>::multQt(int glMPCnum, const Scalar *beta, Scalar *result) const {
	auto &mpc = this->mpc;
	int iMPC = globalToLocalMPC[glMPCnum];
	for (int i = 0; i < mpc[iMPC]->nterms; ++i) {
		int dof = c_dsa->locate(mpc[iMPC]->terms[i].nnum,
		                        (1 << mpc[iMPC]->terms[i].dofnum));
		if (dof < 0) continue;
		result[0] += mpc[iMPC]->terms[i].coef * beta[dof];
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::setMpcRhs(Scalar *interfvec, double _t, int flag) {
	auto &mpc = this->mpc;
	// set the rhs of inequality mpcs to the geometric gap and reset the rhs of the equality mpcs to the original rhs
	// if flag = 0 then interfvec = C*u, used in nonlinear analyses to update LMPCs and tied surfaces
	// if flag = 1 then interfvec = C*(u-u_n), used in nonlinear analyses to update piecewise linear contact surfaces
	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		if ((mpc[locMpcNb]->getSource() != mpc::ContactSurfaces && flag == 0) ||
		    (mpc[locMpcNb]->getSource() == mpc::ContactSurfaces && sinfo.piecewise_contact && flag == 1)) {
			mpc[locMpcNb]->rhs = mpc[locMpcNb]->original_rhs - interfvec[scomm->mapT(SComm::mpc, i)];
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::updateMpcRhs(Scalar *interfvec) {
	auto &mpc = this->mpc;
	bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool) * numMPC);
	for (int i = 0; i < numMPC; ++i) mpcFlag[i] = true;
	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		if (mpcFlag[locMpcNb] == true) {
			mpc[locMpcNb]->rhs = mpc[locMpcNb]->original_rhs + interfvec[scomm->mapT(SComm::mpc, i)];
			mpcFlag[locMpcNb] = false;
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::printMpcStatus() {
	auto &mpc = this->mpc;
	for (int i = 0; i < numMPC; ++i) {
		if (mpc[i]->type == 1) {
			std::cerr << (mpc[i]->active ? 'o' : 'x');
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::computeContactPressure(Scalar *globStress, Scalar *globWeight) {
	auto &mpc = this->mpc;
	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		if (mpc[locMpcNb]->type == 1) { // inequality constraint
			for (int j = 0; j < mpc[locMpcNb]->nterms; ++j) {
				int node = mpc[locMpcNb]->terms[j].nnum;
				int glNode = (domain->outFlag) ? domain->nodeTable[glNums[node]] - 1 : glNums[node];
				globStress[glNode] += std::abs(localLambda[scomm->mapT(SComm::mpc, i)]);
				globWeight[glNode] += 1.0;
			}
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::computeLocalContactPressure(Scalar *stress, Scalar *weight) {
	auto &mpc = this->mpc;
	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int locMpcNb = scomm->mpcNb(i);
		if (mpc[locMpcNb]->type == 1) { // inequality constraint
			for (int j = 0; j < mpc[locMpcNb]->nterms; ++j) {
				int node = mpc[locMpcNb]->terms[j].nnum;
				stress[node] += std::abs(localLambda[scomm->mapT(SComm::mpc, i)]);
				weight[node] += 1.0;
			}
		}
	}
}


template<class Scalar>
void
GenSubDomain<Scalar>::getConstraintMultipliers(std::map<std::pair<int, int>, double> &mu, std::vector<double> &lambda) {
	auto &mpc = this->mpc;
	bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool) * numMPC);
	for (int i = 0; i < numMPC; ++i) mpcFlag[i] = true;

	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int l = scomm->mpcNb(i);
		if (mpcFlag[l]) {
			if (mpc[l]->getSource() == mpc::ContactSurfaces) {
				mu[mpc[l]->id] = (localLambda) ? localLambda[scomm->mapT(SComm::mpc, i)] : 0;
			} else {
				lambda.push_back((localLambda) ? localLambda[scomm->mapT(SComm::mpc, i)] : 0);
			}
			mpcFlag[l] = false;
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::getLocalMultipliers(std::vector<double> &lambda) {
	lambda.resize(numMPC);

	for (int i = 0; i < scomm->lenT(SComm::mpc); ++i) {
		int l = scomm->mpcNb(i);
		lambda[l] = (localLambda) ? localLambda[scomm->mapT(SComm::mpc, i)] : 0;
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::initialize() {
	MPCsparse = 0;
	corotators = 0;
	mpcForces = 0;
	M = 0;
	Muc = 0;
	bcx_scalar = 0;
	a = 0;
	b = 0;
	P = 0;
	Q = 0;
	rebuildPade = true;  // pade
	numC_deriv = 0;
	C_deriv = 0;
	Cuc_deriv = 0;
	numK_deriv = 0;
	K_deriv = 0;
	Kuc_deriv = 0;
	num_K_arubber = 0;
	K_arubber_l = 0;
	K_arubber_m = 0;
	Kuc_arubber_l = 0;
	Kuc_arubber_m = 0;
#ifdef HB_COUPLED_PRECOND
	precNodeToNode = 0;
#endif
	mpcStatus = 0;
	mpcStatus1 = 0;
	mpcStatus2 = 0;
	l2g = 0;
}

template<class Scalar>
GenSubDomain<Scalar>::~GenSubDomain() {

	if (mpcForces) {
		delete[] mpcForces;
	}
	if (M) {
		delete M;
	}
	if (Muc) {
		delete Muc;
	}
	// don't delete boundaryDOFs, rbms
	if (bcx_scalar) {
		delete[] bcx_scalar;
		bcx_scalar = 0;
	}
	// pade
	if (a) {
		for (int i = 0; i < ia; ++i) delete a[i];
		delete[] a;
	}
	if (b) {
		for (int i = 0; i < ib; ++i) delete b[i];
		delete[] b;
	}
	if (P) delete P;
	if (Q) delete Q;
	if (C_deriv) {
		for (int i = 0; i < numC_deriv; ++i) delete C_deriv[i];
		delete[] C_deriv;
	}
	if (Cuc_deriv) {
		for (int i = 0; i < numC_deriv; ++i) delete Cuc_deriv[i];
		delete[] Cuc_deriv;
	}
	if (K_deriv) {
		for (int i = 0; i < numK_deriv; ++i) delete K_deriv[i];
		delete[] K_deriv;
	}
	if (Kuc_deriv) {
		for (int i = 0; i < numK_deriv; ++i) delete Kuc_deriv[i];
		delete[] Kuc_deriv;
	}
	if (K_arubber_l) {
		for (int i = 0; i < num_K_arubber; ++i) delete K_arubber_l[i];
		delete[] K_arubber_l;
	}
	if (Kuc_arubber_l) {
		for (int i = 0; i < num_K_arubber; ++i) delete Kuc_arubber_l[i];
		delete[] Kuc_arubber_l;
	}
	if (K_arubber_m) {
		for (int i = 0; i < num_K_arubber; ++i) delete K_arubber_m[i];
		delete[] K_arubber_m;
	}
	if (Kuc_arubber_m) {
		for (int i = 0; i < num_K_arubber; ++i) delete Kuc_arubber_m[i];
		delete[] Kuc_arubber_m;
	}
#ifdef HB_COUPLED_PRECOND
	if(isMixedSub && precNodeToNode) {
		delete precNodeToNode;
		precNodeToNode = 0;
	}
#endif
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

	if (l2g) delete[] l2g;
}

template<class Scalar>
void
GenSubDomain<Scalar>::deleteMPCs() {
	this->mpc.clear();
	if (localToGlobalMPC) {
		delete[] localToGlobalMPC;
		localToGlobalMPC = 0;
	}
	this->deleteG();
	numMPC = 0;
	scomm->deleteTypeSpecificList(SComm::mpc);

	scomm->deleteTypeSpecificList(SComm::all);
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
	if (mpcMaster) {
		delete[] mpcMaster;
		mpcMaster = 0;
	}
	if (mpcToDof) {
		delete mpcToDof;
		mpcToDof = 0;
	}
	if (localMpcToMpc) {
		delete localMpcToMpc;
		localMpcToMpc = 0;
	}

}

template<class Scalar>
void
GenSubDomain<Scalar>::extractMPCs(int glNumMPC, ResizeArray<LMPCons *> &lmpc) {
	auto &mpc = this->mpc;
	// PASS 1 count number of local mpcs
	numMPC = 0;
	for (int iMPC = 0; iMPC < glNumMPC; ++iMPC) {
		if (lmpc[iMPC]->isPrimalMPC()) continue;
		for (int i = 0; i < lmpc[iMPC]->nterms; ++i) {
			if (globalToLocal(lmpc[iMPC]->terms[i].nnum) > -1) {
				numMPC++;
				break;
			}
		}
	}
	if (numMPC == 0) return;

	mpc.clear();
	mpc.resize(numMPC); // use specific class for subdomain lmpcs
	localToGlobalMPC = new int[numMPC];

	// PASS 2: Get the mpc values
	numMPC = 0;
	int count = 0;
	for (int iMPC = 0; iMPC < glNumMPC; ++iMPC) {
		if (lmpc[iMPC]->isPrimalMPC()) continue;
		int used = 0;
		for (int i = 0; i < lmpc[iMPC]->nterms; ++i) {
			if (globalToLocal(lmpc[iMPC]->terms[i].nnum) > -1) {
				//if(c_dsa->locate((lmpc[iMPC]->terms)[i].nnum, (1 << (lmpc[iMPC]->terms)[i].dofnum)) < 0) continue; // XXXX
				if (mpc[numMPC] == 0) {
					Scalar rhs = lmpc[iMPC]->template getRhs<Scalar>();
					GenLMPCTerm<Scalar> term0 = lmpc[iMPC]->template getTerm<Scalar>(i);
					if (lmpc[iMPC]->isBoundaryMPC()) {
						if (lmpc[iMPC]->psub == subNumber)
							term0.coef = /*(solInfo().solvercntl->fetiInfo.c_normalize) ? 0.707106781 :*/ 1.0;
						else if (lmpc[iMPC]->nsub == subNumber)
							term0.coef = /*(solInfo().solvercntl->fetiInfo.c_normalize) ? -0.707106781 :*/ -1.0;
					}
					mpc[numMPC] = std::make_unique<SubLMPCons<Scalar>>(lmpc[iMPC]->lmpcnum, rhs, term0, lmpc[iMPC]->nterms, i);
					mpc[numMPC]->type = lmpc[iMPC]->type; // this is to be phased out
					mpc[numMPC]->setType(lmpc[iMPC]->getType());
					mpc[numMPC]->setSource(lmpc[iMPC]->getSource());
					mpc[numMPC]->id = lmpc[iMPC]->id;
					used = 1;
				} else {
					GenLMPCTerm<Scalar> term = lmpc[iMPC]->template getTerm<Scalar>(i);
					mpc[numMPC]->addterm(term, i);
				}
			}
		}
		if (used == 1) {
			localToGlobalMPC[numMPC] = count;
			numMPC++;
		}
		count++;
	}
	globalToLocalMPC.initialize(numMPC, localToGlobalMPC);
	//globalToLocalMPC.print();

#ifdef DEBUG_MPC
																															std::cerr << "DUAL MPCs: \n";
  for(int iMPC = 0; iMPC < numMPC; ++iMPC) mpc[iMPC]->print();
#endif
}

template<class Scalar>
void
GenSubDomain<Scalar>::printLMPC() {
	auto &mpc = this->mpc;
	std::cerr << "sub = " << subNumber << ", numMPC = " << numMPC << std::endl;
	for (int iMPC = 0; iMPC < numMPC; ++iMPC) {
		std::cerr << "local mpc # " << iMPC << ", global mpc # " << localToGlobalMPC[iMPC] << std::endl;
		mpc[iMPC]->print();
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::extractMPCs_primal(int glNumMPC, ResizeArray<LMPCons *> &lmpc) {
	auto &mpc_primal = this->mpc_primal;
	// PASS 1 count number of local mpcs
	numMPC_primal = 0;
	for (int iMPC = 0; iMPC < glNumMPC; ++iMPC) {
		if (!lmpc[iMPC]->isPrimalMPC()) continue;
		for (int i = 0; i < lmpc[iMPC]->nterms; ++i) {
			if (globalToLocal(lmpc[iMPC]->terms[i].nnum) > -1) {
				numMPC_primal++;
				break;
			}
		}
	}
	if (numMPC_primal == 0) return;

	mpc_primal.resize(numMPC_primal); // use specific class for subdomain lmpcs
	localToGlobalMPC_primal = new int[numMPC_primal];

	// PASS 2: Get the mpc values
	numMPC_primal = 0;
	int count = 0;
	for (int iMPC = 0; iMPC < glNumMPC; ++iMPC) {
		if (!lmpc[iMPC]->isPrimalMPC()) continue;
		int used = 0;
		for (int i = 0; i < lmpc[iMPC]->nterms; ++i) {
			if (globalToLocal(lmpc[iMPC]->terms[i].nnum) > -1) {
				if (mpc_primal[numMPC_primal] == 0) {
					Scalar rhs = lmpc[iMPC]->template getRhs<Scalar>();
					GenLMPCTerm<Scalar> term0 = lmpc[iMPC]->template getTerm<Scalar>(i);
					if (lmpc[iMPC]->isBoundaryMPC()) {
						if (lmpc[iMPC]->psub == subNumber)
							term0.coef = /*(solInfo().solvercntl->fetiInfo.c_normalize) ? 0.707106781 :*/ 1.0;
						else if (lmpc[iMPC]->nsub == subNumber)
							term0.coef = /*(solInfo().solvercntl->fetiInfo.c_normalize) ? -0.707106781 :*/ -1.0;
					}
					mpc_primal[numMPC_primal] = std::make_unique<SubLMPCons<Scalar>>(numMPC_primal, rhs, term0, lmpc[iMPC]->nterms,
					                                                   i);
					mpc_primal[numMPC_primal]->type = lmpc[iMPC]->type;
					mpc_primal[numMPC_primal]->setType(lmpc[iMPC]->getType());
					mpc_primal[numMPC_primal]->setSource(lmpc[iMPC]->getSource());
					mpc_primal[numMPC_primal]->id = lmpc[iMPC]->id;
					used = 1;
				} else {
					GenLMPCTerm<Scalar> term = lmpc[iMPC]->template getTerm<Scalar>(i);
					mpc_primal[numMPC_primal]->addterm(term, i);
				}
			}
		}
		if (used == 1) {
			localToGlobalMPC_primal[numMPC_primal] = count;
			numMPC_primal++;
		}
		count++;
	}
	globalToLocalMPC_primal.initialize(numMPC_primal, localToGlobalMPC_primal);
	//globalToLocalMPC_primal.print();

#ifdef DEBUG_MPC
																															std::cerr << "PRIMAL MPCs: \n";
  for(int iMPC = 0; iMPC < numMPC_primal; ++iMPC) mpc_primal[iMPC]->print();
#endif
}

template<class Scalar>
void
GenSubDomain<Scalar>::locateMpcDofs() {
	auto &mpc = this->mpc;
	auto &mpc_primal = this->mpc_primal;
	for (int i = 0; i < numMPC; ++i)
		for (int k = 0; k < mpc[i]->nterms; ++k) {
			(mpc[i]->terms)[k].dof = dsa->locate((mpc[i]->terms)[k].nnum, (1 << (mpc[i]->terms)[k].dofnum));
			(mpc[i]->terms)[k].cdof = c_dsa->locate((mpc[i]->terms)[k].nnum, (1 << (mpc[i]->terms)[k].dofnum));
			if (cc_dsa)
				(mpc[i]->terms)[k].ccdof = cc_dsa->locate((mpc[i]->terms)[k].nnum, (1 << (mpc[i]->terms)[k].dofnum));
		}
	for (int i = 0; i < numMPC_primal; ++i)
		for (int k = 0; k < mpc_primal[i]->nterms; ++k) {
			(mpc_primal[i]->terms)[k].dof = dsa->locate((mpc_primal[i]->terms)[k].nnum,
			                                            (1 << (mpc_primal[i]->terms)[k].dofnum));
			(mpc_primal[i]->terms)[k].cdof = c_dsa->locate((mpc_primal[i]->terms)[k].nnum,
			                                               (1 << (mpc_primal[i]->terms)[k].dofnum));
			if (cc_dsa)
				(mpc_primal[i]->terms)[k].ccdof = cc_dsa->locate((mpc_primal[i]->terms)[k].nnum,
				                                                 (1 << (mpc_primal[i]->terms)[k].dofnum));
		}
}



template<class Scalar>
void
GenSubDomain<Scalar>::constraintProduct(int num_vect, const double *R[], Scalar **V, int trans) {
	auto &mpc = this->mpc;
	int i, n, iMPC;

	if (trans) {
		for (iMPC = 0; iMPC < numMPC; ++iMPC) {
			for (n = 0; n < num_vect; ++n) {
				for (i = 0; i < mpc[iMPC]->nterms; ++i) {
					int dof = c_dsa->locate(mpc[iMPC]->terms[i].nnum,
					                        (1 << mpc[iMPC]->terms[i].dofnum));
					if (dof < 0) continue;
					V[n][dof] += mpc[iMPC]->terms[i].coef * R[n][dof];
				}
			}
		}
	} else {
		for (iMPC = 0; iMPC < numMPC; ++iMPC) {
			for (i = 0; i < mpc[iMPC]->nterms; ++i) {
				int dof = c_dsa->locate(mpc[iMPC]->terms[i].nnum,
				                        (1 << mpc[iMPC]->terms[i].dofnum));
				if (dof < 0) continue;
				for (n = 0; n < num_vect; ++n)
					V[iMPC][n] += mpc[iMPC]->terms[i].coef * R[iMPC][dof];
			}
		}
	}
}



template<>
void
GenSubDomain<double>::dualConstraintProjection(std::vector<std::map<int,double> > &W,
                                               Rom::DistrVecBasis &CtW,
                                               Eigen::Matrix<double,Eigen::Dynamic,1> &WtRhs,
                                               int startCol, int blockCols){

	for(int k = 0; k < blockCols; ++k){ // loop through each dual vector
		bool *mpcFlag =  (bool *) dbg_alloca(sizeof(bool)*numMPC);
		for(int i = 0; i < numMPC; ++i) mpcFlag[i] = true;
		// get subvector for target and dual vector  k
		std::map<int,double> &basisVec = W[startCol+k];
		StackVector target(CtW[k].subData(localSubNumber), CtW[k].subLen(localSubNumber)); target.zero();
		for(int l = 0; l < scomm->lenT(SComm::mpc); ++l) { // loop over all MPCs
			int i = scomm->mpcNb(l);
			if(!mpcFlag[i]) continue;
			if(this->mpc[i]->getSource() == mpc::ContactSurfaces) { // check if contact constraint
				// get slave node from mpc and extract corresponding row from dual basis
				std::map<int, double>::iterator it1 = basisVec.find(this->mpc[i]->id.second); // get correct row of dual basis
				// if, the node is in this vector, get multiplier slot associated with this slave node, othewise skip to next one
				if(it1 == basisVec.end()){ continue;}
				double rowVal = it1->second;
				//fprintf(stderr,"slave Node: %d, rowVal %3.2e, rhs: %3.2e\n",mpc[i]->id.second, rowVal, mpc[i]->rhs);
				int sNode = this->mpc[i]->id.second; //global node id
				const int pnId = globalToLocal(sNode);
				if(pnId > 0) {
					WtRhs(k) += this->mpc[i]->rhs*rowVal;
				}
				for(int j = 0; j < this->mpc[i]->nterms; ++j) { // number of nterms associated with this multiplier
					int dof = c_dsa->locate(this->mpc[i]->terms[j].nnum, (1 << this->mpc[i]->terms[j].dofnum));
					//fprintf(stderr,"mpc[%d]->terms[%d].nnum: %d, .dofnum: %d , dof: %d\n",i,j,mpc[i]->terms[j].nnum, mpc[i]->terms[j].dofnum,dof);
					if(dof < 0) continue;
					target[dof] += this->mpc[i]->terms[j].coef*rowVal;
				}
			}
			mpcFlag[i] = false;
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::addConstraintForces(std::map<std::pair<int, int>, double> &mu, std::vector<double> &lambda,
                                          GenVector<Scalar> &f) {
	auto &mpc = this->mpc;

	bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool) * numMPC);
	for (int i = 0; i < numMPC; ++i) mpcFlag[i] = true;

	std::vector<double>::iterator it2 = lambda.begin();
	for (int l = 0; l < scomm->lenT(SComm::mpc); ++l) { // loop over all MPCs
		int i = scomm->mpcNb(l);
		if (!mpcFlag[i]) continue;
		if (mpc[i]->getSource() == mpc::ContactSurfaces) { // check if contact constraint
			std::map<std::pair<int, int>, double>::iterator it1 = mu.find(mpc[i]->id); // get multiplier value
			if (it1 != mu.end()) {
				for (int j = 0; j < mpc[i]->nterms; ++j) { // number of nterms associated with this multiplier
					int dof = c_dsa->locate(mpc[i]->terms[j].nnum, (1 << mpc[i]->terms[j].dofnum));
					if (dof < 0) continue;
					f[dof] += mpc[i]->terms[j].coef * it1->second;
				}
			}
		} else {
			if (it2 != lambda.end()) {
				for (int j = 0; j < mpc[i]->nterms; ++j) {
					int dof = c_dsa->locate(mpc[i]->terms[j].nnum, (1 << mpc[i]->terms[j].dofnum));
					if (dof < 0) continue;
					f[dof] += mpc[i]->terms[j].coef * (*it2);
				}
				it2++;
			}
		}
		mpcFlag[i] = false;
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::addCConstraintForces(std::map<std::pair<int, int>, double> &mu, std::vector<double> &lambda,
                                           GenVector<Scalar> &fc, double s) {
	auto &mpc = this->mpc;

	bool *mpcFlag = (bool *) dbg_alloca(sizeof(bool) * numMPC);
	for (int i = 0; i < numMPC; ++i) mpcFlag[i] = true;

	std::vector<double>::iterator it2 = lambda.begin();
	for (int l = 0; l < scomm->lenT(SComm::mpc); ++l) {
		int i = scomm->mpcNb(l);
		if (!mpcFlag[i]) continue;
		if (mpc[i]->getSource() == mpc::ContactSurfaces) { // contact
			std::map<std::pair<int, int>, double>::iterator it1 = mu.find(mpc[i]->id);
			if (it1 != mu.end()) {
				for (int j = 0; j < mpc[i]->nterms; ++j) {
					int dof = dsa->locate(mpc[i]->terms[j].nnum, (1 << mpc[i]->terms[j].dofnum));
					if (dof < 0) continue;
					int cdof = c_dsa->invRCN(dof);
					if (cdof >= 0)
						fc[cdof] += s * mpc[i]->terms[j].coef * it1->second;
				}
			}
		} else {
			if (it2 != lambda.end()) {
				for (int j = 0; j < mpc[i]->nterms; ++j) {
					int dof = dsa->locate(mpc[i]->terms[j].nnum, (1 << mpc[i]->terms[j].dofnum));
					if (dof < 0) continue;
					int cdof = c_dsa->invRCN(dof);
					if (cdof >= 0)
						fc[cdof] += s * mpc[i]->terms[j].coef * (*it2);
				}
				it2++;
			}
		}
		mpcFlag[i] = false;
	}
}

template<class Scalar>
template<class Scalar1>
void
GenSubDomain<Scalar>::dispatchNodalData(FSCommPattern<Scalar> *pat, DistVec<Scalar1> *vec) {
	Scalar1 *v = vec->subData(localSubNumber);
	for (int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		FSSubRecInfo<Scalar> sInfo = pat->getSendBuffer(subNumber, scomm->subNums[iSub]);
		for (int iNode = 0; iNode < scomm->sharedNodes->num(iSub); ++iNode)
			sInfo.data[iNode] = v[(*scomm->sharedNodes)[iSub][iNode]];
	}
}

template<class Scalar>
template<class Scalar1>
void
GenSubDomain<Scalar>::addNodalData(FSCommPattern<Scalar> *pat, DistVec<Scalar1> *vec) {
	Scalar1 *v = vec->subData(localSubNumber);
	for (int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		FSSubRecInfo<Scalar> rInfo = pat->recData(scomm->subNums[iSub], subNumber);
		for (int iNode = 0; iNode < scomm->sharedNodes->num(iSub); ++iNode)
			v[(*scomm->sharedNodes)[iSub][iNode]] += rInfo.data[iNode];
	}
}

template<>
template<>
void
GenSubDomain<DComplex>::dispatchNodalData(FSCommPattern<DComplex> *pat, DistVec<double> *vec) {
	double *v = vec->subData(localSubNumber);
	for (int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		FSSubRecInfo<DComplex> sInfo = pat->getSendBuffer(subNumber, scomm->subNums[iSub]);
		for (int iNode = 0; iNode < scomm->sharedNodes->num(iSub); ++iNode)
			sInfo.data[iNode] = v[(*scomm->sharedNodes)[iSub][iNode]];
	}
}

template<>
template<>
void
GenSubDomain<DComplex>::addNodalData(FSCommPattern<DComplex> *pat, DistVec<double> *vec) {
	double *v = vec->subData(localSubNumber);
	for (int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		FSSubRecInfo<DComplex> rInfo = pat->recData(scomm->subNums[iSub], subNumber);
		for (int iNode = 0; iNode < scomm->sharedNodes->num(iSub); ++iNode)
			v[(*scomm->sharedNodes)[iSub][iNode]] += ScalarTypes::Real(rInfo.data[iNode]);
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::dispatchInterfaceGeomState(FSCommPattern<double> *pat, GeomState *geomState) {
	for (int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		FSSubRecInfo<double> sInfo = pat->getSendBuffer(subNumber, scomm->subNums[iSub]);
		for (int iNode = 0; iNode < scomm->sharedNodes->num(iSub); ++iNode) {
			if ((*scomm->sharedNodes)[iSub][iNode] < geomState->numNodesFixed()) {
				NodeState &ns = (*geomState)[(*scomm->sharedNodes)[iSub][iNode]];
				sInfo.data[16 * iNode] = 1; // status
				sInfo.data[16 * iNode + 1] = ns.x;
				sInfo.data[16 * iNode + 2] = ns.y;
				sInfo.data[16 * iNode + 3] = ns.z;
				sInfo.data[16 * iNode + 4] = ns.R[0][0];
				sInfo.data[16 * iNode + 5] = ns.R[0][1];
				sInfo.data[16 * iNode + 6] = ns.R[0][2];
				sInfo.data[16 * iNode + 7] = ns.R[1][0];
				sInfo.data[16 * iNode + 8] = ns.R[1][1];
				sInfo.data[16 * iNode + 9] = ns.R[1][2];
				sInfo.data[16 * iNode + 10] = ns.R[2][0];
				sInfo.data[16 * iNode + 11] = ns.R[2][1];
				sInfo.data[16 * iNode + 12] = ns.R[2][2];
				sInfo.data[16 * iNode + 13] = ns.theta[0];
				sInfo.data[16 * iNode + 14] = ns.theta[1];
				sInfo.data[16 * iNode + 15] = ns.theta[2];
			} else {
				for (int j = 0; j < 16; ++j) sInfo.data[16 * iNode + j] = 0;
			}
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::collectInterfaceGeomState(FSCommPattern<double> *pat, GeomState *geomState) {
	for (int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		FSSubRecInfo<double> rInfo = pat->recData(scomm->subNums[iSub], subNumber);
		for (int iNode = 0; iNode < scomm->sharedNodes->num(iSub); ++iNode) {
			if ((*scomm->sharedNodes)[iSub][iNode] >= geomState->numNodesFixed() && rInfo.data[16 * iNode] != 0) {
				NodeState &ns = (*geomState)[(*scomm->sharedNodes)[iSub][iNode]];
				ns.x = rInfo.data[16 * iNode + 1];
				ns.y = rInfo.data[16 * iNode + 2];
				ns.z = rInfo.data[16 * iNode + 3];
				ns.R[0][0] = rInfo.data[16 * iNode + 4];
				ns.R[0][1] = rInfo.data[16 * iNode + 5];
				ns.R[0][2] = rInfo.data[16 * iNode + 6];
				ns.R[1][0] = rInfo.data[16 * iNode + 7];
				ns.R[1][1] = rInfo.data[16 * iNode + 8];
				ns.R[1][2] = rInfo.data[16 * iNode + 9];
				ns.R[2][0] = rInfo.data[16 * iNode + 10];
				ns.R[2][1] = rInfo.data[16 * iNode + 11];
				ns.R[2][2] = rInfo.data[16 * iNode + 12];
				ns.theta[0] = rInfo.data[16 * iNode + 13];
				ns.theta[1] = rInfo.data[16 * iNode + 14];
				ns.theta[2] = rInfo.data[16 * iNode + 15];
			}
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::dispatchInterfaceGeomStateDynam(FSCommPattern<double> *pat, GeomState *geomState) {
	for (int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		FSSubRecInfo<double> sInfo = pat->getSendBuffer(subNumber, scomm->subNums[iSub]);
		for (int iNode = 0; iNode < scomm->sharedNodes->num(iSub); ++iNode) {
			if ((*scomm->sharedNodes)[iSub][iNode] < geomState->numNodesFixed()) {
				NodeState &ns = (*geomState)[(*scomm->sharedNodes)[iSub][iNode]];
				sInfo.data[28 * iNode] = 1; // status
				sInfo.data[28 * iNode + 1] = ns.x;
				sInfo.data[28 * iNode + 2] = ns.y;
				sInfo.data[28 * iNode + 3] = ns.z;
				sInfo.data[28 * iNode + 4] = ns.R[0][0];
				sInfo.data[28 * iNode + 5] = ns.R[0][1];
				sInfo.data[28 * iNode + 6] = ns.R[0][2];
				sInfo.data[28 * iNode + 7] = ns.R[1][0];
				sInfo.data[28 * iNode + 8] = ns.R[1][1];
				sInfo.data[28 * iNode + 9] = ns.R[1][2];
				sInfo.data[28 * iNode + 10] = ns.R[2][0];
				sInfo.data[28 * iNode + 11] = ns.R[2][1];
				sInfo.data[28 * iNode + 12] = ns.R[2][2];
				sInfo.data[28 * iNode + 13] = ns.theta[0];
				sInfo.data[28 * iNode + 14] = ns.theta[1];
				sInfo.data[28 * iNode + 15] = ns.theta[2];
				sInfo.data[28 * iNode + 16] = ns.v[0];
				sInfo.data[28 * iNode + 17] = ns.v[1];
				sInfo.data[28 * iNode + 18] = ns.v[2];
				sInfo.data[28 * iNode + 19] = ns.v[3];
				sInfo.data[28 * iNode + 20] = ns.v[4];
				sInfo.data[28 * iNode + 21] = ns.v[5];
				sInfo.data[28 * iNode + 22] = ns.a[0];
				sInfo.data[28 * iNode + 23] = ns.a[1];
				sInfo.data[28 * iNode + 24] = ns.a[2];
				sInfo.data[28 * iNode + 25] = ns.a[3];
				sInfo.data[28 * iNode + 26] = ns.a[4];
				sInfo.data[28 * iNode + 27] = ns.a[5];
			} else {
				for (int j = 0; j < 28; ++j) sInfo.data[28 * iNode + j] = 0;
			}
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::collectInterfaceGeomStateDynam(FSCommPattern<double> *pat, GeomState *geomState) {
	for (int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		FSSubRecInfo<double> rInfo = pat->recData(scomm->subNums[iSub], subNumber);
		for (int iNode = 0; iNode < scomm->sharedNodes->num(iSub); ++iNode) {
			if ((*scomm->sharedNodes)[iSub][iNode] >= geomState->numNodesFixed() && rInfo.data[28 * iNode] != 0) {
				NodeState &ns = (*geomState)[(*scomm->sharedNodes)[iSub][iNode]];
				ns.x = rInfo.data[28 * iNode + 1];
				ns.y = rInfo.data[28 * iNode + 2];
				ns.z = rInfo.data[28 * iNode + 3];
				ns.R[0][0] = rInfo.data[28 * iNode + 4];
				ns.R[0][1] = rInfo.data[28 * iNode + 5];
				ns.R[0][2] = rInfo.data[28 * iNode + 6];
				ns.R[1][0] = rInfo.data[28 * iNode + 7];
				ns.R[1][1] = rInfo.data[28 * iNode + 8];
				ns.R[1][2] = rInfo.data[28 * iNode + 9];
				ns.R[2][0] = rInfo.data[28 * iNode + 10];
				ns.R[2][1] = rInfo.data[28 * iNode + 11];
				ns.R[2][2] = rInfo.data[28 * iNode + 12];
				ns.theta[0] = rInfo.data[28 * iNode + 13];
				ns.theta[1] = rInfo.data[28 * iNode + 14];
				ns.theta[2] = rInfo.data[28 * iNode + 15];
				ns.v[0] = rInfo.data[28 * iNode + 16];
				ns.v[1] = rInfo.data[28 * iNode + 17];
				ns.v[2] = rInfo.data[28 * iNode + 18];
				ns.v[3] = rInfo.data[28 * iNode + 19];
				ns.v[4] = rInfo.data[28 * iNode + 20];
				ns.v[5] = rInfo.data[28 * iNode + 21];
				ns.a[0] = rInfo.data[28 * iNode + 22];
				ns.a[1] = rInfo.data[28 * iNode + 23];
				ns.a[2] = rInfo.data[28 * iNode + 24];
				ns.a[3] = rInfo.data[28 * iNode + 25];
				ns.a[4] = rInfo.data[28 * iNode + 26];
				ns.a[5] = rInfo.data[28 * iNode + 27];
			}
		}
	}
}

template<class Scalar>
void
GenSubDomain<Scalar>::dispatchInterfaceNodalInertiaTensors(FSCommPattern<double> *pat) {
#ifdef USE_EIGEN3
	for (int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		FSSubRecInfo<double> sInfo = pat->getSendBuffer(subNumber, scomm->subNums[iSub]);
		for (int iNode = 0; iNode < scomm->sharedNodes->num(iSub); ++iNode) {
			for (int j = 0; j < 3; ++j)
				for (int k = 0; k < 3; ++k)
					sInfo.data[9 * iNode + 3 * j + k] = Jn[(*scomm->sharedNodes)[iSub][iNode]](j, k);
		}
	}
#endif
}

template<class Scalar>
void
GenSubDomain<Scalar>::collectInterfaceNodalInertiaTensors(FSCommPattern<double> *pat) {
#ifdef USE_EIGEN3
	for (int iSub = 0; iSub < scomm->numNeighb; ++iSub) {
		FSSubRecInfo<double> rInfo = pat->recData(scomm->subNums[iSub], subNumber);
		for (int iNode = 0; iNode < scomm->sharedNodes->num(iSub); ++iNode) {
			for (int j = 0; j < 3; ++j)
				for (int k = 0; k < 3; ++k)
					Jn[(*scomm->sharedNodes)[iSub][iNode]](j, k) += rInfo.data[9 * iNode + 3 * j + k];
		}
	}
#endif
}

template<>
double GenSubDomain<double>::Bcx(int i);

template<>
DComplex GenSubDomain<DComplex>::Bcx(int i);

template<class Scalar>
void
GenSubDomain<Scalar>::multM(Scalar *localrhs, GenStackVector<Scalar> **u, int k) {
// RT: 11/17/08
//  double omega2 = (isFluid(0) && !sinfo.isCoupled) ? packedEset[0]->helmCoef() : geoSource->shiftVal();
	double omega2 = geoSource->shiftVal();
	double omega = sqrt(omega2);

	// assemble localvec for MatVec product
	Scalar *localvec = (Scalar *) dbg_alloca(sizeof(Scalar) * c_dsa->size());
	if (u == 0) {
		for (int i = 0; i < c_dsa->size(); ++i)
			localrhs[i] = 0.0;
		makeFreqSweepLoad(localrhs, k, omega);
		return;
	}
	for (int i = 0; i < c_dsa->size(); ++i) {
		localvec[i] = double(k) * (double(k - 1) * (*u[k - 1])[i] + 2.0 * omega * (*u[k])[i]);
	}
	M->mult(localvec, localrhs);


	if (C_deriv) {
		for (int j = 0; j <= k - 1; ++j) {
			if (C_deriv[k - j - 1]) {
				double ckj = DCombination(k, j);
				for (int i = 0; i < c_dsa->size(); ++i) localvec[i] = -ckj * (*u[j + 1])[i];
				C_deriv[k - j - 1]->multAdd(localvec, localrhs);
			}
		}
	}
	if (K_deriv) {
		for (int j = 0; j <= k - 1; ++j) {
			if (K_deriv[k - j]) {
				double ckj = DCombination(k, j);
				for (int i = 0; i < c_dsa->size(); ++i) localvec[i] = -ckj * (*u[j + 1])[i];
				K_deriv[k - j]->multAdd(localvec, localrhs);
			}
		}
	}

	makeFreqSweepLoad(localrhs, k, omega);  // this adds the residual effects of any prescribed displacements
	// and other loads that have non-zero kth derivative wrt omega
}

template<class Scalar>
void
GenSubDomain<Scalar>::makeFreqSweepLoad(Scalar *d, int iRHS, double omega) {
	auto &mpc = this->mpc;

	int numUncon = c_dsa->size();
	GenStackVector<Scalar> force(numUncon, d);

	// Compute Right Hand Side Force = Fext + Fgravity + Fnh + Fpressure
	buildFreqSweepRHSForce<Scalar>(force, Muc, Cuc_deriv, Kuc_deriv, iRHS, omega);

	for (int i = 0; i < numMPC; ++i) mpc[i]->rhs = 0.0; // set mpc rhs to zero for higher-order derivative solves
}


template<class Scalar>
void
GenSubDomain<Scalar>::multMCoupled1(Scalar *localrhs, GenStackVector<Scalar> **u, int k,
                                    FSCommPattern<Scalar> *wiPat) {
	auto &localw = this->localw;
	double omega2 = geoSource->shiftVal();
	double omega = sqrt(omega2);

	// assemble localvec for MatVec product
	Scalar *localvec = (Scalar *) dbg_alloca(sizeof(Scalar) * c_dsa->size());
	for (int i = 0; i < c_dsa->size(); ++i) {
		localvec[i] = double(k) * (double(k - 1) * (*u[k - 1])[i] + 2.0 * omega * (*u[k])[i]);
	}

	if (numWIdof) {
		int i, j;
		for (i = 0; i < numWIdof; ++i) { localw[i] = this->localw_copy[i] = 0; }

		for (i = 0; i < numWInodes; ++i) {
			//DofSet thisDofSet = wetInterfaceDofs[i]^DofSet(DofSet::Helm);
			DofSet thisDofSet = wetInterfaceDofs[i] ^(wetInterfaceDofs[i] & DofSet(DofSet::Helm)); // unmark Helm
			int nd = thisDofSet.count();
			int cdofs[6], dofs[6];
			c_dsa->number(wetInterfaceNodes[i], thisDofSet, cdofs);
			dsa->number(wetInterfaceNodes[i], thisDofSet, dofs);
			for (j = 0; j < nd; ++j) {
				localw[wetInterfaceMap[dofs[j]]] = localvec[cdofs[j]];
			}
		}

		// compute C uw to send to neighbors
		for (i = 0; i < scomm->numT(SComm::fsi); ++i) {
			if (subNumber != scomm->neighbT(SComm::fsi, i)) {
				FSSubRecInfo<Scalar> sInfo = wiPat->getSendBuffer(subNumber, scomm->neighbT(SComm::fsi, i));
				for (j = 0; j < numNeighbWIdof[i]; ++j) sInfo.data[j] = 0.0;
				this->neighbKww->multAdd(localw.data(), sInfo.data.data(), glToLocalWImap, neighbGlToLocalWImap[i], true);
			} else {
				this->neighbKww->multAdd(localw.data(), this->localw_copy.data(), glToLocalWImap, true);
			}
		}
	}

	M->mult(localvec, localrhs);

	if (C_deriv) {
		for (int j = 0; j <= k - 1; ++j) {
			if (C_deriv[k - j - 1]) {
				double ckj = DCombination(k, j);
				for (int i = 0; i < c_dsa->size(); ++i) localvec[i] = -ckj * (*u[j + 1])[i];
				C_deriv[k - j - 1]->multAdd(localvec, localrhs);
			}
		}
	}

	makeFreqSweepLoad(localrhs, k, omega);
}

template<class Scalar>
void
GenSubDomain<Scalar>::multMCoupled2(Scalar *localrhs, FSCommPattern<Scalar> *wiPat) {
	int i, j;
	if (numWIdof) {
		for (i = 0; i < scomm->numT(SComm::fsi); ++i) {
			if (subNumber != scomm->neighbT(SComm::fsi, i)) {
				FSSubRecInfo<Scalar> rInfo = wiPat->recData(scomm->neighbT(SComm::fsi, i), subNumber);
				for (j = 0; j < numWIdof; ++j) this->localw_copy[j] += rInfo.data[j] / this->wweight[j];
			}
		}

		// put localw_copy into localrhs
		double omega2 = geoSource->shiftVal();
		for (i = 0; i < numWInodes; ++i) {
			//DofSet thisDofSet = wetInterfaceDofs[i]&DofSet(DofSet::Helm);
			//if(thisDofSet.count()) {
			if (wetInterfaceDofs[i].contains(DofSet::Helm)) {
				int thisNode = wetInterfaceNodes[i];
				int cdof = c_dsa->locate(thisNode, DofSet::Helm);
				int dof = dsa->locate(thisNode, DofSet::Helm);
				localrhs[cdof] -= this->localw_copy[wetInterfaceMap[dof]] / omega2;
			}
		}
	}
}


template<class Scalar>
void
GenSubDomain<Scalar>::multWCAWE(Scalar *localrhs, GenStackVector<Scalar> **u, Scalar *pU, Scalar *pb, int maxRHS,
                                int iRHS) {
	double omega2 = geoSource->shiftVal();
	double omega = sqrt(omega2);

	Scalar *localvec = (Scalar *) dbg_alloca(sizeof(Scalar) * c_dsa->size());

//  pu = -Z{2}*W(:,ii-1); Z{2} = (-2.0*freqn * M+i*C)/1.0;
	for (int i = 0; i < c_dsa->size(); ++i) {
		localvec[i] = 2.0 * omega * (*u[iRHS - 1])[i];
	}
	M->mult(localvec, localrhs);

	if (C_deriv) {
		if (C_deriv[0]) {
			for (int i = 0; i < c_dsa->size(); ++i) localvec[i] = -(*u[iRHS - 1])[i];
			C_deriv[0]->multAdd(localvec, localrhs);
		}
	}

//  pu = pu + b{j+1}*PP(1,ii-j);
	double factorial = 1.0;
	for (int j = 1; j <= iRHS; j++) {
		factorial *= double(j);
		for (int i = 0; i < c_dsa->size(); ++i) localvec[i] = 0.0;
		makeFreqSweepLoad(localvec, j, omega);
		for (int i = 0; i < c_dsa->size(); ++i)
			localrhs[i] += pb[j - 1] * localvec[i] / factorial;
	}

//  pu = pu - Z{j+1}*W(:,1:(ii-j))*PP(:,ii-j);
//  assume j > 2 all 0 (not so for rubber), and for j=2: Z{3} = -M;
	if (iRHS > 1)
		for (int j = 2; j <= 2; j++) {
			for (int i = 0; i < c_dsa->size(); ++i)
				localvec[i] = 0.0;
			for (int k = 0; k < iRHS + 1 - j; k++) {
				for (int i = 0; i < c_dsa->size(); ++i)
					localvec[i] += pU[k + (j - 1) * (iRHS + 1)] * (*u[k])[i];
			}
			M->multAdd(localvec, localrhs);
		}

}

template<class Scalar>
void
GenSubDomain<Scalar>::makeBcx_scalar() {
	int numdofs = dsa->size();
	bcx_scalar = new Scalar[numdofs];
	for (int i = 0; i < numdofs; ++i) bcx_scalar[i] = Bcx(i);
}

template<class Scalar>
void
GenSubDomain<Scalar>::pade(GenStackVector<Scalar> *sol, GenStackVector<Scalar> **u, double *h, double x) {
	int len = sol->size();
	int i, j, k, n, r;
	// general N-point pade extrapolation
	if (rebuildPade) {
		// first time, allocate storage
		int l = domain->solInfo().getSweepParams()->padeL;
		int m = domain->solInfo().getSweepParams()->padeM;
		int padeN = domain->solInfo().getSweepParams()->padeN;
		int nRHS = domain->solInfo().getSweepParams()->nFreqSweepRHS;
		if (subNumber == 0)
			fprintf(stderr, " ... Computing %d-point Pade coefficients (l = %d, m = %d) ... \n", padeN, l, m);
		ia = l + 1;
		ib = m + 1;
		// allocate storage first time only
		ia = l + 1;
		if (a == 0) {
			a = new GenVector<Scalar> *[ia];
			for (i = 0; i < ia; ++i) a[i] = new GenVector<Scalar>(len);
			b = new GenVector<Scalar> *[ib];
			for (i = 0; i < ib; ++i) b[i] = new GenVector<Scalar>(len);
			P = new GenVector<Scalar>(len);
			Q = new GenVector<Scalar>(len);
		}
		// compute P, Q coefficients by solving system Ax = b
		int dim = l + m + 1;
		GenFullM<Scalar> A(dim); // LHS
		Scalar *v = (Scalar *) dbg_alloca(sizeof(Scalar) * dim);
		for (i = 0; i < sol->size(); ++i) {
			A.zero();
			// assemble
			for (j = 0; j < padeN; ++j) {
				int offset = j * (nRHS + 1) + 1;
				int I;
				for (n = 0; n < nRHS; ++n) {
					if ((I = j * nRHS + n) >= dim) break;
					// fill left block (a coefficients) of RHS
					for (k = n; k <= l; ++k)
						A[I][k] = -double(DFactorial(k) / DFactorial(k - n) * pow(h[j], k - n));
					// fill right block (b coefficients) of RHS
					for (r = 0; r <= n; ++r)
						for (k = r; k <= m; ++k)
							if (k > 0)
								A[I][l + k] += double(
									DCombination(n, r) * DFactorial(k) / DFactorial(k - r) * pow(h[j], k - r)) *
								               (*u[offset + n - r])[i];
					// fill LHS
					v[I] = -(*u[offset + n])[i];
				}
			}
			// factorize
			A.factor();
			// solve
			A.reSolve(v);
			// extract coefficients
			for (j = 0; j <= l; ++j) (*a[j])[i] = v[j];
			(*b[0])[i] = 1.0; // b0 = 1
			for (j = 1; j <= m; ++j) (*b[j])[i] = v[l + j];
		}
		rebuildPade = false;
	}
	// assemble P, Q and solution P/Q
	P->zero();
	Q->zero();
	for (i = 0; i < ia; ++i) P->linAdd(pow(x, i), *a[i]);
	for (i = 0; i < ib; ++i) Q->linAdd(pow(x, i), *b[i]);
	for (j = 0; j < sol->size(); ++j) (*sol)[j] = (*P)[j] / (*Q)[j];
}

template<class Scalar>
void
GenSubDomain<Scalar>::bmpcQualify(std::vector<LMPCons *> *bmpcs, int *pstatus, int *nstatus) {
	for (int i = 0; i < bmpcs->size(); ++i) {
		LMPCons *bmpc = (*bmpcs)[i];
		if (bmpc->psub == subNumber) {
			int ccdof = getCCDSA()->locate(globalToLocal((bmpc->terms)[0].nnum), (1 << (bmpc->terms)[0].dofnum));
			pstatus[i] = (ccdof > -1) ? 1 : 0;
		}
		if (bmpc->nsub == subNumber) {
			int ccdof = getCCDSA()->locate(globalToLocal((bmpc->terms)[0].nnum), (1 << (bmpc->terms)[0].dofnum));
			nstatus[i] = (ccdof > -1) ? 1 : 0;
		}
	}
}

template<class Scalar>
void GenSubDomain<Scalar>::mergeElemProps(double *props,
                                          double *weights,
                                          int propType) {
#ifdef DISTRIBUTED
	int *nodeNumbers = (int *) dbg_alloca(sizeof(int) * maxNumNodes);
#endif
	for (int iele = 0; iele < numele; ++iele) {
		const StructProp *sProp = packedEset[iele]->getProperty();
		double eleattr = 0.0;
		switch (propType) {
			case YOUNG:
				eleattr = sProp->E;
				break;
			case MDENS:
				eleattr = sProp->rho;
				break;
			case THICK:
				eleattr = sProp->eh;
				break;
			default:
				assert(0);
		}

		int NodesPerElement = elemToNode->num(iele);
#ifdef DISTRIBUTED
		packedEset[iele]->nodes(nodeNumbers);
#endif
		for (int k = 0; k < NodesPerElement; ++k) {
#ifdef DISTRIBUTED
			int glbNodeNum = nodeNumbers[k];
#else
																																	int locNodeNum = (*elemToNode)[iele][k];
	  int glbNodeNum = this->localToGlobal(locNodeNum);
#endif
			// not thread safe!!!
			props[glbNodeNum] += eleattr;
			weights[glbNodeNum] += 1.0;
		}
	}
	return;
}

template<class Scalar>
void
GenSubDomain<Scalar>::addSommer(SommerElement *ele) {
	ele->dom = this;
	sommer[numSommer++] = ele;
	//if(sinfo.ATDARBFlag != -2.0) packedEset.elemadd(numele++,ele);  // XDEBUG
	if (sinfo.ATDARBFlag != -2.0) {
		ele->renum(glToLocalNode);
		packedEset.elemadd(numele++, ele);
	}
}

#include "LOpsImpl.h"

template
class GenSubDomain<double>;

template
class GenSubDomain<std::complex<double>>;

template
void
GenSubDomain<std::complex<double>>::addNodalData<std::complex<double>>(FSCommPattern<std::complex<double>> *,
                                                                       NewVec::DistVec<std::complex<double>> *);

template
void
GenSubDomain<std::complex<double>>::dispatchNodalData<std::complex<double>>(FSCommPattern<std::complex<double>> *pat,
                                                                            DistVec<std::complex<double>> *vec);

template
void
GenSubDomain<std::complex<double>>::dispatchNodalData<double>(FSCommPattern<std::complex<double>> *pat,
                                                              DistVec<double> *vec);

template
void
GenSubDomain<double>::dispatchNodalData<double>(FSCommPattern<double> *pat, DistVec<double> *vec);

template
void
GenSubDomain<double>::addNodalData<double>(FSCommPattern<double> *, NewVec::DistVec<double> *);
