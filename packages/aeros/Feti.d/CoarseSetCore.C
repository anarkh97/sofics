#include <complex>
#include <Feti.d/CoarseSet.h>
#include <Math.d/matrix.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/BLAS.h>

template<>
void
GenCoarseSet<std::complex<double>>::addMPCContrib(int iMPC, GenSparseMatrix<std::complex<double>> *coarseMat,
                EqNumberer* eqNumber, int cOffset, int gOffset, int mpcOffset,
                double *locR, GenSolver<std::complex<double>> *s)
{
  int numMPCs = sd->numMPCs();
  if(numMPCs > 0) 
    fprintf(stderr, "WARNING: GenCoarseSet<Scalar>::addMPCContrib(...) not implemented \n");
}

template<>
void
GenCoarseSet<double>::addMPCContrib(int iMPC, GenSparseMatrix<double> *coarseMat,
                                    EqNumberer *eqNumber, int cOffset, int gOffset, int mpcOffset,
                                    double *locR, GenSolver<double> *s)
{
	int numMPCs = sd->numMPCs();

	// Find the maximum matrix allocation size
	int maxRight = (numMPCs > numGs) ? numMPCs : numGs;
	maxRight     = (numBCs > maxRight) ? numBCs : maxRight;

	// Find the maximum matrix allocation size
	int maxLeft  =  maxRight;
	int jSub;
	for(jSub =0; jSub < numNeighb; ++jSub) {
		int nGleft = neighbNumRBMs[jSub];
		if(nGleft > maxLeft) maxLeft = nGleft;
		int nCleft = (neighbCSubIndex) ?
		             neighbCSubIndex[jSub+1] - neighbCSubIndex[jSub] : 0;
		if(nCleft > maxLeft) maxLeft = nCleft;
	}
	double stackMem[maxLeft*numMPCs];

	// This sub-domain does not have a contribution
	if(numMPCs == 0) return;

	int myFGcol  = eqNumber->firstdof(myNum);
	int myCCol   = eqNumber->firstdof(myNum+cOffset);
	int myGCol   = eqNumber->firstdof(myNum+gOffset);

	// KHP: note the mpc numbering has to be done carefully.
	//      this number should be done per MPC as a sub-domain may have
	//      contributions to only a few non-consecutive MPCs
	int myMPCRow = eqNumber->firstdof(mpcOffset);
	// cerr << "myFGcol " << myFGcol << "  myGCol  " << myGCol << endl;

	// add QtKQ
	sd->getQtKQ(iMPC,stackMem);
	GenStackFullM<double> qtKq(numMPCs, 1, stackMem);
#ifdef DEBUG_MPC
	qtKq.print("--- QtKQ ---","QtKQ");
#endif
	int jMPC;
	GenFullM<double> qtkq(1,1);
	for(jMPC=0; jMPC<numMPCs; ++jMPC) {
		qtkq[0][0] = qtKq[jMPC][0];
#ifdef DEBUG_MPC
		fprintf(stderr,"QtKQ: Where am I adding this %d %d\n",
           eqNumber->firstdof(mpcOffset-iMPC)+sd->localToGlobalMPC[jMPC], myMPCRow);
#endif
		coarseMat->add(qtkq,
		               eqNumber->firstdof(mpcOffset-iMPC)+sd->getGlobalMPCIndex(jMPC),
		               myMPCRow);
	}

	// add QtKBtG
	if(numGs > 0) {
		GenStackFullM<double> GtBKQ(numGs,1,stackMem);
		GtBKQ.zero(); // NOTE: Never remove this!
		getGtMult(sd->getQtKpBt()+sd->getLocalMPCIndex(iMPC)*sd->interfLen(), GtBKQ.data());
#ifdef DEBUG_MPC
		GtBKQ.print("--- GtBKQ ---","GtBKQ");
   fprintf(stderr,"GtBKQ: Where am I adding this %d %d\n",myFGcol, myMPCRow);
#endif
		coarseMat->add(GtBKQ, myFGcol, myMPCRow);
	}

	// Now loop on the neighbors
	for(jSub =0; jSub < numNeighb; ++jSub) {
		int leftFGcol = eqNumber->firstdof(neighbs[jSub]);
		int nGleft    = neighbNumRBMs[jSub];
		GenStackFullM<double> GtBKQ(nGleft, 1, stackMem);
		if(neighbNumRBMs[jSub] > 0 ) {
			Tgemm('T', 'N', 1, nGleft, subSize(jSub), -1.0,
			      sd->getQtKpBt()+subOffset[jSub]+sd->getLocalMPCIndex(iMPC)*sd->interfLen()
				,gSize, neighbGs[jSub], leadingDimGs[jSub],
				  0.0, GtBKQ.data(), 1);
#ifdef DEBUG_MPC
			GtBKQ.print("---- Neighbors Gt (BKQ) ----","Gt^jBKQ");
#endif
			coarseMat->add(GtBKQ, leftFGcol, myMPCRow);
		}
	}

	// add QtKBtC
	if(numBCs > 0) {
		GenStackFullM<double> CtBKQ(numBCs, 1, stackMem);
		getCtMult(sd->getQtKpBt()+sd->getLocalMPCIndex(iMPC)*sd->interfLen(), CtBKQ.data());
		int i;
		for(i=0; i<numBCs; ++i)
			CtBKQ[i][0] = -CtBKQ[i][0];
		coarseMat->add(CtBKQ, myCCol, myMPCRow);
#ifdef DEBUG_MPC
		CtBKQ.print("---- CtBKQ ---","CtBKQ");
#endif
	}

	// Now loop on the neighbors
	double cornerSign = 1.0;
	for(jSub =0; jSub < numNeighb; ++jSub) {
		int nCleft = (neighbCSubIndex) ?
		             neighbCSubIndex[jSub+1] - neighbCSubIndex[jSub] : 0;
		if(nCleft == 0) continue;

		int (*leftCs)[3] = neighbCs + neighbCSubIndex[jSub];
		if(nCleft > 0) {
			int leftFCcol = eqNumber->firstdof(neighbs[jSub]+cOffset) +
			            neighbCs[neighbCSubIndex[jSub]][1];
			GenStackFullM<double> CtBKQ(nCleft, 1, stackMem);
			const double *BKQ = sd->getQtKpBt();
			int i;
			for(i = 0; i < nCleft; ++i)
				CtBKQ[i][0] = -cornerSign*BKQ[sd->getLocalMPCIndex(iMPC)*gSize + leftCs[i][0] ];
			coarseMat->add(CtBKQ, leftFCcol, myMPCRow);
		}
	}

	// add QtR
	if((numGs > 0) && (isDynamic==0)) {
		GenStackFullM<double> QtR(1, numGs, stackMem);
		QtR.zero();
		sd->multQt(iMPC, locR, numGs, QtR.data());
		GenFullM<double> RtQ = QtR.transpose();
		coarseMat->add(RtQ, myGCol, myMPCRow);
	}
}
