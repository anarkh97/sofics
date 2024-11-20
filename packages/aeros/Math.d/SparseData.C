#include <Utils.d/dbg_alloca.h>
#include <Math.d/SparseMatrix.h>
#include <Utils.d/Connectivity.h>
#include <Driver.d/Mpc.h>

SparseData::~SparseData()
{
}

void
SparseData::clean_up()
{
	unconstrNum.clear();
	unconstrNum.shrink_to_fit();
	constrndNum.clear();
	constrndNum.shrink_to_fit();
	xunonz.clear();
	xunonz.shrink_to_fit();
	rowu.clear();
	rowu.shrink_to_fit();
	colu.clear();
	colu.shrink_to_fit();
}

SparseData::SparseData(const Connectivity *con, const DofSetArray *dsa, const int *bc)
{
	neq = dsa->size(); // Total number of equations
	int numNodes = std::min( dsa->numNodes(), con->csize() ); // number of nodes

	// We build a temporary dof to Node table.
	unconstrNum.resize(neq);
	constrndNum.resize(neq);

	int cn = 0, un = 0;
	for (int i=0; i<neq; ++i) {
		if(bc[i] == BCFIXED) {
			unconstrNum[i] = -1;
			constrndNum[i] = cn;
			cn = cn+1;
		}
		else {
			unconstrNum[i] = un;
			constrndNum[i] = -1;
			un = un+1;
		}
	}

	numConstrained = cn;

	int *numFixed    = (int * ) dbg_alloca(sizeof(int)*numNodes);
	int *numNotFixed = (int * ) dbg_alloca(sizeof(int)*numNodes);

	for (int i=0; i<numNodes; ++i) {
		numFixed[i] = 0;
		int myFirstDof = dsa->firstdof(i);
		int myNumDofs  = dsa->weight(i);
		for (int j=0; j<myNumDofs; ++j)
			if(constrndNum[myFirstDof+j] >= 0) numFixed[i] += 1;
		numNotFixed[i] = myNumDofs - numFixed[i];
	}

	xunonz.resize(numConstrained+1);

	for (int i=0; i<numNodes; ++i) {
		if(numFixed[i] == 0) continue;
		int nentries=0;
		for (int j=0; j < con->num(i); ++j)
			nentries += numNotFixed[(*con)[i][j]];
		int myFirstDof = dsa->firstdof(i);
		int myNumDofs  = dsa->weight(i);
		for (int j=0; j<myNumDofs; ++j)
			if(constrndNum[myFirstDof+j] >= 0)
				xunonz[constrndNum[myFirstDof+j]] = nentries;
	}

	int count = 0;
	for (int i=0; i<numConstrained; ++i) {
		int temp = xunonz[i];
		xunonz[i] = count;
		count += temp;
	}
	// count = # of entries in kuc

	xunonz[numConstrained] = count;

	rowu.resize(count);

	for (int i=0; i<numNodes; ++i) {
		if(numFixed[i] == 0) continue;
		int iFirstDof = dsa->firstdof(i);
		int nentries = 0;
		for (int j =0; j < con->num(i); ++j) {
			int jNode = (*con)[i][j];
			if(numNotFixed[jNode] == 0) continue;
			int jFirstDof = dsa->firstdof(jNode);
			for (int k = 0; k < dsa->weight(jNode); ++k) {
				if(unconstrNum[jFirstDof + k] >= 0) {
					for (int l = 0; l < dsa->weight(i); ++l)
						if(constrndNum[iFirstDof+l] >= 0)
							rowu[ xunonz [constrndNum[iFirstDof+l]] +nentries] = unconstrNum[jFirstDof + k];
					nentries++;
				}
			}
		}
	}

}

SparseData::SparseData(const Connectivity *con, const DofSetArray *dsa,
                       const DofSetArray *c_dsa)
{
	neq = dsa->size();

	int numNodes = c_dsa->numNodes();
	if(con->csize() < numNodes) numNodes = con->csize();

	numUncon       = c_dsa->size();  // number of unconstrained dof
	unconstrNum = c_dsa->getUnconstrNum();
	if(numUncon == 0) return;
	numConstrained = neq - numUncon; // number of constrained dof
	constrndNum = c_dsa->getConstrndNum();
	if(numConstrained == 0) return;

	std::vector<int> numFixed(numNodes);
	std::vector<int> numNotFixed(numNodes);

	for (int i=0; i<numNodes; ++i) {
		numFixed[i] = 0;
		int myFirstDof = dsa->firstdof(i);
		int myNumDofs  = dsa->weight(i);
		for (int j = 0; j < myNumDofs; ++j)
			if(constrndNum[myFirstDof+j] >= 0) numFixed[i] += 1;
		numNotFixed[i] = myNumDofs - numFixed[i];
	}

	xunonz.resize(numConstrained+1);

	for (int i=0; i<numNodes; ++i) {
		if(numFixed[i] == 0) continue;
		int nentries=0;
		for (int j=0; j<con->num(i); ++j)
			nentries += numNotFixed[(*con)[i][j]];
		int myFirstDof = dsa->firstdof(i);
		int myNumDofs  = dsa->weight(i);
		for (int j=0; j<myNumDofs; ++j)
			if(constrndNum[myFirstDof+j] >= 0)
				xunonz[constrndNum[myFirstDof+j]] = nentries;
	}

	int count = 0;
	for (int i=0; i<numConstrained; ++i) {
		int temp = xunonz[i];
		xunonz[i] = count;
		count += temp;
	}
	xunonz[numConstrained] = count;

	rowu.resize(count);

	for (int i=0; i<numNodes; ++i) {
		if(numFixed[i] == 0) continue;
		int iFirstDof = dsa->firstdof(i);
		int nentries=0;
		for (int j=0; j<con->num(i); ++j) {
			int jNode = (*con)[i][j];
			if(numNotFixed[jNode] == 0) continue;
			int jFirstDof = dsa->firstdof(jNode);
			for (int k=0; k<dsa->weight(jNode); ++k) {
				if(unconstrNum[jFirstDof + k] >= 0) {
					for (int l = 0; l < dsa->weight(i); ++l)
						if(constrndNum[iFirstDof+l] >= 0)
							rowu[ xunonz [constrndNum[iFirstDof+l]] +nentries] = unconstrNum[jFirstDof + k];
					nentries++;
				}
			}
		}
	}

}

SparseData::SparseData(const Connectivity *con, const DofSetArray *dsa, const int *glBoundMap, const int *glInternalMap)
{

	int i,j,k,l,iDof;

	// Get the global length
	neq = dsa->size();
	constrndNum.assign(glBoundMap, glBoundMap + neq);
	unconstrNum .assign(glInternalMap, glInternalMap + neq);
	// Get the number of nodes
	int numNodes = dsa->numNodes();
	if(con->csize() < numNodes) numNodes = con->csize();

	// Compute number of Boundary dofs and number of Internal dofs
	int boundaryLen = 0;
	int internalLen = 0;
	for(iDof = 0; iDof < neq; ++iDof) {
		if(glBoundMap[iDof] >= 0 ) boundaryLen += 1;
		if(glInternalMap[iDof] >= 0) internalLen += 1;
	}

	xunonz.resize(boundaryLen+1);
	numConstrained = boundaryLen;
	numUncon       = internalLen; //HB
	std::vector<int> numBound   (numNodes);
	std::vector<int> numInternal(numNodes);

	for(i=0; i<numNodes; ++i) {
		numBound[i]   = 0;
		numInternal[i] = 0;
		int myFirstDof = dsa->firstdof(i);
		int myNumDofs  = dsa->weight(i);
		for(j = 0; j < myNumDofs; ++j){
			if(glBoundMap[myFirstDof+j] >= 0) numBound[i] += 1;
			if(glInternalMap[myFirstDof+j] >= 0) numInternal[i] += 1;
		}
	}

	for(i=0; i<numNodes; ++i) {
		if(numBound[i] == 0) continue;
		int nentries=0;
		for(j=0; j<con->num(i); ++j)
			nentries += numInternal[(*con)[i][j]];
		int myFirstDof = dsa->firstdof(i);
		int myNumDofs  = dsa->weight(i);
		for(j=0; j<myNumDofs; ++j)
			if(glBoundMap[myFirstDof+j] >= 0)
				xunonz[glBoundMap[myFirstDof+j]] = nentries;
	}

	int count = 0;
	for(i=0; i<boundaryLen; ++i) {
		int temp = xunonz[i];
		xunonz[i] = count;
		count += temp;
	}
	xunonz[boundaryLen] = count;

	rowu.resize(count);
	for(i=0; i<numNodes; ++i) {
		if(numBound[i] == 0) continue;
		int iFirstDof = dsa->firstdof(i);
		int nentries=0;
		for(j=0; j<con->num(i); ++j) {
			int jNode = (*con)[i][j];
			if(numInternal[jNode] == 0) continue;
			int jFirstDof = dsa->firstdof(jNode);
			for(k=0; k<dsa->weight(jNode); ++k) {
				if(glInternalMap[jFirstDof + k] >= 0) {
					for(l = 0; l < dsa->weight(i); ++l)
						if(glBoundMap[iFirstDof+l] >= 0)
							rowu[ xunonz [glBoundMap[iFirstDof+l]] +nentries]
								= glInternalMap[jFirstDof + k];
					nentries++;
				}
			}
		}
	}
}

// used for Mumps & Spooles (1st constructor)
// Constructor for DBSparseMatrix data structures
SparseData::SparseData(const EqNumberer *_dsa, const Connectivity *cn, const int *rCN, int expand, int make_colu)
{
	int i, j, k, thisNode;

	neq = _dsa->size();

// We build a temporary dof to Node table.
	int* dofToN = new int[neq]; //HB
	int numNodes = _dsa->numNodes();

	int *firstDOF   = (int *) dbg_alloca(sizeof(int)*numNodes);
	int *nodeWeight = (int *) dbg_alloca(sizeof(int)*numNodes);

	numUncon  = 0;
	unconstrNum.resize(neq);
	for(i = 0; i < neq; ++i) {
		unconstrNum[i] = (rCN) ? rCN[i] : i;
		if(unconstrNum[i] >= 0) numUncon++;
	}

	for(i=0; i < numNodes; ++i) {
		int fdof = _dsa->firstdof(i);
		int ndof = _dsa->weight(i);
		nodeWeight[i] = 0;
		for(j=0; j<ndof; ++j) {
			dofToN[j+fdof] = i;
			if(unconstrNum[j+fdof] >= 0) {
				nodeWeight[i]++;
			}
		}
	}

	for(i=0; i < numNodes; ++i)
		firstDOF[i] = -1;

	for(i = 0; i < neq; ++i) {
		int lDof = unconstrNum[i];
		if(lDof >= 0) {
			int lNode = dofToN[i];
			if(firstDOF[lNode] < 0 || lDof < firstDOF[lNode])
				firstDOF[lNode]  = lDof;
		}
	}

	// We assume that any Dof that has been used has a diagonal term.
	xunonz.resize(numUncon+1);
	xunonz[0] = 1;

	for(i=0; i<neq; ++i) {
		k = unconstrNum[i];     //  or k = unconstrndNum[i]
		if(k == -1) continue;
		xunonz[k+1]  = xunonz[k];
		thisNode = dofToN[i];
		for(j=0; j<cn->num(thisNode); ++j) {
			int jnode = (*cn)[thisNode][j];
			int fjdof = firstDOF[jnode];
			if(fjdof < 0 || fjdof > k) continue;
			if(nodeWeight[jnode] + fjdof <= k)
				xunonz[k+1] += nodeWeight[jnode];
			else
				xunonz[k+1] += k - fjdof + 1;
		}
	}

// Allocate memory for rowu (row numbers)
	rowu.resize(xunonz[numUncon]);
	for(i=0; i<neq; ++i) {
		int m = unconstrNum[i];
		if(m == -1) continue;
		int numFound = 0;
		thisNode = dofToN[i];
		for(j=0; j<cn->num(thisNode); ++j) {
			int jnode = (*cn)[thisNode][j];
			int fjdof = firstDOF[jnode];
			if(fjdof < 0 || fjdof > m) continue;
			for(k=0; k<nodeWeight[jnode]; ++k)
				if(fjdof + k < m) rowu[xunonz[m]-1+numFound++] = fjdof + k + 1;
		}
		rowu[xunonz[m]-1+numFound++] = m + 1; // This is for the diagonal terms.
	}

/* KAS.
   Construct colu : column numbers for each element - used for Mumps
   same size as rowu. use xunonz index to get column number.*/
	if(make_colu) {
		colu.resize(xunonz[numUncon]);
		k     = 0;
		for (i = 0; i < numUncon; i++) {
			// number of off-diagonal elements, will repeat column index numOffDiag times in colu
			int numOffDiag   = xunonz[i+1] - xunonz[i];   // xunonz is of size numUncon +1
			for (j = 0; j < numOffDiag; j++) {
				colu[k++]  = i+1;
			}
		}
	}
	if(dofToN) delete [] dofToN; //HB
}

// used by 2nd Mumps constructor (with make_colu = 1) and 2nd Spooles constructor
SparseData::SparseData(const DofSetArray *_dsa, const DofSetArray *c_dsa,
                       const Connectivity *cn, int expand, int make_colu, bool unsym)
{

	int i, j, k, thisNode;

	neq = _dsa->size(); // Total number of equations
	numUncon = c_dsa->size(); // Number of unconstrained dof (equations)

	// We build a temporary dof to Node table.
	int* dofToN = new int[neq]; //HB
	int numNodes = _dsa->numNodes();

	for(i=0; i < numNodes; ++i) {
		int fdof = _dsa->firstdof(i);
		int ndof = _dsa->weight(i);
		for(j=0; j<ndof; ++j)
			dofToN[j+fdof] = i;
	}

	unconstrNum = c_dsa->getUnconstrNum();

	// We assume that any Dof that has been used has a diagonal term.
	// Allocate memory for the diagonal location pointers
	xunonz.resize(numUncon+1);

	// Set the diagonal location pointers
	xunonz[0] = 1;
	for(i=0; i<neq; ++i) {
		k = unconstrNum[i];
		if(k == -1) continue;
		xunonz[k+1] = xunonz[k];
		thisNode = dofToN[i];
		for(j=0; j<cn->num(thisNode); ++j) {
			int jnode = (*cn)[thisNode][j];
			int fjdof = c_dsa->firstdof(jnode);
			if(unsym) {
				if(fjdof < 0) continue;
				xunonz[k+1] += c_dsa->weight(jnode);
			}
			else {
				if(fjdof < 0 || fjdof > k) continue;
				if(c_dsa->weight(jnode) + fjdof <= k)
					xunonz[k+1] += c_dsa->weight(jnode);
				else
					xunonz[k+1] += k - fjdof + 1;
			}
		}
	}

	// Allocate memory for rowu (row numbers)
	rowu.resize(xunonz[numUncon]);

	// Set the row numbers
	for(i=0; i<neq; ++i) {
		int m = unconstrNum[i];
		if(m == -1) continue;
		int numFound = 0;
		thisNode = dofToN[i];
		for(j=0; j<cn->num(thisNode); ++j) {
			int jnode = (*cn)[thisNode][j];
			int fjdof = c_dsa->firstdof(jnode);
			if(unsym) {
				if(fjdof < 0) continue;
				for(k=0; k<c_dsa->weight(jnode); ++k)
					rowu[xunonz[m]-1+numFound++] = fjdof + k + 1;
			}
			else {
				if(fjdof < 0 || fjdof > m) continue;
				for(k=0; k<c_dsa->weight(jnode); ++k)
					if(fjdof + k < m) rowu[xunonz[m]-1+numFound++] = fjdof + k + 1;
			}
		}
		if(!unsym) rowu[xunonz[m]-1+numFound++] = m + 1; // This is for the diagonal terms.
	}

	delete [] dofToN; //HB

	if(expand && numUncon > 0) { // PJSA
		int k,j;
		std::vector<int> new_xunonz(numUncon+1);
		int* counter = new int[numUncon+1]; //HB
		for(k=0; k < numUncon+1; k++) {
			new_xunonz[k] = xunonz[k];
			counter[k] = 0;
		}
		for(k=0; k < numUncon; k++)
			for(j=xunonz[k]-1; j < xunonz[k+1]-1; j++) {
				if(rowu[j]-1 != k)
					counter[rowu[j]]++;
			}
		for(k=1; k < numUncon; k++)
			counter[k] += counter[k-1];
		for(k=1; k < numUncon; k++)
			new_xunonz[k] += counter[k];
		new_xunonz[numUncon] += counter[numUncon-1];

		for(k=0; k < numUncon; k++)
			counter[k] = 0;

		std::vector<int> new_rowu(new_xunonz[numUncon]);

//  map the other side to rowu and unonz
		for(k=0; k < numUncon; k++) {
			for(j=xunonz[k]-1; j < xunonz[k+1]-1; j++) {
				new_rowu [new_xunonz[k]-1+counter[k]] = rowu[j];
				counter[k]++;
				if(rowu[j]-1 != k) {
					new_rowu [ new_xunonz[rowu[j]-1]-1+counter[rowu[j]-1]] = k+1;
					counter[rowu[j]-1]++;
				}
			}
		}

		rowu   = std::move(new_rowu);
		xunonz = std::move(new_xunonz);
		delete [] counter; //HB
	}


/* KAS.
   Construct colu : column numbers for each element - used for Mumps
   same size as rowu. use xunonz index to get column number.*/
	if(make_colu) {
		colu.resize(xunonz[numUncon]);
		k     = 0;
		for (i = 0; i < numUncon; i++) {
			// number of off-diagonal elements, will repeat column index numOffDiag times in colu
			int numOffDiag   = xunonz[i+1] - xunonz[i];   // xunonz is of siez numUncon +1
			for (j = 0; j < numOffDiag; j++) {
				colu[k++]  = i+1;
			}
		}
	}
}


SparseData::SparseData(const DofSetArray *_dsa, const int *glInternalMap,
                       const Connectivity *cn, int expand)
{
	int i, j, k, thisNode;

	neq  =  _dsa->size(); // Total number of equations

	// Compute number of Boundary dofs and number of Internal dofs
	int internalLen = 0;
	int iDof;
	for(iDof = 0; iDof < neq; ++iDof) {
		if(glInternalMap[iDof] >= 0) internalLen += 1;
	}
	numUncon = internalLen;

	int numNodes = _dsa->numNodes();
	int *firstDOF   = (int *) dbg_alloca(sizeof(int)*numNodes);
	int *nodeWeight = (int *) dbg_alloca(sizeof(int)*numNodes);

// We build a temporary dof to Node table.
	int* dofToN = new int[neq]; //HB

	unconstrNum.assign(glInternalMap, glInternalMap+neq); // TODO remove the const cast.

	for(i=0; i < numNodes; ++i) {
		int fdof = _dsa->firstdof(i);
		int ndof = _dsa->weight(i);
		nodeWeight[i] = 0;
		for(j=0; j<ndof; ++j) {
			dofToN[j+fdof] = i;
			if(unconstrNum[j+fdof] >= 0) {
				nodeWeight[i]++;
			}
		}
	}

	for(i=0; i < numNodes; ++i)
		firstDOF[i] = -1;

	for(i = 0; i < neq; ++i) {
		int lDof = unconstrNum[i];
		if(lDof >= 0) {
			int lNode = dofToN[i];
			if(firstDOF[lNode] < 0 || lDof < firstDOF[lNode])
				firstDOF[lNode]  = lDof;
		}
	}

	// We assume that any Dof that has been used has a diagonal term.
	// Allocate memory for the diagonal location pointers
	xunonz.resize(numUncon+1);

	// Set the diagonal location pointers
	xunonz[0] = 1;
	for(i=0; i<neq; ++i) {
		k = unconstrNum[i];
		if(k == -1) continue;
		xunonz[k+1] = xunonz[k];
		thisNode = dofToN[i];
		for(j=0; j<cn->num(thisNode); ++j) {
			int jnode = (*cn)[thisNode][j];
			int fjdof = firstDOF[jnode];
			if(fjdof < 0 || fjdof > k) continue;
			if(nodeWeight[jnode] + fjdof <= k)
				xunonz[k+1] += nodeWeight[jnode];
			else
				xunonz[k+1] += k - fjdof + 1;
		}
	}

	// Allocate memory for rowu (row numbers)
	rowu.resize(xunonz[numUncon]);

	// Set the row numbers
	for(i=0; i<neq; ++i) {
		int m = unconstrNum[i];
		if(m == -1) continue;
		int numFound = 0;
		thisNode = dofToN[i];
		for(j=0; j<cn->num(thisNode); ++j) {
			int jnode = (*cn)[thisNode][j];
			int fjdof = firstDOF[jnode];
			if(fjdof < 0 || fjdof > m) continue;
			for(k=0; k<nodeWeight[jnode]; ++k) {
				if(fjdof + k < m) rowu[xunonz[m]-1+numFound++] = fjdof + k + 1;
			}
		}
		rowu[xunonz[m]-1+numFound++] = m + 1; // This is for the diagonal terms.
	}
	if(dofToN) delete [] dofToN; //HB

	if(expand) {
		int k,j;
		std::vector<int> new_xunonz(numUncon+1);
		std::vector<int> counter(numUncon+1); //HB

		for(k=0; k < numUncon+1; k++) {
			new_xunonz[k] = xunonz[k];
			counter[k] = 0;
		}

		for(k=0; k < numUncon; k++)
			for(j=xunonz[k]-1; j < xunonz[k+1]-1; j++) {
				if(rowu[j]-1 != k)
					counter[rowu[j]]++;
			}

		for(k=1; k < numUncon; k++)
			counter[k] += counter[k-1];

		for(k=1; k < numUncon; k++)
			new_xunonz[k] += counter[k];

		new_xunonz[numUncon] += counter[numUncon-1];

		for(k=0; k < numUncon; k++)
			counter[k] = 0;

		std::vector<int> new_rowu(new_xunonz[numUncon]);

		// map the other side to rowu and unonz
		for (k=0; k < numUncon; k++) {
			for(j=xunonz[k]-1; j < xunonz[k+1]-1; j++) {
				new_rowu [new_xunonz[k]-1+counter[k]] = rowu[j];
				counter[k]++;
				if(rowu[j]-1 != k) {
					new_rowu [ new_xunonz[rowu[j]-1]-1+counter[rowu[j]-1]] = k+1;
					counter[rowu[j]-1]++;
				}
			}
		}

		rowu   = std::move(new_rowu);
		xunonz = std::move(new_xunonz);
	}

}

#include <algorithm> // PJSA: for sgi intel

SparseData::SparseData(const EqNumberer *eqn, const Connectivity *icn)
{
	int i, j, k, l;
	Connectivity ccn{*icn};
	ccn.sortTargets(); // PJSA: for sgi intel ML: Do we still need this?
	Connectivity *cn = &ccn;

	numUncon = neq  = eqn->size(); // Total number of equations

	int numNodes = cn->csize();
	std::vector<int> nodeCWeight(numNodes);

	for(i = 0; i < numNodes; ++i) {
		nodeCWeight[i] = 0;
		for(j = 0; j < cn->num(i); ++j)
			nodeCWeight[i] += eqn->weight((*cn)[i][j]);
	}

	// Allocate the column pointer
	xunonz.resize(neq+1);

	for(i = 0; i < numNodes; ++i) {
		int fDof = eqn->firstdof(i);
		int nDof = eqn->weight(i);
		for(j = 0; j < nDof; ++j)
			xunonz[fDof+j] = nodeCWeight[i];
	}

	int count = 0;
	for(i = 0; i < neq; ++i) {
		int tmpc = count;
		count += xunonz[i];
		xunonz[i] = tmpc;
	}
	xunonz[neq] = count;

	rowu.resize(count);

	for(i = 0; i < numNodes; ++i) {
		int fDof = eqn->firstdof(i);
		int nDof = eqn->weight(i);
		int index = 0;
		for(k = 0; k < cn->num(i); ++k) {
			int nodeNum = (*cn)[i][k];
			int fkDof = eqn->firstdof(nodeNum);
			int nkDof = eqn->weight(nodeNum);
			for(l = 0; l < nkDof; ++l, ++index)
				for(j = 0; j < nDof; ++j)
					rowu[ xunonz[fDof+j] + index] = fkDof+l;
		}
	}
}

SparseData::SparseData(const Connectivity *cn, const EqNumberer *eqn, double trbm,
                       bool expand)
{
	neq  = eqn->size(); // Total number of equations
	numUncon = eqn->size();

	unconstrNum.resize(neq);
	constrndNum.resize(neq);

	for (int i=0; i<neq; ++i) {
		unconstrNum[i] = i;
		constrndNum[i] = -1;
	}

	// We build a temporary dof to Node table.
	int numNodes = eqn->numNodes();
	if( cn->csize() < numNodes ) numNodes = cn->csize();

	std::vector<int> dofToN(neq);

	for (int i=0; i < numNodes; ++i) {
		int fdof = eqn->firstdof(i);
		int ndof = eqn->weight(i);
		for (int j=0; j<ndof; ++j)
			dofToN[j+fdof] = i;
	}

	// We assume that any Dof that has been used has a diagonal term.
	// Allocate memory for the diagonal location pointers
	xunonz.resize(neq+1);

	// Set the diagonal location pointers
	xunonz[0] = 1;

	for (int i=0; i<neq; ++i) {
		xunonz[i+1] = xunonz[i];
		auto thisNode = dofToN[i];
		for (int j=0; j<cn->num(thisNode); ++j) {
			int jnode = (*cn)[thisNode][j];
			int fjdof = eqn->firstdof(jnode);
			if(fjdof < 0 || fjdof > i) continue;
			if( eqn->weight(jnode) + fjdof <= i )
				xunonz[i+1] += eqn->weight(jnode);
			else
				xunonz[i+1] += i - fjdof + 1;
		}
	}

	// Allocate memory for rowu (row numbers)
	rowu.resize(xunonz[neq]-1);

	// Set the row numbers
	for (int i=0; i<neq; ++i) {
		int numFound = 0;
		auto thisNode = dofToN[i];
		for (int j=0; j<cn->num(thisNode); ++j) {
			int jnode = (*cn)[thisNode][j];
			int fjdof = eqn->firstdof(jnode);

			if(fjdof < 0 || fjdof > i) continue;

			for (int k=0; k<eqn->weight(jnode); ++k)
				if(fjdof + k < i) rowu[xunonz[i]-1+numFound++] = fjdof + k + 1;
		}
		rowu[xunonz[i]-1+numFound++] = i + 1; // This is for the diagonal terms.
	}

	if(expand) {
		std::vector<int> new_xunonz(numUncon+1);
		std::vector<int> counter(numUncon+1); //HB

		for (int k=0; k < numUncon+1; k++)  {
			new_xunonz[k] = xunonz[k];
			counter[k] = 0;
		}

		for (int k=0; k < numUncon; k++)
			for (int j=xunonz[k]-1; j < xunonz[k+1]-1; j++) {
				if(rowu[j]-1 != k)
					counter[rowu[j]]++;
			}

		for (int k=1; k < numUncon; k++)
			counter[k] += counter[k-1];

		for (int k=1; k < numUncon; k++)
			new_xunonz[k] += counter[k];

		new_xunonz[numUncon] += counter[numUncon-1];

		for (int k=0; k < numUncon; k++)
			counter[k] = 0;

		std::vector<int> new_rowu(new_xunonz[numUncon]-1);

		// map the other side to rowu and unonz
		for (int k=0; k < numUncon; k++) {
			for (int j=xunonz[k]-1; j < xunonz[k+1]-1; j++) {
				new_rowu [new_xunonz[k]-1+counter[k]] = rowu[j];
				counter[k]++;

				if(rowu[j]-1 != k) {
					new_rowu [ new_xunonz[rowu[j]-1]-1+counter[rowu[j]-1]] = k+1;
					counter[rowu[j]-1]++;
				}
			}
		}

		rowu   = std::move(new_rowu);
		xunonz = std::move(new_xunonz);
	}

}


SparseData::SparseData(LMPCons **mpc, int numMPC, const DofSetArray *c_dsa)
{
	numConstrained = numMPC;
	numUncon = c_dsa->size();
	xunonz.resize(numConstrained+1);

	int i;
	xunonz[0] = 0;
	for(i=1; i<numConstrained+1; ++i) {
		xunonz[i] = xunonz[i-1] + mpc[i-1]->nterms;
	}

	rowu.resize(xunonz[numConstrained]);

	int mstart, mstop, m;
	for(i=0; i<numConstrained; ++i) {
		mstart = xunonz[i];
		mstop  = xunonz[i+1];
		for(m=mstart; m<mstop; ++m) {
			rowu[m] = c_dsa->locate(mpc[i]->terms[m].nnum,
			                        (1 << mpc[i]->terms[m].dofnum));
		}
	}
}


SparseData::SparseData(int num, const int *xyzCount, const int *xyzList)
{
	numConstrained = num;
	xunonz.resize(numConstrained+1);

	xunonz[0] = 0;
	int i;
	for(i=0; i<numConstrained; ++i)
		xunonz[i+1] = xunonz[i] + xyzCount[i];

	rowu.assign(xyzList, xyzList+xunonz[num]);
}

//KHP:to store local G as a rectangular sparse matrix

SparseData::SparseData(int numInterface,
                       const int *glbmap, int numModes, int ldm)
{
	numConstrained = numModes;
	xunonz.resize(numConstrained+1);
	int i;
	xunonz[0] = 0;
	for(i=1; i<numConstrained+1; ++i) {
		xunonz[i] = xunonz[i-1] + numInterface;
	}
	int count = numInterface*numConstrained;
	rowu.resize(count);

	int mstart, mstop, m;
	for(i=0; i<numConstrained; ++i) {
		mstart = xunonz[i];
		mstop  = xunonz[i+1];
		for(m=mstart; m<mstop; ++m) {
			rowu[m] = glbmap[m-mstart];
		}
	}
}

