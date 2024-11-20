#include <cstdlib>
#include <Utils.d/dbg_alloca.h>

#include <Utils.d/linkfc.h>
#include <Utils.d/Connectivity.h>
#include <Math.d/Skyline.d/SkyMatrix.h>
#include <Math.d/Skyline.d/utility.h>
#include <Solvers.d/Rbm.h>

extern "C"      {
void _FORTRAN(slacol)( int*, int*, int&, int&, int& );
}


SkyData::SkyData()
{
   initialize();
}

void
SkyData::initialize()
{
   lacol = 0;
   pivot = 0;
   seqid = 0;
   dlp = 0;
   myRCN = 0;
   rowColNum = 0;
   neq = 0;
   numUncon = 0;
   nzem = 0;
   rbm = 0;
   isTRBM = 0;
   TOLERANCE = 1.0E-6;
   myRbm = 0;
}

SkyData::~SkyData()
{
   if(lacol) { delete [] lacol; lacol=0; }
   if(pivot) { delete [] pivot; pivot=0; }
   if(seqid) { delete [] seqid; seqid=0; }
   if(dlp) { delete [] dlp; dlp=0; }
   if(myRCN && rowColNum) { delete [] rowColNum; rowColNum=0; }
   if(myRbm && rbm) { delete rbm; rbm=0; }
}

SkyData::SkyData(const Connectivity *cn, const DofSetArray *c_dsa, double trbm,
                 Rbm *rigid)
{
  initialize();
  myRbm = 0;

  int i,j,k,n;

  // ... SET TOLERANCE
  TOLERANCE = trbm;

  // ... GET NUMBER OF UNCONSTRAINED DOF
  numUncon = c_dsa->size();

  // ... GET NUMBER OF NODES
  int numNodes = c_dsa->numNodes();
  if( cn->csize() < numNodes ) numNodes = cn->csize();

  // ... GET THE UNCONSTRAINED NUMBER ARRAY
  rowColNum = c_dsa->getUnconstrNum().data();

  if(numUncon==0) return;

  // Allocate enough space for each of the int arrays in this constructor
  dlp = new int[numUncon];

  // Initialize first element of diagonal location pointers to one.
  dlp[0] = 1;

  // Build the dof to node (dofToN) table from constrained dsa.
  int *dofToN  = (int *) dbg_alloca(sizeof(int)*numUncon);

  for(i = 0; i < numNodes; ++i) {
    int fdof = c_dsa->firstdof(i);
    int ndof = c_dsa->weight(i);
    for(j = 0; j < ndof; ++j)
      dofToN[j + fdof] = i;
  }

  int compNum, numComp;

  // Get the number of components
  if(rigid) 
    numComp = rigid->numComponents();
  else
    numComp = 1;

  // Looping over Equations to find minimum Dof connection.
  for(compNum = 0; compNum < numComp; ++compNum) {
     int firstDof, lastDof;
     firstDof = (rigid && numComp > 1) ? rigid->firstDof(compNum)   : 0;
     lastDof  = (rigid && numComp > 1) ? rigid->firstDof(compNum+1) : numUncon;
     int offset = firstDof;
     if(firstDof == 0) firstDof = 1;
     for(n = firstDof; n < lastDof; ++n) {
       int thisNode = dofToN[n];
       int minDof   = n;
       for(i = 0; i < cn->num(thisNode); ++i) {
         k = c_dsa->firstdof((*cn)[thisNode][i]); 
         if(k >= 0 && k < minDof) minDof = k;
       }
       minDof = minDof - ((minDof-offset) % 4);
       dlp[n] = dlp[n-1] +  n - minDof + 1;
       if(dlp[n] < dlp[n-1]) {
         fprintf(stderr,"*** ERROR: Integer Overflow, exiting\n"); 
         fflush(stderr);
         exit(-1);
       }
     }
  }

  int MAXLAC, AVELAC;

  // ... ALLOCATE LAST ACTIVE COLUMN (LACOL) VECTOR
  lacol = new int[numUncon];

  // ... CONSTRUCT LACOL
  _FORTRAN(slacol)( dlp, lacol, numUncon, MAXLAC, AVELAC );

  // ... ALLOCATE SPACE FOR pivot ARRAY
  pivot = new int [numUncon];

  // INITIALIZE rbm
  rbm = rigid;

  // SET NUMBER OF RBM
  if(rbm) { 
    nzem = rbm->numRBM();
    seqid = new int[numUncon];
  }
}

SkyData::SkyData(const Connectivity *cn, const EqNumberer *_dsa, double trbm, const int *bc)
{
  initialize();

  TOLERANCE = trbm;

  int i, j, n, k;
  neq = _dsa->size();
  int numNodes = _dsa->numNodes();
  if( cn->csize() < numNodes ) numNodes = cn->csize();

// Count & renumber unconstrained dofs & mark constrained dofs
  myRCN = 1;
  auto build_rowColNum = new int[neq];
  makeUnconstrainedNum(neq, bc, build_rowColNum, numUncon);
  rowColNum = build_rowColNum;

// Build the dofTON in unconstrained numbering
  int *dofToN = (int *) dbg_alloca(sizeof(int)*numUncon);

  for(i = 0; i < numNodes; ++i ) {
    int fdof = _dsa->firstdof(i);
    int ndof = _dsa->weight(i);
    for(j = 0; j < ndof; ++j ) {
      if(rowColNum[j+fdof] >= 0) dofToN[rowColNum[j + fdof]] = i;
    }
  }

  int *uncFirstDof = (int *) dbg_alloca(sizeof(int)*numNodes);

  for(i=0; i < numNodes; ++i)
    uncFirstDof[i] = -1;

  for(i=0; i < numUncon; ++i) {
    int nodeNum = dofToN[i];
    if(uncFirstDof[ nodeNum ] < 0)
       uncFirstDof[ nodeNum ] = i;
  }

// Allocate memory for diagonal location pointer

  dlp        = new int[numUncon];

// Initialize first element of diagonal location pointers to one.

  dlp[0] = 1;

// Loop over Equations to find minimum Dof connection.

  //ML Debuging
  int maxdist = 0;
  for(n = 1; n < numUncon; ++n) {
    int thisNode = dofToN[n];
    int minDof   = n;
    for(i = 0; i < cn->num(thisNode); ++i) {
      k = uncFirstDof[(*cn)[thisNode][i] ];
      if( k >= 0 && k < minDof) minDof = k;
    }
    if(n-minDof > maxdist) maxdist = n-minDof;
    minDof = minDof - (minDof % 4);
    dlp[n] = dlp[n-1] + n - minDof + 1;
    if(dlp[n] < dlp[n-1]) {
      fprintf(stderr,"DLP OVERFLOW\n"); fflush(stderr);
      exit(-1);
    }
  }

// Create lacol vector
   int MAXLAC,AVELAC;
   lacol = new int[numUncon];
   _FORTRAN(slacol)( dlp, lacol, numUncon, MAXLAC, AVELAC );

   pivot = new int[numUncon];
}

SkyData::SkyData(const Connectivity *cn, const EqNumberer *eqn, double trbm)
{
  initialize();

  TOLERANCE = trbm;

  int i, j, n, k;
  neq = eqn->size();
  int numNodes = eqn->numNodes();
  if( cn->csize() < numNodes ) numNodes = cn->csize();
  
  int maxdist = 0;
  numUncon = eqn->size();
  myRCN = 1;
  auto build_rowColNum = new int[numUncon];
  for(i=0; i<numUncon; ++i)
	  build_rowColNum[i] = i;

  rowColNum = build_rowColNum;
  if(numUncon==0) return;

  // Allocate memory for diagonal location pointer
  dlp        = new int[numUncon];
  // For each DOF, compute in dlp[dof] the column height
  for(i=0; i < numNodes; ++i) {
     int minDof = neq;
     for(j = 0; j < cn->num(i); ++j ) {
       k = eqn->firstdof((*cn)[i][j]);
       if(k < 0) continue;
       if(k < minDof) minDof = k;
     }
     minDof = minDof - (minDof % 4);  // PJSA: this seemed to be causing bug
     int fDof = eqn->firstdof(i);
     for(j = 0; j < eqn->weight(i); ++j) {
       n = fDof+j;
       if(n<0) continue; // PJSA
       if(n-minDof > maxdist) maxdist = n-minDof;
       dlp[n] = n - minDof + 1;
     }
  }

  // Now sum the column heights to get dlp
  for(i = 1; i < numUncon; ++i) {
     dlp[i] += dlp[i-1];
     if( dlp[n] < 0) {
       fprintf(stderr,"DLP OVERFLOW\n"); fflush(stderr);
       exit(-1);
     }
  }

 // Create lacol vector
  int MAXLAC,AVELAC;
  lacol = new int[numUncon];
  _FORTRAN(slacol)( dlp, lacol, numUncon, MAXLAC, AVELAC );

  pivot = new int[numUncon];
}

SkyData::SkyData(const EqNumberer *_dsa, const Connectivity *cn, double trbm, const int *rCN )
{
  initialize();

  TOLERANCE = trbm;

  int i, j, n, k;
  neq = _dsa->size();
  int numNodes = _dsa->numNodes();
  if(cn->csize() < numNodes) numNodes = cn->csize();

// Count & renumber unconstrained dofs & mark constrained dofs
	numUncon  = 0;
	myRCN = 1;
	auto build_rowColNum = new int[neq];
	for(i = 0; i < neq; ++i) {
		build_rowColNum[i] = rCN[i];
		if(build_rowColNum[i] >= 0) numUncon++;
	}
	rowColNum = build_rowColNum;

// Build the dofTON in unconstrained numbering
  int *dofToN = (int *) dbg_alloca(sizeof(int)*numUncon);

  for(i = 0; i < numNodes; ++i) {
    int fdof = _dsa->firstdof(i);
    int ndof = _dsa->weight(i);
    for(j = 0; j < ndof; ++j ) {
      if(rowColNum[j+fdof] >= 0) dofToN[rowColNum[j + fdof]] = i;
    }
  }

  int *uncFirstDof = (int *) dbg_alloca(sizeof(int)*numNodes);

  for(i=0; i < numNodes; ++i)
    uncFirstDof[i] = -1;

  for(i=0; i < numUncon; ++i) {
    int nodeNum = dofToN[i];
    if(uncFirstDof[ nodeNum ] < 0)
       uncFirstDof[ nodeNum ] = i;
   }

// Allocate memory for diagonal location pointer
  dlp        = new int[numUncon];

// Initialize first element of diagonal location pointers to one.
  dlp[0] = 1;

// Loop over Equations to find minimum Dof connection.
  for(n = 1; n < numUncon; ++n) {
    int thisNode = dofToN[n];
    int minDof   = n;
    for(i = 0; i < cn->num(thisNode); ++i) {
      k = uncFirstDof[(*cn)[thisNode][i] ];
      if( k >= 0 && k < minDof) minDof = k;
    }
    minDof = minDof - (minDof % 4);
    dlp[n] = dlp[n-1] + n - minDof + 1;
    if(dlp[n] < dlp[n-1]) {
      fprintf(stderr,"DLP OVERFLOW\n"); fflush(stderr);
      exit(-1);
    }
  }

// Create lacol vector
   int MAXLAC,AVELAC;
   lacol = new int[numUncon];
   _FORTRAN(slacol)( dlp, lacol, numUncon, MAXLAC, AVELAC );

   pivot = new int[numUncon];
}

SkyData::SkyData(int n, double tolerance)
{
   initialize();

   TOLERANCE = tolerance;
   
   numUncon = n;

   lacol     = new int[n];
   pivot     = new int[n];
   dlp       = new int[n];
   myRCN = 1;
   rowColNum = new int[n];
   if(numUncon==0) return;

   dlp[0] = 1;
   for(int i=1; i< n; i++) {
     dlp[i] = dlp[i-1] + i + 1;
   }

   int MAXLAC, AVELAC;

   _FORTRAN(slacol)(dlp, lacol, numUncon, MAXLAC, AVELAC );
}

