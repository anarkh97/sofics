#if defined(USE_SCALAPACK) && defined(USE_EIGEN3)
#include "DistrNonnegativeMatrixFactorization.h"

#include <Comm.d/Communicator.h>
#include <Rom.d/SparseSolvers.d/ScalaLH.d/Plh.h>
#include <Rom.d/NmfSolvers.d/ScalaPQN.d/Nmf.h>
#include <Utils.d/linkfc.h>
#include <Timers.d/GetTime.h>

#include "DistrNodeDof6Buffer.h"
#include "DistrVecNodeDof6Conversion.h"
#include "DistrSvdOrthogonalization.h"

#include <iostream>
#include <map>
#include <algorithm>

#include <mpi.h>

extern "C" {
  // Context & cpu topology management
  int Csys2blacs_handle(MPI_Comm comm);

  void Cfree_blacs_system_handle(int handle);
  void Cblacs_gridinit(int *ictxt, char *order, int nprow, int npcol);

  void Cblacs_gridinfo(int ictxt, int *nprow, int *npcol, int *myrow, int *mycol);
  void Cblacs_gridexit(int ictxt);

  // Index mapping 
  int _FORTRAN(indxg2l)(const int *indxglob, const int *nb, const int *iproc_dummy, const int *isrcproc_dummy, const int *nprocs);
  int _FORTRAN(indxg2p)(const int *indxglob, const int *nb, const int *iproc_dummy, const int *isrcproc, const int *nprocs);
}

extern int verboseFlag;

namespace Rom {

double t1=0,t2=0,t3=0,t4=0,t5=0,t6=0;

DistrNonnegativeMatrixFactorization
::DistrNonnegativeMatrixFactorization(Communicator * comm, int rowCount, int colCount, int localRows, int basisDimension,
                                      int blockSize, int maxIter, double tol, int method, int nsub, int pqnNumInnerIter, double pqnAlpha) : 
  communicator_(comm),
  rowCount_(rowCount),
  colCount_(colCount),
  localRows_(localRows),
  basisDimension_(basisDimension),
  blockSize_(blockSize),
  maxIter_(maxIter),
  tol_(tol),
  method_(method),
  nsub_(nsub),
  matrixBuffer_(localRows,colCount),
  basisBuffer_(localRows,basisDimension),
  pqnNumInnerIter_(pqnNumInnerIter),
  pqnAlpha_(pqnAlpha),
  energy_(0.0)
{
}

int 
DistrNonnegativeMatrixFactorization
::energySVD(double energy, std::vector<int> rows, std::vector<int> cols)
{
  // initialize distributed svd
  DistrSvdOrthogonalization solver(communicator_, communicator_->numCPUs(), 1);
  const int blockSize = domain->solInfo().svdBlockSize; // default: 64
  {
   solver.blockSizeIs(blockSize);
  }

  // allocate space
  const int localLength = rows.size(); 
  {
    const int maxLocalLength = communicator_->globalMax(localLength);
    const int globalProbSize = ((maxLocalLength/blockSize+1)*solver.rowCpus()+1)*blockSize+1;
    solver.problemSizeIs(globalProbSize, cols.size());
  }

  // read in non-zero rows
  int colCounter = 0;
  for(std::vector<int>::iterator colit = cols.begin(); colit != cols.end(); colit++) {
    double *vecBuffer = solver.matrixColBuffer(colCounter); colCounter++; 
    int rowCounter = 0;
    for(std::vector<int>::iterator rowit = rows.begin(); rowit != rows.end(); rowit++) {
      vecBuffer[rowCounter] = matrixBuffer_(*rowit,*colit); rowCounter++;
    }
  }

  solver.solve();

  int rowCount = communicator_->globalSum((int)rows.size());  
  int numSV = std::min(rowCount,(int)cols.size());
  std::vector<double> toto(numSV+1);
  toto[numSV] = 0;

  // use singular values to determine basis size
  for (int iVec = numSV-1; iVec >= 0; --iVec) {
    toto[iVec] = toto[iVec+1]+solver.singularValue(iVec); // running sum
  }

  bool reset = true;
  for (int iVec = 0; iVec < numSV; ++iVec) {
    double en = toto[iVec]/toto[0];
    if(en < energy && reset){
      numSV = iVec+1;
      reset = false;
    }
    if(communicator_->myID() == 0) 
      std::cout << iVec+1 << " " << solver.singularValue(iVec) << " " << en << std::endl;
  }

  return numSV;
 
}

void
DistrNonnegativeMatrixFactorization::solve()
{
  // data is partitioned row-wise
  // determine which cols which are non-zero
  int *buffer = new int[colCount_];
  for(int i=0; i<colCount_; ++i) buffer[i] = (matrixBuffer_.col(i).array() == 0).all() ? 0 : 1;
  communicator_->globalMax(colCount_, buffer);

  // count the number of non-zero cols
  std::vector<int> cols;
  for(int i=0; i<colCount_; ++i) if(buffer[i]) cols.push_back(i);
  int colCount = cols.size();
  
  // print to screen
  if(communicator_->myID() == 0) std::cerr << "X has " << colCount_ << " columns of which " << colCount << " are non-zero\n";
  delete [] buffer;

  // determine which rows are non-zero
  std::vector<int> rows;
  for(int i=0; i<localRows_; ++i) if(!(matrixBuffer_.row(i).array() == 0).all()) rows.push_back(i);
  
  // perform a global sum to determine total number of non-zero rows
  int rowCount = communicator_->globalSum((int)rows.size());
  if(communicator_->myID() == 0) std::cerr << "X has " << rowCount_ << " rows of which " << rowCount << " are non-zero\n";

  if(energy_ > 0.0){
    if(communicator_->myID() == 0 ) std::cout << "Computing Singular Value cutoff" << std::endl;

    basisDimension_ = energySVD(energy_,rows,cols);
    basisBuffer_.resize(localRows_,basisDimension_);
    std::cout << "Basis Dimension is " << basisDimension_ << std::endl;

  }

  // 1. Construct matrix A
  int blacsHandle = Csys2blacs_handle(*communicator_->getCommunicator());
  int context = blacsHandle;
  char order[] = "R";
  int rowCpus = communicator_->numCPUs();
  int colCpus = 1;
  Cblacs_gridinit(&context, order, rowCpus, colCpus);
  SCDoubleMatrix A(context, rowCount, colCount, blockSize_, blockSize_, *communicator_->getCommunicator());

  // 2. Copy non-zero rows and columns from matrixBuffer into A
  int myrow, mycol;
  int dummy, zero=0;
  Cblacs_gridinfo(context, &rowCpus, &colCpus, &myrow, &mycol);
  double *row = new double[std::max(colCount,basisDimension_)];
/*for(int j=1; j<=rowCount_; j++) {
    int p = _FORTRAN(indxg2p)(&j, &blockSize_, &dummy, &zero, &rowCpus);
    if(myrow == p) {
      int lj = _FORTRAN(indxg2l)(&j, &blockSize_, NULL, NULL, &rowCpus) - 1;
      for(int k=0; k<colCount_; ++k) row[k] = -matrixBuffer_(lj,k);
      A.setMatrixRow(j, row);
    }
  }*/
  // make a local context with just 1 cpu
  int localBlacsHandle = Csys2blacs_handle(MPI_COMM_SELF);
  int localContext = localBlacsHandle;
  Cblacs_gridinit(&localContext, order, 1, 1);
  // make a vector to store one non-zero row of A
  SCDoubleMatrix Ai(localContext, 1, colCount, blockSize_, blockSize_, MPI_COMM_SELF);
  int isNonZero;
  std::map<int,int> localrows;
  // loop over the rows, check if non-zero and if so then copy/redistribute to A
  // also, store mapping from global non-zero row numbers to local row numbers
  for(int j=1,i=1; j<=rowCount_; j++) {
    int p = _FORTRAN(indxg2p)(&j, &blockSize_, &dummy, &zero, &rowCpus);
    if(myrow == p) {
      int lj = _FORTRAN(indxg2l)(&j, &blockSize_, NULL, NULL, &rowCpus) - 1;
      std::vector<int>::iterator it = std::find(rows.begin(), rows.end(), lj);
      isNonZero = (it != rows.end()) ? 1 : 0;
      communicator_->broadcast(1, &isNonZero, p);
      if(isNonZero) {
        for(int k=0; k<colCount; ++k) row[k] = -matrixBuffer_(*it,cols[k]);
        Ai.setMatrixRow(1, row);
        Ai.copyRedist(1, colCount, 1, 1, A, i, 1, A.getContext()); // send and receive
        localrows[i] = lj;
        i++;
      }
    }
    else {
      communicator_->broadcast(1, &isNonZero, p);
      if(isNonZero) {
        SCDoubleMatrix::copyRedist(1, colCount, A, i, 1, A.getContext()); // receive only
        i++;
      }
    }
  }
  Cblacs_gridexit(localContext);
  Cfree_blacs_system_handle(localBlacsHandle);

  // Construct matrices to store factors W and H
  SCDoubleMatrix W(context, rowCount, basisDimension_, blockSize_, blockSize_, *communicator_->getCommunicator());
  SCDoubleMatrix Htranspose(context, colCount, basisDimension_, blockSize_, blockSize_, *communicator_->getCommunicator());

  // Initialize W with random positive entries between 0 and 1
  W.initRandom(1,0.,1.);

  switch(method_) {
    default: case 1 : {
      if(communicator_->myID() == 0) {
        std::cout << "Using Alternating Least Squares NMF method" << std::endl;
        std::cout << "   method_ = " << method_ << std::endl;
      }
      // Alternating non-negative least squares iteration loop
      for(int i=0; i<maxIter_; ++i) {

        // make a copy of W to use for checking the stopping criteria
        SCDoubleMatrix W_copy(W);

        // solve: min ||WH-A|| s.t. H >= 0
        solveNNLS_MRHS(W, A, Htranspose, 0);

        // solve: min ||H^TW^T-A^T|| s.t. W >= 0
        solveNNLS_MRHS(Htranspose, A, W, 1);

        // compute residual and check stopping criteria
        SCDoubleMatrix Err(A); W.multiply(Htranspose, Err, 'N', 'T', -1.0, 1.0); // Err = A-W*H;
        double res = Err.froNorm()/A.froNorm();
        W.add(W_copy, 'N', rowCount, basisDimension_, 1.0, -1.0); // W_copy = W-W_copy
        double inc = W_copy.froNorm();
        if(communicator_->myID() == 0)
          std::cout << "iteration = " << i+1 << ", rel. residual = " << res << ", solution incr. = " << inc << std::endl;
        if(inc < tol_) break;
      }
    } break;

    case 2 : {
      if(communicator_->myID() == 0)
        std::cout << "ERROR: Greedy method is not implemented in DistrNonnegativeMatrixFactorization.C\n";
      exit(-1);
    } break;

    case 3 : {
      // NMF based on PQN
      if(communicator_->myID() == 0) {
        std::cout << "Using NMF based on Projected Quasi-Newton" << std::endl;
      }
      Htranspose.initRandom(2,0.,1.);
      //Htranspose.zero();
      //A.initRandom(2,0.,1.);
      Nmf solver = Nmf(A, W, Htranspose);
      solver.setMaxIter(maxIter_);
      solver.setNumInnerIter(pqnNumInnerIter_);
      solver.setAlpha(pqnAlpha_);
      solver.setTol(tol_);
      solver.summary();
      solver.solve();
      solver.printTimes(true);
    } break;

  }

  // copy W into basisBuffer_
/*for(int j=1; j<=rowCount_; j++) {
    int p = _FORTRAN(indxg2p)(&j, &blockSize_, &dummy, &zero, &rowCpus);
    if(myrow == p) {
      int lj = _FORTRAN(indxg2l)(&j, &blockSize_, NULL, NULL, &rowCpus) - 1;
      W.getMatrixRow(j, row, 'U');
      for(int k=0; k<basisDimension_; ++k) basisBuffer_(lj,k) = row[k];
    }
  }*/
  basisBuffer_.setZero();
  for(int i=1; i<=rowCount; i++) {
    W.getMatrixRow(i, row, 'A');
    std::map<int,int>::iterator it = localrows.find(i);
    if(it != localrows.end()) {
      for(int k=0; k<basisDimension_; ++k) basisBuffer_(it->second,k) = row[k];
    }
  }

  delete [] row;

  Cblacs_gridexit(context);
  Cfree_blacs_system_handle(blacsHandle);

/*for(int i=0; i<communicator_->numCPUs(); ++i) { if(i==communicator_->myID()) {
  std::cerr << "t1 = " << t1/1000 << ", t2 = " << t2/1000 << ", t3 = " << t3/1000 << ", t4 = " << t4/1000 
            << ", t5 = " << t5/1000 << ", t6 = " << t6/1000 << std::endl;
  } communicator_->sync(); }*/
}

void
DistrNonnegativeMatrixFactorization::solveNNLS_MRHS(SCDoubleMatrix &A, SCDoubleMatrix &B, SCDoubleMatrix &X, int flag, int SSCflag)
{
  // if flag = 0 then solve: min ||AX^T-B|| s.t. X >= 0
  // if flag = 1 then solve: min ||AX^T-B^T|| s.t. X >= 0
  if(verboseFlag && communicator_->myID() == 0) {
    std::cerr << "\nm = " << A.getNumberOfRows() 
              << ", n = " << A.getNumberOfCols()
              << ", k = " << X.getNumberOfRows() << std::endl;
  }

  // number of sub-matrices in column-wise partition of B if flag == 0, or row-wise partition of B if flag == 1
  const int nsub = (nsub_ <= 0 || nsub_ > communicator_->numCPUs()) ? communicator_->numCPUs() : nsub_;

  // create new communicator and context if necessary
  MPI_Comm comm;
  int color, blacsHandle, context, rowCpus, colCpus;
  if(nsub > 1) {
    color = communicator_->myID()%nsub;
    MPI_Comm_split(*communicator_->getCommunicator(), color+1, 0, &comm);
    blacsHandle = Csys2blacs_handle(comm);
    context = blacsHandle;
    char order[] = "R";
    MPI_Comm_size(comm, &rowCpus);
    colCpus = 1;
    Cblacs_gridinit(&context, order, rowCpus, colCpus);
  }
  else {
    color = 0;
    comm = *communicator_->getCommunicator();
    context = A.getContext();
    rowCpus = A.getNumberOfProcsRow();
    colCpus = A.getNumberOfProcsCol();
  }

  Plh solver(A.getNumberOfRows(), A.getNumberOfCols());
  solver.setContext(context, rowCpus, colCpus, comm);
  solver.setQProcGrid(rowCpus, colCpus);
  solver.setABlockSize(blockSize_, blockSize_);
  solver.setQBlockSize(blockSize_, blockSize_);
  solver.setColumnScaling(); 
  solver.init();
  MPI_Barrier(comm); t1 -= getTime();
  for(int k = 0; k < nsub; ++k) {
    if(k != color) A.copyRedist(A.getNumberOfRows(), A.getNumberOfCols(), 1, 1, A.getContext());
    else solver.setMatrix(A);
  }
  t1 += getTime();
  solver.setVerbose(0);
  solver.setMaxIterRatio(3);
  if(domain->solInfo().solverTypeSpnnls == 5)
    solver.setOrthogonalMatchingPursuit();

  SCDoubleMatrix &b = solver.getRhsVector(); 
  SCDoubleMatrix &x = solver.getSolutionVector();

  std::vector<int> columns; 

  // copy/redistribute from B to subB, if necessary
  SCDoubleMatrix *subB, *subX; // this is for setting the correct constraints for sparse subspace clustering
  int nrhs;
  if(nsub == 1) { // if there is only 1 subdomain, B stays intact
    nrhs = X.getNumberOfRows();
    subX = &X;
    subB = &B;
    if(SSCflag) {
      for(int col = 1; col <= A.getNumberOfCols(); col++)
        columns.push_back(col); // column ordering doesn't change
    }
  } else {        // otherwise, B is divided up among the processes
    nrhs = X.getNumberOfRows()/nsub + ((color < X.getNumberOfRows()%nsub) ? 1 : 0);
    subX = new SCDoubleMatrix(context, nrhs, X.getNumberOfCols(), blockSize_, blockSize_, comm);
    MPI_Barrier(comm); t5 -= getTime();
    if(flag == 0) {
      subB   = new SCDoubleMatrix(context, B.getNumberOfRows(), nrhs, blockSize_, blockSize_, comm);
      for(int k=0,ja=1,n; k<nsub; ++k, ja+=n) {
        n = X.getNumberOfRows()/nsub + ((k < X.getNumberOfRows()%nsub) ? 1 : 0);
        if(k != color) B.copyRedist(A.getNumberOfRows(), n, 1, ja, A.getContext());              // send only
        else           B.copyRedist(A.getNumberOfRows(), n, 1, ja, *subB, 1, 1, A.getContext()); // send and receive
        if(SSCflag && k == color){
          for(int locCol = ja; locCol < ja+n; ++locCol)
            columns.push_back(locCol); // get local column numbers 
        }
      }
    }
    else {
      subB   = new SCDoubleMatrix(context, nrhs, B.getNumberOfCols(), blockSize_, blockSize_, comm);
      for(int k=0,ia=1,m; k<nsub; ++k, ia+=m) {
        m = X.getNumberOfRows()/nsub + ((k < X.getNumberOfRows()%nsub) ? 1 : 0);
        if(k != color) B.copyRedist(m, A.getNumberOfRows(), ia, 1, A.getContext());              // send only
        else           B.copyRedist(m, A.getNumberOfRows(), ia, 1, *subB, 1, 1, A.getContext()); // send and receive
      }
    }
    t5 += getTime();
  }

  int iter = 0;
  bool hotStart = domain->solInfo().hotstartSample;
  for(int i=1; i<=nrhs; i++) {
  
    t2 -= getTime();
    // copy ith column/row of matrix subB into vector b
    if(flag==0) subB->add(b, 'N', A.getNumberOfRows(), 1, 1.0, 0.0, 1, i, 1, 1);
    else        subB->add(b, 'T', A.getNumberOfRows(), 1, 1.0, 0.0, i, 1, 1, 1);
    t2 += getTime();

    if(hotStart && i > 1) solver.hotStart(); // use pre-computed solution
   
    // solve: min ||Ax-b|| s.t. x >= 0
    t3 -= getTime();
    if(SSCflag){ // don't let a snapshot select itself
      if(communicator_->myID() == 0) {
        fprintf(stderr,"\r");
        fprintf(stderr,"%3.2f percent done",100.0*double(i)/double(nrhs));
      }
      solver.setConstraint(columns[i-1]); // dont let clustering solver select itself
      solver.setRtol(tol_);
    } else {
      solver.setRtol(1e-16);
    }
    solver.solve();
    iter += solver.getIter();
    t3 += getTime();
    // copy x to ith row of subX
    t4 -= getTime();
    x.add(*subX, 'N', 1, A.getNumberOfCols(), 1.0, 0.0, 1, 1, i, 1);
    t4 += getTime();
  }

  if(SSCflag && communicator_->myID() == 0) fprintf(stderr,"\n");

  if(verboseFlag && communicator_->myID() == 0) {
    std::cerr << "Total number of iterations = " << iter << std::endl;
  }
  //solver.printTimes(true);

  // copy/redistribute from subX to X, if necessary
  if(nsub > 1) {
    MPI_Barrier(comm); t6 -= getTime();
    for(int k=0,ib=1,m; k<nsub; ++k,ib+=m) {
      m = X.getNumberOfRows()/nsub + ((k < X.getNumberOfRows()%nsub) ? 1 : 0);
      if(k != color) SCDoubleMatrix::copyRedist(m, A.getNumberOfCols(), X, ib, 1, A.getContext()); // receive only
      else           subX->copyRedist(m, A.getNumberOfCols(), 1, 1, X, ib, 1, A.getContext());     // send and receive
    }
    t6 += getTime();
    delete subB;
    delete subX;
    Cblacs_gridexit(context);
    Cfree_blacs_system_handle(blacsHandle);
  }
}

} // end namespace Rom
#endif
