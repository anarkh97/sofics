#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <Utils.d/dbg_alloca.h>

#include <Driver.d/Domain.h>
#include <Driver.d/Dynam.h>
#include <Problems.d/EigenDescr.h>

#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Math.d/NBSparseMatrix.h>
#include <Math.d/EiSparseMatrix.h>
#include <Math.d/Vector.h>
#include <Math.d/VectorSet.h>
#include <Math.d/AddedMassMatrix.h>

#include <Solvers.d/Solver.h>
#include <Solvers.d/Rbm.h>
#include <Timers.d/StaticTimers.h>
#include <Timers.d/GetTime.h>
#include <Corotational.d/GeomState.h>
#include <Utils.d/DistHelper.h>
#include <Utils.d/ModeData.h>

extern ModeData modeData;

int
SingleDomainEigen::solVecSize()
{
 // ... returns number of unconstrained dof
 return domain->numUncon();
}

int
SingleDomainEigen::solVecInfo()
{
 // ... returns number of unconstrained dof
 return domain->numUncon();
}

void
SingleDomainEigen::preProcess()
{

 // Allocate space for the Static Timers
 times = new StaticTimers;

 times->preProcess -= getTime();
 // Makes renumbering, connectivities and dofsets
 domain->preProcessing();

 // Total number of dof
 int numdof  = domain->numdof();

 times->makeBCs -= getTime();
 int *bc  = new int[numdof];
 bcx = new double[numdof];

 // ... make boundary conditions
 domain->make_bc(bc,bcx);
 times->makeBCs += getTime();

 // ... make constrained dof set array
 times->makeDOFs -= getTime();
 domain->make_constrainedDSA();

 // ... construct all dofs
 domain->makeAllDOFs();
 times->makeDOFs += getTime();

 if(domain->solInfo().HEV) {
   domain->makeAllDOFsFluid();
 }
 
 // if we have initial displacements, we have to consider
 // the nonlinear tangent stiffness matrix instead of the
 // linear stiffness matrix. This is due to the prestress.

 geomKelArray = 0;
 if( domain->numInitDisp6() > 0 && domain->solInfo().gepsFlg == 1) {

   if(domain->solInfo().buckling) {
     fprintf(stderr," ... Buckling Analysis              ...\n");
     times->kelArrayTime -= getTime();
     domain->createKelArray(geomKelArray);
     times->kelArrayTime += getTime();
   }

   FullSquareMatrix *dummy=0;
   domain->computeGeometricPreStress(allCorot, geomState, kelArray, times,
                                     geomKelArray, dummy);

   if(domain->solInfo().buckling == 1) {
     //delete [] kelArray;
     kelArray = 0; 
   }
 }

 times->preProcess += getTime();
}

SDEigenPostProcessor *
SingleDomainEigen::getPostProcessor()
{
 return new SDEigenPostProcessor(domain, times, bcx);
}

void
SingleDomainEigen::buildEigOps( DynamMat &dMat )
{
 AllOps<double> allOps;

 if(domain->solInfo().addedMass == 2) allOps.M =
                                          new AddedMassMatrix<double, Domain>(domain->getNodeToNode(), domain->getDSA(),
                                                                              domain->getCDSA(),
                                                                              domain,
                                                                              &Domain::multCV, &Domain::trMultCV);
#ifdef USE_EIGEN3
 else if(domain->solInfo().printMatLab)
   allOps.M = domain->constructEiSparse<double>();
#endif
 else
   allOps.M = domain->constructDBSparseMatrix<double>();
#ifdef USE_EIGEN3
 if(domain->solInfo().printMatLab)
   allOps.K = domain->constructEiSparse<double>();
 else
#endif
 allOps.K = domain->constructDBSparseMatrix<double>();

 // ... Construct geometric rigid body modes if necessary
 if(domain->solInfo().rbmflg) { 
   dMat.rigidBodyModes = domain->constructRbm();
 }

 // ... Construct "thermal rigid body mode" if necessary
 else if(domain->solInfo().hzemFlag) {
   dMat.rigidBodyModes = domain->constructHzem();
 }

 // ... Construct "sloshing rigid body mode" if necessary
 else if(domain->solInfo().slzemFlag) { 
   dMat.rigidBodyModes = domain->constructSlzem();
 }

 // build stiffness and mass matrices
 melArray = (domain->solInfo().arpack_mode == 4 && domain->solInfo().buckling) ? geomKelArray : 0;
 domain->buildOps<double>(allOps, 1.0, 0.0, 0.0, dMat.rigidBodyModes, kelArray, melArray);
 dMat.dynMat  = allOps.sysSolver;
 dMat.M       = allOps.M;
#ifdef USE_EIGEN3
 if(domain->solInfo().printMatLab) {
   allOps.M->printSparse(std::string(domain->solInfo().printMatLabFile)+".mass");
   allOps.K->printSparse(std::string(domain->solInfo().printMatLabFile)+".stiffness");
   if(domain->solInfo().printMatLabExit) return;
 }
#endif

 if(domain->solInfo().addedMass == 2) ((AddedMassMatrix<double, Domain>*)allOps.M)->setFluidSolver(domain->Mff);

 if(domain->solInfo().explicitK)
   dMat.refK    = allOps.K;

 // If we are doing a pre-stress computation.
 if(geomKelArray && domain->solInfo().gepsFlg == 1) {
   fprintf(stderr," ... Assembling Geometric Stiffness ...\n");
   Connectivity *allDOFs = domain->getAllDOFs();
   dMat.M->zeroAll();
   int iele;
   int size = sizeof(double)*domain->maxNumDOF()*domain->maxNumDOF();
   double *karray = (double *) dbg_alloca(size);
   for(iele=0; iele<domain->numElements(); ++iele) {
     if(domain->solInfo().arpack_mode == 4) {
       // this is an unnecessary recalculation of Element::stiffness but it will do for now
       FullSquareMatrix kel  = domain->getElementSet()[iele]->stiffness(domain->getNodes(), karray);
       dMat.M->add(kel,(*domain->getAllDOFs())[iele].data());
     }
     else
       dMat.M->add(geomKelArray[iele],(*allDOFs)[iele].data());
   }
 }
 else if(domain->solInfo().arpack_mode == 4) { // we may also want to use mode for for eigen (non-buckling)
                                               // if the mass matrix is indefinite
   Connectivity *allDOFs = domain->getAllDOFs();
   dMat.M->zeroAll();
   int iele;
   int size = sizeof(double)*domain->maxNumDOF()*domain->maxNumDOF();
   double *karray = (double *) dbg_alloca(size);
   for(iele=0; iele<domain->numElements(); ++iele) {
     // this is an unnecessary recalculation of Element::stiffness but it will do for now
     FullSquareMatrix kel  = domain->getElementSet()[iele]->stiffness(domain->getNodes(), karray);
     dMat.M->add(kel,(*domain->getAllDOFs())[iele].data());
   }
 }
}

void
SingleDomainEigen::reBuild( DynamMat &dMat )
{
 AllOps<double> allOps;
 dMat.M->zeroAll();
 allOps.M    = dMat.M ;
 allOps.sysSolver = dMat.dynMat;

 // rebuild stiffness and mass matrices
 // watch: no rigid body modes assumed

 domain->rebuildOps<double>(allOps, 1.0, 0.0, 0.0, (Rbm *) NULL, kelArray, melArray);
}

int SingleDomainEigen::getNumEigen()
{
  return domain->solInfo().nEig;
}

int SingleDomainEigen::getEigenSolverType()
{
  return domain->solInfo().eigenSolverType;
}

int SingleDomainEigen::getEigenSolverSubType()
{
  return domain->solInfo().eigenSolverSubType;
}

void
SingleDomainEigen::getSubSpaceInfo(int& subspacesize, int& maxIter, 
                                   double& tolEig, double& tolJac, bool &explicitK )
{
 subspacesize = domain->solInfo().subspaceSize;
 tolEig       = domain->solInfo().tolEig;
 tolJac       = domain->solInfo().tolJac;

 int numIter1 = domain->solInfo().maxitEig;
 int numIter2 = 10*subspacesize;

 maxIter = std::max(numIter1, numIter2);

 if (numIter1 > 0) maxIter = numIter1;

 explicitK = domain->solInfo().explicitK;
}

void
SingleDomainEigen::error(int subspacesize, int numRbm)
{
 fprintf(stderr,"*** subspace size must be bigger than number of rbm   ***\n");
 fprintf(stderr,"subspace size = %d\n",subspacesize);
 fprintf(stderr,"number of rbm = %d\n",numRbm);
}

#ifdef WINDOWS
#define drand48 rand
#endif

void
SingleDomainEigen::initQ(Vector* Q, int size)
{
  int numdof = domain->numUncon();
  int i, j;
  if(size == 0) {
    fprintf(stderr," *** ERROR: Division by zero within initQ routine. \n");
    return;
  }

  Vector t(numdof);

  // This has been changed to uses random vectors for the initial subspace.
  // The old way (commented out) could lead to a linear dependence between
  // RBMs and initial vectors.
  // WARNING: does drand48() exist on all machines?
  for(i=0; i<size; ++i) {
    t.zero();
    for(j=0; j<numdof; ++j)
      t[j] = drand48();
    Q[i] = t;
  }
}

void
SingleDomainEigen::printTimers(Solver *solver)
{
 long memoryUsed = solver->size();
 double solveTime  = solver->getSolutionTime();

 times->printStaticTimers(solveTime, memoryUsed, domain);

  if(domain->solInfo().doEigSweep && domain->solInfo().massFlag) {
    double mass = domain->computeStructureMass();
    filePrint(stderr," --------------------------------------\n");
    filePrint(stderr," ... Structure mass = %e  ...\n",mass);
    filePrint(stderr," --------------------------------------\n");
  }
}

void
SingleDomainEigen::convertModeDataToVecSet(VectorSet& vModeData)
{
  DofSetArray *cdsa = domain->getCDSA();
  for(int j=0; j<domain->numNodes(); ++j) { // loop over nodes
    for(int k=0; k<6; ++k) { // loop over dofs
      int dof = cdsa->locate(modeData.nodes[j], 1 << k);
      if(dof >= 0)
        for(int i=0; i<modeData.numModes; ++i) {
          vModeData[i][dof] = modeData.modes[i][j][k];
        }
    }
  }
}

void
SDEigenPostProcessor::eigenOutput(Vector& eigenValues, VectorSet& eigenVectors, int convEig)
{
 times->output -= getTime(); 
 domain->eigenOutput(eigenValues,eigenVectors,bcx,convEig);
 times->output += getTime(); 
}

#ifdef USE_EIGEN3
void
SDEigenPostProcessor::eigenQROutput(Eigen::MatrixXd& XVectors, Eigen::MatrixXd& QVectors, Eigen::MatrixXd& RVectors)
{
 times->output -= getTime();
 domain->eigenQROutput(XVectors, QVectors,RVectors);
 times->output += getTime();
}
#endif
