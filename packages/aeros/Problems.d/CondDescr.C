#include <cstdio>
#include <Utils.d/dbg_alloca.h>

#include <Problems.d/CondDescr.h>
#include <Driver.d/Domain.h>
#include <Driver.d/Dynam.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>

int
SingleDomainCond::solVecInfo()
{
 // ... returns number of unconstrained dof
 return domain->numUncon();
}

void
SingleDomainCond::preProcess()
{
 // ... Makes renumbering, connectivities and dofsets
 domain->preProcessing();

 // ... Total number of dof
 int numdof = domain->numdof();

 // ... allocate arrays for boundary conditions
 int    *bc  = (int *)    dbg_alloca(sizeof(int)*numdof);
 double *bcx = (double *) dbg_alloca(sizeof(double)*numdof);

 // ... make boundary conditions
 domain->make_bc(bc, bcx);

 // ... make constrained dof set array
 domain->make_constrainedDSA();

 // ... construct all dofs
 domain->makeAllDOFs();

}

SDCondPostProcessor *
SingleDomainCond::getPostProcessor()
{
 return new SDCondPostProcessor(domain);
}

DynamMat
SingleDomainCond::buildCondOps()
{
 AllOps<double> allOps;

 // ... Construct a sparse matrix K
 allOps.K = domain->constructDBSparseMatrix<double>();

 // ... Build K and a Solver
 domain->buildOps<double>(allOps, 1.0, 0.0, 0.0);

 DynamMat dMat;

 dMat.dynMat  = allOps.sysSolver;
 dMat.K       = allOps.K;

 return dMat;
}

void
SingleDomainCond::getConditionInfo(double& tolerance, int& maxit, int& numdof)
{
 tolerance = domain->solInfo().condNumTolerance;
 maxit     = domain->solInfo().condNumMaxit;
 numdof    = domain->numUncon();
}

void
SDCondPostProcessor::conditionOutput(double eigmin, double eigmax)
{
 fprintf(stderr, " ... Lowest  Eigen Value = %9.5e ...\n",eigmin);
 fprintf(stderr, " ... Highest Eigen Value = %9.5e ...\n",eigmax);

 // ... PRINT THE CONDITION NUMBER
 fprintf(stderr, " ... Condition Number    = %9.5e ...\n",eigmax/eigmin);
 fflush(stderr);
}
