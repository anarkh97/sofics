#include <Utils.d/dbg_alloca.h>
#include <cstdio>

#include <Corotational.d/Corotator.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/GeomNLSolver.h>
#include <Math.d/CuCSparse.h>
#include <Solvers.d/Solver.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>

GeomNLSolver::GeomNLSolver(int _numele, Connectivity* _allDOFs,
                           int _maxNumDof, CoordSet &_nodes,int numdofs, 
                           SolverInfo &_sinfo, FullSquareMatrix *_kelArray,
                           DofSetArray *_c_dsa, Vector &f,
 		           Corotator **_allCorot, Solver *_solver)

 : sinfo(_sinfo), force(f), residual(numdofs), pVec(numdofs), 
   arcLenResid(numdofs), elementInternalForce(_maxNumDof), nodes(_nodes)
{
 numele    = _numele;		// number of elements
 allCorot  = _allCorot;		// all Element Corotators
 kelArray  = _kelArray;		// element stiffness matrix array
 allDOFs   = _allDOFs;		// element dof array
 maxNumDof = _maxNumDof;	// maximum number of dof per element
 c_dsa     = _c_dsa;		// constrained dof set array
 solver    = _solver;		// specified solver to be used


 // Need a better formula for this # of iterations
 // or it should be input from the file.
 // i.e. maximum number of iterations for newton
 //      maximum number of iterations for extended newton
 //      should be separate

 // maximum # of iterations for extended newton
 // maxit = 1 / sinfo.getNLInfo().dlambda;

 // Compute norm of external force vector
 forceNorm = force.squareNorm();
 //HB 
 fprintf(stderr,"--- forceNorm = %e\n",forceNorm);

 double lfactor = sinfo.getNLInfo().lfactor;
 lmin = sinfo.getNLInfo().dlambda / lfactor;
 lmax = sinfo.getNLInfo().dlambda * lfactor;
 
}

GeomNLSolver::GeomNLSolver(int _numele, Domain &domain, CoordSet &_nodes, 
			   int numdofs, SolverInfo &_sinfo,
                           DofSetArray *_c_dsa, Vector &f, 
			   Corotator **_allCorot) 
 : sinfo(_sinfo), force(f), residual(numdofs), pVec(numdofs), 
   arcLenResid(numdofs), elementInternalForce(domain.maxNumDOF()), nodes(_nodes)
{
 numele   = _numele;
 allCorot = _allCorot;

 domain.createKelArray(kelArray);

 allDOFs   = domain.getAllDOFs();
 maxNumDof = domain.maxNumDOF();

 c_dsa = _c_dsa;

 // Maximum # of iterations for extended newton
 //maxit = 1 / sinfo.getNLInfo().dlambda;

 // Compute first norm based on input force
 forceNorm = force.squareNorm();
 //HB
 fprintf(stderr,"--- forceNorm = %e\n",forceNorm);
  
}

void
GeomNLSolver::reBuild(int iteration)
{
   if( iteration % sinfo.getNLInfo().updateK == 0 ) solver->reBuild(kelArray);
}

void
GeomNLSolver::getStiffAndForce(GeomState &u, Vector &residual)
{
/*******************************************************************
 *
 * Purpose :
 *  Compute element tangential stiffness and element internal force
 *  and to assemble them.
 *
 * Input :
 *  u        : current node geometric state
 *  nodes    : undeformed coordinate set of nodes
 *
 * Output :
 *  residual : residual vector = force - p 
 *  kelArray : array of element tangential stiffness matrices 
 *             in current configuration
 *
 *****************************************************************/
   int dofNum, idof;

   int iele;
   for(iele = 0; iele < numele; ++iele) {

      // Get updated stiffness and element internal force vector
      allCorot[iele]->getStiffAndForce(u, nodes, kelArray[iele], 
                      elementInternalForce.data(), 0, 0);

      // Symmetrize element stiffness matrix
      if(!domain->solInfo().getNLInfo().unsymmetric) kelArray[iele].symmetrize();

      // Assemble element internal force into residual vector
      // r(lambda,v) = force(lambda) - p(v)

      for(idof = 0; idof < kelArray[iele].dim(); ++idof) {
         dofNum = c_dsa->getRCN((*allDOFs)[iele][idof]);
         if(dofNum >= 0)
            residual[dofNum] -= elementInternalForce[idof];
      }
   }

}

int
GeomNLSolver::checkConvergence(int iter, double residualNorm, double firstNorm)
{
 double relativeNorm = residualNorm / firstNorm;
 fprintf(stderr,"   #%d relative norm = %15.6e norm = %15.6e\n",
                iter,relativeNorm,residualNorm);
 return (relativeNorm <= sinfo.getNLInfo().tolRes) ? 1 : 0;
}

void
GeomNLSolver::newton(GeomState &u, const double lambda)
{

// ***************************************************************
// *                                                             *
// *  Purpose:  To solve r = f - p = 0 using Newton's method     *
// *                                                             *
// *    Input:                                                   *
// *          u = current node geometric state                   *
// *     lambda = current value of control parameter             *
// *                                                             *
// *   Output:                                                   *
// *          u = updated node geometric state   		 *
// *                                                             *
// ***************************************************************


 fprintf(stderr,"--- Newton Method with lambda = %f\n",lambda);

 double firstNorm = (lambda*lambda)*forceNorm;

 //HB
 fprintf(stderr,"--- firstNorm = %e\n",firstNorm);

 if( firstNorm == 0.0 ) return;

 // MAIN ITERATION LOOP
 int maximumIterations = sinfo.getNLInfo().maxiter;
 int iter;
 for(iter = 0; iter < maximumIterations; ++iter) {

   // Set residual = lambda*force
   residual.linC( lambda, force);

   // Compute tangent stiffness matrix and residual force vector
   getStiffAndForce(u, residual); 

   // Compute residual norm
   double residualNorm = residual.squareNorm();

   // rebuild tangent stiffness matrix 
   reBuild(iter);

   // Solve system Kv = residual, Overwriting residual with solution
   solver->reSolve(residual);

   // Update geometric node state
   u.update(residual);

   int converged = checkConvergence(iter, residualNorm, firstNorm);
   
   if (converged) break;
 }

}

void
GeomNLSolver::extendedNewton(GeomState &u, Vector &dU, double &lambda, 
         double deltaLambda, double &deltaS, double w, int &numExtIter)
{
  // Check this formula
  //double dlambda = (deltaS-dU*dU)/(w*deltaS);
  //if(dlambda < lmin) dlambda = lmin;
  //if(dlambda > lmax) dlambda = lmax;
  //lambda += dlambda;

  // This alternative constant formula allows me to apply a larger load.
  lambda += deltaLambda;

  // Compute first norm
  double firstNorm = (lambda*lambda)*forceNorm;
  //HB
  fprintf(stderr,"--- firstNorm = %e\n",firstNorm);

  if(firstNorm == 0.0) return;

  for(numExtIter = 0; numExtIter < 50; ++numExtIter) {

    // Compute residual = lambda*force
    residual.linC(lambda, force);

    getStiffAndForce(u, residual);

    reBuild(numExtIter);

    arcLenResid = residual;

    solver->reSolve(residual);

    solver->solve(force, pVec);

    double mu = -1.0*( dU * residual ) / (w*deltaLambda + dU * pVec);

    // Compute: pVec = mu*pVec + residual
    pVec.linC(mu, pVec, 1.0, residual);

    // Update geometrical state
    u.update(pVec);

    // Update lambda
    lambda += mu;

    // Compute arcLenResid = arcLenResid + mu*force
    //HB
    double residualNorm = arcLenResid.squareNorm();
    arcLenResid.linC(arcLenResid, mu, force);

    //double residualNorm = arcLenResid.squareNorm();
    // KHP:
    // double pnorm = pVec.squareNorm();
    int converged = checkConvergence(numExtIter, residualNorm, firstNorm);


    if (converged) break;
  }
}

void
Domain::arcLength()
{
 // ... CONSTRUCT CONNECTIVITIES, DOFSETARRAY AND RENUMBERING
 preProcessing();

 // ... SET UP BOUNDARY CONDITIONS
 int numdofs = numdof();

 int *bc = (int *) dbg_alloca(sizeof(int)*numdofs);
 double *bcx = new double[numdofs];

 make_bc(bc, bcx);

 // Currently, boundary condition values are not being used
 // i.e. nonhomogeneous dirichlet boundary conditions

 make_constrainedDSA(bc);

 // NOW SET NUMDOFS EQUAL TO THE NUMBER OF UNCONSTRAINED DOF
 numdofs = numUncon();
 
 // MAKE ALL OF THE ELEMENTS DOF ARRAYS
 makeAllDOFs();

 // ALLOCATE MEMORY FOR EXTERNAL LOAD VECTOR
 Vector force(numdofs,0.0);
 Vector dummyAeroForce(numdofs,0.0);

 // ... ALLOCATE MEMORY FOR ARRAY OF POINTERS TO COROTATORS
 Corotator **allCorot = new Corotator *[numele];

 // ... CREATE THE ARRAY OF POINTERS TO COROTATORS
 createCorotators(allCorot);

 // ... CREATE THE ARRAY OF ELEMENT STIFFNESS MATRICES
 FullSquareMatrix *kelArray;
 createKelArray(kelArray);

 // ... BUILD THE SOLVER
 AllOps<double> allOps;
 allOps.Kuc = constructCuCSparse<double>();

 double Kcoef = 1.0;
 double Mcoef = 0.0;
 double Ccoef = 0.0;

 buildOps<double>(allOps, Kcoef, Mcoef, Ccoef, (Rbm *) NULL, (FullSquareMatrix *) NULL,
                  (FullSquareMatrix *) NULL, (FullSquareMatrix *) NULL, true);

 solver = allOps.sysSolver;

 // ... BUILDRHSFORCE SUMS GRAVITY FORCE & USER-DEFINED FORCE.
 buildRHSForce(force, allOps.Kuc);

 // Add Kuc to GeomNLSolver class:
 GeomNLSolver nlSolver(numele,  allDOFs, maxNumDOFs, nodes, numdofs, 
		        sinfo, kelArray,      c_dsa, force, allCorot, solver );

 // ... INITIALIZE THE REFERENCE GEOMETRY STATE
 GeomState u0(*dsa, *c_dsa, nodes);

 // ... INITIALIZE GeomState u = u0
 GeomState u(u0);

 // ... Output initial configuration
 postProcessing(&u0, force, dummyAeroForce, 0.0, 1, 0, 0, allCorot);

 // ... DEFINE deltaLambda0
 double deltaLambda0 = sinfo.getNLInfo().dlambda;

 // ... DEFINE minimum delta Lambda and maximum delta Lambda
 // lmin <= deltaLamba <= lmax
 double lfactor = sinfo.getNLInfo().lfactor;
 double lmin = deltaLambda0 / lfactor;
 double lmax = deltaLambda0 * lfactor;

 // ... DEFINE initial lambda0
 double lambda0 = 0.0;

 // ... COMPUTE FIRST LOAD CONTROL PARAMETER
 double lambda = lambda0 + deltaLambda0;

 // ... CALL NEWTON FOR FIRST SOLUTION
 nlSolver.newton( u, lambda );

 // ... Declare Vector dU
 Vector dU(numdofs);

 // ... Define initial dU = u - u0
 u.diff(u0, dU);

 // ... Define deltaS	(arc length distance)
 double deltaS = dU.squareNorm();

 // ... COMPUTE W = SMOOTHING PARAMETER
 double w = (deltaS*deltaS) / (deltaLambda0*deltaLambda0);

 // ... DEFINE NU = load control parameter multiplier 
 double nu = 1.0;

 // ... DEFINE MAXIMUM NUMBER OF TRAJECTORY ITERATIONS
 int maxit = 2000;
 // HB
 maxit = 1;

// ---------------------------------------------------

 double deltaLambda;
 int numExtIter;

 // ... COMPUTE TRAJECTORY (EQUILIBRIUM) PATH 
 int iter;
 for(iter=0; iter<maxit; ++iter) {

   postProcessing(&u, force, dummyAeroForce, lambda, 1, 0, 0, allCorot);

   // Compute: dU = nu*(u - u0);
   u.diff(u0, dU);
   dU *= nu;

   // ... Set Geometry State u0 = u
   u0 = u;

   // ... COMPUTE NEXT DELTALAMBDA
   deltaLambda = nu*(lambda - lambda0); 

   // ... CHECK MAGNITUDE OF DELTALAMBDA
   if(deltaLambda > lmax ) deltaLambda = lmax;
   if(deltaLambda < lmin ) deltaLambda = lmin;

   // ... SET lambda0 = lambda
   lambda0 = lambda;

   // ... COMPUTE DELTAS
   deltaS *= nu;

   fprintf(stderr,"Equilibrium #%d deltaLambda = %10.6f lambda = %10.6f\n",
                  iter+1, deltaLambda, lambda);

   // ... CALL EXTENDED NEWTON
   nlSolver.extendedNewton(u,dU,lambda,deltaLambda,deltaS,w,numExtIter);

   if( std::abs(lambda) >= 1.0 ) break;

   // ... DETERMINE CONTROL PARAMETER BASED ON # OF ITERATIONS 
   // ... IN EXTENDED NEWTON 
   nu = 1.0;
   //if(numExtIter < 4) nu = 2.0;
   //if(numExtIter > 6) nu = 0.5;
   //fprintf(stderr,"    nu = %e\n", nu);
   //if(iter>20) nu = 0.5; 
 }
 fprintf(stderr,"Equilibrium #%d deltaLambda = %10.6f lambda = %10.6f\n",
                iter+1, deltaLambda, lambda);

 // ... DEFINE INITIAL GUESS
 u.interp((1.0 - lambda0)/(lambda - lambda0), u, u0 );

 // ... CALL NEWTON FOR FINAL SOLUTION
 nlSolver.newton(u);

 // CALL POST PROCESSING OF DISPLACEMENTS
 postProcessing(&u, force, dummyAeroForce, 1.0, 1, 0, 0, allCorot);

 // deallocate memory
 delete [] allCorot;

}
