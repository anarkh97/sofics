#include <Math.d/Vector.h>

template < class ConditionOps, class PostProcessor, class ProblemDescriptor >
void
CondSolver < ConditionOps, PostProcessor, ProblemDescriptor >::solve()
{
 using std::abs;
 int i,j;
 fprintf(stderr, " ... Computing the Condition Number ...\n");

 // ... Build Connectivities, dof set, etc.
 probDesc->preProcess();

 // ... Get the Condition Number post processor
 PostProcessor *postProcessor = probDesc->getPostProcessor();

 // ... Get the necessary matrices
 ConditionOps dMat = probDesc->buildCondOps();

 // ... Get the required input data
 double relTol;
 int ndof, maxIte;
 probDesc->getConditionInfo(relTol,maxIte,ndof);      

 double eigmax, preeigmax = 0.0; 
 double eigmin, preeigmin = 0.0;
 //int maxIte = ndof * 50;

 Vector v(probDesc->solVecInfo());
 Vector z(probDesc->solVecInfo());

// Power method.
// Starts from an arbitray array.

 for (i=0; i< ndof; i++)
    v[i] = (double) (i+1) / (double) ndof;

// Power method iteration loop

 for (i=0; i<maxIte; i++) {
    dMat.K->mult(v,z);

/*
 for (j=0; j< ndof; j++)
    z[j] /= dMat.M->diag(j);
*/

// Normalize

  // compute absolute maximum value of z
  double zmax = z.absMax();

  eigmax = zmax;

  // scalar multiplication: v = (1/zmax)*z
  v = (1.0/zmax)*z;

  if ( abs(eigmax - preeigmax) < relTol*abs(preeigmax) ) break;

  preeigmax = eigmax;
  }

  // Inverse power method.
  // Starts from an arbitray array.

  for (i=0; i< ndof; i++)
     v[i] = (double) (i+1) / (double) ndof;

// Inverse power method iteration loop

  for (i=0; i<maxIte; i++) {

     for (j=0; j< ndof; j++)
        z[j] = v[j]; // * dMat.M->diag(j);

     dMat.dynMat->reSolve(z);

// Normalize

     double zmax = z[0];
     for (j=1; j<ndof; j++)
        if (abs(z[j])>zmax) zmax = abs(z[j]);

     eigmin = 1.0 / zmax;

     // scalar multiplication: v = (1/zmax)*z
     v = (1.0/zmax)*z;

     if ( abs(eigmin - preeigmin) < relTol*abs(preeigmin) ) break;

     preeigmin = eigmin;
  }

  postProcessor->conditionOutput(eigmin,eigmax);

}
