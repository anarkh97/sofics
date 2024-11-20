#include <Corotational.d/PhantomCorotator.h>
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>

void
PhantomCorotator::getStiffAndForce(GeomState &, CoordSet &, FullSquareMatrix &k, double *, double, double)
{
  k.zero();
}

void 
PhantomCorotator::formGeometricStiffness(GeomState&, CoordSet&, FullSquareMatrix &k, double*) 
{ 
  k.zero();
}

double* 
PhantomCorotator::getOriginalStiffness() 
{ 
  return 0; 
}

void
PhantomCorotator::extractDeformations(GeomState&, CoordSet&, double*, int &nlflag) 
{
  nlflag = 1;
}

void
PhantomCorotator::extractDeformations(GeomState&, CoordSet&, DComplex*, int &nlflag)
{
  nlflag = 1;
}

double
PhantomCorotator::getElementEnergy(GeomState&, CoordSet&)
{
  return 0;
}

void
PhantomCorotator::extractRigidBodyMotion(GeomState&, CoordSet&, double*)
{
}

