#include <Corotational.d/Corotator.h>
#include <Corotational.d/GeomState.h>
#include <Math.d/FullSquareMatrix.h>
#include <iostream>

void 
Corotator::formGeometricStiffness(GeomState& gs, CoordSet& cs, FullSquareMatrix& Kg, double* f) 
{ 
  // default implementation
  GeomState gs0(cs); // TODO 
  double t = 0, dt = 0;
  int n = Kg.dim();

  FullSquareMatrix K(n);
  getStiffAndForce(gs0, cs, K, f, dt, t);

  FullSquareMatrix Kt(n);
  getStiffAndForce(gs, cs, Kt, f, dt, t);  

  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j)
      Kg[i][j] = Kt[i][j] - K[i][j];
}

double* 
Corotator::getOriginalStiffness() 
{ 
  std::cerr << " *** WARNING: Corotator::getOriginalStiffness() is not implemented\n";
  return 0; 
}

void
Corotator::extractDeformations(GeomState&, CoordSet&, double*, int&) 
{
  std::cerr << " *** WARNING: Corotator::extractDeformations(GeomState&, CoordSet&, double*, int&) is not implemented\n"; 
}

void
Corotator::extractDeformations(GeomState&, CoordSet&, DComplex*, int&)
{
  std::cerr << " *** WARNING: Corotator::extractDeformations(GeomState&, CoordSet&, DComplex*, int&) is not implemented\n";
}

void
Corotator::getNLVonMises(Vector&, Vector&, GeomState&, GeomState*, CoordSet&, int, int, double, double, int, int)
{
  std::cerr << " *** WARNING: Corotator::getNLVonMises(Vector&, Vector&, GeomState&, GeomState*, CoordSet&, int, ...) is not implemented\n";
}

void
Corotator::getNLVonMises(ComplexVector&, Vector&, GeomState&, GeomState*, CoordSet&, int, int, double, double, int, int)
{
  std::cerr << " *** WARNING: Corotator::getNLVonMises(ComplexVector&, Vector&, GeomState&, GeomState*, CoordSet&, int, ...) is not implemented\n";
}

void
Corotator::getNLAllStress(FullM&, Vector&, GeomState&, GeomState*, CoordSet&, int, int, int)
{
  std::cerr << " *** WARNING: Corotator::getNLAllStress(FullM&, Vector&, GeomState&, GeomState*, CoordSet&, int, ...) is not implemented\n";
}

void
Corotator::getNLAllStress(FullMC&, Vector&, GeomState&, GeomState*, CoordSet&, int, int, int)
{
  std::cerr << " *** WARNING: Corotator::getNLAllStress(FullMC&, Vector&, GeomState&, GeomState*, CoordSet&, int, ...) is not implemented\n";
}

double
Corotator::getElementEnergy(GeomState&, CoordSet&)
{
  std::cerr << " *** WARNING: Corotator::getElementEnergy(GeomState&, CoordSet&) is not implemented\n";
  return 0.0;
}

void
Corotator::extractRigidBodyMotion(GeomState&, CoordSet&, double*)
{
  std::cerr << " *** WARNING: Corotator::extractRigidBodyMotion(GeomState&, CoordSet&, double*) is not implemented\n";
}

void
Corotator::getInternalForceThicknessSensitivity(GeomState*, GeomState&, CoordSet&, Vector&, double, double)
{
  std::cerr << " *** WARNING: Corotator::getInternalForceThicknessSensitivity() is not implemented\n";
}

void
Corotator::getInternalForceNodalCoordinateSensitivity(GeomState*, GeomState&, CoordSet&, Vector*&, double, double)
{
  std::cerr << " *** WARNING: Corotator::getInternalForceNodalCoordinateSensitivity() is not implemented\n";
}

void
Corotator::extractDeformationsDisplacementSensitivity(GeomState&, CoordSet&, double*)
{
  std::cerr << " *** WARNING: Corotator::extractDeformationsDisplacementSensitivity() is not implemented\n";
}

void
Corotator::getNLVonMisesNodalCoordinateSensitivity(Vector&, Vector&, GeomState&, GeomState*, CoordSet&, int, int, double, double, int, int)
{
  std::cerr << " *** WARNING: Corotator::getNLVonMisesNodalCoordinateSensitivity() is not implemented\n";
}

void
Corotator::getNLVonMisesThicknessSensitivity(Vector&, Vector&, GeomState&, GeomState*, CoordSet&, int, int, double, double, int, int)
{
  std::cerr << " *** WARNING: Corotator::getNLVonMisesThicknessSensitivity() is not implemented\n";
}

void
Corotator::getNLVonMisesDisplacementSensitivity(Vector&, Vector&, GeomState&, GeomState*, CoordSet&, int, int, double, double, int, int)
{
  std::cerr << " *** WARNING: Corotator::getNLVonMisesDisplacementSensitivity() is not implemented\n";
}
