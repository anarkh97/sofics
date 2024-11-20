// ----------------------------------------------------------------
// HB - 08/13/03
// ----------------------------------------------------------------
// Standard C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>

// FEM headers
#include <Element.d/Element.h>
#include <Math.d/matrix.h>

#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/MortarElement.d/MortarTri6.d/DualMortarTri6.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
DualMortarTri6::DualMortarTri6() 
{
  Initialize();
} 

DualMortarTri6::DualMortarTri6(FaceElement* FaceElem)
{
  Initialize();
  SetPtrMasterFace(FaceElem);
}

DualMortarTri6::DualMortarTri6(double _area, FaceElement* FaceElem)
{
  Initialize();
  SetArea(_area);
  SetPtrMasterFace(FaceElem);
}

DualMortarTri6::DualMortarTri6(FaceElement* FaceElem, CoordSet &cs)
{
  Initialize();
  SetPtrMasterFace(FaceElem);
  ComputeDualCoeffs(cs);
}

DualMortarTri6::DualMortarTri6(double _area, FaceElement* FaceElem, CoordSet &cs)
{
  Initialize();
  SetArea(_area);
  SetPtrMasterFace(FaceElem);
  ComputeDualCoeffs(cs);
}

// -----------------------------------------------------------------------------------------------------
//                                            DESTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
DualMortarTri6::~DualMortarTri6()
{
  if(Alpha) { delete Alpha; Alpha = 0; }
}

// -----------------------------------------------------------------------------------------------------
//                                   INITILIZATION & CLEAN/CLEAR METHODS
// -----------------------------------------------------------------------------------------------------
void
DualMortarTri6::Initialize()
{
  Alpha = 0;
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------
int
DualMortarTri6::nNodes() { return 6; }

int
DualMortarTri6::nMortarShapeFct() { return 6; }

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
DualMortarTri6::ComputeDualCoeffs(CoordSet &cs)
{
  if(Alpha==0) Alpha = new FullM(6);
  else { delete Alpha; Alpha = new FullM(6); }
  Alpha->zero();

  // 1) compute scalar mass matrix (M) of supporting face element
  FaceElement* FaceElem = GetPtrMasterFace();
  FullM M = FaceElem->ScalarMass(cs, 1.0, 2);
  // 2) compute std shape function integrals (b) over supporting face element
  double ShapeIntg[6];
  FaceElem->IntegrateShapeFcts(ShapeIntg, cs, 1.0, 2);
  // 3) solve M.alpha=b to get the dual shape fcts coeffs alpha
  M.factor();
  for(int j=0; j<6; j++) {
    (*Alpha)[j][j] = ShapeIntg[j];
    M.reSolve((*Alpha)[j]); 
  }
}

#if !defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1300
// !!! ASSUME THAT DUAL MORTAR COEFFS Alpha HAVE ALREADY BEEN COMPUTED !!!
template<>
void
DualMortarTri6::GetShapeFctVal(double* Shape, double* m)
{
  double StdShape[6];
  GetPtrMasterFace()->GetShapeFctVal(StdShape, m);
  for(int i=0; i<6; ++i) {
    Shape[i] = 0;
    for(int j=0; j<6; j++)
      Shape[i] += (*Alpha)[i][j]*StdShape[j]; 
  }
}
#endif

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
DualMortarTri6::GetShapeFctVal(double* Shape, double* m)
{ 
  GetShapeFctVal<double>(Shape, m); 
}
