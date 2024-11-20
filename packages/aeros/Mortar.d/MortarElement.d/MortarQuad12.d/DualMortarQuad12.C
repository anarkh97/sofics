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
#include <Mortar.d/MortarElement.d/MortarQuad12.d/DualMortarQuad12.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
DualMortarQuad12::DualMortarQuad12() 
{
  Initialize();
} 

DualMortarQuad12::DualMortarQuad12(FaceElement* FaceElem)
{
  Initialize();
  SetPtrMasterFace(FaceElem);
}

DualMortarQuad12::DualMortarQuad12(double _area, FaceElement* FaceElem)
{
  Initialize();
  SetArea(_area);
  SetPtrMasterFace(FaceElem);
}

DualMortarQuad12::DualMortarQuad12(FaceElement* FaceElem, CoordSet &cs)
{
  Initialize();
  SetPtrMasterFace(FaceElem);
  ComputeDualCoeffs(cs);
}

DualMortarQuad12::DualMortarQuad12(double _area, FaceElement* FaceElem, CoordSet &cs)
{
  Initialize();
  SetArea(_area);
  SetPtrMasterFace(FaceElem);
  ComputeDualCoeffs(cs);
}

// -----------------------------------------------------------------------------------------------------
//                                            DESTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
DualMortarQuad12::~DualMortarQuad12()
{
  if(Alpha) { delete Alpha; Alpha = 0; }
}

// -----------------------------------------------------------------------------------------------------
//                                   INITILIZATION & CLEAN/CLEAR METHODS
// -----------------------------------------------------------------------------------------------------
void
DualMortarQuad12::Initialize()
{
  Alpha = 0;
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------
int
DualMortarQuad12::nNodes() { return 12; }

int
DualMortarQuad12::nMortarShapeFct() { return 12; }

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
void
DualMortarQuad12::ComputeDualCoeffs(CoordSet &cs)
{
  if(Alpha==0) Alpha = new FullM(12);
  else { delete Alpha; Alpha = new FullM(12); }
  Alpha->zero();

  // 1) compute scalar mass matrix (M) of supporting face element
  FaceElement* FaceElem = GetPtrMasterFace();
  FullM M = FaceElem->ScalarMass(cs, 1.0, 3);
  // 2) compute std shape function integrals (b) over supporting face element
  double ShapeIntg[12];
  FaceElem->IntegrateShapeFcts(ShapeIntg, cs, 1.0, 3);
  // 3) solve M.alpha=b to get the dual shape fcts coeffs alpha
  M.factor();
  for(int j=0; j<12; j++) {
    (*Alpha)[j][j] = ShapeIntg[j];
    M.reSolve((*Alpha)[j]); 
  }
}

#if !defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1300
// !!! ASSUME THAT DUAL MORTAR COEFFS Alpha HAVE ALREADY BEEN COMPUTED !!!
template<>
void
DualMortarQuad12::GetShapeFctVal(double* Shape, double* m)
{
  double StdShape[12];
  GetPtrMasterFace()->GetShapeFctVal(StdShape, m);
  for(int i=0; i<12; ++i) {
    Shape[i] = 0;
    for(int j=0; j<12; j++)
      Shape[i] += (*Alpha)[i][j]*StdShape[j]; 
  }
}
#endif

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
DualMortarQuad12::GetShapeFctVal(double* Shape, double* m)
{ 
  GetShapeFctVal<double>(Shape, m); 
}
