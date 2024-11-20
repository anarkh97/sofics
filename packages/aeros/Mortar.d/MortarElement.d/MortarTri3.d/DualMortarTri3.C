// ----------------------------------------------------------------
// HB - 08/26/03
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
#include <Mortar.d/MortarElement.d/MortarTri3.d/DualMortarTri3.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
DualMortarTri3::DualMortarTri3() 
{
  //Alpha = 0;
} 

DualMortarTri3::DualMortarTri3(FaceElement* FaceElem)
{
  //Alpha = 0;
  SetPtrMasterFace(FaceElem);
}

DualMortarTri3::DualMortarTri3(double area_, FaceElement* FaceElem)
{
  //Alpha = 0;
  SetArea(area_);
  SetPtrMasterFace(FaceElem);
}

DualMortarTri3::DualMortarTri3(FaceElement* FaceElem, CoordSet &cs)
{
  //Alpha = 0;
  SetPtrMasterFace(FaceElem);
  //ComputeDualCoeffs(cs);
}

DualMortarTri3::DualMortarTri3(double area_, FaceElement* FaceElem, CoordSet &cs)
{
  //Alpha = 0;
  SetArea(area_);
  SetPtrMasterFace(FaceElem);
  //ComputeDualCoeffs(cs);
}

// -----------------------------------------------------------------------------------------------------
//                                            DESTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
DualMortarTri3::~DualMortarTri3()
{
  //if(Alpha) { delete Alpha; Alpha = 0; };
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------
int
DualMortarTri3::nNodes() { return 3; }

int
DualMortarTri3::nMortarShapeFct() { return 3; }

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
/*
void
DualMortarTri3::ComputeDualCoeffs(CoordSet &cs)
{
  if(Alpha==0) Alpha = new FullM(4);
  else{ delete Alpha; Alpha = new FullM(4); }
  Alpha->zero();

  // 1) compute scalar mass matrix (M) of supporting face element
  FaceElement* FaceElem = GetPtrMasterFace();
  FullM M = FaceElem->ScalarMass(cs, 1.0, 2);
  // 2) compute std shape function integrals (b) over supporting face element
  double ShapeIntg[4];
  FaceElem->IntegrateShapeFcts(ShapeIntg, cs, 1.0, 2);
  // 3) solve M.alpha=b to get the dual shape fcts coeffs alpha
  M.factor();
  for(int j=0; j<4; j++) {
    (*Alpha)[j][j] = ShapeIntg[j];
    M.reSolve((*Alpha)[j]); 
  }
}
*/

#if !defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1300
template<>
void
DualMortarTri3::GetShapeFctVal(double* Shape, double* m)
{
  double StdShape[3];
  GetPtrMasterFace()->GetShapeFctVal(StdShape, m);

  Shape[0] = 3.*StdShape[0] -    StdShape[1] -    StdShape[2];
  Shape[1] =  - StdShape[0] + 3.*StdShape[1] -    StdShape[2];
  Shape[2] =  - StdShape[0] -    StdShape[1] + 3.*StdShape[2];
}
#endif

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
DualMortarTri3::GetShapeFctVal(double* Shape, double* m)
{ 
  GetShapeFctVal<double>(Shape, m); 
}
