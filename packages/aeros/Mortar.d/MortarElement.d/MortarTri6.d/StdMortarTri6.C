// ----------------------------------------------------------------
// HB - 08/25/03
// ----------------------------------------------------------------

// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>

// FEM headers
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/MortarElement.d/MortarTri6.d/StdMortarTri6.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
StdMortarTri6::StdMortarTri6() {} 
// call the default base class constructor MortarElement::MortarElement() 

StdMortarTri6::StdMortarTri6(FaceElement* FaceElem)
{
  SetPtrMasterFace(FaceElem);
}

StdMortarTri6::StdMortarTri6(double _area, FaceElement* FaceElem)
{
  SetArea(_area);
  SetPtrMasterFace(FaceElem);
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------
int
StdMortarTri6::nNodes() { return 6; }

int
StdMortarTri6::nMortarShapeFct() { return 6; }

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
#if !defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1300
template<>
void
StdMortarTri6::GetShapeFctVal(double* Shape, double* m)
{
  GetPtrMasterFace()->GetShapeFctVal(Shape, m);
}

template<>
void
StdMortarTri6::GetdShapeFct(double* dShapex, double* dShapey, double* m)
{
  GetPtrMasterFace()->GetdShapeFct(dShapex, dShapey, m);
}

template<>
void
StdMortarTri6::Getd2ShapeFct(double* d2Shapex, double* d2Shapey, double* d2Shapexy, double* m)
{
  GetPtrMasterFace()->Getd2ShapeFct(d2Shapex, d2Shapey, d2Shapexy, m);
}
#endif

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
StdMortarTri6::GetShapeFctVal(double* Shape, double* m)
{
  GetShapeFctVal<double>(Shape, m); 
}

void
StdMortarTri6::GetdShapeFct(double* dShapex, double* dShapey, double* m)
{
  GetdShapeFct<double>(dShapex, dShapey, m);
}

void
StdMortarTri6::Getd2ShapeFct(double* d2Shapex, double* d2Shapey, double* d2Shapexy, double* m)
{
  Getd2ShapeFct<double>(d2Shapex, d2Shapey, d2Shapexy, m);
}
