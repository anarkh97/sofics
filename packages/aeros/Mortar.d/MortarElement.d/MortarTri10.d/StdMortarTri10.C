// ----------------------------------------------------------------
// HB - 05/24/05
// ----------------------------------------------------------------

// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>

// FEM headers
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/MortarElement.d/MortarTri10.d/StdMortarTri10.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
StdMortarTri10::StdMortarTri10() {} 
// call the default base class constructor MortarElement::MortarElement() 

StdMortarTri10::StdMortarTri10(FaceElement* FaceElem)
{
  SetPtrMasterFace(FaceElem);
}

StdMortarTri10::StdMortarTri10(double _area, FaceElement* FaceElem)
{
  SetArea(_area);
  SetPtrMasterFace(FaceElem);
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------
int
StdMortarTri10::nNodes() { return 10; }

int
StdMortarTri10::nMortarShapeFct() { return 10; }

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
#if !defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1300
template<>
void
StdMortarTri10::GetShapeFctVal(double* Shape, double* m)
{
  GetPtrMasterFace()->GetShapeFctVal(Shape, m);
}

template<>
void
StdMortarTri10::GetdShapeFct(double* dShapex, double* dShapey, double* m)
{
  GetPtrMasterFace()->GetdShapeFct(dShapex, dShapey, m);
}

template<>
void
StdMortarTri10::Getd2ShapeFct(double* d2Shapex, double* d2Shapey, double* d2Shapexy, double* m)
{
  GetPtrMasterFace()->Getd2ShapeFct(d2Shapex, d2Shapey, d2Shapexy, m);
}
#endif

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
StdMortarTri10::GetShapeFctVal(double* Shape, double* m)
{
  GetShapeFctVal<double>(Shape, m); 
}

void
StdMortarTri10::GetdShapeFct(double* dShapex, double* dShapey, double* m)
{
  GetdShapeFct<double>(dShapex, dShapey, m);
}

void
StdMortarTri10::Getd2ShapeFct(double* d2Shapex, double* d2Shapey, double* d2Shapexy, double* m)
{
  Getd2ShapeFct<double>(d2Shapex, d2Shapey, d2Shapexy, m);
}
