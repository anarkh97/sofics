// ----------------------------------------------------------------
// HB - 06/09/03
// ----------------------------------------------------------------

// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>

// FEM headers
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/MortarElement.d/MortarQuad12.d/StdMortarQuad12.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
StdMortarQuad12::StdMortarQuad12() {} 
// call the default base class constructor MortarElement::MortarElement() 

StdMortarQuad12::StdMortarQuad12(FaceElement* FaceElem)
{
  SetPtrMasterFace(FaceElem);
}

StdMortarQuad12::StdMortarQuad12(double area_, FaceElement* FaceElem)
{
  SetArea(area_);
  SetPtrMasterFace(FaceElem);
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------
int
StdMortarQuad12::nNodes() { return(12); }

int
StdMortarQuad12::nMortarShapeFct() { return(12); }

// -----------------------------------------------------------------------------------------------------
//                                      MAPPING & SHAPE FUNCTION METHODS
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS
// -------------
#if !defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1300
template<>
void
StdMortarQuad12::GetShapeFctVal(double* Shape, double* m)
{
  GetPtrMasterFace()->GetShapeFctVal(Shape, m);
}

template<>
void
StdMortarQuad12::GetdShapeFct(double* dShapex, double* dShapey, double* m)
{
  GetPtrMasterFace()->GetdShapeFct(dShapex, dShapey, m);
}

template<>
void
StdMortarQuad12::Getd2ShapeFct(double* d2Shapex, double* d2Shapey, double* d2Shapexy, double* m)
{
  GetPtrMasterFace()->Getd2ShapeFct(d2Shapex, d2Shapey, d2Shapexy, m);
}
#endif

// ---------------------------------
// IMPLEMENTATION OF VIRTUAL METHODS
// ---------------------------------
void
StdMortarQuad12::GetShapeFctVal(double* Shape, double* m)
{
  GetShapeFctVal<double>(Shape, m); 
}

void
StdMortarQuad12::GetdShapeFct(double* dShapex, double* dShapey, double* m)
{
  GetdShapeFct<double>(dShapex, dShapey, m);
}

void
StdMortarQuad12::Getd2ShapeFct(double* d2Shapex, double* d2Shapey, double* d2Shapexy, double* m)
{
  Getd2ShapeFct<double>(d2Shapex, d2Shapey, d2Shapexy, m);
}
