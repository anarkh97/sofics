
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Mortar.d/MortarElement.d/MortarElement.h>

// -----------------------------------------------------------------------------------------------------
//                                            CONSTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
MortarElement::MortarElement() { Area = 0; MasterFace = 0; } 

MortarElement::MortarElement(FaceElement* PtrMasterFace)
{
  Area       = 0;
  MasterFace = PtrMasterFace;
}

MortarElement::MortarElement(double _area, FaceElement* PtrMasterFace)
{
  Area       = _area;
  MasterFace = PtrMasterFace;
}

// -----------------------------------------------------------------------------------------------------
//                                            DESCTRUCTOR METHODS
// -----------------------------------------------------------------------------------------------------
MortarElement::~MortarElement()
{
  Area       = 0;
  MasterFace = 0;
}

// -----------------------------------------------------------------------------------------------------
//                                    INITIALIZATION & CLEAR/CLEAN METHODS
// -----------------------------------------------------------------------------------------------------
void
MortarElement::Initialize()
{ 
  Area       = 0;
  MasterFace = 0;
}

// -----------------------------------------------------------------------------------------------------
//                                            SET METHODS
// -----------------------------------------------------------------------------------------------------
void
MortarElement::SetArea(double _area) { Area = _area; }

void
MortarElement::SetPtrMasterFace(FaceElement* PtrMasterFace)
{
  MasterFace = PtrMasterFace;
}

// -----------------------------------------------------------------------------------------------------
//                                            GET METHODS
// -----------------------------------------------------------------------------------------------------
// LOCAL METHODS 
// ------------- 
FaceElement*
MortarElement::GetPtrMasterFace() { return(MasterFace); }

int
MortarElement::nMasterFaceNodes() 
{ 
  if(MasterFace) return(MasterFace->nNodes());
  else return(0);
}

// --------------- 
// VIRTUAL METHODS
// --------------- 
int
MortarElement::nNodes() 
{
  fprintf(stderr," *** WARNING: base class MortarElement::nNodes() should NOT be called \n");
  fprintf(stderr," *** WARNING: currently returns the number of nodes of the MASTER face \n");
  return(nMasterFaceNodes());
}

int 
MortarElement::nMortarShapeFct()
{
  fprintf(stderr," *** WARNING: base class MortarElement::nMortarShapeFct() should NOT be called \n");
  fprintf(stderr," *** WARNING: currently returns the number of nodes of the MASTER face \n");
  return(nMasterFaceNodes());
}

// -----------------------------------------------------------------------------------------------------
//                                   MAPPING & SHAPE FUNCTION METHODS 
// -----------------------------------------------------------------------------------------------------
// VIRTUAL METHODS
// ---------------
void
MortarElement::GetShapeFctVal(double* Shape, double* m)
{
  fprintf(stderr," *** WARNING: base class MortarElement::GetShapeFctVal() is not implemented.\n");
} 

void
MortarElement::GetdShapeFct(double* dShapex, double* dShapey, double* m)
{
  fprintf(stderr," *** WARNING: base class MortarElement::GetdShapeFct() is not implemented.\n");
}

void
MortarElement::Getd2ShapeFct(double* d2Shapex, double* d2Shapey, double* d2Shapexy, double* m)
{
  fprintf(stderr," *** WARNING: base class MortarElement::Getd2ShapeFct() is not implemented.\n");
}
