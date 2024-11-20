// ------------------------------------------------
// HB - 08-25-3003 
// Last modif: 09/01/2003
// ------------------------------------------------
// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>

// FEM headers
#include <Mortar.d/FaceElement.d/FaceElement.h>

#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/MortarElement.d/MortarQuad4.d/StdMortarQuad4.h>
#include <Mortar.d/MortarElement.d/MortarQuad4.d/DualMortarQuad4.h>

#include <Mortar.d/MortarElement.d/MortarQuad8.d/StdMortarQuad8.h>
#include <Mortar.d/MortarElement.d/MortarQuad8.d/DualMortarQuad8.h>

#include <Mortar.d/MortarElement.d/MortarQuad9.d/StdMortarQuad9.h>
#include <Mortar.d/MortarElement.d/MortarQuad9.d/DualMortarQuad9.h>

#include <Mortar.d/MortarElement.d/MortarQuad12.d/StdMortarQuad12.h>
#include <Mortar.d/MortarElement.d/MortarQuad12.d/DualMortarQuad12.h>

#include <Mortar.d/MortarElement.d/MortarTri3.d/StdMortarTri3.h>
#include <Mortar.d/MortarElement.d/MortarTri3.d/DualMortarTri3.h>

#include <Mortar.d/MortarElement.d/MortarTri6.d/StdMortarTri6.h>
#include <Mortar.d/MortarElement.d/MortarTri6.d/DualMortarTri6.h>

#include <Mortar.d/MortarElement.d/MortarTri10.d/StdMortarTri10.h>
#include <Mortar.d/MortarElement.d/MortarTri10.d/DualMortarTri10.h>

MortarElement*
CreateMortarElement(FaceElement* FaceElem, CoordSet& cs, bool DualFlag=false)
{
   MortarElement* MortarElem = 0;

   int FaceElemType = FaceElem->GetFaceElemType(); 

   switch(FaceElemType)
   {
     case FaceElement::QUADFACEL4:
     {  
       if(DualFlag) MortarElem = new DualMortarQuad4(FaceElem, cs);
       else         MortarElem = new StdMortarQuad4(FaceElem);
       break;
     }

     case FaceElement::QUADFACEQ8:
     {
       if(DualFlag) MortarElem = new DualMortarQuad8(FaceElem, cs);
       else         MortarElem = new StdMortarQuad8(FaceElem);
       break;
     }

     case FaceElement::QUADFACEQ9:
     {
       if(DualFlag) MortarElem = new DualMortarQuad9(FaceElem, cs);
       else         MortarElem = new StdMortarQuad9(FaceElem);
       break;
     }

     case FaceElement::QUADFACEC12:
     {
       if(DualFlag) MortarElem = new DualMortarQuad12(FaceElem, cs);
       else         MortarElem = new StdMortarQuad12(FaceElem);
       break;
     }

     case FaceElement::TRIFACEL3:
     {  
       if(DualFlag) MortarElem = new DualMortarTri3(FaceElem, cs);
       else         MortarElem = new StdMortarTri3(FaceElem);
       break;
     }

     case FaceElement::TRIFACEQ6:
     {
       if(DualFlag) MortarElem = new DualMortarTri6(FaceElem, cs);
       else         MortarElem = new StdMortarTri6(FaceElem);
       break;
     }

     case FaceElement::TRIFACEC10:
     {
       if(DualFlag) MortarElem = new DualMortarTri10(FaceElem, cs);
       else         MortarElem = new StdMortarTri10(FaceElem);
       break;
     }
                                                                                                                                                                                 
     default:
       fprintf(stderr," *** WARNING: In CreateMortarElement: Face element type %32d" 
                      "is NOT SUPPORTED.\n",FaceElemType); 
       exit(-1);
   }

  return(MortarElem);
}

