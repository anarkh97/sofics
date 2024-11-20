// ------------------------------------------------
// HB - 07/01/03
// Last modif: 09/01/2003 
// ------------------------------------------------
// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>

// FEM headers
#include <Utils.d/DistHelper.h>
#include <Mortar.d/FaceElement.d/FaceElemSet.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>

#include <Mortar.d/FaceElement.d/FaceQuad4.d/FaceQuad4.h>
#include <Mortar.d/FaceElement.d/FaceQuad8.d/FaceQuad8.h>
#include <Mortar.d/FaceElement.d/FaceQuad9.d/FaceQuad9.h>
#include <Mortar.d/FaceElement.d/FaceQuad12.d/FaceQuad12.h>
#include <Mortar.d/FaceElement.d/FaceTri3.d/FaceTri3.h>
#include <Mortar.d/FaceElement.d/FaceTri6.d/FaceTri6.h>
#include <Mortar.d/FaceElement.d/FaceTri10.d/FaceTri10.h>

// ------------------------------------------------
// For input file:
// etype = 1 -> Quad4  face element
//       = 2 -> Quad8  face element
//       = 3 -> Tri3   face element
//       = 4 -> Tri6   face element
//       = 5 -> Quad9  face element
//       = 6 -> Quad12 face element
//       = 7 -> Tri10  face element
// ------------------------------------------------

void
FaceElemSet::elemadd(int num, int etype, int nnodes, int* nodes)
{
   FaceElement* ele = 0;

   switch(etype) 
   {
     case FaceElement::QUADFACEL4: 
       ele = new (ba) FaceQuad4(nodes);
       break;
     case FaceElement::QUADFACEQ8: 
       ele = new (ba) FaceQuad8(nodes);
       break;
     case FaceElement::TRIFACEL3: 
       ele = new (ba) FaceTri3(nodes);
       break;
     case FaceElement::TRIFACEQ6: 
       ele = new (ba) FaceTri6(nodes);
       break;
     case FaceElement::QUADFACEQ9: 
       ele = new (ba) FaceQuad9(nodes);
       break;
     case FaceElement::QUADFACEC12: 
       ele = new (ba) FaceQuad12(nodes);
       break;
     case FaceElement::TRIFACEC10: 
       ele = new (ba) FaceTri10(nodes);
       break;
     default:
       filePrint(stderr," *** ERROR: Face element Type %2d is NOT supported. Abort.\n", etype);
       exit(-1);
       return;
   }

   if(ele) elemadd(num, ele);
}
