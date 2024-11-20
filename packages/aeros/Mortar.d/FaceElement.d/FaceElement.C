// ------------------------------------------------------------
// HB -  08/11/03
// ------------------------------------------------------------
// Std C/C++ lib
#include <cstdio>
#include <cstdlib>
#include <iostream>

// FEM headers
#include <Utils.d/resize_array.h>
#include <Utils.d/Connectivity.h>

#include <Mortar.d/FaceElement.d/FaceElement.h>

// -----------------------------------------------------------------------------------------------------
//                                        MAPPING & SHAPE FCT METHODS 
// -----------------------------------------------------------------------------------------------------
double* 
FaceElement::ViewRefCoords()
{
  fprintf(stderr," *** ERROR: method ViewRefCoords() NOT implemented for face el. type %d\n. Abort",
  GetFaceElemType());
  exit(1);
  return NULL;
}

/*
double* 
FaceElement::ViewRefCoords()
{
 switch(GetFaceElemType())
 {
   case FaceElement::QUADFACEL4:
     return(&FaceQuad4::RefCoords);

   case FaceElement::QUADFACEQ8:
     return(&FaceQuad8::RefCoords);

   case FaceElement::QUADFACEQ9:
     return(&FaceQuad9::RefCoords);

   case FaceElement::QUADFACEC12:
     return(&FaceQuad12::RefCoords);

   case FaceElement::TRIFACEL3:
     return(&FaceTri3::RefCoords);

   case FaceElement::TRIFACEQ6:
     return(&FaceTri6::RefCoords);

   case FaceElement::TRIFACEC10:
     return(&FaceTri10::RefCoords);

   default:
     fprintf(stderr," *** ERROR: method ViewRefCoords() NOT implemented for face el. type %d\n. Abort",
     GetFaceElemType()));
     exit(-1);
     return;
  }
}
*/

int
FaceElement::findEle(Connectivity *nodeToElem, int *eleTouch,
                     int *eleCount, int myNum, int *fnId)
{
  int *nn = new int[nNodes()];
  GetNodes(nn);
  for(int i = 0; i<nNodes(); i++) {
    for(int iele = 0; iele < nodeToElem->num(fnId[nn[i]]); iele++) {
      int eleNum = (*nodeToElem)[fnId[nn[i]]][iele];
      if (eleTouch [eleNum] != myNum) {
        eleTouch [eleNum] = myNum;
        eleCount [eleNum] = 1;
      }
      else {
        eleCount[eleNum]++;
        if(eleCount[eleNum] == nNodes()) {
          delete [] nn;
          return eleNum;
        }
      }
    }
  }
  delete [] nn;
  return -1;
}

