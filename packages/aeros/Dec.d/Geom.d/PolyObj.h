#ifndef _POLYOBJ_H
#define _POLYOBJ_H

#include <Utils.d/resize_array.h>

class PolyObj {

   ResizeArray<int> indices;
   ResizeArray<int> nodeNum;
   float (*xyz)[3];
   int currentStep;
   int numSteps;
   int numPoly;

 public :

   PolyObj();

   void  addTri(int, int, int);
   void  addQuad(int, int, int, int);
   int   get_num_steps();
   void  set_step(int);
   float getTime(int);
   void  setCoord(float (*_xyz)[3]) { xyz = _xyz; }
};

#endif
