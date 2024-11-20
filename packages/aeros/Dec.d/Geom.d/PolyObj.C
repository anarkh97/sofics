#include "PolyObj.h"

PolyObj::PolyObj() : indices(0), nodeNum(0)
{
  numPoly = 0;
  indices[0] = 0; // first polygon starts at the begining
}

void
PolyObj::addTri(int num1, int num2, int num3)
{
  nodeNum[indices[numPoly]] = num1;
  nodeNum[indices[numPoly]+1] = num2;
  nodeNum[indices[numPoly]+2] = num3;

  indices[numPoly+1] = indices[numPoly]+3;
  numPoly++;
}

void
PolyObj::addQuad(int num1, int num2, int num3, int num4)
{
  nodeNum[indices[numPoly]] = num1;
  nodeNum[indices[numPoly]+1] = num2;
  nodeNum[indices[numPoly]+2] = num3;
  nodeNum[indices[numPoly]+3] = num4;

//  indices[++numPoly] = numPoly+4;
  indices[numPoly+1] = indices[numPoly]+4;
  numPoly++;
}

int
PolyObj::get_num_steps()
{
  return numSteps;
}

void
PolyObj::set_step(int _step)
{
  currentStep = _step-1;
}

/*
float
PolyObj::getTime(int v)
{
  if(v == 0) return 0;
  if(v > sr.get_num_steps())
     return 0.0;
  else
     return sr.getTimes()[v-1];
}
*/

#ifdef USE_OPENGL
#include <GL/gl.h>

void
PolyObj::GL_draw()
{
  glDisable(GL_LIGHTING);
  for(int i = 0; i < numPoly; ++i) {
    glBegin(GL_POLYGON); 
    for(int j = indices[i]; j < indices[i+1]; j++) {
//    for(int j = 0; j < 100; ++j){
      glVertex3fv(xyz[nodeNum[j]]);
    }
    glEnd();
  }
 glEnable(GL_LIGHTING);
}

#endif
