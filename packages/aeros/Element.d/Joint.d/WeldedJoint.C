#ifdef USE_EIGEN3
#include <Element.d/Joint.d/WeldedJoint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/CommonPointConstraint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/ParallelAxesConstraint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/RotationBlockerConstraint.h>

WeldedJoint::WeldedJoint(int* _nn)
 : SuperElement(true), elemframe(NULL)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 3;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new CommonPointConstraint(nnloc);
  subElems[1] = new ParallelAxesConstraint(nnloc);
  subElems[2] = new RotationBlockerConstraint(nnloc, 2, 1);
}

WeldedJoint::~WeldedJoint()
{
  if(elemframe && myframe) {
    delete [] elemframe;
  }
}

void
WeldedJoint::setFrame(EFrame *_elemframe)
{
  elemframe = _elemframe;
  myframe = false;
}

void
WeldedJoint::buildFrame(CoordSet &cs)
{
  if(!elemframe) {
    // if no element frame is specified, generate and use default frame
    elemframe = reinterpret_cast<double (*)[3][3]>(new double[9]);
    myframe = true;
    (*elemframe)[0][0] = 1; (*elemframe)[0][1] = 0; (*elemframe)[0][2] = 0;
    (*elemframe)[1][0] = 0; (*elemframe)[1][1] = 1; (*elemframe)[1][2] = 0;
    (*elemframe)[2][0] = 0; (*elemframe)[2][1] = 0; (*elemframe)[2][2] = 1;
  }
  SuperElement::setFrame(elemframe);
  SuperElement::buildFrame(cs);
}

int 
WeldedJoint::getTopNumber() const
{ 
  return 106; 
}
#endif
