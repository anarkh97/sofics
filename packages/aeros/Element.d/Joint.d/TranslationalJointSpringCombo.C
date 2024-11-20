#ifdef USE_EIGEN3
#include <Element.d/Joint.d/TranslationalJointSpringCombo.h>
#include <Element.d/Joint.d/TranslationalJoint.h>
#include <Element.d/Joint.d/NonlinearTranslationalSpring.h>

TranslationalJointSpringCombo::TranslationalJointSpringCombo(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 4;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new TranslationalJoint(nnloc);          // ↓ propIndex
  subElems[1] = new NonlinearTranslationalSpring(nnloc, 0, 0);
  subElems[2] = new NonlinearTranslationalSpring(nnloc, 1, 1);
  subElems[3] = new NonlinearTranslationalSpring(nnloc, 2, 2);
}

int 
TranslationalJointSpringCombo::getTopNumber() const
{ 
  return 106; 
}
#endif
