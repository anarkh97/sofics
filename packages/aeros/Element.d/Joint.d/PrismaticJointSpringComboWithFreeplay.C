#ifdef USE_EIGEN3
#include <Element.d/Joint.d/PrismaticJointSpringComboWithFreeplay.h>
#include <Element.d/Joint.d/PrismaticJointSpringCombo.h>
#include <Element.d/Joint.d/NonlinearTranslationalSpring.h>

PrismaticJointSpringComboWithFreeplay::PrismaticJointSpringComboWithFreeplay(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 3;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new PrismaticJointSpringCombo(nnloc);   // ↓ propIndex
  subElems[1] = new NonlinearTranslationalSpring(nnloc, 0, 0, 1, 1);
  subElems[2] = new NonlinearTranslationalSpring(nnloc, 0, 0, 1, 2);
}

int 
PrismaticJointSpringComboWithFreeplay::getTopNumber() const
{ 
  return 106; 
}
#endif
