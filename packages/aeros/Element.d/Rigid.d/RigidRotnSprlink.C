#ifdef USE_EIGEN3
#include <Element.d/Rigid.d/RigidRotnSprlink.h>
#include <Element.d/MpcElement.d/MpcElement.h>

RigidRotnSprlink::RigidRotnSprlink(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
}

void
RigidRotnSprlink::setProp(StructProp* _prop, bool _myProp)
{
  nSubElems = int(_prop->kx != 0.0) + int(_prop->ky != 0.0) + int(_prop->kz != 0.0);
  subElems = new Element * [nSubElems];
  int indices[2] = { 0, 1 };
  int count = 0;
  if(_prop->kx != 0.0) {
    subElems[count++] = new MpcElement(2, DofSet::Xrot, indices);
  }
  if(_prop->ky != 0.0) {
    subElems[count++] = new MpcElement(2, DofSet::Yrot, indices);
  }
  if(_prop->kz != 0.0) {
    subElems[count++] = new MpcElement(2, DofSet::Zrot, indices);
  }
  for(int i=0; i<nSubElems; ++i) {
    LMPCons* mpc = dynamic_cast<LMPCons*>(subElems[i]);
    mpc->terms[0].coef.r_value = 1.0;
    mpc->terms[1].coef.r_value = -1.0;
  }

  SuperElement::setProp(_prop, _myProp);
}
#endif
