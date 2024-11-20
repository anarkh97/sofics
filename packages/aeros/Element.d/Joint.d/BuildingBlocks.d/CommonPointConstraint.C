#ifdef USE_EIGEN3
#include <Element.d/Joint.d/BuildingBlocks.d/CommonPointConstraint.h>
#include <Element.d/MpcElement.d/MpcElement.h>

CommonPointConstraint::CommonPointConstraint(int* _nn)
 : SuperElement(true)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 3;
  subElems = new Element * [nSubElems];
  int nnloc[2] = { 0, 1 };
  subElems[0] = new MpcElement(2, DofSet::Xdisp, nnloc);
  subElems[1] = new MpcElement(2, DofSet::Ydisp, nnloc);
  subElems[2] = new MpcElement(2, DofSet::Zdisp, nnloc);
  for(int i=0; i<3; ++i) {
    LMPCons* mpc = dynamic_cast<LMPCons*>(subElems[i]);
    mpc->terms[0].coef.r_value = 1.0;
    mpc->terms[1].coef.r_value = -1.0;
  }
}

void
CommonPointConstraint::buildFrame(CoordSet& cs)
{
  double d[3] = { cs.getNode(nn[0]).x - cs.getNode(nn[1]).x,
                  cs.getNode(nn[0]).y - cs.getNode(nn[1]).y,
                  cs.getNode(nn[0]).z - cs.getNode(nn[1]).z };
  for(int i=0; i<3; ++i) {
    LMPCons* mpc = dynamic_cast<LMPCons*>(subElems[i]);
    mpc->rhs.r_value = -d[i];
    mpc->original_rhs = mpc->rhs;
  }
  SuperElement::buildFrame(cs);
}

int 
CommonPointConstraint::getTopNumber() const
{ 
  return 101; 
}
#endif
