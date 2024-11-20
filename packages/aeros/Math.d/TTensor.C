#include <Math.d/TTensor.h>
template <>
void
SymTensor<SymTensor<double,2>,2>::dblContractInto(const Tensor &b, Tensor *res) const
{
  SymTensor<double,2> *resStress = dynamic_cast<SymTensor<double,2> *>(res);
  if(resStress != 0) {
    const SymTensor<double,2> &bb = dynamic_cast<const SymTensor<double,2> &>(b);
    for(int i = 0; i < 3; ++i)
      (*resStress)[i] = (*this)[i][0]*bb[0]+(*this)[i][1]*bb[1]+2*(*this)[i][2]*bb[2];
  } else {
    fprintf(stderr, "Tried to double contract a Sym 4th order tensor with something unknown\n");
  }
}
