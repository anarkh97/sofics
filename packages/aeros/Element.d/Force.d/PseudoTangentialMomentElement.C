#ifdef USE_EIGEN3
#include <Element.d/Force.d/PseudoTangentialMomentElement.h>
#include <Element.d/Force.d/PotentialFunctionElementImpl.h>
#include <Corotational.d/GeomState.h>

const DofSet PseudoTangentialMomentElement::NODALDOFS[1] = { DofSet::XYZrot };

PseudoTangentialMomentElement::PseudoTangentialMomentElement(int* _nn)
 : PotentialFunctionElement<Simo::PseudoTangentialMomentWorkFunction>(1, const_cast<DofSet*>(NODALDOFS),  _nn)
{
  C0 = 0;
}

PseudoTangentialMomentElement::~PseudoTangentialMomentElement()
{
  if(C0) delete [] C0;
}

void
PseudoTangentialMomentElement::setFrame(EFrame *elemframe)
{
  C0 = new double[3][3];
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      C0[i][j] = (*elemframe)[i][j];
}

void
PseudoTangentialMomentElement::getConstants(const CoordSet& cs, Eigen::Array<double,16,1>& sconst, Eigen::Array<int,0,1>&, const GeomState* gs,
                                            double t) const
{
  if(!C0) { std::cerr << "ERROR: Element 149 (pseudo tangential moment) requires an EFRAME\n"; exit(-1); }

  if(gs) {
    sconst << prop->F0,
              C0[0][0], C0[0][1], C0[0][2],
              C0[1][0], C0[1][1], C0[1][2],
              (*gs)[nn[0]].R[0][0], (*gs)[nn[0]].R[0][1], (*gs)[nn[0]].R[0][2],
              (*gs)[nn[0]].R[1][0], (*gs)[nn[0]].R[1][1], (*gs)[nn[0]].R[1][2],
              (*gs)[nn[0]].R[2][0], (*gs)[nn[0]].R[2][1], (*gs)[nn[0]].R[2][2];
  }
  else {
    sconst << prop->Mx,
              C0[0][0], C0[0][1], C0[0][2],
              C0[1][0], C0[1][1], C0[1][2],
              1, 0, 0,
              0, 1, 0,
              0, 0, 1;
  }
}
#endif
