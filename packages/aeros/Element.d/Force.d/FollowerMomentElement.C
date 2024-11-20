#ifdef USE_EIGEN3
#include <Element.d/Force.d/FollowerMomentElement.h>
#include <Element.d/Joint.d/ElementaryFunction.h>

const DofSet FollowerMomentElement::NODALINPUTDOFS[1] = { DofSet::XYZrot };
const DofSet FollowerMomentElement::NODALOUTPUTDOFS[1] = { DofSet::XYZrot };

FollowerMomentElement::FollowerMomentElement(int* _nn)
 : ForceFunctionElement<Simo::FollowerMomentForceFunction>(1, const_cast<DofSet*>(NODALINPUTDOFS),
                                                           const_cast<DofSet*>(NODALOUTPUTDOFS), _nn)
{
  C0 = NULL;
}

FollowerMomentElement::~FollowerMomentElement()
{
  if(C0) delete C0;
}

void
FollowerMomentElement::setFrame(EFrame *elemframe)
{
  C0 = new Eigen::Matrix3d;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      C0->coeffRef(i,j) = (*elemframe)[i][j];
}

void
FollowerMomentElement::getConstants(const CoordSet& cs, Eigen::Array<double,12,1>& sconst, Eigen::Array<int,0,1>&,
                                    const GeomState* gs, double t) const
{
  Eigen::Vector3d M;
  if(!C0) {
    // in this case the components of M are given in the global frame
    M << prop->Mx, prop->My, prop->Mz;
  }
  else {
    // in this case M is defined as M0*l.cross(k), where l and k are the second and third axes of the element frame
    ElementaryFunction f(prop->funtype, prop->amplitude, prop->offset, prop->c1, prop->c2, prop->c3, prop->c4);
    double M0 = f(t);
    M = M0*(C0->row(1).cross(C0->row(2))).normalized();
  }

  if(gs) {
    sconst << M[0], M[1], M[2],
              (*gs)[nn[0]].R[0][0], (*gs)[nn[0]].R[0][1], (*gs)[nn[0]].R[0][2],
              (*gs)[nn[0]].R[1][0], (*gs)[nn[0]].R[1][1], (*gs)[nn[0]].R[1][2],
              (*gs)[nn[0]].R[2][0], (*gs)[nn[0]].R[2][1], (*gs)[nn[0]].R[2][2];
  }
  else {
    sconst << M[0], M[1], M[2],
              1, 0, 0,
              0, 1, 0,
              0, 0, 1;
  }
}
#endif
