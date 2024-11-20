#ifdef USE_EIGEN3
#include <Element.d/Force.d/FollowerForceElement.h>
#include <Element.d/Joint.d/ElementaryFunction.h>

const DofSet FollowerForceElement::NODALINPUTDOFS[1] = { DofSet::XYZrot };
const DofSet FollowerForceElement::NODALOUTPUTDOFS[1] = { DofSet::XYZdisp };

FollowerForceElement::FollowerForceElement(int* _nn)
 : ForceFunctionElement<Simo::FollowerForceFunction>(1, const_cast<DofSet*>(NODALINPUTDOFS),
                                                     const_cast<DofSet*>(NODALOUTPUTDOFS), _nn)
{
  C0 = NULL;
}

FollowerForceElement::~FollowerForceElement()
{
  if(C0) delete C0;
}

void
FollowerForceElement::setFrame(EFrame *elemframe)
{
  C0 = new Eigen::Matrix3d;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      C0->coeffRef(i,j) = (*elemframe)[i][j];
}

void
FollowerForceElement::getConstants(const CoordSet& cs, Eigen::Array<double,12,1>& sconst, Eigen::Array<int,0,1>&,
                                   const GeomState* gs, double t)const
{
  Eigen::Vector3d F;
  if(!C0) {
    // in this case the components of F are given in the global frame
    F << prop->Fx, prop->Fy, prop->Fz;
  }
  else {
    // in this case F is defined as F0*l.cross(k), where l and k are the second and third axes of the element frame
    ElementaryFunction f(prop->funtype, prop->amplitude, prop->offset, prop->c1, prop->c2, prop->c3, prop->c4);
    double F0 = f(t);
    F = F0*(C0->row(1).cross(C0->row(2))).normalized();
  }

  if(gs) {
    sconst << F[0], F[1], F[2],
              (*gs)[nn[0]].R[0][0], (*gs)[nn[0]].R[0][1], (*gs)[nn[0]].R[0][2],
              (*gs)[nn[0]].R[1][0], (*gs)[nn[0]].R[1][1], (*gs)[nn[0]].R[1][2],
              (*gs)[nn[0]].R[2][0], (*gs)[nn[0]].R[2][1], (*gs)[nn[0]].R[2][2];
  }
  else {
    sconst << F[0], F[1], F[2],
              1, 0, 0,
              0, 1, 0,
              0, 0, 1;
  }
}
#endif
