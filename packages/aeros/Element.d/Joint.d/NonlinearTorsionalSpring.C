#ifdef USE_EIGEN3
#include <Element.d/Joint.d/NonlinearTorsionalSpring.h>

NonlinearTorsionalSpring::NonlinearTorsionalSpring(int* _nn, int _axis1, int _axis2, int _propIndex, int _type, int _ieqtype)
 : AngleType1ConstraintElement(_nn, _axis1, _axis2, M_PI/2, _type, _ieqtype)
{
  m_axis1 = _axis1;
  m_axis2 = _axis2;
  propIndex = _propIndex;
  offset2 = 0;
  quadrant = 0;
}

void
NonlinearTorsionalSpring::setProp(StructProp *p, bool _myProp)
{
  StructProp *prop = (_myProp) ? p : new StructProp(*p); 

  const double k[3] = { p->k1, p->k2, p->k3 };
  const int &i = propIndex;
  if(type == 1 && ieqtype == 1) {
    offset += p->freeplay[i].ul;
    prop->penalty = (p->freeplay[i].uz-p->freeplay[i].dz)*k[i];
  }
  else if(type == 1 && ieqtype == 2) {
    offset += p->freeplay[i].ll;
    prop->penalty = (p->freeplay[i].lz-p->freeplay[i].dz)*k[i];
  }
  else if(type == 0) {
    prop->penalty = p->freeplay[i].dz*k[i];
  }
  prop->lagrangeMult = false;

  AngleType1ConstraintElement::setProp(prop, true);
}

void 
NonlinearTorsionalSpring::update(GeomState *refState, GeomState& gState, CoordSet& cs, double t)
{
  // internal states
  if(numStates() > 0) {
    updateStates((GeomState *) NULL, gState, cs, 0.0);
    axis1 = (quadrant == 0 || quadrant == 2) ? m_axis1 : m_axis2, axis2 = m_axis2;
  }

  AngleType1ConstraintElement::update(refState, gState, cs, t);
}

int
NonlinearTorsionalSpring::numStates()
{
  // TODO: consider reparametrization for inequalities
  return (type == 0) ? 3 : 0;
}

void
NonlinearTorsionalSpring::initStates(double *statenp)
{
  if(numStates() == 0) return;
  statenp[0] = offset;
  statenp[1] = offset2;
  statenp[2] = quadrant;
}

void 
NonlinearTorsionalSpring::updateStates(GeomState *, GeomState &gState, CoordSet &, double dt)
{
  if(numStates() == 0) return;
  // TODO: consider if it is better to update the state from the reference state (i.e. the last converged solution)
  //       rather than the current newton iteration, as we do for plasticity
  double *statenp = gState.getElemState(getGlNum()) + stateOffset;
  offset = statenp[0];
  offset2 = statenp[1];
  quadrant = (int)statenp[2];

  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > C0(&AngleType1ConstraintElement::C0[0][0]),
                                                         R1(&gState[nn[0]].R[0][0]),
                                                         R2(&gState[nn[1]].R[0][0]);

  axis1 = (quadrant == 0 || quadrant == 2) ? m_axis1 : m_axis2;
  axis2 = m_axis2;

  Eigen::Matrix<double,1,3> a = C0.row(axis1)*R1.transpose(), b = C0.row(axis2)*R2.transpose();
  double cth = a.dot(b);

  // reparameterization
  if(cth > 0.707106781 || cth < -0.707106781) {
    if(axis1 != axis2) {
      //std::cerr << "reparameterizing #1, cth = " << cth << ", offset = " << offset << std::endl;
      if(quadrant == 0 && cth > 0) { 
        //std::cerr << "going from quadrant 0 to quadrant 1\n";
        quadrant = 1; offset2 += M_PI/2; 
      }
      else if(quadrant == 0 && cth < 0) { 
        //std::cerr << "going from quadrant 0 to quadrant 3\n";
        quadrant = 3; offset2 -= M_PI/2;
      }
      else if(quadrant == 2 && cth < 0) {
        //std::cerr << "going from quadrant 2 to quadrant 3\n";
        quadrant = 3; offset2 += M_PI/2;
      }
      else if(quadrant == 2 && cth > 0) {
        //std::cerr << "going from quadrant 2 to quadrant 1\n";
        quadrant = 1; offset2 -= M_PI/2;
      }

      if(cth > 0) offset = M_PI/2 - offset2;
      else        offset = M_PI/2 + offset2;
    }
    else {
      //std::cerr << "reparameterizing #2, cth = " << cth << ", offset = " << offset << std::endl;
      if(quadrant == 1 && cth < 0) {
        //std::cerr << "going from quadrant 1 to quadrant 2\n";
        quadrant = 2; offset2 += M_PI/2;
      }
      else if(quadrant == 1 && cth > 0) {
        //std::cerr << "going from quadrant 1 to quadrant 0\n";
        quadrant = 0; offset2 -= M_PI/2;
      }
      else if(quadrant == 3 && cth > 0) {
        //std::cerr << "going from quadrant 3 to quadrant 0\n";
        quadrant = 0; offset2 += M_PI/2;
      }
      else if(quadrant == 3 && cth < 0) {
        //std::cerr << "going from quadrant 3 to quadrant 2\n";
        quadrant = 2; offset2 -= M_PI/2;
      }

      if(cth > 0) offset = M_PI/2 + offset2;
      else        offset = M_PI/2 - offset2;
    }
  }

  statenp[0] = offset;
  statenp[1] = offset2;
  statenp[2] = quadrant;
}
#endif
