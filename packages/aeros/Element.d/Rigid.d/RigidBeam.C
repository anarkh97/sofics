#ifdef USE_EIGEN3
#include <Element.d/Rigid.d/RigidBeam.h>
#include <Element.d/Joint.d/BuildingBlocks.d/ConstantDistanceConstraint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/StraightLinePointFollowerConstraint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/RotationBlockerConstraint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/ParallelAxesConstraint.h>
#include <Element.d/MpcElement.d/DotType2ConstraintElement.h>
#include <Element.d/DiscreteMass.d/DiscreteMass6Dof.h>
#include <Corotational.d/utilities.h>

RigidBeam::RigidBeam(int* _nn, int _variant)
 : SuperElement(true),
   elemframe(0),
   variant(_variant)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 4; 
  subElems = new Element * [nSubElems];
  int indices[2] = { 0, 1 };

  switch(variant) {
    case 0: {
      subElems[0] = new ConstantDistanceConstraint(indices);
      subElems[1] = new StraightLinePointFollowerConstraint(indices);
      subElems[2] = new RotationBlockerConstraint(indices, 2, 1);
      subElems[3] = new ParallelAxesConstraint(indices);
    } break;
    case 1: {
      subElems[0] = new DotType2ConstraintElement(indices, 0);
      subElems[1] = new StraightLinePointFollowerConstraint(indices);
      subElems[2] = new RotationBlockerConstraint(indices, 2, 1);
      subElems[3] = new ParallelAxesConstraint(indices);
    } break;
  }  
}

void
RigidBeam::buildFrame(CoordSet& cs)
{
  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  getLength(cs, length);
  if(elemframe != 0) {
    if(length == 0) {
      if(variant != 1) {
        std::cerr << " *** ERROR: Massless rigid beam type 132 element number " << getGlNum()+1 << " has zero length,\n"; 
        std::cerr << " ***        use massless rigid beam type 133 or welded joint type 119 instead.\n";
        exit(-1);
      }
    }
    else {
      EFrame &theFrame = *elemframe;
      theFrame[0][0] = nd2.x-nd1.x;
      theFrame[0][1] = nd2.y-nd1.y;
      theFrame[0][2] = nd2.z-nd1.z;
      normalize(theFrame[0]);
      crossprod(theFrame[0],theFrame[1],theFrame[2]);
      normalize(theFrame[2]);
      crossprod(theFrame[2],theFrame[0],theFrame[1]);
    }
  }
  else {
    if(length == 0) {
      if(variant == 1) {
        c0[0][0] = 1.0;
        c0[0][1] = 0.0;
        c0[0][2] = 0.0;
        c0[1][0] = 0.0;
        c0[1][1] = 1.0;
        c0[1][2] = 0.0;
        c0[2][0] = 0.0;
        c0[2][1] = 0.0;
        c0[2][2] = 1.0;
      }
      else {
        std::cerr << " *** ERROR: Massless rigid beam type 132 element number " << getGlNum()+1 << " has zero length,\n";         
        std::cerr << " ***        use massless rigid beam type 133 or welded joint type 119 instead.\n";
        exit(-1);
      }
    }
    else {
      c0[0][0] = nd2.x-nd1.x;
      c0[0][1] = nd2.y-nd1.y;
      c0[0][2] = nd2.z-nd1.z;
      normalize(c0[0]);
      double N1 = sqrt( c0[0][0]*c0[0][0] + c0[0][1]*c0[0][1] );
      double N2 = sqrt( c0[0][0]*c0[0][0] + c0[0][2]*c0[0][2] );

      if (N1 > N2) {
        c0[1][0] = -c0[0][1]/N1;
        c0[1][1] = c0[0][0]/N1;
        c0[1][2] = 0.0;
      }
      else {
        c0[1][0] = c0[0][2]/N2;
        c0[1][1] = 0.0;
        c0[1][2] = -c0[0][0]/N2;
      }

      c0[2][0] = c0[0][1] * c0[1][2] - c0[0][2] * c0[1][1];
      c0[2][1] = c0[0][2] * c0[1][0] - c0[0][0] * c0[1][2];
      c0[2][2] = c0[0][0] * c0[1][1] - c0[0][1] * c0[1][0];
    }

    elemframe = &c0;
  }

  if(variant == 1) {
    dynamic_cast<DotType2ConstraintElement*>(subElems[0])->setConstantTerm(length);
  }

  SuperElement::setFrame(elemframe);
  SuperElement::buildFrame(cs);
}

void
RigidBeam::getLength(CoordSet& cs, double &length)
{
  // Returns length of element

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);

  double x[2], y[2], z[2];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

  double dx = x[1] - x[0];
  double dy = y[1] - y[0];
  double dz = z[1] - z[0];

  length = sqrt(dx*dx + dy*dy + dz*dz);
}

RigidBeamWithMass::RigidBeamWithMass(int* _nn, int _variant)
 : SuperElement(true),
   elemframe(0),
   variant(_variant)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 6; 
  subElems = new Element * [nSubElems];
  int indices[2] = { 0, 1 };

  switch(variant) {
    case 0: {
      subElems[0] = new ConstantDistanceConstraint(indices);
      subElems[1] = new StraightLinePointFollowerConstraint(indices);
      subElems[2] = new RotationBlockerConstraint(indices, 2, 1);
      subElems[3] = new ParallelAxesConstraint(indices);
      subElems[4] = new DiscreteMass6Dof(&(indices[0]));
      subElems[5] = new DiscreteMass6Dof(&(indices[1]));
    } break;
    case 1: {
      subElems[0] = new DotType2ConstraintElement(indices, 0);
      subElems[1] = new StraightLinePointFollowerConstraint(indices);
      subElems[2] = new RotationBlockerConstraint(indices, 2, 1);
      subElems[3] = new ParallelAxesConstraint(indices);
      subElems[4] = new DiscreteMass6Dof(&(indices[0]));
      subElems[5] = new DiscreteMass6Dof(&(indices[1]));
    } break;
  }  
}

void
RigidBeamWithMass::buildFrame(CoordSet& cs)
{
  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  getLength(cs, length);
  if(elemframe != 0) {
    if(length == 0) {
      if(variant != 1) {
        std::cerr << " *** ERROR: Rigid beam type 66 element number " << getGlNum()+1 << " has zero length,\n"; 
        std::cerr << " ***        use rigid beam type 106 or welded joint type 119 instead.\n";
        exit(-1);
      }
    }
    else {
      EFrame &theFrame = *elemframe;
      theFrame[0][0] = nd2.x-nd1.x;
      theFrame[0][1] = nd2.y-nd1.y;
      theFrame[0][2] = nd2.z-nd1.z;
      normalize(theFrame[0]);
      crossprod(theFrame[0],theFrame[1],theFrame[2]);
      normalize(theFrame[2]);
      crossprod(theFrame[2],theFrame[0],theFrame[1]);
    }
  }
  else {
    if(length == 0) {
      if(variant == 1) {
        c0[0][0] = 1.0;
        c0[0][1] = 0.0;
        c0[0][2] = 0.0;
        c0[1][0] = 0.0;
        c0[1][1] = 1.0;
        c0[1][2] = 0.0;
        c0[2][0] = 0.0;
        c0[2][1] = 0.0;
        c0[2][2] = 1.0;
      }
      else {
        std::cerr << " *** ERROR: Rigid beam type 66 element number " << getGlNum()+1 << " has zero length,\n";         
        std::cerr << " ***        use rigid beam type 106 or welded joint type 119 instead.\n";
        exit(-1);
      }
    }
    else {
      c0[0][0] = nd2.x-nd1.x;
      c0[0][1] = nd2.y-nd1.y;
      c0[0][2] = nd2.z-nd1.z;
      normalize(c0[0]);
      double N1 = sqrt( c0[0][0]*c0[0][0] + c0[0][1]*c0[0][1] );
      double N2 = sqrt( c0[0][0]*c0[0][0] + c0[0][2]*c0[0][2] );

      if (N1 > N2) {
        c0[1][0] = -c0[0][1]/N1;
        c0[1][1] = c0[0][0]/N1;
        c0[1][2] = 0.0;
      }
      else {
        c0[1][0] = c0[0][2]/N2;
        c0[1][1] = 0.0;
        c0[1][2] = -c0[0][0]/N2;
      }

      c0[2][0] = c0[0][1] * c0[1][2] - c0[0][2] * c0[1][1];
      c0[2][1] = c0[0][2] * c0[1][0] - c0[0][0] * c0[1][2];
      c0[2][2] = c0[0][0] * c0[1][1] - c0[0][1] * c0[1][0];
    }

    elemframe = &c0;
  }

  if(variant == 1) {
    dynamic_cast<DotType2ConstraintElement*>(subElems[0])->setConstantTerm(length);
  }

  SuperElement::setFrame(elemframe);
  SuperElement::buildFrame(cs);
}

void
RigidBeamWithMass::setProp(StructProp *_prop, bool _myProp)
{
  if(myProp && prop) delete prop;

  prop = _prop;
  myProp = _myProp;
  subElems[0]->setProp(prop, false);
  subElems[1]->setProp(prop, false);
  subElems[2]->setProp(prop, false);
  subElems[3]->setProp(prop, false);

  StructProp *p1 = new StructProp(*prop);
  StructProp *p2 = new StructProp(*prop);
  // compute the mass and inertia for each half of the beam, assuming the shape
  // is a solid cylinder with radius = sqrt(A/pi)
  // The inertia tensor and offset are defined in the element frame
  double l, m, r, Ixx, Iyy, Izz;
  l = length/2;
  m = prop->rho*prop->A*l;
  Ixx = 0.5*m*prop->A/M_PI;
  Iyy = Izz = 0.5*Ixx + 1/12.*m*l*l;
  p1->rho = p2->rho = m;
  p1->Ixx = p2->Ixx = Ixx;
  p1->Iyy = p2->Iyy = Iyy;
  p1->Izz = p2->Izz = Izz;
  p1->cx = l/2;
  p2->cx = -l/2;
  p1->cy = p2->cy = 0;
  p1->cz = p2->cz = 0;
  subElems[4]->setProp(p1, true);
  subElems[5]->setProp(p2, true);
  SuperElement::makeAllDOFs();
}

void
RigidBeamWithMass::getLength(CoordSet& cs, double &length)
{
  // Returns length of element

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);

  double x[2], y[2], z[2];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

  double dx = x[1] - x[0];
  double dy = y[1] - y[0];
  double dz = z[1] - z[0];

  length = sqrt(dx*dx + dy*dy + dz*dz);
}

double
RigidBeamWithMass::computeStabilityTimeStep(FullSquareMatrix &K, FullSquareMatrix &M, CoordSet &cs, GeomState *gs,
                                            double stable_tol, int stable_maxit)
{
  if(length == 0) return std::numeric_limits<double>::infinity();
  else return Element::computeStabilityTimeStep(K, M, cs, gs, stable_tol, stable_maxit);
}

#endif
