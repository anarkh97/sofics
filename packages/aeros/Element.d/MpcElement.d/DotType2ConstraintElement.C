#ifdef USE_EIGEN3
#include <Element.d/MpcElement.d/DotType2ConstraintElement.h>
#include <Corotational.d/GeomState.h>
#include <iostream>

const DofSet DotType2ConstraintElement::NODALDOFS[2] = { DofSet::XYZdisp | DofSet::XYZrot, DofSet::XYZdisp };

DotType2ConstraintElement::DotType2ConstraintElement(int* _nn, int _axis, int _type, int _ieqtype)
 : ConstraintFunctionElement<Simo::DotType2ConstraintFunction>(2, const_cast<DofSet*>(NODALDOFS), _nn, _type)
{
  axis = _axis;
  C0 = 0;
  d0 = 0;
  ieqtype = _ieqtype;
}

DotType2ConstraintElement::~DotType2ConstraintElement()
{
  if(C0) delete [] C0;
}

void
DotType2ConstraintElement::getConstants(const CoordSet & cs,
                                        Eigen::Array<double,7,1>& sconst, Eigen::Array<int,1,1>& iconst,
                                        const GeomState *gs) const
{
  // gs is the current configuration for eulerian description of rotations
  // if gs == NULL the current configuration is identical to ihe undeformed configuration
  if(gs == NULL) {
    Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > C0(&DotType2ConstraintElement::C0[0][0]);
    sconst << cs[nn[1]]->x - cs[nn[0]]->x, cs[nn[1]]->y - cs[nn[0]]->y, cs[nn[1]]->z - cs[nn[0]]->z,
              C0(axis,0), C0(axis,1), C0(axis,2), d0;
    iconst << ieqtype;
  }
  else {
    Eigen::Map<const Eigen::Matrix<double,3,3,Eigen::RowMajor> > C0(&DotType2ConstraintElement::C0[0][0]),
                                                           R1(&(*gs)[nn[0]].R[0][0]);

    Eigen::Matrix<double,1,3> b0 = C0.row(axis)*R1.transpose();
    sconst << cs[nn[1]]->x - cs[nn[0]]->x, cs[nn[1]]->y - cs[nn[0]]->y, cs[nn[1]]->z - cs[nn[0]]->z,
              b0[0], b0[1], b0[2], d0;
    iconst << ieqtype;
  }
}

void
DotType2ConstraintElement::setFrame(EFrame *elemframe)
{
  C0 = new double[3][3];
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      C0[i][j] = (*elemframe)[i][j];
}

void 
DotType2ConstraintElement::buildFrame(CoordSet& cs)
{
  if(!C0) {
    std::cerr << " *** ERROR: A reference frame needs to be specified using EFRAMES for element " << glNum+1 << ". Exiting...\n";
    exit(-1);
  }

  ConstraintFunctionElement<Simo::DotType2ConstraintFunction>::buildFrame(cs);
}

double
DotType2ConstraintElement::getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double)
{
  return 0.;
}

double
DotType2ConstraintElement::getAccelerationConstraintRhs(GeomState* refState, GeomState& gState, CoordSet& cs, double t)
{
  return MpcElement::getAccelerationConstraintRhs(refState, gState, cs, t);
}

#endif
