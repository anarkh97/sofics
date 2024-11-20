#ifdef USE_EIGEN3
#include <Element.d/MpcElement.d/DotType3ConstraintElement.h>
#include <Corotational.d/GeomState.h>
#include <iostream>

DotType3ConstraintElement::DotType3ConstraintElement(int* _nn, int _axis)
		: ConstraintFunctionElement<Simo::DotType3ConstraintFunction>(2, DofSet::XYZdisp | DofSet::XYZrot, _nn, 0)
{
	axis = _axis;
	C0 = 0;
	d0 = 0;
}

DotType3ConstraintElement::~DotType3ConstraintElement()
{
	if(C0) delete [] C0;
}

void
DotType3ConstraintElement::getConstants(const CoordSet & cs,
                                        Eigen::Array<double,10,1>& sconst, Eigen::Array<int,0,1>&,
                                        const GeomState *gs) const
{
	// gs is the current configuration for eulerian description of rotations
	// if gs == NULL the current configuration is identical to ihe undeformed configuration
	if(gs == NULL) {
		Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > C0(&DotType3ConstraintElement::C0[0][0]);
		sconst << cs[nn[1]]->x - cs[nn[0]]->x, cs[nn[1]]->y - cs[nn[0]]->y, cs[nn[1]]->z - cs[nn[0]]->z,
				C0(axis,0), C0(axis,1), C0(axis,2),
				C0(axis,0), C0(axis,1), C0(axis,2), d0;
	}
	else {
		Eigen::Map<const Eigen::Matrix<double,3,3,Eigen::RowMajor> > C0(&DotType3ConstraintElement::C0[0][0]),
				R1(&(*gs)[nn[0]].R[0][0]),
				R2(&(*gs)[nn[1]].R[0][0]);
		Eigen::Matrix<double,1,3> b0 = C0.row(axis)*R1.transpose(),
				c0 = C0.row(axis)*R2.transpose();
		sconst << cs[nn[1]]->x - cs[nn[0]]->x, cs[nn[1]]->y - cs[nn[0]]->y, cs[nn[1]]->z - cs[nn[0]]->z,
				b0[0], b0[1], b0[2],
				c0[0], c0[1], c0[2], d0;
	}
}

void
DotType3ConstraintElement::setFrame(EFrame *elemframe)
{
	C0 = new double[3][3];
	for(int i = 0; i < 3; ++i)
		for(int j = 0; j < 3; ++j)
			C0[i][j] = (*elemframe)[i][j];
}

void
DotType3ConstraintElement::buildFrame(CoordSet& cs)
{
	if(!C0) {
		std::cerr << " *** ERROR: A reference frame needs to be specified using EFRAMES for element " << glNum+1 << ". Exiting...\n";
		exit(-1);
	}

	ConstraintFunctionElement<Simo::DotType3ConstraintFunction>::buildFrame(cs);
}

double
DotType3ConstraintElement::getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double)
{
	return 0.;
}

double
DotType3ConstraintElement::getAccelerationConstraintRhs(GeomState* refState, GeomState& gState, CoordSet& cs, double t)
{
	return MpcElement::getAccelerationConstraintRhs(refState, gState, cs, t);
}

#endif
