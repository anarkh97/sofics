#ifdef USE_EIGEN3
#include <Element.d/MpcElement.d/AngleType1ConstraintElement.h>
#include <Corotational.d/GeomState.h>
#include <iostream>

AngleType1ConstraintElement::AngleType1ConstraintElement(int* _nn, int _axis1, int _axis2, double _offset, int _type, int _ieqtype)
		: ConstraintFunctionElement<Simo::AngleType1ConstraintFunction>(2, DofSet::XYZrot, _nn, _type)
{
	C0 = 0;
	axis1 = _axis1;
	axis2 = _axis2;
	offset = _offset;
	ieqtype = _ieqtype;
}

AngleType1ConstraintElement::~AngleType1ConstraintElement()
{
	if(C0) delete [] C0;
}

void
AngleType1ConstraintElement::getConstants(const CoordSet &,
                                          Eigen::Array<double,7,1>& sconst, Eigen::Array<int,1,1>& iconst,
                                          const GeomState *gs) const
{
	// gs is the current configuration for eulerian description of rotations
	// if gs == NULL the current configuration is identical to ihe undeformed configuration
	if(gs == NULL) {
		Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > C0(&AngleType1ConstraintElement::C0[0][0]);
		sconst << C0(axis1,0), C0(axis1,1), C0(axis1,2), C0(axis2,0), C0(axis2,1), C0(axis2,2), offset;
		iconst << ieqtype;
	}
	else {
		Eigen::Map<const Eigen::Matrix<double,3,3,Eigen::RowMajor> > C0(&AngleType1ConstraintElement::C0[0][0]),
				R1(&(*gs)[nn[0]].R[0][0]),
				R2(&(*gs)[nn[1]].R[0][0]);
		// rotate specified axes to their position in the specified (either reference or current) configuration
		Eigen::Matrix<double,1,3> a0 = C0.row(axis1)*R1.transpose(), b0 = C0.row(axis2)*R2.transpose();
		sconst << a0[0], a0[1], a0[2], b0[0], b0[1], b0[2], offset;
		iconst << ieqtype;
	}
}

void
AngleType1ConstraintElement::setFrame(EFrame *elemframe)
{
	C0 = new double[3][3];
	for(int i = 0; i < 3; ++i)
		for(int j = 0; j < 3; ++j)
			C0[i][j] = (*elemframe)[i][j];
}

void
AngleType1ConstraintElement::buildFrame(CoordSet& cs)
{
	if(!C0) {
		std::cerr << " *** ERROR: A reference frame needs to be specified using EFRAMES for element " << glNum+1 << ". Exiting...\n";
		exit(-1);
	}

	ConstraintFunctionElement<Simo::AngleType1ConstraintFunction>::buildFrame(cs);
}

double
AngleType1ConstraintElement::getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double)
{
	return 0.;
}

double
AngleType1ConstraintElement::getAccelerationConstraintRhs(GeomState* refState, GeomState& gState, CoordSet& cs, double t)
{
	return MpcElement::getAccelerationConstraintRhs(refState, gState, cs, t);
}

#endif
