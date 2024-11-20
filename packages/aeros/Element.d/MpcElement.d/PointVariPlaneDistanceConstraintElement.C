#ifdef USE_EIGEN3
#include <Element.d/MpcElement.d/PointVariPlaneDistanceConstraintElement.h>

PointVariPlaneDistanceConstraintElement::PointVariPlaneDistanceConstraintElement(int* _nn)
		: ConstraintFunctionElement<Simo::PointVariPlaneDistanceConstraintFunction>(4, DofSet::XYZdisp, _nn, 0)
{
}

void
PointVariPlaneDistanceConstraintElement::getConstants(const CoordSet & cs,
                                                      Eigen::Array<double,17,1>& sconst, Eigen::Array<int,1,1>& iconst,
                                                      const GeomState*) const
{
	// note: StructProps::relop = -1 --> f(x) <= 0
	//                          = 0  --> f(x) = 0
	//                          = +1 --> f(x) >= 0 in this case ConstraintFunction::operator() returns -f(x)
	//                                             and we enforce -f(x) <= 0
	this->type = (prop) ? int(prop->relop != 0) : 0;

	sconst << cs[nn[0]]->x, cs[nn[0]]->y, cs[nn[0]]->z, cs[nn[1]]->x, cs[nn[1]]->y, cs[nn[1]]->z,
			cs[nn[2]]->x, cs[nn[2]]->y, cs[nn[2]]->z, cs[nn[3]]->x, cs[nn[3]]->y, cs[nn[3]]->z,
			(prop) ? prop->amplitude : 0, (prop) ? prop->omega : 0, (prop) ? prop->phase : 0,
			(prop) ? prop->B : 1, (prop) ? prop->C : 0;
	iconst << ((prop) ? int(prop->relop == 1) : 0);
}
#endif
