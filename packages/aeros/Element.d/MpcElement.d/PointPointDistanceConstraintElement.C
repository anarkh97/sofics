#ifdef USE_EIGEN3
#include <Element.d/MpcElement.d/PointPointDistanceConstraintElement.h>

PointPointDistanceConstraintElement::PointPointDistanceConstraintElement(int* _nn)
		: ConstraintFunctionElement<Simo::PointPointDistanceConstraintFunction>(1, DofSet::XYZdisp, _nn, 0)
{
}

void
PointPointDistanceConstraintElement::setFrame(EFrame *elemframe)
{
	for(int i = 0; i < 3; ++i) {
		x1[i] = (*elemframe)[0][i];
	}
}

void
PointPointDistanceConstraintElement::getConstants(const CoordSet & cs,
                                                  Eigen::Array<double,11,1>& sconst, Eigen::Array<int,1,1>& iconst,
                                                  const GeomState*) const
{
	// note: StructProps::relop = -1 --> f(x) <= 0
	//                          = 0  --> f(x) = 0
	//                          = +1 --> f(x) >= 0 in this case ConstraintFunction::operator() returns -f(x)
	//                                             and we enforce -f(x) <= 0
	this->type = (prop) ? int(prop->relop != 0) : 0;

	sconst << cs[nn[0]]->x, cs[nn[0]]->y, cs[nn[0]]->z, x1[0], x1[1], x1[2],
			(prop) ? prop->amplitude : 0, (prop) ? prop->omega : 0, (prop) ? prop->phase : 0,
			(prop) ? prop->B : 1, (prop) ? prop->C : 0;
	iconst << ((prop) ? int(prop->relop == 1) : 0);
}
#endif
