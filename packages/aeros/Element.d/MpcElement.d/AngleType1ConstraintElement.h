#ifndef _ANGLETYPE1CONSTRAINTELEMENT_H_
#define _ANGLETYPE1CONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/AngleType1ConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class AngleType1ConstraintElement : public ConstraintFunctionElement<Simo::AngleType1ConstraintFunction>
{
protected:
	double (*C0)[3]; // initial frame (axes stored row-wise)
	int axis1, axis2;
	double offset;
	int ieqtype; // 1: c(x) <= 0, 2: c(x) >= 0 (equivalently -c(x) <= 0) ... only applies for inequality constraints (i.e. type = 1)

public:
	AngleType1ConstraintElement(int*, int, int, double = M_PI/2, int = 0, int = 1);
	~AngleType1ConstraintElement() override;
	void buildFrame(CoordSet&) override;
	void setFrame(EFrame *) override;
	double getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double) override;
	double getAccelerationConstraintRhs(GeomState*, GeomState&, CoordSet&, double) override;

protected:
	void getConstants(const CoordSet &,
	                  Eigen::Array<double,7,1>& sconst, Eigen::Array<int,1,1>&,
	                  const GeomState *gs = nullptr) const override;
};

#endif
