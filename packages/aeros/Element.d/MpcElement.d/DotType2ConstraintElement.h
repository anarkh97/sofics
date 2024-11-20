#ifndef _DOTTYPE2CONSTRAINTELEMENT_H_
#define _DOTTYPE2CONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/DotType2ConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class DotType2ConstraintElement : public ConstraintFunctionElement<Simo::DotType2ConstraintFunction>
{
protected:
	double (*C0)[3]; // initial frame (axes stored row-wise)
	int axis;
	double d0;
	int ieqtype; // 1: c(x) <= 0, 2: c(x) >= 0 (equivalently -c(x) <= 0) ... only applies for inequality constraints (i.e. type = 1)

public:
	DotType2ConstraintElement(int*, int, int = 0, int = 1);
	~DotType2ConstraintElement() override;
	void setFrame(EFrame *) override;
	void buildFrame(CoordSet&) override;
	static const DofSet NODALDOFS[2];
	void setConstantTerm(double _d0) { d0 = _d0; }
	double getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double) override;
	double getAccelerationConstraintRhs(GeomState*, GeomState&, CoordSet&, double) override;

protected:
	void getConstants(const CoordSet & cs,
	                  Eigen::Array<double,7,1>& sconst, Eigen::Array<int,1,1>&, const GeomState *gs) const override;
};

#endif
