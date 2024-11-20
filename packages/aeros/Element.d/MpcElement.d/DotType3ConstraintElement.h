#ifndef _DOTTYPE3CONSTRAINTELEMENT_H_
#define _DOTTYPE3CONSTRAINTELEMENT_H_

#include <Element.d/Function.d/Constraint.d/DotType3ConstraintFunction.h>
#include <Element.d/MpcElement.d/ConstraintFunctionElement.h>

class DotType3ConstraintElement : public ConstraintFunctionElement<Simo::DotType3ConstraintFunction>
{
protected:
	double (*C0)[3]; // initial frame (axes stored row-wise)
	int axis;
	double d0;

public:
	DotType3ConstraintElement(int*, int);
	~DotType3ConstraintElement();
	Category getCategory() const override { return Category::Structural; }
	void setFrame(EFrame *) override;
	void buildFrame(CoordSet&) override;
	void setConstantTerm(double _d0) { d0 = _d0; }
	double getVelocityConstraintRhs(GeomState*, GeomState&, CoordSet&, double) override;
	double getAccelerationConstraintRhs(GeomState*, GeomState&, CoordSet&, double) override;

protected:
	void getConstants(const CoordSet & cs, Eigen::Array<double,10,1>& sconst, Eigen::Array<int,0,1>&,
	                  const GeomState *gs = nullptr) const override;
};

#endif
