#ifndef _FOLLOWERFORCEELEMENT_H_
#define _FOLLOWERFORCEELEMENT_H_

#include <Element.d/Function.d/ExternalForce.d/FollowerForceFunction.h>
#include <Element.d/Force.d/ForceFunctionElement.h>

class FollowerForceElement : public ForceFunctionElement<Simo::FollowerForceFunction>
{
	Eigen::Matrix3d *C0; // initial frame (axes stored row-wise)

public:
	static const DofSet NODALINPUTDOFS[1];
	static const DofSet NODALOUTPUTDOFS[1];
	explicit FollowerForceElement(int* _nn);
	~FollowerForceElement() override;

	int getElementType() const override { return 140; }
	Category getCategory() const override { return Category::Structural; }

	void setFrame(EFrame *) override;
	bool hasRot() const override { return true; }

protected:
	void getConstants(const CoordSet& cs, Eigen::Array<double,12,1>& sconst, Eigen::Array<int,0,1>&,
					  const GeomState* = nullptr, double = 0) const override;
};

#endif
