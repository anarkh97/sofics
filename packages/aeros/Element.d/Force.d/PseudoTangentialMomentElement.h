#ifndef _PSEUDOTANGENTIALMOMENTELEMENT_H_
#define _PSEUDOTANGENTIALMOMENTELEMENT_H_

#include <Element.d/Function.d/ExternalWork.d/PseudoTangentialMomentWorkFunction.h>
#include <Element.d/Force.d/PotentialFunctionElement.h>

class PseudoTangentialMomentElement : public PotentialFunctionElement<Simo::PseudoTangentialMomentWorkFunction>
{
    double (*C0)[3]; // initial frame (axes stored row-wise)

public:
    static const DofSet NODALDOFS[1];
    PseudoTangentialMomentElement(int* _nn);
    ~PseudoTangentialMomentElement();

	int getElementType() const override { return 149; }
	Category getCategory() const override { return Category::Structural; }
	void setFrame(EFrame *) override;
    bool hasRot() const override { return true; }

protected:
    void getConstants(const CoordSet& cs, Eigen::Array<double,16,1>& sconst, Eigen::Array<int,0,1>&,
                      const GeomState* = nullptr,
                      double = 0) const override;
};

#endif
