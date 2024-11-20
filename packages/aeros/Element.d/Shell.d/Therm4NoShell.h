#ifndef _THERM4NOSHELL_H_
#define _THERM4NOSHELL_H_

#include <Element.d/SuperElement.h>
#include <Element.d/Element.h>

class Therm4NoShell : public SuperElement
{
public:
	explicit Therm4NoShell(int *nodenums);

	int getElementType() const override { return 4646; }
	Element* clone() override;
	int getTopNumber() const override;
	PrioInfo examine(int sub, MultiFront *) override;
	bool hasRot() const override {return true;}
	// aero functions
	void computeTemp(CoordSet &cs, State &state, double[2], double *res) override;
	void getFlFlux(double[2], double *flF, double *res) override;
};

#endif
