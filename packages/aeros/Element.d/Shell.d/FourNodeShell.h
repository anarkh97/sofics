#ifndef _FOURNODESHELL_H_
#define _FOURNODESHELL_H_

#include <Element.d/SuperElement.h>

class FourNodeShell : public SuperElement
{
public:
	explicit FourNodeShell(int *nodenums);

	int getElementType() const override { return 88; }
	Element* clone() override;
	int getTopNumber() const override;
	bool isShell() const override { return true; }

	int nDecFaces() const override { return 1;}
	int getDecFace(int iFace, int *fn) override { for(int i=0; i<4; i++) fn[i] = nn[i]; return 4; }

	int getFace(int iFace, int *fn) override { return getDecFace(iFace,fn); }

	// aero functions
	void computeDisp(CoordSet &cs, State &state, const InterpPoint &ip, double *res, GeomState *gs) override;
	void getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF, double *res, GeomState *gs) override;
	bool hasRot() const override { return true; }
	PrioInfo examine(int sub, MultiFront *mf) override;

};

#endif
