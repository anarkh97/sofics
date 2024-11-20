#ifndef _COMPO4NODESHELL_H_
#define _COMPO4NODESHELL_H_

#include <Element.d/SuperElement.h>

class Compo4NodeShell : public SuperElement
{
  public:
    Compo4NodeShell(int *nodenums);

	int getElementType() const override { return 2020; }
    Element* clone() override;
    int getTopNumber() const override;
    bool isShell() const override { return true; }

    // aero functions
    void computeDisp(CoordSet &cs, State &state, const InterpPoint &ip, double *res, GeomState *gs) override;
    void getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF, double *res, GeomState *gs) override;
    bool hasRot() const override { return true; }

    // Routines for the decomposer
    PrioInfo examine(int sub, MultiFront *) override;
    int nDecFaces() const override { return 1; }
    int getDecFace(int iFace, int *fn) override { for(int i=0; i<4; i++) fn[i] = nn[i]; return 4; }

    int getFace(int iFace, int *fn) override { return getDecFace(iFace,fn); }
};

#endif
