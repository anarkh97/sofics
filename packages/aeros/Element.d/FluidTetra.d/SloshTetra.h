#ifndef _SLOSHTETRA_H_
#define _SLOSHTETRA_H_

#include <Element.d/Element.h>
#include <cstdio>

class SloshTetra: public Element {
private:
    int nn[4];
public:
    explicit SloshTetra(int*);

	int getElementType() const override { return 311; }
    Category getCategory() const override { return Category::Fluid; }
    Element *clone() override;

    void renum(const int *) override;
    void renum(EleRenumMap&) override;

    FullSquareMatrix stiffness(const CoordSet&, double *kel, int flg) const override;
    FullSquareMatrix massMatrix(const CoordSet&,double *mel, int cmflg) const override;
    bool isSloshingElement() override { return true; }
    double getMass(const CoordSet& cs) const override;
    void markDofs(DofSetArray &) const override;
    int* dofs(DofSetArray &, int *p) const override;
    int numDofs() const override;

    int             numNodes() const override;
    int * nodes(int *) const override;
    int getTopNumber() const override;
    PrioInfo examine(int sub, MultiFront *) override {
        fprintf(stderr,"SloshTetra.h: PrioInfo examine is commented in Dec.d/ElemMFCheck.C\n");
        return *(new PrioInfo);
    };

private:

    double volume(double dOmega) const {return dOmega/6.0;}
    double computeMetrics(const CoordSet &, double gN[4][3]) const;
    void buildTetraStiff(double m[4][4], double gN[4][3], double dOmega) const;

};

#endif
