#ifndef _SLOSHTRIANGLEFS_H_
#define _SLOSHTRIANGLEFS_H_

#include <Element.d/Element.h>
#include <cstdio>

class SloshTriangleFS: public Element {
private:
    int nn[3];
public:
    SloshTriangleFS(int*);

	int getElementType() const override { return 312; }
    Category getCategory() const override { return Category::Structural; }
    Element *clone() override;

    void renum(const int *) override;
    void renum(EleRenumMap&) override;

    FullSquareMatrix stiffness(const CoordSet& cs, double *d, int flg = 1) const override;
    FullSquareMatrix massMatrix(const CoordSet& cs, double *mel, int cmflg=1) const override;
    double getArea(const CoordSet&) const;
    bool isSloshingElement() override { return true; }

    void markDofs(DofSetArray &) const override;
    int* dofs(DofSetArray &, int *p) const override;
    int numDofs() const override;

    int numNodes() const override;
    int * nodes(int *) const override;
    int getTopNumber() const override;

    PrioInfo examine(int sub, MultiFront *) override {
        fprintf(stderr,"SloshTriangleFS.h: PrioInfo examine is commented in Dec.d/ElemFSCheck.C");
        return *(new PrioInfo);
    };
};

#endif
