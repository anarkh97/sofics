#ifndef _PENTAHEDRAL_H_
#define _PENTAHEDRAL_H_

#include <Element.d/Element.h>

class Pentahedral: public Element
{
    int nn[6];
    double *cCoefs;
    double *cFrame;
    NLMaterial *mat;

public:
    explicit Pentahedral(int*);
    ~Pentahedral() override;

	int getElementType() const override { return 24; }
    Category getCategory() const override { return Category::Structural; }
    Element *clone() override;

    void renum(const int *) override;
    void renum(EleRenumMap&) override;

    FullSquareMatrix stiffness(const CoordSet&, double *kel, int flg) const override;
    FullSquareMatrix massMatrix(const CoordSet&,double *mel, int cmflg) const override;
    double getMass(const CoordSet& cs) const override;

    void getGravityForce(CoordSet&, double *gravity, Vector&, int gravflg, GeomState *gs) override;
    void getThermalForce(CoordSet &cs, Vector &ndTemps, Vector &force, int glflag, GeomState *gs) override;

    void getVonMises(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
                     int surface, double *ndTemps, double ylayer, double zlayer, int avgnum) override;

    void getAllStress(FullM &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
                      int surface, double *ndTemps) override;

    void markDofs(DofSetArray &) const override;
    int* dofs(DofSetArray &, int *p) const override;
    int numDofs() const override;

    int numNodes() const override;
    int* nodes(int *) const override;

    int getTopNumber() const override;

    PrioInfo examine(int sub, MultiFront *) override;
    int nDecFaces() const override { return 5; }
    int getDecFace(int iFace, int *fn) override;

    int getFace(int iFace, int *fn) override { return getDecFace(iFace, fn); }

    void setCompositeData(int _type, int nlays, double *lData, double *coefs, double *frame) override
    { cCoefs = coefs; cFrame = frame; }

    double* setCompositeData2(int, int, double*, double*, CoordSet&, double) override
    { fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
                     "              for Pentahedral el.\n"); return (double *) 0;
    }
    void getCFrame(CoordSet &cs, double cFrame[3][3]) const override;

    void getVonMisesAniso(Vector &stress, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd,
                          int surface, double *ndTemps, double ylayer, double zlayer, int avgnum);

    void getAllStressAniso(FullM &stress, Vector &weight, CoordSet &cs,
                           Vector &elDisp, int strInd, int surface, double *ndTemps);

    void setMaterial(NLMaterial *) override;
    int numStates() override;
    void initStates(double *st) override;
    Corotator *getCorotator(CoordSet &cs, double *kel, int, int) override;
};

#endif
