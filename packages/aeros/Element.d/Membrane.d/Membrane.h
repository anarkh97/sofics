#ifndef _MEMBRANE_H_
#define _MEMBRANE_H_

#include <Element.d/Element.h>

class Membrane : public Element
{
private:
    int nn[3];
public:
    Membrane(int*);

	int getElementType() const override { return 19; }
    Category getCategory() const override { return Category::Structural; }
    Element *clone() override;

    void renum(const int *) override;
    void renum(EleRenumMap&) override;

    FullSquareMatrix stiffness(const CoordSet&, double *d, int flg=1) const override;
    FullSquareMatrix massMatrix(const CoordSet&, double *mel, int cmflg) const override;
    double getMass(const CoordSet& cs) const override;
    void getGravityForce(CoordSet&, double *gravity, Vector&, int gravflg,
                         GeomState *gs) override;
    void getVonMises(Vector& stress, Vector& weight, CoordSet& cs,
                     Vector& elDisp, int strInd, int surface=0,
                     double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum=0) override;
    void getAllStress(FullM& stress, Vector& weight, CoordSet& cs,
                      Vector& elDisp, int strInd, int surface=0,
                      double *ndTemps=0) override;

    void markDofs(DofSetArray &) const override;
    int* dofs(DofSetArray &, int *p) const override;
    int numDofs() const override;

    int numNodes() const override;
    int* nodes(int *) const override;
    int getTopNumber() const override;
    Corotator* getCorotator(CoordSet &cs, double *kel, int fitAlg, int) override;

    // Routines for the decomposer
    PrioInfo examine(int sub, MultiFront *) override;

    // Miscellaneous
    bool hasRot() const override { return true; }
    int  getMassType() { return 0; } // lumped only

#ifdef USE_EIGEN3
    // NEW STRUCTOPT
    double getMassThicknessSensitivity(CoordSet& cs) override;
    void getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet&, double *gravityAcceleration,
                                             int senMethod = 1);
    void getGravityForceThicknessSensitivity(CoordSet&, double *gravity, int senMethod, Vector&, int gravflg,
                                             GeomState *gs = 0);
    void getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration, int senMethod,
                                                   GenFullM<double> &dGfdx, int gravflg, GeomState *gs = 0);
    void getStiffnessThicknessSensitivity(CoordSet& cs, FullSquareMatrix &dStiffdThick, int flg = 1,
                                          int senMethod = 0);
    void getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs, int senMethod = 0);
    void getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs,
                                               Vector &elDisp, int strInd, int surface, int senMethod = 1,
                                               double *ndTemps = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);
    void getVonMisesThicknessSensitivity(Vector &dStdThick, Vector &weight, CoordSet &cs, Vector &elDisp,
                                         int strInd, int surface, int senMethod = 1, double *ndTemps = 0, int avgnum = 1,
                                         double ylayer = 0, double zlayer = 0);
    void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *dDispDisp,
                                            CoordSet &cs, Vector &elDisp, int strInd, int surface, int senMethod = 1,
                                            double *ndTemps = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);
#endif
};
#endif
