#ifndef _BELYTSCHKOTSAYSHELL_H_
#define _BELYTSCHKOTSAYSHELL_H_

#include <Element.d/Element.h>
#include <Corotational.d/Corotator.h>
#include <Material.d/ElastoPlasticPlaneStressMaterial.h>

class GeomState;
class MultiFront;
class NLMaterial;
class ExpMat;

class BelytschkoTsayShell : virtual public Element, public Corotator
{
public:
    static double t1, t2, t3, t4, t5, t6, t7;
protected:
    // TODO most of this should belong to element property (and therefore be shared to save memory)
    int nn[4];
    int optdmg; // damage model type (0 for no damage, 1 for lematire damage model, 2 for linear softening with scaling)
    int opthgc; // hourglass control (1 for perturbation type hourglass control)
    int opttrc; // bc option (0 for pressure, 1 for traction, -1 for neither)
    int optdmp; // damping (0/1 for damping off/on)
    int optcor[2]; // warping and/or shear correction on/off
    int optprj; // rotation projection 0: off, 1: drilling only, 2: coupled drilling and rigid body modes
    double prmhgc[10]; // hourglass control parameters
    double prmdmp[10]; // damping control parameters
    int ngqpt[3]; // ngqpt[0] = gq rule for regular element
    // ngqpt[1] = gq rule for enriched element
    // ngqpt[2] = gq rule for through thickness
    int ngqpt4; // gq rule for bc or cohesive force integration
    int nndof;
    int ndime;
    int nnode;

    int mgaus[3];
    int mgqpt[2];
    double *gqpoin3;
    double *gqweigt3;
    double *evar1; // effective strain and damage
    double *evar2; // effective stress
    double *evoit1; // voight form of hourglass control stress
    double *evoit2; // voight form of local cauchy stress
    double *evoit3; // strain (local)
    ExpMat *expmat;
    bool myMat;
    ElastoPlasticPlaneStressMaterial **mat;
    PressureBCond *pbc;

public:
    BelytschkoTsayShell(int*);
    ~BelytschkoTsayShell();

	int getElementType() const override { return 16; }
    Category getCategory() const override { return Category::Structural; }
    void setProp(StructProp *p, bool _myProp) override;
    void setMaterial(NLMaterial *) override;
    void setPressure(PressureBCond *_pbc) override;
    PressureBCond* getPressure() override;
    Element *clone() override;

    void renum(const int *) override;
    void renum(EleRenumMap&) override;

    FullSquareMatrix stiffness(const CoordSet&, double* d, int flg) const override;
    FullSquareMatrix massMatrix(const CoordSet&, double* mel, int cmflg) const override;
    double getMass(const CoordSet& cs) const override;
    double getMassThicknessSensitivity(CoordSet& cs) override;

    void getGravityForce(CoordSet&, double* gravity, Vector&, int gravflg,
                         GeomState *gs) override;
    void getGravityForceThicknessSensitivity(CoordSet&, double* gravity, Vector&, int gravflg,
                                             GeomState *gs) override;
    void getVonMises(Vector& stress, Vector& weight, CoordSet& cs,
                     Vector& elDisp,  int strInd, int surface,
                     double *ndTemps, double ylayer,
                     double zlayer, int avgnum) override;

    void markDofs(DofSetArray&) const override;
    int* dofs(DofSetArray&, int* p) const override;
    int numDofs() const override;

    int numNodes() const override;
    int *nodes(int *) const override;
    Corotator *getCorotator(CoordSet&, double*, int , int) override;
    void getStiffAndForce(GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double) override;
    void getInternalForce(GeomState&, CoordSet&, FullSquareMatrix&, double*, double, double) override;

    // (AN): Local implementation for calculating dissipated energy
    double getDissipatedEnergy(GeomState &, CoordSet &) override;

    bool checkElementDeletion(GeomState &) override;
    void extractDeformations(GeomState &geomState, CoordSet &cs, double *vld, int &nlflag) override;

    void computeDisp(CoordSet&, State&, const InterpPoint&, double*,
                     GeomState*) override;
    void getFlLoad(CoordSet&, const InterpPoint&, double*, double *,
                   GeomState*) override;

    int getTopNumber() const override;
    void computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState* gs, int cflg, double t) override;

    void getThermalForce(CoordSet& cs, Vector& ndTemps, Vector &elThermalForce,
                         int glfag, GeomState* gs) override;

    bool isShell() const override { return true; }

    int getMassType() const override { return 0; } // lumped only

    bool hasRot() const override { return true; }

    PrioInfo examine(int sub, MultiFront* mf) override;

    void writeHistory(int fn) const override;
    void readHistory(int fn) override;

    double computeStabilityTimeStep(FullSquareMatrix &K, FullSquareMatrix &M, CoordSet &cs,
                                    GeomState *gs, double stable_tol, int stable_maxit) override;

private:
    void Elefintbt1(double delt, double *ecord, double *edisp, double *evelo,
                    double trac[3], double tmftval, double *efint);

};
#endif

