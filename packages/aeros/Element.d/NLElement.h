#ifndef _NLELEMENT_H_
#define _NLELEMENT_H_
#include <Element.d/Element.h>
#include <cstdio>

// Declaration of the Material Non linear element
class MatNLElement : public Element {
public:
    MatNLElement() {}

    Category getCategory() const override { return Category::Structural; }

    virtual void getStiffAndForce(Node *nodes, double *disp,
                                  double *state, FullSquareMatrix &kTan,
                                  double *force) {
        fprintf(stderr, "MatNLElement::getStiffAndForce is being called on an element "
                        "for which it is not defined\n");
    }
    virtual void integrate(Node *nodes, double *dispn, double *staten,
                           double *dispnp, double *statenp,
                           FullSquareMatrix &kTan,
                           double *force, double dt, double *temps) {
        int nst = numStates();
        int i;
        for(i = 0; i < nst; ++i)
            statenp[i] = staten[i];
        updateStates(nodes, dispn, dispnp, statenp, temps);
        getStiffAndForce(nodes, dispnp, statenp, kTan, force);
    }
    virtual void getInternalForce(Node *nodes, double *disp,
                                  double *state, double *force) {
        fprintf(stderr, "MatNLElement::getInternalForce is being called on an element "
                        "for which it is not defined\n");
    }
    virtual void integrate(Node *nodes, double *dispn, double *staten,
                           double *dispnp, double *statenp,
                           double *force, double dt, double *temps) {
        int nst = numStates();
        int i;
        for(i = 0; i < nst; ++i)
            statenp[i] = staten[i];
        updateStates(nodes, dispn, dispnp, statenp, temps);
        getInternalForce(nodes, dispnp, statenp, force);
    }

    virtual void updateStates(Node *nodes, double *un, double *unp, double *statenp, double *temps, double dt=0) {
        fprintf(stderr, "MatNLElement::updateStates is being called on an element "
                        "for which it is not defined\n");
    }
    virtual void getStrainTens(Node *nodes, double *dispnp, double (*result)[9], int avgnum) {
        fprintf(stderr, "MatNLElement::getStrainTens is being called on an element "
                        "for which it is not defined\n");
    }
    virtual void getVonMisesStrain(Node *nodes, double *dispnp, double *result, int avgnum) {
        fprintf(stderr, "MatNLElement::getVonMisesStrain is being called on an element "
                        "for which it is not defined\n");
    }
    virtual void getStressTens(Node *nodes, double *dispn, double *staten,
                               double *dispnp, double *statenp, double (*result)[9], int avgnum,
                               double *temps) {
        fprintf(stderr, "MatNLElement::getStressTens is being called on an element "
                        "for which it is not defined\n");
    }
    virtual void getVonMisesStress(Node *nodes, double *dispn, double *staten,
                                   double *dispnp, double *statenp, double *result, int avgnum,
                                   double *temps) {
        fprintf(stderr, "MatNLElement::getVonMisesStress is being called on an element "
                        "for which it is not defined\n");
    }
    virtual void getEquivPlasticStrain(double *statenp, double *result, int avgnum) {
        fprintf(stderr, "MatNLElement::getEquivPlasticStrain is being called on an element "
                        "for which it is not defined\n");
    }
    virtual void getBackStressTens(double *statenp, double (*result)[9], int avgnum) {
        fprintf(stderr, "MatNLElement::getBackStressTens is being called on an element "
                        "for which it is not defined\n");
    }
    virtual void getPlasticStrainTens(double *statenp, double (*result)[9], int avgnum) {
        fprintf(stderr, "MatNLElement::getPlasticStrainTens is being called on an element "
                        "for which it is not defined\n");
    }
    virtual void getDamage(double *statenp, double *result, int avgnum) {
        fprintf(stderr, "MatNLElement::getDamage is being called on an element "
                        "for which it is not defined\n");
    }

    virtual double getStrainEnergy(Node *nodes, double *dispnp, double *state, double *temps) {
        fprintf(stderr, "MatNLElement::getStrainEnergy is being called on an element "
                        "for which it is not defined\n");
        return 0.0;
    }

    virtual double getDissipatedEnergy(Node *nodes, double *state) {
        fprintf(stderr, "MatNLElement::getDissipatedEnergy is being called on an element "
                        "for which it is not defined\n");
        return 0.0;
    }

    virtual int getNumGaussPoints() { return 0; }

    virtual bool checkFailure(double *statenp) { return false; }
};

#endif
