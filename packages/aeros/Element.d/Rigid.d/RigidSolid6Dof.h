#ifndef _RIGIDSOLID6DOF_H_
#define _RIGIDSOLID6DOF_H_

#include <Element.d/SuperElement.h>

class RigidSolid6Dof : public SuperElement
{
  public:
    RigidSolid6Dof(int, int*);

	int getElementType() const override { return 74; }
    int getTopNumber() const override;
    int numTopNodes() const override;
    bool isRigidElement() const override { return true; }
    bool hasRot() const override { return true; }
    PrioInfo examine(int sub, MultiFront*) override;
};

#endif

