#ifndef _RIGIDSPRING_H_
#define _RIGIDSPRING_H_

#include <Element.d/SuperElement.h>

class RigidSpring : public SuperElement
{
  public:
    RigidSpring(int*);

	int getElementType() const override { return 67; }
    int getTopNumber() const override;
    bool isRigidElement() const override { return true; }
    bool isSpring() const override { return true; }
    bool hasRot() const override { return true; }
    PrioInfo examine(int sub, MultiFront*) override;
};

#endif
