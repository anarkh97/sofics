#ifndef _FOURNODEMEMBRANE_H_
#define _FOURNODEMEMBRANE_H_

#include <Element.d/SuperElement.h>

class FourNodeMembrane : public SuperElement 
{
  public:
    FourNodeMembrane(int *nodenums);

	int getElementType() const override { return 87; }
    int getTopNumber() const override;
    PrioInfo examine(int sub, MultiFront *mf) override;
};

#endif
