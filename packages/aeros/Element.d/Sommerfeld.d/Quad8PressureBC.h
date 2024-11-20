#ifndef _QUAD8PRESSUREBC_H_
#define _QUAD8PRESSUREBC_H_

#include <Element.d/Sommerfeld.d/PressureElement.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Mortar.d/FaceElement.d/FaceQuad8.d/FaceQuad8.h>

typedef PressureElement<FaceQuad8, GaussLegendre2d, 2, 3, 18> Quad8PressureBC;

#endif
