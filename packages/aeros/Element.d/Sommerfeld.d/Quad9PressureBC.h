#ifndef _QUAD9PRESSUREBC_H_
#define _QUAD9PRESSUREBC_H_

#include <Element.d/Sommerfeld.d/PressureElement.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Mortar.d/FaceElement.d/FaceQuad9.d/FaceQuad9.h>

typedef PressureElement<FaceQuad9, GaussLegendre2d, 2, 3, 19> Quad9PressureBC;

#endif
