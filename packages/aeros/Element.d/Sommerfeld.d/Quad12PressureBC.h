#ifndef _QUAD12PRESSUREBC_H_
#define _QUAD12PRESSUREBC_H_

#include <Element.d/Sommerfeld.d/PressureElement.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Mortar.d/FaceElement.d/FaceQuad12.d/FaceQuad12.h>

typedef PressureElement<FaceQuad12, GaussLegendre2d, 3, 4, 20> Quad12PressureBC;

#endif
