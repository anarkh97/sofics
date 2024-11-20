#ifndef _QUADPRESSUREBC_H_
#define _QUADPRESSUREBC_H_

#include <Element.d/Sommerfeld.d/PressureElement.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Mortar.d/FaceElement.d/FaceQuad4.d/FaceQuad4.h>

typedef PressureElement<FaceQuad4, GaussLegendre2d, 2, 3, 16> QuadPressureBC;

#endif
