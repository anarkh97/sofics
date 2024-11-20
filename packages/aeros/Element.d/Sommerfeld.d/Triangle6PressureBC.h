#ifndef _TRIANGLE6PRESSUREBC_H_
#define _TRIANGLE6PRESSUREBC_H_

#include <Element.d/Sommerfeld.d/PressureElement.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Mortar.d/FaceElement.d/FaceTri6.d/FaceTri6.h>

typedef PressureElement<FaceTri6, TriangleQuadratureRule<double>, 3, 6, 17> Triangle6PressureBC;

#endif
