#ifndef _TRIANGLE10PRESSUREBC_H_
#define _TRIANGLE10PRESSUREBC_H_

#include <Element.d/Sommerfeld.d/PressureElement.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Mortar.d/FaceElement.d/FaceTri10.d/FaceTri10.h>

typedef PressureElement<FaceTri10, TriangleQuadratureRule<double>, 6, 7, 21> Triangle10PressureBC;

#endif
