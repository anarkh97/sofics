#ifndef _TRIANGLEPRESSUREBC_H_
#define _TRIANGLEPRESSUREBC_H_

#include <Element.d/Sommerfeld.d/PressureElement.h>
#include <Element.d/Function.d/QuadratureRule.h>
#include <Mortar.d/FaceElement.d/FaceTri3.d/FaceTri3.h>

typedef PressureElement<FaceTri3, TriangleQuadratureRule<double>, 1, 3, 15> TrianglePressureBC;

#endif
