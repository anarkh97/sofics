// ---------------------------------------------------------------- 
// HB - 05/13/05
// ---------------------------------------------------------------- 
#ifndef _STDMORTARQUAD8_H_
#define _STDMORTARQUAD8_H_

#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/FaceElement.d/FaceQuad8.d/FaceQuad8.h>

class FaceElement;

class StdMortarQuad8: public MortarElement {

  public:
        // Constructors
        // ~~~~~~~~~~~~
        StdMortarQuad8();
        StdMortarQuad8(FaceElement*);
        StdMortarQuad8(double, FaceElement*);  
 
        // Get methods
        // ~~~~~~~~~~~
        // -> implementation of virtual methods
        int nNodes() override;
        int nMortarShapeFct() override;
        bool GetDualFlag() override { return false; }

        // Shape fct methods
        // ~~~~~~~~~~~~~~~~~
        // -> local methods
        template<typename Scalar>
          void GetShapeFctVal(Scalar* Shape, Scalar* m);
        template<typename Scalar>
          void GetdShapeFct(Scalar* dShapex, Scalar* dShapey, Scalar* m);
        template<typename Scalar>
          void Getd2ShapeFct(Scalar* d2Shapex, Scalar* d2Shapey, Scalar* d2Shapexy, Scalar* m);

        // -> implementation of virtual methods
        void GetShapeFctVal(double* Shape, double* m) override;
        void GetdShapeFct(double* dShapex, double* dShapey, double* m) override;
        void Getd2ShapeFct(double* d2Shapex, double* d2Shapey, double* d2Shapexy, double* m) override;
};

template<typename Scalar>
void
StdMortarQuad8::GetShapeFctVal(Scalar* Shape, Scalar* m)
{
  static_cast<FaceQuad8*>(GetPtrMasterFace())->GetShapeFctVal(Shape, m);
}

template<typename Scalar>
void
StdMortarQuad8::GetdShapeFct(Scalar* dShapex, Scalar* dShapey, Scalar* m)
{
  static_cast<FaceQuad8*>(GetPtrMasterFace())->GetdShapeFct(dShapex, dShapey, m);
}

template<typename Scalar>
void
StdMortarQuad8::Getd2ShapeFct(Scalar* d2Shapex, Scalar* d2Shapey, Scalar* d2Shapexy, Scalar* m)
{
  static_cast<FaceQuad8*>(GetPtrMasterFace())->Getd2ShapeFct(d2Shapex, d2Shapey, d2Shapexy, m);
}

#endif
