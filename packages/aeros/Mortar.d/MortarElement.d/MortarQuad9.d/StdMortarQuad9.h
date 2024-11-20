// ---------------------------------------------------------------- 
// HB - 05/13/05
// ---------------------------------------------------------------- 
#ifndef _STDMORTARQUAD9_H_
#define _STDMORTARQUAD9_H_

#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/FaceElement.d/FaceQuad9.d/FaceQuad9.h>

class FaceElement;

class StdMortarQuad9: public MortarElement {

  public:
        // Constructors
        // ~~~~~~~~~~~~
        StdMortarQuad9();
        StdMortarQuad9(FaceElement*);
        StdMortarQuad9(double, FaceElement*);  

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
StdMortarQuad9::GetShapeFctVal(Scalar* Shape, Scalar* m)
{
  static_cast<FaceQuad9*>(GetPtrMasterFace())->GetShapeFctVal(Shape, m);
}

template<typename Scalar>
void
StdMortarQuad9::GetdShapeFct(Scalar* dShapex, Scalar* dShapey, Scalar* m)
{
  static_cast<FaceQuad9*>(GetPtrMasterFace())->GetdShapeFct(dShapex, dShapey, m);
}

template<typename Scalar>
void
StdMortarQuad9::Getd2ShapeFct(Scalar* d2Shapex, Scalar* d2Shapey, Scalar* d2Shapexy, Scalar* m)
{
  static_cast<FaceQuad9*>(GetPtrMasterFace())->Getd2ShapeFct(d2Shapex, d2Shapey, d2Shapexy, m);
}

#endif
