// ---------------------------------------------------------------- 
// HB - 05/13/05
// ---------------------------------------------------------------- 
#ifndef _STDMORTARTRI6_H_
#define _STDMORTARTRI6_H_

#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/FaceElement.d/FaceTri6.d/FaceTri6.h>

class FaceElement;

class StdMortarTri6: public MortarElement {

  public:
        // Constructors
        // ~~~~~~~~~~~~
        StdMortarTri6();
        StdMortarTri6(FaceElement*);
        StdMortarTri6(double, FaceElement*);  
 
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
StdMortarTri6::GetShapeFctVal(Scalar* Shape, Scalar* m)
{
  static_cast<FaceTri6*>(GetPtrMasterFace())->GetShapeFctVal(Shape, m);
}

template<typename Scalar>
void
StdMortarTri6::GetdShapeFct(Scalar* dShapex, Scalar* dShapey, Scalar* m)
{
  static_cast<FaceTri6*>(GetPtrMasterFace())->GetdShapeFct(dShapex, dShapey, m);
}

template<typename Scalar>
void
StdMortarTri6::Getd2ShapeFct(Scalar* d2Shapex, Scalar* d2Shapey, Scalar* d2Shapexy, Scalar* m)
{
  static_cast<FaceTri6*>(GetPtrMasterFace())->Getd2ShapeFct(d2Shapex, d2Shapey, d2Shapexy, m);
}

#endif
