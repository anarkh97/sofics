// ---------------------------------------------------------------- 
// HB - 08/25/03
// ---------------------------------------------------------------- 
#ifndef _STDMORTARTRI3_H_
#define _STDMORTARTRI3_H_

#include <Mortar.d/MortarElement.d/MortarElement.h>
#include <Mortar.d/FaceElement.d/FaceTri3.d/FaceTri3.h>

class FaceElement;

class StdMortarTri3: public MortarElement {

  public:
        // Constructors
        // ~~~~~~~~~~~~
        StdMortarTri3();
        StdMortarTri3(FaceElement*);
        StdMortarTri3(double, FaceElement*);  

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
StdMortarTri3::GetShapeFctVal(Scalar* Shape, Scalar* m)
{
  static_cast<FaceTri3*>(GetPtrMasterFace())->GetShapeFctVal(Shape, m);
}

template<typename Scalar>
void
StdMortarTri3::GetdShapeFct(Scalar* dShapex, Scalar* dShapey, Scalar* m)
{
  static_cast<FaceTri3*>(GetPtrMasterFace())->GetdShapeFct(dShapex, dShapey, m);
}

template<typename Scalar>
void
StdMortarTri3::Getd2ShapeFct(Scalar* d2Shapex, Scalar* d2Shapey, Scalar* d2Shapexy, Scalar* m)
{
  static_cast<FaceTri3*>(GetPtrMasterFace())->Getd2ShapeFct(d2Shapex, d2Shapey, d2Shapexy, m);
}

#endif
