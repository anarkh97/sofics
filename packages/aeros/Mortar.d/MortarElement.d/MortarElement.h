// ---------------------------------------------------------------- 
// HB - 06/09/03
// HB - Modified 02/28/04
// ---------------------------------------------------------------- 
#ifndef _MORTARELEMENT_H_ 
#define _MORTARELEMENT_H_ 

#include <iostream>

class FaceElement;

class MortarElement{
  private:
        double Area;

        FaceElement* MasterFace; // ptr to supporting face element 

  public:
        // Constructors
        // ~~~~~~~~~~~~
        MortarElement(); 
        MortarElement(FaceElement*);
        MortarElement(double, FaceElement*); 
        
        // Destructors
        // ~~~~~~~~~~~
        virtual ~MortarElement(); 
 
        // Initialize & clear/clean methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        void Initialize();

        // Set methods
        // ~~~~~~~~~~~
        // -> local
        void SetArea(double);
        void SetPtrMasterFace(FaceElement*);

        // Get methods
        // ~~~~~~~~~~~
        // -> local
        int nMasterFaceNodes();   
        FaceElement* GetPtrMasterFace(); 
        
        // -> interface 
        virtual int nNodes() ;
        virtual int nMortarShapeFct();
        virtual bool GetDualFlag() = 0;

        // Shape fct methods
        // ~~~~~~~~~~~~~~~~~
        // -> interface 
        virtual void GetShapeFctVal(double* Shape, double* m);
        virtual void GetdShapeFct(double* dShapex, double* dShapey, double* m);
        virtual void Getd2ShapeFct(double* d2Shapex, double* d2Shapey, double* d2Shapexy, double* m);
};

#endif
