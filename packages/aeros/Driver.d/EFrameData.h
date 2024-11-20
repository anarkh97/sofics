#ifndef _EFRAME_DATA_H_
#define _EFRAME_DATA_H_

#include <Element.d/Element.h>

// Element Frame Data structure
struct EFrameData {
    int elnum;
    EFrame frame;
    EFrameData *next;
};

// Node Frame Data structure
struct NFrameData {
    int elnum;
    EFrame frame;
    double origin[3];
    enum FrameType { Rectangular=0, Cylindrical, Spherical } type;
    NFrameData *next;
    template<typename Scalar>
      void transformVector3(Scalar *data);
    template<typename Scalar>
      void invTransformVector3(Scalar *data);
    template<typename Scalar>
      void transformVector6(Scalar *data);
    template<typename Scalar>
      void invTransformVector6(Scalar *data);
    template<typename Scalar>
      void transformMatrix3(Scalar *data);
    template<typename Scalar>
      void invTransformMatrix3(Scalar *data);
    template<typename Scalar>
      void transformMatrix6(Scalar *data);
    template<typename Scalar>
      void invTransformMatrix6(Scalar *data);
    template<typename Scalar>
      void transformSymMatrix3(Scalar *data);
    template<typename Scalar>
      void invTransformSymMatrix3(Scalar *data);
};

#endif
