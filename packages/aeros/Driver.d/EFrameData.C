#ifndef _EFRAME_DATA_C_
#define _EFRAME_DATA_C_

#include <Driver.d/EFrameData.h>
#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif

template<typename Scalar>
inline void
NFrameData::transformVector3(Scalar *data)
{
  // transform 3x1 vector from basic to DOF_FRM coordinates
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&frame[0][0]);
  Eigen::Map<Eigen::Matrix<Scalar,3,1> > v(data);

  v = T*v;
#endif
}

template<typename Scalar>
inline void
NFrameData::invTransformVector3(Scalar *data)
{
  // transform 3x1 vector from DOF_FRM to basic coordinates
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&frame[0][0]);
  Eigen::Map<Eigen::Matrix<Scalar,3,1> > v(data);

  v = T.transpose()*v;
#endif
}

template<typename Scalar>
inline void
NFrameData::transformVector6(Scalar *data)
{
  // transform 6x1 vector from basic to DOF_FRM coordinates
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&frame[0][0]);
  Eigen::Map<Eigen::Matrix<Scalar,3,2,Eigen::ColMajor> > v(data);

  v = T*v;
#endif
}

template<typename Scalar>
inline void
NFrameData::invTransformVector6(Scalar *data)
{
  // transform 6x1 vector from DOF_FRM to basic coordinates
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&frame[0][0]);
  Eigen::Map<Eigen::Matrix<Scalar,3,2,Eigen::ColMajor> > v(data);

  v = T.transpose()*v;
#endif
}

template<typename Scalar>
inline void
NFrameData::transformMatrix3(Scalar *data)
{
  // transform 3x3 matrix from basic to DOF_FRM coordinates
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&frame[0][0]);
  Eigen::Map<Eigen::Matrix<Scalar,3,3,Eigen::RowMajor> > M(data);

  M = T*M*T.transpose();
#endif
}

template<typename Scalar>
inline void
NFrameData::invTransformMatrix3(Scalar *data)
{
#ifdef USE_EIGEN3
  // transform 3x3 matrix from DOF_FRM to basic coordinate
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&frame[0][0]);
  Eigen::Map<Eigen::Matrix<Scalar,3,3,Eigen::RowMajor> > M(data);

  M = T.transpose()*M*T;
#endif
}

template<typename Scalar>
inline void
NFrameData::transformMatrix6(Scalar *data)
{
  // transform 6x6 matrix from basic to DOF_FRM coordinates
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&frame[0][0]);
  Eigen::Map<Eigen::Matrix<Scalar,6,6,Eigen::RowMajor> > M(data);

  Eigen::Matrix<double,6,6> T6;
  T6 << T, Eigen::Matrix3d::Zero(),
        Eigen::Matrix3d::Zero(), T;

  M = T6*M*T6.transpose();
#endif
}

template<typename Scalar>
inline void
NFrameData::invTransformMatrix6(Scalar *data)
{
#ifdef USE_EIGEN3
  // transform 6x6 matrix from DOF_FRM to basic coordinate
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&frame[0][0]);
  Eigen::Map<Eigen::Matrix<Scalar,6,6,Eigen::RowMajor> > M(data);

  Eigen::Matrix<double,6,6> T6;
  T6 << T, Eigen::Matrix3d::Zero(),
        Eigen::Matrix3d::Zero(), T;

  M = T6.transpose()*M*T6;
#endif
}

template<typename Scalar>
inline void
NFrameData::transformSymMatrix3(Scalar *data)
{
  // transform 3x3 symmetric matrix from basic to DOF_FRM coordinates
  // data = { xx, yy, zz, xy, yz, xz }
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&frame[0][0]);

  Eigen::Matrix<Scalar,3,3> M;
  M << data[0], data[3], data[5],
       data[3], data[1], data[4],
       data[5], data[4], data[2];

  M = T*M*T.transpose();

  data[0] = M(0,0);
  data[1] = M(1,1);
  data[2] = M(2,2);
  data[3] = M(0,1);
  data[4] = M(1,2);
  data[5] = M(0,2);
#endif
}

template<typename Scalar>
inline void
NFrameData::invTransformSymMatrix3(Scalar *data)
{
#ifdef USE_EIGEN3
  // transform 3x3 symmetric matrix from DOF_FRM to basic coordinate
  // data = { xx, yy, zz, xy, yz, xz }

  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > T(&frame[0][0]);

  Eigen::Matrix<Scalar,3,3> M;
  M << data[0], data[3], data[5],
       data[3], data[1], data[4],
       data[5], data[4], data[2];

  M = T.transpose()*M*T;

  data[0] = M(0,0);
  data[1] = M(1,1);
  data[2] = M(2,2);
  data[3] = M(0,1);
  data[4] = M(1,2);
  data[5] = M(0,2);
#endif
}

#endif
