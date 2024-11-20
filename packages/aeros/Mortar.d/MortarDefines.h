#ifndef _MORTAR_DEFINES_H_
#define _MORTAR_DEFINES_H_

#ifdef USE_EIGEN3
#  include <Eigen/Dense>
#  include <Eigen/Sparse>
#  include <unsupported/Eigen/AutoDiff>
#endif

#ifdef USE_SACADO
#  include <Sacado.hpp>
#  include <Sacado_trad.hpp>
#  include <Sacado_tradvec.hpp>
#endif

// note: if this is defined, make sure that COMPUTE_CENTROID_AND_LOCAL_EDGE_COORDS is also defined in Acme.d/search.d/ContactFaceFaceSearch.C
#define USE_ACME_CENTROID
//#define NO_SPECIAL_CASES

template<typename Scalar>
struct NodeTemplate {
  Scalar x,y,z;
  NodeTemplate(const Scalar& _x, const Scalar& _y, const Scalar& _z) : x(_x), y(_y), z(_z) {}
};

// helper functions to provide common interface to different AD types
template<typename Scalar>
inline Scalar GetActiveScalarValue(const Scalar& s) 
{ 
  return s; 
}

#ifdef USE_EIGEN3
template<typename Scalar, int MAX_NBDER>
inline Scalar GetActiveScalarValue(const Eigen::AutoDiffScalar<Eigen::Matrix<Scalar,MAX_NBDER,1> >& s)
{
  return s.value();
}

template<typename Scalar, int MAX_NBDER>
inline Scalar GetActiveScalarValue(const Eigen::AutoDiffScalar<Eigen::Matrix<Eigen::AutoDiffScalar<Eigen::Matrix<Scalar,MAX_NBDER,1> >,MAX_NBDER,1> >& s)
{
  return s.value().value();
}
#endif

#ifdef USE_SACADO
template<typename Scalar, int MAX_NBDER>
inline Scalar GetActiveScalarValue(const Sacado::Fad::SFad<Scalar,MAX_NBDER>& s)
{
  return s.val();
}

template<typename Scalar>
inline Scalar GetActiveScalarValue(const Sacado::RadVec::ADvar<Scalar>& s)
{
  return s.val();
}

template<typename Scalar, int MAX_NBDER>
inline Scalar GetActiveScalarValue(const Sacado::Rad::ADvar<Sacado::Fad::SFad<Scalar, MAX_NBDER> >& s)
{
  return s.val().val();
}

template<typename Scalar, int MAX_NBDER>
inline Scalar GetActiveScalarValue(const Sacado::RadVec::ADvar<Sacado::Fad::SFad<Scalar, MAX_NBDER> >& s)
{
  return s.val().val();
}

#ifdef USE_EIGEN3
template<typename Scalar, int MAX_NBDER>
inline Scalar GetActiveScalarValue(const Sacado::RadVec::ADvar<Eigen::AutoDiffScalar<Eigen::Matrix<Scalar, MAX_NBDER,1> > >& s)
{
  return s.val().value();
}
#endif
#endif

#endif
