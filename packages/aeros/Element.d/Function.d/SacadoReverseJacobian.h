#ifndef _SACADOREVERSEJACOBIAN_H_
#define _SACADOREVERSEJACOBIAN_H_

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>

#ifdef USE_SACADO
#include <Sacado.hpp>
#include <Sacado_trad.hpp>
#include <Sacado_tradvec.hpp>

namespace Eigen { namespace internal {

template<typename Scalar, int NumDeriv>
struct scalar_product_traits<Sacado::Rad::ADvar<Sacado::Fad::SFad<Scalar, NumDeriv> >, Scalar>
{
  enum { Defined = 1 };
  typedef Sacado::Rad::ADvar<Sacado::Fad::SFad<Scalar, NumDeriv> > ReturnType;
};

template<typename Scalar, int NumDeriv>
struct scalar_product_traits<Scalar, Sacado::Rad::ADvar<Sacado::Fad::SFad<Scalar, NumDeriv> > >
{
  enum { Defined = 1 };
  typedef Sacado::Rad::ADvar<Sacado::Fad::SFad<Scalar, NumDeriv> > ReturnType;
};

template<typename Scalar, int NumDeriv>
struct scalar_product_traits<Sacado::Rad::ADvar<Eigen::AutoDiffScalar<Eigen::Matrix<Scalar,NumDeriv,1> > >, Scalar>
{
  enum { Defined = 1 };
  typedef Sacado::Rad::ADvar<Eigen::AutoDiffScalar<Eigen::Matrix<Scalar,NumDeriv,1> > > ReturnType;
};

template<typename Scalar, int NumDeriv>
struct scalar_product_traits<Scalar, Sacado::Rad::ADvar<Eigen::AutoDiffScalar<Eigen::Matrix<Scalar,NumDeriv,1> > > >
{
  enum { Defined = 1 };
  typedef Sacado::Rad::ADvar<Eigen::AutoDiffScalar<Eigen::Matrix<Scalar,NumDeriv,1> > > ReturnType;
};

template<typename Scalar>
struct scalar_product_traits<Sacado::RadVec::ADvar<Scalar>, Scalar>
{ 
  enum { Defined = 1 };
  typedef Sacado::RadVec::ADvar<Scalar> ReturnType;
};

template<typename Scalar>
struct scalar_product_traits<Scalar, Sacado::RadVec::ADvar<Scalar> >
{ 
  enum { Defined = 1 };
  typedef Sacado::RadVec::ADvar<Scalar> ReturnType;
};

template<typename Scalar, int NumDeriv>
struct scalar_product_traits<Sacado::RadVec::ADvar<Eigen::AutoDiffScalar<Eigen::Matrix<Scalar, NumDeriv, 1> > >, Scalar>
{ 
  enum { Defined = 1 };
  typedef Sacado::RadVec::ADvar<Eigen::AutoDiffScalar<Eigen::Matrix<Scalar, NumDeriv, 1> > > ReturnType;
};

template<typename Scalar, int NumDeriv>
struct scalar_product_traits<Scalar, Sacado::RadVec::ADvar<Eigen::AutoDiffScalar<Eigen::Matrix<Scalar, NumDeriv, 1> > > >
{
  enum { Defined = 1 };
  typedef Sacado::RadVec::ADvar<Eigen::AutoDiffScalar<Eigen::Matrix<Scalar, NumDeriv, 1> > > ReturnType;
};

template<typename Scalar, int NumDeriv>
struct scalar_product_traits<Eigen::AutoDiffScalar<Eigen::Matrix<Sacado::RadVec::ADvar<Scalar>, NumDeriv, 1> >, Scalar>
{
  enum { Defined = 1 };
  typedef Eigen::AutoDiffScalar<Eigen::Matrix<Sacado::RadVec::ADvar<Scalar>, NumDeriv, 1> > ReturnType;
};

template<typename Scalar, int NumDeriv>
struct scalar_product_traits<Scalar, Eigen::AutoDiffScalar<Eigen::Matrix<Sacado::RadVec::ADvar<Scalar>, NumDeriv, 1> > >
{
  enum { Defined = 1 };
  typedef Eigen::AutoDiffScalar<Eigen::Matrix<Sacado::RadVec::ADvar<Scalar>, NumDeriv, 1> > ReturnType;
};

} }
#endif

template<typename Functor> class SacadoReverseJacobian : public Functor
{
#ifdef USE_SACADO
  typedef Sacado::RadVec::ADvar<typename Functor::Scalar> ActiveScalar;
#endif
public:
  SacadoReverseJacobian() : Functor() {}
  SacadoReverseJacobian(const Functor& f) : Functor(f) {}

  // forward constructors
  template<typename T0>
  SacadoReverseJacobian(const T0& a0) : Functor(a0) {}
  template<typename T0, typename T1>
  SacadoReverseJacobian(const T0& a0, const T1& a1) : Functor(a0, a1) {}
  template<typename T0, typename T1, typename T2>
  SacadoReverseJacobian(const T0& a0, const T1& a1, const T2& a2) : Functor(a0, a1, a2) {}

  typedef typename Functor::Scalar Scalar;
  typedef typename Functor::InputType InputType;
  typedef typename Functor::JacobianType ValueType;
#ifdef USE_SACADO
  typedef Eigen::Matrix<ActiveScalar, InputType::SizeAtCompileTime, 1> ActiveInput;
  typedef Eigen::Matrix<ActiveScalar, Functor::ValueType::SizeAtCompileTime, 1> ActiveValue;
#endif

  template<typename T>
  int operator() (const Eigen::Matrix<T,InputType::SizeAtCompileTime,1>& x,
                  Eigen::Matrix<T,Functor::ValueType::SizeAtCompileTime,InputType::SizeAtCompileTime>& jac) const
  {
#ifdef USE_SACADO
// note: the 2nd condition is a hack to get code to build with buggy version of icpc 12
#if defined(_OPENMP) && (!defined(__INTEL_COMPILER) || __INTEL_COMPILER < 1200 || __INTEL_COMPILER > 1210)
    #pragma omp critical
    {
#endif
    ActiveInput ax = x.template cast<ActiveScalar>();
    ActiveValue av;

    Functor::operator()(ax, av);

    for (int i=0; i<jac.rows(); i++)
    {
      Sacado::RadVec::ADvar<typename Functor::Scalar>::Outvar_Gradcomp(av[i]);
      for (int j=0; j<jac.cols(); j++)
        jac.coeffRef(i,j) = ax[j].adj();
    }

    Sacado::RadVec::ADvar<typename Functor::Scalar>::aval_reset();
#if defined(_OPENMP) && (!defined(__INTEL_COMPILER) || __INTEL_COMPILER < 1200 || __INTEL_COMPILER > 1210)
    }
#endif
#else
    std::cerr << "Error: USE_SACADO is not defined in SacadoReverseJacobian::operator()\n";
    exit(-1);
#endif
    return 0;
  }
};

template<typename Functor> class SacadoHessian : public Functor
{
#ifdef USE_SACADO
  typedef Sacado::Rad::ADvar<Sacado::Fad::SFad<typename Functor::Scalar,Functor::InputType::SizeAtCompileTime> > ActiveScalar;
#endif
public:

  SacadoHessian() : Functor() {}
  SacadoHessian(const Functor& f) : Functor(f) {}

  // forward constructors
  template<typename T0>
  SacadoHessian(const T0& a0) : Functor(a0) {}
  template<typename T0, typename T1>
  SacadoHessian(const T0& a0, const T1& a1) : Functor(a0, a1) {}
  template<typename T0, typename T1, typename T2>
  SacadoHessian(const T0& a0, const T1& a1, const T2& a2) : Functor(a0, a1, a2) {}

  typedef typename Functor::Scalar Scalar;
  typedef typename Functor::InputType InputType;
  typedef typename Functor::JacobianType ValueType;
  typedef typename Functor::HessianType JacobianType;
#ifdef USE_SACADO
  typedef Eigen::Matrix<ActiveScalar, InputType::SizeAtCompileTime, 1> ActiveInput;
  typedef Eigen::Matrix<ActiveScalar, 1, Functor::ValueType::SizeAtCompileTime> ActiveValue;
#endif

  template<typename T>
  int operator() (const Eigen::Matrix<T,InputType::SizeAtCompileTime,1>& x, Eigen::Matrix<T,InputType::SizeAtCompileTime,InputType::SizeAtCompileTime>& hes) const
  {
#ifdef USE_SACADO
    ActiveInput ax = x.template cast<ActiveScalar>();
    ActiveValue av;

    for (int j=0; j<hes.cols(); j++)
      ax[j] = Sacado::Fad::SFad<Scalar,InputType::SizeAtCompileTime>(hes.cols(), j, x[j]);

    Functor::operator()(ax, av);

    Sacado::Rad::ADvar< Sacado::Fad::SFad<Scalar,InputType::SizeAtCompileTime> >::Gradcomp();
    for (int i=0; i<hes.rows(); i++)
    {
      for (int j=0; j<hes.cols(); j++)
        hes.coeffRef(i,j) = ax[i].adj().dx(j);
    }

    Sacado::Rad::ADvar< Sacado::Fad::SFad<Scalar,InputType::SizeAtCompileTime> >::aval_reset();
#else
    std::cerr << "Error: USE_SACADO is not defined in SacadoHessian::operator()\n";
    exit(-1);
#endif
    return 0;
  }

};

#endif
#endif
