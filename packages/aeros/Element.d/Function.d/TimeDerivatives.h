#if !defined(_TIMEDERIVATIVES_H_) && defined(USE_EIGEN3)
#define _TIMEDERIVATIVES_H_

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Sparse>

#include <Element.d/Function.d/SacadoReverseJacobian.h>

namespace Simo {

// wrapper "Functor" to support automatic and numerical differentiation of spatio-temporal
// matrix valued function of a matrix w.r.t the temporal coordinate, t
template<typename _Scalar, template <typename S> class FunctionTemplate>
class TemporalView
{
  public:
    typedef _Scalar Scalar;
    enum {
      InputsAtCompileTime = 1,
      ValuesAtCompileTime = FunctionTemplate<Scalar>::NumberOfValues
    };

  private:
    const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                       FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& sconst;
    const Eigen::Array<int,
                       FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& iconst;
    const Eigen::Matrix<Scalar, FunctionTemplate<Scalar>::InputNumberOfRows, FunctionTemplate<Scalar>::InputNumberOfColumns>& q;

  public:
    TemporalView(const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType, FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst,
                 const Eigen::Array<int, FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& _iconst,
                 const Eigen::Matrix<Scalar, FunctionTemplate<Scalar>::InputNumberOfRows, FunctionTemplate<Scalar>::InputNumberOfColumns>& _q) 
     : sconst(_sconst), iconst(_iconst), q(_q) {}

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& t,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>& y) const 
    {
      // evaluate y = f(t)
      FunctionTemplate<T> f(sconst, iconst);
      assign_coherent(f(q.template cast<T>(), t[0]), y);

      return 1;
    }

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& t,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>* f) const 
    { 
      return (*this)(t,*f);
    }

    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,InputsAtCompileTime> HessianType;

    int inputs() const { return InputsAtCompileTime; }
    int values() const { return ValuesAtCompileTime; }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename _Scalar, template <typename S> class FunctionTemplate, int Options=0>
class FirstPartialTimeDerivative
{
  public:
    typedef _Scalar Scalar;
    typedef typename FunctionTemplate<Scalar>::ScalarConstantType ScalarConstantType;
    enum { InputNumberOfRows           = FunctionTemplate<Scalar>::InputNumberOfRows,
           InputNumberOfColumns        = FunctionTemplate<Scalar>::InputNumberOfColumns,
           NumberOfScalarConstants     = FunctionTemplate<Scalar>::NumberOfScalarConstants,
           NumberOfIntegerConstants    = FunctionTemplate<Scalar>::NumberOfIntegerConstants,
           NumberOfGeneralizedCoordinates  = InputNumberOfRows*InputNumberOfColumns,
           NumberOfValues              = FunctionTemplate<Scalar>::NumberOfValues
    };
    typedef typename FunctionTemplate<Scalar>::ReturnType ReturnType;

  protected:
    const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                       FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& sconst;
    const Eigen::Array<int,
                       FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& iconst;

  public:
    FirstPartialTimeDerivative(const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                               FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const
                               Eigen::Array<int, FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& _iconst)
     : sconst(_sconst), iconst(_iconst) {}

    typename FunctionTemplate<Scalar>::ReturnType
    operator() (const Eigen::Matrix<Scalar,InputNumberOfRows,InputNumberOfColumns>& q, Scalar t) {

      typename FunctionTemplate<Scalar>::ReturnType ret;
#ifndef AEROS_NO_AD
      TemporalView<Scalar, FunctionTemplate> f(sconst, iconst, q);
      Eigen::AutoDiffJacobian<TemporalView<Scalar,FunctionTemplate> > dfdt(f);
      Eigen::Matrix<Scalar,FunctionTemplate<Scalar>::NumberOfValues,1> fdot;
      Eigen::Matrix<Scalar,FunctionTemplate<Scalar>::NumberOfValues,1> y;

      Eigen::Matrix<Scalar,1,1> _t;
      assign_coherent(t, _t);
      dfdt(_t, &y, &fdot);
      assign_coherent(fdot, ret);
#else
      std::cerr << "Error: AEROS_NO_AD is defined in FirstPartialTimeDerivative::operator()\n";
      exit(-1);
#endif
      return ret;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename _Scalar, template <typename S> class FunctionTemplate>
class FirstPartialTimeDerivative<_Scalar, FunctionTemplate, 1>
{
  public:
    typedef _Scalar Scalar;
    typedef typename FunctionTemplate<Scalar>::ScalarConstantType ScalarConstantType;
    enum { InputNumberOfRows           = FunctionTemplate<Scalar>::InputNumberOfRows,
           InputNumberOfColumns        = FunctionTemplate<Scalar>::InputNumberOfColumns,
           NumberOfScalarConstants     = FunctionTemplate<Scalar>::NumberOfScalarConstants,
           NumberOfIntegerConstants    = FunctionTemplate<Scalar>::NumberOfIntegerConstants,
           NumberOfGeneralizedCoordinates  = InputNumberOfRows*InputNumberOfColumns,
           NumberOfValues              = FunctionTemplate<Scalar>::NumberOfValues
    };
    typedef typename FunctionTemplate<Scalar>::ReturnType ReturnType;

  protected:
    const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                       FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& sconst;
    const Eigen::Array<int,
                       FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& iconst;

  public:
    FirstPartialTimeDerivative(const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                               FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const
                               Eigen::Array<int, FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& _iconst)
     : sconst(_sconst), iconst(_iconst) {}

    typename FunctionTemplate<Scalar>::ReturnType
    operator() (const Eigen::Matrix<Scalar,InputNumberOfRows,InputNumberOfColumns>& q, Scalar t) {

      typename FunctionTemplate<Scalar>::ReturnType ret;
#ifndef AEROS_NO_AD
      TemporalView<Scalar, FunctionTemplate> f(sconst, iconst, q);
      SacadoReverseJacobian<TemporalView<Scalar,FunctionTemplate> > dfdt(f);
      Eigen::Matrix<Scalar,FunctionTemplate<Scalar>::NumberOfValues,1> fdot;

      Eigen::Matrix<Scalar,1,1> _t;
      assign_coherent(t, _t);
      dfdt(_t, fdot);
      assign_coherent(fdot, ret);
#else
      std::cerr << "Error: AEROS_NO_AD is defined in FirstPartialTimeDerivative::operator()\n";
      exit(-1);
#endif
      return ret;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#if ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && defined(HAS_CXX11_TEMPLATE_ALIAS)
template<typename _Scalar, template <typename S> class FunctionTemplate, int Options=0>
class SecondPartialTimeDerivative
{
  public:
    typedef _Scalar Scalar;
    typedef typename FunctionTemplate<Scalar>::ScalarConstantType ScalarConstantType;
    enum { InputNumberOfRows           = FunctionTemplate<Scalar>::InputNumberOfRows,
           InputNumberOfColumns        = FunctionTemplate<Scalar>::InputNumberOfColumns,
           NumberOfScalarConstants     = FunctionTemplate<Scalar>::NumberOfScalarConstants,
           NumberOfIntegerConstants    = FunctionTemplate<Scalar>::NumberOfIntegerConstants,
           NumberOfGeneralizedCoordinates  = InputNumberOfRows*InputNumberOfColumns,
           NumberOfValues              = FunctionTemplate<Scalar>::NumberOfValues
    };
    typedef typename FunctionTemplate<Scalar>::ReturnType ReturnType;

  protected:
    const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                       FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& sconst;
    const Eigen::Array<int,
                       FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& iconst;

    template <typename S>
    using FirstPartialTimeDerivativeFunctionTemplate = FirstPartialTimeDerivative<S, FunctionTemplate>;

  public:
    SecondPartialTimeDerivative(const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                                FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const
                                Eigen::Array<int, FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& _iconst)
     : sconst(_sconst), iconst(_iconst) {}

    typename FunctionTemplate<Scalar>::ReturnType
    operator() (const Eigen::Matrix<Scalar,InputNumberOfRows,InputNumberOfColumns>& q, Scalar t) {

      FirstPartialTimeDerivative<Scalar, FirstPartialTimeDerivativeFunctionTemplate, 1> d2fdt2(sconst, iconst);
      return d2fdt2(q,t);
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
#endif

} // namespace Simo
#endif
