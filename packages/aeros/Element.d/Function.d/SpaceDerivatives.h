#if !defined(_SPACEDERIVATIVES_H_) && defined(USE_EIGEN3)
#define _SPACEDERIVATIVES_H_

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Sparse>

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/SacadoReverseJacobian.h>

namespace Simo {

template<typename A, typename B>
struct assign_coherent_impl {
  static void run(const A& a, B& b) { b = a; }
};

template<typename A, typename B>
void assign_coherent(const A& a, B& b)
{
  assign_coherent_impl<A,B>::run(a, b);
}

template<typename Scalar, int Options, int MaxRows, int MaxCols>
struct assign_coherent_impl<Eigen::Matrix<Scalar, 1, 1, Options, MaxRows, MaxCols>, Scalar> {
  typedef Eigen::Matrix<Scalar, 1, 1, Options, MaxRows, MaxCols> A;
  typedef Scalar B;
  static void run(const A& a, B& b) { b = a[0]; }
};

template<typename Scalar, int Options, int MaxRows, int MaxCols>
struct assign_coherent_impl<Scalar, Eigen::Matrix<Scalar, 1, 1, Options, MaxRows, MaxCols> > {
  typedef Scalar A;
  typedef Eigen::Matrix<Scalar, 1, 1, Options, MaxRows, MaxCols> B;
  static void run(const A& a, B& b) { b[0] = a; }
};

template<typename Scalar, int A_Rows, int A_Cols, int A_Options, int A_MaxRows, int A_MaxCols,
                          int B_Rows, int B_Cols, int B_Options, int B_MaxRows, int B_MaxCols>
struct assign_coherent_impl<Eigen::Matrix<Scalar, A_Rows, A_Cols, A_Options, A_MaxRows, A_MaxCols>, 
                            Eigen::Matrix<Scalar, B_Rows, B_Cols, B_Options, B_MaxRows, B_MaxCols> > {
  typedef Eigen::Matrix<Scalar, A_Rows, A_Cols, A_Options, A_MaxRows, A_MaxCols> A;
  typedef Eigen::Matrix<Scalar, B_Rows, B_Cols, B_Options, B_MaxRows, B_MaxCols> B;
  static void run(const A& a, B& b) { b = Eigen::Map<B>(const_cast<Scalar*>(a.data())); }
};

template<typename Scalar, int A_Cols, int A_Options, int A_MaxRows, int A_MaxCols,
                          int B_Rows, int B_Cols, int B_Options, int B_MaxRows, int B_MaxCols,
                          int C_Options, int C_MaxRows, int C_MaxCols>
struct assign_coherent_impl<Eigen::Matrix<Scalar, B_Rows*B_Cols, A_Cols, A_Options, A_MaxRows, A_MaxCols>,
                            Eigen::Array<Eigen::Matrix<Scalar, B_Rows, B_Cols, B_Options, B_MaxRows, B_MaxCols>, 1, A_Cols, C_Options, C_MaxRows, C_MaxCols> > { 
  typedef Eigen::Matrix<Scalar, B_Rows*B_Cols, A_Cols, A_Options, A_MaxRows, A_MaxCols> A;
  typedef Eigen::Array<Eigen::Matrix<Scalar, B_Rows, B_Cols, B_Options, B_MaxRows, B_MaxCols>, 1, A_Cols, C_Options, C_MaxRows, C_MaxCols> B;
  static void run(const A& a, B& b) {
    for(int i=0; i<A_Cols; ++i) assign_coherent<Eigen::Matrix<Scalar, B_Rows*B_Cols,1>,
                                                Eigen::Matrix<Scalar, B_Rows, B_Cols, B_Options, B_MaxRows, B_MaxCols> >(a.col(i), b[i]);
  }
};

// wrapper "Functor" to support automatic and numerical differentiation of spatio-temporal
// matrix valued function of a matrix w.r.t spatial coordinates, q
template<typename _Scalar, template <typename S> class FunctionTemplate>
class SpatialView
{
  public:
    typedef _Scalar Scalar;
    enum {
      InputsAtCompileTime = FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
      ValuesAtCompileTime = FunctionTemplate<Scalar>::NumberOfValues
    };

  protected:
    const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                       FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& sconst;
    const Eigen::Array<int,
                       FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& iconst;
    Scalar t;

  public:
    SpatialView(const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const 
                Eigen::Array<int, FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>&
                _iconst, Scalar _t = 0) 
     : sconst(_sconst), iconst(_iconst), t(_t) {}

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& _q,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>& y) const
    {
      Eigen::Matrix<T,FunctionTemplate<T>::InputNumberOfRows,FunctionTemplate<T>::InputNumberOfColumns> q;
      assign_coherent(_q, q);

      // evaluate y = f(q)
      FunctionTemplate<T> f(sconst, iconst);
      assign_coherent(f(q,static_cast<T>(t)), y);

      return 1;
    }

    template<typename T>
    int operator() (const Eigen::Matrix<T,InputsAtCompileTime,1>& q,
                    Eigen::Matrix<T,ValuesAtCompileTime,1>* y) const 
    { 
      return (*this)(q,*y);
    }

    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;
    typedef Eigen::Array<Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime>,InputsAtCompileTime,1> HessianType;

    int inputs() const { return InputsAtCompileTime; }
    int values() const { return ValuesAtCompileTime; }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename _Scalar, template <typename S> class FunctionTemplate, int Options=0>
class FirstPartialSpaceDerivatives
{
  public:
    typedef _Scalar Scalar;
    typedef typename FunctionTemplate<Scalar>::ScalarConstantType ScalarConstantType;
    enum { InputNumberOfRows               = FunctionTemplate<Scalar>::InputNumberOfRows,
           InputNumberOfColumns            = FunctionTemplate<Scalar>::InputNumberOfColumns,
           NumberOfScalarConstants         = FunctionTemplate<Scalar>::NumberOfScalarConstants,
           NumberOfIntegerConstants        = FunctionTemplate<Scalar>::NumberOfIntegerConstants,
           NumberOfGeneralizedCoordinates  = InputNumberOfRows*InputNumberOfColumns,
           NumberOfValues                  = FunctionTemplate<Scalar>::NumberOfValues*NumberOfGeneralizedCoordinates
    };
    typedef Eigen::Array<typename FunctionTemplate<Scalar>::ReturnType,InputNumberOfColumns,InputNumberOfRows> ReturnType;

  protected:
    const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                       FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& sconst;
    const Eigen::Array<int,
                       FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& iconst;

  public:
    FirstPartialSpaceDerivatives(const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                                 FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const
                                 Eigen::Array<int, FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& _iconst)
     : sconst(_sconst), iconst(_iconst) {}

    typename FunctionTemplate<Scalar>::JacobianType
    operator() (const Eigen::Matrix<Scalar,InputNumberOfRows,InputNumberOfColumns>& q, Scalar t) {

      typename FunctionTemplate<Scalar>::JacobianType ret;
#ifndef AEROS_NO_AD
      SpatialView<Scalar, FunctionTemplate> f(sconst, iconst, t);
      Eigen::AutoDiffJacobian<SpatialView<Scalar, FunctionTemplate> > dfdq(f);
      Eigen::Matrix<Scalar,FunctionTemplate<Scalar>::NumberOfValues,
                    FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates> J;
      Eigen::Matrix<Scalar,FunctionTemplate<Scalar>::NumberOfValues,1> y;

      Eigen::Matrix<Scalar,FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,1> _q;
      assign_coherent(q, _q);
      dfdq(_q, &y, &J);
      assign_coherent(J, ret);
#else
      std::cerr << "Error: AEROS_NO_AD is defined in FirstPartialSpaceDerivatives::operator()\n";
      exit(-1);
#endif
      return ret;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename _Scalar, template <typename S> class FunctionTemplate>
class FirstPartialSpaceDerivatives<_Scalar, FunctionTemplate, 1>
{
  public:
    typedef _Scalar Scalar;
    typedef typename FunctionTemplate<Scalar>::ScalarConstantType ScalarConstantType;
    enum { InputNumberOfRows           = FunctionTemplate<Scalar>::InputNumberOfRows,
           InputNumberOfColumns        = FunctionTemplate<Scalar>::InputNumberOfColumns,
           NumberOfScalarConstants     = FunctionTemplate<Scalar>::NumberOfScalarConstants,
           NumberOfIntegerConstants    = FunctionTemplate<Scalar>::NumberOfIntegerConstants,
           NumberOfGeneralizedCoordinates  = InputNumberOfRows*InputNumberOfColumns,
           NumberOfValues              = FunctionTemplate<Scalar>::NumberOfValues*NumberOfGeneralizedCoordinates
    };
    typedef Eigen::Array<typename FunctionTemplate<Scalar>::ReturnType,InputNumberOfColumns,InputNumberOfRows> ReturnType;

  protected:
    const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                       FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& sconst;
    const Eigen::Array<int,
                       FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& iconst;

  public:
    FirstPartialSpaceDerivatives(const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                                 FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const
                                 Eigen::Array<int, FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& _iconst)
     : sconst(_sconst), iconst(_iconst) {}

    typename FunctionTemplate<Scalar>::JacobianType
    operator() (const Eigen::Matrix<Scalar,InputNumberOfRows,InputNumberOfColumns>& q, Scalar t) {

      typename FunctionTemplate<Scalar>::JacobianType ret;
#ifndef AEROS_NO_AD
      SpatialView<Scalar, FunctionTemplate> f(sconst, iconst, t);
      SacadoReverseJacobian<SpatialView<Scalar, FunctionTemplate> > dfdq(f);
      Eigen::Matrix<Scalar,FunctionTemplate<Scalar>::NumberOfValues,
                    FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates> J;

      Eigen::Matrix<Scalar,FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,1> _q;
      assign_coherent(q, _q);
      dfdq(_q, J);
      assign_coherent(J, ret);
#else
      std::cerr << "Error: AEROS_NO_AD is defined in FirstPartialSpaceDerivatives::operator()\n";
      exit(-1);
#endif
      return ret;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// Jacobian: the M*N matrix of the first partial space derivatives of a M*1 vector valued function of a N*1 (column) vector
template<typename _Scalar, template <typename S> class VectorValuedFunctionTemplate, int Options=0>
class Jacobian : MatrixValuedFunction<VectorValuedFunctionTemplate<_Scalar>::NumberOfGeneralizedCoordinates,
                                      VectorValuedFunctionTemplate<_Scalar>::NumberOfValues,
                                      VectorValuedFunctionTemplate<_Scalar>::NumberOfGeneralizedCoordinates,
                                      typename VectorValuedFunctionTemplate<_Scalar>::Scalar,
                                      VectorValuedFunctionTemplate<_Scalar>::NumberOfScalarConstants,
                                      VectorValuedFunctionTemplate<_Scalar>::NumberOfIntegerConstants,
                                      typename VectorValuedFunctionTemplate<_Scalar>::ScalarConstantType>
{
    typedef MatrixValuedFunction<VectorValuedFunctionTemplate<_Scalar>::NumberOfGeneralizedCoordinates,
                                 VectorValuedFunctionTemplate<_Scalar>::NumberOfValues,
                                 VectorValuedFunctionTemplate<_Scalar>::NumberOfGeneralizedCoordinates,
                                 typename VectorValuedFunctionTemplate<_Scalar>::Scalar,
                                 VectorValuedFunctionTemplate<_Scalar>::NumberOfScalarConstants,
                                 VectorValuedFunctionTemplate<_Scalar>::NumberOfIntegerConstants,
                                 typename VectorValuedFunctionTemplate<_Scalar>::ScalarConstantType> MatrixValuedFunctionBase;

  public:
    typedef typename MatrixValuedFunctionBase::Scalar Scalar;
    typedef typename MatrixValuedFunctionBase::ScalarConstantType ScalarConstantType;
    enum { InputNumberOfRows        = MatrixValuedFunctionBase::InputNumberOfRows,
           InputNumberOfColumns     = MatrixValuedFunctionBase::InputNumberOfColumns,
           NumberOfGeneralizedCoordinates = MatrixValuedFunctionBase::NumberOfGeneralizedCoordinates,
           NumberOfValuesPerColumn  = MatrixValuedFunctionBase::NumberOfValuesPerColumn,
           NumberOfValuesPerRow     = MatrixValuedFunctionBase::NumberOfValuesPerRow,
           NumberOfScalarConstants  = MatrixValuedFunctionBase::NumberOfScalarConstants,
           NumberOfIntegerConstants = MatrixValuedFunctionBase::NumberOfIntegerConstants,
           NumberOfValues           = NumberOfValuesPerRow*NumberOfValuesPerColumn
    };
    typedef typename Eigen::Matrix<Scalar,NumberOfValuesPerColumn,NumberOfValuesPerRow> ReturnType;

  protected:
    const Eigen::Array<typename VectorValuedFunctionTemplate<Scalar>::ScalarConstantType,
                       VectorValuedFunctionTemplate<Scalar>::NumberOfScalarConstants,1>& sconst;
    const Eigen::Array<int,
                       VectorValuedFunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& iconst;

  public:
    Jacobian(const Eigen::Array<typename VectorValuedFunctionTemplate<Scalar>::ScalarConstantType,
                  VectorValuedFunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const
                  Eigen::Array<int, VectorValuedFunctionTemplate<Scalar>::NumberOfIntegerConstants,1> & _iconst)
     : sconst(_sconst), iconst(_iconst) {}

    Eigen::Matrix<Scalar,NumberOfValuesPerColumn,NumberOfValuesPerRow>
    operator() (const Eigen::Matrix<Scalar,NumberOfGeneralizedCoordinates,1>& q, Scalar t)
    {
      FirstPartialSpaceDerivatives<Scalar, VectorValuedFunctionTemplate, Options> J(sconst, iconst);
      return J(q,t);
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
 
template<typename Scalar, template <typename S> class FunctionTemplate, int Options=0>
class JacobianVectorProduct : public VectorValuedFunction<FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
                                                          FunctionTemplate<Scalar>::NumberOfValues,
                                                          Scalar,
                                                          FunctionTemplate<Scalar>::NumberOfScalarConstants+
                                                          FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
                                                          FunctionTemplate<Scalar>::NumberOfIntegerConstants,
                                                          typename FunctionTemplate<Scalar>::ScalarConstantType>
{
  typedef VectorValuedFunction<FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
                               FunctionTemplate<Scalar>::NumberOfValues,
                               Scalar,
                               FunctionTemplate<Scalar>::NumberOfScalarConstants+
                               FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
                               FunctionTemplate<Scalar>::NumberOfIntegerConstants,
                               typename FunctionTemplate<Scalar>::ScalarConstantType> Base;

    const Eigen::Array<int, Base::NumberOfIntegerConstants, 1>& iconst;
    Eigen::Array<typename Base::ScalarConstantType, FunctionTemplate<Scalar>::NumberOfScalarConstants, 1> sconst;
    Eigen::Array<typename Base::ScalarConstantType, FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates, 1> V;

  public:
    JacobianVectorProduct(const Eigen::Array<typename Base::ScalarConstantType, Base::NumberOfScalarConstants, 1>& _sconst,
                          const Eigen::Array<int, Base::NumberOfIntegerConstants, 1> & _iconst)
     : iconst(_iconst)
    {
      sconst = _sconst.template head<FunctionTemplate<Scalar>::NumberOfScalarConstants>();
      V = _sconst.template tail<FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates>();
    }

    typename Base::ReturnType
    operator() (const Eigen::Matrix<Scalar,Base::NumberOfGeneralizedCoordinates,1>& q, Scalar t)
    {
      typename Base::ReturnType ret;
      Jacobian<Scalar, FunctionTemplate, Options> J(sconst, iconst);
      Eigen::Matrix<Scalar,1,Base::NumberOfValues> JV = J(q,t)*V.matrix().template cast<Scalar>();
      assign_coherent(JV, ret);
      return ret;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// Hessian: the N*N matrix of the second partial space derivatives of a scalar valued function of a N*1 (column) vector
template<typename Scalar, template <typename S> class FunctionTemplate>
class Hessian : public MatrixValuedFunction<FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
                                            FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
                                            FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
                                            Scalar,
                                            FunctionTemplate<Scalar>::NumberOfScalarConstants,
                                            FunctionTemplate<Scalar>::NumberOfIntegerConstants,
                                            typename FunctionTemplate<Scalar>::ScalarConstantType>
{
    const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                       FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& sconst;
    const Eigen::Array<int,
                       FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& iconst;

  public:
    Hessian(const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
            FunctionTemplate<Scalar>::NumberOfScalarConstants,1>& _sconst, const
            Eigen::Array<int, FunctionTemplate<Scalar>::NumberOfIntegerConstants,1> & _iconst)
     : sconst(_sconst), iconst(_iconst) {}

    Eigen::Matrix<Scalar,FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
                  FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates>
    operator() (const Eigen::Matrix<Scalar,FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,1>& q, Scalar t)
    {
      Eigen::Matrix<Scalar,FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
                    FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates> H;
#ifndef AEROS_NO_AD
      SpatialView<Scalar, FunctionTemplate> V(sconst, iconst, static_cast<Scalar>(t));
      SacadoHessian<SpatialView<Scalar, FunctionTemplate> > d2Vdq2(V);
      d2Vdq2(q, H);
#else
      std::cerr << "Error: AEROS_NO_AD is defined in Hessian::operator()\n";
      exit(-1);
#endif
      return H;
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace Simo
#endif
