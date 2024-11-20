#if !defined(_CONTRACTIONVIEW_H_) && defined(USE_EIGEN3)
#define _CONTRACTIONVIEW_H_

namespace Simo {

// Scalar-valued function view of the contraction of a vector-valued function with a vector
// The gradient of this function is equal to the Jacobian of the vector-valued function contracted with the vector, etc.
template<typename Scalar, template <typename S> class FunctionTemplate>
class VectorVectorContractionView : ScalarValuedFunction<FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
                                                         Scalar,
                                                         FunctionTemplate<Scalar>::NumberOfScalarConstants+FunctionTemplate<Scalar>::NumberOfValues,
                                                         FunctionTemplate<Scalar>::NumberOfIntegerConstants,
                                                         typename FunctionTemplate<Scalar>::ScalarConstantType>
{
  public:
    typedef typename FunctionTemplate<Scalar>::ScalarConstantType ScalarConstantType;
    enum { InputNumberOfRows              = FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
           InputNumberOfColumns           = 1,
           NumberOfGeneralizedCoordinates = FunctionTemplate<Scalar>::NumberOfGeneralizedCoordinates,
           NumberOfValues                 = 1,
           NumberOfScalarConstants        = FunctionTemplate<Scalar>::NumberOfScalarConstants+FunctionTemplate<Scalar>::NumberOfValues,
           NumberOfIntegerConstants       = FunctionTemplate<Scalar>::NumberOfIntegerConstants
    };
    using ScalarValuedFunction<InputNumberOfRows,Scalar,NumberOfScalarConstants,NumberOfIntegerConstants,ScalarConstantType>::ReturnType;
    using ScalarValuedFunction<InputNumberOfRows,Scalar,NumberOfScalarConstants,NumberOfIntegerConstants,ScalarConstantType>::JacobianType;

  protected:
    Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                 FunctionTemplate<Scalar>::NumberOfScalarConstants,1> sconst;
    const Eigen::Array<int,
                       FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& iconst;
    Eigen::Matrix<typename FunctionTemplate<Scalar>::ScalarConstantType,
                  FunctionTemplate<Scalar>::NumberOfValues,1> V;

  public:
    VectorVectorContractionView(const Eigen::Array<typename FunctionTemplate<Scalar>::ScalarConstantType,
                                FunctionTemplate<Scalar>::NumberOfScalarConstants+FunctionTemplate<Scalar>::NumberOfValues,1>& _sconst,
                                const Eigen::Array<int, FunctionTemplate<Scalar>::NumberOfIntegerConstants,1>& _iconst) 
     : iconst(_iconst)
    {
      for(int i=0; i<FunctionTemplate<Scalar>::NumberOfScalarConstants; ++i) sconst[i] = _sconst[i];
      for(int i=0; i<FunctionTemplate<Scalar>::NumberOfValues; ++i) V[i] = _sconst[FunctionTemplate<Scalar>::NumberOfScalarConstants+i];
    }

    Scalar operator() (const Eigen::Matrix<Scalar,NumberOfGeneralizedCoordinates,1>& q, Scalar t)
    {
      // evaluate y = f(q)*V
      FunctionTemplate<Scalar> f(sconst, iconst);
      return f(q, t).dot(V.template cast<Scalar>());
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}

#endif
