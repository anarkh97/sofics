#if !defined(_SPATIOTEMPORALFUNCTION_H_) && defined(USE_EIGEN3)
#define _SPATIOTEMPORALFUNCTION_H_

#include <Eigen/Core>

// Scalar-valued spatio-temporal function, takes (q,t) as input where q is a vector (spatial) and t is a scalar (temporal)
template<int _NumberOfGeneralizedCoordinates,
         typename _Scalar,
         int _NumberOfScalarConstants = 0,
         int _NumberOfIntegerConstants = 0,
         typename _ScalarConstantType = double>
class ScalarValuedFunction {

  public:
    typedef _Scalar Scalar;
    typedef _ScalarConstantType ScalarConstantType;
    enum { InputNumberOfRows              = _NumberOfGeneralizedCoordinates,
           InputNumberOfColumns           = 1,
           NumberOfGeneralizedCoordinates = _NumberOfGeneralizedCoordinates,
           NumberOfValues                 = 1,
           NumberOfScalarConstants        = _NumberOfScalarConstants,
           NumberOfIntegerConstants       = _NumberOfIntegerConstants
    };
    typedef Scalar ReturnType;
    typedef Eigen::Matrix<Scalar,1,NumberOfGeneralizedCoordinates> JacobianType;

    virtual Scalar operator() (const Eigen::Matrix<Scalar,NumberOfGeneralizedCoordinates,1>& q, Scalar t) = 0;
};

// Vector-valued spatio-temporal function, takes (q,t) as input where q is a vector (spatial) and t is a scalar (temporal)
template<int _NumberOfGeneralizedCoordinates,
         int _NumberOfValues,
         typename _Scalar,
         int _NumberOfScalarConstants = 0,
         int _NumberOfIntegerConstants = 0,
         typename _ScalarConstantType = double>
class VectorValuedFunction {

  public:
    typedef _Scalar Scalar;
    typedef _ScalarConstantType ScalarConstantType;
    enum { InputNumberOfRows              = _NumberOfGeneralizedCoordinates,
           InputNumberOfColumns           = 1,
           NumberOfGeneralizedCoordinates = _NumberOfGeneralizedCoordinates,
           NumberOfValues                 = _NumberOfValues,
           NumberOfScalarConstants        = _NumberOfScalarConstants,
           NumberOfIntegerConstants       = _NumberOfIntegerConstants
    };
    typedef Eigen::Matrix<Scalar,NumberOfValues,1> ReturnType;
    typedef Eigen::Matrix<Scalar,NumberOfValues,NumberOfGeneralizedCoordinates> JacobianType;

    virtual ReturnType operator() (const Eigen::Matrix<Scalar,NumberOfGeneralizedCoordinates,1>& q, Scalar t) = 0;
};

// Matrix-valued spatio-temporal function, takes (q,t) as input where q is a vector (spatial) and t is a scalar (temporal)
template<int _NumberOfGeneralizedCoordinates,
         int _NumberOfValuesPerColumn,
         int _NumberOfValuesPerRow,
         typename _Scalar,
         int _NumberOfScalarConstants = 0,
         int _NumberOfIntegerConstants = 0,
         typename _ScalarConstantType = double>
class MatrixValuedFunction {

  public:
    typedef _Scalar Scalar;
    typedef _ScalarConstantType ScalarConstantType;
    enum { InputNumberOfRows              = _NumberOfGeneralizedCoordinates,
           InputNumberOfColumns           = 1,
           NumberOfGeneralizedCoordinates = _NumberOfGeneralizedCoordinates,
           NumberOfValues                 = _NumberOfValuesPerColumn*_NumberOfValuesPerRow,
           NumberOfValuesPerColumn        = _NumberOfValuesPerColumn,
           NumberOfValuesPerRow           = _NumberOfValuesPerRow,
           NumberOfScalarConstants        = _NumberOfScalarConstants,
           NumberOfIntegerConstants       = _NumberOfIntegerConstants
    };
    typedef Eigen::Matrix<Scalar,NumberOfValuesPerColumn,NumberOfValuesPerRow> ReturnType;
    typedef Eigen::Array<ReturnType,1,NumberOfGeneralizedCoordinates> JacobianType;

    virtual ReturnType operator() (const Eigen::Matrix<Scalar,NumberOfGeneralizedCoordinates,1>& q, Scalar t) = 0;
};
#endif
