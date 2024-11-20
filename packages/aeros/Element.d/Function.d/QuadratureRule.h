#ifndef _QUADRATURERULE_H_
#define _QUADRATURERULE_H_

// (univariate) quadrature rule
template<typename T>
class QuadratureRule
{
  public:
    virtual void getAbscissaAndWeight(int n, int i, T& xi, T& weight) = 0;
};


template<typename T>
class GaussLegendre : public QuadratureRule<T>
{
  public:
    void getAbscissaAndWeight(int n, int i, T& xi, T& weight);
};

template<typename T>
class CompositeTrapezoidalRule : public QuadratureRule<T>
{
  public:
    void getAbscissaAndWeight(int n, int i, T& xi, T& weight);
};

template<typename T>
class CompositeMidpointRule : public QuadratureRule<T>
{
  public:
    void getAbscissaAndWeight(int n, int i, T& xi, T& weight);
};

template<typename T>
class ClosedNewtonCotes : public QuadratureRule<T>
{
  public:
    void getAbscissaAndWeight(int n, int i, T& xi, T& weight);
};

template<typename T>
class OpenNewtonCotes : public QuadratureRule<T>
{
  public:
    void getAbscissaAndWeight(int n, int i, T& xi, T& weight);
};

template<typename T>
class CurtisClenshaw : public QuadratureRule<T>
{
  public:
    void getAbscissaAndWeight(int n, int i, T& xi, T& weight);
};

// multivariate quadrature rule
template<typename T, int dim, typename VecType = T[dim]>
class CubatureRule
{
  public:
    virtual int getN() = 0;
    virtual void getAbscissaAndWeight(int i, VecType& xi, T& weight) = 0;
};

template<typename T, template<typename S> class QuadratureRule, int dim, typename VecType = T[dim]>
class RepeatedQuadratureRule : public CubatureRule<T, dim, VecType>,
                               private QuadratureRule<T>
{
    int n[dim]; // number of integration points in each dimension
  public:
    // Isotropic constructor
    RepeatedQuadratureRule(int n);
    // Anisotropic constructor
    RepeatedQuadratureRule(int n[dim]);

    int getN();
    void getAbscissaAndWeight(int i, VecType& xi, T& weight);

};

// for convenience
typedef RepeatedQuadratureRule<double, GaussLegendre, 3> GaussLegendre3d;
typedef RepeatedQuadratureRule<double, GaussLegendre, 2> GaussLegendre2d;
typedef RepeatedQuadratureRule<double, GaussLegendre, 1> GaussLegendre1d;

template<typename T, typename VecType = T[2]>
class TriangleQuadratureRule : public CubatureRule<T, 2, VecType>
{
    int n;
  public:
    TriangleQuadratureRule(int n);

    int getN();
    void getAbscissaAndWeight(int i, VecType& xi, T& weight);
};

#ifdef _TEMPLATE_FIX_
  #include <Element.d/Function.d/QuadratureRule.C>
#endif
#endif
