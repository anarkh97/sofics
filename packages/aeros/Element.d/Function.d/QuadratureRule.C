#include <iostream>

template<typename T>
void
GaussLegendre<T>::getAbscissaAndWeight(int n, int i, T& xi, T& weight) 
{
//
//     GaussLegendre::getAbscissaAndWeight resturns the absissae and weight factors of
//     the n-th Gauss-Legendre integration rule over the
//     segment (-1,+1).
//
//     Input arguments:
//
//       n         Number of points in the integration rule (1 to 5)
//       i         Index of sample point (0 to n-1)
//
//     Outputs arguments:
//
//       xi        Absissa of sample point (zero of Legendre polynomial)
//       weight    Weight factor
//
  switch(n) {
    case 1 :
      switch(i) {
        case 0 : 
          xi =      0;
          weight =  2;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in GaussLegendre::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
    case 2 :
      switch(i) {
        case 0 :
          xi =     -1/std::sqrt(3);
          weight =  1;
          break;
        case 1 : 
          xi =      1/std::sqrt(3);
          weight =  1;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in GaussLegendre::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
    case 3 :
      switch(i) {
        case 0 :
          xi =     -std::sqrt(15)/5;
          weight =  5/9.;
          break;
        case 1 :
          xi =      0;
          weight =  8/9.;
          break;
        case 2 :
          xi =      std::sqrt(15)/5;
          weight =  5/9.;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in GaussLegendre::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
    case 4 :
      switch(i) {
        case 0 :
          xi =      sqrt((3-2*sqrt(6/5.))/7);
          weight =  (18+sqrt(30))/36;
          break;
        case 1 : 
          xi =     -sqrt((3-2*sqrt(6/5.))/7);
          weight =  (18+sqrt(30))/36;
          break;
        case 2 :
          xi =      sqrt((3+2*sqrt(6/5.))/7);
          weight =  (18-sqrt(30))/36;
          break;
        case 3 :
          xi =     -sqrt((3+2*sqrt(6/5.))/7);
          weight =  (18-sqrt(30))/36;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in GaussLegendre::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
    case 5 :
      switch(i) {
        case 0 :
          xi =      0;
          weight =  128/225.;
          break;
        case 1 :
          xi =      sqrt(5-2*sqrt(10/7.))/3;
          weight =  (322+13*sqrt(70))/900;
          break;
        case 2 :
          xi =      -sqrt(5-2*sqrt(10/7.))/3;
          weight =  (322+13*sqrt(70))/900;
          break;
        case 3 :
          xi =      sqrt(5+2*sqrt(10/7.))/3;
          weight =  (322-13*sqrt(70))/900;
          break;
        case 4 :
          xi =      -sqrt(5+2*sqrt(10/7.))/3;
          weight =  (322-13*sqrt(70))/900;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in GaussLegendre::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
   case 10 :
      switch(i) {
        case 0 :
          xi =      -0.973906528517171720;
          weight =  0.066671344308688138;
          break;
        case 1 :
          xi =      -0.865063366688984511;
          weight =  0.149451349150580593;
          break;
        case 2 :
          xi =      -0.679409568299024406;
          weight =  0.219086362515982044;
          break;
        case 3 :
          xi =      -0.433395394129247191;
          weight =  0.269266719309996355;
          break;
        case 4 :
          xi =      -0.148874338981631211;
          weight =  0.295524224714752870;
          break;
        case 5 :
          xi =      0.148874338981631211;
          weight =  0.295524224714752870;
          break;
        case 6 :
          xi =      0.433395394129247191;
          weight =  0.269266719309996355;
          break;
        case 7 :
          xi =      0.679409568299024406;
          weight =  0.219086362515982044;
          break;
        case 8 :
          xi =      0.865063366688984511;
          weight =  0.149451349150580593;
          break;
        case 9 :
          xi =      0.973906528517171720;
          weight =  0.066671344308688138;
          break;

        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in GaussLegendre::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
    default :
      std::cerr << " *** ERROR: number of points = " << n << " not supported for GaussLegendre quadrature rule\n";
      exit(-1);
  }
}

template<typename T>
void
CompositeTrapezoidalRule<T>::getAbscissaAndWeight(int n, int i, T& xi, T& weight)
{
  // TODO this can be generalized see: http://en.wikipedia.org/wiki/Trapezoidal_rule
  switch(n) {
    case 2 : 
      switch(i) {
        case 0 :
          xi =      -1;
          weight =  1;
          break;
        case 1:
          xi =      1;
          weight =  1;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in CompositeTrapezoidalRule::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
    case 3 :
      switch(i) {
        case 0 :
          xi =      -1;
          weight =  0.5;
          break;
        case 1:
          xi =      0;
          weight =  1;
          break;
        case 2:
          xi =      1;
          weight =  0.5;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in CompositeTrapezoidalRule::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
    case 4 :
      switch(i) {
        case 0 :
          xi =      -1;
          weight =  1/3.;
          break;
        case 1:
          xi =      -1/3.;
          weight =  2/3.;
          break;
        case 2:
          xi =      1/3.;
          weight =  2/3.;
          break;
        case 3 :
          xi =      1;
          weight =  1/3.;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in CompositeTrapezoidalRule::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;

    default :
      std::cerr << " *** ERROR: number of points = " << n << " not supported for composite trapezoidal quadrature rule\n";
      exit(-1);
  }
}

template<typename T>
void
CompositeMidpointRule<T>::getAbscissaAndWeight(int n, int i, T& xi, T& weight)
{
  switch(n) {
    case 1 :
      switch(i) {
        case 0 :
          xi =      0;
          weight =  2;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in CompositeMidpointRule::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
    case 2 :
      switch(i) {
        case 0 :
          xi =      -0.5;
          weight =  1;
          break;
        case 1:
          xi =      0.5;
          weight =  1;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in CompositeMidpointRule::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
    case 3 :
      switch(i) {
        case 0 :
          xi =      -2/3.;
          weight =  2/3.;
          break;
        case 1:
          xi =      0;
          weight =  2/3.;
          break;
        case 2:
          xi =      2/3.;
          weight =  2/3.;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in CompositeMidpointRule::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;

    default :
      std::cerr << " *** ERROR: number of points = " << n << " not supported for composite midpoint quadrature rule\n";
      exit(-1);
  }
}

template<typename T>
void
ClosedNewtonCotes<T>::getAbscissaAndWeight(int n, int i, T& xi, T& weight)
{
  switch(n) {
    case 2 : // trapezoidal rule
      switch(i) {
        case 0 :
          xi =      -1;
          weight =  1;
          break;
        case 1:
          xi =      1;
          weight =  1;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in ClosedNewtonCotes::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
    case 3 : // simpson's rule
      switch(i) {
        case 0 :
          xi =      -1;
          weight =  1/3.;
          break;
        case 1:
          xi =      0;
          weight =  4/3.;
          break;
        case 2:
          xi =      1;
          weight =  1/3.;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in ClosedNewtonCotes::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
    case 4 : // simpson's 3/8 rule
      switch(i) {
        case 0 : 
          xi =      -1; 
          weight =  1/4.;
          break;
        case 1:
          xi =      -1/3.;
          weight =  3/4.;
          break;
        case 2:
          xi =      1/3.;
          weight =  3/4.;
          break;
        case 3 : 
          xi =      1;  
          weight =  1/4.;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in ClosedNewtonCotes::getAbscissaAndWeight\n";
          exit(-1);
      }   
      break;
    case 5 : // Booles's rule
      switch(i) {
        case 0 :
          xi =      -1;
          weight =  7/45.;
          break;
        case 1:
          xi =      -0.5;
          weight =  32/45.;
          break;
        case 2:
          xi =      0;
          weight =  12/45.;
          break;
        case 3 :
          xi =      0.5;
          weight =  32/45.;
          break;
        case 4 :
          xi =      1;
          weight =  7/45.;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in ClosedNewtonCotes::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;

    default :
      std::cerr << " *** ERROR: number of points = " << n << " not supported for closed Newton Cotes quadrature rule\n";
      exit(-1);
  }
}

template<typename T>
void
OpenNewtonCotes<T>::getAbscissaAndWeight(int n, int i, T& xi, T& weight)
{
  switch(n) {
    case 1 : // midpoint rule
      switch(i) {
        case 0 : 
          xi =      0;  
          weight =  2;  
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in OpenNewtonCotes::getAbscissaAndWeight\n";
          exit(-1);
      }   
      break;
    case 2 : // unnamed rule
      switch(i) {
        case 0 :
          xi =      -1/3.;
          weight =  1;
          break;
        case 1:
          xi =      1/3.;
          weight =  1;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in OpenNewtonCotes::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
    case 3 : // Milne's rule
      switch(i) {
        case 0 :
          xi =      -0.5;
          weight =  4/3.;
          break;
        case 1:
          xi =      0;
          weight =  -2/3.;
          break;
        case 2:
          xi =      0.5;
          weight =  4/3.;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in OpenNewtonCotes::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
    case 4 : // unnamed rule
      switch(i) {
        case 0 : 
          xi =      -0.6; 
          weight =  11/12.;
          break;
        case 1:
          xi =      -0.2;
          weight =  1/12.;
          break;
        case 2:
          xi =      0.2;
          weight =  1/12.;
          break;
        case 3 : 
          xi =      0.6;  
          weight =  11/12.;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in OpenNewtonCotes::getAbscissaAndWeight\n";
          exit(-1);
      }   
      break;

    default :
      std::cerr << " *** ERROR: number of points = " << n << " not supported for open Newton Cotes quadrature rule\n";
      exit(-1);
  }
}

template<typename T>
void
CurtisClenshaw<T>::getAbscissaAndWeight(int n, int i, T& xi, T& weight)
{
  switch(n) {
    case 1 :
      switch(i) {
        case 0 : 
          xi =      0.0000000000000001;  
          weight =  2;  
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in CurtisClenshaw::getAbscissaAndWeight\n";
          exit(-1);
      }   
      break;
    case 2 :
      switch(i) {
        case 0 :
          xi =      -1;
          weight =  0.0000000000000002;
          break;
        case 1:
          xi =      0;
          weight =  2;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in CurtisClenshaw::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
    case 3 :
      switch(i) {
        case 0 :
          xi =      -1;
          weight =  0.3333333333333335;
          break;
        case 1:
          xi =      0.0000000000000001;
          weight =  1.3333333333333330;
          break;
        case 2:
          xi =      1;
          weight =  0.3333333333333333;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in CurtisClenshaw::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;
    case 4 :
      switch(i) {
        case 0 : 
          xi =      -1; 
          weight =  0.3333333333333335;
          break;
        case 1:
          xi =      -0.7071067811865475;
          weight =  0.0000000000000002;
          break;
        case 2:
          xi =      0.0000000000000001;
          weight =  1.3333333333333330;
          break;
        case 3 : 
          xi =      1;  
          weight =  0.3333333333333334;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in CurtisClenshaw::getAbscissaAndWeight\n";
          exit(-1);
      }   
      break;
    case 5 :
      switch(i) {
        case 0 :
          xi =      -1;
          weight =  0.0666666666666669;
          break;
        case 1:
          xi =     -0.7071067811865475; 
          weight =  0.5333333333333333;
          break;
        case 2:
          xi =     0.0000000000000001; 
          weight =  0.7999999999999997;
          break;
        case 3 :
          xi =      0.7071067811865475;
          weight =  0.5333333333333332;
          break;
        case 4 :
          xi =      1;
          weight =  0.0666666666666666;
          break;
        default :
          std::cerr << " *** ERROR: index of sample point " << i << " outside interval [0,n) in CurtisClenshaw::getAbscissaAndWeight\n";
          exit(-1);
      }
      break;

    default :
      std::cerr << " *** ERROR: number of points = " << n << " not supported for Curtis-Clenshaw quadrature rule\n";
      exit(-1);
  }
}

template<typename T, template<typename S> class QuadratureRule, int dim, typename VecType>
RepeatedQuadratureRule<T,QuadratureRule,dim,VecType>::RepeatedQuadratureRule(int nq)
{
  for(int j = 0; j < dim; ++j)
    n[j] = nq;
}

template<typename T, template<typename S> class QuadratureRule, int dim, typename VecType>
RepeatedQuadratureRule<T,QuadratureRule,dim,VecType>::RepeatedQuadratureRule(int _n[dim])
{
  for(int j = 0; j < dim; ++j)
    n[j] = _n[j];
}

template<typename T, template<typename S> class QuadratureRule, int dim, typename VecType>
int
RepeatedQuadratureRule<T,QuadratureRule,dim,VecType>::getN()
{
  int N = 1;
  for(int j = 0; j < dim; ++j) N *= n[j];
  return N;
}

template<typename T, template<typename S> class QuadratureRule, int dim, typename VecType>
void
RepeatedQuadratureRule<T,QuadratureRule,dim,VecType>::getAbscissaAndWeight(int i, VecType& xi, T& weight)
{
  T wj;
  int kk = getN();
  weight = 1;
  for(int j = 0; j < dim; ++j) {
    kk /= n[j];
    QuadratureRule<T>::getAbscissaAndWeight(n[j], (i/kk)%n[j], xi[j], wj);
    weight *= wj;
  }
}

template<typename T, typename VecType>
TriangleQuadratureRule<T,VecType>::TriangleQuadratureRule(int _n)
{
  n = _n;
}

template<typename T, typename VecType>
int
TriangleQuadratureRule<T,VecType>::getN()
{
  return n;
}

extern void getGaussPtOnTriangle(int, int, double&, double&, double&, double&); // Mortar.d/Divers.d/Divers.C

template<typename T, typename VecType>
void
TriangleQuadratureRule<T,VecType>::getAbscissaAndWeight(int i, VecType& xi, T& weight)
{
  T zeta;
  getGaussPtOnTriangle(n, i+1, xi[0], xi[1], zeta, weight);
}

