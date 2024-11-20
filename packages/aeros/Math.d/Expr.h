#ifndef _EXPR_H_
#define _EXPR_H_
#include <complex>
//------------------------------------------------------------------------------

template<class A, class B>
class ProdRes {
  public:
    typedef B ResType;
};


template<>
class ProdRes<complex<double>, double> {
  public:
    typedef complex<double> ResType;
};

//------------------------------------------------------------------------------

template<class T, class Scalar, class IType = typename T::InfoType>
class Expr {

public:

  typedef IType InfoType;
  InfoType inf;
  T x;

  Expr(T v) : inf(v.info()),  x(v) {}
  Expr(T v, IType i) : inf(i), x(v) { }

  Scalar operator[] (int i) const { return x[i]; }
  IType info() const { return inf; }

};

//------------------------------------------------------------------------------
template<class T1, class T2, class Scalar, class IType = typename T1::InfoType>
class Sum {
public:
  typedef IType InfoType;
private:
  T1 a;
  T2 b;
  InfoType len;

public:

  Sum(T1 aa, T2 bb, InfoType l) : a(aa), b(bb), len(l) { }

  Scalar operator[](int i) const { return a[i]+b[i]; }
  InfoType info() const { return len; }

};
//------------------------------------------------------------------------------

template<class T1, class T2, class Scalar>
inline
Expr<Sum<T1, T2, Scalar>, Scalar>
operator+(const Expr<T1, Scalar> &x1, const Expr<T2, Scalar> &x2)
{

  return Expr<Sum<T1, T2, Scalar>, Scalar>
    ( Sum<T1, T2, Scalar>(x1.x, x2.x, x1.info()) );

}

//------------------------------------------------------------------------------

template<class T1, class T2, class Scalar, class IType = typename T1::InfoType>
class Diff {
public:
  typedef IType InfoType;
private:
  T1 a;
  T2 b;
  InfoType len;

public:

  Diff(T1 aa, T2 bb, InfoType l) : a(aa), b(bb), len(l) { }

  Scalar operator[](int i) const { return a[i]-b[i]; }
  InfoType info() const { return len; }

};
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------

template<class T1, class T2, class Scalar>
inline
Expr<Diff<T1, T2, Scalar>, Scalar>
operator-(const Expr<T1, Scalar> &x1, const Expr<T2, Scalar> &x2)
{

  return Expr<Diff<T1, T2, Scalar>, Scalar>
    ( Diff<T1, T2, Scalar>(x1.x, x2.x, x1.info()) );

}
template<class T, class Scalar, class Res, class IType = typename T::InfoType>
class OuterProd {

public:
  typedef IType InfoType;
private:
  Scalar y;
  T a;
  InfoType len;

public:

  OuterProd(Scalar yy, T aa, InfoType l) : y(yy), a(aa), len(l) { }

  Res operator[](int i) const { return y*a[i]; }
  InfoType info() const { return len; }

};

//------------------------------------------------------------------------------

template<class T, class S2>
inline
Expr<OuterProd<T, double, typename ProdRes<double, S2>::ResType>,
    typename ProdRes<double,S2>::ResType > operator*(double y, const Expr<T, S2> &x)
{

  return Expr<OuterProd<T, double, typename ProdRes<double,S2>::ResType>,
         typename ProdRes<double,S2>::ResType>
    ( OuterProd<T, double, typename ProdRes<double,S2>::ResType>(y, x.x, x.info()) );
}

//------------------------------------------------------------------------------

template<class T, class S2>
inline
Expr<OuterProd<T, complex<double>, typename ProdRes<complex<double>, S2>::ResType>,
    typename ProdRes<complex<double>,S2>::ResType > operator*(complex<double> y, const Expr<T, S2> &x)
{

  return Expr<OuterProd<T, complex<double>, typename ProdRes<complex<double>,S2>::ResType>,
         typename ProdRes<complex<double>,S2>::ResType>
    ( OuterProd<T, complex<double>, typename ProdRes<complex<double>,S2>::ResType>(y, x.x, x.info()) );
}

template<class T1, class T2, class Scalar, class OpApp, 
   class IType = typename T1::InfoType>
class BinOp {
  public:
    typedef IType InfoType;
  private:
    T1 a;
    T2 b;
  public:
    BinOp(T1 aa, T2 bb) : a(aa), b(bb) {}
    Scalar operator[](int i) const { return OpApp::apply(a[i],b[i]); }
    InfoType info() { return a.info(); }
};

#define binaryOpDec(op, applicator, tmpl, type) \
  template<class Scalar tmpl> \
  inline \
  Expr< \
  BinOp<const type &, const type &, Scalar, applicator, typename type::InfoType> \
    , Scalar > \
  operator op (const type &v1, const type &v2) { \
    return Expr<BinOp<const type &, const type &, Scalar, applicator, typename type::InfoType>,\
     bool > \
     ( BinOp<const type &, const type &, Scalar, applicator, typename type::InfoType>(v1,v2)); \
    }\
  template<class T, class Scalar tmpl> \
  inline \
  Expr< \
  BinOp<const type &, T, Scalar, applicator, typename type::InfoType>, \
    Scalar > \
  operator op (const type &v1, const Expr<T,Scalar> &v2) { \
    return Expr<BinOp<const type &, T, Scalar, applicator, typename type::InfoType>, \
     bool > \
     ( BinOp<const type &, T, Scalar, applicator, typename type::InfoType> \
    (v1,v2.x)); \
    } 
    /*
  template<class T, class Scalar tmpl> \
  inline \
  Expr< \
  BinOp<T, const type &, Scalar, \
    applicator>, Scalar > \
  operator op (const Expr<T,Scalar> &v1, const type &v2) { \
    return Expr<BinOp<const Expr<T, Scalar> &, const type &, Scalar, \
    applicator>, Scalar > \
     ( BinOp<const Expr<T, Scalar> &, const type &, Scalar, \
    applicator>(v1,v2)); \
    }
*/
#endif
