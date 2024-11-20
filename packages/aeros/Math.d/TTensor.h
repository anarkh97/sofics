#ifndef _TTENSOR_H_
#define _TTENSOR_H_
#include <Utils.d/NodeSpaceArray.h>
#include <cstdio>

template <class T, int d>
class CTD { // Chained Tensor Data
  T v[d];
 public:
  const T &operator[] (int i) const { return v[i]; }
  T &operator[] (int i) { return v[i]; }
};

template <class S, int d>
class SimpleTensor : public Tensor, public CTD<S,d> {
  public:
    SimpleTensor() {}
    void dblContractInto(const Tensor &, Tensor *) const ;
    SimpleTensor<S,d> &operator+=(const SimpleTensor<S,d> &b) {
      for(int i=0; i < d; ++i)
        (*this)[i] += b[i];
      return *this;
    }
    SimpleTensor<S,d> operator+(const SimpleTensor<S,d> &b) {
      SimpleTensor<S,d> res;
      for(int i=0; i < d; ++i)
        res[i] = (*this)[i] + b[i];
      return res;
    }
    SimpleTensor &operator = (const S &v) {
      for(int i=0; i < d; ++i)
        (*this)[i] = v;
      return *this;
    }
};

template <class S, int d>
class SymTensor :  public Tensor, public CTD<S,(d*(d+1))/2> {
  public:
    SymTensor<S,d> &operator= (const DoubleContraction &dc)
      {
       dc.assignTo(this);
       return *this;
      }
  void dblContractInto(const Tensor &, Tensor *) const;
};

template <class S, int d>
class SimpleTensor<SymTensor<S,2>,d> : public Tensor, public CTD<SymTensor<S,2>,d> {
  public:
    void dblContractInto(const Tensor &, Tensor *) const;
};

typedef SymTensor<double,2> Stress2D;
typedef SimpleTensor<Stress2D, 9> Stress2DDeriv9;
typedef SimpleTensor<SimpleTensor<double,3>,2> Grad2D;
typedef SimpleTensor<Grad2D, 9> Grad2DDeriv9;

inline
double dblContract(const SymTensor<double,2> &a, const SymTensor<double,2> &b) {
  return a[0]*b[0]+a[1]*b[1]+2*a[2]*b[2];
  }
  
inline
SymTensor<double,2> dblContract(const SymTensor<SymTensor<double,2>,2> &a, SymTensor<double,2> &b)
{
  SymTensor<double,2> res;
  res[0] = dblContract(a[0],b);
  res[1] = dblContract(a[1],b);
  res[2] = dblContract(a[2],b);
  return res; 
}


template<int n, int dim >
SimpleTensor<SymTensor<double,dim>, n>
dblContractTransp(SymTensor<SymTensor<double,dim>,dim> &a,
                  SimpleTensor<SymTensor<double,dim>, n> &b) {
  SimpleTensor<SymTensor<double,dim>, n> res;
  for(int i=0; i < n; ++i) {
     res[i] = dblContract(a,b[i]);
  }
  return res;
}


template <>
void
SymTensor<SymTensor<double,2>,2>::dblContractInto(const Tensor &b, Tensor *res) const;

template<class S, int n>
void
SimpleTensor<SymTensor<S,2>,n>::dblContractInto(const Tensor &b, Tensor *res) const {
  Tensor_d2s0 *resMat = dynamic_cast<Tensor_d2s0 *>(res);
  if(resMat != 0) {
    const SimpleTensor<SymTensor<S,2>, n> &bb = 
               dynamic_cast<const SimpleTensor<SymTensor<S,2>, n> &>(b);
    for(int i = 0; i < n; ++i)
      for(int j = 0; j < n; ++j) {
        (*resMat)[i*n+j] = 
          (*this)[i][0]*bb[j][0]+(*this)[i][1]*bb[j][1]+2*(*this)[i][2]*bb[j][2];
      }
  } else {
    fprintf(stderr, "Tried to double contract a Sym 3th order tensor into something unknown\n");
  }   
}

template <class S, int n>
SimpleTensor<S,n> operator||(const SymTensor<S,2> &a, const SimpleTensor<SymTensor<S,2>,n> &b)
{
 SimpleTensor<S,n> res;
 for(int i=0; i < n; ++i) {
     res[i] = dblContract(a,b[i]);
  }
  return res;
}

template<class S, int n>
SimpleTensor<S,n> operator*(double a, const SimpleTensor<S,n> &b)
{
 SimpleTensor<S,n> res;
 for(int i=0; i < n; ++i) {
     res[i] = a*b[i];
 }
 return res;
}

template<class S, int n, int d>
SimpleTensor<SimpleTensor<S,n>,n> operator||(const SimpleTensor<SimpleTensor<SymTensor<S,d>, n>, n> &a,
                                             SymTensor<S,d> &b)
{
    SimpleTensor<SimpleTensor<S,n>,n> res;
    for(int i = 0; i < n; ++i)
      for(int j = 0; j < n; ++j)
        res[i][j] = dblContract(a[i][j], b);
    return res;
}

template <class S, int d>
void SymTensor<S,d>::dblContractInto(const Tensor &, Tensor *) const
{
  fprintf(stderr, "ERROR: Calling dblContractInto where it is not defined for SymTensor\n");
}
template <class S, int d>
void SimpleTensor<S,d>::dblContractInto(const Tensor &, Tensor *) const
{
  fprintf(stderr, "ERROR: Calling dblContractInto where it is not defined for SimpleTensor\n");
}
#endif
