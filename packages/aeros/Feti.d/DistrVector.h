#ifndef _DISTRVECTOR_H_
#define _DISTRVECTOR_H_

#include <cstdio>
#include <iostream>
#include <Driver.d/Communicator.h>
#include <Utils.d/MyComplex.h>

#include <Eigen/Dense>

#include <Math.d/Vector.h>

struct DistrInfo {

   int len;
   union {
     int numDom;
     int numLocSub;
   };
   union {
     int *domLen;
     int *subLen;
   };
   int *subOffset;
   // For parallel operations grouped by thread:
   int numLocThreads;
   int *threadOffset;
   int *threadLen;
   bool *masterFlag;
         
   FSCommunicator *com;
   
   DistrInfo(int i);
   DistrInfo() { initialize(); };
   ~DistrInfo(); 
   void setMasterFlag();
   void setMasterFlag(bool *_masterFlag) { if(masterFlag) delete [] masterFlag; masterFlag = _masterFlag; }
   void computeOffsets();
   int totLen() const { return len; }
   int masterLen() const;
   int *getMasterFlag(int i) const { return 0; }
   void recomputeOffsets();
   bool operator==(const DistrInfo& other) const;
   bool operator!=(const DistrInfo& other) const;
 private:
   void initialize();
};

//------------------------------------------------------------------------------

#include <Math.d/Expr.h>

template<class Scalar> class GenPartialDistrVector;

template <typename Scalar>
using VecRef = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>;
template <typename Scalar>
using ConstVecRef = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>;

template<class Scalar>
class GenDistrVector {
protected:
	bool myMemory;
	int len;		//!< \brief entire length of the vector
	int numDom;		//!< \brief number of domains
	Scalar *v;		//!< \brief entire vector data
	std::vector<Scalar *> subV;	//!< \brief pointers to each domains sub-vector
	std::vector<int> subVLen;	//!< \brief length of each domains sub-vector
	int nT; //!< \brief Number of threads
	std::vector<int> thLen;     //!< \brief lengths per thread
	std::vector<Scalar*> thV;       // each thread's subvector
	std::vector<int> subVOffset;
	std::vector<int> thOffset;
	bool *masterFlag;
	bool infoFlag;
	DistrInfo const * inf;
	Scalar *partial;
public:
    GenDistrVector() : myMemory(false), len(0), numDom(0), v(NULL), nT(0),
                       masterFlag(NULL),
                       infoFlag(false), inf(new DistrInfo), partial(NULL) {}
    GenDistrVector(const DistrInfo &dinfo);
    GenDistrVector(const GenDistrVector<Scalar> &v);
    GenDistrVector(const DistrInfo &dinfo, Scalar *, bool myMemory = true);
    virtual ~GenDistrVector();
    void initialize();
    void zero();
    void clean_up();
    int size() const { return len; }
    void resize(const DistrInfo &dinfo); // no-op if the sizes match, otherwise data is lost
    void conservativeResize(const DistrInfo &dinfo); // resizing with data preservation
    int num() const { return numDom; }
    Scalar &operator[](int i) { return v[i]; }
    Scalar operator[](int i) const { return v[i]; } 
    Scalar operator*(const GenDistrVector &) const;
    Scalar operator^(const GenDistrVector &) const;
    double norm();
    double infNorm();
    double sqNorm() const { return ScalarTypes::norm((*this) * (*this)); }
    GenDistrVector &operator=(const GenDistrVector<Scalar> &);
    GenDistrVector &operator=(Scalar c);
    GenDistrVector &operator*=(Scalar c);
    GenDistrVector &operator/=(Scalar c);
    GenDistrVector &operator+=(GenDistrVector<Scalar> &);
    GenDistrVector &operator+=(const GenDistrVector<Scalar> &);
    GenDistrVector &operator-=(GenDistrVector<Scalar> &);
    GenDistrVector &operator/=(GenDistrVector<Scalar> &);
    GenDistrVector &linAdd(GenDistrVector<Scalar> &);
    GenDistrVector &linAdd(Scalar, GenDistrVector<Scalar> &);
    GenDistrVector &linAdd_inv(Scalar, GenDistrVector<Scalar> &);
    GenDistrVector &linAdd(Scalar, GenPartialDistrVector<Scalar> &);
    GenDistrVector &linAdd(Scalar, GenDistrVector<Scalar> &, Scalar, GenDistrVector<Scalar> &);
    GenDistrVector &linC(const GenDistrVector<Scalar> &, Scalar);
    GenDistrVector &linC(const GenDistrVector<Scalar> &, Scalar, const GenDistrVector<Scalar> &);
    GenDistrVector &linC(Scalar, const GenDistrVector<Scalar> &, Scalar, const GenDistrVector<Scalar> &);
    GenDistrVector &linC(Scalar, const GenDistrVector<Scalar> &, Scalar, const GenDistrVector<Scalar> &, Scalar, const GenDistrVector<Scalar> &);
    GenDistrVector &swap(GenDistrVector<Scalar> &);
    template <class T>
      GenDistrVector &operator=(const Expr<T,Scalar> &);
    template <class T>
      GenDistrVector &operator+=(const Expr<T,Scalar> &);
    template <class T>
      GenDistrVector &operator-=(const Expr<T,Scalar> &);

    void updateBlock(int ii, Scalar c, GenDistrVector<Scalar> &) {std::cerr << "GenDistrVector::updateBlock not implemented" << std::endl;}
    void copyBlock(GenDistrVector<Scalar> &, int ii) {std::cerr << "GenDistrVector::copyBlock not implemented" << std::endl;}
    void copy(const Scalar *v) {std::cerr << "GenDistrVector::copy(const Scalar *v) not implemented" << std::endl;}
    void copy(const Scalar v)  {std::cerr << "GenDistrVector::copy(const Scalar v) not implemented" << std::endl;}
    void addBlockSqr(int ii, Scalar c, GenDistrVector<Scalar> &);
    void computeSqrt();
    void computeRealz(int ii, Scalar c, GenDistrVector<Scalar> &) {std::cerr << "GenDistrVector::computeRealz not implemented" << std::endl;}
    void setn(int _n) {};      
    GenDistrVector<Scalar>&  getBlock(int iblock) { std::cerr << "GenDistrVector::getBlock not implemented" << std::endl; 
                                                    return *(new GenDistrVector<Scalar>()); }
	///\brief Obtain a view to the part of the vector for one subdomain.
	VectorView<Scalar> subVec(int iSub) { return VectorView<Scalar>{subData(iSub), subLen(iSub), 1}; }
	///\copydoc
	VectorView<const Scalar> subVec(int iSub) const { return VectorView<const Scalar>{subData(iSub), subLen(iSub), 1}; }
    void negate();
    Scalar *data() const      { return v;          }
    Scalar *subData(int i)    { return subV[i];    }
    const Scalar *subData(int i) const { return subV[i];    }
    int subLen(int i) const     { return subVLen[i]; }
    int subOffset(int i) const  { return subVOffset[i]; }
    int numThreads() const      { return nT; }
    int threadLen(int i) const  { return thLen[i];   }
    int threadOffset(int i) const { return thOffset[i]; }
    Scalar *threadData(int i) { return thV[i];     }
	const Scalar *threadData(int i) const { return thV[i]; }
    bool *threadMasterFlag(int i) { return masterFlag +thOffset[i]; }
    bool *subMasterFlag(int i) { return masterFlag +subVOffset[i]; }
    Scalar ident();
    void print();
    void printNonZeroTerms();
    void printAll();
    void initRand();
    void doubleUp(int, Scalar, GenDistrVector<Scalar> * , GenDistrVector<Scalar> * , 
                  GenDistrVector<Scalar> *);
    void tripleUp(int, Scalar, GenDistrVector<Scalar> * , GenDistrVector<Scalar> * ,
                  GenDistrVector<Scalar> *, GenDistrVector<Scalar> * , GenDistrVector<Scalar> *);
    virtual Scalar getPartial(int iSub) { return 1.0; }

    Scalar sum() {
       Scalar x =0;
       for(int i = 0; i < len; ++i) x+= v[i];
       return x;
    }

   void scaleBlock(int k, Scalar s) { std::cerr << "Error : GenDistrVector::scaleBlock not implemented " << std::endl; } 

   typedef const DistrInfo &InfoType;
   const DistrInfo &info() const { return *inf; }
};

template<class Scalar>
Scalar
dot_ignore_master_flag(const GenDistrVector<Scalar> &, const GenDistrVector<Scalar> &);

template<class Scalar>
class GenStackDistVector : public GenDistrVector<Scalar> {
   public:
      GenStackDistVector(const DistrInfo &dinfo, Scalar *v)
               : GenDistrVector<Scalar>(dinfo, v, false) { }
      virtual ~GenStackDistVector() { }
};

template<class Scalar>
class GenPartialDistrVector : public GenDistrVector<Scalar> {
     // PJSA 1-22-07 distributed vector for which some subdomains have zero subvectors
     // can speed up some operations like * by skipping these subvectors
   public:
     GenPartialDistrVector(const DistrInfo &dinfo)
        : GenDistrVector<Scalar>(dinfo) { this->partial = new Scalar[this->numDom]; for(int i=0; i<this->numDom; ++i) this->partial[i] = 1.0; }
     virtual ~GenPartialDistrVector() { delete [] this->partial; }

     Scalar operator * (const GenDistrVector<Scalar> &) const;
     void computePartial();
     Scalar getPartial(int iSub) { return this->partial[iSub]; }
};

typedef GenDistrVector<double> DistrVector;
typedef GenDistrVector<DComplex> ComplexDistrVector;
typedef GenStackDistVector<double> StackDistVector;
typedef GenStackDistVector<DComplex> ComplexStackDistVector;

//-----------------------------------------------------------------------------
template<class T1, class T2, class Scalar>
Scalar operator,(const Expr<T1,Scalar> &,const Expr<T2,Scalar> &);
//-----------------------------------------------------------------------------
template<class T2, class Scalar>
Scalar operator,(const GenDistrVector<Scalar> &,const Expr<T2,Scalar> &);
//-----------------------------------------------------------------------------
template<class T1,class Scalar>
Scalar operator,(const Expr<T1,Scalar> &,const GenDistrVector<Scalar> &);
//-----------------------------------------------------------------------------
template<class Scalar>
Scalar operator,(const GenDistrVector<Scalar> &,const GenDistrVector<Scalar> &);
//------------------------------------------------------------------------------
template<class Scalar>
double norm(const GenDistrVector<Scalar> &);
//------------------------------------------------------------------------------
template<class T1, class Scalar>
double norm(const Expr<T1,Scalar,const DistrInfo&> &);
//------------------------------------------------------------------------------

template<class Scalar>
inline
Expr<
 Sum<const GenDistrVector<Scalar> &, const GenDistrVector<Scalar> &, Scalar,
    typename GenDistrVector<Scalar>::InfoType>
    , Scalar>
operator+(const GenDistrVector<Scalar> &v1, const GenDistrVector<Scalar> &v2)
{

  return Expr<
    Sum<const GenDistrVector<Scalar> &, const GenDistrVector<Scalar> &, Scalar,
     typename GenDistrVector<Scalar>::InfoType> , Scalar >
    ( Sum<const GenDistrVector<Scalar> &, const GenDistrVector<Scalar> &,
      Scalar, typename GenDistrVector<Scalar>::InfoType>(v1, v2, v1.info()) );

}

//------------------------------------------------------------------------------

template<class T, class Scalar>
inline
Expr<Sum<T, Scalar *, Scalar>, Scalar>
operator+(const Expr<T, Scalar> &x, const GenDistrVector<Scalar> &v)
{

  return Expr<Sum<T, Scalar *, Scalar>, Scalar>
    ( Sum<T, Scalar *, Scalar>(x.x, v.data(), v.info()) );

}

//------------------------------------------------------------------------------

template<class T, class Scalar>
inline
Expr<Sum<const GenDistrVector<Scalar> &, T, Scalar, 
   typename GenDistrVector<Scalar>::InfoType>, Scalar>
operator+(const GenDistrVector<Scalar> &v, const Expr<T, Scalar> &x)
{

  return Expr<Sum<const GenDistrVector<Scalar> &, T, Scalar, 
   typename GenDistrVector<Scalar>::InfoType>
   , Scalar>
    ( Sum<const GenDistrVector<Scalar> &, T, Scalar, 
   typename GenDistrVector<Scalar>::InfoType>(v.data(), x.x, v.info()) );

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
Expr<
 Diff<const GenDistrVector<Scalar> &, const GenDistrVector<Scalar> &, Scalar,
    typename GenDistrVector<Scalar>::InfoType>
    , Scalar>
operator-(const GenDistrVector<Scalar> &v1, const GenDistrVector<Scalar> &v2)
{

  return Expr<
    Diff<const GenDistrVector<Scalar> &, const GenDistrVector<Scalar> &, Scalar,
     typename GenDistrVector<Scalar>::InfoType> , Scalar >
    ( Diff<const GenDistrVector<Scalar> &, const GenDistrVector<Scalar> &,
      Scalar, typename GenDistrVector<Scalar>::InfoType>(v1, v2, v1.info()) );

}

//------------------------------------------------------------------------------

template<class T, class Scalar>
inline
Expr<Diff<T, const GenDistrVector<Scalar> &, Scalar>, Scalar>
operator-(const Expr<T, Scalar> &x, const GenDistrVector<Scalar> &v)
{

  return Expr<Diff<T, const GenDistrVector<Scalar> &, Scalar>, Scalar>
    ( Diff<T, const GenDistrVector<Scalar> &, Scalar>(x.x, v, v.info()) );

}

//------------------------------------------------------------------------------

template<class T, class Scalar>
inline
Expr<Diff<const GenDistrVector<Scalar> &, T, Scalar, 
   typename GenDistrVector<Scalar>::InfoType>, Scalar>
operator-(const GenDistrVector<Scalar> &v, const Expr<T, Scalar> &x)
{

  return Expr<Diff<const GenDistrVector<Scalar> &, T, Scalar, 
   typename GenDistrVector<Scalar>::InfoType>
   , Scalar>
    ( Diff<const GenDistrVector<Scalar> &, T, Scalar, 
   typename GenDistrVector<Scalar>::InfoType>(v.data(), x.x, v.info()) );

}

//------------------------------------------------------------------------------

template<class Scalar, class Res>
inline
Expr<OuterProd<const GenDistrVector<Res> &, Scalar, 
       typename ProdRes<Scalar,Res>::ResType, 
       typename GenDistrVector<Res>::InfoType>,
     typename ProdRes<Scalar,Res>::ResType> operator*(Scalar y, const
     GenDistrVector<Res> &v)
{

  return Expr<
     OuterProd<const GenDistrVector<Res> &, Scalar, 
       typename ProdRes<Scalar,Res>::ResType, 
       typename GenDistrVector<Res>::InfoType>,
         typename ProdRes<Scalar,Res>::ResType>
    ( OuterProd<const GenDistrVector<Res> &, Scalar, 
       typename ProdRes<Scalar,Res>::ResType,
       typename GenDistrVector<Res>::InfoType>(y,
    v, v.info()) );

}

class CompIsL {
  public:
    template <class A, class B>
    static bool apply(const A &a, const B&b) { return a < b; }
};

class CompIsLE {
  public:
    template <class A, class B>
    static bool apply(const A &a, const B&b) { return a <= b; }
};

class CompIsG {
  public:
    template <class A, class B>
    static bool apply(const A &a, const B&b) { return a > b; }
};

class CompIsGE {
  public:
    template <class A, class B>
    static bool apply(const A &a, const B&b) { return a >= b; }
};

class CompIsE {
  public:
    template <class A, class B>
    static bool apply(const A &a, const B&b) { return a == b; }
};

class BoolAnd {
  public:
    template <class A, class B>
    static bool apply(const A &a, const B&b) { return a && b; }
};

class BoolOr {
  public:
    template <class A, class B>
    static bool apply(const A &a, const B&b) { return a || b; }
};

//-----------------------------------------------------------------------------

template<class T, class IType = typename T::InfoType>
class BoolNot {
public:
  typedef IType InfoType;
private:
  T a;
  InfoType len;
public:
  BoolNot(T t, IType it) : a(t), len(it) {}
  bool operator[](int i) const { return !a[i]; }
  InfoType info() const { return len; }
};

#ifdef _TEMPLATE_FIX_
#include <Feti.d/DistrVector.C>
#endif

#endif
