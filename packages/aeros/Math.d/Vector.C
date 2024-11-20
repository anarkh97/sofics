#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <algorithm>

template<class Scalar>
GenVector<Scalar>::GenVector(const GenVector<Scalar> &v2)
{    
 myMemory = true;
 len = v2.len;

 if(v2.d) {
   d = new Scalar[len];
   copy(v2.d);
 }
 else
   d = 0;

 n = v2.n;
}

template<class Scalar>
GenVector<Scalar>::GenVector(Scalar *v, int l, bool _myMemory)
{
 len = l;
 myMemory = _myMemory;

 if(myMemory) {
   d = new Scalar[len];
   copy(v);
 } else {
   d = v;
 }
}

template<class Scalar>
GenVector<Scalar>::GenVector(int l, Scalar *v, bool _myMemory)
{
 len = l;
 myMemory = _myMemory;

 if(myMemory) {
   d = new Scalar[len];
   copy(v);
 } else {
   d = v;
 }
}

template<class Scalar>
GenVector<Scalar>::GenVector(int l, Scalar initialValue)
{
 myMemory = true;
 len = l;
 d   = new Scalar[len];
 copy(initialValue);
}

template<class Scalar>
void
GenVector<Scalar>::copy(const GenVector<Scalar> &v2)
{
 int i;
 for(i=0; i<len; ++i)
   d[i] = v2[i];

 n = v2.n;
}

template<class Scalar>
void
GenVector<Scalar>::copy(const Scalar *v2)
{
 int i;
 for(i=0; i<len; ++i)
   d[i] = v2[i];
}

template<class Scalar>
void
GenVector<Scalar>::copy(const Scalar value)
{
 int i;
 for(i=0; i<len; ++i)
   d[i] = value;
}

template<class Scalar>
void
GenVector<Scalar>::setData(const GenVector<Scalar> &v1)
{
 if(d) { delete [] d; }
 
 len = v1.len; 
 d = new Scalar[len];
 for(int i=0; i<len; ++i)
   d[i] = v1[i];
}

template<class Scalar>
void
GenVector<Scalar>::reset(int newlen, Scalar initialvalue)
{
 if(myMemory && d) { delete [] d; }

 len = newlen;
 d = new Scalar[len];
 myMemory = true;
 for(int i=0; i<len; ++i)
   d[i] = initialvalue;
}

template<class Scalar>
void
GenVector<Scalar>::resize(int newlen)
{
 if(len == newlen) return;
 if(myMemory && d) { delete [] d; }
 len = newlen;
 d = new Scalar[len];
 myMemory = true;
}

template<class Scalar>
void
GenVector<Scalar>::conservativeResize(int newlen)
{
 if(len == newlen) return;
 Scalar *newd = new Scalar[newlen];
 for(int i=0; i<std::min(len,newlen); ++i)
   newd[i] = d[i];
 if(newlen > len) {
   for(int i=len; i<newlen; ++i) newd[i] = 0;
 }
 if(myMemory && d) delete [] d;
 myMemory = true;
 len = newlen;
 d  = newd;
}
                                                                                                                   
template<class Scalar>
void
GenVector<Scalar>::putIn(Scalar *array, int position, int num)
{
  if (num>len) fprintf(stderr,"Incompatible length in GenVector<Scalar> putIn\n");
  for (int i=0; i<num; i++)
    array[position+i]=d[i];
}

template<class Scalar>
void
GenVector<Scalar>::getFrom(Scalar *array, int position, int numdata)
{
  if ( numdata > len )
    fprintf(stderr,"Incompatible length of GenVector<Scalar>s in GenVector<Scalar> getFrom\n");

  for (int i=0; i<numdata; i++)
    d[i]=array[position+i];
}

template<class Scalar>
void
GenVector<Scalar>::getDataFrom(Scalar *array, int num)
{
  if (num!=len) fprintf(stderr,"Incompatible length in GenVector<Scalar> getDataFrom\n");
  for (int i=0; i<num; i++)
    d[i]=array[i];
}

template<class Scalar>
void
GenVector<Scalar>::addDataFrom(double *array, int num)
{
  if (num!=len) fprintf(stderr,"Incompatible length in GenVector<Scalar> addDataFrom\n");
  for (int i=0; i<num; i++)
    d[i]+=array[i];
}

template<class Scalar>
GenVector<Scalar>
GenVector<Scalar>::operator+(const GenVector<Scalar> &v2)
{
 if(len != v2.len) {
   fprintf(stderr,"Incompatible length of GenVector<Scalar>s in GenVector<Scalar> addition\n");
   return GenVector<Scalar>(1);
 }

 GenVector<Scalar> res(len);

 int i;
 for(i=0; i<len; ++i)
   res[i] = (*this)[i] + v2[i];

 return res;
}

template<class Scalar>
GenVector<Scalar>
GenVector<Scalar>::operator-(const GenVector<Scalar> &v2)
{
 if(len != v2.len) {
   fprintf(stderr,"Incompatible length of GenVector<Scalar>s in GenVector<Scalar> subtraction\n");
   return GenVector<Scalar>(1);
 }
 GenVector<Scalar> res(len);

 int i;
 for(i=0; i<len; ++i)
   res[i] = (*this)[i] - v2[i];

 return res;
}

template<class Scalar>
void
GenVector<Scalar>::operator*=(const Scalar c)
{
 if(c == 1.0) return;

 int i;
 for(i=0; i< len; ++i)
    d[i] *= c;
}

template<class Scalar>
void
GenVector<Scalar>::operator/=(const Scalar c)
{
 if(c == 1.0) return;

 int i;
 for(i=0; i< len; ++i)
    d[i] /= c;
}

template<class Scalar>
void
GenVector<Scalar>::operator+=(const GenVector<Scalar> &v2)
{
 if(len != v2.len) {
   fprintf(stderr,"Incompatible length of GenVector<Scalar>s in operator +=.\n");
   exit(-1);
 }

 int i;
 for(i=0; i<len; ++i)
    d[i] += v2[i];
}

template<class Scalar>
void
GenVector<Scalar>::operator-=(const GenVector<Scalar> &v2)
{
 if(len != v2.len) {
   fprintf(stderr,"Incompatible length of GenVector<Scalar>s in operator -=.\n");
   exit(-1);
 }

 int i;
 for(i=0; i< len; ++i)
    d[i] -= v2[i];
}

template<class Scalar>
void
GenVector<Scalar>::operator/=(const GenVector<Scalar> &v2)
{
 if(len != v2.len) {
   fprintf(stderr,"Incompatible length of GenVector<Scalar>s in operator /=.\n");
   exit(-1);
 }

 int i;
 for(i=0; i< len; ++i)
   d[i] /= v2[i];

}

// Mathematical dot product
/*
                  T
      (x, y) = {x} {y} 
*/

template<class Scalar>
Scalar
GenVector<Scalar>::operator*(const GenVector<Scalar> &v2) const
{
 int i;
 if(len != v2.len) {
   fprintf(stderr,"Incompatible length of GenVector<Scalar>s in dot product (operator *) \n");
   exit(-1);
 }
 Scalar answer = 0.0;
 for(i=0; i < len; ++i)
   answer += d[i]*ScalarTypes::conj(v2.d[i]);
 return answer;
}

template<class Scalar>
Scalar
GenVector<Scalar>::operator^(const GenVector<Scalar> &v2)
{
 int i;
 if(len != v2.len) {
   fprintf(stderr,"Incompatible length of GenVector<Scalar>s in dot product (operator ^)\n");
   exit(-1);
 }
 Scalar answer = 0.0;
 for(i=0; i < len; ++i)
   answer += d[i]*ScalarTypes::conj(v2.d[i]);
 return answer;
}


template<class Scalar>
GenVector<Scalar>
operator*(const Scalar c,const GenVector<Scalar> &v)
{
 int sz = v.size();
 GenVector<Scalar> res(sz);
 int i;
 for(i=0; i<sz; ++i)
    res[i] = c*v[i];
 return res;
}


template<class Scalar>
GenVector<Scalar>
GenVector<Scalar>::operator/(const GenVector<Scalar> &v2)
{
 if(len != v2.len) {
   fprintf(stderr,"Incompatible length of GenVector<Scalar>s in GenVector<Scalar> operator /\n");
   return GenVector<Scalar>(1);
 }

 GenVector<Scalar> res(len);

 int i;
 for(i=0; i<len; ++i) {
   if(v2[i] == 0.0) {
     fprintf(stderr," *** ERROR: Division by zero in GenVector<Scalar> operator/ \n");
     break;
   }
   res[i] = (*this)[i]/v2[i];
 }
 return res;
}

template<class Scalar>
GenVector<Scalar> &
GenVector<Scalar>::operator=(const GenVector<Scalar> &v2)
{
 // check lengths of GenVector<Scalar>s first
 if(len != v2.len) {
   if(myMemory) delete [] d;
   d   = new Scalar[v2.len];
   len = v2.len;
   myMemory = true;
 }

 copy( v2.d );

 return *this;
}

template<class Scalar>
GenVector<Scalar> &
GenVector<Scalar>::operator=(const Scalar c)
{
 copy(c);
 return *this;
}

template<class Scalar>
GenVector<Scalar> &
GenVector<Scalar>::operator=(const Scalar *data)
{
 copy(data);
 return *this;
}

template<class Scalar>
void
GenVector<Scalar>::zero()
{
  for(int i=0; i < len; ++i) ScalarTypes::initScalar(d[i], 0.0, 0.0);
}

template<class Scalar>
void
GenVector<Scalar>::linC(const GenVector<Scalar> &v1, Scalar c)
{

// linear combination where v = c*v1
// where c = constant

 if(len != v1.len) {
   delete [] d;
   d   = new Scalar[v1.len];
   len = v1.len;
 }

 int i;
 for(i=0; i < len; ++i)
   d[i] = c*v1.d[i];

}

template<class Scalar>
void
GenVector<Scalar>::linC(Scalar c, const GenVector<Scalar> &v1)
{
// linear combination where v = c*v1
// where c = constant

 if(len != v1.len) {
   delete [] d;
   d   = new Scalar[v1.len];
   len = v1.len;
 }

 int i;
 for(i=0; i < len; ++i)
   d[i] = c*v1.d[i];

}

template<class Scalar>
void
GenVector<Scalar>::linC(const GenVector<Scalar> &v1, Scalar c, const GenVector<Scalar> &v2)
{
// linear combination: v1 + c*v2
// where c = constant

 if(v1.len != v2.len) {
   fprintf(stderr,"Incompatible length of GenVector<Scalar>s in linC\n");
   exit(-1);
 }
 if(len != v2.len) {
   delete [] d;
   d = new Scalar[v2.len];
   len = v2.len;
 }
 int i;
 for(i=0; i < len; ++i)
   d[i] = v1.d[i] + c*v2.d[i];
}

template<class Scalar>
void
GenVector<Scalar>::linC(Scalar c1, const GenVector<Scalar> &v1, Scalar c2, const GenVector<Scalar> &v2)
{

// linear combination: c1*v1 + c2*v2
// where c1 and c2 are constants

 if(v1.len != v2.len) {
   fprintf(stderr,"Incompatible length of GenVector<Scalar>s in linC\n");
   exit(-1);
 }
 if(len != v2.len) {
   delete [] d;
   d = new Scalar[v2.len];
   len = v2.len;
 }

 int i;
 for(i=0; i < len; ++i)
   d[i] = c1*v1.d[i] + c2*v2.d[i];
}

template<class Scalar>
void
GenVector<Scalar>::linC(Scalar c1, const GenVector<Scalar> &v1, Scalar c2, const GenVector<Scalar> &v2,
                        Scalar c3, const GenVector<Scalar> &v3)
{
  for(int i = 0; i < len; ++i)
    d[i] = c1*v1.d[i] + c2*v2.d[i] + c3*v3.d[i];
}

template<class Scalar>
void
GenVector<Scalar>::swap(GenVector<Scalar> &v1)
{
 std::swap(len, v1.len);
 std::swap(d, v1.d);
 std::swap(myMemory, v1.myMemory);
 std::swap(n, v1.n);
}

template<class Scalar>
void 
GenVector<Scalar>::updateBlock(int ii, Scalar c, GenVector<Scalar> &v)
{
 for (int i=0; i<n; ++i)   d[ii*n+i]=d[ii*n+i]+c*v[i];
}


template<class Scalar>
void
GenVector<Scalar>::copyBlock(GenVector<Scalar> &v, int ii) 
{
 for (int i=0; i<n; ++i) d[i] = v[ii*n+i];
}


template<class Scalar>
void 
GenVector<Scalar>::addBlockSqr(int ii, Scalar c, GenVector<Scalar> &v)
{
 for (int i=0; i<n; ++i) d[i] = d[i] + c*ScalarTypes::sqNorm(v[ii*n+i]);
}


template<class Scalar>
void
GenVector<Scalar>::computeSqrt()
{
 for (int i=0; i<n; ++i) d[i] = ScalarTypes::sqrt(d[i]);
}


template<class Scalar>
void 
GenVector<Scalar>::computeRealz(int ii, Scalar c, GenVector<Scalar> &v)
{
 for (int i=0; i<n; ++i)   d[i]=d[i]+c*v[ii*n+i];
}


template<class Scalar>
void
GenVector<Scalar>::print(const char *msg, const char* msg2)
{
 if(msg) std::cerr << msg << " ";
 for(int i=0; i<len; ++i) std::cerr << std::setprecision(10) << d[i] << " ";
 std::cerr << std::endl;
}

template<class Scalar>
GenVector<Scalar>
GenVector<Scalar>::cross(GenVector<Scalar> &v2)
{
 if(len != v2.len) {
   fprintf(stderr,"Incompatible length of GenVector<Scalar>s in GenVector<Scalar> product.\n");
   return GenVector<Scalar>(1);
 }
 GenVector<Scalar> v3(len);

 v3[0] = d[1]*v2[2] - d[2]*v2[1];
 v3[1] = d[2]*v2[0] - d[0]*v2[2];
 v3[2] = d[0]*v2[1] - d[1]*v2[0];

 return v3;
}

template<class Scalar>
double
GenVector<Scalar>::magnitude()
{
 double res = 0.0;
 for(int i = 0; i < len; ++i)
   res += ScalarTypes::sqNorm(d[i]);
 return sqrt( res );
}

template<class Scalar>
double
GenVector<Scalar>::absMax()
{
 double max = ScalarTypes::norm((*this)[0]);
 int i;
 for(i = 0; i<len; ++i)
     if(ScalarTypes::norm((*this)[i]) > max)
        max = ScalarTypes::norm((*this)[i]);

 return max;

}

template<class Scalar>
double
GenVector<Scalar>::squareNorm() const
{
 double res = 0.0;
 for(int i = 0; i < len; ++i)
   res += ScalarTypes::sqNorm(d[i]);
 return res;
}

template<class Scalar>
double
GenVector<Scalar>::norm() const
{
 return sqrt( squareNorm() );
}

template<class Scalar>
void
GenVector<Scalar>::linAdd(Scalar c, const GenVector<Scalar> &v2)
{
 // linear addition: v1 = v1 + c*v2 
 if(len != v2.len) {
   fprintf(stderr,"Incompatible length of GenVector<Scalar>s in linAdd\n");
   delete [] d;
   d   = new Scalar[v2.len];
   len = v2.len;
 }

 int i;
 for(i=0; i < len; ++i)
   d[i] += c*v2.d[i]; 
}

template<class Scalar>
void
GenVector<Scalar>::linAdd(Scalar c1, const GenVector<Scalar> &v1, Scalar c2, const GenVector<Scalar> &v2)
{
// linear addition:  this + c1*v1 + c2*v2
// where c1 and c2 are constants

 if(v1.len != v2.len) {
   fprintf(stderr,"Incompatible length of GenVector<Scalar>s in linAdd\n");
   exit(-1);
 }
 if(len != v2.len) {
   delete [] d;
   d   = new Scalar[v2.len];
   len = v2.len;
 }
 int i;
 for(i=0; i < len; ++i)
   d[i] += c1*v1.d[i] + c2*v2.d[i];
}

template<class Scalar>
void
GenVector<Scalar>::mult(GenVector<Scalar> &v)
{
 int i;
 for(i=0; i < len; ++i)
   d[i] *= v[i];
}

template<class Scalar>
void
GenVector<Scalar>::diff(GenVector<Scalar> &v1, GenVector<Scalar> &v2)
{
 // Check lengths of v1 and v2
 if(v1.len != v2.len) {
   fprintf(stderr,"Incompatible length of GenVector<Scalar>s in difference\n");
   exit(-1);
 }

 if(len != v2.len) {
   delete [] d;
   d   = new Scalar[v2.len];
   len = v2.len;
 }

 int i;
 for(i=0; i<len; ++i)
   d[i] = v1[i] - v2[i];
}

template<class Scalar>
GenVector<Scalar>::GenVector(const GenVector<Scalar> &v, int nr, int sr)
{
 myMemory = true;
 len = nr;
 d = new Scalar[len] ;

 int i;
 for(i=0; i < len; ++i)
    d[i] = v[i+sr] ;
}

template<class Scalar>
void
GenVector<Scalar>::add(GenVector<Scalar> &vec, int fRow)
{
  int size = vec.size();

  int irow;
  for(irow = 0; irow < size; ++irow) {
    d[fRow+irow] += vec[irow];
  }
}

template<class Scalar>
void
GenVector<Scalar>::add(Scalar *v)
{
  int irow;
  for(irow = 0; irow < len; ++irow) {
    d[irow] += v[irow];
  }
}

template<class Scalar>
void
GenVector<Scalar>::add(gsl::span<Scalar> v)
{
    int irow;
    for(irow = 0; irow < len; ++irow) {
        d[irow] += v[irow];
    }
}

template<class Scalar>
void
GenVector<Scalar>::add(GenVector<Scalar> &vec, int *rows)
{
 int size = vec.size();

  int irow;
  for(irow = 0; irow < size; ++irow) {
    d[rows[irow]] += vec[irow];
  }
}

template<class Scalar>
void
GenVector<Scalar>::insertData(Scalar *v)
{
  int irow;
  for(irow = 0; irow < len; ++irow) {
    d[irow] = v[irow];
  }
}

template<class Scalar>
void
GenVector<Scalar>::insertData(gsl::span<Scalar> v)
{
    for(int irow = 0; irow < len; ++irow) {
        d[irow] = v[irow];
    }
}


template<class Scalar>
void
GenVector<Scalar>::subtract(GenVector<Scalar> &vec, int fRow)
{
  int size = vec.size();

  int irow;
  for(irow = 0; irow < size; ++irow) {
    d[fRow+irow] -= vec[irow];
  }
}



template<class Scalar>
GenVector<Scalar>::GenVector(const GenVector<Scalar> &v, int nr, int *rows)
{
 myMemory = true;
 len = nr;
 d = new Scalar[len] ;

 int i;
 for(i=0; i < len; ++i)
    d[i] = v[rows[i]] ;
}


template<class Scalar>
void
GenVector<Scalar>::zeroAll()
{
  for(int i=0; i < len; ++i) ScalarTypes::initScalar(d[i], 0.0, 0.0);
}

template<class Scalar>
void
GenVector<Scalar>::vectorToTensor(Tensor_d1s0 &t) {
if ( len != t.getSize())
 {
	fprintf(stderr,"Wrong VectorToTensor Conversion, dimension mismatch\n");
	exit(-1);
 }

for (int i = 0; i<len; i++)
  t[i]=d[i]; 
}


template<class Scalar>
GenVector<Scalar>&
GenVector<Scalar>::getBlock(int k) 
{ 
 GenStackVector<Scalar> *v = new GenStackVector<Scalar>(n, d + n*k); 
 return *v; 
}  
