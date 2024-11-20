#include <iostream>
#include <cstdio>
#include <cassert>
#include <algorithm>

#include <Math.d/FullSquareMatrix.h>

/*
extern "C"      {
   void _FORTRAN(dsyev)(const char &, const char &, const int &,
                        double *, const int &, double *, double *,
                        const int &, int &);
   //use zheev if complex matrix
   // for more functions see dsyevr and zheevr
}

inline void Tdsyev(const char &a, const char &b, const int &c,
                   double *d, const int &e, double *f, double *g,
                   const int &h, int &i)
{ _FORTRAN(dsyev)(a,b,c,d,e,f,g,h,i); }

*/
#include <Utils.d/MyComplex.h>
#include <vector>

template<class Scalar>
GenFullSquareMatrix<Scalar>::GenFullSquareMatrix()
{
 size = 0; value = 0; myval=1; length = 0;
}


template<class Scalar>
GenFullSquareMatrix<Scalar>::GenFullSquareMatrix(int i, Scalar *l)
{
  size = i;
  length = size*size;
  if(l) {
    value=l;
    myval=0;
  }
  else {
    value=new Scalar[size*size];
    myval=1;
  }
}

template<class Scalar>
GenFullSquareMatrix<Scalar>::GenFullSquareMatrix(GenFullSquareMatrix<Scalar> &m, Scalar s)
{
  size = m.dim();
  length = size*size;
  value= new Scalar[length];
  myval=1;
  int p = 0;
  for(int i=0; i<m.numRow(); ++i)
    for(int j=0; j<m.numCol(); ++j)
      value[p++] = m[i][j]*s;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::copy(const GenFullSquareMatrix<Scalar> &m)
{
  if(size != m.dim()) {
    size = m.dim();
    length = size*size;
    if(myval && value) delete [] value;
    value= new Scalar[length];
    myval=1;
  }
  int p = 0;
  for(int i=0; i<m.numRow(); ++i)
    for(int j=0; j<m.numCol(); ++j)
      value[p++] = m[i][j];
}

template<class Scalar>
GenFullSquareMatrix<Scalar>::~GenFullSquareMatrix()
{
  if(value && myval) delete [] value;
  value=0;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::setSize(int i)
{
  if(value && myval) delete [] value;
  size = i;
  length = size*size;
  if(myval) value = new Scalar[length];
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::changeSize(int i, int numMax)
{
 if (i>numMax) fprintf(stderr," May be a Problem in GenFullSquareMatrix<Scalar>::changeSize ");
 else {
   size = i;
   length = size*size;
 }
}


template<class Scalar>
void
GenFullSquareMatrix<Scalar>::reSize(int newSize)
{
  if (newSize == size) {
    return;
  }
  
  // Create updated data
  int newLength = newSize * newSize;
  Scalar * newValue = new Scalar[newLength];
  int copySize = std::min(newSize, size);

  for (int i = 0; i < copySize; ++i) {
    const Scalar * rowBegin = value + i * size;
    Scalar * newRowBegin = newValue + i * newSize;
    std::copy(rowBegin, rowBegin + copySize, newRowBegin);
  } 

  // Commit changes
  size = newSize;
  length = newLength;
  if (myval) {
    delete[] value;
  }
  value = newValue;
  myval = 1;
}

//template<class Scalar>
//GenFullSquareMatrix<Scalar>
//GenFullSquareMatrix<Scalar>::operator + (const GenFullSquareMatrix<Scalar> &M2)
//{
//  GenFullSquareMatrix<Scalar> res(size);
//  for(int i = 0; i < length; ++i)
//    res.setData(i, value[i] + M2.value[i]);
//  return res;
//}
//
//template<class Scalar>
//GenFullSquareMatrix<Scalar>
//GenFullSquareMatrix<Scalar>::operator - (const GenFullSquareMatrix<Scalar> &M2)
//{
//  GenFullSquareMatrix<Scalar> res(size);
//  for(int i = 0; i < length; ++i)
//    res.setData(i, value[i] - M2.value[i]);
//  return res;
//}
//
//template<class Scalar>
//GenFullSquareMatrix<Scalar>
//GenFullSquareMatrix<Scalar>::operator * (Scalar s)
//{
//  GenFullSquareMatrix<Scalar> res(size);
//  for(int i = 0; i < length; ++i)
//    res.setData(i, s*value[i]);
//  return res;
//}
//
//template<class Scalar>
//GenFullSquareMatrix<Scalar>
//GenFullSquareMatrix<Scalar>::operator / (Scalar s)
//{
//  GenFullSquareMatrix<Scalar> res(size);
//  for(int i = 0; i < length; ++i)
//    res.setData(i, value[i]/s);
//  return res;
//}

template<class Scalar>
GenFullSquareMatrix<Scalar> &
GenFullSquareMatrix<Scalar>::operator += (const GenFullSquareMatrix<Scalar> &M2)
{
  for(int i = 0; i < length; ++i)
    value[i] += M2.value[i];
  return *this;
}

template<class Scalar>
GenFullSquareMatrix<Scalar> &
GenFullSquareMatrix<Scalar>::operator -= (const GenFullSquareMatrix<Scalar> &M2)
{
  for(int i = 0; i < length; ++i)
    value[i] -= M2.value[i];
  return *this;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::multiply(GenFullSquareMatrix<Scalar> &M2, GenFullSquareMatrix<Scalar> &result)
{
  for(int i = 0; i < size; ++i) {
    for(int j = 0; j < size; ++j) {
      result[i][j] = 0.0;
      for(int k = 0; k < size; ++k) {
        result[i][j] += (*this)[i][k]*M2[k][j];
      }
    }
  }
}


template<class Scalar>
GenFullSquareMatrix<Scalar> &
GenFullSquareMatrix<Scalar>::operator += (const Tensor_d2s0 &tens)
{
  for(int i = 0; i < length; ++i)
    value[i] += tens[i];
  return *this;
}

template<class Scalar>
GenFullSquareMatrix<Scalar> &
GenFullSquareMatrix<Scalar>::operator *= (Scalar s)
{
  for(int i = 0; i < length; ++i)
    value[i] *= s;
  return *this;
}

template<class Scalar>
GenFullSquareMatrix<Scalar> &
GenFullSquareMatrix<Scalar>::operator /= (Scalar s)
{
  for(int i = 0; i < length; ++i)
    value[i] /= s;
  return *this;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::zero()
{
  for(int i=0; i<length; ++i)
    value[i] = 0.0;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::unitary()
{
  for(int i=0; i<length; ++i)
    value[i] = 1.0;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::unitaryDiag()
{
  int i;
  for(i=0; i<length; ++i)
     value[i] = 0.0;
  for(i=0; i<size; ++i)
    (*this)[i][i] = 1.0;
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::symmetrize()
{
 int i, j;
 for(i=0; i<size; ++i)
   for(j=0;  j<i; ++j)
     (*this)[j][i] = (*this)[i][j] = 0.5*((*this)[j][i] + (*this)[i][j]);
}

// ... THIS FUNCTION PRINTS A FULLSQUARE MATRIX
// ... IN MATHEMATICA MATRIX INPUT FORM
template<class Scalar>
void
GenFullSquareMatrix<Scalar>::print(const char *msg, const char *msg2)
{
 int mathematica = 0;

 if(mathematica) {
   if(*msg) printf("%s\n",msg);
   int i,j;
   printf("{\n");
   for(i=0; i<size; ++i) {
     printf("{ ");
     for(j=0; j<size; ++j) {
       if(j==size-1)
         printf("%e ",(*this)[i][j]);
       else
         printf("%e, ",(*this)[i][j]);
     }
     if(i==size-1)
       printf("}");
     else
       printf("},");
     printf("\n");
   }
   printf("}\n");
 } else {
   if(*msg) fprintf(stderr,"%s\n",msg);
   int i,j;
   for(i=0; i<size; ++i) {
     for(j=0; j<size; ++j)
       //fprintf(stderr,"%s(%d,%d)=%3.2e, ",msg2,i+1,j+1,(*this)[i][j]);
       fprintf(stderr,"%20.16e,",(*this)[i][j]);
     fprintf(stderr,"\n");
   }
 }
}

template<>
inline
void
GenFullSquareMatrix<complex<double> >::print(const char*, const char*) {
  fprintf(stderr, "GenFullSquareMatrix<complex<double> >::print not implemented\n");
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::printDiagonals()
{
  int i;
  for(i=0; i<size; ++i)
    fprintf(stderr,"% e\n",std::real((*this)[i][i]));
}


template<class Scalar>
void
GenFullSquareMatrix<Scalar>::copy(Scalar *d)
{
  // WARNING : No check is performed on the size of d
  for(int i=0; i<length; ++i)
    value[i] = d[i];
}
template<class Scalar>
void
GenFullSquareMatrix<Scalar>::add(GenFullSquareMatrix<Scalar> &m, int *rc)
{
  for(int i=0; i<m.dim(); ++i)
    for(int j=0; j<m.dim(); ++j)
      (*this)[rc[i]][rc[j]] += m[i][j];
}

template<class Scalar>
void
GenFullSquareMatrix<Scalar>::multiply(GenFullSquareMatrix<Scalar> &res, double d)
{
  for(int i = 0; i < length; ++i)
    res.value[i] = value[i] * d;
}

template<>
void
GenFullSquareMatrix<double>::eigenVals(double* eigenV)
{
	//Scalar work[6*size];
	std::vector<double> copyValue{ value, value+size*size };
	//use the lapack routine dsyev
	char a = 'N';//'N' for eigenVlaues alone and 'V' for eigenvalues + eigenvectors
	char b = 'U';//store upper triangle
	int c = size;//order of the matrix
	double *d = copyValue.data();
	int e = size;
	double *f = eigenV;
	double *g = new double[6*size];
	int h = 6*size;
	int i = 0;
//  _FORTRAN(dsyev)(a, b, c, d, e, f, g, h, i);
	//other function can be used: zheev for complex or dsyevr & zheevr for faster methods
	Tdsyev(a, b, c, d, e, f, g, h, i);
	delete[] g;
}

template<>
void
GenFullSquareMatrix<double>::eigenV(double* eigenV)
{
  // value will change to store the eigenvectors
  //use the lapack routine dsyev
  char a = 'V';//'N' for eigenVlaues alone and 'V' for eigenvalues + eigenvectors
  char b = 'U';//store upper triangle
  int c = size;//order of the matrix
  double *d = value;
  int e = size;
  double *f = eigenV;
  double *g = new double[6*size];
  int h = 6*size;
  int i = 0;
//  _FORTRAN(dsyev)(a, b, c, d, e, f, g, h, i);
  //other function can be used: zheev for complex or dsyevr & zheevr for faster methods
  Tdsyev(a, b, c, d, e, f, g, h, i);
  delete[] g;
}

/*
template<class Scalar>
void
GenFullSquareMatrix<Scalar>::invert(GenFullSquareMatrix<Scalar> GFSM)
{
//see sommerElement.C
}
*/

#include <Utils.d/linkfc.h>
#include "BLAS.h"

template<> template<>
inline void GenFullSquareMatrix<DComplex>::multiply(GenVector<DComplex>& a, GenVector<DComplex>& b, double c, GenFullSquareMatrix<DComplex>::TransposeFlag transposed)
{
  DComplex alpha = c;
  DComplex beta  = 1.0;
  int one = 1;
  char trans = (transposed == TRANSPOSED) ? 'N' : 'T';
  _FORTRAN(zgemv)(trans, size, size, alpha, value, size, a.data(), one, beta, b.data(), one);
  return;
}

template<> template<>
inline void GenFullSquareMatrix<DComplex>::multiply(GenVector<DComplex>& a, GenVector<DComplex>& b, DComplex c, GenFullSquareMatrix<DComplex>::TransposeFlag transposed)
{
  DComplex beta  = 1.0;
  int one = 1;
  char trans = (transposed == TRANSPOSED) ? 'N' : 'T';
  _FORTRAN(zgemv)(trans, size, size, c, value, size, a.data(), one, beta, b.data(), one);
  return;
}

template<> template<>
inline void GenFullSquareMatrix<double>::multiply(GenVector<double>& a, GenVector<double>& b, double c, GenFullSquareMatrix<double>::TransposeFlag transposed)
{
  double beta  = 1.0;
  int one = 1;
  char trans = (transposed == TRANSPOSED) ? 'N' : 'T';
  _FORTRAN(dgemv)(trans, size, size, c, value, size, a.data(), one, beta, b.data(), one);
  return;
}

template<> template<>
inline void GenFullSquareMatrix<double>::multiply(GenVector<DComplex>& a, GenVector<DComplex>& b, double c, GenFullSquareMatrix<double>::TransposeFlag transposed)
{
  double beta  = 1.0;
  int two = 2;
  char trans = (transposed == TRANSPOSED) ? 'N' : 'T';
  _FORTRAN(dgemv)(trans, size, size, c, value, size,
		  reinterpret_cast<double*>(a.data()), two, beta,
		  reinterpret_cast<double*>(b.data()), two);
  _FORTRAN(dgemv)(trans, size, size, c, value, size,
		  reinterpret_cast<double*>(a.data())+1, two, beta,
		  reinterpret_cast<double*>(b.data())+1, two);
  return;
}

template<> template<>
inline void GenFullSquareMatrix<double>::multiply(GenVector<double>& a, GenVector<double>& b, DComplex c, GenFullSquareMatrix<double>::TransposeFlag transposed)
{
  assert(0);
}

template <typename Scalar>
inline
void
GenFullSquareMatrix<Scalar>::multAdd(const GenVector<Scalar> &a, GenVector<Scalar> &b) const {
  const_cast<GenFullSquareMatrix<Scalar> *>(this)->multiply(const_cast<GenVector<Scalar> &>(a), b, 1.0, NORMAL);
}

template <typename Scalar>
inline
void
GenFullSquareMatrix<Scalar>::mult(const GenVector<Scalar> &a, GenVector<Scalar> &b) const {
  b.zero();
  this->multAdd(a, b);
}

template <typename Scalar>
void
GenFullSquareMatrix<Scalar>::linAdd(const Scalar &coef, const GenFullSquareMatrix<Scalar> &other) {
  assert(this->length == other.length);
  for (int i = 0; i < this->length; ++i) {
    this->value[i] += coef * other.value[i];
  }
}

template class GenFullSquareMatrix<double>;
template class GenFullSquareMatrix<std::complex<double>>;
