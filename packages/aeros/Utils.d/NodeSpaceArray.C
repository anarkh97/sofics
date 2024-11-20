#include <Utils.d/NodeSpaceArray.h>
#include <cstdio>
#include <cstdlib>

#define UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
#define NO_DYNAMIC_CAST_IN_NODESPACE_ARRAY_C

Contraction
Tensor::operator | (const Tensor &b) const
{
  return Contraction(*this, b);
}

DoubleContraction 
Tensor::operator || (const Tensor &b) const
{
  return DoubleContraction(*this, b);
}

void 
Tensor::dblContractInto(const Tensor &, Tensor *) const 
{
  fprintf(stderr," It seems that the operands don't fit the operator ||...   \n");
  exit(-1);
}

void 
Tensor::splContractInto(const Tensor &, Tensor *) const
{
  fprintf(stderr," It seems that the operands don't fit the operator |...   \n");
  exit(-1);
}

void
Tensor::dblContractWith(const Tensor_d0s4_Ss12s34 &, Tensor *) const
{
  fprintf(stderr," It seems that the operands don't fit the operator ||...   \n");
  exit(-1);
}

void
Tensor::dblContractWith(const Tensor_d0s4_Ss12s34_diag &, Tensor *) const
{
  fprintf(stderr," It seems that the operands don't fit the operator ||...   \n");
  exit(-1);
}

void
Tensor::dblContractWith(const Tensor_d0s4 &, Tensor *) const
{
  fprintf(stderr," It seems that the operands don't fit the operator ||...   \n");
  exit(-1);
}

void
Tensor::splContractWith(const Tensor_d0s2 &, Tensor *) const
{
  fprintf(stderr," It seems that the operands don't fit the operator |...   \n");
  exit(-1);
}

void
Tensor::setZero()
{
  fprintf(stderr," Tensor::setZero() is not implemented ...\n");
  exit(-1);
}

Tensor_d0s2
Tensor_d0s2::operator + (const Tensor_d0s2 &tens) const
{
  Tensor_d0s2 t;
  for (int i = 0; i < 9; i++)
    t[i] = v[i] + tens[i];
  return t;
}

Tensor_d0s2
Tensor_d0s2::operator * (double scal)
{
  Tensor_d0s2 t;
  for (int i = 0; i < 9; i++)
    t[i] = v[i] * scal;
  return t;
}

Tensor_d0s2
operator * (double scal, const Tensor_d0s2 &tens)
{
  Tensor_d0s2 t;
  for (int i = 0; i < 9; i++)
    t[i] = tens[i] * scal;
  return t;
}

Tensor_d0s2 &
Tensor_d0s2::operator = (const Tensor_d0s2 &t)
{
  for(int i = 0; i < 9; ++i)
    v[i] = t.v[i];
  return *this;
}

Tensor_d0s2 &
Tensor_d0s2::operator = (const Tensor_d0s2_Ss12 &t)
{
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      v[3*i+j] = t(i,j);
  return *this;
}

/*
Tensor_d0s2
Tensor_d0s2::operator | (const Tensor_d0s2 &tens) const 
{
  Tensor_d0s2 t;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        t[3*i+j] += v[3*i+k] * tens[3*k+j];
  return t;
}

Tensor_d1s2_full
Tensor_d0s2::operator | (const Tensor_d1s2_full &tens) const
{
  int size = tens.getSize();
  Tensor_d1s2_full t(size);
  for (int m = 0; m < size; m++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
          t[m][3*i+j] += v[3*i+k] * tens[m][3*k+j];
  return t;
}

Tensor_d1s2_full
Tensor_d0s2::operator | (const Tensor_d1s2_sparse &tens) const
{
  int size = tens.getSize();
  Tensor_d1s2_full t(size);
  for (int n = 0; n < size; n++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        t[n][3*i+j] = v[3*i+(n%3)] * tens[(n-(n%3))/3][3*(n%3)+j];
  return t;
}
*/

void
Tensor_d0s2::dblContractWith(const Tensor_d0s4 &tens, Tensor *result) const
{
  Tensor_d0s2 &t = static_cast<Tensor_d0s2 &>(*result);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      t[3*i+j] = 0.0;
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          t[3*i+j] +=  tens[i*27+j*9+k*3+l] * v[3*k+l];
    }
}

void
Tensor_d0s2::splContractWith(const Tensor_d0s2 &tens, Tensor *result) const
{
  Tensor_d0s2 &t = static_cast<Tensor_d0s2 &>(*result);
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  t[0] = tens[0]*v[0] + tens[1]*v[3] + tens[2]*v[6];
  t[1] = tens[0]*v[1] + tens[1]*v[4] + tens[2]*v[7];
  t[2] = tens[0]*v[2] + tens[1]*v[5] + tens[2]*v[8];
  t[3] = tens[3]*v[0] + tens[4]*v[3] + tens[5]*v[6];
  t[4] = tens[3]*v[1] + tens[4]*v[4] + tens[5]*v[7];
  t[5] = tens[3]*v[2] + tens[4]*v[5] + tens[5]*v[8];
  t[6] = tens[6]*v[0] + tens[7]*v[3] + tens[8]*v[6];
  t[7] = tens[6]*v[1] + tens[7]*v[4] + tens[8]*v[7];
  t[8] = tens[6]*v[2] + tens[7]*v[5] + tens[8]*v[8];
#else
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      t[3*i+j] = 0.0;
      for (int k = 0; k < 3; k++)
        t[3*i+j] += tens[3*i+k] * v[3*k+j];
    }
#endif
}

void
Tensor_d0s2::splContractInto(const Tensor &b, Tensor *result) const
{
#ifdef NO_DYNAMIC_CAST_IN_NODESPACE_ARRAY_C
  b.splContractWith(*this, result);
#else
  const Tensor_d0s2 *tens = dynamic_cast<const Tensor_d0s2 *>(&b);
  if(tens)
  {
    Tensor_d0s2 &t = static_cast<Tensor_d0s2 &>(*result);  
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        t[3*i+j] = 0.0;
        for (int k = 0; k < 3; k++)
          t[3*i+j] += v[3*i+k] * (*tens)[3*k+j];
      }
  }
  else
  {
    const Tensor_d1s2_sparse *tens1 = dynamic_cast<const Tensor_d1s2_sparse *>(&b);
    if(tens1)
    {
      Tensor_d1s2_full &t = static_cast<Tensor_d1s2_full &>(*result); 
      int size = tens1->getSize();
      for (int n = 0; n < size; n++)
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            t[n][3*i+j] = v[3*i+(n%3)] * (*tens1)[(n-(n%3))/3][3*(n%3)+j];
    }
    else
    {
      const Tensor_d1s2_full &tens2 = static_cast<const Tensor_d1s2_full &>(b);
      Tensor_d1s2_full &t = static_cast<Tensor_d1s2_full &>(*result);
      int size = tens2.getSize();
      for (int m = 0; m < size; m++)
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++) {
            t[m][3*i+j] = 0.0;
            for (int k = 0; k < 3; k++)
              t[m][3*i+j] += v[3*i+k] * tens2[m][3*k+j];
          }
    }
  }
#endif
}

void
Tensor_d0s2::getDeterminant(double &det)
{ 
  det = v[0]*v[4]*v[8]+v[3]*v[7]*v[2]+v[1]*v[5]*v[6] - (v[6]*v[4]*v[2]+v[3]*v[1]*v[8]+v[0]*v[7]*v[5]);
}

void
Tensor_d0s2::getInverse(Tensor_d0s2 &t)
{
  double d;
  getDeterminant(d);
  if (d == 0.0) {
    fprintf(stderr," Caution : This Matrix is singular  \n");
    exit(-1);
  }

  double detinv = 1./d;
  t[0] = detinv*(v[4]*v[8]-v[5]*v[7]);
  t[1] = detinv*(-v[1]*v[8]+v[2]*v[7]);
  t[2] = detinv*(v[1]*v[5]-v[2]*v[4]);
  
  t[3] = detinv*(v[6]*v[5]-v[3]*v[8]);
  t[4] = detinv*(-v[6]*v[2]+v[0]*v[8]);
  t[5] = detinv*(v[3]*v[2]-v[0]*v[5]);
  
  t[6] = detinv*(-v[6]*v[4]+v[3]*v[7]);
  t[7] = detinv*(v[6]*v[1]-v[0]*v[7]);
  t[8] = detinv*(-v[3]*v[1]+v[0]*v[4]);
}

double
Tensor_d0s2::getInverseAndDeterminant(Tensor_d0s2 &t)
{
  double d;
  getDeterminant(d);
  if (d == 0.0) {
    fprintf(stderr," Caution : This Matrix is singular  \n");
    exit(-1);
  }

  double detinv = 1./d;
  t[0] = detinv*(v[4]*v[8]-v[5]*v[7]);
  t[1] = detinv*(-v[1]*v[8]+v[2]*v[7]);
  t[2] = detinv*(v[1]*v[5]-v[2]*v[4]);

  t[3] = detinv*(v[6]*v[5]-v[3]*v[8]);
  t[4] = detinv*(-v[6]*v[2]+v[0]*v[8]);
  t[5] = detinv*(v[3]*v[2]-v[0]*v[5]);

  t[6] = detinv*(-v[6]*v[4]+v[3]*v[7]);
  t[7] = detinv*(v[6]*v[1]-v[0]*v[7]);
  t[8] = detinv*(-v[3]*v[1]+v[0]*v[4]);

  return d;
}

void 
Tensor_d0s2::getTranspose(Tensor_d0s2 &t) const
{
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  t[0] = v[0]; t[1] = v[3]; t[2] = v[6];
  t[3] = v[1]; t[4] = v[4]; t[5] = v[7];
  t[6] = v[2]; t[7] = v[5]; t[8] = v[8];
#else
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      t[3*i+j] = v[3*j+i];
#endif
}

void 
Tensor_d0s2::convertToSym(Tensor_d0s2_Ss12 &t)
{
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  t[0] = v[0]; t[1] = v[1]; t[2] = v[2];
               t[3] = v[4]; t[4] = v[5];
                            t[5] = v[8];
#else
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
      t[i*(5-i)/2+j] = v[3*i+j];
#endif
} 

void
Tensor_d0s2::dblContractInto(const Tensor &b, Tensor *result) const
{
  const Tensor_d1s2_full &tens = static_cast<const Tensor_d1s2_full &>(b);
  Tensor_d1s0 &t = static_cast<Tensor_d1s0 &>(*result);
  int size = tens.getSize();
  for (int m = 0; m < size; m++) {
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
    t[m] = v[0]*tens[m][0] + v[1]*tens[m][1] + v[2]*tens[m][2] +
           v[3]*tens[m][3] + v[4]*tens[m][4] + v[5]*tens[m][5] +
           v[6]*tens[m][6] + v[7]*tens[m][7] + v[8]*tens[m][8];
#else
    t[m] = 0.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        t[m] += v[3*i+j] * tens[m][3*i+j];
#endif
  }
}

Tensor_d1s0::Tensor_d1s0(int _size)
{
  size = _size;
  v = new double[size];
  for (int i = 0; i < size; i++) 
    v[i] = 0;
}

Tensor_d1s0
Tensor_d1s0::operator + (const Tensor_d1s0 &tens) const
{
  Tensor_d1s0 t(size);
  for (int k = 0; k < size; k++)
    t[k] = v[k] + tens[k];
  return t;
}

Tensor_d1s0 
operator * (double d, const Tensor_d1s0 &tens)
{
  Tensor_d1s0 t(tens.size);
  for (int k = 0; k < t.size; k++)
    t[k] = d * tens[k];
  return t;
}

Tensor_d1s0 &
Tensor_d1s0::operator = (const Tensor_d1s0 &t)
{
  if (size != t.size) {
    if(v) delete[] v;
    size = t.size;
    v = new double[size];
  }
  for (int i = 0; i < size; ++i)
    v[i] = t.v[i];
  return *this;
}

Tensor_d1s0::Tensor_d1s0(const Tensor_d1s0 &t)
{
  size = t.size;
  v = new double[size];
  for (int i = 0; i < size; ++i)
    v[i] = t.v[i];
}

Tensor_d1s1::Tensor_d1s1(int _size)
{
  size = _size;
  v = new Tensor_d0s1[size];
}

Tensor_d1s1 &
Tensor_d1s1::operator = (const Tensor_d1s1 &t)
{
  if (size != t.size) {
    if(v) delete[] v;
    size = t.size;
    v = new Tensor_d0s1[size];
  }
  for(int i = 0; i < size; ++i)
    v[i] = t.v[i];
  return *this;
}

Tensor_d1s1::Tensor_d1s1(const Tensor_d1s1 &t)
{
  size = t.size;
  v = new Tensor_d0s1[size];
  for(int i = 0; i < size; ++i)
    v[i] = t.v[i];
}

Tensor_d1s2_full::Tensor_d1s2_full(int _size)
{
  size = _size;
  v = new Tensor_d0s2[size];
}

Tensor_d1s2_full::Tensor_d1s2_full(const Tensor_d1s2_full &t)
{
  size = t.size;
  v = new Tensor_d0s2[size];
  for(int i = 0; i < size; ++i)
    v[i] = t.v[i];
}

Tensor_d1s2_full
Tensor_d1s2_full::operator + (const Tensor_d1s2_full &tens) const
{
  Tensor_d1s2_full t(size);
  for (int k = 0; k < size; k++)
    for (int i = 0; i < 9; i++)
      t[k][i] = v[k][i] + tens[k][i];
  return t;
}

Tensor_d1s2_full
Tensor_d1s2_full::operator * (double scal)
{
  Tensor_d1s2_full t(size);
  for (int k = 0; k < size; k++)
    for (int i = 0; i < 9; i++)
      t[k][i] = v[k][i] * scal;
  return t;
}

Tensor_d1s2_full
operator * (double scal, const Tensor_d1s2_full &tens)
{
  int size = tens.size;
  Tensor_d1s2_full t(size);
  for (int k = 0; k < size; k++)
    for (int i = 0; i < 9; i++)
      t[k][i] = tens[k][i] * scal;
  return t;
};

/*
Tensor_d2s2
Tensor_d1s2_full::operator | (const Tensor_d1s2_full &tens) const
{
  Tensor_d2s2 t(size);
  for (int m = 0; m < size; m++)
    for (int n = 0; n < size; n++)
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          for (int k = 0; k < 3; k++)
            t[size*m+n][3*i+j] += v[m][3*i+k] * tens[n][3*k+j];
  return t;
}

Tensor_d1s2_full
Tensor_d1s2_full::operator | (const Tensor_d0s2 &tens) const
{
  Tensor_d1s2_full t(size);
  for (int n = 0; n < size; n++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
          t[n][3*i+j] += v[n][3*i+k] * tens[3*k+j];
  return t;
}
*/

void
Tensor_d1s2_full::dblContractWith(const Tensor_d0s4 &tens, Tensor *result) const
{
  Tensor_d1s2_full &t = static_cast<Tensor_d1s2_full &>(*result);
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  for (int m = 0; m < size; m++) {
    t[m][0] = tens[0 ]*v[m][0] + tens[1 ]*v[m][1] + tens[2 ]*v[m][2] +
              tens[3 ]*v[m][3] + tens[4 ]*v[m][4] + tens[5 ]*v[m][5] +
              tens[6 ]*v[m][6] + tens[7 ]*v[m][7] + tens[8 ]*v[m][8];
    t[m][1] = tens[9 ]*v[m][0] + tens[10]*v[m][1] + tens[11]*v[m][2] +
              tens[12]*v[m][3] + tens[13]*v[m][4] + tens[14]*v[m][5] +
              tens[15]*v[m][6] + tens[16]*v[m][7] + tens[17]*v[m][8];
    t[m][2] = tens[18]*v[m][0] + tens[19]*v[m][1] + tens[20]*v[m][2] +
              tens[21]*v[m][3] + tens[22]*v[m][4] + tens[23]*v[m][5] +
              tens[24]*v[m][6] + tens[25]*v[m][7] + tens[26]*v[m][8];
    t[m][3] = tens[27]*v[m][0] + tens[28]*v[m][1] + tens[29]*v[m][2] +
              tens[30]*v[m][3] + tens[31]*v[m][4] + tens[32]*v[m][5] +
              tens[33]*v[m][6] + tens[34]*v[m][7] + tens[35]*v[m][8];
    t[m][4] = tens[36]*v[m][0] + tens[37]*v[m][1] + tens[38]*v[m][2] +
              tens[39]*v[m][3] + tens[40]*v[m][4] + tens[41]*v[m][5] +
              tens[42]*v[m][6] + tens[43]*v[m][7] + tens[44]*v[m][8];
    t[m][5] = tens[45]*v[m][0] + tens[46]*v[m][1] + tens[47]*v[m][2] +
              tens[48]*v[m][3] + tens[49]*v[m][4] + tens[50]*v[m][5] +
              tens[51]*v[m][6] + tens[52]*v[m][7] + tens[53]*v[m][8];
    t[m][6] = tens[54]*v[m][0] + tens[55]*v[m][1] + tens[56]*v[m][2] +
              tens[57]*v[m][3] + tens[58]*v[m][4] + tens[59]*v[m][5] +
              tens[60]*v[m][6] + tens[61]*v[m][7] + tens[62]*v[m][8];
    t[m][7] = tens[63]*v[m][0] + tens[64]*v[m][1] + tens[65]*v[m][2] +
              tens[66]*v[m][3] + tens[67]*v[m][4] + tens[68]*v[m][5] +
              tens[69]*v[m][6] + tens[70]*v[m][7] + tens[71]*v[m][8];
    t[m][8] = tens[72]*v[m][0] + tens[73]*v[m][1] + tens[74]*v[m][2] +
              tens[75]*v[m][3] + tens[76]*v[m][4] + tens[77]*v[m][5] +
              tens[78]*v[m][6] + tens[79]*v[m][7] + tens[80]*v[m][8];
  }
#else
  for (int m = 0; m < size; m++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        t[m][3*i+j] = 0.0;
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
            t[m][3*i+j] += tens[i*27+j*9+k*3+l] * v[m][3*k+l];
      }
#endif
}

void
Tensor_d1s2_full::splContractWith(const Tensor_d0s2 &tens, Tensor *result) const
{
  Tensor_d1s2_full &t = static_cast<Tensor_d1s2_full &>(*result);
  for (int m = 0; m < size; m++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        t[m][3*i+j] = 0.0;
        for (int k = 0; k < 3; k++)
          t[m][3*i+j] += tens[3*i+k] * v[m][3*k+j];
      }
}

void
Tensor_d1s2_full::splContractInto(const Tensor &b, Tensor *result) const
{
  const Tensor_d1s2_full *tens = dynamic_cast<const Tensor_d1s2_full *>(&b);
  if(tens)
  {
    Tensor_d2s2 &t = static_cast<Tensor_d2s2 &>(*result);   
    for (int m = 0; m < size; m++)
      for (int n = 0; n < size; n++)
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++) {
            t[size*m+n][3*i+j] = 0.0;
	    for (int k = 0; k < 3; k++)
              t[size*m+n][3*i+j] += v[m][3*i+k]*(*tens)[n][3*k+j];
          }
  }
  else
  {
    const Tensor_d0s2 &tens1 = static_cast<const Tensor_d0s2 &>(b);
    Tensor_d1s2_full &t = static_cast<Tensor_d1s2_full &>(*result);
    for (int n = 0; n < size; n++)
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          t[n][3*i+j] = 0.0;
          for (int k = 0; k < 3; k++)
            t[n][3*i+j] += v[n][3*i+k] * tens1[3*k+j];
        }
  }
}

void
Tensor_d1s2_full::dblContractInto(const Tensor &b, Tensor *result) const
{
  const Tensor_d1s2_full &tens = static_cast<const Tensor_d1s2_full &>(b);
  Tensor_d2s0 &t = static_cast<Tensor_d2s0 &>(*result);
  for (int m = 0; m < size; m++)
    for (int n = 0; n < size; n++) {
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
      t[n+size*m] = v[m][0]*tens[n][0] + v[m][1]*tens[n][1] + v[m][2]*tens[n][2] +
                    v[m][3]*tens[n][3] + v[m][4]*tens[n][4] + v[m][5]*tens[n][5] +
                    v[m][6]*tens[n][6] + v[m][7]*tens[n][7] + v[m][8]*tens[n][8];
#else
      t[n+size*m] = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
          t[n+size*m] += v[m][3*i+j] * tens[n][3*i+j];
      }
#endif
    }                               
}

Tensor_d0s2 
Tensor_d1s2_full::operator % (const Tensor_d1s0 &tens) const
{
  Tensor_d0s2 t;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int m = 0; m < size; m++)
        t[3*i+j] += v[m][3*i+j] * tens[m];
  return t;
}

Tensor_d1s2_full &
Tensor_d1s2_full::operator = (const Tensor_d1s2_full &t)
{
  if(size != t.size) {
    if(v) delete[] v;
    size = t.size;
    v = new Tensor_d0s2[size];
  }
  for(int i = 0; i < size; ++i)
    v[i] = t.v[i];
  return *this;
}

Tensor_d1s2_full &
Tensor_d1s2_full::operator = (const Tensor_d1s2_Ss23 &t)
{
  if(size != t.getSize()) {
    if(v) delete[] v;
    size = t.getSize();
    v = new Tensor_d0s2[size];
  }
  for(int i = 0; i < size; ++i)
    v[i] = t[i];
  return *this;
}

Tensor_d1s2_full &
Tensor_d1s2_full::operator = (const Tensor_d1s2_sparse &t)
{
  if(size != t.getSize()) {
    if(v) delete[] v;
    size = t.getSize();
    v = new Tensor_d0s2[size];
  }
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  for(int i = 0, j = 0, k = 1, l = 2; i < size/3; ++i, j+=3, k+=3, l+=3) {
    v[j][0] = t[i][0]; v[j][1] = t[i][1]; v[j][2] = t[i][2];
    v[j][3] = 0;       v[j][4] = 0;       v[j][5] = 0;
    v[j][6] = 0;       v[j][7] = 0;       v[j][8] = 0;

    v[k][0] = 0;       v[k][1] = 0;       v[k][2] = 0;
    v[k][3] = t[i][3]; v[k][4] = t[i][4]; v[k][5] = t[i][5];
    v[k][6] = 0;       v[k][7] = 0;       v[k][8] = 0;

    v[l][0] = 0;       v[l][1] = 0;       v[l][2] = 0;
    v[l][3] = 0;       v[l][4] = 0;       v[l][5] = 0;
    v[l][6] = t[i][6]; v[l][7] = t[i][7]; v[l][8] = t[i][8];
  }
#else
  for(int i = 0; i < size; ++i)
    for(int j = 0; j < 3; ++j)
      for(int k = 0; k < 3; ++k)
        (*this)(i,j,k) = t(i,j,k);
#endif

  return *this;
}

void
Tensor_d1s2_full::getSpaceTranspose(Tensor_d1s2_full &t) const
{
  for (int k = 0; k < size; k++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        t[k][3*i+j] = v[k][3*j+i];
}

void 
Tensor_d1s2_full::convertToSym(Tensor_d1s2_Ss23 &t)
{
  for (int k = 0; k < size; k++)   
    for (int i = 0; i < 3; i++)
      for (int j = i; j < 3; j++)
        t[k][i*(5-i)/2+j] = v[k][3*i+j];
} 

Tensor_d1s2_Ss23
Tensor_d1s2_full::symPart() const
{
  Tensor_d1s2_Ss23 t(size);
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  for (int m = 0; m < size; m++) {
    t[m][0] = v[m][0];
    t[m][1] = 0.5*(v[m][1] + v[m][3]);
    t[m][2] = 0.5*(v[m][2] + v[m][6]);
    t[m][3] = v[m][4];
    t[m][4] = 0.5*(v[m][5] + v[m][7]);
    t[m][5] = v[m][8];
  }
#else
  for (int m = 0; m < size; m++)
    for (int i = 0; i < 3; i++)
      for (int j = i; j < 3; j++)
         t[m][i*(5-i)/2+j] = (1./2)*(v[m][3*i+j] + v[m][3*j+i]);
#endif
  return t;
};

Tensor_d1s2_sparse::Tensor_d1s2_sparse(int _size)
{
  size = _size;
  v = new Tensor_d0s2[size/3];
}

Tensor_d1s2_sparse::Tensor_d1s2_sparse(const Tensor_d1s2_sparse &t)
{
  size = t.size;
  v = new Tensor_d0s2[size/3];
  for(int i = 0; i < (size/3); ++i)
    v[i] = t.v[i];
}

Tensor_d1s2_sparse
Tensor_d1s2_sparse::operator + (const Tensor_d1s2_sparse &tens) const
{
  Tensor_d1s2_sparse t(size);
  for (int k = 0; k < size/3; k++)
    for (int i = 0; i < 9; i++)
      t[k][i] = v[k][i] + tens[k][i];
  return t;
}

Tensor_d1s2_sparse &
Tensor_d1s2_sparse::operator = (const Tensor_d1s2_sparse &t)
{
  if(size != t.size) {
    if(v) delete[] v;
    size = t.size;
    v = new Tensor_d0s2[size/3];
  }
  for(int i = 0; i < size/3; ++i)
    v[i] = t.v[i];
  return *this;
}

Tensor_d0s2 
Tensor_d1s2_sparse::operator % (const Tensor_d1s0 &b) const
{
  Tensor_d0s2 t;
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  switch(size) {
    case 12: { // 4-node tetra
      t[0] = v[0][0]*b[0]+v[1][0]*b[3]+v[2][0]*b[6]+v[3][0]*b[9 ];
      t[1] = v[0][1]*b[0]+v[1][1]*b[3]+v[2][1]*b[6]+v[3][1]*b[9 ];
      t[2] = v[0][2]*b[0]+v[1][2]*b[3]+v[2][2]*b[6]+v[3][2]*b[9 ];
      t[3] = v[0][3]*b[1]+v[1][3]*b[4]+v[2][3]*b[7]+v[3][3]*b[10];
      t[4] = v[0][4]*b[1]+v[1][4]*b[4]+v[2][4]*b[7]+v[3][4]*b[10];
      t[5] = v[0][5]*b[1]+v[1][5]*b[4]+v[2][5]*b[7]+v[3][5]*b[10];
      t[6] = v[0][6]*b[2]+v[1][6]*b[5]+v[2][6]*b[8]+v[3][6]*b[11];
      t[7] = v[0][7]*b[2]+v[1][7]*b[5]+v[2][7]*b[8]+v[3][7]*b[11];
      t[8] = v[0][8]*b[2]+v[1][8]*b[5]+v[2][8]*b[8]+v[3][8]*b[11];
    } break;
    case 24: { // 8-node hexa
      t[0] = v[0][0]*b[0]+v[1][0]*b[3]+v[2][0]*b[6]+v[3][0]*b[9 ]+v[4][0]*b[12]+v[5][0]*b[15]+v[6][0]*b[18]+v[7][0]*b[21];
      t[1] = v[0][1]*b[0]+v[1][1]*b[3]+v[2][1]*b[6]+v[3][1]*b[9 ]+v[4][1]*b[12]+v[5][1]*b[15]+v[6][1]*b[18]+v[7][1]*b[21];
      t[2] = v[0][2]*b[0]+v[1][2]*b[3]+v[2][2]*b[6]+v[3][2]*b[9 ]+v[4][2]*b[12]+v[5][2]*b[15]+v[6][2]*b[18]+v[7][2]*b[21];
      t[3] = v[0][3]*b[1]+v[1][3]*b[4]+v[2][3]*b[7]+v[3][3]*b[10]+v[4][3]*b[13]+v[5][3]*b[16]+v[6][3]*b[19]+v[7][3]*b[22];
      t[4] = v[0][4]*b[1]+v[1][4]*b[4]+v[2][4]*b[7]+v[3][4]*b[10]+v[4][4]*b[13]+v[5][4]*b[16]+v[6][4]*b[19]+v[7][4]*b[22];
      t[5] = v[0][5]*b[1]+v[1][5]*b[4]+v[2][5]*b[7]+v[3][5]*b[10]+v[4][5]*b[13]+v[5][5]*b[16]+v[6][5]*b[19]+v[7][5]*b[22];
      t[6] = v[0][6]*b[2]+v[1][6]*b[5]+v[2][6]*b[8]+v[3][6]*b[11]+v[4][6]*b[14]+v[5][6]*b[17]+v[6][6]*b[20]+v[7][6]*b[23];
      t[7] = v[0][7]*b[2]+v[1][7]*b[5]+v[2][7]*b[8]+v[3][7]*b[11]+v[4][7]*b[14]+v[5][7]*b[17]+v[6][7]*b[20]+v[7][7]*b[23];
      t[8] = v[0][8]*b[2]+v[1][8]*b[5]+v[2][8]*b[8]+v[3][8]*b[11]+v[4][8]*b[14]+v[5][8]*b[17]+v[6][8]*b[20]+v[7][8]*b[23];
    } break;
    case 30: { // 10-node tetra
      t[0] = v[0][0]*b[0]+v[1][0]*b[3]+v[2][0]*b[6]+v[3][0]*b[9 ]+v[4][0]*b[12]+v[5][0]*b[15]+v[6][0]*b[18]+v[7][0]*b[21]+v[8][0]*b[24]+v[9][0]*b[27];
      t[1] = v[0][1]*b[0]+v[1][1]*b[3]+v[2][1]*b[6]+v[3][1]*b[9 ]+v[4][1]*b[12]+v[5][1]*b[15]+v[6][1]*b[18]+v[7][1]*b[21]+v[8][1]*b[24]+v[9][1]*b[27];
      t[2] = v[0][2]*b[0]+v[1][2]*b[3]+v[2][2]*b[6]+v[3][2]*b[9 ]+v[4][2]*b[12]+v[5][2]*b[15]+v[6][2]*b[18]+v[7][2]*b[21]+v[8][2]*b[24]+v[9][2]*b[27];
      t[3] = v[0][3]*b[1]+v[1][3]*b[4]+v[2][3]*b[7]+v[3][3]*b[10]+v[4][3]*b[13]+v[5][3]*b[16]+v[6][3]*b[19]+v[7][3]*b[22]+v[8][3]*b[25]+v[9][3]*b[28];
      t[4] = v[0][4]*b[1]+v[1][4]*b[4]+v[2][4]*b[7]+v[3][4]*b[10]+v[4][4]*b[13]+v[5][4]*b[16]+v[6][4]*b[19]+v[7][4]*b[22]+v[8][4]*b[25]+v[9][4]*b[28];
      t[5] = v[0][5]*b[1]+v[1][5]*b[4]+v[2][5]*b[7]+v[3][5]*b[10]+v[4][5]*b[13]+v[5][5]*b[16]+v[6][5]*b[19]+v[7][5]*b[22]+v[8][5]*b[25]+v[9][5]*b[28];
      t[6] = v[0][6]*b[2]+v[1][6]*b[5]+v[2][6]*b[8]+v[3][6]*b[11]+v[4][6]*b[14]+v[5][6]*b[17]+v[6][6]*b[20]+v[7][6]*b[23]+v[8][6]*b[26]+v[9][6]*b[29];
      t[7] = v[0][7]*b[2]+v[1][7]*b[5]+v[2][7]*b[8]+v[3][7]*b[11]+v[4][7]*b[14]+v[5][7]*b[17]+v[6][7]*b[20]+v[7][7]*b[23]+v[8][7]*b[26]+v[9][7]*b[29];
      t[8] = v[0][8]*b[2]+v[1][8]*b[5]+v[2][8]*b[8]+v[3][8]*b[11]+v[4][8]*b[14]+v[5][8]*b[17]+v[6][8]*b[20]+v[7][8]*b[23]+v[8][8]*b[26]+v[9][8]*b[29];
    } break;
    default: {
      for (int m = 0, i = 0; m < size; m += 3, i++) {
        t[0] += v[i][0]*b[m];
        t[1] += v[i][1]*b[m];
        t[2] += v[i][2]*b[m];
        t[3] += v[i][3]*b[m+1];
        t[4] += v[i][4]*b[m+1];
        t[5] += v[i][5]*b[m+1];
        t[6] += v[i][6]*b[m+2];
        t[7] += v[i][7]*b[m+2];
        t[8] += v[i][8]*b[m+2];
      }
    } break;
  }
#else
  for (int j = 0; j < 3; j++)
    for (int m = 0; m < size; m++) {
      t[3*(m%3)+j] += v[(m-(m%3))/3][3*(m%3)+j] * b[m];
    }
#endif
  return t;
}

/*
Tensor_d1s2_sparse
Tensor_d1s2_sparse::operator | (const Tensor_d0s2 &tens) const
{
  Tensor_d1s2_sparse t(size);
  for (int n = 0; n < size; n++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        t[(n-(n%3))/3][3*(n%3)+i] += v[(n-(n%3))/3][3*(n%3)+j] * tens[3*j+i];        
  return t;
}
*/

void
Tensor_d1s2_sparse::splContractWith(const Tensor_d0s2 &tens, Tensor *result) const
{
  Tensor_d1s2_full &t = static_cast<Tensor_d1s2_full &>(*result);
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  for (int n = 0, k = 0; n < size; n += 3, k++) {
    t[n  ][0] = tens[0]*v[k][0]; t[n  ][1] = tens[0]*v[k][1]; t[n  ][2] = tens[0]*v[k][2];
    t[n  ][3] = tens[3]*v[k][0]; t[n  ][4] = tens[3]*v[k][1]; t[n  ][5] = tens[3]*v[k][2];
    t[n  ][6] = tens[6]*v[k][0]; t[n  ][7] = tens[6]*v[k][1]; t[n  ][8] = tens[6]*v[k][2];
    t[n+1][0] = tens[1]*v[k][3]; t[n+1][1] = tens[1]*v[k][4]; t[n+1][2] = tens[1]*v[k][5];
    t[n+1][3] = tens[4]*v[k][3]; t[n+1][4] = tens[4]*v[k][4]; t[n+1][5] = tens[4]*v[k][5];
    t[n+1][6] = tens[7]*v[k][3]; t[n+1][7] = tens[7]*v[k][4]; t[n+1][8] = tens[7]*v[k][5];
    t[n+2][0] = tens[2]*v[k][6]; t[n+2][1] = tens[2]*v[k][7]; t[n+2][2] = tens[2]*v[k][8];
    t[n+2][3] = tens[5]*v[k][6]; t[n+2][4] = tens[5]*v[k][7]; t[n+2][5] = tens[5]*v[k][8];
    t[n+2][6] = tens[8]*v[k][6]; t[n+2][7] = tens[8]*v[k][7]; t[n+2][8] = tens[8]*v[k][8];
  }
#else
  for (int n = 0; n < size; n++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        t[n][3*i+j] = tens[3*i+(n%3)] * v[(n-(n%3))/3][3*(n%3)+j];
#endif
}

void
Tensor_d1s2_sparse::splContractInto(const Tensor &b, Tensor *result) const
{
  const Tensor_d0s2 &tens = static_cast<const Tensor_d0s2 &>(b);
  Tensor_d1s2_sparse &t = static_cast<Tensor_d1s2_sparse &>(*result);
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  for (int n = 0, i = 0; n < size; n += 3, i++) {
    t[i][0] = v[i][0]*tens[0] + v[i][1]*tens[3] + v[i][2]*tens[6];
    t[i][1] = v[i][0]*tens[1] + v[i][1]*tens[4] + v[i][2]*tens[7];
    t[i][2] = v[i][0]*tens[2] + v[i][1]*tens[5] + v[i][2]*tens[8];
    t[i][3] = v[i][3]*tens[0] + v[i][4]*tens[3] + v[i][5]*tens[6];
    t[i][4] = v[i][3]*tens[1] + v[i][4]*tens[4] + v[i][5]*tens[7];
    t[i][5] = v[i][3]*tens[2] + v[i][4]*tens[5] + v[i][5]*tens[8];
    t[i][6] = v[i][6]*tens[0] + v[i][7]*tens[3] + v[i][8]*tens[6];
    t[i][7] = v[i][6]*tens[1] + v[i][7]*tens[4] + v[i][8]*tens[7];
    t[i][8] = v[i][6]*tens[2] + v[i][7]*tens[5] + v[i][8]*tens[8];
  }
#else
  for (int n = 0; n < size; n++)
    for (int i = 0; i < 3; i++) {
      t[(n-(n%3))/3][3*(n%3)+i] = 0.0;
      for (int j = 0; j < 3; j++)
        t[(n-(n%3))/3][3*(n%3)+i] += v[(n-(n%3))/3][3*(n%3)+j] * tens[3*j+i];
    }
#endif
}

Tensor_d1s2_Ss23 
Tensor_d1s2_sparse::symPart() const
{
  Tensor_d1s2_Ss23 t(size);
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  for (int n = 0, k = 0; n < size; n += 3, k++) {
    t[n  ][0] = v[k][0];
    t[n  ][1] = 0.5*v[k][1];
    t[n  ][2] = 0.5*v[k][2];
    t[n+1][1] = 0.5*v[k][3];
    t[n+1][3] = v[k][4];
    t[n+1][4] = 0.5*v[k][5];
    t[n+2][2] = 0.5*v[k][6];
    t[n+2][4] = 0.5*v[k][7];
    t[n+2][5] = v[k][8];
  }
#else
  for (int n = 0; n < size; n++) {
    for (int i = (n%3); i < 3; i++) {
      t[n][(n%3)*(5-(n%3))/2+i] += (1./2) * v[(n-(n%3))/3][3*(n%3)+i]; // t[n](I,J) += 1/2*v[N](I,J) where I=n%3 and J=i
    }
    for (int j = 0; j < (n%3)+1; j++) {
      t[n][j*(5-j)/2+(n%3)] += (1./2) * v[(n-(n%3))/3][3*(n%3)+j];     // t[n](I,J) += 1/2*v[N](J,I) where I=j and J=n%3
    }
  }
#endif
  return t;
}

void
Tensor_d1s2_sparse::getSymSquare(Tensor_d2s2 &t) const
{
  std::cerr << "Tensor_d1s2_sparse::getSymSquare(Tensor_d2s2 &t) is not implemented\n";
}

void 
Tensor_d1s2_sparse::getSymSquare(Tensor_d2s2_Sd12s34_dense &t) const
{
  std::cerr << "Tensor_d1s2_sparse::getSymSquare(Tensor_d2s2_Sd12s34_dense &t) is not implemented\n";
}

void 
Tensor_d1s2_sparse::getSymSquare(Tensor_d2s2_Sd12s34_sparse &t) const
{
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  for(int m = 0, k = 0, i = 0; m < size; m += 3, i++) {
    for (int n = m, j = i; n < size; n += 3, k++, j++) {
      t[k][0] = v[i][0]*v[j][0];
      t[k][1] = 0.5*(v[i][1]*v[j][0] + v[i][0]*v[j][1]);
      t[k][2] = 0.5*(v[i][2]*v[j][0] + v[i][0]*v[j][2]);
      t[k][3] = v[i][1]*v[j][1];
      t[k][4] = 0.5*(v[i][2]*v[j][1] + v[i][1]*v[j][2]);
      t[k][5] = v[i][2]*v[j][2];
    }
    for (int n = m+1, j = i; n < size; n += 3, k++, j++) {
      t[k][0] = v[i][3]*v[j][3];
      t[k][1] = 0.5*(v[i][4]*v[j][3] + v[i][3]*v[j][4]);
      t[k][2] = 0.5*(v[i][5]*v[j][3] + v[i][3]*v[j][5]);
      t[k][3] = v[i][4]*v[j][4];
      t[k][4] = 0.5*(v[i][5]*v[j][4] + v[i][4]*v[j][5]);
      t[k][5] = v[i][5]*v[j][5];
    }
    for (int n = m+2, j = i; n < size; n += 3, k++, j++) {
      t[k][0] = v[i][6]*v[j][6];
      t[k][1] = 0.5*(v[i][7]*v[j][6] + v[i][6]*v[j][7]);
      t[k][2] = 0.5*(v[i][8]*v[j][6] + v[i][6]*v[j][8]);
      t[k][3] = v[i][7]*v[j][7];
      t[k][4] = 0.5*(v[i][8]*v[j][7] + v[i][7]*v[j][8]);
      t[k][5] = v[i][8]*v[j][8];
    }
  }
#else
  for (int m = 0; m < size; m++)
    for (int n = m; n < size; n++)
      if (((n-m)%3) == 0) {
        for (int i = 0; i < 3; i++)
          for (int j = i; j < 3; j++) {
            t[(size*(size+3)-(size-m+m%3)*(size-m-m%3+3)+2*(n-m))/6][i*(5-i)/2+j] 
              = 0.5*(v[(m-(m%3))/3][3*(m%3)+j] * v[(n-(n%3))/3][3*(n%3)+i] + v[(m-(m%3))/3][3*(m%3)+i] * v[(n-(n%3))/3][3*(n%3)+j]);
          }
      }	           
#endif
}

Tensor_d2s0::Tensor_d2s0(int _size, bool isNull)
{
  size = _size; 
  if(isNull) 
    v = 0;
  else {
    int len = size*size;
    v = new double[len]; 
    //for (int i = 0; i < len; i++)
    //  v[i] = 0.;
  }
}

Tensor_d2s0::Tensor_d2s0(const Tensor_d2s0 &t)
{
  size = t.size;
  int len = size*size;
  v = new double[len];
  for(int i = 0; i < len; i++)
    v[i] = t.v[i];
}

Tensor_d2s0
Tensor_d2s0::operator + (const Tensor_d2s0 &tens) const
{
  if(tens.v == 0)
    return *this;
  if(v == 0)
    return tens;
  Tensor_d2s0 t(size);
  int len = size*size;
  for (int i = 0; i < len; i++)
    t[i] = v[i]+tens[i];
  return t;
}

Tensor_d2s0 
operator * (double scal, const Tensor_d2s0 &tens)
{
  int size = tens.getSize();
  Tensor_d2s0 t(size);
  int len = size*size;
  for(int i = 0; i < len; ++i)
    t[i] = tens[i]*scal;
  return t;
};

Tensor_d2s0 &
Tensor_d2s0::operator = (const Tensor_d2s0 &t)
{
  int len = size*size;
  if(size != t.size) {
    if(v) delete[] v;
    size = t.size;
    v = new double[len];
  }
  for(int i = 0; i < len; ++i)
    v[i] = t.v[i];
  return *this;
}

/*
Tensor_d2s0_null::Tensor_d2s0_null(const Tensor_d2s0_null &t)
{
  size = t.getSize();
  v = 0;
}
*/

Tensor_d2s2::Tensor_d2s2(int _size)
{
  size=_size;
  v = new Tensor_d0s2[size*size];
}

Tensor_d2s2::Tensor_d2s2(const Tensor_d2s2 &t)
{
  size = t.size;
  v = new Tensor_d0s2[size*size];
  for(int i = 0; i < size*size; ++i)
    v[i] = t.v[i];
}

Tensor_d2s2
Tensor_d2s2::operator + (const Tensor_d2s2 &tens) const
{
  Tensor_d2s2 t(size);
  for (int k = 0; k < size*size; k++)
    for (int i = 0; i < 9; i++)
      t[k][i] = v[k][i] + tens[k][i];
  return t;
}

Tensor_d2s2
Tensor_d2s2::operator * (double scal)
{
  Tensor_d2s2 t(size);
  for (int k = 0; k < size*size; k++)
    for (int i = 0; i < 9; i++)
      t[k][i] = v[k][i] * scal;
  return t;
}

Tensor_d2s2 
operator * (double scal, const Tensor_d2s2 &tens)
{
  int size = tens.size;
  Tensor_d2s2 t(size);
  for (int k = 0; k < size*size; k++)
    for (int i = 0; i < 9; i++)
      t[k][i] = tens[k][i] * scal;
  return t;
};

Tensor_d2s2 &
Tensor_d2s2::operator = (const Tensor_d2s2 &t)
{
  if(size != t.size) {
    if(v) delete[] v;
    size = t.size;
    v = new Tensor_d0s2[size*size];
  }
  for(int i = 0; i < size*size; ++i)
    v[i] = t.v[i];
  return *this;
}

void
Tensor_d2s2::getSpaceTranspose(Tensor_d2s2 &t)
{
  for (int k = 0; k < size*size; k++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        t[k][3*i+j] = v[k][3*j+i];
}

void 
Tensor_d2s2::convertToSym(Tensor_d2s2_Sd12s34_dense &t)
{
  for (int i = 0; i < size; i++)
    for (int j = i; j < size; j++)
      for (int k = 0; k < 3; k++)
        for (int l = k; l < 3; l++)
          t[i*(2*size-i-1)/2+j][k*(5-k)/2+l] = v[i*size+j][3*k+l];
}     

#include <typeinfo>
void
Tensor_d2s2::dblContractInto(const Tensor &b, Tensor *result) const
{
  try {
    const Tensor_d0s2 &tens = dynamic_cast<const Tensor_d0s2 &>(b);
    Tensor_d2s0 &t = static_cast<Tensor_d2s0 &>(*result);
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++) {
        t[i*size+j] = 0.0;
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
            t[i*size+j] += tens[3*k+l] * v[i*size+j][3*k+l];
      }
    return;
  }
  catch(std::bad_cast e) {}

  try {
    const Tensor_d0s2_Ss12 &tens = dynamic_cast<const Tensor_d0s2_Ss12 &>(b);
    Tensor_d2s0 &t = static_cast<Tensor_d2s0 &>(*result);
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++) {
        t[i*size+j] = 0.0;
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < k; l++)
            t[i*size+j] += tens[l*(5-l)/2+k] * v[i*size+j][3*k+l]; // t(i,j) += tens(l,k)*v(i,j,k,l)
          for (int l = k; l < 3; l++)
            t[i*size+j] += tens[k*(5-k)/2+l] * v[i*size+j][3*k+l];
        }
      }
  }
  catch(std::bad_cast e) {}
}

Tensor_d2s2_Sd12s34_null::Tensor_d2s2_Sd12s34_null(const Tensor_d2s2_Sd12s34_null &t)
{
  size = t.size;
  v = new Tensor_d0s2_Ss12[size*(size+1)/2];
}

Tensor_d2s2_Sd12s34_null &
Tensor_d2s2_Sd12s34_null::operator = (const Tensor_d2s2_Sd12s34_null &t)
{
  int len = size*(size+1)/2;
  if(size != t.size) {
    if(v) delete[] v;
    size = t.size;
    v = new Tensor_d0s2_Ss12[len];
  }
  for (int i = 0; i < len; ++i)
    v[i] = t.v[i];
  return *this;
}

Tensor_d2s2_Sd12s34_dense::Tensor_d2s2_Sd12s34_dense(int _size)
{
  size = _size;
  v = new Tensor_d0s2_Ss12[size*(size+1)/2];
}

Tensor_d2s2_Sd12s34_dense::Tensor_d2s2_Sd12s34_dense(const Tensor_d2s2_Sd12s34_dense &t)
{
  int len = size*(size+1)/2;
  size = t.size;
  v = new Tensor_d0s2_Ss12[len];
  for (int i = 0; i < len; ++i)
    v[i] = t.v[i];
}

Tensor_d2s2_Sd12s34_dense &
Tensor_d2s2_Sd12s34_dense::operator = (const Tensor_d2s2_Sd12s34_dense &t)
{
  int len = size*(size+1)/2;
  if(size != t.size) {
    if(v) delete[] v;
    size = t.size;
    v = new Tensor_d0s2_Ss12[len];
  }
  for (int i = 0; i < len; ++i)
    v[i] = t.v[i];
  return *this;
}

Tensor_d2s2_Sd12s34_sparse::Tensor_d2s2_Sd12s34_sparse(int _size)
{
  size=_size;
  v = new Tensor_d0s2_Ss12[size*(size+3)/6];
}

Tensor_d2s2_Sd12s34_sparse::Tensor_d2s2_Sd12s34_sparse(const Tensor_d2s2_Sd12s34_sparse &t)
{
  int len = size*(size+3)/6;
  size = t.size;
  v = new Tensor_d0s2_Ss12[len];
  for (int i = 0; i < len; ++i)
    v[i] = t.v[i];
}

Tensor_d2s2_Sd12s34_sparse &
Tensor_d2s2_Sd12s34_sparse::operator = (const Tensor_d2s2_Sd12s34_sparse &t)
{
  int len = size*(size+3)/6;
  if(size != t.size) {
    if(v) delete[] v;
    size = t.size;
    v = new Tensor_d0s2_Ss12[len];
  }
  for(int i = 0; i < len; ++i)
    v[i] = t.v[i];
  return *this;
}

/*
Tensor_d2s0
Tensor_d2s2_Sd12s34_sparse::operator || (const Tensor_d0s2_Ss12 &tens) const
{
  Tensor_d2s0 t(size);
  for (int i = 0; i < size; i++) {
    for (int k = 0; k < 3; k++) {
      t[(size+1)*i] += -tens[k*(5-k)/2+k]*v[i][k*(5-k)/2+k];
      for (int n = k; n < 3; n++)                       
        t[(size+1)*i] += 2*tens[k*(5-k)/2+n] * v[i][k*(5-k)/2+n];
    }
  }
  return t;
}
*/

void
Tensor_d2s2_Sd12s34_sparse::dblContractInto(const Tensor &b, Tensor *result) const
{
  const Tensor_d0s2_Ss12 &tens = static_cast<const Tensor_d0s2_Ss12 &>(b);
  Tensor_d2s0 &t = static_cast<Tensor_d2s0 &>(*result);

#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  int size2 = size*size;
  for (int i = 0; i < size2; ++i) t[i] = 0;
  for (int i = 0, k = 0; i < size; i++) {
    for (int j = i; j < size; j += 3, k++) {
      t[j*size+i] = t[i*size+j] = tens[0]*v[k][0] + tens[3]*v[k][3] + tens[5]*v[k][5] + 2*(tens[1]*v[k][1] + tens[2]*v[k][2] + tens[4]*v[k][4]);
    }
  }
#else
  int i,j;
  for (i = 0; i < size; i++)
    for (j = i; j < size; j++) {
      t[i*size+j] = 0;
      if(((j-i)%3) == 0) {
        for (int k = 0; k < 3; k++) {       
          t[i*size+j] += -tens[k*(5-k)/2+k] * v[(size*(size+3)-(size-i+i%3)*(size-i-i%3+3)+2*(j-i))/6][k*(5-k)/2+k];
          for (int n = k; n < 3; n++)                       
            t[i*size+j] += 2*tens[k*(5-k)/2+n] * v[(size*(size+3)-(size-i+i%3)*(size-i-i%3+3)+2*(j-i))/6][k*(5-k)/2+n];
        }
      }
      t[j*size+i] = t[i*size+j];           
    }
#endif
}

/*
Tensor_d2s0
Tensor_d2s2_Sd12s34_null::operator ||(const Tensor_d0s2_Ss12 &tens) const
{
  Tensor_d2s0 t(size,true);
  return t;
}
*/

void
Tensor_d2s2_Sd12s34_null::dblContractInto(const Tensor &b, Tensor *result) const
{
  Tensor_d2s0 &t = static_cast<Tensor_d2s0 &>(*result);
  int len = size*size;
  for(int i = 0; i < len; ++i)
    t[i] = 0.0;
}

/*
Tensor_d2s0
Tensor_d2s2_Sd12s34_dense::operator ||(const Tensor_d0s2_Ss12 &tens) const
{
  Tensor_d2s0 t(size);
  for (int i = 0; i < size; i++)
    for (int j = i; j < size; j++) {
      for (int k = 0; k < 3; k++) {
        t[j+size*i] += -tens[k*(5-k)/2+k] * v[i*(2*size-i-1)/2+j][k*(5-k)/2+k];
        for (int n = k; n < 3; n++)                       
          t[j+size*i] += 2*tens[k*(5-k)/2+n] * v[i*(2*size-i-1)/2+j][k*(5-k)/2+n];
      }
    }
    for (int l = 1; l < size; l++)
      for (int m = 0; m < l; m++)
        t[size*l+m] = t[size*m+l];
  return t;
}
*/

void
Tensor_d2s2_Sd12s34_dense::dblContractInto(const Tensor &b, Tensor *result) const
{
  const Tensor_d0s2_Ss12 &tens = static_cast<const Tensor_d0s2_Ss12 &>(b);
  Tensor_d2s0 &t = static_cast<Tensor_d2s0 &>(*result);
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  for (int i = 0, k = 0; i < size; i++)
    for (int j = i; j < size; j++, k++)
      t[j*size+i] = t[i*size+j] = tens[0]*v[k][0] + tens[3]*v[k][3] + tens[5]*v[k][5]
                             + 2*(tens[1]*v[k][1] + tens[2]*v[k][2] + tens[4]*v[k][4]);
#else
  for (int i = 0; i < size; i++)
    for (int j = i; j < size; j++) {
      t[i*size+j] = 0.0;
      for (int k = 0; k < 3; k++) {
        t[i*size+j] += -tens[k*(5-k)/2+k] * v[i*(2*size-i-1)/2+j][k*(5-k)/2+k];
        for (int n = k; n < 3; n++)                       
          t[i*size+j] += 2 * tens[k*(5-k)/2+n] * v[i*(2*size-i-1)/2+j][k*(5-k)/2+n];
      }
    }
  for (int l = 1; l < size; l++)
    for (int m = 0; m < l; m++)
      t[l*size+m] = t[m*size+l];
#endif
}

/*
Tensor_d1s0 
Tensor_d0s2_Ss12::operator || (const Tensor_d1s2_Ss23 &tens) const
{
  int size = tens.getSize();
  Tensor_d1s0 t(size);
  for (int i = 0; i < size; i++)
    t[i]i = v[0] * tens[i][0]
          + v[3] * tens[i][3]
          + v[5] * tens[i][5]
          + 2 * v[1] * tens[i][1]
          + 2 * v[2] * tens[i][2]
          + 2 * v[4] * tens[i][4];
  return t;
}
*/

void
Tensor_d0s2_Ss12::dblContractWith(const Tensor_d0s4_Ss12s34 &tens, Tensor *result) const
{
  Tensor_d0s2_Ss12 &t = static_cast<Tensor_d0s2_Ss12 &>(*result);
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  t[0] = tens[0][0]*v[0] + tens[0][3]*v[3] + tens[0][5]*v[5] + 2*(tens[0][1]*v[1] + tens[0][2]*v[2] + tens[0][4]*v[4]);
  t[1] = tens[1][0]*v[0] + tens[1][3]*v[3] + tens[1][5]*v[5] + 2*(tens[1][1]*v[1] + tens[1][2]*v[2] + tens[1][4]*v[4]);
  t[2] = tens[2][0]*v[0] + tens[2][3]*v[3] + tens[2][5]*v[5] + 2*(tens[2][1]*v[1] + tens[2][2]*v[2] + tens[2][4]*v[4]);
  t[3] = tens[3][0]*v[0] + tens[3][3]*v[3] + tens[3][5]*v[5] + 2*(tens[3][1]*v[1] + tens[3][2]*v[2] + tens[3][4]*v[4]);
  t[4] = tens[4][0]*v[0] + tens[4][3]*v[3] + tens[4][5]*v[5] + 2*(tens[4][1]*v[1] + tens[4][2]*v[2] + tens[4][4]*v[4]);
  t[5] = tens[5][0]*v[0] + tens[5][3]*v[3] + tens[5][5]*v[5] + 2*(tens[5][1]*v[1] + tens[5][2]*v[2] + tens[5][4]*v[4]);
#else
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++) {
      t[i*(5-i)/2+j] = 0.0;
      for (int k = 0; k < 3; k++) {
        t[i*(5-i)/2+j] += -tens[i*(5-i)/2+j][k*(5-k)/2+k] * v[k*(5-k)/2+k];
        for (int l = k; l < 3; l++)
          t[i*(5-i)/2+j] += 2.0 * tens[i*(5-i)/2+j][k*(5-k)/2+l] * v[k*(5-k)/2+l];
      }
    }
#endif
}

void
Tensor_d0s2_Ss12::dblContractInto(const Tensor &b, Tensor *result) const
{
  const Tensor_d1s2_Ss23 &tens = static_cast<const Tensor_d1s2_Ss23 &>(b);
  Tensor_d1s0 &t = static_cast<Tensor_d1s0 &>(*result);
  int size = tens.getSize();
  for (int i = 0; i < size; i++)
    t[i] = v[0]*tens[i][0] + v[3]*tens[i][3] + v[5]*tens[i][5]
         + 2.0*(v[1]*tens[i][1] + v[2]*tens[i][2] + v[4]*tens[i][4]);
}

void
Tensor_d0s2_Ss12::buildTensorOf(double *state)
{
  for(int i = 0; i < 6; ++i)
    v[i] = state[i];
}

Tensor_d0s2_Ss12
Tensor_d0s2_Ss12::operator + (const Tensor_d0s2_Ss12 &tens) const
{
  Tensor_d0s2_Ss12 t;
  for (int i = 0; i < 6; i++)
    t.v[i] = v[i] + tens[i];
  return t;
}

Tensor_d0s2_Ss12
Tensor_d0s2_Ss12::operator - (const Tensor_d0s2_Ss12 &tens) const
{
  Tensor_d0s2_Ss12 t;
  for (int i = 0; i < 6; i++)
    t.v[i] = v[i] - tens[i];
  return t;
}

Tensor_d0s2_Ss12 &
Tensor_d0s2_Ss12::operator = (const Tensor_d0s2_Ss12 &t)
{ 
  for(int i = 0; i < 6; ++i)
    v[i] = t.v[i];
  return *this;
}

Tensor_d0s2_Ss12
operator * (double scal, const Tensor_d0s2_Ss12 &tens)
{
  Tensor_d0s2_Ss12 t;
  for (int i = 0; i < 6; i++)
    t[i] = tens[i] * scal;
  return t;
};

void 
Tensor_d0s2_Ss12::getDeviation(Tensor_d0s2_Ss12 &t)
{
  double tracethird = (1./3)*(v[0]+v[3]+v[5]);
  t[0] = v[0]- tracethird;
  t[3] = v[3]- tracethird;
  t[5] = v[5]- tracethird;
  t[1] = v[1];
  t[2] = v[2];
  t[4] = v[4];
}

double 
Tensor_d0s2_Ss12::innerProduct()
{
  // return value: tr(A*A)
  return ((v[0]*v[0]+v[3]*v[3]+v[5]*v[5])+2*(v[1]*v[1]+v[2]*v[2]+v[4]*v[4]));
}

double
Tensor_d0s2_Ss12::secondInvariant()
{
  // return value: the second invariant defined as 0.5*((trA)^2 - tr(A*A))
  return v[0]*v[3]+v[3]*v[5]+v[0]*v[5] - v[1]*v[1] - v[2]*v[2] - v[4]*v[4];
}

double 
Tensor_d0s2_Ss12::getTrace()
{
  return v[0]+v[3]+v[5];
}

void
Tensor_d0s2_Ss12::addSymPart(const Tensor_d0s2 &t)
{
  v[0] += t[0];
  v[1] += 0.5*(t[1]+t[3]);
  v[2] += 0.5*(t[2]+t[6]);
  v[3] += t[4];
  v[4] += 0.5*(t[5]+t[7]);
  v[5] += t[8];
}

void
Tensor_d0s4::print() const
{
  for(int i = 0; i < 81; ++i) std::cerr << v[i] << " "; std::cerr << std::endl;
}

void
Tensor_d0s4::dblContractInto(const Tensor &b, Tensor *result) const
{
#ifdef NO_DYNAMIC_CAST_IN_NODESPACE_ARRAY_C
  b.dblContractWith(*this, result);
#else
  const Tensor_d0s2 *tens = dynamic_cast<const Tensor_d0s2 *>( &b);
  if(tens)
  {
    Tensor_d0s2 &t = static_cast<Tensor_d0s2 &>(*result);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        t[3*i+j] = 0.0;
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
            t[3*i+j] +=  v[i*27+j*9+k*3+l] * (*tens)[3*k+l];
      }
  }
  else
  {
    const Tensor_d1s2_full &tens = static_cast<const Tensor_d1s2_full &>(b);
    Tensor_d1s2_full &t = static_cast<Tensor_d1s2_full &>(*result);
    int size = tens.getSize();
    for (int m = 0; m < size; m++)
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          t[m][3*i+j] = 0.0;
          for (int k = 0; k < 3; k++) 
            for (int l = 0; l < 3; l++)
              t[m][3*i+j] += v[i*27+j*9+k*3+l] * tens[m][3*k+l];
        }
  }
#endif
}

/*
Tensor_d1s2_Ss23
Tensor_d0s4_Ss12s34::operator || (const Tensor_d1s2_Ss23 &tens) const
{
  int size = tens.getSize();
  Tensor_d1s2_Ss23 t(size);
  for (int m = 0; m < size; m++)
    for (int i = 0; i < 3; i++)
      for (int j = i; j < 3; j++)      
        for (int k = 0; k < 3; k++) {
          t[m][i*(5-i)/2+j] += -v[(i*(5-i)/2+j)][k*(5-k)/2+k] * tens[m][k*(5-k)/2+k];
          for (int l = k; l < 3; l++)
            t[m][i*(5-i)/2+j] += 2*v[(i*(5-i)/2+j)][k*(5-k)/2+l] * tens[m][k*(5-k)/2+l];
        }
  return t;
}

Tensor_d0s2_Ss12
Tensor_d0s4_Ss12s34::operator || (const Tensor_d0s2_Ss12 &tens) const
{
  Tensor_d0s2_Ss12 t;
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)      
      for (int k = 0; k < 3; k++) {
        for (int l = k; l < 3; l++)
          t[i*(5-i)/2+j] += v[(i*(5-i)/2+j)][k*(5-k)/2+l] * tens[k*(5-k)/2+l];
      }      
  return t;
}
*/

void
Tensor_d0s4_Ss12s34::print() const
{
  int I, J, K, L;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j) {
      if(j >= i) { J = j; I = i; } else { J = i; I = j; }
      for(int k = 0; k < 3; ++k)
        for(int l = 0; l < 3; ++l) {
          if(l >= k) { L = l; K = k; } else { L = k; K = l; }
          std::cerr << v[(I*(5-I)/2+J)][K*(5-K)/2+L] << " ";
        }
    }
  std::cerr << std::endl;
}

void
Tensor_d0s4_Ss12s34::dblContractInto(const Tensor &b, Tensor *result) const
{
#ifdef NO_DYNAMIC_CAST_IN_NODESPACE_ARRAY_C
  b.dblContractWith(*this, result);
#else
  const Tensor_d0s2_Ss12 *tens = dynamic_cast<const Tensor_d0s2_Ss12 *>( &b);
  if(tens)
  {
    Tensor_d0s2_Ss12 &t = static_cast<Tensor_d0s2_Ss12 &>(*result);   
    for (int i = 0; i < 3; i++)
      for (int j = i; j < 3; j++) {
        t[i*(5-i)/2+j] = 0.0;     
        for (int k = 0; k < 3; k++) {
          t[i*(5-i)/2+j] += -v[i*(5-i)/2+j][k*(5-k)/2+k] * (*tens)[k*(5-k)/2+k];
          for (int l = k; l < 3; l++)
            t[i*(5-i)/2+j] += 2.0 * v[i*(5-i)/2+j][k*(5-k)/2+l] * (*tens)[k*(5-k)/2+l];
        } 
      }     
  }
  else
  {
    const Tensor_d1s2_Ss23 &tens = static_cast<const Tensor_d1s2_Ss23 &>(b);
    Tensor_d1s2_Ss23 &t = static_cast<Tensor_d1s2_Ss23 &>(*result);
    int size = tens.getSize();
    for (int m = 0; m < size; m++)
      for (int i = 0; i < 3; i++)
        for (int j = i; j < 3; j++){
          t[m][i*(5-i)/2+j] = 0.0;     
          for (int k = 0; k < 3; k++) {
            t[m][i*(5-i)/2+j] += -v[(i*(5-i)/2+j)][k*(5-k)/2+k] * tens[m][k*(5-k)/2+k];
            for (int l = k; l < 3; l++)
              t[m][i*(5-i)/2+j] += 2.0*v[(i*(5-i)/2+j)][k*(5-k)/2+l] * tens[m][k*(5-k)/2+l];
          }
        }
  }
#endif
}

Tensor_d1s2_Ss23 
Tensor_d1s2_Ss23::operator + (const Tensor_d1s2_Ss23 &tens) const
{
  Tensor_d1s2_Ss23 t(size);
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  for (int k = 0; k < size; k++) {
    t.v[k][0] = v[k][0] + tens[k][0];
    t.v[k][1] = v[k][1] + tens[k][1];
    t.v[k][2] = v[k][2] + tens[k][2];
    t.v[k][3] = v[k][3] + tens[k][3];
    t.v[k][4] = v[k][4] + tens[k][4];
    t.v[k][5] = v[k][5] + tens[k][5];
  }
#else
  for (int k = 0; k < size; k++)
    for (int i = 0; i < 6; i++)
       t.v[k][i] = v[k][i] + tens[k][i];
#endif
  return t;
}

Tensor_d1s2_Ss23 &
Tensor_d1s2_Ss23::operator = (const Tensor_d1s2_Ss23 &t)
{
  if(size != t.size) {
    if(v) delete[] v;
    size = t.size;
    v = new Tensor_d0s2_Ss12[size];
  }
  for (int i = 0; i < size; ++i)
    v[i] = t.v[i];
  return *this;
}

Tensor_d1s2_Ss23::Tensor_d1s2_Ss23(int _size)
{
  size = _size;
  v = new Tensor_d0s2_Ss12[size];
}

Tensor_d1s2_Ss23::Tensor_d1s2_Ss23(const Tensor_d1s2_Ss23 &t)
{
  size = t.getSize();
  v = new Tensor_d0s2_Ss12[size];
  for (int i = 0; i < size; ++i)
    v[i] = t.v[i];
}

/*
Tensor_d2s0
Tensor_d1s2_Ss23::operator || (const Tensor_d1s2_Ss23 &tens) const
{
  Tensor_d2s0 t(size);
  for (int m = 0; m < size; m++)
    for (int n = m; n < size; n++)
      for (int i = 0; i < 3; i++) {
        t[n+size*m] += -v[m][i*(5-i)/2+i] * tens[n][i*(5-i)/2+i];
        for (int j = i; j < 3; j++)
          t[n+size*m] += 2 * v[m][i*(5-i)/2+j] * tens[n][i*(5-i)/2+j];
      }
  for (int o = 1; o < size; o++)
    for (int p = 0; p < o; p++)
      t[p+size*o] = t[o+size*p];
  return t;
}
*/

void
Tensor_d1s2_Ss23::dblContractWith(const Tensor_d0s4_Ss12s34 &tens, Tensor *result) const
{
  Tensor_d1s2_Ss23 &t = static_cast<Tensor_d1s2_Ss23 &>(*result);
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  for (int m = 0; m < size; m++) {
    t[m][0] = tens[0][0]*v[m][0] + tens[0][3]*v[m][3] + tens[0][5]*v[m][5] + 2*(tens[0][1]*v[m][1] + tens[0][2]*v[m][2] + tens[0][4]*v[m][4]);
    t[m][1] = tens[1][0]*v[m][0] + tens[1][3]*v[m][3] + tens[1][5]*v[m][5] + 2*(tens[1][1]*v[m][1] + tens[1][2]*v[m][2] + tens[1][4]*v[m][4]);
    t[m][2] = tens[2][0]*v[m][0] + tens[2][3]*v[m][3] + tens[2][5]*v[m][5] + 2*(tens[2][1]*v[m][1] + tens[2][2]*v[m][2] + tens[2][4]*v[m][4]);
    t[m][3] = tens[3][0]*v[m][0] + tens[3][3]*v[m][3] + tens[3][5]*v[m][5] + 2*(tens[3][1]*v[m][1] + tens[3][2]*v[m][2] + tens[3][4]*v[m][4]);
    t[m][4] = tens[4][0]*v[m][0] + tens[4][3]*v[m][3] + tens[4][5]*v[m][5] + 2*(tens[4][1]*v[m][1] + tens[4][2]*v[m][2] + tens[4][4]*v[m][4]);
    t[m][5] = tens[5][0]*v[m][0] + tens[5][3]*v[m][3] + tens[5][5]*v[m][5] + 2*(tens[5][1]*v[m][1] + tens[5][2]*v[m][2] + tens[5][4]*v[m][4]);
  }
#else
  for (int m = 0; m < size; m++)
    for (int i = 0; i < 3; i++)
      for (int j = i; j < 3; j++){
        t[m][i*(5-i)/2+j] = 0.0;
        for (int k = 0; k < 3; k++) {
          t[m][i*(5-i)/2+j] += -tens[(i*(5-i)/2+j)][k*(5-k)/2+k] * v[m][k*(5-k)/2+k];
          for (int l = k; l < 3; l++)
            t[m][i*(5-i)/2+j] += 2.0*tens[(i*(5-i)/2+j)][k*(5-k)/2+l] * v[m][k*(5-k)/2+l];
        }
      }
#endif
}

void
Tensor_d1s2_Ss23::dblContractInto(const Tensor &b, Tensor *result) const
{
  const Tensor_d1s2_Ss23 &tens = static_cast<const Tensor_d1s2_Ss23 &>(b);
  Tensor_d2s0 &t = static_cast<Tensor_d2s0 &>(*result);
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  for (int m = 0; m < size; m++)
    for (int n = m; n < size; n++) {
      t[m+size*n] = t[n+size*m] = v[m][0]*tens[n][0] + v[m][3]*tens[n][3] + v[m][5]*tens[n][5] +
                               2*(v[m][1]*tens[n][1] + v[m][2]*tens[n][2] + v[m][4]*tens[n][4]);
    }
#else
  for (int m = 0; m < size; m++)
    for (int n = m; n < size; n++) {
      t[n+size*m] = 0.0;
      for (int i = 0; i < 3; i++) {
        t[n+size*m] += -v[m][i*(5-i)/2+i] * tens[n][i*(5-i)/2+i];
        for (int j = i; j < 3; j++)
          t[n+size*m] += 2 * v[m][i*(5-i)/2+j] * tens[n][i*(5-i)/2+j];
      }
    } 
  for (int o = 1; o < size; o++)
    for (int p = 0; p < o; p++)
      t[p+size*o] = t[o+size*p];
#endif
}

void
Tensor_d1s2_Ss23::addSymPart(const Tensor_d1s2_sparse &t)
{
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  for (int n = 0, k = 0; n < size; n += 3, k++) {
    v[n  ][0] += t[k][0];
    v[n  ][1] += 0.5*t[k][1];
    v[n  ][2] += 0.5*t[k][2];
    v[n+1][1] += 0.5*t[k][3];
    v[n+1][3] += t[k][4];
    v[n+1][4] += 0.5*t[k][5];
    v[n+2][2] += 0.5*t[k][6];
    v[n+2][4] += 0.5*t[k][7];
    v[n+2][5] += t[k][8];
  }
#else
  for (int n = 0; n < size; n++) {
    for (int i = (n%3); i < 3; i++) {
      v[n][(n%3)*(5-(n%3))/2+i] += (1./2) * t[(n-(n%3))/3][3*(n%3)+i]; // t[n](I,J) += 1/2*v[N](I,J) where I=n%3 and J=i
    }
    for (int j = 0; j < (n%3)+1; j++) {
      v[n][j*(5-j)/2+(n%3)] += (1./2) * t[(n-(n%3))/3][3*(n%3)+j];     // t[n](I,J) += 1/2*v[N](J,I) where I=j and J=n%3
    }
  }
#endif
}

void
Tensor_d1s2_Ss23::addSymPart(const Tensor_d1s2_full &t)
{
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  for (int m = 0; m < size; m++) {
    v[m][0] += t[m][0];
    v[m][1] += 0.5*(t[m][1] + t[m][3]);
    v[m][2] += 0.5*(t[m][2] + t[m][6]);
    v[m][3] += t[m][4];
    v[m][4] += 0.5*(t[m][5] + t[m][7]);
    v[m][5] += t[m][8];
  }
#else
  for (int m = 0; m < size; m++)
    for (int i = 0; i < 3; i++)
      for (int j = i; j < 3; j++)
         v[m][i*(5-i)/2+j] += (1./2)*(t[m][3*i+j] + t[m][3*j+i]);
#endif
}

void
Tensor_d1s2_Ss23::assignSymPart(const Tensor_d1s2_sparse &s, const Tensor_d1s2_full &t)
{
#ifdef UNROLL_LOOPS_IN_NODESPACE_ARRAY_C
  for (int n = 0, k = 0; n < size; n += 3, k++) {
    v[n  ][0] = s[k][0] + t[n][0];
    v[n  ][1] = 0.5*(s[k][1] + t[n][1] + t[n][3]);
    v[n  ][2] = 0.5*(s[k][2] + t[n][2] + t[n][6]);
    v[n  ][3] = t[n][4];
    v[n  ][4] = 0.5*(t[n][5] + t[n][7]);
    v[n  ][5] = t[n][8];

    v[n+1][0] = t[n+1][0];
    v[n+1][1] = 0.5*(s[k][3] + t[n+1][1] + t[n+1][3]);
    v[n+1][2] = 0.5*(t[n+1][2] + t[n+1][6]);
    v[n+1][3] = s[k][4] + t[n+1][4];
    v[n+1][4] = 0.5*(s[k][5] + t[n+1][5] + t[n+1][7]);
    v[n+1][5] = t[n+1][8];

    v[n+2][0] = t[n+2][0];
    v[n+2][1] = 0.5*(t[n+2][1] + t[n+2][3]);
    v[n+2][2] = 0.5*(s[k][6] + t[n+2][2] + t[n+2][6]);
    v[n+2][3] = t[n+2][4];
    v[n+2][4] = 0.5*(s[k][7] + t[n+2][5] + t[n+2][7]);
    v[n+2][5] = s[k][8] + t[n+2][8];
  }
#else
  setZero();
  addSymPart(s);
  addSymPart(t);
#endif
}

void
Tensor_d0s2_Ss12_diag::dblContractWith(const Tensor_d0s4_Ss12s34_diag &tens, Tensor *result) const
{
  Tensor_d0s2_Ss12_diag &t = static_cast<Tensor_d0s2_Ss12_diag &>(*result);
  t[0] = tens[0][0]*v[0] + tens[0][1]*v[1] + tens[0][2]*v[2];
  t[1] = tens[1][0]*v[0] + tens[1][1]*v[1] + tens[1][2]*v[2];
  t[2] = tens[2][0]*v[0] + tens[2][1]*v[1] + tens[2][2]*v[2];
}

void
Tensor_d0s2_Ss12_diag::dblContractInto(const Tensor &b, Tensor *result) const
{
  const Tensor_d1s2_Ss23_diag &tens = static_cast<const Tensor_d1s2_Ss23_diag &>(b);
  Tensor_d1s0 &t = static_cast<Tensor_d1s0 &>(*result);
  int size = tens.getSize();
  for (int i = 0; i < size; i++)
    t[i] = v[0]*tens[i][0] + v[1]*tens[i][1] + v[2]*tens[i][2];
}

void
Tensor_d0s2_Ss12_diag::buildTensorOf(double *state)
{
  for(int i = 0; i < 3; ++i)
    v[i] = state[i];
}

Tensor_d0s2_Ss12_diag
Tensor_d0s2_Ss12_diag::operator + (const Tensor_d0s2_Ss12_diag &tens) const
{
  Tensor_d0s2_Ss12_diag t;
  for (int i = 0; i < 3; i++)
    t.v[i] = v[i] + tens[i];
  return t;
}

Tensor_d0s2_Ss12_diag
Tensor_d0s2_Ss12_diag::operator - (const Tensor_d0s2_Ss12_diag &tens) const
{
  Tensor_d0s2_Ss12_diag t;
  for (int i = 0; i < 3; i++)
    t.v[i] = v[i] - tens[i];
  return t;
}

Tensor_d0s2_Ss12_diag &
Tensor_d0s2_Ss12_diag::operator = (const Tensor_d0s2_Ss12_diag &t)
{ 
  for(int i = 0; i < 3; ++i)
    v[i] = t.v[i];
  return *this;
}

Tensor_d0s2_Ss12_diag
operator * (double scal, const Tensor_d0s2_Ss12_diag &tens)
{
  Tensor_d0s2_Ss12_diag t;
  for (int i = 0; i < 3; i++)
    t[i] = tens[i] * scal;
  return t;
};

void 
Tensor_d0s2_Ss12_diag::getDeviation(Tensor_d0s2_Ss12_diag &t)
{
  double tracethird = (1./3)*(v[0]+v[1]+v[2]);
  t[0] = v[0] - tracethird;
  t[1] = v[1] - tracethird;
  t[2] = v[2] - tracethird;
}

double 
Tensor_d0s2_Ss12_diag::innerProduct()
{
  // return value: tr(A*A)
  return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}

double
Tensor_d0s2_Ss12_diag::secondInvariant()
{
  // return value: the second invariant defined as 0.5*((trA)^2 - tr(A*A))
  return v[0]*v[1]+v[1]*v[2]+v[0]*v[2];
}

double 
Tensor_d0s2_Ss12_diag::getTrace()
{
  return v[0]+v[1]+v[2];
}

Tensor_d1s2_Ss23_diag 
Tensor_d1s2_Ss23_diag::operator + (const Tensor_d1s2_Ss23_diag &tens) const
{
  Tensor_d1s2_Ss23_diag t(size);
  for (int k = 0; k < size; k++) {
    t.v[k][0] = v[k][0] + tens[k][0];
    t.v[k][1] = v[k][1] + tens[k][1];
    t.v[k][2] = v[k][2] + tens[k][2];
  }
  return t;
}

Tensor_d1s2_Ss23_diag &
Tensor_d1s2_Ss23_diag::operator = (const Tensor_d1s2_Ss23_diag &t)
{
  if(size != t.size) {
    if(v) delete[] v;
    size = t.size;
    v = new Tensor_d0s2_Ss12_diag[size];
  }
  for (int i = 0; i < size; ++i)
    v[i] = t.v[i];
  return *this;
}

Tensor_d1s2_Ss23_diag::Tensor_d1s2_Ss23_diag(int _size)
{
  size = _size;
  v = new Tensor_d0s2_Ss12_diag[size];
}

Tensor_d1s2_Ss23_diag::Tensor_d1s2_Ss23_diag(const Tensor_d1s2_Ss23_diag &t)
{
  size = t.getSize();
  v = new Tensor_d0s2_Ss12_diag[size];
  for (int i = 0; i < size; ++i)
    v[i] = t.v[i];
}

void
Tensor_d1s2_Ss23_diag::dblContractWith(const Tensor_d0s4_Ss12s34_diag &tens, Tensor *result) const
{
  Tensor_d1s2_Ss23_diag &t = static_cast<Tensor_d1s2_Ss23_diag &>(*result);
  for (int m = 0; m < size; m++) {
    t[m][0] = tens[0][0]*v[m][0] + tens[0][1]*v[m][1] + tens[0][2]*v[m][2];
    t[m][1] = tens[1][0]*v[m][0] + tens[1][1]*v[m][1] + tens[1][2]*v[m][2];
    t[m][2] = tens[2][0]*v[m][0] + tens[2][1]*v[m][1] + tens[2][2]*v[m][2];
  }
}

void
Tensor_d1s2_Ss23_diag::dblContractInto(const Tensor &b, Tensor *result) const
{
  const Tensor_d1s2_Ss23_diag &tens = static_cast<const Tensor_d1s2_Ss23_diag &>(b);
  Tensor_d2s0 &t = static_cast<Tensor_d2s0 &>(*result);
  for (int m = 0; m < size; m++)
    for (int n = m; n < size; n++) {
      t[m+size*n] = t[n+size*m] = v[m][0]*tens[n][0] + v[m][1]*tens[n][1] + v[m][2]*tens[n][2];
    }
}

void
Tensor_d0s4_Ss12s34_diag::dblContractInto(const Tensor &b, Tensor *result) const
{
  b.dblContractWith(*this, result);
}

Tensor_d2s2_Sd12s34_dense_diag::Tensor_d2s2_Sd12s34_dense_diag(int _size)
{
  size = _size;
  v = new Tensor_d0s2_Ss12_diag[size*(size+1)/2];
}

Tensor_d2s2_Sd12s34_dense_diag::Tensor_d2s2_Sd12s34_dense_diag(const Tensor_d2s2_Sd12s34_dense_diag &t)
{
  int len = size*(size+1)/2;
  size = t.size;
  v = new Tensor_d0s2_Ss12_diag[len];
  for (int i = 0; i < len; ++i)
    v[i] = t.v[i];
}

Tensor_d2s2_Sd12s34_dense_diag &
Tensor_d2s2_Sd12s34_dense_diag::operator = (const Tensor_d2s2_Sd12s34_dense_diag &t)
{
  int len = size*(size+1)/2;
  if(size != t.size) {
    if(v) delete[] v;
    size = t.size;
    v = new Tensor_d0s2_Ss12_diag[len];
  }
  for (int i = 0; i < len; ++i)
    v[i] = t.v[i];
  return *this;
}

void
Tensor_d2s2_Sd12s34_dense_diag::dblContractInto(const Tensor &b, Tensor *result) const
{
  const Tensor_d0s2_Ss12_diag &tens = static_cast<const Tensor_d0s2_Ss12_diag &>(b);
  Tensor_d2s0 &t = static_cast<Tensor_d2s0 &>(*result);
  for (int i = 0, k = 0; i < size; i++)
    for (int j = i; j < size; j++, k++)
      t[j*size+i] = t[i*size+j] = tens[0]*v[k][0] + tens[1]*v[k][1] + tens[2]*v[k][2];
}

