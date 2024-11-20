#include <cstdio>
#include <Math.d/FullRectMatrix.h>

FullRectMatrix::FullRectMatrix()
{
 sizem = 0; sizen = 0; value = 0; myval=0;
}


FullRectMatrix::FullRectMatrix(int i, int j,  double *l)
{
  sizem = i; 
  sizen = j; 
  if (l) { 
     value=l; 
     myval=0; 
  } 
  else { 
     value=new double [sizem*sizen]; 
     myval=1; 
   } 
}

FullRectMatrix::~FullRectMatrix()
{
// RT 
  if (myval) delete [] value; 
  value=0; 
}


FullRectMatrix
FullRectMatrix::operator *= (double v)
{
 int length = sizem*sizen;
 int i;
 for(i = 0; i < length; ++i)
    value[i] *= v;
 return *this;
}


void
FullRectMatrix::zero()
{
 int length = sizem*sizen;
 int i;
 for(i=0; i<length; ++i)
   value[i] = 0.0;
}
