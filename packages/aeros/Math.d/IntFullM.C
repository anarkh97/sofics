#include <cstdio>
#include <Math.d/IntFullM.h>

IntFullM::IntFullM()
{
  nrow    = 0; 
  ncolumn = 0;
  v = new int[1];
}

IntFullM::IntFullM(int nr)
{
  nrow    = nr;
  ncolumn = nr;
  v = new int [nrow*ncolumn];
}

IntFullM::IntFullM(int nr, int nc)
{
  nrow    = nr;
  ncolumn = nc;
  v = new int [nrow*ncolumn];
}

IntFullM::~IntFullM()
{
  if(v) delete [] v;
}

void
IntFullM::zero()
{
 int i;
 for(i=0; i<nrow*ncolumn; ++i)
   v[i] = 0;
}

void
IntFullM::print()
{
  int i,j;
  for (i=0; i< nrow; i++) {
    for (j=0; j< ncolumn; j++)
      fprintf(stderr,"%4d  ",(*this)[i][j]);
    fprintf(stderr,"\n");

  }
}

