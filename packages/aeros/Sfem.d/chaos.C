#include <cmath>
#include <cstdlib>
#include <vector>
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Sfem.d/chaos.h>

/* defining the function nterm */

int nterm(int x,int y)
{

  int fac = 1;
  int nterm = 1;

  for (int i=0; i<x; i++)
    {
      if (i>0)
        fac = fac * i;
      nterm = nterm*(y+i);
    }
  fac = fac*x;
  if (fac>0)
    nterm = nterm/fac;
  else
    nterm = 1;

  return nterm;
}


Chaos::Chaos(int nd, int od) {
  ndim = nd;
  order = od; 
}
int Chaos::tnterms(){
  int x = 0;
  for (int i=1; i<=order; i++)
    x += nterm(i, ndim);
  return x;
}

void Chaos::HermiteToPowers(int &mr, GenFullM<double> &icoef){

  const int maxorder = mr;

  for (int i=0; i<=maxorder; i++)
    {
      for (int j=0; j<=maxorder; j++)
        {
          icoef[i][j] = 0;
        }
    }
  icoef [1][1] = 1;
  icoef [2][2] = 1;
 
  for (int k=2; k<=maxorder; k++)
    {
      for (int i=1; i<=maxorder; i++)
        {icoef[k][i] = icoef[k-1][i-1];
        }
      for (int i=0; i<=maxorder; i++)
        {icoef[k][i] = icoef[k][i]-icoef[k-2][i]*(k-2);
        }
    }
}

void Chaos::ScHermiteToPowers(int &mr,double & sc, GenFullM<double> &icoef){
  
  const int maxorder = mr;

  for (int i=0; i<=maxorder; i++)
    {
      for (int j=0; j<=maxorder; j++)
        {   
          icoef[i][j] = 0;
        }
    }
  icoef [1][1] = 1;
  icoef [2][2] = 1;
    
  for (int k=2; k<=maxorder; k++)
    {
      for (int i=1; i<=maxorder; i++)
        {icoef[k][i] = icoef[k-1][i-1];
        }
      for (int i=0; i<=maxorder; i++)
        {icoef[k][i] = icoef[k][i]-icoef[k-2][i]*(k-2);
        }
    }
  for(int i=1; i<=maxorder; i++)
    {
     for(int j=1; j<=maxorder; j++) 
       {
        if(j<i)
         icoef[i][j] *= pow(sc,(i-j));
       }
    }
}

void Chaos::mult1dHermite(GenVector<double> &A, GenVector<double> &B, GenVector<double> &C){
  int noi = A.size();
  int noj = B.size();
  int nok = C.size(); 
 
  // arrays containing the coef. of each of the polynomials 
   
  for(int i=0; i<nok; i++)
    C[i] = 0;
 
  for (int i=0; i<noi; i++)
    {
      for(int j=0; j<noj; j++)
        {
          C[i+j] += A[i]*B[j];
        }
    }
}
void Chaos::OrthoOrder(GenFullM<double> &D){

  int nterms;
  int  maxterms = nterm(order, ndim);

  //int c[maxterms+1][order+1];
  //int c1[maxterms+1][ndim+1];
  std::vector< std::vector<int> > c(maxterms+1, std::vector<int>(order+1));
  std::vector< std::vector<int> > c1(maxterms+1, std::vector<int>(ndim+1));

  
  int count = 0;
  int count1 = 0;

  for (int iorder=1; iorder<=order; iorder++)
    {
      nterms = nterm(iorder, ndim);

      for (int i=1; i<=iorder; i++)
	{
	  for(int j=1; j<=nterms; j++)
	    {
	      c[j][i] = 0;
	    }
	}

      int iterm = 1;
      count++;

      for(int i=1; i<=iorder; i++)
	c[iterm][i] = 1;


      for(int i=1; i<=ndim; i++)
	{
	  c1[iterm][i] = 0;
	  for(int j=1; j<=iorder; j++)
	    {
	      if (c[iterm][j]== i)
		c1[iterm][i]+=1;
	    }
	}

  
      while (iterm < nterms)
	{
          iterm++;
	  int iswitch = 0;
	  count++;

	  for(int i=iorder; i>=1; i--)
	    {
	      if(((c[iterm-1][i]+1)<=ndim)&&(iswitch==0))
		{
		  iswitch = 1;
		  c[iterm][i] = c[iterm-1][i] + 1;
		  
		  if(i < iorder)
		    {
		      for(int j=i+1; j<=iorder; j++)
			c[iterm][j] = c[iterm][i];
		    }
		  
		  if (i>1)
		    {
		      for(int j=1; j<=i-1; j++)
			c[iterm][j] = c[iterm-1][j];
		    }
		}
	    }

  
      
	  for (int i=1;  i<=ndim; i++)
	    {
	      c1[iterm][i] = 0;
	      for(int j=1; j<=iorder; j++)
		{
		  if( c[iterm][j] == i)
		    c1[iterm][i]+=1;
		}
	    }      
	}
      
      
      for (int it = 1; it<=nterms; it++)
	{
          count1++;
	  for(int j=1; j<=ndim; j++)
	    {
               D[count1][j] = c1[it][j];	   
	    }
	}
    }
}

