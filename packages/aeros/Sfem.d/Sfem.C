#include <Sfem.d/Sfem.h>
#include <Sfem.d/cijk.h>
#include <cmath>

void Sfem::setOrder(int _output_order) 
{
  output_order = _output_order;
}
 
void Sfem::computeLP()
{
  if(Gauss) L=ndim+1;             
  else {                           // if(GaussSqr)
    L=0;
    for (int i=0; i<=2; i++) {
      L = L + no_terms(i, ndim);
    }
  }
 
  nonzindex = new int[L];

  if(Gauss) for (int i=0; i<L; i++) nonzindex[i]=i; 
  else {                           // if(GaussSqr)
    if(ndim==1) {
      nonzindex[0] = 0;
      nonzindex[1] = -1;
      nonzindex[2] = 1;
    }
    else {
      for (int i=0; i<L; i++) nonzindex[i]=-1;
      nonzindex[0] = 0;

      int tempos = ndim+1; 
      int count = 1; 
      for (int i=0; i<ndim; i++) {
        for (int j=i; j<ndim; j++) {
          if(j==i) {
           nonzindex[tempos]=count;
           count=count+1;
          }
          tempos = tempos+1;
        }
      }
    }

  }

  P=0;
  for (int i=0; i<=output_order; i++) {
    P = P + no_terms(i, ndim);
  }

}


int Sfem::no_terms(int x,int y)
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
    nterm = 1; /// assuming P >= L, else should be 0

  return nterm;
}

void Sfem::build_psisqr()
{
  Cijk secmom(ndim, output_order, P);
  for (int i=0;i<P;i++) {
    psisqr[i] = secmom.expectation(0,i,i);
  }
}

void Sfem::computeNnzBlocks(double* bln)
{
 build_psisqr();
 double totalbln = 0;
 for (int i=0; i<P; i++)  totalbln = totalbln + bln[i]*psisqr[i];
 totalbln = sqrt(totalbln); // YYY DG \|.\|_(\Omega \times D), recall bln[i] are already squared

 nnzblkindex = new int[P];
 for (int i=0; i<P; i++) {
   if(sqrt(bln[i]*psisqr[i])/totalbln > 0.001) nnzblkindex[i]=1; // YYY DG Decide the cut-off
   else nnzblkindex[i]=0;
 }
}

void Sfem::init_XiPsi()
{
  xi=new double[ndim];
  psi=new double[P];
  psisqr=new double[P];
}

void Sfem::genXi(int seed)
{
// generate and store in memory xi for a particular seed

 zufalli_(seed);
 normalen_(ndim,xi);


/* if(seed==0) {
   readfile.open("xifile",ios::in);

   std::ifstream readintegfile2("integparamfile",ios::in);
   int junk;
   readintegfile2 >> junk;
   readintegfile2 >> nosamp_deletelater;
   readintegfile2.close();
 }
 for(int i=0;i<ndim;i++)  readfile >> xi[i];

 cerr << "xi read from file = " << endl;
 for(int i=0;i<ndim;i++)  cerr << xi[i] << " ";
 cerr << endl;

 if(seed==nosamp_deletelater-1) readfile.close();
*/

}

void Sfem::genXiPsi(int seed)
{
 genXi(seed);
 makealpha();
 build_psi();
}

struct klotz0_1_ {
    double buff[607];
    int ptr;
};
                                                                                                                                             
#define klotz0_1 (*(struct klotz0_1_ *) &klotz0_)
                                                                                                                                             
struct klotz1_1_ {
    double xbuff[1024];
    int first, xptr;
};
                                                                                                                                             
#define klotz1_1 (*(struct klotz1_1_ *) &klotz1_)
                                                                                                                                             
// Initialized data
struct {
    int fill_1[1214];
    int e_2;
    } klotz0_ = { {0}, 0 };
                                                                                                                                             
struct {
    double fill_1[1024];
    int e_2[2];
    double e_3;
    } klotz1_ = { {0}, 0, 0, 0. };


int Sfem::zufall_(int nn, double *a)
{
    int buffsz = 607;

    int left, aptr, bptr, aptr0, i, k, q;
    double t;
    int vl, qq, k273, k607, kptr;

/* portable lagged Fibonacci series uniform random number */
/* generator with "lags" -273 und -607: */
/* W.P. Petersen, IPS, ETH Zuerich, 19 Mar. 92 */

    aptr = 0;

L1:

    if (nn <= 0) {
	return 0;
    }

/* factor nn = q*607 + r */

    q = (nn - 1) / 607;
    left = buffsz - klotz0_1.ptr;

    if (q <= 1) {

/* only one or fewer full segments */

	if (nn < left) {
            kptr = klotz0_1.ptr;
	    for (i = 0; i < nn; ++i) {
		a[i + aptr] = klotz0_1.buff[kptr + i];
	    }
	    klotz0_1.ptr += nn;
	    return 0;
	} else {
            kptr = klotz0_1.ptr;
//#pragma _CRI ivdep
	    for (i = 0; i < left; ++i) {
		a[i + aptr] = klotz0_1.buff[kptr + i];
	    }
	    klotz0_1.ptr = 0;
	    aptr += left;
	    nn -= left;
/*  buff -> buff case */
	    vl = 273;
	    k273 = 334;
	    k607 = 0;
	    for (k = 0; k < 3; ++k) {
//#pragma _CRI ivdep
		for (i = 0; i < vl; ++i) {
		   t = klotz0_1.buff[k273+i]+klotz0_1.buff[k607+i];
		   klotz0_1.buff[k607+i] = t - (double) ((int) t);
		}
		k607 += vl;
		k273 += vl;
		vl = 167;
		if (k == 0) {
		    k273 = 0;
		}
	    }
	    goto L1;
	}
    } else {

/* more than 1 full segment */

        kptr = klotz0_1.ptr;
//#pragma _CRI ivdep
	for (i = 0; i < left; ++i) {
	    a[i + aptr] = klotz0_1.buff[kptr + i];
	}
	nn -= left;
	klotz0_1.ptr = 0;
	aptr += left;

/* buff -> a(aptr0) */

	vl = 273;
	k273 = 334;
	k607 = 0;
	for (k = 0; k < 3; ++k) {
	    if (k == 0) {
//#pragma _CRI ivdep
		for (i = 0; i < vl; ++i) {
		    t = klotz0_1.buff[k273+i]+klotz0_1.buff[k607+i];
		    a[aptr + i] = t - (double) ((int) t);
		}
		k273 = aptr;
		k607 += vl;
		aptr += vl;
		vl = 167;
	    } else {
//#pragma _CRI ivdep
		for (i = 0; i < vl; ++i) {
		    t = a[k273 + i] + klotz0_1.buff[k607 + i];
		    a[aptr + i] = t - (double) ((int) t);
		}
		k607 += vl;
		k273 += vl;
		aptr += vl;
	    }
	}
	nn += -607;

/* a(aptr-607) -> a(aptr) for last of the q-1 segments */

	aptr0 = aptr - 607;
	vl = 607;

	for (qq = 0; qq < q-2; ++qq) {
	    k273 = aptr0 + 334;
//#pragma _CRI ivdep
	    for (i = 0; i < vl; ++i) {
		t = a[k273 + i] + a[aptr0 + i];
		a[aptr + i] = t - (double) ((int) t);
	    }
	    nn += -607;
	    aptr += vl;
	    aptr0 += vl;
	}

/* a(aptr0) -> buff, last segment before residual */

	vl = 273;
	k273 = aptr0 + 334;
	k607 = aptr0;
	bptr = 0;
	for (k = 0; k < 3; ++k) {
	    if (k == 0) {
//#pragma _CRI ivdep
		for (i = 0; i < vl; ++i) {
		    t = a[k273 + i] + a[k607 + i];
		    klotz0_1.buff[bptr + i] = t - (double) ((int) t);
		}
		k273 = 0;
		k607 += vl;
		bptr += vl;
		vl = 167;
	    } else {
//#pragma _CRI ivdep
		for (i = 0; i < vl; ++i) {
		    t = klotz0_1.buff[k273 + i] + a[k607 + i];
		    klotz0_1.buff[bptr + i] = t - (double) ((int) t);
		}
		k607 += vl;
		k273 += vl;
		bptr += vl;
	    }
	}
	goto L1;
    }
} /* zufall_ */


int Sfem::zufalli_(int seed_local)
{
    /* Initialized data */

    int kl = 9373;
    int ij = 1802;

    /* Local variables */
    int i, j, k, l, m;
    double s, t;
    int ii, jj;


/*  generates initial seed buffer by linear congruential */
/*  method. Taken from Marsaglia, FSU report FSU-SCRI-87-50 */
/*  variable seed should be 0 < seed <31328 */

    if (seed_local != 0) {
	ij = seed_local;
    }

    i = ij / 177 % 177 + 2;
    j = ij % 177 + 2;
    k = kl / 169 % 178 + 1;
    l = kl % 169;
    for (ii = 0; ii < 607; ++ii) {
	s = 0.;
	t = .5;
	for (jj = 1; jj <= 24; ++jj) {
	    m = i * j % 179 * k % 179;
	    i = j;
	    j = k;
	    k = m;
	    l = (l * 53 + 1) % 169;
	    if (l * m % 64 >= 32) {
		s += t;
	    }
	    t *= (double).5;
	}
	klotz0_1.buff[ii] = s;
    }
    return 0;
} /* zufalli_ */


int Sfem::normalen_(int nn, double *x)
{
    int buffsz = 1024;

    /* Local variables */
    int left, i, ptr, kptr;
//    extern int normal00_();
/* Box-Muller method for Gaussian random numbers */

    if (nn <= 0) {
	return 0;
    }
    if (klotz1_1.first == 0) {
	normal00_();
	klotz1_1.first = 1;
    }
    ptr = 0;

L1:
    left = buffsz - klotz1_1.xptr;
    if (nn < left) {
        kptr = klotz1_1.xptr;
//#pragma _CRI ivdep
	for (i = 0; i < nn; ++i) {
	    x[i + ptr] = klotz1_1.xbuff[kptr + i];
	}
	klotz1_1.xptr += nn;
	return 0;
    } else {
        kptr = klotz1_1.xptr;
//#pragma _CRI ivdep
	for (i = 0; i < left; ++i) {
	    x[i + ptr] = klotz1_1.xbuff[kptr + i];
	}
	klotz1_1.xptr = 0;
	ptr += left;
	nn -= left;
	normal00_();
	goto L1;
    }
} /* normalen_ */


int Sfem::normal00_()
{
    int i;
    double twopi, r1, r2, t1, t2;
 //   extern int zufall_();

    twopi = 6.2831853071795862;
    zufall_(1024, klotz1_1.xbuff);
//#pragma _CRI ivdep
    for (i = 0; i < 1023; i += 2) {
	r1 = twopi * klotz1_1.xbuff[i];
	t1 = cos(r1);
	t2 = sin(r1);
	r2 = sqrt(-2.*(log(1. - klotz1_1.xbuff[i+1])));
	klotz1_1.xbuff[i]   = t1 * r2;
	klotz1_1.xbuff[i+1] = t2 * r2;
    }

    return 0;
} /* normal00_ */


void Sfem::build_psi()
{
 double **pceach;
 
 pceach = new double*[ndim];
 for (int i=0; i<ndim; ++i) {
   pceach[i] = new double[output_order+1];
   for (int j=0; j<output_order+1; ++j) pceach[i][j]=0;
 }
 
 for (int i=0; i<ndim; ++i) {
   pceach[i][0]=1;
   pceach[i][1]=xi[i];
 }
 if (output_order > 1) {
   for (int j=2; j<output_order+1; ++j) {
     for (int i=0; i<ndim; ++i) {
       pceach[i][j]=xi[i]*pceach[i][j-1]-(j-1)*pceach[i][j-2];
     }
   } 
 }

 double temp_prod;
 for (int j=0; j<P; ++j) {
   temp_prod=1;
   for (int i=0; i<ndim; i++) {
     temp_prod=temp_prod*pceach[i][alphamat[i][j]];
   }
    psi[j]=temp_prod;
 }
  for (int i=0; i<ndim; i++) delete [] pceach[i];
  delete [] pceach;
}

void Sfem::makealpha()
{
 if(alphamat==0) {
   alphamat=new int*[ndim];
   for (int i=0; i<ndim; ++i) alphamat[i] = new int[P];
 }

 for (int i=0; i<ndim; ++i) {
   for (int j=0; j<P; ++j) {
     if(j==i+1)  
       alphamat[i][j]=1;
     else
       alphamat[i][j]=0;
   }
 }

 int **parray;
 parray = new int*[ndim];
 for (int i=0; i<ndim; ++i) parray[i] = new int[output_order];
 for (int i=0; i<ndim; ++i) {
   parray[i][0]=1;
   for (int j=1; j<output_order; ++j) parray[i][j]=0;
 }

 int L_cap;
 int sum;
 int P_cap=ndim;
 for (int k=1; k<output_order; ++k) {
   L_cap=P_cap;
   for (int i=0; i<ndim; ++i) {
     sum=0;
     for (int m=i; m<ndim; ++m) {
       sum=sum+parray[m][k-1];
     }
     parray[i][k]=sum;
   }

  for (int j=0; j<ndim; ++j) {
     for (int m=L_cap-parray[j][k]; m<L_cap; ++m) {
       P_cap=P_cap+1;
       for (int i=0; i<ndim; ++i) {
         alphamat[i][P_cap]=alphamat[i][m+1];
       }
       alphamat[j][P_cap]=alphamat[j][P_cap]+1;
     }
   }
 }

 for (int i=0; i<ndim; i++) delete [] parray[i];
 delete [] parray;
    
}


int Sfem::nchooser(int n, int r)
{
  int p;
  int q;
  if (r<n-r) {
    p=r;
    q=n-r;
  }
  else {
    p=n-r;
    q=r;
  }

  int num=1;
  for (int i=q+1;i<=n;i++)  num=num*i;
  int den = 1;
  for (int i=2; i<=p;i++)  den = den* i;
  int nd =num/den;
  return nd;
}


