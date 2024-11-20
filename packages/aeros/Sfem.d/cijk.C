#include <cmath>
#include <cstdlib>
#include <vector>
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Sfem.d/cijk.h>

Cijk::Cijk(int nd, int od, int P) : PC1(nd, od)
{
  ndim = nd; 
  order = od;
  maxorder = 11;
  maxnterms = P;
//  if(P != PC1.tnterms()) cerr << " *** WARNING Cijk, P = " << P << ", PC1.tnterms() = " << PC1.tnterms() << endl;
//  cerr << " ndim is  = " << ndim << " order is  " << od << endl;
}

double Cijk::expectation(int i, int j, int k)
{
  //Chaos PC1(ndim, order);

  //maxnterms = PC1.tnterms();

  if (i>maxnterms)
    printf ("ERROR index exceedes specified Expansion\n");
  if (j>maxnterms)
    printf ("ERROR index exceedes specified Expansion\n");
  if (k>maxnterms)
    printf ("ERROR index exceedes specified Expansion\n");
  
  GenFullM<double> icoef(maxorder+1);
  GenFullM<double> Orthorder(maxnterms+1, ndim+1);

  PC1.HermiteToPowers(maxorder, icoef);
  PC1.OrthoOrder(Orthorder);


  //int A[ndim+1];
  //int B[ndim+1];
  //int C[ndim+1];
  std::vector<int> A(ndim+1);
  std::vector<int> B(ndim+1);
  std::vector<int> C(ndim+1);

  for(int id=1; id<=ndim; id++)
    {
      A[id] = int(Orthorder[i][id]);
      B[id] = int(Orthorder[j][id]);
      C[id] = int(Orthorder[k][id]);
      if (i==0)
	A[id] = 0;
      if (j==0)
	B[id] =0;
      if (k==0)
	C[id] = 0;
    }

  double prod = 1;

  //double Mu[ndim+1];
  //double Sigma[ndim+1];
  std::vector<double> Mu(ndim+1);
  std::vector<double> Sigma(ndim+1);

  for(int in=1; in<=ndim; in++)
    {
      Mu[in] = 0;
      Sigma[in] = 1;
    }

  for(int dm=1; dm<=ndim; dm++)
    {
      double m = Mu[dm]; 
      double s = Sigma[dm];

      double Momt[40];
      Momt[0] = 1;
      Momt[1] = m;
      Momt[2] = pow(m,2) + pow(s,2);
      for(int mg = 3; mg<40; mg++)
	Momt[mg] = m*Momt[mg-1]+((mg-1)*pow(s,2)*Momt[mg-2]);

      GenVector<double> A1(maxorder);
      GenVector<double> B1(maxorder); 
      GenVector<double> C1(maxorder);

      GenVector<double> P1(2*maxorder);
      GenVector<double> P2(3*maxorder);

      for(int in=1; in<=maxorder; in++)
	{
	  A1[in-1] = (icoef[A[dm]+1][in]);
	  B1[in-1] = (icoef[B[dm]+1][in]);
	  C1[in-1] = (icoef[C[dm]+1][in]);
	}

      PC1.mult1dHermite(A1,B1,P1);
      PC1.mult1dHermite(P1,C1,P2);

      double p1=0;
      for(int ig=0; ig<3*maxorder; ig++)
	p1 += P2[ig]*Momt[ig];

      prod *= p1;
    } 

  return prod;

}
