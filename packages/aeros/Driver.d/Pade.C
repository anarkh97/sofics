#include <Utils.d/DistHelper.h>

template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::pade(VecType *sol, VecType **u, double *h, double x)
{
  // note u0[k+1] is the kth derivate of u0, similarly for uh
  int i, j, k, n, r;
  // PJSA 12-06-04: pade for general case (single or multipoint)
  if(a == 0) {
    // first time, compute P and Q coeffs (a, b)
    int l = domain->solInfo().getSweepParams()->padeL;
    int m = domain->solInfo().getSweepParams()->padeM;
    int padeN = domain->solInfo().getSweepParams()->padeN;
    int nRHS = domain->solInfo().getSweepParams()->nFreqSweepRHS;
    filePrint(stderr, " ... Computing %d-point Pade coefficients (l = %d, m = %d) ... \n", padeN, l, m);
    // allocate storage for P coefficients
    ia = l+1;
    a = new VecType * [ia];
    for(i = 0; i < ia; ++i) a[i] = new VecType(probDesc->solVecInfo());
    // allocate storage for Q coefficients
    ib = m+1;
    b = new VecType * [ib];
    for(i = 0; i < ib; ++i) b[i] = new VecType(probDesc->solVecInfo());
    // compute P, Q coefficients by solving system Ax = b
    int dim = l+m+1;
    GenFullM<Scalar> A(dim); // LHS
    Scalar *v = new Scalar[dim]; // RHS
    for(i=0; i<sol->size(); ++i) {
      A.zero();
      // assemble
      for(j=0; j<padeN; ++j) {
        int offset = j*(nRHS+1)+1;
        int I;
        for(n=0; n<nRHS; ++n) {
          if((I = j*nRHS+n) >= dim) break;
          // fill left block (a coefficients) of RHS
          for(k=n; k<=l; ++k)
            A[I][k] = -double(DFactorial(k)/DFactorial(k-n)*pow(h[j],k-n));
          // fill right block (b coefficients) of RHS
          for(r=0; r<=n; ++r)
            for(k=r; k<=m; ++k)
              if(k>0) A[I][l+k] += double(DCombination(n,r)*DFactorial(k)/DFactorial(k-r)*pow(h[j],k-r))*(*u[offset+n-r])[i];
          // fill RHS
          v[I] = -(*u[offset+n])[i];
        }
      }
      if(domain->solInfo().getSweepParams()->pade_pivot) { // factor and solve with full pivoting
        double tolerance = domain->solInfo().getSweepParams()->pade_tol; // default is 1.0e-12
        A.Factor(tolerance);
        A.ReSolve(v);
      }
      else { // factor and solve without pivoting
        A.factor();
        A.reSolve(v);
      }
      // extract coefficients
      for(j=0; j<=l; ++j) (*a[j])[i] = v[j];
      (*b[0])[i] = 1.0; // b0 = 1
      for(j=1; j<=m; ++j) (*b[j])[i] = v[l+j];
    }
    delete [] v;
    // allocate P, Q (will be reused)
    P = new VecType(probDesc->solVecInfo());
    Q = new VecType(probDesc->solVecInfo());
  }
                                                                                                                                                    
  // assemble P, Q and solution P/Q
  P->zero();
  Q->zero();
  for(i = 0; i < ia; ++i) P->linAdd(pow(x,i),*a[i]);
  for(i = 0; i < ib; ++i) Q->linAdd(pow(x,i),*b[i]);
  for(j=0; j<sol->size(); ++j) (*sol)[j] = (*P)[j] / (*Q)[j];  // this should be added to DistrVector and Vector classes
}


extern "C" {
void _FORTRAN(pade)(double *c, int &ic, double *a, int &ia, double *b, int &ib,
                    int &l, int &m, double &x, double *w, double *w1, double *w2, int *ik, int &iw);
void _FORTRAN(zpade)(DComplex *c, int &ic, DComplex *a, int &ia, DComplex *b, int &ib,
                    int &l, int &m, double &x, DComplex *w, DComplex *w1, DComplex *w2, int *ik, int &iw);
}
inline void Tpade(double *c, int &ic, double *a, int &ia, double *b, int &ib,
                    int &l, int &m, double &x, double *w, double *w1, double *w2, int *ik, int &iw)
{ 
  _FORTRAN(pade)(c,ic,a,ia,b,ib,l,m,x,w,w1,w2,ik,iw); 
}
inline void Tpade(DComplex *c, int &ic, DComplex *a, int &ia, DComplex *b, int &ib,
                  int &l, int &m, double &x, DComplex *w, DComplex *w1, DComplex *w2, int *ik, int &iw)
{ 
  _FORTRAN(zpade)(c,ic,a,ia,b,ib,l,m,x,w,w1,w2,ik,iw);
}

template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::pade1(VecType *sol, VecType **sol_prev, double x)
{
  int i, j;
  // PJSA 10-29-04: 1 point pade extrapolation using Charbel's fortran code see Driver.d/pade.f
  if(a == 0) {
    filePrint(stderr," ... Computing 1-point Pade coefficients ... \n");
    // first time, compute P and Q coeffs (a, b)
    int l = domain->solInfo().getSweepParams()->padeL;
    int m = domain->solInfo().getSweepParams()->padeM;
    // allocate storage for P coefficients
    ia = l+1;
    a = new VecType * [ia];
    for(i = 0; i < ia; ++i) a[i] = new VecType(probDesc->solVecInfo());
    // allocate storage for Q coefficients
    ib = m+1;
    b = new VecType * [ib];
    for(i = 0; i < ib; ++i) b[i] = new VecType(probDesc->solVecInfo());
    // allocate storage for series coefficients
    int ic = l+m+1;
    VecType **c = new VecType * [ic];
    for(i = 0; i < ic; ++i) c[i] = new VecType(probDesc->solVecInfo());
    // allocate temporary storage
    int iw = m;
    Scalar *w = (Scalar *) dbg_alloca(iw*iw*sizeof(Scalar));
    Scalar *w1 = (Scalar *) dbg_alloca(iw*sizeof(Scalar));
    Scalar *w2 = (Scalar *) dbg_alloca(iw*sizeof(Scalar));
    int *ik = (int *) dbg_alloca(iw*sizeof(int));
    // fill series coefficients: F(x0+deltax) = F(x0) + 1/1! F'(x0) deltax + 1/2! F''(x0) deltax^2 + ...
    for(i=0; i<ic; ++i) {
      c[i]->zero();
      c[i]->linAdd(1.0/DFactorial(i), *sol_prev[i+1]);
    }
    // compute P, Q coefficients
    Scalar *ai = (Scalar *) dbg_alloca(ia*sizeof(Scalar));
    Scalar *bi = (Scalar *) dbg_alloca(ib*sizeof(Scalar));
    Scalar *ci = (Scalar *) dbg_alloca(ic*sizeof(Scalar));
    for(i=0; i<sol->size(); ++i) {
      for(j = 0; j < ic; ++j) ci[j] = (*(c[j]))[i];  // extract series coef for dof i
      Tpade(ci,ic,ai,ia,bi,ib,l,m,x,w,w1,w2,ik,iw);
      for(j = 0; j < ia; ++j) (*(a[j]))[i] = ai[j];  // insert P coef for dof i
      for(j = 0; j < ib; ++j) (*(b[j]))[i] = bi[j];  // insert Q coef for dof i
    }
    // delete series coefs
    for(i=0; i<ic; ++i) delete c[i];
    delete [] c;
    // allocate P, Q (will be reused)
    P = new VecType(probDesc->solVecInfo());
    Q = new VecType(probDesc->solVecInfo());
  }

  // assemble P, Q and solution P/Q
  P->zero();
  Q->zero();
  for(i = 0; i < ia; ++i) P->linAdd(pow(x,i),*a[i]);
  for(i = 0; i < ib; ++i) Q->linAdd(pow(x,i),*b[i]);
  for(j=0; j<sol->size(); ++j) (*sol)[j] = (*P)[j] / (*Q)[j];  // this should be added to DistrVector and Vector classes
}

template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::fourier(VecType *sol, VecType **u, double *h, double x)
{
  int padeN = domain->solInfo().getSweepParams()->padeN; 
  double x0 = (h[0]+h[padeN-1])/2;
  double X = x-x0;
  double L = (h[padeN-1]-h[0])*2.0; // the 2.0 factor makes L the "effective" period
  ComplexD ii=ComplexD(0.0, 1.0);
  DComplex tpiol = 2.0*PI*ii/L;
  int l = domain->solInfo().getSweepParams()->padeL;
  int m = domain->solInfo().getSweepParams()->padeM;
  // note u0[k+1] is the kth derivate of u0, similarly for uh
  int i, j, k, n, r;
  // PJSA 12-06-04: fourier-pade for general case (single or multipoint)
  if(ca == 0) {
    // first time, compute P and Q coeffs (a, b)
    double *H = new double [domain->solInfo().getSweepParams()->padeN];
    for(i=0; i<padeN; ++i) H[i] = h[i]-x0;
    int nRHS = domain->solInfo().getSweepParams()->nFreqSweepRHS;
    filePrint(stderr, " ... Computing %d-point Fourier-Pade coefficients (l = %d, m = %d) ... \n", padeN, l, m);
    // allocate storage for P coefficients
    ia = l+1;
    ca = new ComplexVecType * [ia];
    for(i = 0; i < ia; ++i) ca[i] = new ComplexVecType(probDesc->solVecInfo());
    // allocate storage for Q coefficients
    ib = m+1;
    cb = new ComplexVecType * [ib];
    for(i = 0; i < ib; ++i) cb[i] = new ComplexVecType(probDesc->solVecInfo());
    // compute P, Q coefficients by solving system Ax = b
    int dim = l+m+1;
    GenFullM<DComplex> A(dim); // LHS
    DComplex *v = new DComplex[dim]; // RHS
    for(i=0; i<sol->size(); ++i) { 
      A.zero();
       // assemble
      for(j=0; j<padeN; ++j) {
        int offset = j*(nRHS+1)+1;
        int row,col;
        for(n=0; n<nRHS; ++n) {
          col = 0;
          if((row = j*nRHS+n) >= dim) break;
          // fill left block (a coefficients) of RHS
          for(k=-l/2; k<=(l+1)/2; ++k) {
            col = k + l/2;
            A[row][col] = -pow(tpiol*double(k),n)*exp(tpiol*double(k)*H[j]); 
          }
          // fill right block (b coefficients) of RHS
          for(r=0; r<=n; ++r)
            for(k=-m/2; k<=(m+1)/2; ++k) {
              col = l + (k+m/2);
              if(k<0) col += 1;
              if(k!=0) {
                A[row][col] += double(DCombination(n,r))*(*u[offset+n-r])[i]*pow(tpiol*double(k),r)*exp(tpiol*double(k)*H[j]);
              }
            }
          // fill RHS
          v[row] = -(*u[offset+n])[i];
        }
      }
      if(domain->solInfo().getSweepParams()->pade_pivot) { // factor and solve with full pivoting
        double tolerance = domain->solInfo().getSweepParams()->pade_tol; // default is 1.0e-12
        A.Factor(tolerance); 
        A.ReSolve(v);
      }
      else { // factor and solve without pivoting
        A.factor(); 
        A.reSolve(v);
      }

      // extract coefficients
      for(j=0; j<=l; ++j) (*ca[j])[i] = v[j];
      for(j=0; j<m/2; ++j) (*cb[j])[i] = v[l+j+1];
      (*cb[m/2])[i] = 1.0;
      for(j=m/2+1; j<=m; ++j) (*cb[j])[i] = v[l+j];
    }
    delete [] v;
    // allocate P, Q (will be reused)
    cP = new ComplexVecType(probDesc->solVecInfo());
    cQ = new ComplexVecType(probDesc->solVecInfo());
    delete [] H;
  }
                                                                                                                                      
  // assemble P, Q and solution P/Q
  cP->zero();
  cQ->zero();
  i=0;
  for(k=-l/2; k<=(l+1)/2; ++k) cP->linAdd(exp(k*X*tpiol),*ca[i++]);
  i = 0;
  for(k=-m/2; k<=(m+1)/2; ++k) cQ->linAdd(exp(k*X*tpiol),*cb[i++]);
  for(j=0; j<sol->size(); ++j) ScalarTypes::copy((*sol)[j], (*cP)[j]/(*cQ)[j]);  
}

