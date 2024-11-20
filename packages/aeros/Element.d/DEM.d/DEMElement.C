#include <cstdio>
#include <cstdlib>
#include <alloca.h>
#include <cmath>
#include <ctime>
#include <complex>

#include <Element.d/Helm.d/IsoParamUtils.h>
#include <Element.d/Helm.d/GaussRules.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>

#include <Element.d/DEM.d/DEMElement.h>
#include <Driver.d/Domain.h>
extern Domain *domain;
//#else
//
//#ifdef F_NEEDS_UNDSC
//#define _FORTRAN(a) a ## _
//#else
//#define _FORTRAN(a) a
//#endif
//
//#ifdef _STANDARD_C_PLUS_PLUS
// #include <complex>
// using std::complex;
//#else
//#include <complex.h>
//#endif
//#ifdef COMPLEX_NON_TEMPLATE
//typedef complex ComplexD;
//typedef complex DComplex;
//#else
//typedef complex<double> ComplexD;
//typedef complex<double> DComplex;
//#endif
//
//#include "IsoParamUtils.h"
//#include "GaussRules.h"
//#endif
#include "DEMElement.h"

#include <Driver.d/GeoSource.h>
extern GeoSource *geoSource;

using namespace std;

DEMInterfaceElement::DEMInterfaceElement(DEMElement *_deme, DEMElement *_deme2, int _fi) {
 deme = _deme;
 deme2 = _deme2;
 fi = _fi;
}

#ifndef SANDIA
// --------------------------------------------------------------------------
// Element functions

void DEMElement::renum(const int *table) {
 fprintf(stderr,"DEMElement::renum is not implemented.\n");
}


void DEMElement::renum(EleRenumMap& table) {
 fprintf(stderr,"DEMElement::renum is not implemented.\n");
}


void DEMInterfaceElement::renum(const int *table) {
 fprintf(stderr,"DEMInterfaceElement::renum is not implemented.\n");
}

void DEMInterfaceElement::renum(EleRenumMap& table) {
 fprintf(stderr,"DEMInterfaceElement::renum is not implemented.\n");
}


int* DEMElement::nodes(int *p) const {

 if (p == 0) p = new int[numNodes()];
 int i,ii=0;
 if (!dgmFlag())
   for(i=0;i<nGeomNodes();i++) p[ii++] = nn[i];
 for(i=0;i<nFaces();i++) {
   if (lm[i] != 0) {
     int nlm = lm[i]->nDofs();
     for(int j=0;j<nlm;j++) p[ii++] = lm[i]->getNodeOffset()+j;
   }
 }
 if (!condensedFlag()) for(i=0;i<nEnrichmentDofs();i++) p[ii++] = ne+i;

 return p;
}


int* DEMInterfaceElement::nodes(int *p) const {

 if (p == 0) p = new int[numNodes()];

 int ii=0;
 if (!deme->dgmFlag()) 
   for(int i=0;i<deme->nGeomNodes();i++) p[ii++] = deme->nn[i];
 for(int i=0;i<deme->nEnrichmentDofs();i++) p[ii++] = deme->ne+i;
 if (!deme2->dgmFlag()) 
   for(int i=0;i<deme2->nGeomNodes();i++) p[ii++] = deme2->nn[i];
 for(int i=0;i<deme2->nEnrichmentDofs();i++) p[ii++] = deme2->ne+i;

 return p;
}


int* DEMElement::dofs(DofSetArray &dsa, int *p) const  {

 if (p == 0) p = new int[numDofs()];
 int ii=0;
 if (!dgmFlag()) {
   for(int i=0;i<nGeomNodes();i++) { 
     dsa.number(nn[i],polyDofType(),p+ii);
     ii += polyDofsPerNode();
   }
 }

 for(int i=0;i<nFaces();i++) {
   if (lm[i] != 0) {
     int j;
     for(j=0;j<lm[i]->nDofs();j++) {
        dsa.number(lm[i]->getNodeOffset()+j, DofSet::LagrangeE, p+ii++);
       
     }
   }
 }
 
 if (!condensedFlag()) 
// RT: Temporarily put DofSet::Helm, but should create some other new name
   for(int i=0;i<nEnrichmentDofs();i++) dsa.number(ne+i,DofSet::Helm,p+ii++);
 return p;
}


int* DEMInterfaceElement::dofs(DofSetArray &dsa, int *p) const  {

 if (p == 0) p = new int[numDofs()];
 int ii=0;

 if (!deme->dgmFlag()) {
   for(int i=0;i<deme->nGeomNodes();i++) {
     dsa.number(deme->nn[i], deme->polyDofType(),p+ii);
     ii += deme->polyDofsPerNode();
   }
 }

// RT: Temporarily put DofSet::Helm, but should create some other new name
 for(int i=0;i<deme->nEnrichmentDofs();i++)
   dsa.number(deme->ne+i,DofSet::Helm,p+ii++);

 if (!deme2->dgmFlag()) {
   for(int i=0;i<deme2->nGeomNodes();i++) {
     dsa.number(deme2->nn[i], deme2->polyDofType(),p+ii);
     ii += deme2->polyDofsPerNode();
   }
 }

// RT: Temporarily put DofSet::Helm, but should create some other new name
 for(int i=0;i<deme2->nEnrichmentDofs();i++)
   dsa.number(deme2->ne+i,DofSet::Helm,p+ii++);

 return p;
}


void DEMElement::markDofs(DofSetArray &dsa) const {

 if (!dgmFlag())
   for(int i=0;i<nGeomNodes();i++) { dsa.mark(nn[i],polyDofType());
 }

 for(int i=0;i<nFaces();i++) {
   if (lm[i] != 0) {
     int j;
     for(j=0;j<lm[i]->nDofs();j++) {
        dsa.mark(lm[i]->getNodeOffset()+j, DofSet::LagrangeE);
     }
   }
 }

 if (!condensedFlag()) 
// RT: Temporarily put DofSet::Helm, but should create some other new name
   for(int i=0;i<nEnrichmentDofs();i++) dsa.mark(ne+i,DofSet::Helm);
}


void DEMInterfaceElement::markDofs(DofSetArray &dsa) const {

 if (!deme->dgmFlag())
   for(int i=0;i<deme->nGeomNodes();i++) dsa.mark(deme->nn[i],deme->polyDofType());

// RT: Temporarily put DofSet::Helm, but should create some other new name
 for(int i=0;i<deme->nEnrichmentDofs();i++) dsa.mark(deme->ne+i,DofSet::Helm);

 if (!deme2->dgmFlag())
   for(int i=0;i<deme2->nGeomNodes();i++) dsa.mark(deme2->nn[i],deme2->polyDofType());

// RT: Temporarily put DofSet::Helm, but should create some other new name
 for(int i=0;i<deme2->nEnrichmentDofs();i++) dsa.mark(deme2->ne+i,DofSet::Helm);
}


FullSquareMatrix DEMElement::massMatrix(const CoordSet &cs, double *K, int fl) const {
 fprintf(stderr,"DEMElement::massMatrix not implemented.\n");
 FullSquareMatrix ret(0,K);
 return ret;
}


FullSquareMatrix DEMElement::stiffness(const CoordSet &cs, double *K, int flg ) const {
 fprintf(stderr,"DEMElement::massMatrix not implemented.\n");
 FullSquareMatrix ret(0,K);
 return ret;
}


int DEMElement::nodesE(int no) {
 ne = no;
 return condensedFlag()?0:nEnrichmentDofs();
}

double DEMElement::getOmega() const { return geoSource->omega(); }

double *DEMElement::getWaveDirection() const { return domain->getWaveDirection(); }

double DEMElement::getSpeedOfSound() const
{
    return geoSource->omega()/prop->kappaHelm;
}

void
DEMElement::getNodalCoord(int n, int *nn, double *xyz) const
{
    (domain->getNodes()).getCoordinates(nn,n,xyz,xyz+n,xyz+2*n);
}

void DEMElement::nodalSolution(CoordSet &cs, complex<double> *sol,
                          complex<double> (*nodalSol)[8]) {

 complex<double> *a;
 if (condensedFlag()) {
   int n = nPolynomialDofs()+nLagrangeDofs()+nEnrichmentDofs();
   int nc, *condensed, *kept;
   condensedDofs(nc,condensed,kept);
   if (!storeMatricesFlag()) {
     complex<double> *kk;
     kk = new complex<double>[n*n];
     createM(kk);
     complex<double> *A = new complex<double>[nc*nc];
     createCondensedMatrices(kk,A);
     delete[] A;
     delete[] kk;
   }
   a = new complex<double>[nc]; 
   complex<double>* c = new complex<double>[n-nc]; 
   for(int i=0;i<n-nc;i++) c[i] = sol[kept[i]];
   delete[] kept;
   delete[] condensed;
   staticCondensationSolution(nc, n-nc, demm.C, demm.B,
                              c, demm.f, a);
   delete[] c;
   if (!storeMatricesFlag()) {
     delete[] demm.B;
     delete[] demm.C;
   }
 } else {
   a = sol + nPolynomialDofs()+nLagrangeDofs();
 }
 int j=0;
 for(int i=0;i<nGeomNodes();i++) {
   double xyz[3];
   cs.getCoordinates(nn+i,1,xyz,xyz+1,xyz+2);
   createSol(xyz, a, nodalSol[i]);
   if (!dgmFlag()) {
      if (polyDofType()& DofSet::Xdisp) { 
        nodalSol[i][1] += sol[j];
        j++;
      }
      if (polyDofType()& DofSet::Ydisp)  {
        nodalSol[i][2] += sol[j];
        j++;
      }
      if (polyDofType()& DofSet::Zdisp) {
        nodalSol[i][3] += sol[j];
        j++;
      }
      if (polyDofType()& DofSet::Helm) {
        nodalSol[i][0] += sol[j];
        j++;
      }
   }
 }
 if (condensedFlag()) delete[] a;
}


// -----------------------------------------------------------------------
#endif


int DEMCoreElement::nLagrangeDofs() const {

 int nl = 0;
 for(int i=0;i<nFaces();i++)
   if (lm[i] != 0) nl += lm[i]->nDofs();
 return nl;
}


void DEMCoreElement::condensedDofs(int &nc, int *&condensed, int *&kept) {
 int n = nPolynomialDofs()+nLagrangeDofs()+nEnrichmentDofs();
 nc = nEnrichmentDofs();
 condensed = new int[nc];
 kept = new int[n-nc];
 int i;
 for(i=0;i<nc;i++)
   condensed[i] = nPolynomialDofs()+nLagrangeDofs()+i;
 for(i=0;i<n-nc;i++)
   kept[i] = i;
}


int DEMCoreElement::createCondensedMatrices(complex<double>* kk, 
                                         complex<double>* A) {

 int nc, *condensed, *kept;
 int n = nPolynomialDofs()+nLagrangeDofs()+nEnrichmentDofs();
 condensedDofs(nc,condensed,kept);
 demm.B = new complex<double>[nc*(n-nc)];
 demm.C = new complex<double>[nc*nc];
 int i,j;
 for(i=0;i<nc;i++) for(j=0;j<nc;j++)
   demm.C[i+j*nc] = kk[condensed[i]+n*condensed[j]];
 for(i=0;i<n-nc;i++) for(j=0;j<n-nc;j++)
   A[i+j*(n-nc)] = kk[kept[i]+n*kept[j]];
 for(i=0;i<n-nc;i++) for(j=0;j<nc;j++)
   demm.B[j+i*nc] = kk[kept[i]+n*condensed[j]];
 delete[] condensed;
 delete[] kept;

 staticCondensationLHS(nc,n-nc,demm.C,demm.B,A);
 return nc;
}


void DEMCoreElement::systemMatrix(complex<double> *K) {

 int n = nPolynomialDofs()+nLagrangeDofs()+nEnrichmentDofs();
 if (!condensedFlag()) {
   createM(K);
 } else {
   complex<double> *kk;
   kk = new complex<double>[n*n];

   createM(kk);
   int nc = createCondensedMatrices(kk,K);
   delete[] kk;
   
   if (!storeMatricesFlag()) {
     delete[] demm.B;
     delete[] demm.C;
   }
 }
}


void DEMCoreInterfaceElement::systemMatrix(complex<double> *K) {
 deme->interfMatrix(fi,deme2,K);
}


void DEMCoreElement::interfMatrix(int fi, DEMElement* deme2, complex<double> *K) {
 fprintf(stderr,"DEMElement::interfMatrix not implemented.\n");
}

void DEMCoreElement::enrichmentF(double *x, complex<double> *f) {
 fprintf(stderr,"DEMElement::enrichmentF is not implemented.\n");
}

void DEMCoreElement::polynomialF(double *x, double *f) {
 fprintf(stderr,"DEMElement::polynomialF is not implemented.\n");
}


void DEMCoreElement::systemRHS(complex<double> *v) {

 if (!condensedFlag()) {
   createRHS(v);
 } else {
   int n = nPolynomialDofs()+nLagrangeDofs()+nEnrichmentDofs();
   complex<double> *vv = new complex<double>[n];
   createRHS(vv);
   int nc, *condensed, *kept;
   condensedDofs(nc,condensed,kept);
   if (!storeMatricesFlag()) {
     complex<double> *kk;
     kk = new complex<double>[n*n];
     createM(kk);
     complex<double> *A = new complex<double>[nc*nc];
     createCondensedMatrices(kk,A);
     delete[] A;
     delete[] kk;
   }

   if (demm.f==0)
     demm.f = new complex<double>[nc];
   int i;
   for(i=0;i<nc;i++) demm.f[i] = vv[condensed[i]];
   for(i=0;i<n-nc;i++) v[i] = vv[kept[i]];
   delete[] condensed;
   delete[] kept;
   delete[] vv;
   staticCondensationRHS(nc, n-nc, 1, demm.C, demm.B, demm.f, v);
   if (!storeMatricesFlag()) {
     delete[] demm.B;
     delete[] demm.C;
   }
 }
}


extern "C" {
void _FORTRAN(zgesvx)( const char &fact, const char & trans, const int &n,
                       const int &nrhs, ComplexD *a, const int &lda,
                       ComplexD *af, const int &ldaf, int *ipiv,
                       char &equed, double *r, double *c,
                       ComplexD *b, const int &ldb,
                       ComplexD *x, const int &ldx,
                       double &rcond,
                       double *ferr, double *berr,
                       ComplexD  *work, double *rwork, int &info);
}

extern "C" {
void _FORTRAN(zgesv)
(const int &n, const int &nrhs,
                            ComplexD *a, const int &lda, int *ipiv,
                            ComplexD *b, const int &ldb, int &info);
}

extern "C" {
void _FORTRAN(zgetrf)
(const int &m, const int &n, ComplexD *a, const int &lda, int *ipiv, int &info);
}

extern "C" {
void _FORTRAN(zgetrs) (const char & trans, const int &n, const int &nrhs,
                            ComplexD *a, const int &lda, int *ipiv,
                            ComplexD *b, const int &ldb, int &info);
}



//#define RT_USE_ALLOCA
#define LAFVERSION 1

#if LAFVERSION == 0
void DEMElement::staticCondensationLHS(int ne, int nl, complex<double> *kee,
                         complex<double> *kel, complex<double> *kll) {
#ifdef RT_USE_ALLOCA
 complex<double> *af = (complex<double>*)alloca(sizeof(complex<double>)*ne*ne);
 complex<double> *x =  (complex<double>*)alloca(sizeof(complex<double>)*ne*nl);
#else
 complex<double> *af = new complex<double>[ne*ne];
 complex<double> *x = new complex<double>[ne*nl]; 
#endif

 int info;
 int lda = ne;
 int ldaf = ne;
 int ldb = ne;
 int ldx = ne;
#ifdef RT_USE_ALLOCA
 double *r = (double*)alloca(sizeof(double)*ne);
 double *c = (double*)alloca(sizeof(double)*ne);
 int *ipiv = (int*)alloca(sizeof(int)*ne);
 double *ferr = (double*)alloca(sizeof(double)*nl);
 double *berr = (double*)alloca(sizeof(double)*nl);
 complex<double> *work = (complex<double>*)alloca(sizeof(complex<double>)*2*ne);
 double *rwork = (double*)alloca(sizeof(double)*2*ne);
#else
 double *r = new double[ne];
 double *c = new double[ne];
 int *ipiv = new int[ne];
 double *ferr = new double[nl];
 double *berr = new double[nl];
 complex<double> *work =  new complex<double>[2*ne];
 double *rwork = new double[2*ne];
#endif
// double *r = 0, *c = 0;
 char fact = 'N';
 char trans = 'N';
 char equed = 'N';
 double rcond;
 _FORTRAN(zgesvx)(fact, trans, ne, nl, kee, lda, af, ldaf, ipiv,
                  equed, r, c, kel, ldb, x, ldx,
                  rcond, ferr, berr, work, rwork, info);
//  fprintf(stderr, " rcond = %e \n", rcond);
 if (rcond<1e-15) fprintf(stderr,
   "Nearly singular matrix: rcond = %e in DEMElement::staticCondensationLHS\n",
    rcond);
 if (info!=0) fprintf(stderr,
   "Singular matrix in DEMElement::staticCondensationLHS\n");

 char transa = 'T';
 char transb = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int ldc = nl;
 _FORTRAN(zgemm)( transa, transb, nl, nl, ne, alpha , kel, lda, x, ldx,
                  beta, kll, ldc);
#ifndef RT_USE_ALLOCA
 delete[] af;
 delete[] x;
 delete[] r;
 delete[] c;
 delete[] ipiv;
 delete[] ferr;
 delete[] berr;
 delete[] work;
 delete[] rwork;
#endif
}

#elif LAFVERSION == 1
void DEMCoreElement::staticCondensationLHS(int ne, int nl, complex<double> *kee,
                         complex<double> *kel, complex<double> *kll) {

 ipiv = new int[ne];

 int info;
 int lda = ne;
 _FORTRAN(zgetrf)(ne, ne, kee, lda, ipiv, info);
 if (info>0) fprintf(stderr,
   "Singular matrix in DEMElement::staticCondensationLHS\n");

 char trans = 'N';
 int ldb = ne;
 complex<double> *b = (complex<double>*)alloca(sizeof(complex<double>)*ne*nl);
 for(int i=0;i<ne*nl;i++) b[i] = kel[i];
 _FORTRAN(zgetrs)(trans, ne, nl, kee, lda, ipiv, b, ldb, info);
 if (info!=0) fprintf(stderr,
   "Error by zgetrs in DEMElement::staticCondensationLHS\n");

 char transa = 'T';
 char transb = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int ldc = nl;
 _FORTRAN(zgemm)( transa, transb, nl, nl, ne, alpha , kel, lda, b, ldb,
                  beta, kll, ldc);
}
#else
void DEMCoreElement::staticCondensationLHS(int ne, int nl, complex<double> *kee,
                         complex<double> *kel, complex<double> *kll) {

 int *ipiv = (int*)alloca(sizeof(int)*ne);
 complex<double> *a = (complex<double>*)alloca(sizeof(complex<double>)*ne*ne);
 int i;
 for(i=0;i<ne*ne;i++) a[i] = kee[i];
 complex<double> *b = (complex<double>*)alloca(sizeof(complex<double>)*ne*nl);
 for(i=0;i<ne*nl;i++) b[i] = kel[i];

 int info;
 int lda = ne;
 int ldb = ne;
 _FORTRAN(zgesv)(ne, nl, a, lda, ipiv, b, ldb, info);

 char transa = 'T';
 char transb = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int ldc = nl;
 _FORTRAN(zgemm)( transa, transb, nl, nl, ne, alpha , kel, lda, b, ldb,
                  beta, kll, ldc);
}
#endif


#if LAFVERSION == 0
void DEMCoreElement::staticCondensationRHS(int ne, int nl, int nr,
                         complex<double> *kee, complex<double> *kel,
                         complex<double> *re, complex<double>* rl) {

 complex<double> *af = (complex<double>*)alloca(sizeof(complex<double>)*ne*ne);
 complex<double> *x = (complex<double>*)alloca(sizeof(complex<double>)*ne*nr);

 int info;
 int lda = ne;
 int ldaf = ne;
 int ldb = ne;
 int ldx = ne;
 // double *r = (double*)alloca(sizeof(double)*ne);
 // double *c = (double*)alloca(sizeof(double)*ne);
 double *r = 0, *c = 0;
 char fact = 'N';
 char trans = 'N';
 char equed = 'N';
 double rcond;
 int *ipiv = (int*)alloca(sizeof(int)*ne);
 double *ferr = (double*)alloca(sizeof(double)*nr);
 double *berr = (double*)alloca(sizeof(double)*nr);
 complex<double> *work = (complex<double>*)alloca(sizeof(complex<double>)*2*ne);
 double *rwork = (double*)alloca(sizeof(double)*2*ne);
 _FORTRAN(zgesvx)(fact, trans, ne, nr, kee, lda, af, ldaf, ipiv,
                  equed, r, c, re, ldb, x, ldx,
                  rcond, ferr, berr, work, rwork, info);

 char transa = 'T';
 char transb = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int ldc = nl;
 _FORTRAN(zgemm)( transa, transb, nl, nr, ne, alpha , kel, lda, x, ldx,
                  beta, rl, ldc);
}

#elif LAFVERSION == 1
void DEMCoreElement::staticCondensationRHS(int ne, int nl, int nr,
                         complex<double> *kee, complex<double> *kel,
                         complex<double> *re, complex<double>* rl) {

 int info;
 char trans = 'N';
 int lda = ne;
 int ldb = ne;
 complex<double> *b = (complex<double>*)alloca(sizeof(complex<double>)*ne*nr);
 for(int i=0;i<ne*nr;i++) b[i] = re[i];
 _FORTRAN(zgetrs)(trans, ne, nr, kee, lda, ipiv, b, ldb, info);
 if (info!=0) fprintf(stderr,
   "Error by zgetrs in DEMElement::staticCondensationRHS\n");

 char transa = 'T';
 char transb = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int ldc = nl;
 _FORTRAN(zgemm)( transa, transb, nl, nr, ne, alpha , kel, lda, b, ldb,
                  beta, rl, ldc);
}

#else
void DEMCoreElement::staticCondensationRHS(int ne, int nl, int nr,
                         complex<double> *kee, complex<double> *kel,
                         complex<double> *re, complex<double>* rl) {

 int *ipiv = (int*)alloca(sizeof(int)*ne);
 complex<double> *a = (complex<double>*)alloca(sizeof(complex<double>)*ne*ne);
 int i;
 for(i=0;i<ne*ne;i++) a[i] = kee[i];
 complex<double> *b = (complex<double>*)alloca(sizeof(complex<double>)*ne*nr);
 for(i=0;i<ne*nr;i++) b[i] = re[i];

 int info;
 int lda = ne;
 int ldb = ne;
 _FORTRAN(zgesv)(ne, nr, a, lda, ipiv, b, ldb, info);

 char transa = 'T';
 char transb = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int ldc = nl;
 _FORTRAN(zgemm)( transa, transb, nl, nr, ne, alpha , kel, lda, b, ldb,
                  beta, rl, ldc);
}
#endif


#if LAFVERSION == 0
void DEMCoreElement::staticCondensationSolution(int ne, int nl,
                          complex<double> *kee, complex<double> *kel,
                          complex<double> *pl,
                          complex<double> *re, complex<double> *pe) {
 int i;
 for(i=0;i<ne;i++) pe[i] = re[i];

 char trans = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int incx = 1, incy = 1;
 _FORTRAN(zgemv) ( trans, ne, nl, alpha, kel, ne, pl, incx, beta, pe, incy);

 int nr = 1;

 complex<double> *af = (complex<double>*)alloca(sizeof(complex<double>)*ne*ne);
 complex<double> *x = (complex<double>*)alloca(sizeof(complex<double>)*ne*nr);
 for(int i=0;i<ne;i++) x[i] = pe[i];

 int info;
 int lda = ne;
 int ldaf = ne;
 int ldb = ne;
 int ldx = ne;
 // double *r = (double*)alloca(sizeof(double)*ne);
 // double *c = (double*)alloca(sizeof(double)*ne);
 double *r = 0, *c = 0;
 char fact = 'N';
 trans = 'N';
 char equed = 'N';
 double rcond;
 int *ipiv = (int*)alloca(sizeof(int)*ne);
 double *ferr = (double*)alloca(sizeof(double)*nr);
 double *berr = (double*)alloca(sizeof(double)*nr);
 complex<double> *work = (complex<double>*)alloca(sizeof(complex<double>)*2*ne);
 double *rwork = (double*)alloca(sizeof(double)*2*ne);
 _FORTRAN(zgesvx)(fact, trans, ne, nr, kee, lda, af, ldaf, ipiv,
                  equed, r, c, x, ldb, pe, ldx,
                  rcond, ferr, berr, work, rwork, info);
 if (rcond<1e-15) fprintf(stderr,
   "Nearly singular matrix: rcond = %e in DEMElement::staticCondensationSolution\n",
    rcond);
 if (info!=0) fprintf(stderr,
    "Singular matrix in DEMElement::staticCondensationSolution\n");
}

#elif LAFVERSION == 1
void DEMCoreElement::staticCondensationSolution(int ne, int nl,
                          complex<double> *kee, complex<double> *kel,
                          complex<double> *pl,
                          complex<double> *re, complex<double> *pe) {

 int i;
 for(i=0;i<ne;i++) pe[i] = re[i];

 char trans = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int incx = 1, incy = 1;
 _FORTRAN(zgemv) ( trans, ne, nl, alpha, kel, ne, pl, incx, beta, pe, incy);

 int info;
 int lda = ne;
 int ldb = ne;
 int nr = 1;

 _FORTRAN(zgetrs)(trans, ne, nr, kee, lda, ipiv, pe, ldb, info);
}

#else
void DEMCoreElement::staticCondensationSolution(int ne, int nl,
                          complex<double> *kee, complex<double> *kel,
                          complex<double> *pl,
                          complex<double> *re, complex<double> *pe) {

 int i;
 for(i=0;i<ne;i++) pe[i] = re[i];

 char trans = 'N';
 complex<double> alpha = -1.0;
 complex<double> beta = 1.0;
 int incx = 1, incy = 1;
 _FORTRAN(zgemv) ( trans, ne, nl, alpha, kel, ne, pl, incx, beta, pe, incy);

 int *ipiv = (int*)alloca(sizeof(int)*ne);
 complex<double> *a = (complex<double>*)alloca(sizeof(complex<double>)*ne*ne);
 for(i=0;i<ne*ne;i++) a[i] = kee[i];

 int info;
 int lda = ne;
 int ldb = ne;
 int nr = 1;
 _FORTRAN(zgesv)(ne, nr, a, lda, ipiv, pe, ldb, info);
}
#endif

#ifdef SANDIA

#include "DEMHelm3d.h"

hex8dem::hex8dem() {
 nodal_coords = 0;
 de = 0;
}


int hex8dem::NumEhancementDofs() const { return (de==0)?-1:de->nEnrichmentDofs(); }
int hex8dem::NumConstraintDofs()  const{ return (de==0)?-1:de->nLagrangeDofs(); }
int hex8dem::NumPolynomialDofs()  const{ return (de==0)?-1:de->nPolynomialDofs(); }


int hex8dem::Initialize(bool condense_enrichment_flag, // indicates whether condensation takes place
                   bool dgm_flag, // indicates whether the element has polynomial field
                   bool store_matrices_flag, // indicates whether element matrices should be stored
                   int element_enrichment_type, // type of enrichment
                   int* lagrange_multiplier_type, // type of lagrange multiplier
                   const double* _nodal_coords, // coordinates of element nodes
                   int * connectivity, // global node numbers of element nodes
                   double _omega, // omega = 2*pi*f, where f is the frequency
                   double _c0, // speed of sound
                   double _rho // density of fluid
                   ) {
 omega = _omega;
 c0 = _c0;
 rho = _rho;
 switch (element_enrichment_type) {
   case 1:
      if (dgm_flag) de = new DGMHelm3d_6(8,connectivity); 
      else de = new DEMHelm3d_6(8,connectivity);
      break;
   case 2:
      if (dgm_flag) de = new DGMHelm3d_26(8,connectivity); 
      else de = new DEMHelm3d_26(8,connectivity);
      break;
   case 3:
      if (dgm_flag) de = new DGMHelm3d_56(8,connectivity); 
      else de = new DEMHelm3d_56(8,connectivity);
      break;
   case 4:
      if (dgm_flag) de = new DGMHelm3d_98(8,connectivity); 
      else de = new DEMHelm3d_98(8,connectivity);
      break;
   default:
      fprintf(stderr,"Invalid enrichment type %d.\n",element_enrichment_type); return -1;
 }

 if (!store_matrices_flag) de->storeMatricesF = false;
 if (!condense_enrichment_flag) de->condensedF = false;

 nodal_coords = new double[8*3];
 for(int i=0;i<8*3;i++) nodal_coords[i] = _nodal_coords[i];

 for(int fi=0;fi<6;fi++) {
   DEMLM *lm=0;
   switch (lagrange_multiplier_type[fi]) {
     case 51:  lm = new DGMHelm3d_1_LM(); break;
     case 52:  lm = new DGMHelm3d_4_LM(); break;
     case 53:  lm = new DGMHelm3d_8_LM(); break;
     case 54:  lm = new DGMHelm3d_12_LM(); break;
     default: fprintf(stderr,"Invalid LM type %d.\n",lagrange_multiplier_type[fi]);
   }
   de->lm[fi] = lm;
   if (lm) de->lm[fi]->e1 = de;
 }
 return 0;
}

hex8dem::~hex8dem() {
  if (de)  for(int fi=0;fi<6;fi++) if (de->lm[fi]) delete de->lm[fi];
  if (de) delete de;
  if (nodal_coords) delete[] nodal_coords;
}

void hex8dem::ElemDynamicMtx(complex<double> *matrix) const {
 if (de==0) { matrix=0; return; }
 matrix = new complex<double>[de->numDofs()*de->numDofs()];
 de->systemMatrix(matrix);
}

void hex8dem::ScalarFaceIntegral(int faceid, const double *scalar_field,
                            complex<double> **vector) const {
}


#endif

