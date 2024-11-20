#ifndef _ARUBBERF_H_
#define _ARUBBERF_H_ 


class ARubberF {
public:
 double omega;
 complex<double> *lambda;
 complex<double> *mu;
 double E0;
 double dE;
 double mu0;
 double dmu;
 double eta_E;
 double deta_E;
 double eta_mu;
 double deta_mu;
 int nd;
 ARubberF() {}
 ARubberF(int _nd, double _omega,
                      double _E0, double _dE, double _mu0, double _dmu,
                      double _eta_E, double _deta_E,
                      double _eta_mu, double _deta_mu) {
   nd = _nd; omega = _omega;
   E0 = _E0; dE = _dE; mu0 = _mu0; dmu = _dmu;
   eta_E = _eta_E; deta_E = _deta_E; eta_mu = _eta_mu; deta_mu = _deta_mu;
//   lambda = nu*E/((1.0+nu)*(1.0-2.0*nu));
//   mu = E/(2.0*(1.0+nu));
/*
 complex<double> p[6] = { 3,2,1,-5,-10,1.2};
 complex<double> q[5] = { 4,-2,1,-1.2,2.3};
 double x = 5.0;
 complex<double> r[11];

 rational_derivatives(5,p,4,q,x,10,r);
 for(int i=0;i<10;i++) 
   fprintf(stderr,"%e\n",real(r[i]));

 complex<double> rr[10];
 mult(5,p,4,q,rr);
 for(int i=0;i<10;i++) 
   fprintf(stderr,"%e\n",real(rr[i]));
*/

// lamda = mu(E-2mu) / (3mu-E)
// mu = (mu0+dmu*f)*(1+i*(eta_mu+deta_mu*f));
//
   double f = omega/2.0/M_PI;
  
   complex<double> mu_n[3] = { mu0*complex<double>(1.0,-eta_mu),
                               dmu*complex<double>(1.0,-eta_mu)+
                               mu0*complex<double>(0.0,-deta_mu),
                               dmu*complex<double>(0.0,-deta_mu)};
   complex<double> mu_d[1]= {1.0 };
   mu= new complex<double>[nd+1];
//   rational_derivatives(2,mu_n,0,mu_d,f,nd,mu);
   rational_derivatives(2,mu_n,0,mu_d,omega,nd,mu);

   complex<double> E_n[3] = { E0*complex<double>(1.0,-eta_E),
                               dE*complex<double>(1.0,-eta_E)+
                               E0*complex<double>(0.0,-deta_E),
                               dE*complex<double>(0.0,-deta_E)};


   complex<double> lambda_x[3] = { -2.0*mu_n[0] + E_n[0], -2.0*mu_n[1] + E_n[1],
                                   -2.0*mu_n[2] + E_n[2]};
   complex<double> lambda_n[5];
   mult(2,mu_n,2,lambda_x,lambda_n);
   complex<double> lambda_d[3] = { 3.0*mu_n[0] - E_n[0], 3.0*mu_n[1] - E_n[1],
                                   3.0*mu_n[2] - E_n[2]};
   lambda = new complex<double>[nd+1];
//   rational_derivatives(4,lambda_n,2,lambda_d,f,nd,lambda);
   rational_derivatives(4,lambda_n,2,lambda_d,omega,nd,lambda);
 }

 complex<double> poly(int n_p, complex<double>*p, double x) {
  complex<double> y = p[n_p];
  for(int i=n_p;i>0;i--) {
    y = y*x+p[i-1];
  }
  return y;
 }
 
 void mult(int n_p, complex<double>* p,
           int n_q, complex<double>* q,
           complex<double> *r) {
 
 
 for(int i=0;i<=n_p+n_q;i++) r[i] = 0.0;
 for(int i=0;i<=n_p;i++) for(int j=0;j<=n_q;j++)
   r[i+j] += p[i]*q[j];
 }
 
 ~ARubberF() { 
    delete[] lambda;
    delete[] mu;
 } 
 void rational_derivatives(int n_p, complex<double>* p,
                           int n_q, complex<double>* q,
                           double x, int n, complex<double> *rd) {
 
  complex<double> *pd = new complex<double>[n_p+1];
  complex<double> *qd = new complex<double>[n_q+1];
 
  int nm = (n_p>n_q)? n_p:n_q;
  complex<double> *c = new complex<double>[nm+1];
 
 // Derivatives of p
  for(int i=0;i<=n_p;i++) c[i] = p[i];
  for(int i=0;i<=n_p;i++) {
    pd[i] = poly(n_p-i,c+i,x);
    for(int j=i+1;j<=n_p;j++) 
      c[j] *= double(j-i);
  }
 // Derivatives of q
  for(int i=0;i<=n_q;i++) c[i] = q[i];
  for(int i=0;i<=n_q;i++) {
    qd[i] = poly(n_q-i,c+i,x);
    for(int j=i+1;j<=n_q;j++) 
      c[j] *= double(j-i);
  }
 // Derivatives of r = p/q
  for(int i=0; i<=n; i++) {
    nm = (n_q<i)?n_q:i;
    double binomial = 1.0;
    rd[i] = (i<=n_p)?pd[i]:0.0;
    for(int j=1;j<=nm;j++) {
      binomial *= double(i-j+1)/double(j);
      rd[i] -= binomial*rd[i-j]*qd[j];
    }
    rd[i] /= qd[0];
  }
 
  delete[] pd;
  delete[] qd;
  delete[] c;
 }
 complex<double> d_lambda(int n) {
   return lambda[n];
 }
 complex<double> d_mu(int n) {
   return mu[n];
 }
};

class ARubberShell: public ARubberF {
public:
 complex<double> *c1;
 complex<double> *c2;
 complex<double> *c3;
 ARubberShell(int _nd, double _omega,
                      double _E0, double _dE, double _mu0, double _dmu,
                      double _eta_E, double _deta_E,
                      double _eta_mu, double _deta_mu) {
   nd = _nd; omega = _omega;
   E0 = _E0; dE = _dE; mu0 = _mu0; dmu = _dmu;
   eta_E = _eta_E; deta_E = _deta_E; eta_mu = _eta_mu; deta_mu = _deta_mu;

   double f = omega/2.0/M_PI;

// lamda = mu(E-2mu) / (3mu-E)
// mu = (mu0+dmu*f)*(1+i*(eta_mu+deta_mu*f));
// c1 = E/(1-nu^2) = 4mu^2/(4mu-E)
// c2 = nu*E(1-nu^2) = (2E-4mu)mu/(4mu-E)
// c3 = 2mu
  
   complex<double> mu_n[3] = { mu0*complex<double>(1.0,-eta_mu),
                               dmu*complex<double>(1.0,-eta_mu)+
                               mu0*complex<double>(0.0,-deta_mu),
                               dmu*complex<double>(0.0,-deta_mu)};
   complex<double> mu_d[1]= {1.0 };

   complex<double> E_n[3] = { E0*complex<double>(1.0,-eta_E),
                               dE*complex<double>(1.0,-eta_E)+
                               E0*complex<double>(0.0,-deta_E),
                               dE*complex<double>(0.0,-deta_E)};

   complex<double> c3_n[3] = { 2.0*mu_n[0], 2.0*mu_n[1], 2.0*mu_n[2] };
   complex<double> c3_d[1] = { 1.0 };
   c3 = new complex<double>[nd+1];
   rational_derivatives(2,c3_n,0,c3_d,omega,nd,c3);

   complex<double> c2_x[3] = { -4.0*mu_n[0] + 2.0*E_n[0],
                               -4.0*mu_n[1] + 2.0*E_n[1],
                               -4.0*mu_n[2] + 2.0*E_n[2] };
   complex<double> c2_n[5];
   mult(2,mu_n,2,c2_x,c2_n);
   complex<double> c2_d[3] = { 4.0*mu_n[0] - E_n[0], 4.0*mu_n[1] - E_n[1],
                                   4.0*mu_n[2] - E_n[2]};
   c2 = new complex<double>[nd+1];
   rational_derivatives(4,c2_n,2,c2_d,omega,nd,c2);

   complex<double> c1_x[3] = { 4.0*mu_n[0],
                               4.0*mu_n[1],
                               4.0*mu_n[2] };
   complex<double> c1_n[5];
   mult(2,mu_n,2,c1_x,c1_n);
   complex<double> c1_d[3] = { 4.0*mu_n[0] - E_n[0], 4.0*mu_n[1] - E_n[1],
                                   4.0*mu_n[2] - E_n[2]};
   c1 = new complex<double>[nd+1];
   rational_derivatives(4,c1_n,2,c1_d,omega,nd,c1);

 }
 complex<double> d_c1(int n) { return c1[n]; }
 complex<double> d_c2(int n) { return c2[n]; }
 complex<double> d_c3(int n) { return c3[n]; }
 ~ARubberShell() { 
    delete[] c1;
    delete[] c2;
    delete[] c3;
 } 
};

#endif
