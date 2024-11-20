#include <cmath>
#include <iostream>

#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Element.d/Element.h>
#include <Corotational.d/BarFCorotator.h>
#include <Corotational.d/GeomState.h>
#include <Corotational.d/utilities.h>

BarFCorotator::BarFCorotator(int _n1, int _n2, double _e,
                             double _lambda,
                             double _a0, 
                             int _op, double _h, double _d,
                             double _Uc, double _Uf,
                             int _np, int _Nf, double _dlambda, int _Seed,
                             double _preload, CoordSet& cs)
{
 n1 = _n1; // Node 1
 n2 = _n2; // Node 2
 op = _op; // Material Option
 Ucrit = _Uc; // Stretch where Damage Initiates
 Uf = _Uf; // Failure Stretch for the Yarn
 // Assign Material Values
 if (op == 1)
 {
   // Micro-Scale Simulation
   double h, d, lambda_g, ef;
   int np, Nf;
   ef = _e;
   lambda_g = _lambda;
   h = _h;
   d = _d;
   np = _np;
   Nf = _Nf;
   AssignMicroScaleProp(ef, h, d, np, Nf, lambda_g);
 }
 else if (op == 2)
 {
   // Define Using Distributions
   double lam_mean, ey, de, dlambda, lam_slope, lam_intercept;
   ey = _e;
   de = _d;
   // Calculate Young's Modulus
   em = RndNorm(ey, de);
   // Calculate Damage Parameter
   lam_slope = _lambda;
   lam_intercept = _h;
   dlambda = _dlambda;
   lam_mean = lam_slope*em + lam_intercept;
   lambda = RndNorm(lam_mean, dlambda);
 }
 else
 {
   // Assign directly
   em = _e; // Elastic modulus
   lambda = _lambda; // Damage growth parameter
 }

 a0 = _a0; // Original Cross-sectional area
 preload = _preload;  // Preload in Truss

 // Get original coordinates of bar's nodes
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);

 // Compute original length of bar element
 double dx = node2.x - node1.x;
 double dy = node2.y - node1.y;
 double dz = node2.z - node1.z;
 l0 = sqrt(dx*dx + dy*dy + dz*dz);

 // Initialize Parameters
 damage = 1.0;
 //lambda = 150.0;
  
 if (_Seed == -1)
 {
   // Display Material Values
   std::cerr << "Young's Modulus = " << em << "; Lambda = " << lambda << std::endl;
 }
}

void BarFCorotator::AssignMicroScaleProp(double ef, double h, double d, int np,
                                         int Nf, double lambda_g)
/*******************************************************************
 *
 * Purpose :
 *  Preform a micro-scale level simulation to determine the macro-
 *  scale material properties 
 *
 * Input Variables:
 * ef       : elastic modulus of the fibrils
 * Ucrit    : critical stretch of the fibrils
 * Uf       : final stretch for the yarn (used for total failure)
 * h        : initial height of the yarn
 * d        : standard deviation for the inclined distance of the fibril
 * np       : number of points used in curve fitting the damage parameter lambda
 * Nf       : number of fibrils in a single yarn
 * lambda_g : initial guess for the damage parameter (for newton iteration)
 *
 * Local Variables:
 * em       : Elastic modulus
 * lambda   : damage parameter
 *
 *****************************************************************/
{
  int i, j;
  double delta_step, delta, U, IE_dam;
  double *pU_vec, *pd, *pEps, *pDamage;
  // If np < 1, set to default value of 10
  if (np < 1)
    np = 10;
  pd = new double[Nf];
  pEps = new double[Nf];
  pU_vec = new double[np+2];
  pDamage = new double[np+2];
  // Populate array of fibril inclination values
  for (i = 0; i < Nf; i++)
    *(pd+i) = RndNorm(0.0, d);
  em = 0.0; // Initialize Sum
  for (i = 0; i < Nf; i++)
    em += (ef/Nf)*pow(h,3)/pow(pow(h,2) + pow(*(pd+i),2),1.5);
  // Initialize Failure Flag for the Fibrils
  for (i = 0; i < Nf; i++)
    *(pEps+i) = 1.0;
  delta_step = (Uf - Ucrit)/(np+1); // size of the stretch steps
  // Run a micro-scale simulation - stretch the fibrils by displacing one end of the yarn
  for (i = 0; i < (np+1); i++)
  {
    *(pU_vec+i) = Ucrit + i*delta_step; // library of stretch values
    delta = h*(Ucrit + i*delta_step - 1.0); // displacement of upper endpoint of the yarn
    IE_dam = 0.0; // initialize sum
    for (j = 0; j < Nf; j++)
    {
      // Determine the stretch in the jth fibril
      U = sqrt(pow(h+delta,2) + pow(*(pd+j),2))/sqrt(pow(h,2) + pow(*(pd+j),2));
      if (U > Ucrit)
        *(pEps+j) = 0.0; // Fibril has failed
      IE_dam += (ef/Nf)*pow(h,3)/pow(pow(h,2) + pow(*(pd+i),2),1.5)* *(pEps+j);
    }
    *(pDamage+i) = IE_dam/ef;
  }
  // Set Final Damage and Stretch Value
  *(pU_vec+np+1) = Uf;
  *(pDamage+np+1) = 0.0;
  delete [] pd;
  delete [] pEps;
  // Calculate the Damage Parameter
  MicroCalcLambda(np, pDamage, pU_vec, lambda_g);
  delete [] pU_vec;
  delete [] pDamage;
  return;
}

void
BarFCorotator::MicroCalcLambda(int np, double *pDamage, double *pU_vec, double lambda_g)
/*******************************************************************
 * 
 * Purpose :
 *  Determine the damage parameter by doing a curve-fit
 *
 * Input :
 *  np       : number of points used in the curve fit not including the two end points
 *  pDamage  : damage values, y values for the curve
 *  pU_vec   : stretch values, x values for the curve
 *  Ucrit    : stretch at which damage first initiates
 *  Uf       : stretch at which element fails
 *  lambda_g : initial guess for newton iteration
 *
 * Output :
 *  lambda   : damage parameter
 *
 *****************************************************************/
{
  double a, b, tol, dr, ddr;
  double alpha, alpha2, U;
  double da_dl, d2a_dl2;
  int i, k, k_max;
  lambda = lambda_g; // Initial Guess
  b = Uf - Ucrit;
  tol = 1e-8; // tolerance on newton iteration
  k_max = 2*np; // maximum number of newton steps
  // Netwon steps
  for (k = 0; k < k_max; k++){
    dr = 0; // initialize sum
    ddr = 0; // initialize sum
    for (i = 0; i < (np+2); i++){
      // Input value of U 
      U = *(pU_vec+i);
      // Input value of alpha2 - value I am trying to match
      alpha2 = *(pDamage+i);
      a = U - Ucrit;
      alpha = (exp(-lambda*a) - exp(-lambda*b))/(1 - exp(-lambda*b));
      da_dl = (-a*exp(-lambda*a)+b*exp(-lambda*b))/(1-exp(-lambda*b))-(exp(-lambda*a)-exp(-lambda*b))/pow((1-exp(-lambda*b)),2)*b*exp(-lambda*b);
      d2a_dl2 = (a*a*exp(-lambda*a)-b*b*exp(-lambda*b))/(1-exp(-lambda*b))-2*(-a*exp(-lambda*a)+b*exp(-lambda*b))/pow((1-exp(-lambda*b)),2)*b*exp(-lambda*b)+2*(exp(-lambda*a)-exp(-lambda*b))/pow((1-exp(-lambda*b)),3)*b*b*pow(exp(-lambda*b),2)+(exp(-lambda*a)-exp(-lambda*b))/pow((1-exp(-lambda*b)),2)*b*b*exp(-lambda*b);
      // sum over each point
      dr = dr + 2*(alpha - alpha2)*da_dl;
      ddr = ddr + 2*(da_dl*da_dl) + 2*(alpha - alpha2)*d2a_dl2;
    }
    // calculate new value of lambda
    lambda = lambda - dr/ddr;
    if (fabs(dr) < tol)
      break;
  }
}

double
BarFCorotator::RndNorm(double mean, double stdev)
/*******************************************************************
 * 
 * Purpose :
 *  Return a random number according to a normal (gaussian) distribution
 *
 * Input :
 *  mean  : mean value
 *  stdev : standard deviation
 *
 *****************************************************************/
{
  double U1, U2, V1, V2, X, Y, S = 2.0;
  Y = (double) RAND_MAX;
  while (S >= 1){
    U1 = rand()/Y;
    U2 = rand()/Y;
    V1 = 2.0*U1 - 1.0;
    V2 = 2.0*U2 - 1.0;
    S = pow(V1,2.0) + pow(V2,2.0);
  }
  X = V1*sqrt( (-2.0*log(S))/S );
  return (mean + stdev*X);
}

void
BarFCorotator::getStiffAndForce(GeomState &geomState, CoordSet &cs, 
                                FullSquareMatrix &elK, double *f, double dt, double)
/*******************************************************************
 *
 * Purpose :
 *  Compute Internal force vector and Tangent Stiffness matrix 
 *  for Co-rotated bar element in current configuration.
 *
 * Input Variables:
 *  geomState : current configuration 
 *  cs        : coordinate set, contains reference configuration
 *
 * Local Variables:
 * x0[i][j] : global coordinate component j of node i in reference configuration
 * xn[i][j] : global coordinate component j of node i in current configuration
 * t[i]     : transformation matrix (1st vector only) for current state
 * t0[i]    : transformation matrix (1st vector only) for initial state
 * l0       : original length 
 * ld       : deformed length
 * a0       : Cross-sectional area
 * em       : Elastic modulus
 * e        : Green-Lagrange (GL) strain
 * sigma    : PK2 axial stress
 * p        : axial force in local coordinates
 * preload  : axial preload
 * f0       : internal force vector in initial configuration due to preload
 *
 * Output :
 *  f        : internal force vector in current configuration
 *  elK      : element tangent stiffness matrix
 *
 *****************************************************************/
{
 // Declare local variables 
 int    i;
 double xn[2][3], x0[2][3], t[3], t0[3], f0[6];

 // Get original coordinates of bar's nodes
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);

 // Get current Node State
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];

 // Set coordinates of C0 configuration
 x0[0][0] = node1.x;
 x0[0][1] = node1.y;
 x0[0][2] = node1.z;
 x0[1][0] = node2.x;
 x0[1][1] = node2.y;
 x0[1][2] = node2.z;

 // Set coordinates of Cn configuration 
 xn[0][0] = ns1.x; // xn coordinate of node 1
 xn[0][1] = ns1.y; // yn coordinate of node 1
 xn[0][2] = ns1.z; // zn coordinate of node 1
 xn[1][0] = ns2.x; // xn coordinate of node 2
 xn[1][1] = ns2.y; // yn coordinate of node 2
 xn[1][2] = ns2.z; // zn coordinate of node 2

 // Compute deformed length: ld
 double dx = xn[1][0] - xn[0][0];
 double dy = xn[1][1] - xn[0][1];
 double dz = xn[1][2] - xn[0][2];
 double ld = sqrt(dx*dx + dy*dy + dz*dz);

 // Form transformation tensor: t (1st vector only)
 // for current state 
 t[0] = dx/ld;
 t[1] = dy/ld;
 t[2] = dz/ld;

 // Form transformation tensor: t0 (1st vector only)
 // for initial state 
 dx = x0[1][0] - x0[0][0];
 dy = x0[1][1] - x0[0][1];
 dz = x0[1][2] - x0[0][2];
 t0[0] = dx/l0;
 t0[1] = dy/l0;
 t0[2] = dz/l0;

 // Compute current GL-strain
 //double e = (ld - l0)/l0;
 double e = 0.5*((ld/l0)*(ld/l0) - 1);

 // Compute damage
 double alpha;
 if (e < 0.03) 
   alpha = 1.0;
 else
   alpha = (exp(-lambda*(ld/l0 - Ucrit)) - exp(-lambda*(Uf - Ucrit)))/(1.0 - exp(-lambda*(Uf - Ucrit)));
 //alpha = exp(-lambda*(e - .03));
 // Set upper limit on damage
 if (alpha < damage)
    damage = alpha;
 if (damage < 0.0)
    damage = 0.0;
 // Compute current PK2-stress
 double sigma = damage*em*e;
 // Enforce zero stress on compression
 if (e < 0.0)
    sigma = 0.0;
 // Compute current axial force: p
 double p = (ld/l0)*sigma*a0;

 // Add Preload
 p += preload;

 // Form current internal force: f
 formInternalForce(t, p, f);

 // Form initial internal force: f0
 formInternalForce(t0, preload, f0);

 for(i=0; i<6; ++i)
   f[i] -= f0[i];

 elK.zero();  // this element is only used for nonlinear explicit
 //damage += 1.0;
}

//----------------------------------------------------------------------

void
BarFCorotator::getExternalForce(GeomState &geomState, CoordSet &cs, double *f)
{
 // Declare local variables 
 int    i, j;
 double xn[2][3], t[3];
 
 // Get current Node State
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];
 
 // Set coordinates of Cn configuration 
 xn[0][0] = ns1.x; // xn coordinate of node 1
 xn[0][1] = ns1.y; // yn coordinate of node 1
 xn[0][2] = ns1.z; // zn coordinate of node 1
 xn[1][0] = ns2.x; // xn coordinate of node 2
 xn[1][1] = ns2.y; // yn coordinate of node 2
 xn[1][2] = ns2.z; // zn coordinate of node 2

 // Compute deformed length: ld
 double dx = xn[1][0] - xn[0][0];
 double dy = xn[1][1] - xn[0][1];
 double dz = xn[1][2] - xn[0][2];
 double ld = sqrt(dx*dx + dy*dy + dz*dz);

 // Form transformation tensor: t (1st vector only)
 // for current state 
 t[0] = dx/ld;
 t[1] = dy/ld;
 t[2] = dz/ld;
 
 // build force
 f[3] = t[0]*f[0];
 f[4] = t[1]*f[0];
 f[5] = t[2]*f[0];

 f[0] = -f[3];
 f[1] = -f[4];
 f[2] = -f[5];
}

//----------------------------------------------------------------------

void 
BarFCorotator::getDExternalForceDu(GeomState &geomState, CoordSet &cs, 
                                   FullSquareMatrix &elK, double *locF)
{
 // Declare local variables 
 int    i, j;
 double xn[2][3], t[3];
 
 // Get current Node State
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];
    
 // Set coordinates of Cn configuration 
 xn[0][0] = ns1.x; // xn coordinate of node 1
 xn[0][1] = ns1.y; // yn coordinate of node 1
 xn[0][2] = ns1.z; // zn coordinate of node 1
 xn[1][0] = ns2.x; // xn coordinate of node 2
 xn[1][1] = ns2.y; // yn coordinate of node 2
 xn[1][2] = ns2.z; // zn coordinate of node 2

 // Compute deformed length: ld
 double dx = xn[1][0] - xn[0][0];
 double dy = xn[1][1] - xn[0][1];
 double dz = xn[1][2] - xn[0][2];
 double ld = sqrt(dx*dx + dy*dy + dz*dz);

 // Form transformation tensor: t (1st vector only)
 // for current state 
 t[0] = dx/ld;
 t[1] = dy/ld;
 t[2] = dz/ld;

 double p = -locF[0];
 double dummyK[3][3];
  
 dummyK[0][0] = p*(t[1]*t[1] + t[2]*t[2])/ld;
 dummyK[0][1] = -p*t[0]*t[1]/ld;
 dummyK[0][2] = -p*t[0]*t[2]/ld;
 dummyK[1][0] = -p*t[1]*t[0]/ld;
 dummyK[1][1] = p*(t[0]*t[0] + t[2]*t[2])/ld;
 dummyK[1][2] = -p*t[1]*t[2]/ld;
 dummyK[2][0] = -p*t[2]*t[0]/ld;
 dummyK[2][1] = -p*t[2]*t[1]/ld;
 dummyK[2][2] = p*(t[0]*t[0] + t[1]*t[1])/ld;
 
 // Fill element stiffness matrix
 for(i=0; i<3; i++) {
   for(j=0; j<3; j++) {
      elK[  i][  j] += dummyK[i][j];
      elK[  i][3+j] -= dummyK[i][j];
      elK[3+i][  j] -= dummyK[i][j];
      elK[3+i][3+j] += dummyK[i][j];
   }
 } 

}

void
BarFCorotator::formInternalForce(double t[3], double p, double *f)
/*******************************************************************
 *
 * Purpose :
 *  Compute internal force vector for Co-rotated bar
 *  element in current configuration.
 *
 * Input :
 *  t[i]     : transformation matrix (1st vector only)
 *  p        : axial force in local coordinates
 *
 * Output :
 *  f        : internal force vector in current configuration
 *
 *****************************************************************/
{
  // Compute internal force in local system and store in f
  f[0] =   p;
  f[1] = 0.0;
  f[2] = 0.0;

  // Transform to global coordinate by Fg = T'*Fl and store in f
  // Shortened form, since f[1] = f[2] = 0.0
  f[3] = t[0]*f[0];
  f[4] = t[1]*f[0];
  f[5] = t[2]*f[0];

  f[0] = -f[3];
  f[1] = -f[4];
  f[2] = -f[5];
}

void
BarFCorotator::formGeometricStiffness(GeomState &geomState, CoordSet &cs,
                                      FullSquareMatrix &kg, double *f)
{
/*******************************************************************
 *
 * Purpose :
 *  Compute Geometric stiffness for Co-rotated bar element
 *  in current configuration.
 *
 * Input :
 * geomState: current state of nodes (coordinates, rotation tensors)
 *
 * Local variables:
 *  t        : transformation matrix (1st vector only)
 *  p        : axial force in local coordinates
 *  l0       : undeformed (initial) length of bar
 *  ld       : deformed length
 *  e        : current GL-strain
 *  sigma    : current PK-2 stress in bar
 *  p        : current axial force in bar
 *  a0       : initial area of bar
 *
 * Output :
 *  kg       : geometric tangent stiffness matrix in current configuration
 *
 *****************************************************************/
 // Declare local variables
 int    i, j;
 double xn[2][3],t[3];
 
 // Get current Node State
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];

 // Set coordinates of Cn configuration
 xn[0][0] = ns1.x; // xn coordinate of node 1
 xn[0][1] = ns1.y; // yn coordinate of node 1
 xn[0][2] = ns1.z; // zn coordinate of node 1
 xn[1][0] = ns2.x; // xn coordinate of node 2
 xn[1][1] = ns2.y; // yn coordinate of node 2
 xn[1][2] = ns2.z; // zn coordinate of node 2

 // Compute deformed length: ld
 double dx = xn[1][0] - xn[0][0];
 double dy = xn[1][1] - xn[0][1];
 double dz = xn[1][2] - xn[0][2];
 double ld = sqrt(dx*dx + dy*dy + dz*dz);

 // Form transformation tensor: t (1st vector only)
 // for current state 
 t[0] = dx/ld;
 t[1] = dy/ld;
 t[2] = dz/ld;

 // Compute current GL-strain
 double e = (ld - l0)/l0;

 // Compute current PK2-stress
 double sigma = em*e;

 // Compute current axial force: p
 double p = sigma*a0;

 // Add Preload
 p += preload;

 // Zero stiffness matrix
 for(i=0; i<6; ++i )
   for(j=0; j<6; ++j )
     kg[i][j] = 0.0;

 kg[0][0] += p*(t[1]*t[1] + t[2]*t[2])/ld;
 kg[0][1] += -p*t[0]*t[1]/ld;
 kg[0][2] += -p*t[0]*t[2]/ld;
 kg[1][0] += -p*t[1]*t[0]/ld;
 kg[1][1] += p*(t[0]*t[0] + t[2]*t[2])/ld;
 kg[1][2] += -p*t[1]*t[2]/ld;
 kg[2][0] += -p*t[2]*t[0]/ld;
 kg[2][1] += -p*t[2]*t[1]/ld;
 kg[2][2] += p*(t[0]*t[0] + t[1]*t[1])/ld;


 // Fill element stiffness matrix
 for(i=0; i<3; i++) {
   for(j=0; j<3; j++) {
     kg[  i][3+j] = -kg[i][j];
     kg[3+i][  j] =  kg[i][3+j];
     kg[3+i][3+j] =  kg[i][j];
   }
 }
}

void
BarFCorotator::extractDeformations(GeomState &geomState, CoordSet &cs, 
                                  double *vld, int &nlflag)
{
 // Set Flag to Use Linear Routines for Stress
 nlflag = 1;

 // Get original coordinates of bar's nodes
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);

 double xn[2][3], x0[2][3];

 // Set coordinates of C0 configuration
 x0[0][0] = node1.x;
 x0[0][1] = node1.y;
 x0[0][2] = node1.z;
 x0[1][0] = node2.x;
 x0[1][1] = node2.y;
 x0[1][2] = node2.z;

 // Get current Node State
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];

 // Set coordinates of Cn configuration
 xn[0][0] = ns1.x; // xn coordinate of node 1
 xn[0][1] = ns1.y; // yn coordinate of node 1
 xn[0][2] = ns1.z; // zn coordinate of node 1
 xn[1][0] = ns2.x; // xn coordinate of node 2
 xn[1][1] = ns2.y; // yn coordinate of node 2
 xn[1][2] = ns2.z; // zn coordinate of node 2

 // Compute deformed length: ld
 double dx = xn[1][0] - xn[0][0];
 double dy = xn[1][1] - xn[0][1];
 double dz = xn[1][2] - xn[0][2];
 double ld = sqrt(dx*dx + dy*dy + dz*dz);
 double delta_l = ld-l0;

 dx = x0[1][0] - x0[0][0];
 dy = x0[1][1] - x0[0][1];
 dz = x0[1][2] - x0[0][2];

 // scale dx, dy, and dz by the initial length
 dx /= l0;
 dy /= l0;
 dz /= l0;

 // Transform to global coordinate by  dXg = T'*dXl
 // Shortened form, since def[1] = def[2] = 0.0
 vld[3] = dx*delta_l/2.0;
 vld[4] = dy*delta_l/2.0;
 vld[5] = dz*delta_l/2.0;
 vld[0] = -vld[3];
 vld[1] = -vld[4];
 vld[2] = -vld[5]; 
}

void
BarFCorotator::extractRigidBodyMotion(GeomState &geomState, CoordSet &cs,
                                      double *vlr)
{
}

double
BarFCorotator::getElementEnergy(GeomState &geomState, CoordSet &cs)
{
// Computes Internal Energy of Element in Given State

 // Get original coordinates of bar's nodes
 //Node &node1 = cs.getNode(n1);
 //Node &node2 = cs.getNode(n2);

 // Get current Node State
 NodeState &ns1 = geomState[n1];
 NodeState &ns2 = geomState[n2];

 // Set coordinates of Cn configuration
 double xn[2][3];
 xn[0][0] = ns1.x; // xn coordinate of node 1
 xn[0][1] = ns1.y; // yn coordinate of node 1
 xn[0][2] = ns1.z; // zn coordinate of node 1
 xn[1][0] = ns2.x; // xn coordinate of node 2
 xn[1][1] = ns2.y; // yn coordinate of node 2
 xn[1][2] = ns2.z; // zn coordinate of node 2

 // Compute deformed length: ld
 double dx = xn[1][0] - xn[0][0];
 double dy = xn[1][1] - xn[0][1];
 double dz = xn[1][2] - xn[0][2];
 double ld = sqrt(dx*dx + dy*dy + dz*dz);

 // Compute current GL-strain
 //double e = (ld - l0)/l0;
 double e = 0.5*((ld/l0)*(ld/l0) - 1);

 // Compute damage
 double alpha;
 if (e < 0.03)
    alpha = 1.0;
 else
    alpha = (exp(-lambda*(ld/l0 - Ucrit)) - exp(-lambda*(Uf - Ucrit)))/(1.0 - exp(-lambda*(Uf - Ucrit)));
 // Set upper limit on damage
 if (alpha < damage)
    damage = alpha;
 if (damage < 0.0)
    damage = 0.0;

 // Compute current PK2-stress
 double sigma = damage*em*e;
 // Enforce zero stress on compression
 if (e < 0.0)
    sigma = 0.0;

 // Add Preload????
 sigma += preload/a0;

 // Compute Energy as 1/2 Integral[e*sigma*dV]
 double Energy = 0.5*e*sigma*a0*l0;

 return Energy;
}
