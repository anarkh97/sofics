#include <cmath>
#include <Math.d/FullSquareMatrix.h>
#include <Element.d/Element.h>
#include <Corotational.d/BarThermalCorotator.h>
#include <Corotational.d/TemperatureState.h> 

BarThermalCorotator::BarThermalCorotator(int _n1, int _n2, double __P, 
                                         double _eps, double _sigma, double _Tr, CoordSet& cs)
{
 n1    = _n1;         	// Node 1
 n2    = _n2;   	// Node 2
 P     = __P;    	// Perimeter
 eps   = _eps;          // Emissivity of the body
 sigma = _sigma;        // Stefan's constant
 Tr    = _Tr;           // Temperature of the enclosure receiving the radiation

 // Get original coordinates of bar's nodes
 Node &node1 = cs.getNode(n1);
 Node &node2 = cs.getNode(n2);

 // Compute original length of bar element
 double dx = node2.x - node1.x;
 double dy = node2.y - node1.y;
 double dz = node2.z - node1.z;
 l0 = sqrt(dx*dx + dy*dy + dz*dz);
}

void
BarThermalCorotator::getStiffAndForce(GeomState &ts, CoordSet &cs, 
                               FullSquareMatrix &elK, double *f, double dt, double t)
/*******************************************************************
 *
 * Purpose :
 *  Compute Tangent Stiffness matrix 
 *  for bar thermal element in current configuration.
 *
 * Input Variables:
 *  ts        : current temperature 
 *  cs        : coordinate set, contains reference configuration
 *
 * Local Variables:
 * xn : temperature in the current configuration
 *
 * Output :
 * elK      : element tangent stiffness matrix
 * f        : element internal force
 *
 *****************************************************************/
{
 // Declare local variables 
 int    i, j, k;
 double kt[2][2], ff[2], xn[2];

 // Get current Node State
 NodeState &tn1 = ts[n1];
 NodeState &tn2 = ts[n2];

 // Set temperature of Cn configuration 
 xn[0] = tn1.x; // temperature of node 1
 xn[1] = tn2.x; // temperature of node 2

 // Form tangent stiffness matrix
 formTangentStiffness( xn, l0, P, eps, sigma, kt);

 // Copy tangent stiffness matrix to element K matrix
 for(i=0; i<2; ++i)
   for(j=0; j<2; ++j)
     elK[i][j] = kt[i][j];
 
 // Form internal force
 formInternalForce( xn, l0, P, eps, sigma, Tr, ff);

 // Copy internal force to element f matrix
 for(k=0; k<2; ++k)
   f[k]=ff[k];
}

void
BarThermalCorotator::getInternalForce(GeomState &ts, CoordSet &cs, 
                                      FullSquareMatrix &elK, double *f, double dt, double t)
/*******************************************************************
 *
 * Purpose :
 *  Compute internal force vector
 *  for bar thermal element in current configuration.
 *
 * Input Variables:
 *  ts        : current temperature 
 *  cs        : coordinate set, contains reference configuration
 *
 * Local Variables:
 * xn : temperature in the current configuration
 *
 * Output :
 * f        : element internal force
 *
 *****************************************************************/
{
 // Declare local variables 
 int    i, j, k;
 double ff[2], xn[2];

 // Get current Node State
 NodeState &tn1 = ts[n1];
 NodeState &tn2 = ts[n2];

 // Set temperature of Cn configuration 
 xn[0] = tn1.x; // temperature of node 1
 xn[1] = tn2.x; // temperature of node 2

 // Form internal force
 formInternalForce( xn, l0, P, eps, sigma, Tr, ff);

 // Copy internal force to element f matrix
 for(k=0; k<2; ++k)
   f[k]=ff[k];
}

void
BarThermalCorotator::formInternalForce(double xn[2], double l0, double P, double eps, double sigma, double Tr, double *f) 
/*******************************************************************
 *
 * Purpose :
 *  Compute internal force vector for thermal bar
 *  element in current configuration.
 *
 * Output :
 *  f        : internal force vector in current configuration
 *
 *****************************************************************/
{
     f[0]=0.0;
     f[1]=0.0;

     // Compute internal force in local system and store in f (analytically derived)
     double coeff = (eps*sigma*P*l0)/pow(2.0,6.0);
     f[0] = coeff*(pow(xn[0],4.0)*10.666667+4*pow(xn[0],3.0)*xn[1]*2.133333+6*pow(xn[0],2.0)*pow(xn[1],2.0)*1.066667+4*xn[0]*pow(xn[1],3.0)*1.066667+pow(xn[1],4.0)*2.133333)-pow(2.0,5.0)*coeff*pow(Tr,4.0);
     f[1] = coeff*(pow(xn[0],4.0)*2.133333+4*pow(xn[0],3.0)*xn[1]*1.066667+6*pow(xn[0],2.0)*pow(xn[1],2.0)*1.066667+4*xn[0]*pow(xn[1],3.0)*2.133333+pow(xn[1],4.0)*10.666667)-pow(2.0,5.0)*coeff*pow(Tr,4.0);

}


void
BarThermalCorotator::formTangentStiffness(double xn[2], double l0, double P, double eps, double sigma, double kt[2][2])
/*******************************************************************
 * 
 * Purpose :
 *  Compute tangential stiffness for bar thermal element
 *  in current configuration.
 *
 * Input :
 *  xn       : current nodal temperature vector
 *  l0       : original length of the element
 *  P        : perimeter of the element
 *  eps      : emissivity of the body
 *  
 * Output :
 *  kt       : tangent stiffness matrix in current configuration
 *
 *****************************************************************/
{
     int i, j;
  // Zero stiffness matrix
     for(i=0; i<2; ++i )
       for(j=0; j<2; ++j )
         kt[i][j] = 0.0;

  // Tangent stiffness matrix (analytically derived)
 kt[0][0]=(4*eps*sigma*l0*P/pow(2.0,6.0))*(pow(xn[0],3.0)*10.6666666667+3*pow(xn[0],2.0)*xn[1]*2.1333333333+3*xn[0]*pow(xn[1],2.0)*1.0666666667+pow(xn[1],3.0)*1.0666666667);
 kt[0][1]=(4*eps*sigma*l0*P/pow(2.0,6.0))*(pow(xn[0],3.0)*2.1333333333+3*pow(xn[0],2.0)*xn[1]*1.0666666667+3*xn[0]*pow(xn[1],2.0)*1.0666666667+pow(xn[1],3.0)*1.0666666667);
 kt[1][0]=(4*eps*sigma*l0*P/pow(2.0,6.0))*(pow(xn[0],3.0)*2.1333333333+3*pow(xn[0],2.0)*xn[1]*1.0666666667+3*xn[0]*pow(xn[1],2.0)*1.0666666667+pow(xn[1],3.0)*1.0666666667);
 kt[1][1]=(4*eps*sigma*l0*P/pow(2.0,6.0))*(pow(xn[0],3.0)*1.0666666667+3*pow(xn[0],2.0)*xn[1]*1.0666666667+3*xn[0]*pow(xn[1],2.0)*2.1333333333+pow(xn[1],3.0)*10.6666666667);
}
