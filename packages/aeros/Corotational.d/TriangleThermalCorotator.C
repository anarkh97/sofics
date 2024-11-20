#include <cmath>

#include <Math.d/FullSquareMatrix.h>
#include <Element.d/Element.h>
#include <Corotational.d/TriangleThermalCorotator.h>
#include <Corotational.d/TemperatureState.h> 

TriangleThermalCorotator::TriangleThermalCorotator(int _n1, int _n2, int _n3, double _A, 
                                                   double _eps, double _sigma, double _Tr, CoordSet& cs)
{
 n1    = _n1;         	// Node 1
 n2    = _n2;   	// Node 2
 n3    = _n3;           // Node 3
 A     = _A;    	// Area
 eps   = _eps;          // Emissivity of the body
 sigma = _sigma;	// Stefan's constant
 Tr    = _Tr;           // Temperature of the enclosure receiving the radiation
}

void
TriangleThermalCorotator::getStiffAndForce(GeomState &ts, CoordSet &cs, 
                                           FullSquareMatrix &elK, double *f, double dt, double t)
/*******************************************************************
 *
 * Purpose :
 *  Compute Tangent Stiffness matrix 
 *  for triangle thermal element in current configuration.
 *
 * Input Variables:
 *  ts        : current temperature 
 *  cs        : coordinate set, contains reference configuration
 *
 * Local Variables:
 * xn         : temperature in the current configuration
 *
 * Output :
 * elK      : element tangent stiffness matrix
 *
 *****************************************************************/
{
 // Declare local variables 
 int    i, j;
 double kt[3][3], ff[3], xn[3];

 // Get original coordinates of bar's nodes
 //Node &node1 = cs.getNode(n1);
 //Node &node2 = cs.getNode(n2);
 //Node &node3 = cs.getNode(n3);

 // Get current temperature state
 NodeState &tn1 = ts[n1];
 NodeState &tn2 = ts[n2];
 NodeState &tn3 = ts[n3];

 // Set temperature of Cn configuration 
 xn[0] = tn1.x; // temperature of node 1
 xn[1] = tn2.x; // temperature of node 2
 xn[2] = tn3.x; // temperature of node 3
 // Form tangent stiffness matrix
 formTangentStiffness( xn, A, eps, sigma, kt);

 // Copy tangent stiffness matrix to element K matrix
 for(i=0; i<3; ++i)
   for(j=0; j<3; ++j)
     elK[i][j] = kt[i][j];
 
 // Form internal force
 formInternalForce( xn, A, eps, sigma, Tr, ff);

 // Copy internal force to element f matrix
 for(j=0; j<3; ++j)
   f[j] = ff[j];
}

void
TriangleThermalCorotator::getInternalForce(GeomState &ts, CoordSet &cs, 
                                           FullSquareMatrix &elK, double *f, double dt, double t)
/*******************************************************************
 *
 * Purpose :
 *  Compute internal force vector
 *  for triangle thermal element in current configuration.
 *
 * Input Variables:
 *  ts        : current temperature 
 *  cs        : coordinate set, contains reference configuration
 *
 * Local Variables:
 * xn         : temperature in the current configuration
 *
 * Output :
 * elK      : element tangent stiffness matrix
 *
 *****************************************************************/
{
 // Declare local variables 
 int    i, j;
 double ff[3], xn[3];

 // Get current temperature state
 NodeState &tn1 = ts[n1];
 NodeState &tn2 = ts[n2];
 NodeState &tn3 = ts[n3];

 // Set temperature of Cn configuration 
 xn[0] = tn1.x; // temperature of node 1
 xn[1] = tn2.x; // temperature of node 2
 xn[2] = tn3.x; // temperature of node 3

 // Form internal force
 formInternalForce( xn, A, eps, sigma, Tr, ff);

 // Copy internal force to element f matrix
 for(j=0; j<3; ++j)
   f[j] = ff[j];
}

void
TriangleThermalCorotator::formInternalForce(double xn[3], double A, double eps, double sigma, double Tr, double f[3])
/*******************************************************************
 *
 * Purpose :
 *  Compute internal force vector for thermal triangle
 *  element in current configuration.
 *
 * Output :
 *  f        : internal force vector in current configuration
 *
 *****************************************************************/
{
  f[0]=0;
  f[1]=0;
  f[2]=0;
  
  double coeff = (A*eps*sigma)/630;

  // Compute internal force in local system and store in f (analytically derived)
  f[0]=coeff*(pow(xn[0],4.0)*30+pow(xn[1],4.0)*6+pow(xn[2],4.0)*6+4*pow(xn[0],3.0)*xn[1]*6+4*pow(xn[0],3.0)*xn[2]*6+4*xn[0]*pow(xn[1],3.0)*3+4*pow(xn[1],3.0)*xn[2]*1.5+4*xn[0]*pow(xn[2],3.0)*3+4*xn[1]*pow(xn[2],3.0)*1.5+6*pow(xn[0],2.0)*pow(xn[1],2.0)*3+6*pow(xn[0],2.0)*pow(xn[2],2.0)*3+6*pow(xn[1],2.0)*pow(xn[2],2.0)+12*pow(xn[0],2.0)*xn[1]*xn[2]*1.5+12*xn[0]*pow(xn[1],2.0)*xn[2]+12*xn[0]*xn[1]*pow(xn[2],2.0))-210*coeff*pow(Tr,4.0);
  f[1]=coeff*(pow(xn[0],4.0)*6+pow(xn[1],4.0)*30+pow(xn[2],4.0)*6+4*pow(xn[0],3.0)*xn[1]*3+4*pow(xn[0],3.0)*xn[2]*1.5+4*xn[0]*pow(xn[1],3.0)*6+4*pow(xn[1],3.0)*xn[2]*6+4*xn[0]*pow(xn[2],3.0)*1.5+4*xn[1]*pow(xn[2],3.0)*3+6*pow(xn[0],2.0)*pow(xn[1],2.0)*3+6*pow(xn[0],2.0)*pow(xn[2],2.0)+6*pow(xn[1],2.0)*pow(xn[2],2.0)*3+12*pow(xn[0],2.0)*xn[1]*xn[2]+12*xn[0]*pow(xn[1],2.0)*xn[2]*1.5+12*xn[0]*xn[1]*pow(xn[2],2.0))-210*coeff*pow(Tr,4.0);
  f[2]=coeff*(pow(xn[0],4.0)*6+pow(xn[1],4.0)*6+pow(xn[2],4.0)*30+4*pow(xn[0],3.0)*xn[1]*1.5+4*pow(xn[0],3.0)*xn[2]*3+4*xn[0]*pow(xn[1],3.0)*1.5+4*pow(xn[1],3.0)*xn[2]*3+4*xn[0]*pow(xn[2],3.0)*6+4*xn[1]*pow(xn[2],3.0)*6+6*pow(xn[0],2.0)*pow(xn[1],2.0)+6*pow(xn[0],2.0)*pow(xn[2],2.0)*3+6*pow(xn[1],2.0)*pow(xn[2],2.0)*3+12*pow(xn[0],2.0)*xn[1]*xn[2]+12*xn[0]*pow(xn[1],2.0)*xn[2]+12*xn[0]*xn[1]*pow(xn[2],2.0)*1.5)-210*coeff*pow(Tr,4.0);
}


void
TriangleThermalCorotator::formTangentStiffness(double xn[3], double A, double eps, double sigma, double kt[3][3])
/*******************************************************************
 * 
 * Purpose :
 *  Compute tangential stiffness for triangle thermal element
 *  in current configuration.
 *
 * Input :
 *  xn       : current nodal temperature vector
 *  A        : area of the element
 *  eps      : emissivity of the body
 *  
 * Output :
 *  kt       : tangent stiffness matrix in current configuration
 *
 *****************************************************************/
{
     int i, j;
     double coeff = (4*eps*sigma*A)/630;
  // Zero stiffness matrix
     for(i=0; i<3; ++i )
       for(j=0; j<3; ++j )
         kt[i][j] = 0.0;

  // Tangent stiffness matrix (analytically derived)
 kt[0][0]=coeff*(pow(xn[0],3.0)*30+pow(xn[1],3.0)*3+pow(xn[2],3.0)*3+3*pow(xn[0],2.0)*xn[1]*6+3*pow(xn[0],2.0)*xn[2]*6+3*xn[0]*pow(xn[1],2.0)*3+3*pow(xn[1],2.0)*xn[2]+3*xn[0]*pow(xn[2],2.0)*3+3*xn[1]*pow(xn[2],2.0)+6*xn[0]*xn[1]*xn[2]*1.5);
 kt[0][1]=coeff*(pow(xn[0],3.0)*6+pow(xn[1],3.0)*6+pow(xn[2],3.0)*1.5+3*pow(xn[0],2.0)*xn[1]*3+3*pow(xn[0],2.0)*xn[2]*1.5+3*xn[0]*pow(xn[1],2.0)*3+3*pow(xn[1],2.0)*xn[2]*1.5+3*xn[0]*pow(xn[2],2.0)+3*xn[1]*pow(xn[2],2.0)+6*xn[0]*xn[1]*xn[2]);
 kt[0][2]=coeff*(pow(xn[0],3.0)*6+pow(xn[1],3.0)*1.5+pow(xn[2],3.0)*6+3*pow(xn[0],2.0)*xn[1]*1.5+3*pow(xn[0],2.0)*xn[2]*3+3*xn[0]*pow(xn[1],2.0)+3*pow(xn[1],2.0)*xn[2]+3*xn[0]*pow(xn[2],2.0)*3+3*xn[1]*pow(xn[2],2.0)*1.5+6*xn[0]*xn[1]*xn[2]);
 kt[1][0]=coeff*(pow(xn[0],3.0)*6+pow(xn[1],3.0)*6+pow(xn[2],3.0)*1.5+3*pow(xn[0],2.0)*xn[1]*3+3*pow(xn[0],2.0)*xn[2]*1.5+3*xn[0]*pow(xn[1],2.0)*3+3*pow(xn[1],2.0)*xn[2]*1.5+3*xn[0]*pow(xn[2],2.0)+3*xn[1]*pow(xn[2],2.0)+6*xn[0]*xn[1]*xn[2]);
 kt[1][1]=coeff*(pow(xn[0],3.0)*3+pow(xn[1],3.0)*30+pow(xn[2],3.0)*3+3*pow(xn[0],2.0)*xn[1]*3+3*pow(xn[0],2.0)*xn[2]+3*xn[0]*pow(xn[1],2.0)*6+3*pow(xn[1],2.0)*xn[2]*6+3*xn[0]*pow(xn[2],2.0)+3*xn[1]*pow(xn[2],2.0)*3+6*xn[0]*xn[1]*xn[2]*1.5);
 kt[1][2]=coeff*(pow(xn[0],3.0)*1.5+pow(xn[1],3.0)*6+pow(xn[2],3.0)*6+3*pow(xn[0],2.0)*xn[1]+3*pow(xn[0],2.0)*xn[2]+3*xn[0]*pow(xn[1],2.0)*1.5+3*pow(xn[1],2.0)*xn[2]*3+3*xn[0]*pow(xn[2],2.0)*1.5+3*xn[1]*pow(xn[2],2.0)*3+6*xn[0]*xn[1]*xn[2]);
 kt[2][0]=coeff*(pow(xn[0],3.0)*6+pow(xn[1],3.0)*1.5+pow(xn[2],3.0)*6+3*pow(xn[0],2.0)*xn[1]*1.5+3*pow(xn[0],2.0)*xn[2]*3+3*xn[0]*pow(xn[1],2.0)+3*pow(xn[1],2.0)*xn[2]+3*xn[0]*pow(xn[2],2.0)*3+3*xn[1]*pow(xn[2],2.0)*1.5+6*xn[0]*xn[1]*xn[2]);
 kt[2][1]=coeff*(pow(xn[0],3.0)*1.5+pow(xn[1],3.0)*6+pow(xn[2],3.0)*6+3*pow(xn[0],2.0)*xn[1]+3*pow(xn[0],2.0)*xn[2]+3*xn[0]*pow(xn[1],2.0)*1.5+3*pow(xn[1],2.0)*xn[2]*3+3*xn[0]*pow(xn[2],2.0)*1.5+3*xn[1]*pow(xn[2],2.0)*3+6*xn[0]*xn[1]*xn[2]);
 kt[2][2]=coeff*(pow(xn[0],3.0)*3+pow(xn[1],3.0)*3+pow(xn[2],3.0)*30+3*pow(xn[0],2.0)*xn[1]+3*pow(xn[0],2.0)*xn[2]*3+3*xn[0]*pow(xn[1],2.0)+3*pow(xn[1],2.0)*xn[2]*3+3*xn[0]*pow(xn[2],2.0)*6+3*xn[1]*pow(xn[2],2.0)*6+6*xn[0]*xn[1]*xn[2]*1.5);
}

