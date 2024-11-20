#include <cmath>

#include <Math.d/FullSquareMatrix.h>
#include <Element.d/Element.h>
#include <Corotational.d/QuadThermalCorotator.h>
#include <Corotational.d/TemperatureState.h> 
#include <Corotational.d/utilities.h> 

extern "C"      {
void _FORTRAN(qgauss)(int &, int &, int &, int &,
                      double &,  double &, double &);

void _FORTRAN(q4shpe)(double &, double &, double *, double *,
                      double *, double *, double *, double &);
}

QuadThermalCorotator::QuadThermalCorotator(int _n1, int _n2, int _n3, int _n4,
                                           double _eps, double _sigma, double _Tr, CoordSet& cs)
{
 n1    = _n1;         	// Node 1
 n2    = _n2;   	// Node 2
 n3    = _n3;           // Node 3
 n4    = _n4;           // Node 4 
 eps   = _eps;          // Emissivity of the body
 sigma = _sigma;        // Stefan's constant
 Tr    = _Tr;           // Temperature of the enclosure receiving the radiation
}

void
QuadThermalCorotator::getStiffAndForce(GeomState &ts, CoordSet &cs, 
                                       FullSquareMatrix &elK, double *f, double dt, double t)
/*******************************************************************
 *
 * Purpose :
 *  Compute Tangent Stiffness matrix and internal force 
 *  for quad thermal element in current configuration.
 *
 * Input Variables:
 *  ts        : current temperature 
 *  cs        : coordinate set, contains reference configuration
 *
 * Local Variables:
 * xn[i] : temperature in the current configuration
 *
 * Output :
 * elK      : element tangent stiffness matrix
 * f        : element internal force
 *
 *****************************************************************/
{
 // Declare local variables 
 int    i, j;
 double kt[4][4], ff[4], xn[4];
 double x[4], y[4], z[4];
 double xl[4], yl[4];
 double T1[3], T2[3], T3[3], V[3];

 // Get original coordinates of quad's nodes
 auto &nd1 = cs.getNode(n1);
 auto &nd2 = cs.getNode(n2);
 auto &nd3 = cs.getNode(n3);
 auto &nd4 = cs.getNode(n4);
    
 x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
 x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
 x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
 x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
    
 // Redefine plane coordinates for the quad, with origin at node 1
 double origin[3];
 origin[0] = x[0];
 origin[1] = y[0];
 origin[2] = z[0];

 //  Shift Origin to Node #1
 for (i = 0; i < 4; ++i) {
   x[i] -= origin[0];
   y[i] -= origin[1];
   z[i] -= origin[2];
 }
 // Local X-axis (Node 1->2)
 T1[0] = x[1];
 T1[1] = y[1];
 T1[2] = z[1];
 normalize( T1 );
 // Vector 2 from Node 4->2
 T2[0] = x[1] - x[3];
 T2[1] = y[1] - y[3];
 T2[2] = z[1] - z[3];
 normalize( T2 );
 // Vector 3 from Node 1->3
 T3[0] = x[2];
 T3[1] = y[2];
 T3[2] = z[2];
 normalize( T3 );
 // Perpendicular as cross between v2 and v3
 crossprod( T2, T3, V );
 normalize( V );
 // Local Y-axis as cross between X and V
 crossprod( V, T1, T2 );
 normalize( T2);
 // Local Z-axis as cross between X and Y
 crossprod( T1, T2, T3 );
 normalize( T3);

 // Compute Local "In-plane" Coordinates (required by shape routines)
 for (i = 0; i < 4; ++i) {
   xl[i] = 0.0;
   yl[i] = 0.0;
 }
 for (i = 0; i < 4; ++i) {
   xl[i] = (T1[0]*x[i]) + (T1[1]*y[i]) + (T1[2]*z[i]);
   yl[i] = (T2[0]*x[i]) + (T2[1]*y[i]) + (T2[2]*z[i]);
 }

 // Get current temperature state
 NodeState &tn1 = ts[n1];
 NodeState &tn2 = ts[n2];
 NodeState &tn3 = ts[n3];
 NodeState &tn4 = ts[n4];

 // Set temperature of Cn configuration 
 xn[0] = tn1.x; // temperature of node 1
 xn[1] = tn2.x; // temperature of node 2
 xn[2] = tn3.x; // temperature of node 3
 xn[3] = tn4.x; // temperature of node 4

 // Form tangent stiffness matrix
 formTangentStiffness(xl, yl, xn, eps, sigma, kt); 

 // Copy tangent stiffness matrix to element K matrix
 for(i=0; i<4; ++i)
   for(j=0; j<4; ++j)
     elK[i][j] = kt[i][j];
 
 // Form internal force
 formInternalForce(xl, yl, xn, eps, sigma, Tr, ff);

 // Copy internal force to element f matrix
 for(j=0; j<4; ++j)
   f[j]=ff[j];
}

void
QuadThermalCorotator::getInternalForce(GeomState &ts, CoordSet &cs, 
                                       FullSquareMatrix &elK, double *f, double dt, double t)
/*******************************************************************
 *
 * Purpose :
 *  Compute internal force 
 *  for quad thermal element in current configuration.
 *
 * Input Variables:
 *  ts        : current temperature 
 *  cs        : coordinate set, contains reference configuration
 *
 * Local Variables:
 * xn[i] : temperature in the current configuration
 *
 * Output :
 * f        : element internal force
 *
 *****************************************************************/
{
 // Declare local variables 
 int    i, j;
 double ff[4], xn[4];
 double x[4], y[4], z[4];
 double xl[4], yl[4];
 double T1[3], T2[3], T3[3], V[3];

 // Get original coordinates of quad's nodes
 auto &nd1 = cs.getNode(n1);
 auto &nd2 = cs.getNode(n2);
 auto &nd3 = cs.getNode(n3);
 auto &nd4 = cs.getNode(n4);
    
 x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
 x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
 x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
 x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
    
 // Redefine plane coordinates for the quad, with origin at node 1
 double origin[3];
 origin[0] = x[0];
 origin[1] = y[0];
 origin[2] = z[0];

 //  Shift Origin to Node #1
 for (i = 0; i < 4; ++i) {
   x[i] -= origin[0];
   y[i] -= origin[1];
   z[i] -= origin[2];
 }
 // Local X-axis (Node 1->2)
 T1[0] = x[1];
 T1[1] = y[1];
 T1[2] = z[1];
 normalize( T1 );
 // Vector 2 from Node 4->2
 T2[0] = x[1] - x[3];
 T2[1] = y[1] - y[3];
 T2[2] = z[1] - z[3];
 normalize( T2 );
 // Vector 3 from Node 1->3
 T3[0] = x[2];
 T3[1] = y[2];
 T3[2] = z[2];
 normalize( T3 );
 // Perpendicular as cross between v2 and v3
 crossprod( T2, T3, V );
 normalize( V );
 // Local Y-axis as cross between X and V
 crossprod( V, T1, T2 );
 normalize( T2);
 // Local Z-axis as cross between X and Y
 crossprod( T1, T2, T3 );
 normalize( T3);

 // Compute Local "In-plane" Coordinates (required by shape routines)
 for (i = 0; i < 4; ++i) {
   xl[i] = 0.0;
   yl[i] = 0.0;
 }
 for (i = 0; i < 4; ++i) {
   xl[i] = (T1[0]*x[i]) + (T1[1]*y[i]) + (T1[2]*z[i]);
   yl[i] = (T2[0]*x[i]) + (T2[1]*y[i]) + (T2[2]*z[i]);
 }

 // Get current temperature state
 NodeState &tn1 = ts[n1];
 NodeState &tn2 = ts[n2];
 NodeState &tn3 = ts[n3];
 NodeState &tn4 = ts[n4];

 // Set temperature of Cn configuration 
 xn[0] = tn1.x; // temperature of node 1
 xn[1] = tn2.x; // temperature of node 2
 xn[2] = tn3.x; // temperature of node 3
 xn[3] = tn4.x; // temperature of node 4

 // Form internal force
 formInternalForce(xl, yl, xn, eps, sigma, Tr, ff);

 // Copy internal force to element f matrix
 for(j=0; j<4; ++j)
   f[j]=ff[j];
}

void
QuadThermalCorotator::formInternalForce(double xl[4], double yl[4], double xn[4], double eps, double sigma, double Tr, double f[4])
/*******************************************************************
 *
 * Purpose :
 *  Compute internal force vector for thermal quad
 *  element in current configuration.
 * Input  :
 *  xl, yl   : plane coordinates for quad's nodes
 *  xn       : temperature of quad's nodes
 * Output :
 *  f        : internal force vector in current configuration
 *
 *****************************************************************/
{
  f[0]=0;
  f[1]=0;
  f[2]=0;
  f[3]=0;
  
  int numgauss = 3;
  int i;
 
  // Compute internal force in local system and store in f

  int fortran = 1;  // fortran routines start from index 1
  int pt1, pt2;
  for (pt1 = 0 + fortran; pt1 < numgauss + fortran; pt1++)  {
      for (pt2 = 0 + fortran; pt2 < numgauss + fortran; pt2++)  {

        // get gauss point
        double xi, eta, wt;
        _FORTRAN(qgauss)(numgauss, pt1, numgauss, pt2, xi, eta, wt);

        // compute shape functions
        double shapeFunc[4], shapeGradX[4], shapeGradY[4];
        double detJ;  //det of jacobian

        _FORTRAN(q4shpe)(xi, eta, xl, yl,
                         shapeFunc, shapeGradX, shapeGradY, detJ);
        
        // compute the integrand for the internal force
        double integrandInternal[4];
        for (i = 0; i < 4; ++i)
          integrandInternal[i] = 0;

        double T4power = pow(xn[0]*shapeFunc[0]+xn[1]*shapeFunc[1]+xn[2]*shapeFunc[2]+xn[3]*shapeFunc[3],4.0);

        for (i = 0; i < 4; ++i) 
          integrandInternal[i] = eps*sigma*shapeFunc[i]*(T4power-pow(Tr,4.0));

        // get internal force
        for (i = 0; i < 4; ++i) 
          f[i] += wt*integrandInternal[i]*detJ;
        
      }
  }

}


void
QuadThermalCorotator::formTangentStiffness(double xl[4], double yl[4], double xn[4], double eps, double sigma, double kt[4][4]) 
/*******************************************************************
 * 
 * Purpose :
 *  Compute tangential stiffness for triangle quad element
 *  in current configuration.
 *
 * Input :
 *  xn       : current nodal temperature vector
 *  xl       : local x-coordinates
 *  yl       : local y-coordinates
 *  
 * Output :
 *  kt       : tangent stiffness matrix in current configuration
 *
 *****************************************************************/
{
     int i, j;
     int numgauss = 3;

     for(i=0; i<4; ++i )
       for(j=0; j<4; ++j )
         kt[i][j] = 0.0;

  // Tangent stiffness matrix 
    int fortran = 1;  // fortran routines start from index 1
    int pt1, pt2;
    for (pt1 = 0 + fortran; pt1 < numgauss + fortran; pt1++)  {
      for (pt2 = 0 + fortran; pt2 < numgauss + fortran; pt2++)  {
        // get gauss point
        double xi, eta, wt;
        _FORTRAN(qgauss)(numgauss, pt1, numgauss, pt2, xi, eta, wt);

        // compute shape functions
        double shapeFunc[4], shapeGradX[4], shapeGradY[4];
        double detJ;  //det of jacobian

        _FORTRAN(q4shpe)(xi, eta, xl, yl,
                         shapeFunc, shapeGradX, shapeGradY, detJ);

        // compute integrand for the stiffness matrix
        double integrandStiff[4][4];
        for (i = 0; i < 4; ++i)
          for (j = 0; j < 4; ++j)
            integrandStiff[i][j] = 0;

        double T3power = pow(xn[0]*shapeFunc[0]+xn[1]*shapeFunc[1]+xn[2]*shapeFunc[2]+xn[3]*shapeFunc[3],3.0);
     
        for (i = 0; i < 4; ++i) 
          for (j = 0; j < 4; ++j) 
            integrandStiff[i][j] = 4*eps*sigma*shapeFunc[i]*shapeFunc[j]*T3power;         
          
        // get the tangent stiffness matrix
        for (i = 0; i < 4; ++i) 
          for (j = 0; j < 4; ++j) 
            kt[i][j] += wt*integrandStiff[i][j]*detJ;
      }
    }

}

