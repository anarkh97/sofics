#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Element.d/BulkFluid.d/PentaBulk.h>
#include <Corotational.d/utilities.h>

extern "C"      {
void _FORTRAN(qgauss)(int &, int &, int &, int &,
                      double &,  double &, double &);

void _FORTRAN(q4shpe)(double &, double &, double *, double *,
                      double *, double *, double *, double &);
}


/* First node of the triangle = Bulk node (in the fluid)
   2nd,3rd,4th and 5th nodes = quad's nodes, on the cavity surface */

PentaBulk::PentaBulk(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
        nn[2] = nodenums[2];
        nn[3] = nodenums[3];
        nn[4] = nodenums[4];
}

Element *
PentaBulk::clone()
{
 return new PentaBulk(*this);
}

void
PentaBulk::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
}

void
PentaBulk::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
}

double
PentaBulk::getMass(const CoordSet &cs) const
{
        // The node 1 is at the top of the pyramid, nodes 2-3-4-5 form the base.

        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);
        auto &nd5 = cs.getNode(nn[4]);
        
        double x[5], y[5], z[5];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
        x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;

        // To compute the area of the base

        double T1[3], T2[3], T3[3], N[3], P[3];

        T1[0] = x[2]-x[1];
        T1[1] = y[2]-y[1];
        T1[2] = z[2]-z[1];

        double length = sqrt(pow(T1[0],2)+pow(T1[1],2)+pow(T1[2],2));

        T2[0] = x[4]-x[1];
        T2[1] = y[4]-y[1];
        T2[2] = z[4]-z[1];

        crossprod(T1,T2,N);
        normalize(N);

        crossprod(N,T1,P);
        normalize(P);

        double h = fabs(T2[0]*P[0]+T2[1]*P[1]+T2[2]*P[2]);

        double areabase = length*h;

        // To compute height of pyramid

        T3[0] = x[0]-x[1];
        T3[1] = y[0]-y[1];
        T3[2] = z[0]-z[1]; 

        double hpenta = fabs(T3[0]*N[0]+T3[1]*N[1]+T3[2]*N[2]);

        double volume = (1/3.0)*areabase*hpenta;
  
        double mass = prop->rho*volume;
 
        return mass;

}

FullSquareMatrix
PentaBulk::massMatrix(const CoordSet &cs, double *d, int cmflg) const
{

        FullSquareMatrix massMatrix(5,d);
        massMatrix.zero();

  // ... Only one element in the mass matrix is non zero

        double bulkMass = getMass(cs);
        massMatrix[0][0] = bulkMass*prop->Q; 

        return massMatrix;
}

FullSquareMatrix
PentaBulk::stiffness(const CoordSet &cs, double *Kcv, int flg) const
{

        double x[5], y[5], z[5];
        double xquad[4], yquad[4], zquad[4];
        double xl[4], yl[4];
        double T1[3], T2[3], T3[3], N[3], P[3], V[3];
        int i, j, k;

  // Get original coordinates of base's nodes

        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);
        auto &nd5 = cs.getNode(nn[4]);

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
        x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;

  // Coordinates of the base

        xquad[0] = nd2.x; yquad[0] = nd2.y; zquad[0] = nd2.z;
        xquad[1] = nd3.x; yquad[1] = nd3.y; zquad[1] = nd3.z;
        xquad[2] = nd4.x; yquad[2] = nd4.y; zquad[2] = nd4.z;
        xquad[3] = nd5.x; yquad[3] = nd5.y; zquad[3] = nd5.z;

  // To compute the area of the base

        T1[0] = x[2]-x[1];
        T1[1] = y[2]-y[1];
        T1[2] = z[2]-z[1];

        double length = sqrt(T1[0]*T1[0]+T1[1]*T1[1]+T1[2]*T1[2]);

        T2[0] = x[4]-x[1];
        T2[1] = y[4]-y[1];
        T2[2] = z[4]-z[1];

        crossprod(T1,T2,N);
        normalize(N);

        crossprod(N,T1,P);
        normalize(P);

        double h = fabs(T2[0]*P[0]+T2[1]*P[1]+T2[2]*P[2]);

        double areabase = length*h;
 
  // Redefine plane coordinates for the base, with origin at node 1 of the pyramid base
        double origin[3];
        origin[0] = xquad[0];
        origin[1] = yquad[0];
        origin[2] = zquad[0];


        //  Shift Origin to Node #1
        for (i = 0; i < 4; ++i) {
          xquad[i] -= origin[0];
          yquad[i] -= origin[1];
          zquad[i] -= origin[2];
        }
        // Local X-axis (Node 1->2)
        T1[0] = xquad[1];
        T1[1] = yquad[1];
        T1[2] = zquad[1];
        normalize( T1 );
        // Vector 2 from Node 4->2
        T2[0] = xquad[1] - xquad[3];
        T2[1] = yquad[1] - yquad[3];
        T2[2] = zquad[1] - zquad[3];
        normalize( T2 );
        // Vector 3 from Node 1->3
        T3[0] = xquad[2];
        T3[1] = yquad[2];
        T3[2] = zquad[2];
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
          xl[i] = (T1[0]*xquad[i]) + (T1[1]*yquad[i]) + (T1[2]*zquad[i]);
          yl[i] = (T2[0]*xquad[i]) + (T2[1]*yquad[i]) + (T2[2]*zquad[i]);
        }

  // ... Compute bulk fluid contribution matrix

          FullSquareMatrix ret(5,Kcv);
          ret.zero();

          double integrandStiffness[5][5];
          for (i = 0; i < 5; i++) {
              for (j = 0; j < 5; j++) {
                  integrandStiffness[i][j] = 0;
              }
          }
          integrandStiffness[0][0] = areabase;

          int fortran = 1;  // fortran routines start from index 1
          int pt1, pt2;
          int numgauss = 2;
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
                
                 // compute the integrand for the stiffness matrix

                     // First row and first column

                 for (k = 1; k < 5; k++) {
                     integrandStiffness[0][k] += -shapeFunc[k-1]*detJ*wt;
                     integrandStiffness[k][0] += -shapeFunc[k-1]*detJ*wt;
                 }

                 for (i = 1; i < 5; i++) {
                     for (k = 1; k < 5; k++) { 
                         integrandStiffness[i][k] += shapeFunc[i-1]*shapeFunc[k-1]*detJ*wt;
                     }
                 } 
             }
          }
          for (i = 0; i < 5; i++)
              for (j = 0; j < 5; j++)
                  ret[i][j] = prop->c*integrandStiffness[i][j];
          return ret;
}

int
PentaBulk::numNodes() const
{
        return 5;
}

int*
PentaBulk::nodes(int *p) const
{
        if(p == 0) p = new int[5];
        p[0] = nn[0];
        p[1] = nn[1];
        p[2] = nn[2];
        p[3] = nn[3];
        p[4] = nn[4];
        return p;
}

int
PentaBulk::numDofs() const
{
        return 5;
}

int*
PentaBulk::dofs(DofSetArray &dsa, int *p) const
{
        if(p == 0) p = new int[5];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);
        p[2] = dsa.locate(nn[2],DofSet::Temp);
        p[3] = dsa.locate(nn[3],DofSet::Temp);
        p[4] = dsa.locate(nn[4],DofSet::Temp);

        return p;
}

void
PentaBulk::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 5, DofSet::Temp);
}

int
PentaBulk::getTopNumber() const
{ 
  // Wrong top number - this type of element still needs to be defined in XPost
  return 151;
}
