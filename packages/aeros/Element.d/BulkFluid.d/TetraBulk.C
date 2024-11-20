#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Element.d/BulkFluid.d/TetraBulk.h>
#include <Corotational.d/utilities.h>

extern "C"      {
void    _FORTRAN(mass23)(double*, double*, double*, double&, double*,
                         const int&, double*, double*, const int&,
                         double&, const int&);
void   _FORTRAN(trianarea)(double*, double*, double*, double&);
}


/* First node of the triangle = Bulk node (in the fluid)
   2nd,3rd and 4th nodes = triangle's nodes, on the cavity surface*/

TetraBulk::TetraBulk(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
        nn[2] = nodenums[2];
        nn[3] = nodenums[3];
}

Element *
TetraBulk::clone()
{
 return new TetraBulk(*this);
}

void
TetraBulk::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

void
TetraBulk::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

double
TetraBulk::getMass(const CoordSet &cs) const
{
        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);

        double x[4], y[4], z[4];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

        double T1[3], T2[3], T3[3], V[3];
        
        T1[0] = x[0]-x[3];
        T1[1] = y[0]-y[3];
        T1[2] = z[0]-z[3];

        T2[0] = x[1]-x[3];
        T2[1] = y[1]-y[3];
        T2[2] = z[1]-z[3];

        T3[0] = x[2]-x[3];
        T3[1] = y[2]-y[3];
        T3[2] = z[2]-z[3];

        crossprod(T2,T3,V);

        double volume = (1/6.0)*fabs(T1[0]*V[0]+T1[1]*V[1]+T1[2]*V[2]);
        double mass = prop->rho*volume;
 
        return mass;

}

FullSquareMatrix
TetraBulk::massMatrix(const CoordSet &cs, double *d, int cmflg) const
{

        FullSquareMatrix mass(3,d);
        mass.zero();

  // ... Only one element in the mass matrix is non zero

        double bulkMass = getMass(cs);
        mass[0][0] = bulkMass*prop->Q; 

        return mass;
}

FullSquareMatrix
TetraBulk::stiffness(const CoordSet &cs, double *Kcv, int flg) const
{

  // ... Need the area of the tetrahedral's triangular base, which is on the surface of the cavity
  // ... Base = node 2, 3, 4

        auto &nd1 = cs.getNode(nn[1]);
        auto &nd2 = cs.getNode(nn[2]);
        auto &nd3 = cs.getNode(nn[3]);

        double x[3], y[3], z[3];
        double area;

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  // ... Compute area of triangle

        _FORTRAN(trianarea)(x,y,z,area);

// ... Compute bulk fluid contribution matrix

          FullSquareMatrix ret(4,Kcv);

          double coeff = prop->c*area/12;
          ret[0][0] = 12;
          ret[0][1] = -4;
          ret[0][2] = -4;
          ret[0][3] = -4;
          ret[1][0] = -4; 
          ret[1][1] =  2;
          ret[1][2] =  1;
          ret[1][3] =  1;
          ret[2][0] = -4;
          ret[2][1] =  1;
          ret[2][2] =  2;
          ret[2][3] =  1;
          ret[3][0]=  -4;
          ret[3][1] =  1;
          ret[3][2] =  1;
          ret[3][3] =  2;

          ret *= coeff;

        return ret;
}

int
TetraBulk::numNodes() const
{
        return 4;
}

int*
TetraBulk::nodes(int *p) const
{
        if(p == 0) p = new int[4];
        p[0] = nn[0];
        p[1] = nn[1];
        p[2] = nn[2];
        p[3] = nn[3];
        return p;
}

int
TetraBulk::numDofs() const
{
        return 4;
}

int*
TetraBulk::dofs(DofSetArray &dsa, int *p) const
{
        if(p == 0) p = new int[4];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);
        p[2] = dsa.locate(nn[2],DofSet::Temp);
        p[3] = dsa.locate(nn[3],DofSet::Temp);

        return p;
}

void
TetraBulk::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 4, DofSet::Temp);
}

int
TetraBulk::getTopNumber() const
{
  return 150;
}
