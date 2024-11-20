#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Element.d/BulkFluid.d/TriangleBulk.h>

extern "C"      {
void   _FORTRAN(trianarea)(double*, double*, double*, double&);
}

/* First node of the triangle = Bulk node (in the fluid)
   2nd and 3rd nodes = bar's nodes, on the cavity surface */

TriangleBulk::TriangleBulk(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
        nn[2] = nodenums[2];
}

Element *
TriangleBulk::clone()
{
 return new TriangleBulk(*this);
}

void
TriangleBulk::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

void
TriangleBulk::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

double
TriangleBulk::getMass(const CoordSet &cs) const
{
        double mass = 0;

        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3];
        double area;

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  // ... Compute area of triangle

        _FORTRAN(trianarea)(x,y,z,area);

        mass = prop->rho*prop->eh*area;

        return mass;
}

FullSquareMatrix
TriangleBulk::massMatrix(const CoordSet &cs, double *d, int cmflg) const
{

        FullSquareMatrix mass(3,d);
        mass.zero();

  // Only one element in the mass matrix is non zero

        double bulkMass = getMass(cs);
        mass[0][0] = bulkMass*prop->Q;
        return mass;
}

FullSquareMatrix
TriangleBulk::stiffness(const CoordSet &cs, double *Kcv, int flg) const
{

// ... Compute length of the interface

        auto &nd2 = cs.getNode( nn[1] );
        auto &nd3 = cs.getNode( nn[2] );

        double xi[2], yi[2], zi[2];

        xi[0] = nd2.x; yi[0] = nd2.y; zi[0] = nd2.z;
        xi[1] = nd3.x; yi[1] = nd3.y; zi[1] = nd3.z;

        //double dx = xi[1] - xi[0];
        //double dy = yi[1] - yi[0];
        //double dz = zi[1] - zi[0];

        //double length = sqrt( dx*dx + dy*dy + dz*dz );

// ... Compute bulk fluid contribution matrix

        FullSquareMatrix ret(3,Kcv);

        //double coeff = prop->c*length/6;
        ret[0][0] =  6;
        ret[0][1] = -3;
        ret[0][2] = -3;
        ret[1][0] = -3;
        ret[1][1] =  2;
        ret[1][2] =  1;
        ret[2][0] = -3;
        ret[2][1] =  1;
        ret[2][2] =  2;

        return ret;
}

int
TriangleBulk::numNodes() const
{
        return 3;
}

int*
TriangleBulk::nodes(int *p) const
{
        if(p == 0) p = new int[3];
        p[0] = nn[0];
        p[1] = nn[1];
        p[2] = nn[2];
        return p;
}

int
TriangleBulk::numDofs() const
{
        return 3;
}

int*
TriangleBulk::dofs(DofSetArray &dsa, int *p) const
{
        if(p == 0) p = new int[3];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);
        p[2] = dsa.locate(nn[2],DofSet::Temp);

        return p;
}

void
TriangleBulk::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 3, DofSet::Temp);
}

int
TriangleBulk::getTopNumber() const
{
  return 149;
}
