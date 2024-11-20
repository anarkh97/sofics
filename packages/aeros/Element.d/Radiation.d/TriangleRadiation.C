#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Element.d/Radiation.d/TriangleRadiation.h>
#include <Corotational.d/TriangleThermalCorotator.h>
#include <Corotational.d/GeomState.h>

extern "C"      {
void   _FORTRAN(trianarea)(double*, double*, double*, double&);
}

TriangleRadiation::TriangleRadiation(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
        nn[2] = nodenums[2];
}

TriangleRadiation::~TriangleRadiation()
{
}

Element *
TriangleRadiation::clone()
{
 return new TriangleRadiation(*this);
}

void
TriangleRadiation::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

void
TriangleRadiation::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

double
TriangleRadiation::getMass(const CoordSet&) const
{
 return 0.0;
}

FullSquareMatrix
TriangleRadiation::massMatrix(const CoordSet &cs, double *d, int cmflg) const
{

        FullSquareMatrix mass(3,d);
        mass.zero();
        return mass;
}

FullSquareMatrix
TriangleRadiation::stiffness(const CoordSet &cs, double *Kcv, int flg) const
{

// ... Compute Radiative matrix

        FullSquareMatrix ret(3,Kcv);

        ret[0][0] = 0;
        ret[1][1] = 0;
        ret[2][2] = 0;
        ret[0][1] = 0;
        ret[0][2] = 0;
        ret[1][0] = 0;
        ret[1][2] = 0;
        ret[2][0] = 0;
        ret[2][1] = 0;

        return ret;
}

Corotator *
TriangleRadiation::getCorotator(CoordSet &cs, double* kel, int, int)
{
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


// ... Return thermal corotator 

	return new TriangleThermalCorotator(nn[0], nn[1], nn[2], area, prop->eps, prop->sigma, prop->Tr, cs);
}


int 
TriangleRadiation::numNodes() const
{
        return 3;
}

int*
TriangleRadiation::nodes(int *p) const
{
        if(p == 0) p = new int[3];
        p[0] = nn[0];
        p[1] = nn[1];
        p[2] = nn[2];
        return p;
}

int
TriangleRadiation::numDofs() const
{
        return 3;
}

int*
TriangleRadiation::dofs(DofSetArray &dsa, int *p) const
{
        if(p == 0) p = new int[3];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);
        p[2] = dsa.locate(nn[2],DofSet::Temp);

        return p;
}

void
TriangleRadiation::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 3, DofSet::Temp);
}

int
TriangleRadiation::getTopNumber() const
{
  return 149;
}

void 
TriangleRadiation::computeTemp(CoordSet&, State &, double[2], double*)
{
  fprintf(stderr," *** WARNING: TriangleRadiation::computeTemp is not implemented\n");
}

void 
TriangleRadiation::getFlFlux(double[2], double *, double *)
{
  fprintf(stderr," *** WARNING: TriangleRadiation::getFlFlux is not implemented\n");
}

