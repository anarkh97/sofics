#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Element.d/Radiation.d/QuadRadiation.h>
#include <Corotational.d/QuadThermalCorotator.h>
#include <Corotational.d/GeomState.h>

// Four node quadrilateral

QuadRadiation::QuadRadiation(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
        nn[2] = nodenums[2];
        nn[3] = nodenums[3];
}

QuadRadiation::~QuadRadiation()
{
}

Element *
QuadRadiation::clone()
{
 return new QuadRadiation(*this);
}

void
QuadRadiation::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

void
QuadRadiation::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

double
QuadRadiation::getMass(const CoordSet&) const
{
 return 0.0;
}

FullSquareMatrix
QuadRadiation::massMatrix(const CoordSet &cs, double *d, int cmflg) const
{

        FullSquareMatrix mass(4,d);
        mass.zero();
        return mass;
}

FullSquareMatrix
QuadRadiation::stiffness(const CoordSet &cs, double *Kcv, int flg) const
{

// ... Compute Radiative matrix

        FullSquareMatrix ret(4,Kcv);

        ret[0][0] = 0;
        ret[1][1] = 0;
        ret[2][2] = 0;
        ret[3][3] = 0;
        ret[0][1] = 0;
        ret[0][2] = 0;
        ret[0][3] = 0;
        ret[1][0] = 0;
        ret[1][2] = 0;
        ret[1][3] = 0;
        ret[2][0] = 0;
        ret[2][1] = 0;
        ret[2][2] = 0;
        ret[2][3] = 0;
        ret[3][0] = 0;
        ret[3][1] = 0;
        ret[3][2] = 0;
        ret[3][3] = 0;

        return ret;
}

Corotator *
QuadRadiation::getCorotator(CoordSet &cs, double* kel, int, int)
{

// ... Return thermal corotator 
	
	return new QuadThermalCorotator(nn[0], nn[1], nn[2], nn[3], prop->eps, prop->sigma, prop->Tr, cs);
}


int 
QuadRadiation::numNodes() const
{
        return 4;
}

int*
QuadRadiation::nodes(int *p) const
{
        if(p == 0) p = new int[4];
        p[0] = nn[0];
        p[1] = nn[1];
        p[2] = nn[2];
        p[3] = nn[3];
        return p;
}

int
QuadRadiation::numDofs() const
{
        return 4;
}

int*
QuadRadiation::dofs(DofSetArray &dsa, int *p) const
{
        if(p == 0) p = new int[4];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);
        p[2] = dsa.locate(nn[2],DofSet::Temp);
        p[3] = dsa.locate(nn[3],DofSet::Temp);

        return p;
}

void
QuadRadiation::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 4, DofSet::Temp);
}

int
QuadRadiation::getTopNumber() const
{
  return 148;
}

void 
QuadRadiation::computeTemp(CoordSet&, State &, double[2], double*)
{
  fprintf(stderr," *** WARNING: QuadRadiation::computeTemp is not implemented\n");
}

void 
QuadRadiation::getFlFlux(double[2], double *, double *)
{
  fprintf(stderr," *** WARNING: QuadRadiation::getFlFlux is not implemented\n");
}

