#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Element.d/ContactResistance.d/PentaContact.h>

extern "C"      {
void   _FORTRAN(trianarea)(double*, double*, double*, double&);
}


// Interface with 3-node triangle on each side

PentaContact::PentaContact(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
        nn[2] = nodenums[2];
        nn[3] = nodenums[3];
        nn[4] = nodenums[4];
        nn[5] = nodenums[5];
}

Element *
PentaContact::clone()
{
 return new PentaContact(*this);
}

void
PentaContact::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
  nn[5] = table[nn[5]];
}

void
PentaContact::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
  nn[5] = table[nn[5]];
}

double
PentaContact::getMass(const CoordSet&) const
{
 return 0.0;
}

FullSquareMatrix
PentaContact::massMatrix(const CoordSet &cs, double *d, int cmflg) const
{

        FullSquareMatrix mass(6,d);
        mass.zero();
        return mass;
}

FullSquareMatrix
PentaContact::stiffness(const CoordSet &cs, double *Kcv, int flg) const
{
 
        int i;
        int j;

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

// ... Compute the interface contribution to the stiffness matrix

        FullSquareMatrix ret(6,Kcv);
        ret.zero();

        double coeff = prop->c*area/12;

        for( i = 0 ; i < 3 ; ++i ) {
           for( j = 0 ; j < 3 ; ++j ) {
              if ( i == j ) {
                 ret[i][j]     = coeff*2;
                 ret[i+3][j+3] = coeff*2;
                 ret[i][j+3]   = -coeff*2;
                 ret[i+3][j]   = -coeff*2;
              }
              else {
                 ret[i][j]     =  coeff;
                 ret[i+3][j+3] =  coeff;
                 ret[i][j+3]   = -coeff;
                 ret[i+3][j]   = -coeff;                  
              }
           }
        }
        return ret;
}


int
PentaContact::numNodes() const
{
        return 6;
}

int*
PentaContact::nodes(int *p) const
{
        if(p == 0) p = new int[6];
        p[0] = nn[0];
        p[1] = nn[1];
        p[2] = nn[2];
        p[3] = nn[3];
        p[4] = nn[4];
        p[5] = nn[5];
        return p;
}

int
PentaContact::numDofs() const
{
        return 6;
}

int*
PentaContact::dofs(DofSetArray &dsa, int *p) const
{
        if(p == 0) p = new int[6];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);
        p[2] = dsa.locate(nn[2],DofSet::Temp);
        p[3] = dsa.locate(nn[3],DofSet::Temp);
        p[4] = dsa.locate(nn[4],DofSet::Temp);
        p[5] = dsa.locate(nn[5],DofSet::Temp);

        return p;
}

void
PentaContact::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 6, DofSet::Temp);
}

int
PentaContact::getTopNumber() const
{
  return 124;
}

