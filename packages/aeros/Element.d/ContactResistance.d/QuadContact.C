#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <Element.d/ContactResistance.d/QuadContact.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>

// Interface with 2-node bars on each side

QuadContact::QuadContact(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
        nn[2] = nodenums[2];
        nn[3] = nodenums[3];
}

Element *
QuadContact::clone()
{
        return new QuadContact(*this);
}

void
QuadContact::renum(const int *table)
{
        nn[0] = table[nn[0]];
        nn[1] = table[nn[1]];
        nn[2] = table[nn[2]];
        nn[3] = table[nn[3]];
}

void
QuadContact::renum(EleRenumMap& table)
{
        nn[0] = table[nn[0]];
        nn[1] = table[nn[1]];
        nn[2] = table[nn[2]];
        nn[3] = table[nn[3]];
}

double
QuadContact::getMass(const CoordSet& cs) const
{
        return 0.0;
}

FullSquareMatrix
QuadContact::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{

        FullSquareMatrix elementMassMatrix(4,mel);

// zero the element mass matrix
        elementMassMatrix.zero();

        return elementMassMatrix;
}

FullSquareMatrix
QuadContact::stiffness(const CoordSet &cs, double *Kcv, int flg) const
{
// This is the additional matrix when thermal contact at an interface is present.
// It is added into the conductance matrix.

        auto &nd1 = cs.getNode( nn[0] );
        auto &nd2 = cs.getNode( nn[1] );
        auto &nd3 = cs.getNode( nn[2] );
        auto &nd4 = cs.getNode( nn[3] );

        double x[4], y[4], z[4];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

// Length of the interface

        double dx = x[1] - x[0];
        double dy = y[1] - y[0];
        double dz = z[1] - z[0];

        double length = sqrt( dx*dx + dy*dy + dz*dz );

// Construct contact matrix 

        FullSquareMatrix ret(4,Kcv);
        double coeff = (prop->c*length)/6;

        ret[0][0] =  2.0*coeff;
        ret[0][1] =  1.0*coeff;
        ret[0][2] = -1.0*coeff;
        ret[0][3] = -2.0*coeff;
        ret[1][0] =  1.0*coeff;
        ret[1][1] =  2.0*coeff;
        ret[1][2] = -2.0*coeff;
        ret[1][3] = -1.0*coeff;
        ret[2][0] = -1.0*coeff;
        ret[2][1] = -2.0*coeff;
        ret[2][2] =  2.0*coeff;
        ret[2][3] =  1.0*coeff;
        ret[3][0] = -2.0*coeff;
        ret[3][1] = -1.0*coeff;
        ret[3][2] =  1.0*coeff;
        ret[3][3] =  2.0*coeff;

        return ret;
}


int
QuadContact::numNodes() const
{
        return 4;
}

int*
QuadContact::nodes(int *p) const
{
        if(p == 0) p = new int[4];
        p[0] = nn[0];
        p[1] = nn[1];
        p[2] = nn[2];
        p[3] = nn[3];
        return p;
}

int
QuadContact::numDofs() const
{
        return 4;
}

int*
QuadContact::dofs(DofSetArray &dsa, int *p) const
{
        if(p == 0) p = new int[4];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);
        p[2] = dsa.locate(nn[2],DofSet::Temp);
        p[3] = dsa.locate(nn[3],DofSet::Temp);

        return p;
}

void
QuadContact::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 4, DofSet::Temp);
}

int
QuadContact::getTopNumber() const
{
  return 148;
}
