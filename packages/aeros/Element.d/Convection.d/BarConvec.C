#include <cmath>

#include <Element.d/Convection.d/BarConvec.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>


BarConvec::BarConvec(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
}

Element *
BarConvec::clone()
{
	return new BarConvec(*this);
}

void
BarConvec::renum(const int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
}

void
BarConvec::renum(EleRenumMap& table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
}

double
BarConvec::getMass(const CoordSet& cs) const
{
        return 0.0;
}

FullSquareMatrix
BarConvec::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{

        FullSquareMatrix elementMassMatrix(2,mel);

// zero the element mass matrix

	elementMassMatrix.zero();

        return elementMassMatrix;
}

FullSquareMatrix
BarConvec::stiffness(const CoordSet &cs, double *Kcv, int flg) const
{
// This is the additional matrix when convection is present.
// It is added into the conductance matrix.

        auto &nd1 = cs.getNode( nn[0] );
        auto &nd2 = cs.getNode( nn[1] );

        double x[2], y[2], z[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

	double dx = x[1] - x[0];
	double dy = y[1] - y[0];
	double dz = z[1] - z[0];

	double length = sqrt( dx*dx + dy*dy + dz*dz );

//... BOUNDARY convection for thermal bars ONLY ...
//... k and t are respectively the boundary convective coefficients of node 1 and 2 ...
//... Remember to add a new ATTRIBUTE for the convective elements ...

        double kcvb1 = prop->k*prop->A;
        double kcvb2 = prop->eh*prop->A;

//... LATERAL convection for the thermal bar ...
//... Also used as BOUNDARY convection for the thermal quad ...
//... If used for the quad, P is the depth (profondeur) ...

        double c = prop->c*prop->P*length/6;

//... Construct convective matrix ...

        FullSquareMatrix ret(2,Kcv);

        ret[0][0] = kcvb1 + 2*c;
        ret[1][1] = kcvb2 + 2*c;
        ret[1][0] =  c;
        ret[0][1] =  c;
                                    
        return ret;
}

int
BarConvec::numNodes() const
{
        return 2;
}

int *
BarConvec::nodes(int *p) const
{
        if(p == 0) p = new int[2];
        p[0] = nn[0];
        p[1] = nn[1];
        return p;
}

int
BarConvec::numDofs() const
{
        return 2;
}

int *
BarConvec::dofs(DofSetArray &dsa, int *p) const
{
        if(p == 0) p = new int[2];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);

        return p;
}

void
BarConvec::markDofs(DofSetArray& dsa) const
{
        dsa.mark( nn, 2, DofSet::Temp);
}

int
BarConvec::getTopNumber() const
{
  return 147;//1;
}
