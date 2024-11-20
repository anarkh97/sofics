#include <cmath>

#include        <Element.d/Convection.d/QuadConvec.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>

extern "C"      {
void   _FORTRAN(thermquad3a)(double*, double*, double*, double*, double*);
void   _FORTRAN(convecquad)(double*, double*,  double &, double*, const int&, 
                                        const int&);
}


// Four Node Quad element (Galerkin)

QuadConvec::QuadConvec(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
}

Element *
QuadConvec::clone()
{
 return new QuadConvec(*this);
}

void
QuadConvec::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

void
QuadConvec::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

double
QuadConvec::getMass(const CoordSet&) const
{
 return 0.0;
}

FullSquareMatrix
QuadConvec::massMatrix(const CoordSet &cs, double *d, int cmflg) const
{

        FullSquareMatrix mass(4,d);
        mass.zero();
        return mass;
}

FullSquareMatrix
QuadConvec::stiffness(const CoordSet &cs, double *Kcv, int flg) const
{

	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);
	auto &nd3 = cs.getNode(nn[2]);
	auto &nd4 = cs.getNode(nn[3]);

        int i;
	double x[4], y[4], z[4], Ke[16];
	double xl[4], yl[4];

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
	x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

	const int numgauss = 2;
	const int numdof   = 4;

// ... Get Local coordinates ...

        _FORTRAN(thermquad3a)(x, y, z, xl, yl);

        double c = prop ->c;

        _FORTRAN(convecquad)(xl, yl,c, Ke, numgauss, numdof);

        for (i=0;i<16;i++)
         Kcv[i] = Ke[i];

        FullSquareMatrix ret(4, Kcv);

        return ret;
}

int
QuadConvec::numNodes() const
{
 	return 4;
}

int*
QuadConvec::nodes(int *p) const
{
 	if(p == 0) p = new int[4];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
	return p;
}

int
QuadConvec::numDofs() const
{
 	return 4;
}

int*
QuadConvec::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[4];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);
        p[2] = dsa.locate(nn[2],DofSet::Temp);
        p[3] = dsa.locate(nn[3],DofSet::Temp);

	return p;
}

void
QuadConvec::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 4, DofSet::Temp);
}

int
QuadConvec::getTopNumber() const
{
  return 188;//148;//2;
}
