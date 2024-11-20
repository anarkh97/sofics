#include        <Element.d/Convection.d/TriangleConvec.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>

extern "C"      {
void   _FORTRAN(trianarea)(double*, double*, double*, double&);
}

// Tree node triangle

TriangleConvec::TriangleConvec(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
}

Element *
TriangleConvec::clone()
{
 return new TriangleConvec(*this);
}

void
TriangleConvec::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

void
TriangleConvec::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

double
TriangleConvec::getMass(const CoordSet&) const
{
 return 0.0;
}

FullSquareMatrix
TriangleConvec::massMatrix(const CoordSet &cs, double *d, int cmflg) const
{

        FullSquareMatrix mass(3,d);
        mass.zero();
        return mass;
}

FullSquareMatrix
TriangleConvec::stiffness(const CoordSet &cs, double *Kcv, int flg) const
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
        double c = prop ->c;

// ... Compute Convective matrix

          FullSquareMatrix ret(3,Kcv);

          ret[0][0] = c*area/6;
          ret[1][1] = c*area/6;
          ret[2][2] = c*area/6;
          ret[0][1] = c*area/12;
          ret[0][2] = c*area/12;
          ret[1][0] = c*area/12;
          ret[1][2] = c*area/12;
          ret[2][0] = c*area/12;
          ret[2][1] = c*area/12;

        return ret;
}

int
TriangleConvec::numNodes() const
{
 	return 3;
}

int*
TriangleConvec::nodes(int *p) const
{
 	if(p == 0) p = new int[3];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
	return p;
}

int
TriangleConvec::numDofs() const
{
 	return 3;
}

int*
TriangleConvec::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[3];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);
        p[2] = dsa.locate(nn[2],DofSet::Temp);

	return p;
}

void
TriangleConvec::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 3, DofSet::Temp);
}

int
TriangleConvec::getTopNumber() const
{
  return 149;//4;
}
