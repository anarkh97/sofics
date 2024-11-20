#include	<Element.d/FluidTriangle3.d/SloshTriangleFS.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Math.d/Vector.h>
#include        <Utils.d/dofset.h>
#include        <cmath>
#include        <Element.d/State.h>

SloshTriangleFS::SloshTriangleFS(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
}


Element *
SloshTriangleFS::clone()
{
 return new SloshTriangleFS(*this);
}


void
SloshTriangleFS::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

void
SloshTriangleFS::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}


double
SloshTriangleFS::getArea(const CoordSet& cs) const
{
  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  Vector r1(3), r2(3), r3(3);

  r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
  r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
  r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;

  Vector v1(3), v2(3), v3(3), v4(3), v5(3);

  v1 = r2 - r1;
  v2 = r3 - r1;

  v3 = v1.cross(v2);

  double area = 0.5*v3.magnitude();

  return area;
}

FullSquareMatrix
SloshTriangleFS::massMatrix(const CoordSet &cs,double *mel,int cmflg) const
{
	double area = getArea(cs);

        FullSquareMatrix ret(3,mel);

	ret.zero();

//Calculate entries
	int i;
	int j;

// This is the LUMPED mass 

        //for(i=0; i<3; ++i)
        //  ret[i][i] = area/3;

// This is the CONSISTENT mass

        for(i=0; i<3; ++i) {
           for(j=0; j<3; ++j)
              ret[i][j] = area/12;
        }

        for(i=0; i<3; ++i)
           ret[i][i] = area/6;
        
        return ret;
}

FullSquareMatrix
SloshTriangleFS::stiffness(const CoordSet &cs, double *d, int flg) const
{
        FullSquareMatrix K(3,d);

        K.zero();

        return K;
}

int
SloshTriangleFS::numNodes() const
{
 	return 3;
}

int*
SloshTriangleFS::nodes(int *p) const
{
 	if(p == 0) p = new int[3];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
	return p;
}

int
SloshTriangleFS::numDofs() const
{
 	return 3;
}

int*
SloshTriangleFS::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[3];

        p[0] = dsa.locate(nn[0],DofSet::Potential);
        p[1] = dsa.locate(nn[1],DofSet::Potential);
        p[2] = dsa.locate(nn[2],DofSet::Potential);

	return p;
}

void
SloshTriangleFS::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 3, DofSet::Potential);
}

int
SloshTriangleFS::getTopNumber() const
{
  return 153;//4;
}

