#include	<Element.d/ContactResistance.d/BrickContact.h>
#include        <Math.d/matrix.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>

extern "C"      {
void _FORTRAN(thermquad3a)(double*, double*, double*, double*, double*);

void _FORTRAN(qgauss)(int &, int &, int &, int &,
                      double &,  double &, double &);

void _FORTRAN(q4shpe)(double &, double &, double *, double *,
                      double *, double *, double *, double &);
}

// Interface with 4-node quads on each side

BrickContact::BrickContact(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
	nn[4] = nodenums[4];
	nn[5] = nodenums[5];
	nn[6] = nodenums[6];
	nn[7] = nodenums[7];
}

Element *
BrickContact::clone()
{
 return new BrickContact(*this);
}

void
BrickContact::renum(const int *table)
{
        nn[0] = table[nn[0]];
        nn[1] = table[nn[1]];
        nn[2] = table[nn[2]];
        nn[3] = table[nn[3]];
        nn[4] = table[nn[4]];
        nn[5] = table[nn[5]];
        nn[6] = table[nn[6]];
        nn[7] = table[nn[7]];
}

void
BrickContact::renum(EleRenumMap& table)
{
        nn[0] = table[nn[0]];
        nn[1] = table[nn[1]];
        nn[2] = table[nn[2]];
        nn[3] = table[nn[3]];
        nn[4] = table[nn[4]];
        nn[5] = table[nn[5]];
        nn[6] = table[nn[6]];
        nn[7] = table[nn[7]];
}

double
BrickContact::getMass(const CoordSet& cs) const
{

 return 0.0;

}

FullSquareMatrix
BrickContact::massMatrix(const CoordSet &cs,double *d,int cmflg) const
{
      FullSquareMatrix ret(8,d);
      ret.zero();

      return ret;
}


FullSquareMatrix
BrickContact::stiffness(const CoordSet &cs, double *Ks, int flg) const
{
	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);
	auto &nd3 = cs.getNode(nn[2]);
	auto &nd4 = cs.getNode(nn[3]);
	auto &nd5 = cs.getNode(nn[4]);
	auto &nd6 = cs.getNode(nn[5]);
	auto &nd7 = cs.getNode(nn[6]);
	auto &nd8 = cs.getNode(nn[7]);

	double x[8], y[8], z[8];

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z; 
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
	x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
	x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
	x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;
	x[6] = nd7.x; y[6] = nd7.y; z[6] = nd7.z;
	x[7] = nd8.x; y[7] = nd8.y; z[7] = nd8.z;
     
        double face1x[4], face1y[4], face1z[4];

        face1x[0] = nd1.x; face1y[0] = nd1.y; face1z[0] = nd1.z;   
        face1x[1] = nd2.x; face1y[1] = nd2.y; face1z[1] = nd2.z;   
        face1x[2] = nd3.x; face1y[2] = nd3.y; face1z[2] = nd3.z;   
        face1x[3] = nd4.x; face1y[3] = nd4.y; face1z[3] = nd4.z;   
  
        double xl[4], yl[4];

        // To get the local coordinates of face 1 of the iinterface

        _FORTRAN(thermquad3a)(face1x, face1y, face1z, xl, yl);

        int i,j;

        double h = prop->c;

        int numgauss = 2;

        FullSquareMatrix ret(8,Ks);
        ret.zero();

        int fortran = 1;  // fortran routines start from index 1
        int pt1, pt2;
        for (pt1 = 0 + fortran; pt1 < numgauss + fortran; pt1++)  {
           for (pt2 = 0 + fortran; pt2 < numgauss + fortran; pt2++)  {

              // get gauss point
              double xi, eta, wt;
              _FORTRAN(qgauss)(numgauss, pt1, numgauss, pt2, xi, eta, wt);

              //compute shape functions
              double shapeFunc[4], shapeGradX[4], shapeGradY[4];
              double detJ;  //det of jacobian

              _FORTRAN(q4shpe)(xi, eta, x, y,
                               shapeFunc, shapeGradX, shapeGradY, detJ);

              for (i = 0; i < 4; ++i) { 

                 for (j = 0; j < 4; ++j) 
                    ret[i][j] += h*shapeFunc[i]*shapeFunc[j]*detJ*wt;

                 for (j = 4; j < 8; ++j)
                    ret[i][j] += (-1)*h*shapeFunc[i]*shapeFunc[j-4]*detJ*wt;  

              }

              for (i = 4; i < 8; ++i) {

                 for (j = 0; j < 4; ++j)
                    ret[i][j] += (-1)*h*shapeFunc[i]*shapeFunc[j]*detJ*wt;

                 for (j = 4; j < 8; ++j)
                    ret[i][j] += h*shapeFunc[i]*shapeFunc[j-4]*detJ*wt;
              }
           }
        }

        return ret;
}

int
BrickContact::numNodes() const
{
 	return 8;
}

int*
BrickContact::nodes(int *p) const
{
 	if(p == 0) p = new int[8];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
 	p[4] = nn[4];
 	p[5] = nn[5];
 	p[6] = nn[6];
 	p[7] = nn[7];
	return p;
}

int
BrickContact::numDofs() const
{
 	return 8;
}

int*
BrickContact::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[8];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);
        p[2] = dsa.locate(nn[2],DofSet::Temp);
        p[3] = dsa.locate(nn[3],DofSet::Temp);
        p[4] = dsa.locate(nn[4],DofSet::Temp);
        p[5] = dsa.locate(nn[5],DofSet::Temp);
        p[6] = dsa.locate(nn[6],DofSet::Temp);
        p[7] = dsa.locate(nn[7],DofSet::Temp);

	return p;
}

void
BrickContact::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 8, DofSet::Temp);
}

int
BrickContact::getTopNumber() const
{
  return 151;
}

