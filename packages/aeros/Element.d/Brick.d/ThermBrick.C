#include	<Element.d/Brick.d/ThermBrick.h>
#include        <Math.d/matrix.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>

extern "C"      {
void  _FORTRAN(thermbrik8v)(double*, double*, double*, const int&,
                            double*, const int&, const int& );
void  _FORTRAN(thermbr8mas)(const int&,double&,double*,const int&,
                            double*,double*, double*, double &);
}

ThermBrick::ThermBrick(int* nodenums)
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
ThermBrick::clone()
{
 return new ThermBrick(*this);
}

void
ThermBrick::renum(const int *table)
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
ThermBrick::renum(EleRenumMap& table)
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
ThermBrick::getMass(const CoordSet& cs) const
{

 return 0.0;

}

FullSquareMatrix
ThermBrick::massMatrix(const CoordSet &cs,double *d,int cmflg) const
{
     auto &nd1 = cs.getNode(nn[0]);
     auto &nd2 = cs.getNode(nn[1]);
     auto &nd3 = cs.getNode(nn[2]);
     auto &nd4 = cs.getNode(nn[3]);
     auto &nd5 = cs.getNode(nn[4]);
     auto &nd6 = cs.getNode(nn[5]);
     auto &nd7 = cs.getNode(nn[6]);
     auto &nd8 = cs.getNode(nn[7]);

     double x[8], y[8], z[8], mm[64];

     x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
     x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
     x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
     x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
     x[4] = nd5.x; y[4] = nd5.y; z[4] = nd5.z;
     x[5] = nd6.x; y[5] = nd6.y; z[5] = nd6.z;
     x[6] = nd7.x; y[6] = nd7.y; z[6] = nd7.z;
     x[7] = nd8.x; y[7] = nd8.y; z[7] = nd8.z;

     int i;

     double totmas = 0.0;

     const int numgauss = 2;
     const int numdof   = 8;

/* Lumped mass matrix */

     _FORTRAN(thermbr8mas)(numgauss,prop->rho,mm,numdof,x,y,z,totmas);

      for (i=0;i<64;i++){
       d[i] = prop->Q*mm[i];
      } 

      FullSquareMatrix ret(8,d);
      return ret;
}


FullSquareMatrix
ThermBrick::stiffness(const CoordSet &cs, double *Ks, int flg) const
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

  	const int numgauss = 2;
  	const int numdof   = 8;
	const int outerr   = 6;

        _FORTRAN(thermbrik8v)(x, y, z, numgauss, Ks, numdof, outerr);

        FullSquareMatrix ret(8,Ks);

        double k = prop ->k;

        int i,j;
         for(i=0;i<8;i++)
          for(j=0;j<8;j++) {
           ret[i][j] = k*ret[i][j];
         }
        return ret;
}

int
ThermBrick::numNodes() const
{
 	return 8;
}

int*
ThermBrick::nodes(int *p) const
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
ThermBrick::numDofs() const
{
 	return 8;
}

int*
ThermBrick::dofs(DofSetArray &dsa, int *p) const
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
ThermBrick::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 8, DofSet::Temp);
}

int
ThermBrick::getTopNumber() const
{
  return 117;//151;//3;
}


Corotator *
ThermBrick::getCorotator(CoordSet &cs, double* kel, int, int)
{
  return 0;
}

