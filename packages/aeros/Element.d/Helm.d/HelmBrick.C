#include	<Element.d/Helm.d/HelmBrick.h>
#include        <Math.d/matrix.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>

extern "C"      {
void  _FORTRAN(helmbrik8v)(double*, double*, double*, const int&,
                       double*, const int&, const int& );
void  _FORTRAN(helmbr8mas)(const int&,double*,const int&,double*,double*, double*, double &);
}

HelmBrick::HelmBrick(int* nodenums)
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
HelmBrick::clone()
{
 return new HelmBrick(*this);
}

void
HelmBrick::renum(const int *table)
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
HelmBrick::renum(EleRenumMap& table)
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
HelmBrick::getMass(const CoordSet& cs) const
{

 return 0.0;

}

FullSquareMatrix
HelmBrick::massMatrix(const CoordSet &cs,double *d,int cmflg) const
{
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
        Node nd4 = cs.getNode(nn[3]);
        Node nd5 = cs.getNode(nn[4]);
        Node nd6 = cs.getNode(nn[5]);
        Node nd7 = cs.getNode(nn[6]);
        Node nd8 = cs.getNode(nn[7]);

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

        double totmas;
        _FORTRAN(helmbr8mas)(numgauss,d,numdof,x,y,z,totmas);

        FullSquareMatrix ret(8,d);
        ret /= getProperty()->rho;
        return ret;
}


FullSquareMatrix
HelmBrick::stiffness(const CoordSet &cs, double *Ks, int flg) const
{
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
        Node nd4 = cs.getNode(nn[3]);
        Node nd5 = cs.getNode(nn[4]);
        Node nd6 = cs.getNode(nn[5]);
        Node nd7 = cs.getNode(nn[6]);
        Node nd8 = cs.getNode(nn[7]);

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

        _FORTRAN(helmbrik8v)(x, y, z, numgauss, Ks, numdof, outerr);

        FullSquareMatrix ret(8,Ks);
        ret /= getProperty()->rho;
        return ret;
}


FullSquareMatrix
HelmBrick::acousticm(CoordSet &cs, double *Ks)
{
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
	Node nd4 = cs.getNode(nn[3]);
	Node nd5 = cs.getNode(nn[4]);
	Node nd6 = cs.getNode(nn[5]);
	Node nd7 = cs.getNode(nn[6]);
	Node nd8 = cs.getNode(nn[7]);

	double x[8], y[8], z[8], Kstiff[64];

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

        _FORTRAN(helmbrik8v)(x, y, z, numgauss, Ks, numdof, outerr);

        int i;
        for (i=0;i<64;i++)
          Kstiff[i] = Ks[i];

        double kappa = prop ->kappaHelm;
        double totmas;
        _FORTRAN(helmbr8mas)(numgauss,Ks,numdof,x,y,z,totmas);


	double kk = kappa*kappa;

        for (i=0;i<64;i++)
          Ks[i] = Kstiff[i] - kk*Ks[i];

        FullSquareMatrix ret(8,Ks);
        ret /= getProperty()->rho;
        return ret;
}

int
HelmBrick::numNodes() const
{
 	return 8;
}

int*
HelmBrick::nodes(int *p) const
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
HelmBrick::numDofs() const
{
 	return 8;
}

int*
HelmBrick::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[8];

        p[0] = dsa.locate(nn[0],DofSet::Helm);
        p[1] = dsa.locate(nn[1],DofSet::Helm);
        p[2] = dsa.locate(nn[2],DofSet::Helm);
        p[3] = dsa.locate(nn[3],DofSet::Helm);
        p[4] = dsa.locate(nn[4],DofSet::Helm);
        p[5] = dsa.locate(nn[5],DofSet::Helm);
        p[6] = dsa.locate(nn[6],DofSet::Helm);
        p[7] = dsa.locate(nn[7],DofSet::Helm);

	return p;
}

void
HelmBrick::markDofs(DofSetArray &dsa) const
{
 	dsa.mark(nn[0],DofSet::Helm);
 	dsa.mark(nn[1],DofSet::Helm);
 	dsa.mark(nn[2],DofSet::Helm);
 	dsa.mark(nn[3],DofSet::Helm);
 	dsa.mark(nn[4],DofSet::Helm);
 	dsa.mark(nn[5],DofSet::Helm);
 	dsa.mark(nn[6],DofSet::Helm);
 	dsa.mark(nn[7],DofSet::Helm);
}

int
HelmBrick::getTopNumber() const
{
  return 145;//3;
}



void
HelmBrick::addFaces(PolygonSet *pset)
{
        fprintf(stderr,"HelmBrick::addFaces not implemented.\n");
/*
        pset->addQuad(this,nn[0], nn[1], nn[2], nn[3]);
        pset->addQuad(this,nn[4], nn[5], nn[6], nn[7]); // *
        pset->addQuad(this,nn[3], nn[0], nn[4], nn[7]); // *
        pset->addQuad(this,nn[0], nn[1], nn[5], nn[4]); // *
        pset->addQuad(this,nn[2], nn[1], nn[5], nn[6]);
        pset->addQuad(this,nn[3], nn[2], nn[6], nn[7]);
*/
}

int HelmBrick::getDecFace(int iFace, int *fn) {
 switch(iFace) {
  case 0: fn[0] = nn[0];  fn[1] = nn[1]; fn[2] = nn[2]; fn[3] = nn[3]; break;
  case 1: fn[0] = nn[4];  fn[1] = nn[5]; fn[2] = nn[6]; fn[3] = nn[7]; break;
  case 2: fn[0] = nn[3];  fn[1] = nn[0]; fn[2] = nn[4]; fn[3] = nn[7]; break;
  case 3: fn[0] = nn[0];  fn[1] = nn[1]; fn[2] = nn[5]; fn[3] = nn[4]; break;
  case 4: fn[0] = nn[2];  fn[1] = nn[1]; fn[2] = nn[5]; fn[3] = nn[6]; break;
  default:
  case 5: fn[0] = nn[3];  fn[1] = nn[2]; fn[2] = nn[6]; fn[3] = nn[7]; break;
 }
 return 4;
}

