#include        <Element.d/Helm.d/HelmQuadGal.h>
#include        <Math.d/matrix.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>

extern "C"      {
void _FORTRAN(quad1dofm)( double*, double*, const int&, double*, const int&);
void _FORTRAN(q4d1dofmas)(double*, double*, const int&, double*, const int&);
}


// Helmholtz Quad element (Galerkin)

HelmQuadGal::HelmQuadGal(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
}

Element *
HelmQuadGal::clone()
{
 	return new HelmQuadGal(*this);
}

void
HelmQuadGal::renum(const int *table)
{
  	nn[0] = table[nn[0]];
  	nn[1] = table[nn[1]];
  	nn[2] = table[nn[2]];
  	nn[3] = table[nn[3]];
}

void
HelmQuadGal::renum(EleRenumMap& table)
{
  	nn[0] = table[nn[0]];
  	nn[1] = table[nn[1]];
  	nn[2] = table[nn[2]];
  	nn[3] = table[nn[3]];
}

FullSquareMatrix
HelmQuadGal::massMatrix(const CoordSet &cs, double *d, int cmflg) const
{
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
        Node nd4 = cs.getNode(nn[3]);

        double x[4], y[4];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;
        x[3] = nd4.x; y[3] = nd4.y;

        const int numgauss = 2;
        const int numdof   = 4;

        _FORTRAN(q4d1dofmas)(x, y, numgauss, d, numdof);

        FullSquareMatrix ret(4, d);
        ret /= getProperty()->rho; 
        return ret;
}

FullSquareMatrix
HelmQuadGal::stiffness(const CoordSet &cs, double *Ks, int flg ) const
{
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
        Node nd4 = cs.getNode(nn[3]);

        double x[4], y[4];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;
        x[3] = nd4.x; y[3] = nd4.y;

        const int numgauss = 2;
        const int numdof   = 4;

        _FORTRAN(quad1dofm)(x, y, numgauss, Ks, numdof);

        FullSquareMatrix ret(4, Ks);
        ret /= getProperty()->rho; 
        return ret;
}


FullSquareMatrix
HelmQuadGal::acousticm(CoordSet &cs, double *K)
{
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
	Node nd4 = cs.getNode(nn[3]);

        int i;
	double x[4], y[4], Kstiff[16];

	x[0] = nd1.x; y[0] = nd1.y; 
	x[1] = nd2.x; y[1] = nd2.y;
	x[2] = nd3.x; y[2] = nd3.y; 
	x[3] = nd4.x; y[3] = nd4.y;

        // Calculate stiffness matrix here:
        // The mass term is added locally to the stiffness

	const int numgauss = 2;
	const int numdof   = 4;

        // First get the stiffness  
        _FORTRAN(quad1dofm)(x, y, numgauss, K, numdof);


        // Put the stiffness temporarely in Kstiff
        for (i=0;i<16;i++)
          Kstiff[i] = K[i];

        // Then get the mass 
        // For while the wave number = el. Area
        double kappa = prop ->kappaHelm;
        _FORTRAN(q4d1dofmas)(x, y, numgauss, K, numdof);

        // Now add together [K] + kappa^2 [M] 
	double k2 = kappa*kappa;
        for (i=0;i<16;i++)
            K[i] = Kstiff[i] - k2*K[i];

        FullSquareMatrix ret(4, K);

        ret /= getProperty()->rho; 
        return ret;
}

void
HelmQuadGal::getHelmForce(CoordSet& cs, ComplexVector &vc, ComplexVector &force)
{
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
	Node nd4 = cs.getNode(nn[3]);

        int i;
	double x[4], y[4], mass[16];

	x[0] = nd1.x; y[0] = nd1.y; 
	x[1] = nd2.x; y[1] = nd2.y;
	x[2] = nd3.x; y[2] = nd3.y; 
	x[3] = nd4.x; y[3] = nd4.y;


	const int numgauss = 2;
	const int numdof   = 4;

        _FORTRAN(q4d1dofmas)(x, y, numgauss, mass, numdof);
         force.zero();
         for(i=0;i<16;i++) force[i/4] += mass[i]*vc[i%4];
}


int
HelmQuadGal::numNodes() const
{
 	return 4;
}

int*
HelmQuadGal::nodes(int *p) const
{
 	if(p == 0) p = new int[4];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
	return p;
}

int
HelmQuadGal::numDofs() const
{
 	return 4;
}

int*
HelmQuadGal::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[4];

        dsa.number(nn[0],DofSet::Helm , p);
        dsa.number(nn[1],DofSet::Helm , p+1);
        dsa.number(nn[2],DofSet::Helm , p+2);
        dsa.number(nn[3],DofSet::Helm , p+3);

	return p;
}

void
HelmQuadGal::markDofs(DofSetArray &dsa) const
{
 	dsa.mark(nn[0],DofSet::Helm);
 	dsa.mark(nn[1],DofSet::Helm);
 	dsa.mark(nn[2],DofSet::Helm);
 	dsa.mark(nn[3],DofSet::Helm);
}

int
HelmQuadGal::getTopNumber() const
{
  return 130; // 2
}


void
HelmQuadGal::addFaces(PolygonSet *pset)
{
        fprintf(stderr,"HelmQuadGal::addFaces not implemented.\n");
/*
        pset->addLine(this,nn[0], nn[1]);
        pset->addLine(this,nn[1], nn[2]);
        pset->addLine(this,nn[2], nn[3]);
        pset->addLine(this,nn[3], nn[0]);
*/
}

