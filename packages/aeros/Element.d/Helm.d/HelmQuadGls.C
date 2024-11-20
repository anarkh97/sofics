#include        <Element.d/Helm.d/HelmQuadGls.h>
#include        <Math.d/matrix.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>
#include	<cmath>

extern "C"      {
void  _FORTRAN(quad1dofm)(double*, double*, const int&, double*, const int&);
void  _FORTRAN(q4d1dofmas)(double*, double*, const int&, double*, const int&);
}


// Helmholtz Quad element (Galerkin)

HelmQuadGls::HelmQuadGls(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
}

Element *
HelmQuadGls::clone()
{
	return new HelmQuadGls(*this);
}

void
HelmQuadGls::renum(const int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
	nn[2] = table[nn[2]];
	nn[3] = table[nn[3]];
}

void
HelmQuadGls::renum(EleRenumMap& table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
	nn[2] = table[nn[2]];
	nn[3] = table[nn[3]];
}

FullSquareMatrix
HelmQuadGls::massMatrix(const CoordSet &cs, double *d, int cmflg) const
{
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
        Node nd4 = cs.getNode(nn[3]);

        int i;
        double x[4], y[4];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;
        x[3] = nd4.x; y[3] = nd4.y;

        const int numgauss = 2;
        const int numdof   = 4;

        _FORTRAN(q4d1dofmas)(x, y, numgauss, d, numdof);

        FullSquareMatrix ret(4, d);

        double kappa = prop ->kappaHelm;
        double h = 0.0;
        for(i=0;i<16;i++) h += d[i];
        h = sqrt(h);
        double tauksq = 1.0 - 6.0/(kappa*h*kappa*h)*(1-cos(kappa*h))/(2+cos(kappa*h));
        coef = (-kappa*kappa+kappa*kappa*tauksq);

        ret /= getProperty()->rho; 
        return ret;
}


FullSquareMatrix
HelmQuadGls::stiffness(const CoordSet &cs, double *Ks, int flg) const
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
HelmQuadGls::acousticm(CoordSet &cs,double *K)
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

/* 
//RT: The following implementation is due to Antonini
        // GLS parameters
        // Element diameter
        // Using Harari's formula for rectangles only
        double sx, sy, h_ele;

        sx = sqrt( (x[1]-x[0])*(x[1]-x[0])
                 + (y[1]-y[0])*(y[1]-y[0]));

        sy = sqrt( (x[2]-x[1])*(x[2]-x[1])
                 + (y[2]-y[1])*(y[2]-y[1]));

        h_ele = sqrt(2.0)*sx*sy / sqrt(sx*sx + sy*sy);

        double kk = kappa*kappa;
        double csgls = cos(kappa*h_ele);
        double alphagls = 0.5*(1.0-csgls)/(2.0+csgls);
        double alphah = h_ele*h_ele*kk/12.0; 
        double tau = -1.0/(kk) * (1.0-alphagls/alphah);

        // Now add together [K] + (1.0+kk*tau)*kk* [M] 
        // Note the GLS contribution in the mass term
        coef = (1.0+kk*tau)*kk;
        for (i=0;i<16;i++)
            K[i] = Kstiff[i] - coef*K[i];
       

*/

   double h = 0.0;
   for(i=0;i<16;i++) h += K[i];
   h = sqrt(h);

   double tauksq = 1.0 - 6.0/(kappa*h*kappa*h)*(1-cos(kappa*h))/(2+cos(kappa*h));

   coef = (-kappa*kappa+kappa*kappa*tauksq);
   for (i=0;i<16;i++) K[i] = Kstiff[i] + coef*K[i];
        
        FullSquareMatrix ret(4, K);

        ret /= getProperty()->rho; 
        return ret;
}

int
HelmQuadGls::numNodes() const
{
 	return 4;
}

int*
HelmQuadGls::nodes(int *p) const
{
 	if(p == 0) p = new int[4];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
	return p;
}

int
HelmQuadGls::numDofs() const
{
 	return 4;
}

int*
HelmQuadGls::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[4];


        dsa.number(nn[0],DofSet::Helm , p); 
        dsa.number(nn[1],DofSet::Helm , p+1);
        dsa.number(nn[2],DofSet::Helm , p+2);
        dsa.number(nn[3],DofSet::Helm , p+3);

	return p;
}

void
HelmQuadGls::markDofs(DofSetArray &dsa) const
{
 	dsa.mark(nn[0],DofSet::Helm);
 	dsa.mark(nn[1],DofSet::Helm);
 	dsa.mark(nn[2],DofSet::Helm);
 	dsa.mark(nn[3],DofSet::Helm);
}


int
HelmQuadGls::getTopNumber() const
{
  return 131; //2;
}



void
HelmQuadGls::addFaces(PolygonSet *pset)
{
        fprintf(stderr,"HelmQuadGls::addFaces not implemented.\n");
/*
        pset->addLine(this,nn[0], nn[1]);
        pset->addLine(this,nn[1], nn[2]);
        pset->addLine(this,nn[2], nn[3]);
        pset->addLine(this,nn[3], nn[0]);
*/
}


