#include        <Element.d/Helm.d/HelmQuad8Gal.h>
#include        <Math.d/matrix.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>

extern "C"      {
void    _FORTRAN(quad8stif1)(double*, double*, const int&, double*, const int&);
};
extern "C"	{
void    _FORTRAN(quad8mass1)(double*, double*, const int&, double*, const int&);
};


// Helmholtz Quad element (Galerkin)

HelmQuad8Gal::HelmQuad8Gal(int* nodenums)
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
HelmQuad8Gal::clone()
{
 return new HelmQuad8Gal(*this);
}

void
HelmQuad8Gal::renum(const int *table)
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
HelmQuad8Gal::renum(EleRenumMap& table)
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
HelmQuad8Gal::getMass(const CoordSet&) const
{
 return 0.0;
}


FullSquareMatrix
HelmQuad8Gal::massMatrix(const CoordSet &cs, double *d, int cmflg) const
{
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
        Node nd4 = cs.getNode(nn[3]);
        Node nd5 = cs.getNode(nn[4]);
        Node nd6 = cs.getNode(nn[5]);
        Node nd7 = cs.getNode(nn[6]);
        Node nd8 = cs.getNode(nn[7]);

        double x[8], y[8];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;
        x[3] = nd4.x; y[3] = nd4.y;
        x[4] = nd5.x; y[4] = nd5.y;
        x[5] = nd6.x; y[5] = nd6.y;
        x[6] = nd7.x; y[6] = nd7.y;
        x[7] = nd8.x; y[7] = nd8.y;

        _FORTRAN(quad8mass1)(x, y, 3, d, 8);

        FullSquareMatrix ret(8, d);
        ret /= getProperty()->rho; 
        return ret;
}


FullSquareMatrix
HelmQuad8Gal::stiffness(const CoordSet &cs, double *Ks, int flg) const
{
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
        Node nd4 = cs.getNode(nn[3]);
        Node nd5 = cs.getNode(nn[4]);
        Node nd6 = cs.getNode(nn[5]);
        Node nd7 = cs.getNode(nn[6]);
        Node nd8 = cs.getNode(nn[7]);

        double x[8], y[8];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;
        x[3] = nd4.x; y[3] = nd4.y;
        x[4] = nd5.x; y[4] = nd5.y;
        x[5] = nd6.x; y[5] = nd6.y;
        x[6] = nd7.x; y[6] = nd7.y;
        x[7] = nd8.x; y[7] = nd8.y;

        _FORTRAN(quad8stif1)(x, y, 3, Ks, 8);

        FullSquareMatrix ret(8, Ks);
        ret /= getProperty()->rho; 
        return ret;
}


FullSquareMatrix
HelmQuad8Gal::acousticm(CoordSet &cs,double *K)
{
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
	Node nd4 = cs.getNode(nn[3]);
        Node nd5 = cs.getNode(nn[4]);
        Node nd6 = cs.getNode(nn[5]);
        Node nd7 = cs.getNode(nn[6]);
        Node nd8 = cs.getNode(nn[7]);

        int i;
	double x[8], y[8], Kstiff[64];

	x[0] = nd1.x; y[0] = nd1.y; 
	x[1] = nd2.x; y[1] = nd2.y;
	x[2] = nd3.x; y[2] = nd3.y; 
	x[3] = nd4.x; y[3] = nd4.y;
        x[4] = nd5.x; y[4] = nd5.y;
        x[5] = nd6.x; y[5] = nd6.y;
        x[6] = nd7.x; y[6] = nd7.y;
        x[7] = nd8.x; y[7] = nd8.y;

        // Calculate stiffness matrix here:
        // The mass term is added locally to the stiffness

        // First get the stiffness  
        _FORTRAN(quad8stif1)(x, y, 3, K, 8);


        // Put the stiffness temporarely in Kstiff
        for (i=0;i<64;i++)
          Kstiff[i] = K[i];

        // Then get the mass 
        // For while the wave number = el. Area
        double kappa = prop ->kappaHelm;
        _FORTRAN(quad8mass1)(x, y, 3, K, 8);

        // Now add together [K] - kappa^2 [M] 
        for (i=0;i<64;i++)
            K[i] = Kstiff[i] - (kappa*kappa)*K[i];
       

        FullSquareMatrix ret(8, K);
        ret /= getProperty()->rho; 

        return ret;
}

extern bool useFull;

int
HelmQuad8Gal::numNodes() const
{
  if(useFull)
 	return 8;
  else
    return 4;
}

int*
HelmQuad8Gal::nodes(int *p) const
{
  if(useFull)
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
  else
    {
      	if(p == 0) p = new int[4];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
	return p;
    }
}

int
HelmQuad8Gal::numDofs() const
{
 	return 8;
}

int*
HelmQuad8Gal::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[8];

        dsa.number(nn[0],DofSet::Helm , p);
        dsa.number(nn[1],DofSet::Helm , p+1);
        dsa.number(nn[2],DofSet::Helm , p+2);
        dsa.number(nn[3],DofSet::Helm , p+3);
        dsa.number(nn[4],DofSet::Helm , p+4);
        dsa.number(nn[5],DofSet::Helm , p+5);
        dsa.number(nn[6],DofSet::Helm , p+6);
        dsa.number(nn[7],DofSet::Helm , p+7);

	return p;
}

void
HelmQuad8Gal::markDofs(DofSetArray &dsa) const
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

void
HelmQuad8Gal::addFaces(PolygonSet *pset)
{
        fprintf(stderr,"HelmQuad8Gal::addFaces not implemented.\n");
/*
	pset->addLine(this,nn[0], nn[4]);	
	pset->addLine(this,nn[4], nn[1]);	
	pset->addLine(this,nn[1], nn[5]);	
	pset->addLine(this,nn[5], nn[2]);	
        pset->addLine(this,nn[2], nn[6]);
        pset->addLine(this,nn[6], nn[3]);
        pset->addLine(this,nn[3], nn[7]);
        pset->addLine(this,nn[7], nn[0]);
*/
}
