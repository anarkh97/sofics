#include	<Element.d/Helm.d/HelmTri6Gal.h>
#include        <Math.d/matrix.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Math.d/Vector.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>
#include        <cmath>

extern bool useFull;

extern "C"      {
void    _FORTRAN(trig6stif1)(double*, double*, const int&, double*, const int&);
};
extern "C"      {
void    _FORTRAN(trig6mass1)(double*, double*, const int&, double*, const int&);};


HelmTri6Gal::HelmTri6Gal(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
	nn[4] = nodenums[4];
	nn[5] = nodenums[5];
}

Element *
HelmTri6Gal::clone()
{
 return new HelmTri6Gal(*this);
}

void
HelmTri6Gal::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
  nn[5] = table[nn[5]];

}


void
HelmTri6Gal::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
  nn[4] = table[nn[4]];
  nn[5] = table[nn[5]];

}


double
HelmTri6Gal::getMass(const CoordSet& cs) const
{
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);

        double x[3], y[3];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;

	double area = 0.5*((x[1]*y[2]-x[2]*y[1])+
                           (x[2]*y[0]-x[0]*y[2])+
                           (x[0]*y[1]-x[1]*y[0]));

	double mass = area;

        return mass;
}



FullSquareMatrix
HelmTri6Gal::massMatrix(const CoordSet &cs,double *mel,int cmflg) const
{
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
        Node nd4 = cs.getNode(nn[3]);
        Node nd5 = cs.getNode(nn[4]);
        Node nd6 = cs.getNode(nn[5]);

        double x[6], y[6];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;
        x[3] = nd4.x; y[3] = nd4.y;
        x[4] = nd5.x; y[4] = nd5.y;
        x[5] = nd6.x; y[5] = nd6.y;

        _FORTRAN(trig6mass1)(x, y, 7, mel, 6);

        FullSquareMatrix ret(6, mel);
        ret /= getProperty()->rho; 
        return ret;
}

FullSquareMatrix
HelmTri6Gal::stiffness(const CoordSet &cs, double *Ks, int flg ) const
{
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
        Node nd4 = cs.getNode(nn[3]);
        Node nd5 = cs.getNode(nn[4]);
        Node nd6 = cs.getNode(nn[5]);

        double x[6], y[6];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;
        x[3] = nd4.x; y[3] = nd4.y;
        x[4] = nd5.x; y[4] = nd5.y;
        x[5] = nd6.x; y[5] = nd6.y;

        _FORTRAN(trig6stif1)(x, y, 7, Ks, 6);

        FullSquareMatrix ret(6, Ks);
        ret /= getProperty()->rho; 
        return ret;
}


FullSquareMatrix
HelmTri6Gal::acousticm(CoordSet &cs,double *K)
{
        Node nd1 = cs.getNode(nn[0]);
        Node nd2 = cs.getNode(nn[1]);
        Node nd3 = cs.getNode(nn[2]);
        Node nd4 = cs.getNode(nn[3]);
        Node nd5 = cs.getNode(nn[4]);
        Node nd6 = cs.getNode(nn[5]);

        int i;
        double x[6], y[6], Kstiff[36];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;
        x[3] = nd4.x; y[3] = nd4.y;
        x[4] = nd5.x; y[4] = nd5.y;
        x[5] = nd6.x; y[5] = nd6.y;

        // Calculate stiffness matrix here:
        // The mass term is added locally to the stiffness

        // First get the stiffness
        _FORTRAN(trig6stif1)(x, y, 7, K, 6);


        // Put the stiffness temporarely in Kstiff
        for (i=0;i<36;i++)
          Kstiff[i] = K[i];

        // Then get the mass
        // For while the wave number = el. Area
        double kappa = prop ->kappaHelm;
        _FORTRAN(trig6mass1)(x, y, 7, K, 6);

        // Now add together [K] + kappa^2 [M]
        for (i=0;i<36;i++)
            K[i] = Kstiff[i] - (kappa*kappa)*K[i];

        FullSquareMatrix ret(6, K);

        ret /= getProperty()->rho; 
        return ret;
}

int
HelmTri6Gal::numNodes() const
{
  if(useFull)
    return 6;
  else
    return 3;
}

int*
HelmTri6Gal::nodes(int *p) const
{
  if(useFull)
    {
      
 	if(p == 0) p = new int[6];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
        p[3] = nn[3];
        p[4] = nn[4];
        p[5] = nn[5];
	return p;
    }
  else
    {
	if(p == 0) p = new int[3];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
	return p;
    }
}

int
HelmTri6Gal::numDofs() const
{
 	return 6;
}

int*
HelmTri6Gal::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[6];

	dsa.number(nn[0], DofSet::Helm, p);
	dsa.number(nn[1], DofSet::Helm, p+1);
	dsa.number(nn[2], DofSet::Helm, p+2);
        dsa.number(nn[3], DofSet::Helm, p+3);
        dsa.number(nn[4], DofSet::Helm, p+4);
        dsa.number(nn[5], DofSet::Helm, p+5);

	return p;
}

void
HelmTri6Gal::markDofs(DofSetArray &dsa) const
{
 	dsa.mark(nn[0], DofSet::Helm);
 	dsa.mark(nn[1], DofSet::Helm);
 	dsa.mark(nn[2], DofSet::Helm);
 	dsa.mark(nn[3], DofSet::Helm);
 	dsa.mark(nn[4], DofSet::Helm);
 	dsa.mark(nn[5], DofSet::Helm);
}




void
HelmTri6Gal::addFaces(PolygonSet *pset)
{
        fprintf(stderr,"HelmTri6Gal::addFaces not implemented.\n");
/*
        pset->addLine(this,nn[0], nn[3]);
        pset->addLine(this,nn[3], nn[1]);
        pset->addLine(this,nn[1], nn[4]);
        pset->addLine(this,nn[4], nn[2]);
        pset->addLine(this,nn[2], nn[5]);
        pset->addLine(this,nn[5], nn[0]);
*/
}


void HelmTri6Gal::computedxdxi(CoordSet &cs, int nint, double (*derivatives)[6][2], Matrix22 *dxdxi, double *det) {

 int i,j,k;
 Node nd[6];
 for(i=0;i<6;i++) nd[i] = cs.getNode(nn[i]);

 double coord[6][2];
 for(i=0;i<6;i++) {
   coord[i][0] = nd[i].x;
   coord[i][1] = nd[i].y;
 }

 for(j=0;j<2;j++) for(k=0;k<2;k++) {
   for(i=0;i<nint;i++) {
     dxdxi[i][j][k] = 0.0;
     int m;
     for(m=0;m<6;m++)
       dxdxi[i][j][k] += derivatives[i][m][k] * coord[m][j];
   }
 }

 for(i=0;i<nint;i++) {
   Matrix22 &x = dxdxi[i];
   det[i] = x[0][0]*x[1][1] - x[1][0]*x[0][1];
 }

}



void HelmTri6Gal::getNormalDeriv(CoordSet&cs,ComplexD *uel, int nsc,
                                  int *sc, ComplexD *grad, 
                                  double kappa, double *waveDir) {

 int iVertex[2];
 int i;

 for(i=0;i<3;i++) {
    if (sc[0]==nn[i]) iVertex[0] = i;
    if (sc[2]==nn[i]) iVertex[1] = i;
 }

 if (iVertex[0] > iVertex[1]) { 
   int tmp = iVertex[0]; iVertex[0] = iVertex[1]; iVertex[1] = tmp;
 }

 double r,s;
 if ( (iVertex[0] == 0) && (iVertex[1] == 1) ) {
    r=0.5; s = 0.0;
 } else 
 if ( (iVertex[0] == 1) && (iVertex[1] == 2) ) {
    r=0.5; s = 0.5;
 } else 
 if ( (iVertex[0] == 0) && (iVertex[1] == 2) ) {
    r=0.0; s = 0.5;
 } else { 
   fprintf(stderr, "Error in HelmTri6Gal::getNormalDeriv\n");
 }

 double dN[1][6][2];
 dN[0][0][0] = 1.0 - 4.0*(1.0 - r - s );
 dN[0][1][0] = -1.0 + 4.0*r;
 dN[0][2][0] = 0.0;
 dN[0][3][0] =  -4.0*r + 4.0*(1.0 - r - s );
 dN[0][4][0] = 4.0*s;
 dN[0][5][0] = -4.0*s;

 dN[0][0][1] = 1.0 - 4.0*(1.0 - r - s );
 dN[0][1][1] = 0.0;
 dN[0][2][1] = -1.0 + 4.0*s;
 dN[0][3][1] = -4.0*r;
 dN[0][4][1] = 4.0*r;
 dN[0][5][1] =  -4.0*s + 4.0*(1.0 - r - s );

 double dxdxi[1][2][2], det[1];
 computedxdxi(cs,1,dN,dxdxi,det);
 Matrix22 &x = dxdxi[0];
 double inv[2][2];
 inv[0][0] = x[1][1]/det[0];
 inv[0][1] = -x[0][1]/det[0];
 inv[1][0] = -x[1][0]/det[0];
 inv[1][1] = x[0][0]/det[0];

 grad[0] = grad[1] = grad[2] = ComplexD(0.0,0.0);
 int l,m,n;
 for(l=0;l<6;l++) {
   for(m=0;m<2;m++) {
     for(n=0;n<2;n++)
       grad[m] += inv[n][m] * dN[0][l][n] * uel[l];
   }
 }

 double N[6];
 N[0] = (1.0-r-s)*(2.0*(1.0-r-s)-1.0);
 N[1] = r*(2.0*r-1.0);
 N[2] = s*(2.0*s-1.0);
 N[3] = 4.0*r*(1.0-r-s);
 N[4] = 4.0*r*s;
 N[5] = 4.0*s*(1.0-r-s);

 Node nd[6];
 for(i=0;i<6;i++) nd[i] = cs.getNode(nn[i]);

 double coord[6][3];
 for(i=0;i<6;i++) {
   coord[i][0] = nd[i].x;
   coord[i][1] = nd[i].y;
   coord[i][2] = nd[i].z;
 }
 
 double point[3] = {0.0,0.0,0.0};
 for(m=0;m<2;m++) {
   for(l=0;l<6;l++) {
     point[m] += N[l]*coord[l][m];
   }
 }

 ComplexD tmp = exp(ComplexD(0.0,kappa)*
               (waveDir[0]*point[0]+waveDir[1]*point[1]));
 for(m=0;m<2;m++) {
   grad[m] -= -ComplexD(0.0,1.0)*kappa*waveDir[m]*tmp; 
 }
}
