#include        <cmath>
#include	<Element.d/Helm.d/HelmTri3Gal.h>
#include        <Math.d/matrix.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Math.d/Vector.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>

HelmTri3Gal::HelmTri3Gal(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
}

Element *
HelmTri3Gal::clone()
{
 return new HelmTri3Gal(*this);
}

void
HelmTri3Gal::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

void
HelmTri3Gal::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

double
HelmTri3Gal::getMass(const CoordSet& cs) const
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
HelmTri3Gal::massMatrix(const CoordSet &cs,double *mel,int cmflg) const
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

        FullSquareMatrix M(3,mel);

        double dd = area/6.0;
        double ee = area/12.0;
        if (area<=0.0) fprintf(stderr,"Area less than zero: %f %f\n",x[0],y[0]);

        M[0][0] = dd;
        M[0][1] = ee;
        M[0][2] = ee;

        M[1][0] = ee;
        M[1][1] = dd;
        M[1][2] = ee;

        M[2][0] = ee;
        M[2][1] = ee;
        M[2][2] = dd;

        M /= getProperty()->rho; 
        return M;
}

void
HelmTri3Gal::getHelmForce(CoordSet& cs, ComplexVector &vc, ComplexVector &force)
{
	double mass = getMass(cs);

        int i;
        for(i=0;i<3;i++) force[i] = mass/12.0*(vc[0] + vc[1] + vc[2]);
        for(i=0;i<3;i++) force[i] += mass/12.0*vc[i]; 
}

FullSquareMatrix
HelmTri3Gal::stiffness(const CoordSet &cs, double *Ks, int flg ) const
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

        double x21 = x[1] - x[0];
        double x32 = x[2] - x[1];
        double x13 = x[0] - x[2];

        double y12 = y[0] - y[1];
        double y23 = y[1] - y[2];
        double y31 = y[2] - y[0];

        FullSquareMatrix K(3,Ks);
        K.zero();

        double ke = 0.25/area;

        K[0][0] = ke*(x32*x32 + y23*y23);
        K[0][1] = ke*(x13*x32 + y23*y31);
        K[0][2] = ke*(x21*x32 + y12*y23);

        K[1][0] = K[0][1];
        K[1][1] = ke*(x13*x13 + y31*y31);
        K[1][2] = ke*(x13*x21 + y12*y31);

        K[2][0] = K[0][2];
        K[2][1] = K[1][2];
        K[2][2] = ke*(x21*x21 + y12*y12);

        K /= getProperty()->rho; 
        return K;
}



FullSquareMatrix
HelmTri3Gal::acousticm(CoordSet &cs,double *d)
{
	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);

	double x[3], y[3];
	x[0] = nd1.x; y[0] = nd1.y; 
	x[1] = nd2.x; y[1] = nd2.y;
	x[2] = nd3.x; y[2] = nd3.y; 

        // Mine
	double area = 0.5*((x[1]*y[2]-x[2]*y[1])+
                           (x[2]*y[0]-x[0]*y[2])+
                           (x[0]*y[1]-x[1]*y[0]));

	double x21 = x[1] - x[0];
	double x32 = x[2] - x[1];
	double x13 = x[0] - x[2];

        double y12 = y[0] - y[1];
        double y23 = y[1] - y[2];
        double y31 = y[2] - y[0];

	// Frederic
	/*
	double x12 = x[1] - x[0];
	double x23 = x[2] - x[1];
	double x13 = x[2] - x[0];
	double y12 = y[1] - y[0];
	double y23 = y[2] - y[1];
	double y13 = y[2] - y[0];
	*/

        FullSquareMatrix K(3,d);

	K.zero();

	double ke = 0.25/area;
	//double ke = 0.50/(x12*y13-y12*x13);

	// Mine
	K[0][0] = ke*(x32*x32 + y23*y23);
	K[0][1] = ke*(x13*x32 + y23*y31);
	K[0][2] = ke*(x21*x32 + y12*y23);

	K[1][0] = K[0][1]; 
	K[1][1] = ke*(x13*x13 + y31*y31);
	K[1][2] = ke*(x13*x21 + y12*y31);

	K[2][0] = K[0][2]; 
	K[2][1] = K[1][2]; 
	K[2][2] = ke*(x21*x21 + y12*y12);

        // Frederic
	/*
        K[0][0] = ke*(x23*x23+y23*y23);
        K[0][1] = -ke*(x13*x23+y13*y23);
        K[0][2] = ke*(x12*x23+y12*y23);

        K[1][0] = K[0][1];
        K[1][1] = ke*(x13*x13+y13*y13); 
        K[1][2] = -ke*(x12*x13+y12*y13); 

        K[2][0] = K[0][2];
        K[2][1] = K[1][2];
        K[2][2] = (x12*x12+y12*y12); 
	*/

	FullSquareMatrix M(3,d);

	//double r = x12*y13-y12*x13; // Frederic

	// Mine
	double dd = area/6.0;
	double ee = area/12.0;
        if (area<=0.0) fprintf(stderr,"Area less than zero: %f %f\n",x[0],y[0]);


	// Frederic
	//double dd = r/12.0; 
	//double ee = r/24.0; 

	double kappa = prop ->kappaHelm;

	K[0][0] -= kappa*kappa*dd; 
        K[0][1] -= kappa*kappa*ee; 
        K[0][2] -= kappa*kappa*ee; 

        K[1][0] -= kappa*kappa*ee; 
        K[1][1] -= kappa*kappa*dd; 
        K[1][2] -= kappa*kappa*ee; 

        K[2][0] -= kappa*kappa*ee; 
        K[2][1] -= kappa*kappa*ee; 
        K[2][2] -= kappa*kappa*dd; 

	//fprintf(stderr,"\n\n");
	//K.print();
	//fprintf(stderr,"\n\n");

        K /= getProperty()->rho; 
        return K;
}

int
HelmTri3Gal::numNodes() const
{
 	return 3;
}

int*
HelmTri3Gal::nodes(int *p) const
{
 	if(p == 0) p = new int[3];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
	return p;
}

int
HelmTri3Gal::numDofs() const
{
 	return 3;
}

int*
HelmTri3Gal::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[3];

	dsa.number(nn[0], DofSet::Helm, p);
	dsa.number(nn[1], DofSet::Helm, p+1);
	dsa.number(nn[2], DofSet::Helm, p+2);

	return p;
}

void
HelmTri3Gal::markDofs(DofSetArray &dsa) const
{
 	dsa.mark(nn[0], DofSet::Helm);
 	dsa.mark(nn[1], DofSet::Helm);
 	dsa.mark(nn[2], DofSet::Helm);
}




void
HelmTri3Gal::addFaces(PolygonSet *pset)
{
        fprintf(stderr,"HelmTri3Gal::addFaces not implemented.\n");
/*
        pset->addLine(this,nn[0], nn[1]);
        pset->addLine(this,nn[1], nn[2]);
        pset->addLine(this,nn[2], nn[0]);
*/
}

int
HelmTri3Gal::getTopNumber() const
{
  return 135;//4;
}


void HelmTri3Gal::computedxdxi(CoordSet &cs, int nint, double (*derivatives)[3][2], Matrix22 *dxdxi, double *det) {

 int i,j,k;
 Node nd[3];
 for(i=0;i<3;i++) nd[i] = cs.getNode(nn[i]);

 double coord[3][2];
 for(i=0;i<3;i++) {
   coord[i][0] = nd[i].x;
   coord[i][1] = nd[i].y;
 }

 for(j=0;j<2;j++) for(k=0;k<2;k++) {
   for(i=0;i<nint;i++) {
     dxdxi[i][j][k] = 0.0;
     int m;
     for(m=0;m<3;m++)
       dxdxi[i][j][k] += derivatives[i][m][k] * coord[m][j];
   }
 }

 for(i=0;i<nint;i++) {
   Matrix22 &x = dxdxi[i];
   det[i] = x[0][0]*x[1][1] - x[1][0]*x[0][1];
 }
}



void HelmTri3Gal::getNormalDeriv(CoordSet&cs,ComplexD *uel, int nsc,
                                  int *sc, ComplexD *grad, 
                                  double kappa, double *waveDir) {

 double dN[1][3][2];
 dN[0][0][0] = -1.0;
 dN[0][0][1] = -1.0;
 dN[0][1][0] = 1.0;
 dN[0][1][1] = 0.0;
 dN[0][2][0] = 0.0;
 dN[0][2][1] = 1.0;

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
 for(l=0;l<3;l++) {
   for(m=0;m<2;m++) {
     for(n=0;n<2;n++)
       grad[m] += inv[n][m] * dN[0][l][n] * uel[l];
   }
 }

 int iVertex[2];
 int i;

 for(i=0;i<3;i++) {
    if (sc[0]==nn[i]) iVertex[0] = i;
    if (sc[1]==nn[i]) iVertex[1] = i;
 }

 Node nd[3];
 for(i=0;i<3;i++) nd[i] = cs.getNode(nn[i]);

 double coord[3][3];
 for(i=0;i<3;i++) {
   coord[i][0] = nd[i].x;
   coord[i][1] = nd[i].y;
   coord[i][2] = nd[i].z;
 }

 double point[3] = {0.0,0.0,0.0};
 for(m=0;m<3;m++) {
   for(l=0;l<2;l++) {
     point[m] += coord[iVertex[l]][m];
   }
   point[m] /= 2.0;
 }

 ComplexD tmp = exp(ComplexD(0.0,kappa)*
               (waveDir[0]*point[0]+waveDir[1]*point[1]));
 for(m=0;m<2;m++) {
   grad[m] -= -ComplexD(0.0,1.0)*kappa*waveDir[m]*tmp; 
 }
}
