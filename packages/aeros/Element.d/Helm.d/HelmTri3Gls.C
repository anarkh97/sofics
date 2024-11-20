#include        <cmath>
#include        <Element.d/Helm.d/HelmTri3Gls.h>
#include        <Math.d/matrix.h>
#include	<Math.d/FullSquareMatrix.h>
#include        <Math.d/Vector.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>

HelmTri3Gls::HelmTri3Gls(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
}

Element *
HelmTri3Gls::clone()
{
 return new HelmTri3Gls(*this);
}

void
HelmTri3Gls::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

void
HelmTri3Gls::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

double
HelmTri3Gls::getMass(const CoordSet& cs) const
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
HelmTri3Gls::massMatrix(const CoordSet &cs,double *mel, int cmflg) const
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

        FullSquareMatrix M(3, mel);
        double dd = area/6.0;
        double ee = area/12.0;

        M[0][0] = dd;
        M[0][1] = ee;
        M[0][2] = ee;

        M[1][0] = ee;
        M[1][1] = dd;
        M[1][2] = ee;

        M[2][0] = ee;
        M[2][1] = ee;
        M[2][2] = dd;

        double kappa = prop ->kappaHelm;
        double h = 0.0;
        int i,j;
        for(i=0;i<3;i++) for(j=0;j<3;j++) h += M[i][j];
        h = 2.0*sqrt(h/sqrt(3.0));
        double tauksq = 1.0 - 8.0/(kappa*h*kappa*h)*
                        (1.0-cos(sqrt(3.0)/2.0*kappa*h))/(2+cos(sqrt(3.0)/2.0*kappa*h));
        coef = (-kappa*kappa+kappa*kappa*tauksq);

        M /= getProperty()->rho; 
        return M;
}


FullSquareMatrix
HelmTri3Gls::stiffness(const CoordSet &cs, double *Ks, int flg) const
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
HelmTri3Gls::acousticm(CoordSet &cs, double *d)
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

        FullSquareMatrix K(3,d);

	K.zero();

	double ke = 0.25/area;

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

	double M[3][3];
	double dd = area/6.0;
	double ee = area/12.0;

	double kappa = prop ->kappaHelm;

/*
//RT: The following implementation is due to Antonini

	// GLS parameters
        // Element diameter

	double x_centroid = (x[0]+x[1]+x[2])/3.0;
	double y_centroid = (y[0]+y[1]+y[2])/3.0;

	double d1 = sqrt((x[0]-x_centroid)*(x[0]-x_centroid) 
		        +(y[0]-y_centroid)*(y[0]-y_centroid));

        double d2 = sqrt((x[1]-x_centroid)*(x[1]-x_centroid) 
                        +(y[1]-y_centroid)*(y[1]-y_centroid));

        double d3 = sqrt((x[2]-x_centroid)*(x[2]-x_centroid) 
                        +(y[2]-y_centroid)*(y[2]-y_centroid));

	double h_ele = (d1+d2+d3)/3.0;

        double kk = kappa*kappa;

        double csgls = cos(kappa*h_ele);
        double alphagls = 0.5*(1.0-csgls)/(2.0+csgls);
        double alphah = h_ele*h_ele*kk/12.0;
        double tau = -1.0/(kk) * (1.0-alphagls/alphah);



//        double theta = M_PI/12.0;
//        double f = cos(kappa*h_ele*cos(theta))+2.0*cos(kappa*h_ele/2.0*cos(theta))*cos(sqrt(3.0)/2.0*kappa*h_ele*sin(theta));
//        double tau = -(1.0-8.0/(kappa*kappa*h_ele*h_ele)*(3.0-f)/(3.0+f))/kk;

	dd *= (1.0+kk*tau)*kk;
	ee *= (1.0+kk*tau)*kk;

	K[0][0] -= dd; 
        K[0][1] -= ee; 
        K[0][2] -= ee; 

        K[1][0] -= ee; 
        K[1][1] -= dd; 
        K[1][2] -= ee; 

        K[2][0] -= ee; 
        K[2][1] -= ee; 
        K[2][2] -= dd; 
*/
   M[0][0] = dd; 
   M[0][1] = ee; 
   M[0][2] = ee; 

   M[1][0] = ee; 
   M[1][1] = dd; 
   M[1][2] = ee; 

   M[2][0] = ee; 
   M[2][1] = ee; 
   M[2][2] = dd; 

   double h = 0.0;
   int i,j;
   for(i=0;i<3;i++) for(j=0;j<3;j++) h += M[i][j];
fprintf(stderr,"ssss %f\n",h);
   h = 2.0*sqrt(h/sqrt(3.0));

   double tauksq = 1.0 - 8.0/(kappa*h*kappa*h)*
               (1.0-cos(sqrt(3.0)/2.0*kappa*h))/(2+cos(sqrt(3.0)/2.0*kappa*h));

   coef = (-kappa*kappa+kappa*kappa*tauksq);
   for (i=0;i<3;i++) for(j=0;j<3;j++) K[i][j] += coef*M[i][j];

        K /= getProperty()->rho; 
        return K;
}

int
HelmTri3Gls::numNodes() const
{
 	return 3;
}

int*
HelmTri3Gls::nodes(int *p) const
{
 	if(p == 0) p = new int[3];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
	return p;
}

int
HelmTri3Gls::numDofs() const
{
 	return 3;
}

int*
HelmTri3Gls::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[3];

	dsa.number(nn[0], DofSet::Helm, p);
	dsa.number(nn[1], DofSet::Helm, p+1);
	dsa.number(nn[2], DofSet::Helm, p+2);

	return p;
}

void
HelmTri3Gls::markDofs(DofSetArray &dsa) const
{
 	dsa.mark(nn[0], DofSet::Helm);
 	dsa.mark(nn[1], DofSet::Helm);
 	dsa.mark(nn[2], DofSet::Helm);
}



void
HelmTri3Gls::addFaces(PolygonSet *pset)
{
        fprintf(stderr,"HelmTri3Gls::addFaces not implemented.\n");
/*
        pset->addLine(this,nn[0], nn[1]);
        pset->addLine(this,nn[1], nn[2]);
        pset->addLine(this,nn[2], nn[0]);
*/
}

int
HelmTri3Gls::getTopNumber() const
{
  return 136;//4;
}
