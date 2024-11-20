#include 	<cstdio>
#include 	<cmath>

#include	<Element.d/Helm.d/TetraHelmGLS.h>
#include        <Math.d/matrix.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Driver.d/Domain.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>


TetraHelmGLS::TetraHelmGLS(int* nodenums)
{
	//nn[0] = nodenums[0];
	//nn[1] = nodenums[1];
	//nn[2] = nodenums[2];
	//nn[3] = nodenums[3];


	// Change order to my notation
	nn[0] = nodenums[2];
	nn[1] = nodenums[1];
	nn[2] = nodenums[0];
	nn[3] = nodenums[3];
}

Element *
TetraHelmGLS::clone()
{
	return new TetraHelmGLS(*this);
}

void
TetraHelmGLS::renum(const int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
	nn[2] = table[nn[2]];
	nn[3] = table[nn[3]];
}

void
TetraHelmGLS::renum(EleRenumMap& table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
	nn[2] = table[nn[2]];
	nn[3] = table[nn[3]];
}

int
TetraHelmGLS::getTopNumber() const
{
	return 140;//5;
}


double
TetraHelmGLS::getMass(const CoordSet& cs) const
{
	return 0.0;
}

FullSquareMatrix
TetraHelmGLS::massMatrix(const CoordSet &cs,double *mel,int cmflg) const
{
	FullSquareMatrix sm(4,mel);
	double gN[4][3];
	double dOmega = computeMetrics(cs, gN);

	double TetraMass[4][4];
	buildTetraMass(TetraMass, dOmega);

	int i, j;
	for (i=0;i<4;i++) for(j=0;j<4;j++)
			sm[i][j] = TetraMass[i][j];

	double kappa = prop ->kappaHelm;
	double h = 0.0;
	for(i=0;i<4;i++) for(j=0;j<4;j++) h += TetraMass[i][j];
	h = pow(6*h*sqrt(2.0),1.0/3.0);
	double tauksq = 1.0 - 8.0/(kappa*h*kappa*h)*
	                      (1.0-cos(sqrt(3.0)/2.0*kappa*h))/(2+cos(sqrt(3.0)/2.0*kappa*h));
	coef = (-kappa*kappa+kappa*kappa*tauksq);

	sm /= getProperty()->rho;
	return sm;
}


FullSquareMatrix
TetraHelmGLS::stiffness(const CoordSet &cs, double *Ks, int flg ) const
{
	FullSquareMatrix sm(4,Ks);

	double gN[4][3];
	double dOmega = computeMetrics(cs, gN);

	double TetraStiff[4][4];
	buildTetraStiff(TetraStiff, gN, dOmega);

	int i,j;
	for (i=0;i<4;i++) for(j=0;j<4;j++)
			sm[i][j] = TetraStiff[i][j];

	sm /= getProperty()->rho;
	return sm;
}



FullSquareMatrix
TetraHelmGLS::acousticm(CoordSet &cs, double *d)
{
	// Element stiffness
	FullSquareMatrix sm(4,d);

	// Volume of the tetrahedra:
	//double V = volume();

	double gN[4][3];
	double dOmega = computeMetrics(cs, gN);

	double TetraStiff[4][4];
	buildTetraStiff(TetraStiff, gN, dOmega);

	// Get the TETRAHEDRA mass matrix	
	double TetraMass[4][4];
	buildTetraMass(TetraMass, dOmega);

	// Get the wave number
	double kappa = prop -> kappaHelm;

	double h = 0.0;
	int i,j;
	for(i=0;i<4;i++) for(j=0;j<4;j++) h += TetraMass[i][j];
	h = pow(6*h*sqrt(2.0),1.0/3.0);

	double tauksq = 1.0 - 8.0/(kappa*h*kappa*h)*
	                      (1.0-cos(sqrt(3.0)/2.0*kappa*h))/(2+cos(sqrt(3.0)/2.0*kappa*h));

	coef = (-kappa*kappa+kappa*kappa*tauksq);
	for (i=0;i<4;i++) for(j=0;j<4;j++)
			sm[i][j] = TetraStiff[i][j] + coef*TetraMass[i][j];

	sm /= getProperty()->rho;
	return sm;
}

int
TetraHelmGLS::numNodes() const
{
	return 4;
}

int*
TetraHelmGLS::nodes(int *p) const
{
	if(p == 0) p = new int[4];
	p[0] = nn[0];
	p[1] = nn[1];
	p[2] = nn[2];
	p[3] = nn[3];
	return p;
}

int
TetraHelmGLS::numDofs() const
{
	return 4;
}

int*
TetraHelmGLS::dofs(DofSetArray &dsa, int *p) const
{
	if(p == 0) p = new int[4];

	dsa.number(nn[0],DofSet::Helm,p);
	dsa.number(nn[1],DofSet::Helm,p+1);
	dsa.number(nn[2],DofSet::Helm,p+2);
	dsa.number(nn[3],DofSet::Helm,p+3);

	return p;
}

void
TetraHelmGLS::markDofs(DofSetArray &dsa) const
{

	dsa.mark(nn[0],DofSet::Helm);
	dsa.mark(nn[1],DofSet::Helm);
	dsa.mark(nn[2],DofSet::Helm);
	dsa.mark(nn[3],DofSet::Helm);
}


/*
double
TetraHelmGLS::volume(CoordSet &cs)
{

 Node nd1 = cs.getNode(nn[0]);
 Node nd2 = cs.getNode(nn[1]);
 Node nd3 = cs.getNode(nn[2]);
 Node nd4 = cs.getNode(nn[3]);

 double X0 = nd2.x - nd1.x;     double Y0 = nd3.x - nd1.x;
 double X1 = nd2.y - nd1.y;     double Y1 = nd3.y - nd1.y;
 double X2 = nd2.z - nd1.z;     double Y2 = nd3.z - nd1.z;
 double Z0 = X1*Y2-X2*Y1;
 double Z1 = X2*Y0-X0*Y2;
 double Z2 = X0*Y1-X1*Y0;
 double vol = (Z0*(nd4.x-nd1.x)+Z1*(nd4.y-nd1.y)+Z2*(nd4.z-nd1.z))/6.0;
 return vol;
}
*/

double
TetraHelmGLS::computeMetrics(const CoordSet &cs, double gN[4][3]) const
{

	Node nd1 = cs.getNode(nn[0]);
	Node nd2 = cs.getNode(nn[1]);
	Node nd3 = cs.getNode(nn[2]);
	Node nd4 = cs.getNode(nn[3]);

	double x[4], y[4], z[4];

	// Nodes coordinates
	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
	x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
	double J[3][3];
	double Jinv[3][3];
	double dOmega;
	// Jacobian
	J[0][0] = x[1]-x[0]; J[0][1] = x[2]-x[0]; J[0][2] = x[3]-x[0];
	J[1][0] = y[1]-y[0]; J[1][1] = y[2]-y[0]; J[1][2] = y[3]-y[0];
	J[2][0] = z[1]-z[0]; J[2][1] = z[2]-z[0]; J[2][2] = z[3]-z[0];

	double a,b,c,d,e,f,g,h,i;

	a = J[0][0]; b = J[0][1]; c = J[0][2];
	d = J[1][0]; e = J[1][1]; f = J[1][2];
	g = J[2][0]; h = J[2][1]; i = J[2][2];

	// dOmega = determinant(J)
	dOmega = -a*e*i + a*f*h + d*b*i - d*c*h - g*b*f + g*c*e;

	// compute Jinv = J^-1; Jinv[i][j] =  dxi_i/dx_j
	// Note: Maple was used here
	double t4,t6,t8,t10,t12,t14,t17;

	t4 = a*e;
	t6 = a*f;
	t8 = b*d;
	t10 = c*d;
	t12 = b*g;
	t14 = c*g;
	t17 = 1/dOmega;
	Jinv[0][0] = (-e*i+f*h)*t17;
	Jinv[0][1] = (b*i-c*h)*t17;
	Jinv[0][2] = -(b*f-c*e)*t17;
	Jinv[1][0] = -(-d*i+f*g)*t17;
	Jinv[1][1] = -(a*i-t14)*t17;
	Jinv[1][2] = (t6-t10)*t17;
	Jinv[2][0] = -(d*h-e*g)*t17;
	Jinv[2][1] = (a*h-t12)*t17;
	Jinv[2][2] = -(t4-t8)*t17;

	a = Jinv[0][0]; b = Jinv[0][1]; c = Jinv[0][2];
	d = Jinv[1][0]; e = Jinv[1][1]; f = Jinv[1][2];
	g = Jinv[2][0]; h = Jinv[2][1]; i = Jinv[2][2];


	// Shape function gradients 
	// Note: 1st index = shape function #
	//       2nd index = direction (0=x, 1=y, 2=z) 

	gN[0][0] = -(a+b+c);
	gN[0][1] = -(d+e+f);
	gN[0][2] = -(g+h+i);

	gN[1][0] = a;
	gN[1][1] = d;
	gN[1][2] = g;

	gN[2][0] = b;
	gN[2][1] = e;
	gN[2][2] = h;

	gN[3][0] = c;
	gN[3][1] = f;
	gN[3][2] = i;

	return dOmega;
}


void
TetraHelmGLS::buildTetraStiff(double TetraStiff[4][4], double gN[4][3], double dOmega) const
{
	int i,j,k;
	double dot;

	double V = volume(dOmega);

	// This is the stiffness matrix for the tetrahedra

	for (i=0;i<4;i++)
		for (j=0;j<4;j++) {
			dot = 0.0;
			for (k=0;k<3;k++)
				dot += gN[i][k]*gN[j][k];
			TetraStiff[i][j] = dot*V;
		}

}

void
TetraHelmGLS::buildTetraMass(double TetraMass[4][4], double dOmega) const
{

	// Volume of the tetrahedra:
	double V = volume(dOmega);

	// Exact integration using Maple or Tetrahedra formulas
	double diag_mass=V/10.0, off_diag_mass=V/20.0;

	int i, j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++) {
			TetraMass[i][j] = off_diag_mass;
			TetraMass[i][i] = diag_mass;
		}

}




void
TetraHelmGLS::addFaces(PolygonSet *pset)
{
	fprintf(stderr,"TetraHelmGLS::addFaces not implemented.\n");
/*
        pset->addTri(this,nn[0], nn[1], nn[2]);
        pset->addTri(this,nn[0], nn[1], nn[3]);
        pset->addTri(this,nn[0], nn[3], nn[2]);
        pset->addTri(this,nn[2], nn[3], nn[1]);
*/
}

int TetraHelmGLS::getDecFace(int iFace, int *fn) {
	switch(iFace) {
		case 0: fn[0] = nn[0];  fn[1] = nn[2]; fn[2] = nn[1]; break;
		case 1: fn[0] = nn[0];  fn[1] = nn[1]; fn[2] = nn[3]; break;
		case 2: fn[0] = nn[0];  fn[1] = nn[3]; fn[2] = nn[2]; break;
		default:
		case 3: fn[0] = nn[2];  fn[1] = nn[3]; fn[2] = nn[1]; break;
	}
	return 3;
}

