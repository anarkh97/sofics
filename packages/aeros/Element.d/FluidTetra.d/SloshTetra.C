#include	<Element.d/FluidTetra.d/SloshTetra.h>
#include        <Math.d/matrix.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>
#include 	<cstdio>


SloshTetra::SloshTetra(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
}

Element *
SloshTetra::clone()
{
	return new SloshTetra(*this);
}

void
SloshTetra::renum(const int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
	nn[2] = table[nn[2]];
	nn[3] = table[nn[3]];
}

void
SloshTetra::renum(EleRenumMap& table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
	nn[2] = table[nn[2]];
	nn[3] = table[nn[3]];
}

double
SloshTetra::getMass(const CoordSet& cs) const
{
	if(prop->rho == 0) {
		// PJSA 12/2/2014: note this can happen if no attribute is provided in the input file for this element.
		// In this case a dummy StructProp is assigned (see GeoSource::setUpData)
		std::cerr << " *** WARNING: Zero density in SloshTetra element (type 311). The contribution of this element to\n"
		          << "     the total mass of the system will be neglected.\n";
		return 0;
	}

	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);
	auto &nd3 = cs.getNode(nn[2]);
	auto &nd4 = cs.getNode(nn[3]);

	Vector r1(3), r2(3), r3(3), r4(3);

	r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
	r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
	r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;
	r4[0] = nd4.x; r4[1] = nd4.y; r4[2] = nd4.z;

	Vector v1(3), v2(3), v3(3), v4(3);

	v1 = r2 - r1;
	v2 = r3 - r1;
	v3 = r4 - r1;

	v4 = v2.cross(v3);
	double volume = fabs(v1*v4) / 6.0;

	double mass = volume*prop->rho;

	return mass;
}

FullSquareMatrix
SloshTetra::massMatrix(const CoordSet &cs,double *mel,int cmflg) const
{
	FullSquareMatrix ma(4,mel);

	ma.zero();

	return ma;
}

FullSquareMatrix
SloshTetra::stiffness(const CoordSet &cs, double *d, int flg) const
{
	// Calculate the ELEMENT stiffness matrix here:

	// Element stiffness
	FullSquareMatrix sm(4,d);

	double gN[4][3];
	double dOmega = computeMetrics(cs, gN);

	// Get the TETRAHEDRA stiffness matrix
	double TetraStiff[4][4];
	buildTetraStiff(TetraStiff, gN, 0);

	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++) {
			sm[i][j] = TetraStiff[i][j];
		}

	return sm;
}

int
SloshTetra::numNodes() const
{
	return 4;
}

int*
SloshTetra::nodes(int *p) const
{
	if(p == 0) p = new int[4];
	p[0] = nn[0];
	p[1] = nn[1];
	p[2] = nn[2];
	p[3] = nn[3];
	return p;
}

int
SloshTetra::numDofs() const
{
	return 4;
}

int*
SloshTetra::dofs(DofSetArray &dsa, int *p) const
{
	if(p == 0) p = new int[4];

	p[0] = dsa.locate(nn[0],DofSet::Potential);
	p[1] = dsa.locate(nn[1],DofSet::Potential);
	p[2] = dsa.locate(nn[2],DofSet::Potential);
	p[3] = dsa.locate(nn[3],DofSet::Potential);

	return p;
}

void
SloshTetra::markDofs(DofSetArray &dsa) const
{

	dsa.mark(nn[0],DofSet::Potential);
	dsa.mark(nn[1],DofSet::Potential);
	dsa.mark(nn[2],DofSet::Potential);
	dsa.mark(nn[3],DofSet::Potential);
}


double
SloshTetra::computeMetrics(const CoordSet &cs, double gN[4][3]) const
{

	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);
	auto &nd3 = cs.getNode(nn[2]);
	auto &nd4 = cs.getNode(nn[3]);

	double x[4], y[4], z[4];

	// Nodes coordinates

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
	x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

	// Jacobian
	double J[3][3];
	double Jinv[3][3];
	double dOmega;


	J[0][0] = x[1]-x[0]; J[0][1] = x[2]-x[0]; J[0][2] = x[3]-x[0];
	J[1][0] = y[1]-y[0]; J[1][1] = y[2]-y[0]; J[1][2] = y[3]-y[0];
	J[2][0] = z[1]-z[0]; J[2][1] = z[2]-z[0]; J[2][2] = z[3]-z[0];

	double a,b,c,d,e,f,g,h,i;

	a = J[0][0]; b = J[0][1]; c = J[0][2];
	d = J[1][0]; e = J[1][1]; f = J[1][2];
	g = J[2][0]; h = J[2][1]; i = J[2][2];

	// dOmega = determinant(J)

	//Original code
	//dOmega = -a*e*i + a*f*h + d*b*i - d*c*h - g*b*f + g*c*e;

	//New, EC, 20070713
	dOmega = a*e*i - a*f*h - d*b*i + d*c*h + g*b*f - g*c*e;

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
	//Original in Heat
	//Jinv[0][0] = (-e*i+f*h)*t17;
	//Jinv[0][1] = (b*i-c*h)*t17;
	//Jinv[0][2] = -(b*f-c*e)*t17;
	//Jinv[1][0] = -(-d*i+f*g)*t17;
	//Jinv[1][1] = -(a*i-t14)*t17;
	//Jinv[1][2] = (t6-t10)*t17;
	//Jinv[2][0] = -(d*h-e*g)*t17;
	//Jinv[2][1] = (a*h-t12)*t17;
	//Jinv[2][2] = -(t4-t8)*t17;


	//New, EC, 20070713
	Jinv[0][0] = (e*i-f*h)*t17;
	Jinv[0][1] = (c*h-b*i)*t17;
	Jinv[0][2] = (b*f-c*e)*t17;
	Jinv[1][0] = (f*g-d*i)*t17;
	Jinv[1][1] = (a*i-t14)*t17;
	Jinv[1][2] = (t10-t6)*t17;
	Jinv[2][0] = (d*h-e*g)*t17;
	Jinv[2][1] = (t12-a*h)*t17;
	Jinv[2][2] = (t4-t8)*t17;

	a = Jinv[0][0]; b = Jinv[0][1]; c = Jinv[0][2];
	d = Jinv[1][0]; e = Jinv[1][1]; f = Jinv[1][2];
	g = Jinv[2][0]; h = Jinv[2][1]; i = Jinv[2][2];

/*      int m,n,p;
        FullSquareMatrix eyetest(3,0);

        eyetest.zero();

        for (m=0;m<3;m++)
           for (n=0;n<3;n++)
              for (p=0;p<3;p++)  {
                  eyetest[m][n]+=J[m][p]*Jinv[p][n];
              }

        cout << "dOmega = " << dOmega << endl;
        eyetest.print();
        cout<<endl;
*/
	// Shape function gradients 
	// Note: 1st index = shape function #
	//       2nd index = direction (0=x, 1=y, 2=z) 


	gN[0][0] = -(a+d+g);
	gN[0][1] = -(b+e+h);
	gN[0][2] = -(c+f+i);

	gN[1][0] = a;
	gN[1][1] = b;
	gN[1][2] = c;

	gN[2][0] = d;
	gN[2][1] = e;
	gN[2][2] = f;

	gN[3][0] = g;
	gN[3][1] = h;
	gN[3][2] = i;
	return dOmega;
}


void
SloshTetra::buildTetraStiff(double TetraStiff[4][4], double gN[4][3], double dOmega) const
{
	int i,j,k;
	double dot;

	// This is the stiffness matrix for the tetrahedra

	double v = volume(dOmega);
	for (i=0;i<4;i++)
		for (j=0;j<4;j++) {
			dot = 0.0;
			for (k=0;k<3;k++)
				dot += gN[i][k]*gN[j][k];
			TetraStiff[i][j] = v*dot;
		}
}

int
SloshTetra::getTopNumber() const
{
	return 150;//5;
}

