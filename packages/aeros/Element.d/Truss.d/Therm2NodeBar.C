#include <cmath>

#include <Element.d/Truss.d/Therm2NodeBar.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Element.d/State.h>


Therm2NodeBar::Therm2NodeBar(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
}

Element *
Therm2NodeBar::clone()
{
	return new Therm2NodeBar(*this);
}

void
Therm2NodeBar::renum(const int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
}

void
Therm2NodeBar::renum(EleRenumMap& table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
}

double
Therm2NodeBar::getMass(const CoordSet& cs) const
{
	auto &nd1 = cs.getNode( nn[0] );
	auto &nd2 = cs.getNode( nn[1] );

	double x[2], y[2], z[2];

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

	double dx = x[1] - x[0];
	double dy = y[1] - y[0];
	double dz = z[1] - z[0];

	double length = sqrt( dx*dx + dy*dy + dz*dz );
	double mass = length*prop->A*prop->rho*prop->Q;

	return mass;
}

FullSquareMatrix
Therm2NodeBar::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
	double mass = getMass(cs);
	double massPerNode = 0.5*mass;

	FullSquareMatrix elementMassMatrix(2,mel);

// zero the element mass matrix

	elementMassMatrix.zero();

// set the diagonal elements

	elementMassMatrix[0][0] = massPerNode;
	elementMassMatrix[1][1] = massPerNode;

	return elementMassMatrix;
}

FullSquareMatrix
Therm2NodeBar::stiffness(const CoordSet &cs, double *Ks, int flg) const
{
	auto &nd1 = cs.getNode( nn[0] );
	auto &nd2 = cs.getNode( nn[1] );

	double x[2], y[2], z[2];

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

	double dx = x[1] - x[0];
	double dy = y[1] - y[0];
	double dz = z[1] - z[0];

	double length = sqrt( dx*dx + dy*dy + dz*dz );

	FullSquareMatrix ret(2,Ks);

// conduction contribution
// A is the average cross sectional area

	double k = prop->k*prop->A/length;
	double &k1 = prop->ymin; // temperature dependent heat flux at node 1: q1 = k1*T1
	double &k2 = prop->ymax; // temperature dependent heat flux at node 2: q2 = k2*T2

	ret[0][0] = k-k1;
	ret[1][1] = k-k2;
	ret[1][0] = -k;
	ret[0][1] = -k;

	return ret;
}

void
Therm2NodeBar::getGravityForce(CoordSet& cs, double *, Vector &force, int, GeomState *)
{
	// compute body source term (not gravity force)
	auto &nd1 = cs.getNode( nn[0] );
	auto &nd2 = cs.getNode( nn[1] );

	double x[2], y[2], z[2];

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

	double dx = x[1] - x[0];
	double dy = y[1] - y[0];
	double dz = z[1] - z[0];

	double length = sqrt( dx*dx + dy*dy + dz*dz );

	double &a = prop->Ixx;
	double &b = prop->Iyy;
	double &c = prop->Izz;

	double Q1 = prop->A*(a*x[0] + b*y[0] + c*z[0]);
	double Q2 = prop->A*(a*x[1] + b*y[1] + c*z[1]);

	force[0] = length/6*(2*Q1 + Q2);
	force[1] = length/6*(Q1 + 2*Q2);
}

int
Therm2NodeBar::numNodes() const
{
	return 2;
}

int *
Therm2NodeBar::nodes(int *p) const
{
	if(p == 0) p = new int[2];
	p[0] = nn[0];
	p[1] = nn[1];
	return p;
}

int
Therm2NodeBar::numDofs() const
{
	return 2;
}

int *
Therm2NodeBar::dofs(DofSetArray &dsa, int *p) const
{
	if(p == 0) p = new int[2];

	p[0] = dsa.locate(nn[0],DofSet::Temp);
	p[1] = dsa.locate(nn[1],DofSet::Temp);

	return p;
}

void
Therm2NodeBar::markDofs(DofSetArray& dsa) const
{
	dsa.mark( nn, 2, DofSet::Temp);
}

int
Therm2NodeBar::getTopNumber() const
{
	return 109; // 1
}

void
Therm2NodeBar::computeTemp(CoordSet&cs,
                           State &state, double gp[2], double*tres)
{
	double Temp[2][2];

	state.getTemp(nn[0], Temp[0], Temp[0]+1);
	state.getTemp(nn[1], Temp[1], Temp[1]+1);

	int j;
	for(j=0; j<2; ++j)
		tres[j] = (1-gp[0])*Temp[0][j] + gp[0]*Temp[1][j] ;
}

void
Therm2NodeBar::getFlFlux(double gp[2], double *flF, double *tresF)
{
// Projects a fluid flux contained in flF[0] to all 2 nodes
// Returns tresF

	tresF[0] = (1-gp[0])*flF[0] ;
	tresF[1] = gp[0]    *flF[0] ;
}

void
Therm2NodeBar::trussHeatFluxes(double &trussflux, CoordSet &cs, Vector& elTemp,
                               int hflInd )
{
// Heat flu per area
	// trussflux must be a referenced parameter!!!

	auto &nd1 = cs.getNode( nn[0] );
	auto &nd2 = cs.getNode( nn[1] );

	double x[2], y[2], z[2];

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

	double dx = x[1] - x[0];
	double dy = y[1] - y[0];
	double dz = z[1] - z[0];

	double length = sqrt( dx*dx + dy*dy + dz*dz );
	double k = prop ->k;

	if (hflInd == 0)
		trussflux = -k*(elTemp[1]-elTemp[0])/length;
	else
		trussflux = (elTemp[1]-elTemp[0])/length;
}
