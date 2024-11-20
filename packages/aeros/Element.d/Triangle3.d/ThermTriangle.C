#include	<Element.d/Triangle3.d/ThermTriangle.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Math.d/Vector.h>
#include        <Utils.d/dofset.h>
#include        <cmath>
#include        <Element.d/State.h>

ThermTriangle::ThermTriangle(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
}


Element *
ThermTriangle::clone()
{
 return new ThermTriangle(*this);
}


void
ThermTriangle::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}


void
ThermTriangle::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}


double
ThermTriangle::getMass(const CoordSet& cs) const
{

     auto &nd1 = cs.getNode(nn[0]);
     auto &nd2 = cs.getNode(nn[1]);
     auto &nd3 = cs.getNode(nn[2]);

     double x[3], y[3];
     x[0] = nd1.x; y[0] = nd1.y;
     x[1] = nd2.x; y[1] = nd2.y;
     x[2] = nd3.x; y[2] = nd3.y;

     double area = 0.5*((x[1]*y[2]-x[2]*y[1])+
                        (x[2]*y[0]-x[0]*y[2])+
                        (x[0]*y[1]-x[1]*y[0]));

     double density = prop->rho;
     double t       = prop->eh;
     double cv      = prop->Q;

     double mass = density*cv*t*area;

     return mass;

}

FullSquareMatrix
ThermTriangle::massMatrix(const CoordSet &cs,double *mel,int cmflg) const
{
	double mass = getMass(cs);
	double massPerNode = mass/3.0;

        FullSquareMatrix ret(3,mel);

	ret.zero();

// This is the LUMPED mass 

	int i;
        for(i=0; i<3; ++i)
          ret[i][i] = massPerNode;

        return ret;
}

FullSquareMatrix
ThermTriangle::stiffness(const CoordSet &cs, double *d, int flg) const
{
	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);
	auto &nd3 = cs.getNode(nn[2]);

	double x[3], y[3];
	x[0] = nd1.x; y[0] = nd1.y; 
	x[1] = nd2.x; y[1] = nd2.y;
	x[2] = nd3.x; y[2] = nd3.y; 

	double area2 = ((x[1]*y[2]-x[2]*y[1])+
                        (x[2]*y[0]-x[0]*y[2])+
                        (x[0]*y[1]-x[1]*y[0]));

        double m21 = (y[1]-y[2])/area2;
        double m22 = (y[2]-y[0])/area2;
        double m23 = (y[0]-y[1])/area2;

        double m31 = (x[2]-x[1])/area2;
        double m32 = (x[0]-x[2])/area2;
        double m33 = (x[1]-x[0])/area2;

	double t  = prop->eh;
        double k = prop ->k;

        FullSquareMatrix K(3,d);

	double ke = k*t*area2/2;

        K[0][0] = ke*(m21*m21 + m31*m31);
        K[0][1] = ke*(m21*m22 + m31*m32);
        K[0][2] = ke*(m23*m21 + m33*m31);
        K[1][0] = ke*(m21*m22 + m31*m32);
        K[1][1] = ke*(m22*m22 + m32*m32);
        K[1][2] = ke*(m22*m23 + m32*m33);
        K[2][0] = ke*(m23*m21 + m33*m31);
        K[2][1] = ke*(m22*m23 + m32*m33);
        K[2][2] = ke*(m23*m23 + m33*m33);

        return K;
}

int
ThermTriangle::numNodes() const
{
 	return 3;
}

int*
ThermTriangle::nodes(int *p) const
{
 	if(p == 0) p = new int[3];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
	return p;
}

int
ThermTriangle::numDofs() const
{
 	return 3;
}

int*
ThermTriangle::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[3];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);
        p[2] = dsa.locate(nn[2],DofSet::Temp);

	return p;
}

void
ThermTriangle::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 3, DofSet::Temp);
}

int
ThermTriangle::getTopNumber() const
{
  return 153;//4;
}

void
ThermTriangle::computeTemp(CoordSet&cs,
      State &state, double gp[2], double*tres)
{
// 3 is for the number of nodes, 2 is for temp and its derivative
// with respect to time
 double Temp[3][2];

 state.getTemp(nn[0], Temp[0], Temp[0]+1);
 state.getTemp(nn[1], Temp[1], Temp[1]+1);
 state.getTemp(nn[2], Temp[2], Temp[2]+1);

/*   fprintf(stderr, "TEMP iS : %14.5e\n", Temp[0][0]);
   fprintf(stderr, "TEMP iS : %14.5e\n", Temp[1][0]);
   fprintf(stderr, "TEMP iS : %14.5e\n", Temp[2][0]); */

// tres[0] = temperature
// tres[1] = d(Temperature)/dt

 int j;
 for(j=0; j<2; ++j)
    tres[j] = (1-gp[0]-gp[1])* Temp[0][j] +
                       gp[0] * Temp[1][j] +
                       gp[1] * Temp[2][j] ;

//     fprintf(stderr, "TEMP1 : %14.5e\n",tres[0]);
//     fprintf(stderr, "DTEMP1: %14.5e\n",tres[1]);
}

void
ThermTriangle::getFlFlux(double gp[2], double *flF, double *tresF)
{
// Projects a fluid flux contained in flF[0] to all 3 nodes of triangle
// Returns tresF
// fprintf(stderr, "Gauss Points %f %f\n ", gp[0], gp[1]);

   tresF[0]  = (1-gp[0]-gp[1])* flF[0];
   tresF[1]  = gp[0] * flF[0];
   tresF[2]  = gp[1] * flF[0];

//   fprintf(stderr, "Fluxes are node 1: %f\n", tresF[0]);
//   fprintf(stderr, "Fluxes are node 2: %f\n", tresF[1]);
//   fprintf(stderr, "Fluxes are node 3: %f\n", tresF[2]);
//   fflush(stderr);
}
