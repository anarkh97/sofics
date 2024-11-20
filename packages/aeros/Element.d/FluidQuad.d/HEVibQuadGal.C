#include        <Element.d/FluidQuad.d/HEVibQuadGal.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>
#include        <Element.d/State.h>

extern "C"      {
// Overload thermquad3b for the stiffness operator of the HEV problem
void    _FORTRAN(thermquad3b)(double*, double*, double&, double&, double&, 
                              double*, const int *, const int *);
// void    _FORTRAN(q4d1dofmas)(double*, double*, const int&, double*, const int&);
void    _FORTRAN(slsas2)(char*, double*, double*, double&, double*, double*, 
                         int&, const int&, const int&, const int&);
}


// Four Node Quad element (Galerkin)

HEVibQuadGal::HEVibQuadGal(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
}

Element *
HEVibQuadGal::clone()
{
 return new HEVibQuadGal(*this);
}

void
HEVibQuadGal::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

void
HEVibQuadGal::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

double
HEVibQuadGal::getMass(const CoordSet& cs) const
{
  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);
  auto &nd4 = cs.getNode(nn[3]);

  Vector r1(3), r2(3), r3(3), r4(3);

  r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
  r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
  r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;
  r4[0] = nd4.x; r4[1] = nd4.y; r4[2] = nd4.z;

  Vector v1(3), v2(3), v3(3), v4(3), v5(3);

  v1 = r2 - r1;
  v2 = r3 - r1;
  v3 = r4 - r1;

  v4 = v1.cross(v2);
  v5 = v2.cross(v3);

  double area = 0.5*(v4.magnitude() + v5.magnitude());
  double mass = area*prop->rho*prop->eh;

  return mass;
}

FullSquareMatrix
HEVibQuadGal::massMatrix(const CoordSet &cs, double *d, int cmflg) const
{
        FullSquareMatrix mass(4,d);
        
        mass.zero();
        return mass;
}

FullSquareMatrix
HEVibQuadGal::stiffness(const CoordSet &cs,double *Ks, int) const
{

	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);
	auto &nd3 = cs.getNode(nn[2]);
	auto &nd4 = cs.getNode(nn[3]);

        int i;
	double x[4], y[4], Kstiff[16];

	x[0] = nd1.x; y[0] = nd1.y; 
	x[1] = nd2.x; y[1] = nd2.y;
	x[2] = nd3.x; y[2] = nd3.y; 
	x[3] = nd4.x; y[3] = nd4.y;

        // Calculate stiffness matrix 

	const int numgauss = 2;
	const int numdof   = 4;
        
        double k = 1.;
        double h = prop ->eh;
        double c = 0.;
        double rho = prop ->rho;

        _FORTRAN(thermquad3b)(x, y, k, c, h, Kstiff, &numgauss, &numdof);

        for (i=0;i<16;i++)
         Ks[i] = rho*Kstiff[i];

        FullSquareMatrix ret(4, Ks);
        //ret.print();

        return ret;
}

int
HEVibQuadGal::numNodes() const
{
 	return 4;
}

int*
HEVibQuadGal::nodes(int *p) const
{
 	if(p == 0) p = new int[4];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
	return p;
}

int
HEVibQuadGal::numDofs() const
{
 	return 4;
}

int*
HEVibQuadGal::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[4];

        p[0] = dsa.locate(nn[0],DofSet::Potential);
        p[1] = dsa.locate(nn[1],DofSet::Potential);
        p[2] = dsa.locate(nn[2],DofSet::Potential);
        p[3] = dsa.locate(nn[3],DofSet::Potential);

	return p;
}

void
HEVibQuadGal::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 4, DofSet::Potential);
}

int
HEVibQuadGal::getTopNumber() const
{
  return 110; // 2
}

/*
void
HEVibQuadGal::computeTemp(CoordSet&cs,
      State &state, double gp[2], double*tres)
{
// 4 is for the number of nodes, 2 is for temp and its derivative
// with respect to time
 double Temp[4][2];

 state.getTemp(nn[0], Temp[0], Temp[0]+1);
 state.getTemp(nn[1], Temp[1], Temp[1]+1);
 state.getTemp(nn[2], Temp[2], Temp[2]+1);
 state.getTemp(nn[3], Temp[3], Temp[3]+1);

// fprintf(stderr, "TEMP iS : %14.5e\n", Temp[0][0]);
// fprintf(stderr, "TEMP iS : %14.5e\n", Temp[1][0]);
// fprintf(stderr, "TEMP iS : %14.5e\n", Temp[2][0]);
// fprintf(stderr, "TEMP iS : %14.5e\n", Temp[3][0]);

// tres[0] = temperature
// tres[1] = d(Temperature)/dt

 int j;
 for(j=0; j<2; ++j)
    tres[j] = (1-gp[0])*(1-gp[1])* Temp[0][j] +
              gp[0]*(1-gp[1])    * Temp[1][j] +
              gp[0]*gp[1]        * Temp[2][j] +
              (1-gp[0])*gp[1]    * Temp[3][j]; 
//     fprintf(stderr, "TEMP1 : %14.5e\n",tres[0]);
//     fprintf(stderr, "DTEMP1: %14.5e\n",tres[1]);
}

void
ThermQuadGal::getFlFlux(double gp[2], double *flF, double *tresF)
{
// Projects a fluid flux contained in flF[0] to all 4 nodes of quad
// Returns tresF
// fprintf(stderr, "Gauss Points %f %f\n ", gp[0], gp[1]);

   tresF[0]  = (1-gp[0])*(1-gp[1])* flF[0];
   tresF[1]  = gp[0]*(1-gp[1])    * flF[0];
   tresF[2]  = gp[0]*gp[1]        * flF[0];
   tresF[3]  = (1-gp[0])*gp[1]    * flF[0];

//   fprintf(stderr, "Fluxes are node 1: %f\n", tresF[0]);
//   fprintf(stderr, "Fluxes are node 2: %f\n", tresF[1]);
//   fprintf(stderr, "Fluxes are node 3: %f\n", tresF[2]);
//   fprintf(stderr, "Fluxes are node 4: %f\n", tresF[3]);
//   fflush(stderr); 
}

void
HEVibQuadGal::computeHEVibDisp(Vector& fluidDispHEVib, CoordSet &cs, Vector& elPotHEVib,
                                   int hgInd)

{
// Fluid Displacement
// ... For Z-Direction :
// hgInd are defined in Eigen.C

   if(hgInd ==2) {
      fluidDispHEVib = 0.;
      return;
   }
 
   auto &nd1 = cs.getNode(nn[0]);
   auto &nd2 = cs.getNode(nn[1]);
   auto &nd3 = cs.getNode(nn[2]);
   auto &nd4 = cs.getNode(nn[3]);

   double x[4], y[4];

   x[0] = nd1.x; y[0] = nd1.y;
   x[1] = nd2.x; y[1] = nd2.y;
   x[2] = nd3.x; y[2] = nd3.y;
   x[3] = nd4.x; y[3] = nd4.y;

   int maxgus = 4;
   int maxstr = 6;
   int elm    = 1;
   int numel  = 1;

   //double k = prop ->k;
   double k = 1.;

   double elGradTemp[4][6], elFluidDispHEVib[4][6];

   char ESCM[7] = "DIRECT"; // ... DIRECT Heat Fluxes CALCULATION
   //char ESCM[7] = "EXTRAP";   // ... EXTRAPOLATION FROM GAUSS POINTS

//   elPotHEVib.print();

   _FORTRAN(slsas2)(ESCM, x, y, k, elPotHEVib.data(),
                    (double*)elFluidDispHEVib, maxgus, maxstr, elm, numel);

   //..hgInd=0 : fluidDispHEVib-x     
   //..hgInd=1 : fluidDispHEVib-y    
   //..hgInd=2 : fluidDispHEVib-z = 0.

     fluidDispHEVib[0] = elFluidDispHEVib[0][hgInd];
     fluidDispHEVib[1] = elFluidDispHEVib[1][hgInd];
     fluidDispHEVib[2] = elFluidDispHEVib[2][hgInd];
     fluidDispHEVib[3] = elFluidDispHEVib[3][hgInd];
}
*/
