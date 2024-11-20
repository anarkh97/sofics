#include        <Element.d/ThermQuad.d/ThermQuadGal.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>
#include        <Element.d/State.h>

extern "C"      {
void    _FORTRAN(thermquad3b)(double*, double*, double&, double&, double&, 
                              double*, const int *, const int *);
void    _FORTRAN(q4d1dofmas)(double*, double*, const int&, double*, const int&);
void    _FORTRAN(q4maslumpheat)(double*, double*, const int&, double*, 
                                const int&);
void    _FORTRAN(htsas2)(char*, double*, double*, double&, double*, double*, 
                         double*, const int&, const int&, const int&, const int&);
}


// Four Node Quad element (Galerkin)

ThermQuadGal::ThermQuadGal(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
}

Element *
ThermQuadGal::clone()
{
 return new ThermQuadGal(*this);
}

void
ThermQuadGal::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

void
ThermQuadGal::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

double
ThermQuadGal::getMass(const CoordSet&) const
{
 return 0.0;
}

FullSquareMatrix
ThermQuadGal::massMatrix(const CoordSet &cs, double *d, int cmflg) const
{
        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);

        int i;
        double x[4], y[4], mm[16];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;
        x[3] = nd4.x; y[3] = nd4.y;

        const int numgauss = 2;
        const int numdof   = 4;

        double capacitance = prop->rho*prop->Q*prop->eh;

// Consistent mass

//        _FORTRAN(q4d1dofmas)(x, y, numgauss, mm, numdof);

// Lumped mass

        _FORTRAN(q4maslumpheat)(x,y,numgauss, mm, numdof);

        for (i=0;i<16;i++)
         d[i] = capacitance*mm[i];

        FullSquareMatrix mass(4,d);

        return mass;
}

FullSquareMatrix
ThermQuadGal::stiffness(const CoordSet &cs,double *Ks, int) const
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
        
        double k = prop ->k;
        double h = prop ->eh;
        double c = 0.;

        _FORTRAN(thermquad3b)(x, y, k, c, h, Kstiff, &numgauss, &numdof);

        for (i=0;i<16;i++)
         Ks[i] = Kstiff[i];

        FullSquareMatrix ret(4, Ks);

        return ret;
}

int
ThermQuadGal::numNodes() const
{
 	return 4;
}

int*
ThermQuadGal::nodes(int *p) const
{
 	if(p == 0) p = new int[4];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
	return p;
}

int
ThermQuadGal::numDofs() const
{
 	return 4;
}

int*
ThermQuadGal::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[4];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);
        p[2] = dsa.locate(nn[2],DofSet::Temp);
        p[3] = dsa.locate(nn[3],DofSet::Temp);

	return p;
}

void
ThermQuadGal::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 4, DofSet::Temp);
}

int
ThermQuadGal::getTopNumber() const
{
  return 110; // 2
}

void
ThermQuadGal::computeTemp(CoordSet&cs,
      State &state, double gp[2], double*tres)
{
// 4 is for the number of nodes, 2 is for temp and its derivative
// with respect to time
 double Temp[4][2];

 state.getTemp(nn[0], Temp[0], Temp[0]+1);
 state.getTemp(nn[1], Temp[1], Temp[1]+1);
 state.getTemp(nn[2], Temp[2], Temp[2]+1);
 state.getTemp(nn[3], Temp[3], Temp[3]+1);

/*   fprintf(stderr, "TEMP iS : %14.5e\n", Temp[0][0]);
   fprintf(stderr, "TEMP iS : %14.5e\n", Temp[1][0]);
   fprintf(stderr, "TEMP iS : %14.5e\n", Temp[2][0]);
   fprintf(stderr, "TEMP iS : %14.5e\n", Temp[3][0]); */

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
ThermQuadGal::computeHeatFluxes(Vector& heatflux, CoordSet &cs, Vector& elTemp,
                                   int hgInd)

{
// Heat FLuxes per area
// ... For Z-Direction :
// hgInd are defined in TempDynam.C

   if(hgInd ==2 || hgInd == 5) {
      heatflux = 0.;
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

   double k = prop ->k;

   double elGradTemp[4][6], elHeatFlux[4][6];

   char ESCM[7] = "DIRECT"; // ... DIRECT Heat Fluxes CALCULATION
   //char ESCM[7] = "EXTRAP";   // ... EXTRAPOLATION FROM GAUSS POINTS

//   elTemp.print();

   _FORTRAN(htsas2)(ESCM, x, y, k, elTemp.data(), (double*)elGradTemp, 
                    (double*)elHeatFlux, maxgus, maxstr, elm, numel);

   //..hgInd=0 : heatflux-x             hgInd=3 : gradtemp-x
   //..hgInd=1 : heatflux-y             hgInd=4 : gradtemp-y
   //..hgInd=2 : heatflux-z = 0.        hgInd=5 : gradtemp-z = 0.

    if(hgInd<=2) {
     heatflux[0] = elHeatFlux[0][hgInd];
     heatflux[1] = elHeatFlux[1][hgInd];
     heatflux[2] = elHeatFlux[2][hgInd];
     heatflux[3] = elHeatFlux[3][hgInd];
    }
    if(hgInd>=3) {
     heatflux[0] = elGradTemp[0][hgInd-3];
     heatflux[1] = elGradTemp[1][hgInd-3];
     heatflux[2] = elGradTemp[2][hgInd-3];
     heatflux[3] = elGradTemp[3][hgInd-3];
    }

}
   
