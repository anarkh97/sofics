#include        <Element.d/ThermQuad.d/Therm3DQuad.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>

extern "C"      {
void   _FORTRAN(thermquad3a)(double*, double*, double*, double*, double*);
void   _FORTRAN(thermquad3b)(double*, double*, double &, double &, 
                             double &, double*, const int&, const int&);
void    _FORTRAN(q4d1dofmas)(double*, double*, const int&, double*, const int&);
void    _FORTRAN(htsas2)(char*, double*, double*, double&, double*, double*,
                         double*, const int&, const int&, const int&, 
                          const int&);
}


// Four Node Quad element (Galerkin)

Therm3DQuad::Therm3DQuad(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
}

Element *
Therm3DQuad::clone()
{
 return new Therm3DQuad(*this);
}

void
Therm3DQuad::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

void
Therm3DQuad::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

double
Therm3DQuad::getMass(const CoordSet&) const
{
 return 0.0;
}

FullSquareMatrix
Therm3DQuad::massMatrix(const CoordSet &cs, double *d, int cmflg) const
{
        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);

        int i;
        double x[4], y[4], z[4], mm[16];
        double xl[4], yl[4];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

        const int numgauss = 2;
        const int numdof   = 4;

        _FORTRAN(thermquad3a)(x, y, z, xl, yl);

        double capacitance = prop->rho*prop->Q*prop->eh;

        _FORTRAN(q4d1dofmas)(xl, yl, numgauss, mm, numdof);

        for (i=0;i<16;i++)
         d[i] = capacitance*mm[i];

        FullSquareMatrix mass(4,d);

        return mass;
}

FullSquareMatrix
Therm3DQuad::stiffness(const CoordSet &cs, double *Ks, int flg) const
{

	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);
	auto &nd3 = cs.getNode(nn[2]);
	auto &nd4 = cs.getNode(nn[3]);

        int i;
	double x[4], y[4], z[4], Ke[16];
	double xl[4], yl[4];

	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
	x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

        // Calculate stiffness matrix 

	const int numgauss = 2;
	const int numdof   = 4;

        _FORTRAN(thermquad3a)(x, y, z, xl, yl);

        double h = prop ->eh;
        double k = prop ->k;
        double c = 0.;

        _FORTRAN(thermquad3b)(xl, yl, k, c, h, Ke, numgauss, numdof);

        for (i=0;i<16;i++)
         Ks[i] = Ke[i];

        FullSquareMatrix ret(4, Ks);

        return ret;
}

int
Therm3DQuad::numNodes() const
{
 	return 4;
}

int*
Therm3DQuad::nodes(int *p) const
{
 	if(p == 0) p = new int[4];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
	return p;
}

int
Therm3DQuad::numDofs() const
{
  return 4; 
}

int*
Therm3DQuad::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[4];

        p[0] = dsa.locate(nn[0],DofSet::Temp);
        p[1] = dsa.locate(nn[1],DofSet::Temp);
        p[2] = dsa.locate(nn[2],DofSet::Temp);
        p[3] = dsa.locate(nn[3],DofSet::Temp);

	return p;
}

void
Therm3DQuad::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 4, DofSet::Temp);
}

int
Therm3DQuad::getTopNumber() const
{
  return 103; // 2
}

void
Therm3DQuad::computeHeatFluxes(Vector& heatflux, CoordSet &cs, Vector& elTemp,
                                   int hgInd)
{
// Heat FLuxes per area
// hgInd are defined in TempDynam.C
// ... For Z-Direction :

   if(hgInd ==2 || hgInd == 5) {
      heatflux = 0.;
      return;
   }

   auto &nd1 = cs.getNode(nn[0]);
   auto &nd2 = cs.getNode(nn[1]);
   auto &nd3 = cs.getNode(nn[2]);
   auto &nd4 = cs.getNode(nn[3]);

   double x[4], y[4], z[4];
   double xl[4], yl[4];

   x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
   x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
   x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
   x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

   _FORTRAN(thermquad3a)(x, y, z, xl, yl);

   int maxgus = 4;
   int maxstr = 6;
   int elm    = 1;
   int numel  = 1;

   double k = prop ->k;

   double elGradTemp[4][6], elHeatFlux[4][6];

   char ESCM[7] = "DIRECT"; // ... DIRECT Heat Fluxes CALCULATION
   //char ESCM[7] = "EXTRAP";   // ... EXTRAPOLATION FROM GAUSS POINTS

   _FORTRAN(htsas2)(ESCM, xl, yl, k, elTemp.data(), (double*)elGradTemp,
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

