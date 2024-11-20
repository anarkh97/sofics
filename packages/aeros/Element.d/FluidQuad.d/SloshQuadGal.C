#include        <Element.d/FluidQuad.d/SloshQuadGal.h>
#include        <Math.d/FullSquareMatrix.h>
#include        <Utils.d/dofset.h>
#include        <Utils.d/linkfc.h>
#include        <Element.d/State.h>

extern "C"      {
// Overload thermquad3b for the stiffness operator of the sloshing problem
void    _FORTRAN(thermquad3b)(double*, double*, double&, double&, double&, 
                              double*, const int *, const int *);
void    _FORTRAN(slsas2)(char*, double*, double*, double&, double*, double*, 
                         int&, const int&, const int&, const int&);
}


// Four Node Quad element (Galerkin)

SloshQuadGal::SloshQuadGal(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
}

Element *
SloshQuadGal::clone()
{
 return new SloshQuadGal(*this);
}

void
SloshQuadGal::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

void
SloshQuadGal::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

double
SloshQuadGal::getMass(const CoordSet& cs) const
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
SloshQuadGal::massMatrix(const CoordSet &cs, double *d, int cmflg) const
{
        FullSquareMatrix mass(4,d);
        
        mass.zero();
        return mass;
}

FullSquareMatrix
SloshQuadGal::stiffness(const CoordSet &cs,double *Ks, int) const
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
        if(h == 0) { // PJSA 12/2/2014
          std::cerr << " *** ERROR: SloshQuadGal element (type 301) has zero thickness.\n";
        }
        double c = 0.;

        _FORTRAN(thermquad3b)(x, y, k, c, h, Kstiff, &numgauss, &numdof);

        for (i=0;i<16;i++)
         Ks[i] = Kstiff[i];

        FullSquareMatrix ret(4, Ks);

        return ret;
}

int
SloshQuadGal::numNodes() const
{
 	return 4;
}

int*
SloshQuadGal::nodes(int *p) const
{
 	if(p == 0) p = new int[4];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
	return p;
}

int
SloshQuadGal::numDofs() const
{
 	return 4;
}

int*
SloshQuadGal::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[4];

        p[0] = dsa.locate(nn[0],DofSet::Potential);
        p[1] = dsa.locate(nn[1],DofSet::Potential);
        p[2] = dsa.locate(nn[2],DofSet::Potential);
        p[3] = dsa.locate(nn[3],DofSet::Potential);

	return p;
}

void
SloshQuadGal::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 4, DofSet::Potential);
}

int
SloshQuadGal::getTopNumber() const
{
  return 110; // 2
}

void
SloshQuadGal::computeSloshDisp(Vector& fluidDispSlosh, CoordSet &cs, Vector& elPotSlosh,
                                   int hgInd)

{
// Fluid Displacement
// ... For Z-Direction :
// hgInd are defined in Eigen.C

   if(hgInd ==2) {
      fluidDispSlosh = 0.;
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

   double k = 1.;

   double elFluidDispSlosh[4][6];

   char ESCM[7] = "DIRECT"; // ... DIRECT Heat Fluxes CALCULATION
   //char ESCM[7] = "EXTRAP";   // ... EXTRAPOLATION FROM GAUSS POINTS

   _FORTRAN(slsas2)(ESCM, x, y, k, elPotSlosh.data(),
                    (double*)elFluidDispSlosh, maxgus, maxstr, elm, numel);

   //..hgInd=0 : fluidDispSlosh-x     
   //..hgInd=1 : fluidDispSlosh-y    
   //..hgInd=2 : fluidDispSlosh-z = 0.

     fluidDispSlosh[0] = elFluidDispSlosh[0][hgInd];
     fluidDispSlosh[1] = elFluidDispSlosh[1][hgInd];
     fluidDispSlosh[2] = elFluidDispSlosh[2][hgInd];
     fluidDispSlosh[3] = elFluidDispSlosh[3][hgInd];
}

void
SloshQuadGal::computeSloshDispAll(Vector& fluidDispSlosh, CoordSet &cs, Vector& elPotSlosh)

{
// Fluid Displacement
// ... For Z-Direction :
// hgInd are defined in Eigen.C

//   if(hgInd ==2) {
//      fluidDispSlosh = 0.;
//      return;
//   }
 
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

   double k = 1.;

   double elFluidDispSlosh[4][6];

   char ESCM[7] = "DIRECT"; // ... DIRECT Heat Fluxes CALCULATION
   //char ESCM[7] = "EXTRAP";   // ... EXTRAPOLATION FROM GAUSS POINTS

   _FORTRAN(slsas2)(ESCM, x, y, k, elPotSlosh.data(),
                    (double*)elFluidDispSlosh, maxgus, maxstr, elm, numel);

   //..hgInd=0 : fluidDispSlosh-x     
   //..hgInd=1 : fluidDispSlosh-y    
   //..hgInd=2 : fluidDispSlosh-z = 0.

     fluidDispSlosh[0] = elFluidDispSlosh[0][0];
     fluidDispSlosh[1] = elFluidDispSlosh[0][1];
     fluidDispSlosh[2] = 0.;
     fluidDispSlosh[3] = elFluidDispSlosh[1][0];
     fluidDispSlosh[4] = elFluidDispSlosh[1][1];
     fluidDispSlosh[5] = 0.;
     fluidDispSlosh[6] = elFluidDispSlosh[2][0];
     fluidDispSlosh[7] = elFluidDispSlosh[2][1];
     fluidDispSlosh[8] = 0.;
     fluidDispSlosh[9] = elFluidDispSlosh[3][0];
     fluidDispSlosh[10] = elFluidDispSlosh[3][1];
     fluidDispSlosh[11] = 0.;
}
