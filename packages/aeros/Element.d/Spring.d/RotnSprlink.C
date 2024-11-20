#include <Utils.d/dbg_alloca.h>
#include <cstdio>
#include <Element.d/Spring.d/RotnSprlink.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Corotational.d/SpringCorotator.h>

extern "C" {
 void _FORTRAN(mstf22)(double*, double&, double&, double&, double&,
                       double&, double&, double&, double&,
                       double&, double&, double&, double&);
}

RotnSprlink::RotnSprlink(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
}

Element *
RotnSprlink::clone()
{
  return new RotnSprlink(*this);
}

void
RotnSprlink::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
}

void
RotnSprlink::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
}

FullSquareMatrix
RotnSprlink::massMatrix(const CoordSet&, double* mel, int cmflg) const
{
  FullSquareMatrix elementMassMatrix(numDofs(), mel);

  elementMassMatrix.zero();

  return elementMassMatrix;
}

FullSquareMatrix
RotnSprlink::stiffness(const CoordSet&, double* d, int flg) const
{
  _FORTRAN(mstf22)((double*)d, prop->A, prop->E, prop->nu, prop->rho,
                   prop->c, prop->k, prop->eh, prop->P,
                   prop->Ta, prop->Q, prop->W, prop->Ixx);

  FullSquareMatrix kel(6, d);

  return kel;
}

int
RotnSprlink::numNodes() const
{
  return 2;
}

int *
RotnSprlink::nodes(int *p) const
{
  if(p == 0) p = new int[2];
  p[0] = nn[0];
  p[1] = nn[1];
  return p;
}

int
RotnSprlink::numDofs() const
{
  return 6;
}

int*
RotnSprlink::dofs(DofSetArray& dsa, int* p) const
{
  int ndof = 0;
  if(p == 0) p = new int[numDofs()];

  dsa.number(nn[0], DofSet::Xrot, p+ndof++);
  dsa.number(nn[0], DofSet::Yrot, p+ndof++);
  dsa.number(nn[0], DofSet::Zrot, p+ndof++);

  dsa.number(nn[1], DofSet::Xrot, p+ndof++);
  dsa.number(nn[1], DofSet::Yrot, p+ndof++);
  dsa.number(nn[1], DofSet::Zrot, p+ndof++);

  return p;
}

void
RotnSprlink::markDofs(DofSetArray& dsa) const
{
  dsa.mark(nn[0], DofSet::Xrot);
  dsa.mark(nn[1], DofSet::Xrot);

  dsa.mark(nn[0], DofSet::Yrot);
  dsa.mark(nn[1], DofSet::Yrot);

  dsa.mark(nn[0], DofSet::Zrot);
  dsa.mark(nn[1], DofSet::Zrot);
}

int
RotnSprlink::getTopNumber() const
{
  return 122;
}

Corotator*
RotnSprlink::getCorotator(CoordSet& cs, double* kel, int, int)
{
  int flag = 0;
  FullSquareMatrix myStiff = stiffness(cs, kel, flag);

  return new SpringCorotator(nn, cs, 2, prop->A, prop->E, prop->nu,
                             myStiff);
}

