#include <Utils.d/dbg_alloca.h>
#include <iostream>
#include <Element.d/Spring.d/TransSprlink.h>
#include <Corotational.d/SpringCorotator.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>

extern "C" {
 void _FORTRAN(mstf21)(double*,double&, double&,double&,double&,
                               double&, double&,double&,double&,
                               double&, double&,double&,double&);
} 

TransSprlink::TransSprlink(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
}

Element *
TransSprlink::clone()
{
  return new TransSprlink(*this);
}

void
TransSprlink::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
}

void
TransSprlink::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
}

FullSquareMatrix
TransSprlink::massMatrix(const CoordSet &, double *mel, int cmflg) const
{
  FullSquareMatrix elementMassMatrix(numDofs(), mel);

  elementMassMatrix.zero();

  return elementMassMatrix;
}

FullSquareMatrix
TransSprlink::stiffness(const CoordSet &, double *d, int flg) const
{
  double karray[6*6];

  _FORTRAN(mstf21)((double*)karray, prop->A, prop->E, prop->nu, prop->rho,
                   prop->c, prop->k, prop->eh, prop->P,
                   prop->Ta, prop->Q, prop->W, prop->Ixx);

  //
  // karray contains the spring's stiffness matrix
  //

  FullSquareMatrix ret(numDofs(), d);
  FullSquareMatrix kel(6,karray);

  int i,j;
  for(i=0; i<6; ++i) {
    for(j=0; j<6; ++j) {
      ret[i][j] = kel[i][j];
    }
  }
 
  return ret;
}

int
TransSprlink::numNodes() const
{
  return 2;
}

int *
TransSprlink::nodes(int* p) const
{
  if(p == 0) p = new int[2];
  p[0] = nn[0];
  p[1] = nn[1];
  return p;
}

int
TransSprlink::numDofs() const
{
  return 6;
}

int *
TransSprlink::dofs(DofSetArray& dsa, int* p) const
{
  int ndof = 0;
  if(p == 0) p = new int[numDofs()];

  dsa.number(nn[0], DofSet::Xdisp, p+ndof++);
  dsa.number(nn[0], DofSet::Ydisp, p+ndof++);
  dsa.number(nn[0], DofSet::Zdisp, p+ndof++);

  dsa.number(nn[1], DofSet::Xdisp, p+ndof++);
  dsa.number(nn[1], DofSet::Ydisp, p+ndof++);
  dsa.number(nn[1], DofSet::Zdisp, p+ndof++);

  return p;
}

void
TransSprlink::markDofs(DofSetArray& dsa) const
{
  dsa.mark(nn[0], DofSet::Xdisp);
  dsa.mark(nn[1], DofSet::Xdisp);

  dsa.mark(nn[0], DofSet::Ydisp);
  dsa.mark(nn[1], DofSet::Ydisp);

  dsa.mark(nn[0], DofSet::Zdisp);
  dsa.mark(nn[1], DofSet::Zdisp);
}

int
TransSprlink::getTopNumber() const
{
  return 121;
}

Corotator *
TransSprlink::getCorotator(CoordSet& cs, double* kel, int, int)
{
  int flag = 0;
  FullSquareMatrix myStiff = stiffness(cs, kel, flag);

  return new SpringCorotator(cs, nn, prop->A, prop->E, prop->nu,
                             myStiff);
}

