#include <cstdio>
#include <cmath>
#include <cstdlib>

#include <Element.d/MpcElement.d/FsiElement.h>
#include <Driver.d/Domain.h>
#include <Corotational.d/Corotator.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>

FsiElement::FsiElement(LMPCons* _fsi)
{
  fsi = _fsi;
  nn = new int[fsi->nterms + 1];

  // count the number of nodes. also fill node numbers array
  nnodes = 0;
  int i,j;
  for(i = 0; i < fsi->nterms; i++) {
    int nnum = (fsi->terms)[i].nnum;
    bool found = false;
    for(j=0; j<nnodes; ++j) {
      if(nnum == nn[j]) { found = true; break; } 
    }
    if(!found) nn[nnodes++] = nnum;
  }
  ndofs = fsi->nterms + 1;
  nn[nnodes++] = fsi->lmpcnum; // fluid node
  // cerr << " ndofs = " << ndofs << ", nnodes = " << nnodes << ", nodes: "; for(i=0; i<nnodes; ++i) cerr << nn[i] << " "; cerr << endl;

  renumTable = 0;
  prop = new StructProp();
}

void
FsiElement::renum(const int *table)
{
  int i;
  renumTable = new int[fsi->nterms+1];
  for(i=0; i<nnodes; ++i) nn[i] = table[nn[i]];
  for(i = 0; i < fsi->nterms; i++) {
    int nnum = (fsi->terms)[i].nnum;
    renumTable[i] = table[nnum];
  }
  renumTable[fsi->nterms] = table[nn[nnodes-1]];
}


void
FsiElement::renum(EleRenumMap& table)
{
  int i;
  renumTable = new int[fsi->nterms+1];
  for(i=0; i<nnodes; ++i) nn[i] = table[nn[i]];
  for(i = 0; i < fsi->nterms; i++) {
    int nnum = (fsi->terms)[i].nnum;
    renumTable[i] = table[nnum];
  }
  renumTable[fsi->nterms] = table[nn[nnodes-1]];
}


FullSquareMatrix
FsiElement::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
  FullSquareMatrix ret(ndofs,mel);
  ret.zero();

  for(int i = 0; i < fsi->nterms; i++) {
    if((fsi->terms)[i].isComplex)
      ret[ndofs-1][i] = -(fsi->terms)[i].coef.c_value.real();
    else
      ret[ndofs-1][i] = -(fsi->terms)[i].coef.r_value;
  }

  return ret;
}

FullSquareMatrix
FsiElement::stiffness(const CoordSet &cs, double *k, int flg) const
{
  FullSquareMatrix ret(ndofs,k);
  ret.zero();
  
  for(int i = 0; i < fsi->nterms; i++) {
    if((fsi->terms)[i].isComplex)
      ret[i][ndofs-1] = (fsi->terms)[i].coef.c_value.real();
    else
      ret[i][ndofs-1] = (fsi->terms)[i].coef.r_value;
  }

  return ret;
}

FullSquareMatrix
FsiElement::imagStiffness(CoordSet &cs, double *k, int flg)
{
  FullSquareMatrix ret(ndofs,k);
  ret.zero();

  for(int i = 0; i < fsi->nterms; i++) {
    if((fsi->terms)[i].isComplex)
      ret[i][ndofs-1] = (fsi->terms)[i].coef.c_value.imag();
  }

  return ret;
}

int *
FsiElement::nodes(int *p) const
{
  if(p == 0) p = new int[nnodes];
  int i;
  for(i=0; i<nnodes; ++i) p[i] = nn[i];
 
  // cerr << "nnodes = " << nnodes << ", FsiElement nodes: "; for(i=0; i<nnodes; ++i) cerr << p[i] << " "; cerr << endl; 
  return p;
}

int *
FsiElement::dofs(DofSetArray &dsa, int *p) const
{
  if(p == 0) p = new int[ndofs];

  int i;
  for(i = 0; i < fsi->nterms; i++) {
    int nnum = (fsi->terms)[i].nnum;
    if(renumTable) nnum = renumTable[i];
    int dofnum = (fsi->terms)[i].dofnum;
    dsa.number(nnum, (1 << dofnum), p+i);
  }
  dsa.number(nn[nnodes-1],DofSet::Helm, p+ndofs-1);

  // cerr << "FsiElement dofs: "; for(i=0; i<ndofs; ++i) cerr << p[i] << " "; cerr << endl;
  return p;
}

void
FsiElement::markDofs(DofSetArray &dsa) const
{
  int i;
  for(i = 0; i < fsi->nterms; i++) {
    int nnum = (fsi->terms)[i].nnum;
    if(renumTable) nnum = renumTable[i];
    int dofnum = (fsi->terms)[i].dofnum;
    dsa.mark(nnum, (1 << dofnum));
  }
  dsa.mark(nn[nnodes-1],DofSet::Helm);
}

int
FsiElement::getElementType() const
{
	return 1002;
}