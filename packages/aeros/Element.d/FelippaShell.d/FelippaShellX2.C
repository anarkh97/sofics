#ifdef USE_EIGEN3
#include <Element.d/FelippaShell.d/FelippaShell.h>
#include <Element.d/FelippaShell.d/FelippaShellX2.h>
#include <Hetero.d/FlExchange.h>

// Four node shell comprising two three node shells, 6 dof per node

FelippaShellX2::FelippaShellX2(int *nodenums)
{
  int i,j,k;
  nn = new int[4];
  for(i=0; i<4; ++i) nn[i] = nodenums[i];
 
  nSubElems = 2;
  subElems = new Element * [2];
  subElemNodes = new int * [2];
  subElemDofs = new int * [2];

  subElemNodes[0] = new int[3];
  subElemNodes[0][0] = 0; subElemNodes[0][1] = 1; subElemNodes[0][2] = 3;

  subElemNodes[1] = new int[3];
  subElemNodes[1][0] = 2; subElemNodes[1][1] = 3; subElemNodes[1][2] = 1;

  for(i=0; i<2; ++i) {
    int tmp[3];
    subElemDofs[i] = new int[18];
    for(j=0; j<3; ++j) {
      int nij = subElemNodes[i][j];
      tmp[j] = nodenums[nij]; // global node numbers
      for(k=0;k<6;++k) {
        subElemDofs[i][6*j+k] = 6*nij+k;   
      }
    }
    subElems[i] = new FelippaShell(tmp);
  }
  nnodes = 4;
  ndofs = 24;
}

Element *
FelippaShellX2::clone()
{
  return new FelippaShellX2(*this);
}

int
FelippaShellX2::getTopNumber() const
{
  return 188;
}

void
FelippaShellX2::computeDisp(CoordSet &cs, State &state, const InterpPoint &ip, double *res, GeomState *gs)
{
  const double *gp = ip.xy;
  double gpsum = gp[0] + gp[1];
  int i;
  InterpPoint subip;
  if(gpsum <= 1.) {
    i = 0;
    subip.xy[0] = gp[0];
    subip.xy[1] = gp[1];
  }
  else {
    i = 1;
    subip.xy[0] = 1.0 - gp[0];
    subip.xy[1] = 1.0 - gp[1];
  }

  subElems[i]->computeDisp(cs, state, subip, res);
}

void
FelippaShellX2::getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF,
                         double *res, GeomState *gs)
{
  const double *gp = ip.xy;
  double gpsum = gp[0] + gp[1];
  int i;
  InterpPoint subip;
  if(gpsum <= 1.) {
    i = 0;
    subip.xy[0] = gp[0];
    subip.xy[1] = gp[1];
  }
  else {
    i = 1;
    subip.xy[0] = 1.0 - gp[0];
    subip.xy[1] = 1.0 - gp[1];
  }

  double subres[18];
  subElems[i]->getFlLoad(cs, subip, flF, subres);

  int j;
  for(j=0; j<24; ++j) res[j] = 0.0;
  for(j=0; j<18; ++j) res[subElemDofs[i][j]] = subres[j];
}

#endif
