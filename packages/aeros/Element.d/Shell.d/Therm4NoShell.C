#include <Element.d/Shell.d/Therm3NoShell.h>
#include <Element.d/Shell.d/Therm4NoShell.h>
#include <Hetero.d/FlExchange.h>
#include        <Element.d/Element.h>
// Four node thermal shell comprising two three thermal node shells, 1 dof per node

Therm4NoShell::Therm4NoShell(int *nodenums)
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
    subElemDofs[i] = new int[3];
    for(j=0; j<3; ++j) {
      int nij = subElemNodes[i][j];
      tmp[j] = nodenums[nij]; // global node numbers
      for(k=0;k<1;++k) {
        // subElemDofs[i][6*j+k] = 6*nij+k;  
        subElemDofs[i][1*j+k] = 1*nij+k ;
      }
    }
    // cerr << "SubElem " << i << ", nodes " << tmp[0] << " " << tmp[1] << " " << tmp[2] << ", subElemNodes "
    //      << subElemNodes[i][0] << " " << subElemNodes[i][1] << " " << subElemNodes[i][2] << ", subElemDofs ";
    // for(j=0;j<3;++j) cerr << subElemDofs[i][j] << " "; cerr << endl;
    subElems[i] = new Therm3NoShell(tmp);
  }
  nnodes = 4;
  ndofs = 4;
}

Element *
Therm4NoShell::clone()
{
  return new Therm4NoShell(*this);
}

int
Therm4NoShell::getTopNumber() const
{
  return 4746;//2;
}

// the following aero functions have not been tested

void
Therm4NoShell::computeTemp(CoordSet &cs, State &state, double gp[2], double *res)
{
  double gpsum = gp[0] + gp[1];
  int i;
  double subgp[2];
  if(gpsum <= 1.) {
    i = 0; 
    subgp[0] = gp[0];
    subgp[1] = gp[1];
  }
  else {
    i = 1;
    subgp[0] = 1.0 - gp[0];
    subgp[1] = 1.0 - gp[1];
  }
  subElems[i]->computeTemp(cs, state, subgp, res);
}
  
void
Therm4NoShell::getFlFlux(double gp[2], double *flF, double *res)
{
  double gpsum = gp[0] + gp[1];
  int i;
  double subgp[2];
  if(gpsum <= 1.) {
    i = 0;
    subgp[0] = gp[0];
    subgp[1] = gp[1];
  }
  else {
    i = 1;
    subgp[0] = 1.0 - gp[0];
    subgp[1] = 1.0 - gp[1];
  }
  double subres[3];
  subElems[i]->getFlFlux(subgp, flF, subres);

  int j;
  for(j=0; j<4; ++j) res[j] = 0.0;
  for(j=0; j<3; ++j) res[subElemNodes[i][j]] = subres[j];
}

