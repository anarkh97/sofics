#include <Element.d/Membrane.d/Membrane.h>
#include <Element.d/Membrane.d/FourNodeMembrane.h>

FourNodeMembrane::FourNodeMembrane(int *nodenums)
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
    subElems[i] = new Membrane(tmp);
  }
  nnodes = 4;
  ndofs = 24;
}

int
FourNodeMembrane::getTopNumber() const
{
  return 188;
}

