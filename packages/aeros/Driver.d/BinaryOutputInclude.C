#ifndef _BINARYOUTPUTINCLUDE_C_
#define _BINARYOUTPUTINCLUDE_C_
#include <Driver.d/GeoSource.h>
#include <algorithm>

// ---------------------------------------------------------------------
// binary output routines for double and DComplex
// ---------------------------------------------------------------------

template<class Scalar, int dim>
void GeoSource::writeNodeVectorToFile(SVec<Scalar, dim> &sVec, int glSub,
                int offset, int fileNumber, int iter, int numRes,
                double time, int numComponents, int startComponent, int *glNodeNums)
{
  // this writes numComponents columns of sVec, starting at startComponent
  // for example: displacements -> numComponents = 3, startComponent = 0
  //              temperatures  -> numComponents = 1, startComponent = 6 
  const int numData = numComponents * sVec.size();
  Scalar *data = new Scalar[numData];
  for (int iData = 0; iData < sVec.size(); ++iData) {
    const Scalar *origin = &sVec[iData][startComponent];
    std::copy(origin, origin + numComponents, &data[iData * numComponents]);
  }
    
  writeNodeScalarToFile(data, numData, glSub, offset * numComponents, fileNumber, iter, numRes, 
                        time, numComponents, glNodeNums);

  delete[] data;
}

#endif
