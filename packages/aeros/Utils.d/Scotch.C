#include <cstdio>
#ifdef USE_SCOTCH
#include <stdint.h>
#include "scotch.h"
#endif
#include <Utils.d/Connectivity.h>
#include <iostream>
using namespace std;

Connectivity * 
Connectivity::SCOTCH_graphPart(int partnbr) const
{
  //cerr << "here in Connectivity::SCOTCH_graphPart, partnbr = " << partnbr << ", vertToVert = \n"; print();
#ifdef USE_SCOTCH
  SCOTCH_Graph * grafptr = new SCOTCH_Graph;
  SCOTCH_Num baseval = 0; // the graph base value for index arrays (typically 0 for structures
                          // built from C and 1 for structures built from Fortran).
  SCOTCH_Num vertnbr = size; // number of vertices
  SCOTCH_Num * verttab = new SCOTCH_Num[vertnbr + 1]; // the adjacency index array, of size (vertnbr + 1) if
                                                      // the edge array is compact (that is, if vendtab equals verttab + 1 or NULL),
                                                      // or of size vertnbr else.
  for(size_t i = 0; i < vertnbr + 1; ++i) verttab[i] = pointer[i];
  SCOTCH_Num * vendtab = 0; // vendtab is the adjacency end index array, of size
                            // vertnbr if it is disjoint from verttab.
  SCOTCH_Num * velotab = 0; // the vertex load array, of size vertnbr if it exists.
  SCOTCH_Num * vlbltab = 0; // the vertex label array, of size vertnbr if it exists.
  SCOTCH_Num edgenbr = target.size(); // the number of arcs (that is, twice the number of edges).
  SCOTCH_Num * edgetab = new SCOTCH_Num[edgenbr]; // the adjacency array, of size at least edgenbr (it can be more if the
                                                  // edge array is not compact).
  for(size_t i = 0; i < edgenbr; ++i) edgetab[i] = target[i];
  SCOTCH_Num * edlotab = 0; // the arc load array, of size edgenbr if it exists
  int ierr;

  // NOTE: vendtab, velotab, vlbltab and edlotab arrays are optional, and a NULL
  // pointer can be passed as argument whenever they are not defined.

  ierr = ::SCOTCH_graphBuild(grafptr, baseval, vertnbr, verttab, vendtab, velotab, vlbltab, edgenbr, edgetab, edlotab);
  if(ierr) cerr << "ERROR: SCOTCH_graphBuild unsuccessful\n";

  ierr = ::SCOTCH_graphCheck(grafptr);
  if(ierr) cerr << "ERROR: SCOTCH_graphCheck unsuccessful\n";

  SCOTCH_Strat * straptr = new SCOTCH_Strat; // partitioning strategy
  ierr = ::SCOTCH_stratInit(straptr);
  if(ierr) cerr << "ERROR: SCOTCH_stratInit unsuccessful\n";

  // NOTE: if SCOTCH stratGraphBipart is not called then default strategy will be used
  SCOTCH_Num * parttab = new SCOTCH_Num[vertnbr]; // partition data (returned)
#if defined(_OPENMP)
  #pragma omp critical
#endif
  ierr = ::SCOTCH_graphPart(grafptr, partnbr, straptr, parttab);
  if(ierr) cerr << "ERROR: SCOTCH_graphPart unsuccessful\n";

  int *ptr = new int[vertnbr+1];
  for(int i=0; i<vertnbr+1; ++i) ptr[i] = i;
  Connectivity vertToPart(vertnbr, ptr, parttab);
  Connectivity *partToVert = vertToPart.alloc_reverse();
  //cerr << "partToVert = \n"; partToVert->print();
  delete [] verttab;
  delete [] edgetab;
  return partToVert;
#else
  std::cerr << "USE_SCOTCH is not defined here in Connectivity::SCOTCH_graphPart\n";
  exit(-1);
#endif
}
