#include <Parser.d/AuxDefs.h>
#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>

#include <map>

GeoSource *geoSource = new GeoSource;

std::map<int,double > weightList;    // allows to change the weight of each class of elements dynamically
std::map<int,double> fieldWeightList;

BCList::BCList(int _loadsetid)
{
 maxbc = 32;
 d     = new BCond[maxbc];
 n     = 0;
 loadsetid = _loadsetid;
}

void
BCList::add(BCond&bc)
{
 if(n == maxbc) {
   int nmaxbc = (3*maxbc)/2;
   BCond *nd = new BCond[nmaxbc];
   int i;
   for(i=0; i < maxbc; ++i)
      nd[i] = d[i];
   delete [] d;
   d = nd;
   maxbc = nmaxbc;
 }
 d[n++] = bc;
}

ComplexBCList::ComplexBCList()
{
 maxbc = 32;
 d     = new ComplexBCond[maxbc];
 n     = 0;
}

void
ComplexBCList::add(ComplexBCond &bc)
{
 if(n==maxbc) {
   int nmaxbc = (3*maxbc)/2;
   ComplexBCond *nd = new ComplexBCond[nmaxbc];
   int i;
   for(i=0; i < maxbc; ++i)
      nd[i] = d[i];
   delete[] d;
   d = nd;
   maxbc=nmaxbc;
 }
 d[n++] = bc;
}
