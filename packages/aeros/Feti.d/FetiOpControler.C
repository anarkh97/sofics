#include <Feti.d/Feti.h>
#include <Feti.d/CoarseSet.h>

template<class Scalar>
int
GenFetiOpControler<Scalar>::neighbSubId(int sub, int isub)
{
 return cset[sub].neighbs[isub];
}

template<class Scalar>
int
GenFetiOpControler<Scalar>::numNeighbSubs(int sub)
{
 return cset[sub].numNeighb;
}

template<class Scalar>
int
GenFetiOpControler<Scalar>::index(int sub, int subJ)
{
 int i;
 for(i=0; i < cset[sub].numNeighb; ++i)
   if(subJ == cset[sub].neighbs[i]) return i;
 return -1;
}

