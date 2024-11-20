#include <Utils.d/ParallelTask.h>
#include <Utils.d/Connectivity.h>
#include <Utils.d/dofset.h>
#include <Threads.d/PHelper.h>
#include <Threads.d/Paral.h>
#include <Element.d/Element.h>
#include <Driver.d/Domain.h>


void
ParallelTask::findSizes(int iSub, long *sizes, Connectivity *stoe, 
                        Elemset *allElements )
{
 int numele = stoe->num(iSub);

 // create subdomain element set
 Elemset elemset;

 int iele;
 for(iele=0; iele < numele; ++iele)
    elemset.elemadd(iele, (*allElements)[(*stoe)[iSub][iele]]);

 // create element to node connectivity from packed element set
 Connectivity *eton = new Connectivity( &elemset );

 // create node to element connectivity from element to node connectivity
 Connectivity *ntoe = eton->reverse();

 // create element to element connectivity
 Connectivity *nton = ntoe->transcon(eton);

 // renumber the nodes
 // 0 = none
 // 1 = sloan
 // 2 = rcm
 compStruct renumber = nton->renumByComponent(1);

 // create subdomain dof set array
 DofSetArray *equations = new DofSetArray(nton->csize(),elemset,renumber.renum);

 // Find subdomain profile size
 sizes[iSub] = nton->findProfileSize( equations );

 // Clean up memory
 delete equations;
 delete nton;
 delete ntoe;
 delete eton;

}

long
ParallelTask::findProfileSize(Connectivity *stoe,Elemset &allElements,
                              long *sizes)
{
 execParal(numSub, this, &ParallelTask::findSizes, sizes, stoe, &allElements);

 long totSize = 0;

 int iSub;
 for(iSub=0; iSub<numSub; ++iSub)
   totSize += sizes[iSub];

 return totSize;
}

int
ParallelTask::findProfileSize(Connectivity **nodeToNode, 
                              EqNumberer **eqNums, int *sizes)
{

 execParal(numSub, this, &ParallelTask::findSize, nodeToNode, eqNums, sizes);

 int totSize = 0;

 int iSub;
 for(iSub=0; iSub<numSub; ++iSub)
   totSize += sizes[iSub];

 return totSize;

}

void
ParallelTask::findSize(int iSub, Connectivity **nodeToNode, 
                       EqNumberer **eqNums, int *sizes)
{
 sizes[iSub] = nodeToNode[iSub]->findProfileSize(eqNums[iSub]);
}

