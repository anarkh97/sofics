#include <cstdio>
#include <Element.d/Element.h>
#include <Element.d/Beam.d/EulerBeam.h>
//#include <Element.d/FourNodeQuad.h>
//#include <Element.d/TimoshenkoBeam.h>
//#include <Element.d/Triangle3.h>
#include <Element.d/Penta.d/Pentahedral.h>
#include <Element.d/Tetra.d/Tetrahedral.h>
#include <Element.d/Brick.d/EightNodeBrick.h>
#include <Element.d/Shell.d/ThreeNodeShell.h>
//#include <Element.d/Shell.d/FourNodeShell.h>
#include <Element.d/Truss.d/TwoNodeTruss.h>
#include <Driver.d/MultiFront.h>

/*
PrioInfo
EulerBeam::examine(int sub, MultiFront *mf)
{
 int wn1 = mf->weight(sub, nn[0]);
 int wn2 = mf->weight(sub, nn[1]);
 int cn1 = mf->weight(nn[0]);
 int cn2 = mf->weight(nn[1]);
 int wr1 = mf->rotWeight(sub, nn[0]);
 int wr2 = mf->rotWeight(sub, nn[1]);

 PrioInfo res;
 res.isReady = (wr1 > 0) || (wr2 > 0) ;
 res.priority = -110 + (wn1-cn1-1) + wn2-cn2-1;
 return res;
}
*/

PrioInfo
ThreeNodeShell::examine(int sub, MultiFront *mf)
{
 int wn[3];
 int rw[3];
 int cn[3];
 int i;
 for(i = 0; i < 3; ++i) {
   wn[i] = mf->weight(sub, nn[i]);
   rw[i] = mf->rotWeight(sub, nn[i]);
   cn[i] = mf->weight(nn[i]);
 }
 int nLast = 0, nTouched = 0, nFirst = 0, nRot = 0;
 for(i = 0; i < 3; ++i) {
/*
   nLast += ((cn[i]-wn[i] == 1) ? 1 : 0);
   nTouched += ((wn[i] > 0) ? 1 : 0);
   nRot += ((rw[i] > 0) ? 1 : 0);
   nFirst += ((wn[i] == 0 && cn[i] != 1) ? 1 : 0);
*/
   if(cn[i]-wn[i] == 1) nLast++;
   if (wn[i] > 0) nTouched++;
   if(rw[i] > 0) nRot++;
   if(wn[i] == 0 && cn[i] != 1) nFirst++;
 }
// fprintf(stderr, "Check: %d %d\n", (nLast-nFirst), nTouched);
 PrioInfo res;
 //res.priority = 0;

 res.isReady = (nRot >= 2) || (nTouched > 2);
 if(res.isReady == false) {
 // fprintf(stderr, ".");
  return res;
 }

 if(nTouched == 3) {
   res.priority = -100;
   return res;
 }
 res.priority = -5*(nLast-nFirst) - 2*nTouched;
 res.isReady = true;
 return res;
}


/*
PrioInfo
FourNodeShell::examine(int sub, MultiFront *mf)
{
 int wn[4];
 int rw[4];
 int cn[4];
 int i;
 for(i = 0; i < 4; ++i) {
   wn[i] = mf->weight(sub, nn[i]);
   rw[i] = mf->rotWeight(sub, nn[i]);
   cn[i] = mf->weight(nn[i]);
 }
 int nLast = 0, nTouched = 0, nFirst = 0, nRot = 0;
 for(i = 0; i < 4; ++i) {

//   nLast += ((cn[i]-wn[i] == 1) ? 1 : 0);
//   nTouched += ((wn[i] > 0) ? 1 : 0);
//   nRot += ((rw[i] > 0) ? 1 : 0);
//   nFirst += ((wn[i] == 0 && cn[i] != 1) ? 1 : 0);

   if(cn[i]-wn[i] == 1) nLast++;
   if (wn[i] > 0) nTouched++;
   if(rw[i] > 0) nRot++;
   if(wn[i] == 0 && cn[i] != 1) nFirst++;
 }
// fprintf(stderr, "Check: %d %d\n", (nLast-nFirst), nTouched);
 PrioInfo res;
 //res.priority = 0;

 res.isReady = (nRot >= 2) || (nTouched > 2);
 if(res.isReady == false) {
 // fprintf(stderr, ".");
  return res;
 }

 if(nTouched == 4) {
   res.priority = -100;
   return res;
 }
 res.priority = -5*(nLast-nFirst) - 2*nTouched;
 res.isReady = true;
 return res;
}
*/
