#include <Math.d/Vector.h>
#include <Element.d/State.h>
#include <Driver.d/Domain.h>
#include <Driver.d/EFrameData.h>

extern Domain *domain;

void
State::getDV(int node, double xyz[3], double v[3], bool noPrescribed)
{
 int loc, loc1;

 if((loc = dsa->locate(node,DofSet::Xdisp)) >= 0) { // free
   xyz[0] = disp[loc];
     v[0] = veloc[loc];
 } else if((loc1 = DSA->locate(node,DofSet::Xdisp)) >= 0 && !noPrescribed) { // prescribed
   xyz[0] = bcx[loc1];
     v[0] = vcx[loc1];
 } else { // not defined
   xyz[0] = 0.0;
     v[0] = 0.0;
 }

 if((loc = dsa->locate(node,DofSet::Ydisp)) >= 0) { // free
   xyz[1] = disp[loc];
     v[1] = veloc[loc];
 } else if((loc1 = DSA->locate(node,DofSet::Ydisp)) >= 0 && !noPrescribed) { // prescribed
   xyz[1] = bcx[loc1];
     v[1] = vcx[loc1];
 } else { // not defined
   xyz[1] = 0.0;
     v[1] = 0.0;
 }

 if((loc = dsa->locate(node,DofSet::Zdisp)) >= 0) { // free
   xyz[2] = disp[loc];
     v[2] = veloc[loc];
 } else if((loc1 = DSA->locate(node,DofSet::Zdisp)) >= 0 && !noPrescribed) { // prescribed
   xyz[2] = bcx[loc1];
     v[2] = vcx[loc1];
 } else { // not defined
   xyz[2] = 0.0;
     v[2] = 0.0;
 }

 if(!domain->solInfo().basicDofCoords) {
   if(cs) {
     if(NFrameData *cd = cs->dofFrame(node)) {
       cd->invTransformVector3(xyz);
       cd->invTransformVector3(v);
     }
   }
   else {
     domain->transformVectorInv(xyz, node, false);
     domain->transformVectorInv(v, node, false);
   }
 }
}

void
State::getTemp(int node, double Temp[1], double dTempdt[1])
{
 int tloc, tloc1;

 if((tloc = dsa->locate(node,DofSet::Temp)) >= 0) { // free
   Temp[0] = disp[tloc];
   dTempdt[0] = veloc[tloc];
 } else if((tloc1 = DSA->locate(node,DofSet::Temp)) >= 0) { // prescribed
   Temp[0] = bcx[tloc1];
   dTempdt[0] = 0.0;
 } else { // not defined
   Temp[0] = 0.0;
   dTempdt[0] = 0.0;
 }
}

void
State::getDVRot(int node, double xyz[6], double v[6])
{
 int loc, loc1;

 if((loc = dsa->locate(node,DofSet::Xdisp)) >= 0) { // free
   xyz[0] = disp[loc];
     v[0] = veloc[loc];
 } else if((loc1 = DSA->locate(node,DofSet::Xdisp)) >= 0) { // prescribed
   xyz[0] = bcx[loc1];
     v[0] = vcx[loc1];
 } else { // not defined
   xyz[0] = 0.0;
     v[0] = 0.0;
 }

 if((loc = dsa->locate(node,DofSet::Ydisp)) >= 0) { // free
   xyz[1] = disp[loc];
     v[1] = veloc[loc];
 } else if((loc1 = DSA->locate(node,DofSet::Ydisp)) >= 0) { // prescribed
   xyz[1] = bcx[loc1];
     v[1] = vcx[loc1];
 } else { // not defined
   xyz[1] = 0.0;
     v[1] = 0.0;
 }

 if((loc = dsa->locate(node,DofSet::Zdisp)) >= 0) { // free
   xyz[2] = disp[loc];
     v[2] = veloc[loc];
 } else if((loc1 = DSA->locate(node,DofSet::Zdisp)) >= 0) { // prescribed
   xyz[2] = bcx[loc1];
     v[2] = vcx[loc1];
 } else { // not defined
   xyz[2] = 0.0;
     v[2] = 0.0;
 }

 if((loc = dsa->locate(node,DofSet::Xrot)) >= 0) { // free
   xyz[3] = disp[loc];
     v[3] = veloc[loc];
 } else if((loc1 = DSA->locate(node,DofSet::Xrot)) >= 0) { // prescribed
   xyz[3] = bcx[loc1];
     v[3] = vcx[loc1];
 } else { // not defined
   xyz[3] = 0.0;
     v[3] = 0.0;
 }

 if((loc = dsa->locate(node,DofSet::Yrot)) >= 0) { // free
   xyz[4] = disp[loc];
     v[4] = veloc[loc];
 } else if((loc1 = DSA->locate(node,DofSet::Yrot)) >= 0) { // prescribed
   xyz[4] = bcx[loc1];
     v[4] = vcx[loc1];
 } else { // not defined
   xyz[4] = 0.0;
     v[4] = 0.0;
 }

 if((loc = dsa->locate(node,DofSet::Zrot)) >= 0) { // free
   xyz[5] = disp[loc];
     v[5] = veloc[loc];
 } else if((loc1 = DSA->locate(node,DofSet::Zrot)) >= 0) { // prescribed
   xyz[5] = bcx[loc1];
     v[5] = vcx[loc1];
 } else { // not defined
   xyz[5] = 0.0;
     v[5] = 0.0;
 }

 if(!domain->solInfo().basicDofCoords) {
   if(cs) {
     if(NFrameData *cd = cs->dofFrame(node)) {
       cd->invTransformVector6(xyz);
       cd->invTransformVector6(v);
     }
   }
   else {
     domain->transformVectorInv(xyz, node, true);
     domain->transformVectorInv(v, node, true);
   }
 }
}
