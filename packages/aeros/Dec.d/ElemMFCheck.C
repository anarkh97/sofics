#include <cstdio>
#include <Element.d/Element.h>
#include <Element.d/Beam.d/EulerBeam.h>
#include <Element.d/Quad4.d/FourNodeQuad.h>
#include <Element.d/Beam.d/TimoshenkoBeam.h>
#include <Element.d/Triangle3.d/Triangle3.h>
#include <Element.d/Membrane.d/Membrane.h>
#include <Element.d/Membrane.d/FourNodeMembrane.h>
#include <Element.d/Penta.d/Pentahedral.h>
#include <Element.d/Tetra.d/Tetrahedral.h>
#include <Element.d/Tetra10.d/TenNodeTetrahedral.h>
#include <Element.d/Brick.d/EightNodeBrick.h>
#include <Element.d/Shell.d/ThreeNodeShell.h>
#include <Element.d/CompShell.d/Compo3NodeShell.h>
#include <Element.d/CompShell.d/Compo4NodeShell.h>
#include <Element.d/FelippaShell.d/FelippaShell.h>
#include <Element.d/FelippaShell.d/FelippaShellX2.h>
#include <Element.d/Truss.d/TwoNodeTruss.h>
#include <Element.d/Truss.d/TwoNodeTrussF.h>
#include <Driver.d/MultiFront.h>
#include <Element.d/Shell.d/FourNodeShell.h>

#include <Element.d/Spring.d/TorSpring.h>
#include <Element.d/Spring.d/LinSpring.h>
#include <Element.d/Shear.d/ShearPanel.h>
#include <Element.d/Spring.d/TransSprlink.h>
#include <Element.d/Spring.d/RotnSprlink.h>
#include <Element.d/Triangle3.d/ThermTriangle.h>
#include <Element.d/Brick20.d/Brick20.h>
#include <Element.d/Shell.d/ConnectedTri.h>

#include <Element.d/ThermQuad.d/Therm3DQuad.h>
#include <Element.d/Truss.d/Therm2NodeBar.h>
#include <Element.d/ThermQuad.d/ThermQuadGal.h>
#include <Element.d/Convection.d/BarConvec.h>
#include <Element.d/Convection.d/QuadConvec.h>
#include <Element.d/Convection.d/TriangleConvec.h>
#include <Element.d/Radiation.d/BarRadiation.h>
#include <Element.d/Radiation.d/TriangleRadiation.h>
#include <Element.d/Radiation.d/QuadRadiation.h>
#include <Element.d/ContactResistance.d/QuadContact.h>
#include <Element.d/ContactResistance.d/BrickContact.h>
#include <Element.d/ContactResistance.d/PentaContact.h>
#include <Element.d/BulkFluid.d/TriangleBulk.h>
#include <Element.d/BulkFluid.d/TetraBulk.h>
#include <Element.d/BulkFluid.d/PentaBulk.h>
#include <Element.d/Tetra.d/ThermIsoParamTetra.h>
#include <Element.d/Brick.d/ThermBrick.h>
#include <Element.d/Shell.d/Therm4NoShell.h>

#include <Element.d/Helm.d/HelmQuadGal.h>
#include <Element.d/Helm.d/HelmTri3Gal.h>
#include <Element.d/Helm.d/HelmTri3Gls.h>
#include <Element.d/Helm.d/HelmQuadGls.h>
#include <Element.d/Helm.d/TetraHelmGal.h>
#include <Element.d/Helm.d/TetraHelmGLS.h>
#include <Element.d/Helm.d/HelmBrick.h>
#include <Element.d/Helm.d/HelmBrickGLS.h>
#include <Element.d/Helm.d/Tetra10HelmGal.h>
#include <Element.d/Helm.d/HelmPenta.h>
#include <Element.d/Helm.d/HelmIsoParamTetra.h>
#include <Element.d/Helm.d/LEIsoParamTetra.h>
#include <Element.d/Helm.d/HelmIsoParamHexa.h>
#include <Element.d/Helm.d/ThermIsoParamHexa.h>
#include <Element.d/Helm.d/HelmIsoParamQuad.h>
#include <Element.d/Helm.d/HelmIsoParamTri.h>
#include <Element.d/Helm.d/HelmSpectralIsoParamHexa.h>
#include <Element.d/Helm.d/LEIsoParamHexa.h>
#include <Element.d/Helm.d/HelmLagQuadGal.h>
#include <Element.d/Helm.d/HelmSpectralIsoParamQuad.h>

#include <Element.d/Shell.d/Therm3NoShell.h>

#include <Element.d/Brick32.d/Brick32.h>
#include <Element.d/Penta15.d/Penta15.h>
#include <Element.d/Penta26.d/Penta26.h>
#include <Element.d/Helm.d/HelmBrick32.h>
#include <Element.d/Helm.d/HelmPenta26.h>
#include <Element.d/NonLinearity.d/NLMembrane.h>

#include <Element.d/MpcElement.d/MpcElement.h>
#include <Element.d/MpcElement.d/FsiElement.h>

#include <Element.d/FluidQuad.d/SloshQuadGal.h>
#include <Element.d/FluidQuad.d/BarSloshFS.h>
#include <Element.d/FluidTetra.d/SloshTetra.h>
#include <Element.d/FluidTriangle3.d/SloshTriangleFS.h>
#include <Element.d/FluidTetra.d/HEVibTetra.h>

#ifdef USE_EIGEN3
#include <Element.d/Rigid.d/RigidTwoNodeTruss.h>
#include <Element.d/Rigid.d/RigidBeam.h>
#include <Element.d/Rigid.d/RigidSpring.h>
#include <Element.d/Rigid.d/RigidTransSprlink.h>
#include <Element.d/Rigid.d/RigidRotnSprlink.h>
#include <Element.d/Rigid.d/RigidThreeNodeShell.h>
#include <Element.d/Rigid.d/RigidEightNodeBrick.h>
#include <Element.d/Rigid.d/RigidSolid.h>
#include <Element.d/Rigid.d/RigidSolid6Dof.h>
#include <Element.d/Rigid.d/RigidFourNodeShell.h>

#include <Element.d/Joint.d/WeldedJoint.h>
#include <Element.d/Joint.d/SphericalJoint.h>
#include <Element.d/Joint.d/RevoluteJoint.h>
#include <Element.d/Joint.d/TranslationalJoint.h>
#include <Element.d/Joint.d/UniversalJoint.h>
#include <Element.d/Joint.d/CylindricalJoint.h>
#include <Element.d/Joint.d/PrismaticJoint.h>
#include <Element.d/Joint.d/PinInSlotJoint.h>
#include <Element.d/Joint.d/PlanarJoint.h>

#include <Element.d/Joint.d/RevoluteDriver.h>
#include <Element.d/Joint.d/PrismaticDriver.h>
#include <Element.d/Joint.d/RevoluteActuator.h>
#include <Element.d/Joint.d/PrismaticActuator.h>

#include <Element.d/Joint.d/SphericalJointSpringCombo.h>
#include <Element.d/Joint.d/TranslationalJointSpringCombo.h>
#include <Element.d/Joint.d/UniversalJointSpringCombo.h>
#include <Element.d/Joint.d/RevoluteJointSpringCombo.h>
#include <Element.d/Joint.d/CylindricalJointSpringCombo.h>
#include <Element.d/Joint.d/PrismaticJointSpringCombo.h>
#include <Element.d/Joint.d/PinInSlotJointSpringCombo.h>
#include <Element.d/Joint.d/RevoluteJointSpringComboWithFreeplay.h>
#include <Element.d/Joint.d/PrismaticJointSpringComboWithFreeplay.h>
#endif

#include <Element.d/BelytschkoTsayShell.d/BelytschkoTsayShell.h>

extern bool allowMechanisms;

PrioInfo examineBeam2(int sub, MultiFront *mf, int *nn);
PrioInfo examineTri3Shell(int sub, MultiFront *mf, int *nn);
PrioInfo examineQuad4Shell(int sub, MultiFront *mf, int *nn);

PrioInfo
examineBar2(int sub, MultiFront *mf, int *nn)
{
 if (allowMechanisms) return examineBeam2(sub, mf, nn);

 int wn1 = mf->weight(sub, nn[0]);
 int wn2 = mf->weight(sub, nn[1]);
 int cn1 = mf->weight(nn[0]);
 int cn2 = mf->weight(nn[1]);

 PrioInfo res;
 res.isReady = wn1 > 0 && wn2 > 0;
 if(res.isReady == false) return res;
 res.priority = -110 + (wn1-cn1-1) + wn2-cn2-1;
 // Compiler bug workaround
 res.isReady = true;
 return res;
}

PrioInfo
TwoNodeTruss::examine(int sub, MultiFront *mf)
{
  return examineBar2(sub, mf, nn);
}

PrioInfo
TwoNodeTrussF::examine(int sub, MultiFront *mf)
{
  return examineBar2(sub, mf, nn);
}

PrioInfo
TransSprlink::examine(int sub, MultiFront *mf)
{
  return examineBar2(sub, mf, nn);
}

PrioInfo
RotnSprlink::examine(int sub, MultiFront *mf)
{
  return examineBar2(sub, mf, nn);
}

PrioInfo
BarConvec::examine(int sub, MultiFront *mf)
{
  return examineBar2(sub, mf, nn);
}

PrioInfo
BarRadiation::examine(int sub, MultiFront *mf)
{
  return examineBar2(sub, mf, nn);
}

PrioInfo
Therm2NodeBar::examine(int sub, MultiFront *mf)
{
  return examineBar2(sub, mf, nn);
}

// ADDED FOR SLOSHING PROBLEM, EC, 20070713
/*
PrioInfo
BarSloshFS::examine(int sub, MultiFront *mf)
{
  return examineBar2(sub, mf, nn);
}
*/

PrioInfo
examineBeam2(int sub, MultiFront *mf, int *nn)
{
 int wn1 = mf->weight(sub, nn[0]);
 int wn2 = mf->weight(sub, nn[1]);
 int cn1 = mf->weight(nn[0]);
 int cn2 = mf->weight(nn[1]);
 int wr1 = mf->rotWeight(sub, nn[0]);
 int wr2 = mf->rotWeight(sub, nn[1]);

 PrioInfo res;
 res.isReady = (wr1 > 0) || (wr2 > 0) ;
 if(res.isReady == false) return res;
 res.priority = -110 + (wn1-cn1-1) + wn2-cn2-1;
 // Compiler bug workaround
 res.isReady = true;
 return res;
}

PrioInfo
TimoshenkoBeam::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
EulerBeam::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo 
TorSpring::examine(int sub, MultiFront *mf)
{
 int wn1 = mf->weight(sub, nn[0]);
 int cn1 = mf->weight(nn[0]);

 PrioInfo res;
 res.isReady = wn1 > 0;
 if(res.isReady == false) return res;
 res.priority = -110 + (wn1-cn1-1);
 // Compiler bug workaround
 res.isReady = true;
 return res;
}

PrioInfo
LinSpring::examine(int sub, MultiFront *mf)
{
 int wn1 = mf->weight(sub, nn[0]);
 int cn1 = mf->weight(nn[0]);

 PrioInfo res;
 res.isReady = wn1 > 0;
 if(res.isReady == false) return res;
 res.priority = -110 + (wn1-cn1-1);
 // Compiler bug workaround
 res.isReady = true;
 return res;
}

PrioInfo
examineQuad4(int sub, MultiFront *mf, int *nn)
{
 if(allowMechanisms) return examineQuad4Shell(sub, mf, nn);

 int wn[4];
 int cn[4];
 int i;
 for(i = 0; i < 4; ++i) {
   wn[i] = mf->weight(sub, nn[i]);
   cn[i] = mf->weight(nn[i]);
 }
 PrioInfo res;
 // this ready rule is rather restrictive.. We could simply accept
 // two nodes as being in the subdomain so far;
 res.isReady = (wn[2] > 0 || wn[0] > 0) && (wn[1] > 0 || wn[3] > 0);
 if(!res.isReady) return res;

 if(wn[2] > 0 && wn[0] > 0 && wn[1] > 0 && wn[3] > 0) {
   res.priority = -100;
   return res;
 }
 int nLast = 0, nTouched = 0, nFirst = 0;
 for(i = 0; i < 4; ++i) {
   nLast += (cn[i]-wn[i] == 1) ? 1 : 0;
   nTouched += (wn[i] > 0) ? 1 : 0;
   nFirst += (wn[i] == 0) ? 1 : 0;
 }
 res.priority = -5*(nLast-nFirst) - 2*nTouched;
 // Compiler bug workaround
 res.isReady = true;
 return res;
}

PrioInfo
MpcElement::examine(int sub, MultiFront *mf)
{
 // PJSA 7-31-06
 int count = 0;
 for(int i = 0; i < nNodes; ++i) {
   int wn = mf->weight(sub, nn[i]);
   if(wn > 0) count++;
 }
 PrioInfo res;
 res.isReady = (count > 0);

 if(res.isReady) res.priority = -100;
 return res;
}

PrioInfo
FsiElement::examine(int sub, MultiFront *mf)
{
 // this needs to be better
 int scount = 0, fcount = 0;
 for(int i = 0; i < nnodes; ++i) {
   int wn = mf->weight(sub, nn[i]);
   if((i < nnodes-1) && (wn > 0)) scount++;
   if((i == nnodes-1) && (wn > 0)) fcount++;
 }
 PrioInfo res;
 res.isReady = (fcount || scount);

 // PJSA: seems to minimize unsafe nodes in corner maker by always using maximum priority
 if(res.isReady) res.priority = -100; 
 return res;
}

PrioInfo
examineHex8(int sub, MultiFront *mf, int *nn)
{
 int wn[8];
 int cn[8];
 int i;
 for(i = 0; i < 8; ++i) {
   wn[i] = mf->weight(sub, nn[i]);
   cn[i] = mf->weight(nn[i]);
 }
 PrioInfo res;
 // this ready rule is rather restrictive.. We could simply accept
 // Four nodes as being in the subdomain so far;
 res.isReady =     ((wn[1] && wn[0]) &&
                   ((wn[2] && wn[3]) || (wn[4] && wn[5])))
                || ((wn[5] && wn[6]) &&
                   ((wn[1] && wn[2]) || (wn[7] && wn[4])))
                || ((wn[3] && wn[7]) &&
                   ((wn[0] && wn[4]) || (wn[2] && wn[6])));
 if(!res.isReady) return res;

 if(wn[2] > 0 && wn[0] > 0 && wn[1] > 0 && wn[3] > 0) {
   res.priority = -100;
   return res;
 }
 int nLast = 0, nTouched = 0, nFirst = 0;
 for(i = 0; i < 8; ++i) {
   nLast += (cn[i]-wn[i] == 1) ? 1 : 0;
   nTouched += (wn[i] > 0) ? 1 : 0;
   nFirst += (wn[i] == 0) ? 1 : 0;
 }
 res.priority = -5*(nLast-nFirst) - 2*nTouched;
 // Compiler bug workaround
 res.isReady = true;
 return res;
}


PrioInfo
FourNodeQuad::examine(int sub, MultiFront *mf)
{
  return examineQuad4(sub, mf, nn);
}

PrioInfo
NLMembrane4::examine(int sub, MultiFront *mf)
{
  return examineQuad4(sub, mf, nn);
}

PrioInfo
BelytschkoTsayShell::examine(int sub, MultiFront *mf)
{
  return examineQuad4(sub, mf, nn);
}

/*
PrioInfo Therm4NoShell::examine(int sub, MultiFront *mf)
{
  return examineQuad4(sub, mf, nn);
}
*/

PrioInfo
Therm3DQuad::examine(int sub, MultiFront *mf)
{
  return examineQuad4(sub, mf, nn);
}

PrioInfo
ThermQuadGal::examine(int sub, MultiFront *mf)
{
  return examineQuad4(sub, mf, nn);
}

PrioInfo
QuadConvec::examine(int sub, MultiFront *mf)
{
  return examineQuad4(sub, mf, nn);
}

PrioInfo
QuadRadiation::examine(int sub, MultiFront *mf)
{
  return examineQuad4(sub, mf, nn);
}

PrioInfo
QuadContact::examine(int sub, MultiFront *mf)
{
  return examineQuad4(sub, mf, nn);
}

PrioInfo
HelmQuadGal::examine(int sub, MultiFront *mf)
{
  return examineQuad4(sub, mf, nn);
}

PrioInfo
HelmQuadGls::examine(int sub, MultiFront *mf)
{
  return examineQuad4(sub, mf, nn);
}

// ADDED FOR SLOSHING PROBLEM, EC, 20070713
/*
PrioInfo
SloshQuadGal::examine(int sub, MultiFront *mf)
{
  return examineQuad4(sub, mf, nn);
}
*/

PrioInfo
HelmIsoParamQuad::examine(int sub, MultiFront *mf)
{
  int nSummit[4];

  nSummit[0] = nn[0];
  nSummit[1] = nn[order-1];
  nSummit[2] = nn[order*order-1];
  nSummit[3] = nn[order*order-order];
  return examineQuad4(sub, mf, nSummit);
}


PrioInfo
examineTri3(int sub, MultiFront *mf, int *nn)
{
 if(allowMechanisms) return examineTri3Shell(sub, mf, nn);

 int wn[3];
 int cn[3];
 int i;

 //bool monitor = nn[0] == 22763 && nn[1] == 28696 && nn[2] == 28698;
 for(i = 0; i < 3; ++i) {
   wn[i] = mf->weight(sub, nn[i]);
   cn[i] = mf->weight(nn[i]);
 }
 int nLast = 0, nTouched = 0, nFirst = 0, nRotTouched = 0;
 for(i = 0; i < 3; ++i) {
   nLast += (cn[i]-wn[i] == 1) ? 1 : 0;
   nTouched += (wn[i] > 0) ? 1 : 0;
   nFirst += (wn[i] == 0 && cn[i] != 1) ? 1 : 0;
   nRotTouched += (mf->rotWeight(sub, nn[i]) > 0) ? 1 : 0;
 }

 //if(monitor) fprintf(stderr, "nTouched = %d nFirst  = %d nLast = %d\n", nTouched, nFirst, nLast);
 PrioInfo res;
 res.isReady = (nRotTouched > 1) || (nTouched == 3);
 if(!res.isReady) return res;

 if(nTouched == 3) {
   res.priority = -100;
   return res;
 }
 res.priority = -5*(nLast-nFirst) - 2*nTouched;
 // Compiler bug workaround
 res.isReady = true;
 return res;
}

PrioInfo
Triangle3::examine(int sub, MultiFront *mf)
{
  return examineTri3(sub, mf, nn);
}

PrioInfo
NLMembrane::examine(int sub, MultiFront *mf)
{
  return examineTri3(sub, mf, n);
}

PrioInfo
Membrane::examine(int sub, MultiFront *mf)
{
  return examineTri3(sub, mf, nn);
}

PrioInfo
FourNodeMembrane::examine(int sub, MultiFront *mf)
{
  return examineQuad4(sub, mf, nn);
}

PrioInfo 
TriangleConvec::examine(int sub, MultiFront *mf)
{
  return examineTri3(sub, mf, nn);
}

PrioInfo
TriangleRadiation::examine(int sub, MultiFront *mf)
{
  return examineTri3(sub, mf, nn);
}

PrioInfo
TriangleBulk::examine(int sub, MultiFront *mf)
{
  return examineTri3(sub, mf, nn);
}

PrioInfo
HelmTri3Gal::examine(int sub, MultiFront *mf)
{
  return examineTri3(sub, mf, nn);
}

PrioInfo
HelmTri3Gls::examine(int sub, MultiFront *mf)
{
  return examineTri3(sub, mf, nn);
}

/*
PrioInfo
Therm3NoShell::examine(int sub, MultiFront *mf)
{
  return examineTri3(sub, mf, nn);
}
*/


PrioInfo
ThermTriangle::examine(int sub, MultiFront *mf)
{
  return examineTri3(sub, mf, nn);
}

// ADDED FOR SLOSHING ELEMENT, EC, 20070713
/*
PrioInfo
ThermTriangle::examine(int sub, MultiFront *mf)
{
  return examineTri3(sub, mf, nn);
}
*/


PrioInfo
HelmIsoParamTri::examine(int sub, MultiFront *mf)
{
  int nSummit[3];

  nSummit[0] = nn[0];
  nSummit[1] = nn[order-1];
  nSummit[2] = nn[(order*(order+1))/2-1];
  return examineTri3(sub, mf, nSummit);
}

PrioInfo 
examineTri6(int sub, MultiFront *mf, int *nn)
{
 // NEW & UNTESTED
 int wn[6];
 int cn[6];
 int i;

 //bool monitor = nn[0] == 22763 && nn[1] == 28696 && nn[2] == 28698;
 for(i = 0; i < 6; ++i) {
   wn[i] = mf->weight(sub, nn[i]);
   cn[i] = mf->weight(nn[i]);
 }
 int nLast = 0, nTouched = 0, nFirst = 0;
 for(i = 0; i < 6; ++i) {
   nLast += (cn[i]-wn[i] == 1) ? 1 : 0;
   nTouched += (wn[i] > 0) ? 1 : 0;
   nFirst += (wn[i] == 0 && cn[i] != 1) ? 1 : 0;
   //nRotTouched += (mf->rotWeight(sub, nn[i]) > 0) ? 1 : 0;
 }
 //if(monitor) fprintf(stderr, "nTouched = %d nFirst  = %d nLast = %d\n", nTouched, nFirst, nLast);
 PrioInfo res;
 res.isReady = (nTouched > 0);
 if(!res.isReady) return res;

 if(nTouched > 0) {
   res.priority = -100;
   return res;
 }
 res.priority = -5*(nLast-nFirst) - 2*nTouched;
 // Compiler bug workaround
 res.isReady = true;
 return res;
}

PrioInfo 
examineQuad8(int sub, MultiFront *mf, int *nn)
{
 // NEW & UNTESTED
 int wn[8];
 int cn[8];
 int i;

 //bool monitor = nn[0] == 22763 && nn[1] == 28696 && nn[2] == 28698;
 for(i = 0; i < 8; ++i) {
   wn[i] = mf->weight(sub, nn[i]);
   cn[i] = mf->weight(nn[i]);
 }
 int nLast = 0, nTouched = 0, nFirst = 0;
 for(i = 0; i < 8; ++i) {
   nLast += (cn[i]-wn[i] == 1) ? 1 : 0;
   nTouched += (wn[i] > 0) ? 1 : 0;
   nFirst += (wn[i] == 0 && cn[i] != 1) ? 1 : 0;
   //nRotTouched += (mf->rotWeight(sub, nn[i]) > 0) ? 1 : 0;
 }
 //if(monitor) fprintf(stderr, "nTouched = %d nFirst  = %d nLast = %d\n", nTouched, nFirst, nLast);
 PrioInfo res;
 res.isReady = (nTouched > 0);
 if(!res.isReady) return res;

 if(nTouched > 0) {
   res.priority = -100;
   return res;
 }
 res.priority = -5*(nLast-nFirst) - 2*nTouched;
 // Compiler bug workaround
 res.isReady = true;
 return res;
}

PrioInfo 
examineTet4(int sub, MultiFront *mf, int *nn)
{
 int wn[4];
 int cn[4];
 int i;
 for(i = 0; i < 4; ++i) {
   wn[i] = mf->weight(sub, nn[i]);
   cn[i] = mf->weight(nn[i]);
 }
 int nLast = 0, nTouched = 0, nFirst = 0;
 for(i = 0; i < 4; ++i) {
   nLast += (cn[i]-wn[i] == 1) ? 1 : 0;
   nTouched += (wn[i] > 0) ? 1 : 0;
   nFirst += (wn[i] == 0) ? 1 : 0;
 }
 PrioInfo res;
 res.isReady = nTouched >= 3;
 if(!res.isReady) return res;

 if(wn[2] > 0 && wn[0] > 0 && wn[1] > 0 && wn[3] > 0) {
   res.priority = -100;
   return res;
 }
 res.priority = -5*(nLast-nFirst) - 2*nTouched;
 // Compiler bug workaround
 res.isReady = true;
 return res;
}

PrioInfo
Tetrahedral::examine(int sub, MultiFront *mf)
{
 return examineTet4(sub, mf, nn);
}

PrioInfo
TetraBulk::examine(int sub, MultiFront *mf)
{
 return examineTet4(sub, mf, nn);
}

PrioInfo 
TetraHelmGal::examine(int sub, MultiFront *mf)
{
  return examineTet4(sub, mf, nn);
}

PrioInfo
TetraHelmGLS::examine(int sub, MultiFront *mf)
{
  return examineTet4(sub, mf, nn);
}

PrioInfo
ShearPanel::examine(int sub, MultiFront *mf)
{
  return examineTet4(sub, mf, nn);
}

// ADDED FOR SLOSHING PROBLEM, EC, 20070713
/*
PrioInfo
SloshTetra::examine(int sub, MultiFront *mf)
{
  return examineTet4(sub, mf, nn);
}
*/

// ADDED FOR HEV PROBLEM, EC, 20070815
/*
PrioInfo
HEVibTetra::examine(int sub, MultiFront *mf)
{
  return examineTet4(sub, mf, nn);
}
*/

PrioInfo
examineHex20(int sub, MultiFront *mf, int *nn)
{
 // NEW & UNTESTED
 int wn[20];
 int cn[20];
 int i;

 //bool monitor = nn[0] == 22763 && nn[1] == 28696 && nn[2] == 28698;
 for(i = 0; i < 20; ++i) {
   wn[i] = mf->weight(sub, nn[i]);
   cn[i] = mf->weight(nn[i]);
 }
 int nLast = 0, nTouched = 0, nFirst = 0;
 for(i = 0; i < 20; ++i) {
   nLast += (cn[i]-wn[i] == 1) ? 1 : 0;
   nTouched += (wn[i] > 0) ? 1 : 0;
   nFirst += (wn[i] == 0 && cn[i] != 1) ? 1 : 0;
   //nRotTouched += (mf->rotWeight(sub, nn[i]) > 0) ? 1 : 0;
 }
 //if(monitor) fprintf(stderr, "nTouched = %d nFirst  = %d nLast = %d\n", nTouched, nFirst, nLast);
 PrioInfo res;
 res.isReady =  ((wn[1] && wn[0]) && 
		   ((wn[2] && wn[3]) || (wn[4] && wn[5])))
		|| ((wn[5] && wn[6]) && 
		   ((wn[1] && wn[2]) || (wn[7] && wn[4])))
		|| ((wn[3] && wn[7]) && 
		    ((wn[0] && wn[4]) || (wn[2] && wn[6])));//(nTouched > 3);
 if(!res.isReady) return res;

 if(wn[2] > 0 && wn[0] > 0 && wn[1] > 0 && wn[3] > 0)
   {
   res.priority = -100;
   return res;
 }
 res.priority = -5*(nLast-nFirst) - 2*nTouched;
 // Compiler bug workaround
 res.isReady = true;
 return res;
}

PrioInfo 
Brick20::examine(int sub, MultiFront *mf)
{ 
  return examineHex20(sub, mf, nn);
}

PrioInfo
examineTri3Shell(int sub, MultiFront *mf, int *nn)
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

PrioInfo
ThreeNodeShell::examine(int sub, MultiFront *mf)
{
  return examineTri3Shell(sub, mf, nn);
}

PrioInfo
Compo3NodeShell::examine(int sub, MultiFront *mf)
{
  return examineTri3Shell(sub, mf, nn);
}

#ifdef USE_EIGEN3
PrioInfo
FelippaShell::examine(int sub, MultiFront *mf)
{
  return examineTri3Shell(sub, mf, nn);
}
#endif

PrioInfo
Therm3NoShell::examine(int sub, MultiFront *mf)
{
  return examineTri3Shell(sub, mf, nn);
}


PrioInfo
examineQuad4Shell(int sub, MultiFront *mf, int *nn)
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

 if(nTouched == 4) {
   res.priority = -100;
   return res;
 }
 res.priority = -5*(nLast-nFirst) - 2*nTouched;
 res.isReady = true;
 return res;
}

PrioInfo
FourNodeShell::examine(int sub, MultiFront *mf)
{
  return examineQuad4Shell(sub, mf, nn);
}

PrioInfo
Compo4NodeShell::examine(int sub, MultiFront *mf)
{
  return examineQuad4Shell(sub, mf, nn);
}

#ifdef USE_EIGEN3
PrioInfo
FelippaShellX2::examine(int sub, MultiFront *mf)
{
  return examineQuad4Shell(sub, mf, nn);
}
#endif

PrioInfo Therm4NoShell::examine(int sub, MultiFront *mf)
{
  return examineQuad4Shell(sub, mf, nn);
}

PrioInfo
ConnectedTri::examine(int sub, MultiFront *mf)
{
  return examineQuad4Shell(sub, mf, nn);
}

PrioInfo
examinePenta6(int sub, MultiFront *mf, int *nn)
{
 int wn[6];
 int cn[6];
 int i;
 for(i = 0; i < 6; ++i) {
   wn[i] = mf->weight(sub, nn[i]);
   cn[i] = mf->weight(nn[i]);
 }
 PrioInfo res;
 // this ready rule is rather restrictive.. We could simply accept
 // three nodes as being in the subdomain so far;
 res.isReady = ((wn[0] && wn[1]) && (wn[2] || (wn[3] && wn[4]))) ||
                 ((wn[4] && wn[5]) && (wn[3] || (wn[1] && wn[2]))) ||
                 (wn[5] && wn[3] && wn[2] && wn[0]);
/*
 int cnt = 0;
  for(i = 0; i < 6; ++i)
  if(wn[i] > 0)
    cnt++;
  res.isReady = cnt >= 3;
*/
 if(!res.isReady) return res;

 int nLast = 0, nTouched = 0, nFirst = 0;
 for(i = 0; i < 6; ++i) {
   nLast += (cn[i]-wn[i] == 1) ? 1 : 0;
   nTouched += (wn[i] > 0) ? 1 : 0;
   nFirst += (wn[i] == 0) ? 1 : 0;
 }
 if(nTouched == 6) {
   res.priority = -100;
   // Compiler bug workaround
   res.isReady = true;
   return res;
 }
 res.priority = -5*(nLast-nFirst) - 2*nTouched;
 // Compiler bug workaround
 res.isReady = true;
 return res;
}

PrioInfo
Pentahedral::examine(int sub, MultiFront *mf)
{
  return examinePenta6(sub, mf, nn);
}

PrioInfo
PentaContact::examine(int sub, MultiFront *mf)
{
  return examinePenta6(sub, mf, nn);
}

PrioInfo
HelmPenta::examine(int sub, MultiFront *mf)
{
  return examinePenta6(sub, mf, nn);
}

PrioInfo
Penta15::examine(int sub, MultiFront *mf)
{
  return examinePenta6(sub, mf, nn);
}

PrioInfo
Penta26::examine(int sub, MultiFront *mf)
{
  return examinePenta6(sub, mf, nn);
}

PrioInfo
HelmPenta26::examine(int sub, MultiFront *mf)
{
  return examinePenta6(sub, mf, nn);
}

PrioInfo
examineTet10(int sub, MultiFront *mf, int *nn)
{
 int wn[4];
 int cn[4];
 int i;
 for(i = 0; i < 4; ++i) {
   wn[i] = mf->weight(sub, nn[i]);
   cn[i] = mf->weight(nn[i]);
 }
 int nLast = 0, nTouched = 0, nFirst = 0;
 for(i = 0; i < 4; ++i) {
   nLast += (cn[i]-wn[i] == 1) ? 1 : 0;
   nTouched += (wn[i] > 0) ? 1 : 0;
   nFirst += (wn[i] == 0) ? 1 : 0;
 }
 PrioInfo res;
 res.isReady = nTouched >= 3;
 if(!res.isReady) return res;

 if(wn[2] > 0 && wn[0] > 0 && wn[1] > 0 && wn[3] > 0) {
   res.priority = -100;
   return res;
 }
 res.priority = -5*(nLast-nFirst) - 2*nTouched;
 // Compiler bug workaround
 res.isReady = true;
 return res;
}

PrioInfo
examinePyramid5(int sub, MultiFront *mf, int *nn)
{
 int wn[5];
 int cn[5];
 int i;
 for(i = 0; i < 5; ++i) {
   wn[i] = mf->weight(sub, nn[i]);
   cn[i] = mf->weight(nn[i]);
 }
 int nLast = 0, nTouched = 0, nFirst = 0;
 for(i = 0; i < 5; ++i) {
   nLast += (cn[i]-wn[i] == 1) ? 1 : 0;
   nTouched += (wn[i] > 0) ? 1 : 0;
   nFirst += (wn[i] == 0) ? 1 : 0;
 }
 PrioInfo res;
 res.isReady = (nTouched >= 3);
 if(res.isReady) {
   res.priority = (nTouched == 5) ? -100 : 5*(nLast-nFirst) - 2*nTouched;
 }
 return res;
}

PrioInfo
TenNodeTetrahedral::examine(int sub, MultiFront *mf)
{
  return examineTet10(sub, mf, nn);
}

PrioInfo
Tetra10HelmGal::examine(int sub, MultiFront *mf)
{
  return examineTet10(sub, mf, nn);
}

PrioInfo
EightNodeBrick::examine(int sub, MultiFront *mf)
{
  return examineHex8(sub, mf, nn);
}

PrioInfo
ThermBrick::examine(int sub, MultiFront *mf)
{
  return examineHex8(sub, mf, nn);
}

PrioInfo
BrickContact::examine(int sub, MultiFront *mf)
{
  return examineHex8(sub, mf, nn);
}

PrioInfo
PentaBulk::examine(int sub, MultiFront *mf)
{
  return examinePyramid5(sub, mf, nn);
}

PrioInfo
HelmBrick::examine(int sub, MultiFront *mf)
{
  return examineHex8(sub, mf, nn);
}

PrioInfo
HelmBrickGLS::examine(int sub, MultiFront *mf)
{
  return examineHex8(sub, mf, nn);
}

PrioInfo
Brick32::examine(int sub, MultiFront *mf)
{
  return examineHex8(sub, mf, nn);
}

PrioInfo
HelmBrick32::examine(int sub, MultiFront *mf)
{
  return examineHex8(sub, mf, nn);
}

PrioInfo
HelmIsoParamHexa::examine(int sub, MultiFront *mf)
{
  //find the 8 nodes at the summits of the Hexaedron
  int nSummit[8];

  nSummit[0] = nn[0];
  nSummit[1] = nn[order-1];
  nSummit[2] = nn[order*order-1];
  nSummit[3] = nn[order*order-order];

  nSummit[4] = nn[order*order*order-order*order];
  nSummit[5] = nn[order*order*order-order*order+order-1];
  nSummit[6] = nn[order*order*order-1];
  nSummit[7] = nn[order*order*order-order];

  // this is just copied from Brick32
  return examineHex8(sub, mf, nSummit);
}

PrioInfo
ThermIsoParamHexa::examine(int sub, MultiFront *mf)
{
  //find the 8 nodes at the summits of the Hexaedron
  int nSummit[8];

  nSummit[0] = nn[0];
  nSummit[1] = nn[order-1];
  nSummit[2] = nn[order*order-1];
  nSummit[3] = nn[order*order-order];

  nSummit[4] = nn[order*order*order-order*order];
  nSummit[5] = nn[order*order*order-order*order+order-1];
  nSummit[6] = nn[order*order*order-1];
  nSummit[7] = nn[order*order*order-order];

  // this is just copied from Brick32
  return examineHex8(sub, mf, nSummit);
}


PrioInfo
HelmSpectralIsoParamHexa::examine(int sub, MultiFront *mf)
{
  //find the 8 nodes at the summits of the Hexaedron
  int nSummit[8];

  nSummit[0] = nn[0];
  nSummit[1] = nn[order-1];
  nSummit[2] = nn[order*order-1];
  nSummit[3] = nn[order*order-order];

  nSummit[4] = nn[order*order*order-order*order];
  nSummit[5] = nn[order*order*order-order*order+order-1];
  nSummit[6] = nn[order*order*order-1];
  nSummit[7] = nn[order*order*order-order];

  // this is just copied from Brick32
  return examineHex8(sub, mf, nSummit);
}


PrioInfo
LEIsoParamHexa::examine(int sub, MultiFront *mf)
{
  //find the 8 nodes at the summits of the Hexaedron
  int nSummit[8];

  nSummit[0] = nn[0];
  nSummit[1] = nn[order-1];
  nSummit[2] = nn[order*order-1];
  nSummit[3] = nn[order*order-order];

  nSummit[4] = nn[order*order*order-order*order];
  nSummit[5] = nn[order*order*order-order*order+order-1];
  nSummit[6] = nn[order*order*order-1];
  nSummit[7] = nn[order*order*order-order];

  // this is just copied from Brick32
  return examineHex8(sub, mf, nSummit);
}

PrioInfo
HelmIsoParamTetra::examine(int sub, MultiFront *mf)
{
  int nSummit[4];
  nSummit[0] = nn[0];
  nSummit[1] = nn[order-1];
  nSummit[2] = nn[order*(order+1)/2-1];
  nSummit[3] = nn[(((order+3)*order+2)*order-6)/6];
  return examineTet4(sub, mf, nSummit);
}


PrioInfo
ThermIsoParamTetra::examine(int sub, MultiFront *mf)
{
  int nSummit[4];
  nSummit[0] = nn[0];
  nSummit[1] = nn[order-1];
  nSummit[2] = nn[order*(order+1)/2-1];
  nSummit[3] = nn[(((order+3)*order+2)*order-6)/6];
  return examineTet4(sub, mf, nSummit);
}


PrioInfo
LEIsoParamTetra::examine(int sub, MultiFront *mf)
{
  int nSummit[4];
  nSummit[0] = nn[0];
  nSummit[1] = nn[order-1];
  nSummit[2] = nn[order*(order+1)/2-1];
  return examineTet4(sub, mf, nSummit);
}


PrioInfo
HelmLagQuadGal::examine(int sub, MultiFront *mf)
{
  int nSummit[4];
  nSummit[0] = nn[0];
  nSummit[1] = nn[order-1];
  nSummit[2] = nn[order*(order)-1];
  nSummit[3] = nn[order*(order)-(order)];
  return examineHex8(sub, mf, nSummit);
}

PrioInfo
HelmSpectralIsoParamQuad::examine(int sub, MultiFront *mf)
{
  return examineQuad4(sub, mf, nn);
}

#ifdef USE_EIGEN3
PrioInfo
RigidTwoNodeTruss::examine(int sub, MultiFront *mf)
{
  return examineBar2(sub, mf, nn);
}

PrioInfo
RigidTwoNodeTrussWithMass::examine(int sub, MultiFront *mf)
{
  return examineBar2(sub, mf, nn);
}

PrioInfo
RigidBeam::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
RigidBeamWithMass::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
RigidSpring::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
RigidTransSprlink::examine(int sub, MultiFront *mf)
{
  return examineBar2(sub, mf, nn);
}

PrioInfo
RigidRotnSprlink::examine(int sub, MultiFront *mf)
{
  return examineBar2(sub, mf, nn);
}

PrioInfo
RigidThreeNodeShell::examine(int sub, MultiFront *mf)
{
  return examineTri3Shell(sub, mf, nn);
}

PrioInfo
RigidFourNodeShell::examine(int sub, MultiFront *mf)
{
  return examineQuad4Shell(sub, mf, nn);
}

PrioInfo
RigidEightNodeBrick::examine(int sub, MultiFront *mf)
{
  return examineHex8(sub, mf, nn);
}

PrioInfo
RigidSolid::examine(int sub, MultiFront *mf)
{
  switch(nnodes-numInternalNodes()) {
    case 4 : return examineTet4(sub, mf, nn); // 4-node tetra
    case 6 : return examinePenta6(sub, mf, nn); // 6-node penta
    case 8 : return examineHex8(sub, mf, nn); // 8-node hexa
    case 10 : return examineTet10(sub, mf, nn); // 10-node tetra
    case 15 : return examinePenta6(sub, mf, nn); // 15-node penta
    case 20 : return examineHex20(sub, mf, nn); // 20-node hexa
    case 26 : return examinePenta6(sub, mf, nn); // 26-node penta
    case 32 : return examineHex8(sub, mf, nn); // 32-node hexa
    default : return examineBar2(sub, mf, nn);
  }
}

PrioInfo
RigidSolid6Dof::examine(int sub, MultiFront *mf)
{
  switch(nnodes-numInternalNodes()) {
    case 4 : return examineTet4(sub, mf, nn); // 4-node tetra
    case 6 : return examinePenta6(sub, mf, nn); // 6-node penta
    case 8 : return examineHex8(sub, mf, nn); // 8-node hexa
    case 10 : return examineTet10(sub, mf, nn); // 10-node tetra
    case 15 : return examinePenta6(sub, mf, nn); // 15-node penta
    case 20 : return examineHex20(sub, mf, nn); // 20-node hexa
    case 26 : return examinePenta6(sub, mf, nn); // 26-node penta
    case 32 : return examineHex8(sub, mf, nn); // 32-node hexa
    default : return examineBar2(sub, mf, nn);
  }
}

PrioInfo
WeldedJoint::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
SphericalJoint::examine(int sub, MultiFront *mf)
{
  return examineBar2(sub, mf, nn);
}

PrioInfo
RevoluteJoint::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
TranslationalJoint::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
UniversalJoint::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
CylindricalJoint::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
PrismaticJoint::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
RevoluteDriver::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
PrismaticDriver::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
RevoluteActuator::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
PrismaticActuator::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
PinInSlotJoint::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
PlanarJoint::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
SphericalJointSpringCombo::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
TranslationalJointSpringCombo::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
UniversalJointSpringCombo::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
RevoluteJointSpringCombo::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
CylindricalJointSpringCombo::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
PrismaticJointSpringCombo::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
PinInSlotJointSpringCombo::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
RevoluteJointSpringComboWithFreeplay::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}

PrioInfo
PrismaticJointSpringComboWithFreeplay::examine(int sub, MultiFront *mf)
{
  return examineBeam2(sub, mf, nn);
}
#endif
