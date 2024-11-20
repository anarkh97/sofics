#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <map>

#include <Utils.d/DistHelper.h>

#include <Element.d/Element.h>
#include <Element.d/Truss.d/TwoNodeTruss.h>
#include <Element.d/Truss.d/TwoNodeTrussF.h>
#include <Element.d/Beam.d/EulerBeam.h>
#include <Element.d/Beam.d/TimoshenkoBeam.h>
#include <Element.d/Shell.d/ThreeNodeShell.h>
#include <Element.d/Shell.d/FourNodeShell.h>
#include <Element.d/Shell.d/Therm3NoShell.h>
#include <Element.d/Shell.d/Therm4NoShell.h>
#include <Element.d/Quad4.d/FourNodeQuad.h>
#include <Element.d/Brick.d/EightNodeBrick.h>
#include <Element.d/Tetra.d/Tetrahedral.h>
#include <Element.d/Tetra10.d/TenNodeTetrahedral.h>
#include <Element.d/Penta.d/Pentahedral.h>
#include <Element.d/Membrane.d/Membrane.h>
#include <Element.d/Membrane.d/FourNodeMembrane.h>
#include <Element.d/Spring.d/LinSpring.h>
#include <Element.d/Spring.d/TorSpring.h>
#include <Element.d/Spring.d/TransSprlink.h>
#include <Element.d/Spring.d/RotnSprlink.h>
#include <Element.d/CompShell.d/Compo3NodeShell.h>
#include <Element.d/CompShell.d/Compo4NodeShell.h>
#include <Element.d/FelippaShell.d/FelippaShell.h>
#include <Element.d/FelippaShell.d/FelippaShellX2.h>
#include <Element.d/Triangle3.d/Triangle3.h>
#include <Element.d/Triangle3.d/ThermTriangle.h>
#include <Element.d/ThermQuad.d/ThermQuadGal.h>
#include <Element.d/ThermQuad.d/Therm3DQuad.h>
#include <Element.d/Brick.d/ThermBrick.h>
#include <Element.d/Tetra.d/ThermIsoParamTetra.h>
#include <Element.d/Truss.d/Therm2NodeBar.h>
#include <Element.d/Helm.d/HelmQuadGal.h>
#include <Element.d/Helm.d/Tetra10HelmGal.h>
#include <Element.d/Helm.d/HelmQuadGls.h>
#include <Element.d/Helm.d/TetraHelmGal.h>
#include <Element.d/Helm.d/TetraHelmGLS.h>
#include <Element.d/Helm.d/HelmTri3Gal.h>
#include <Element.d/Helm.d/HelmTri3Gls.h>
#include <Element.d/Helm.d/HelmBrick.h>
#include <Element.d/Helm.d/HelmTri6Gal.h>
#include <Element.d/Helm.d/HelmQuad8Gal.h>
#include <Element.d/Helm.d/HelmLagQuadGal.h>
#include <Element.d/Helm.d/Tetra10HelmGal.h>
#include <Element.d/Helm.d/HelmBrickGLS.h>
#include <Element.d/Helm.d/HelmPenta.h>

#include <Element.d/Helm.d/HelmIsoParamHexa.h>
#include <Element.d/Helm.d/ThermIsoParamHexa.h>
#include <Element.d/Helm.d/HelmSpectralIsoParamHexa.h>
#include <Element.d/Helm.d/HelmIsoParamTetra.h>
#include <Element.d/Helm.d/HelmIsoParamQuad.h>
#include <Element.d/Helm.d/HelmSpectralIsoParamQuad.h>
#include <Element.d/Helm.d/HelmIsoParamTri.h>
#include <Element.d/Helm.d/LEIsoParamQuad.h>
#include <Element.d/Helm.d/LEIsoParamTri.h>
#include <Element.d/Helm.d/LEIsoParamHexa.h>
#include <Element.d/Helm.d/LEIsoParamTetra.h>

#include <Element.d/Shear.d/ShearPanel.h>
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
#include <Element.d/Brick20.d/Brick20.h>
#include <Element.d/NonLinearity.d/NLMembrane.h>
#include <Element.d/Shell.d/ConnectedTri.h>

#include <Element.d/MpcElement.d/MpcElement.h>
#include <Element.d/MpcElement.d/FsiElement.h>
#include <Element.d/MatrixElement.d/MatrixElement.h>

#include <map>
extern std::map<int,double> weightList;
extern std::map<int,double> fieldWeightList;
#ifdef USE_EIGEN3
#include <Element.d/Rigid.d/RigidBeam.h>
#include <Element.d/Rigid.d/RigidSpring.h>
#include <Element.d/Rigid.d/RigidTransSprlink.h>
#include <Element.d/Rigid.d/RigidRotnSprlink.h>
#include <Element.d/Rigid.d/RigidEightNodeBrick.h>
//#include <Element.d/Rigid.d/RigidTetrahedral.h>
#include <Element.d/Rigid.d/RigidSolid.h>
#include <Element.d/Rigid.d/RigidThreeNodeShell.h>
#include <Element.d/Rigid.d/RigidTwoNodeTruss.h>
#include <Element.d/Rigid.d/RigidSolid6Dof.h>
#include <Element.d/Rigid.d/RigidFourNodeShell.h>
#include <Element.d/Joint.d/BuildingBlocks.d/CommonPointConstraint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/ParallelAxesConstraint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/StraightLinePointFollowerConstraint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/RotationBlockerConstraint.h>
#include <Element.d/Joint.d/BuildingBlocks.d/ConstantDistanceConstraint.h>
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
#include <Element.d/MpcElement.d/PointPointDistanceConstraintElement.h>
#include <Element.d/MpcElement.d/PointLineDistanceConstraintElement.h>
#include <Element.d/MpcElement.d/PointPlaneDistanceConstraintElement.h>
#include <Element.d/MpcElement.d/LineLineDistanceConstraintElement.h>
#include <Element.d/MpcElement.d/SegmentSegmentDistanceConstraintElement.h>
#include <Element.d/MpcElement.d/PointVariPointDistanceConstraintElement.h>
#include <Element.d/MpcElement.d/PointVariLineDistanceConstraintElement.h>
#include <Element.d/MpcElement.d/PointVariPlaneDistanceConstraintElement.h>
#include <Element.d/MpcElement.d/PointVariPlaneSegmentDistanceConstraintElement.h>
#include <Element.d/MpcElement.d/PointVariPlaneSegmentDistanceConstraintElement2.h>
#include <Element.d/MpcElement.d/LineVariLineDistanceConstraintElement.h>
#include <Element.d/MpcElement.d/SegVariSegDistanceConstraintElement.h>
#include <Element.d/Joint.d/LinearTranslationalSpring.h>
#include <Element.d/Joint.d/NonlinearTranslationalSpring.h>
#include <Element.d/Joint.d/NonlinearTorsionalSpring.h>
#include <Element.d/Joint.d/SphericalJointSpringCombo.h>
#include <Element.d/Joint.d/TranslationalJointSpringCombo.h>
#include <Element.d/Joint.d/RevoluteJointSpringCombo.h>
#include <Element.d/Joint.d/UniversalJointSpringCombo.h>
#include <Element.d/Joint.d/CylindricalJointSpringCombo.h>
#include <Element.d/Joint.d/PrismaticJointSpringCombo.h>
#include <Element.d/Joint.d/PinInSlotJointSpringCombo.h>
#include <Element.d/Joint.d/RevoluteJointSpringComboWithFreeplay.h>
#include <Element.d/Joint.d/PrismaticJointSpringComboWithFreeplay.h>
#include <Element.d/Force.d/FollowerMomentElement.h>
#include <Element.d/Force.d/FollowerForceElement.h>
#include <Element.d/Force.d/PseudoTangentialMomentElement.h>
#include <Element.d/DiscreteMass.d/DiscreteMass6Dof.h>
#include <Element.d/Force.d/HexaQ1P0.h>
#include <Element.d/Force.d/HexaQ2P0.h>
#include <Element.d/Force.d/HexaQ2P1.h>
#include <Element.d/Force.d/IncompressibleHexaQ1P0.h>
#endif

#include <Element.d/Brick32.d/Brick32.h> 
#include <Element.d/Penta26.d/Penta26.h> 
#include <Element.d/Helm.d/HelmBrick32.h> 
#include <Element.d/Helm.d/HelmPenta26.h> 
#include <Element.d/Penta15.d/Penta15.h> 

#include <Element.d/DEM.d/DEMHelm2d.h>
#include <Element.d/DEM.d/DEMHelm3d.h>
#include <Element.d/DEM.d/DEMLE2d.h>
#include <Element.d/DEM.d/DEMLE3d.h>

#include <Element.d/FluidQuad.d/SloshQuadGal.h> 
#include <Element.d/FluidQuad.d/BarSloshFS.h>
#include <Element.d/FluidTetra.d/SloshTetra.h> 
#include <Element.d/FluidTriangle3.d/SloshTriangleFS.h> 

#include <Element.d/FluidQuad.d/HEVibQuadGal.h> 
#include <Element.d/FluidTetra.d/HEVibTetra.h> 

#include <Element.d/BelytschkoTsayShell.d/BelytschkoTsayShell.h>
#include <Driver.d/Domain.h>
#include <Parser.d/AuxDefs.h>
#include <Utils.d/SolverInfo.h>

extern SolverInfo &solInfo;
extern std::unique_ptr<ElementFactory> elemFact;

void
Elemset::elemadd(int num, int etype, int nnodes, int*n)
{
  Element *ele = elemFact->elemadd(num, etype, nnodes, n, ba);
  elemadd(num, ele);

  etypes.push_back(std::pair<int,int>(etype,num));
}

void
Elemset::setWeights()
{
  for(std::vector<std::pair<int,int> >::iterator it0 = etypes.begin(); it0 != etypes.end(); it0++) {

    int etype = it0->first;
    int num = it0->second;
    Element *ele = elem[num];
  

  }

  etypes.clear();
}

Element*
ElementFactory::elemadd(int num, int etype, int nnodes, int*n, BlockAlloc& ba)
{
   Element *ele;
   switch(etype) 
   {
     case 0:
       ele = new (ba) MatrixElement(nnodes, n);
       break;
     case 1:
       ele = new (ba) TwoNodeTruss(n);
       ele->setCategory(Element::Structural);
       break;
     case 2:
       ele = new (ba) FourNodeQuad(n);
       ele->setCategory(Element::Structural);
       break;
     case 3:
       ele = new (ba) Therm3DQuad(n);
       ele->setCategory(Element::Thermal);
       break;
     case 4:
       ele = new (ba) Triangle3(n);
        ele->setCategory(Element::Structural);
       break;
     case 6: {
       int nn[3] = { n[0], n[1], -1 };
       if(nnodes > 2) nn[2] = n[2];
       ele = new (ba) EulerBeam(nn);
       ele->setCategory(Element::Structural);
       break;
       }
     case 7: {
       int nn[3] = { n[0], n[1], -1 };
       if(nnodes > 2) nn[2] = n[2];
       ele = new (ba) TimoshenkoBeam(nn);
       ele->setCategory(Element::Structural);
       break;
       }
     case 8:
       ele = new (ba) ThreeNodeShell(n);
       ele->setCategory(Element::Structural);
       break;
     case 9:
       ele = new (ba) Therm2NodeBar(n);
       ele->setCategory(Element::Thermal);
       break;
     case 10:
       ele = new (ba) ThermQuadGal(n);
       ele->setCategory(Element::Thermal);
       break;
     case 11:
       ele = new (ba) TorSpring(n);
       ele->setCategory(Element::Structural);
       break;
     case 12:
       ele = new (ba) LinSpring(n);
       ele->setCategory(Element::Structural);
       break;
     case 15: case 1515:
       if(nnodes == 3 || (nnodes == 4 && n[2] == n[3]))
         ele = new (ba) FelippaShell(n);
       else 
         ele = new (ba) FelippaShellX2(n);
       ele->setCategory(Element::Structural);
       break;
     case 16:
       if(nnodes == 3) {
         int n_copy[4] = { n[0], n[1], n[2], n[2] };
         ele = new (ba) BelytschkoTsayShell(n_copy);
       }
       else 
         ele = new (ba) BelytschkoTsayShell(n);
       ele->setCategory(Element::Structural);
       break;
     case 17:
       ele = new (ba) EightNodeBrick(n);
       ele->setCategory(Element::Structural);
       break;
     case 18:
       ele = new (ba) ShearPanel(n);
       ele->setCategory(Element::Structural);
       break;
     case 19:
       ele = new (ba) Membrane(n);
       ele->setCategory(Element::Structural);
       break;
     case 20:
       ele = new (ba) Compo3NodeShell(n);
       ele->setCategory(Element::Structural);
       break;
     case 21:
       ele = new (ba) TransSprlink(n);
       ele->setCategory(Element::Structural);
       break;
     case 22:
       ele = new (ba) RotnSprlink(n);
       ele->setCategory(Element::Structural);
       break;
     case 23:
       ele = new (ba) Tetrahedral(n);
       ele->setCategory(Element::Structural);
       break;
     case 24:
       ele = new (ba) Pentahedral(n);
       ele->setCategory(Element::Structural);
       break;
     case 25:
       ele = new (ba) TenNodeTetrahedral(n);
       ele->setCategory(Element::Structural);
       break;
     case 30:
       ele = new (ba) HelmQuadGal(n);
       ele->setCategory(Element::Acoustic);
       break;
     case 31:
       ele = new (ba) HelmQuadGls(n);
       ele->setCategory(Element::Acoustic);
       break;
     case 32:
       ele = new (ba) HelmQuad8Gal(n);
       ele->setCategory(Element::Acoustic);
       break;
     case 35:
       ele = new (ba) HelmTri3Gal(n);
       ele->setCategory(Element::Acoustic);
       break;
     case 36:
       ele = new (ba) HelmTri3Gls(n);
       ele->setCategory(Element::Acoustic);
       break;
     case 38:
       ele = new (ba) HelmTri6Gal(n);
       ele->setCategory(Element::Acoustic);
       break;
     case 40:
       ele = new (ba) TetraHelmGal(n);
       ele->setCategory(Element::Acoustic);
       break;
     case 41:
       ele = new (ba) TetraHelmGLS(n);
       ele->setCategory(Element::Acoustic);
       break;
     case 42:
       ele = new (ba) Tetra10HelmGal(n);
       ele->setCategory(Element::Acoustic);
       break;
     case 43:
       ele = new (ba) HelmLagQuadGal(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 44:
       ele = new (ba) HelmBrickGLS(n);
       ele->setCategory(Element::Acoustic);
       break;
     case 45:
       ele = new (ba) HelmBrick(n);
       ele->setCategory(Element::Acoustic);
       break;
     case 46:
       ele = new (ba) Therm3NoShell(n);
       ele->setCategory(Element::Thermal);
       break;
     case 47:
       ele = new (ba) BarConvec(n);
       ele->setCategory(Element::Thermal);
       break;
     case 48:
       ele = new (ba) QuadConvec(n);
       ele->setCategory(Element::Thermal);
       break;
     case 49:
       ele = new (ba) TriangleConvec(n);
       ele->setCategory(Element::Thermal);
       break;
     case 50:
       ele = new (ba) ThermIsoParamTetra(nnodes,n);
       ele->setCategory(Element::Thermal);
       break;     
     case 51:
       ele = new (ba) ThermBrick(n);
       ele->setCategory(Element::Thermal);
       break;
     case 53:
       ele = new (ba) ThermTriangle(n);
       ele->setCategory(Element::Thermal);
       break;
     case 56:
       ele = new (ba) BarRadiation(n);
       ele->setCategory(Element::Thermal);
       break;
     case 57:
       ele = new (ba) TriangleRadiation(n);
       ele->setCategory(Element::Thermal);
       break;
     case 58:
       ele = new (ba) QuadRadiation(n);
       ele->setCategory(Element::Thermal);
       break;
     case 63:
       ele = new (ba) HelmLagQuadGal63(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 65:
       ele = new (ba) RigidTwoNodeTrussWithMass(n);
       ele->setCategory(Element::Structural);
       break;
     case 66:
       ele = new (ba) RigidBeamWithMass(n);
       ele->setCategory(Element::Structural);
       break;
     case 67:
       ele = new (ba) RigidSpring(n);
       ele->setCategory(Element::Structural);
       break;
     case 68:
       ele = new (ba) RigidTransSprlink(n);
       ele->setCategory(Element::Structural);
       break;
     case 69:
       ele = new (ba) RigidRotnSprlink(n);
       ele->setCategory(Element::Structural);
       break;
     case 70:
       ele = new (ba) RigidEightNodeBrick(n);
       ele->setCategory(Element::Structural);
       break;
     case 71:
       ele = new (ba) RigidSolid(nnodes,n);
       ele->setCategory(Element::Structural);
       break;
     case 72:
       ele = new (ba) Brick20(n);
       ele->setCategory(Element::Structural);
       break;
     case 73:
       ele = new (ba) RigidThreeNodeShell(n);
       ele->setCategory(Element::Structural);
       break;
     case 74:
       ele = new (ba) RigidSolid6Dof(nnodes,n);
       ele->setCategory(Element::Structural);
       break;
     case 76:
       ele = new (ba) RigidFourNodeShell(n);
       ele->setCategory(Element::Structural);
       break;
     case 77:
       ele = new (ba) PointPointDistanceConstraintElement(n);
       ele->setCategory(Element::Structural);
       break;
     case 78:
       ele = new (ba) PointLineDistanceConstraintElement(n);
       ele->setCategory(Element::Structural);
       break;
     case 79:
       ele = new (ba) PointPlaneDistanceConstraintElement(n);
       ele->setCategory(Element::Structural);
       break;
     case 80:
       ele = new (ba) ConnectedTri(n);
       ele->setCategory(Element::Structural);
       break;
     case 81:
       ele = new (ba) QuadContact(n);  // Interface element
       ele->setCategory(Element::Thermal);
       break;
     case 82:
       ele = new (ba) BrickContact(n);  // Interface element
       ele->setCategory(Element::Thermal);
       break;
     case 83:
       ele = new (ba) PentaContact(n);  // Interface element
       ele->setCategory(Element::Thermal);
       break;
     case 84:
       ele = new (ba) TriangleBulk(n);  // 2D Bulk Fluid Element
       ele->setCategory(Element::Thermal);
       break;
     case 85:
       ele = new (ba) TetraBulk(n);  // 3D Bulk Fluid Element
       ele->setCategory(Element::Thermal);
       break;
     case 86:
       ele = new (ba) PentaBulk(n);  // 3D Bulk Fluid Element
       ele->setCategory(Element::Thermal);
       break;
     case 87:
       ele = new (ba) FourNodeMembrane(n);
       ele->setCategory(Element::Structural);
       break;
     case 88: // Does one need to distinguish this version vs the implementation in backward mode?
       {
         // PJSA: first count number of unique nodes & if there are only 3 make a tri
         int i, j;
         int count = 0;
         bool isUnique[4] = { true, true, true, true };
         for(i=0; i<4; ++i)
           for(j=0; j<i; ++j)
             if((n[i] == n[j]) && isUnique[j]) { count++; isUnique[i] = false; }
         if(count > 1) {
           filePrint(stderr," *** ERROR: FourNodeShell type 88, number %d, has %d duplicate nodes \n",num+1,count);
           exit(-1);
         }
         else if(count == 1) {
           int new_n[3];
           j=0;
           for(i=0; i<4; ++i) if(isUnique[i]) new_n[j++] = n[i];
           //filePrint(stderr," *** WARNING: converting degenerate FourNodeShell to ThreeNodeShell, ele = %d,  nodes: %6d %6d %6d\n",
           //         num+1,new_n[0]+1,new_n[1]+1,new_n[2]+1);
           ele = new (ba) ThreeNodeShell(new_n);
         }
         else {
           ele = new (ba) FourNodeShell(n); // superelement comprising two ThreeNodeShells
         }
         ele->setCategory(Element::Structural);
       }
       break;
     case 90:
       ele = new (ba) HelmPenta(n); 
       ele->setCategory(Element::Acoustic);
       break;
     case 91:
       ele = new (ba) Brick32(n); 
       ele->setCategory(Element::Structural);
       break;
     case 92:
       ele = new (ba) Penta26(n); 
       ele->setCategory(Element::Structural);
       break;
     case 93:
       ele = new (ba) HelmBrick32(n); 
       ele->setCategory(Element::Acoustic);
       break;
     case 94:
       ele = new (ba) HelmPenta26(n); 
       ele->setCategory(Element::Acoustic);
       break;
     case 95:
       ele = new (ba) HelmIsoParamHexa(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 96:
       ele = new (ba) HelmIsoParamTetra(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 97:
       ele = new (ba) Penta15(n);
       ele->setCategory(Element::Structural);
       break;
     case 98:
       ele = new (ba) HelmIsoParamQuad(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 99:
       ele = new (ba) HelmIsoParamTri(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 100:
       ele = new (ba) LEIsoParamQuad(nnodes,n);
       ele->setCategory(Element::Structural);
       break;
     case 101:
       ele = new (ba) LEIsoParamTri(nnodes,n);
       ele->setCategory(Element::Structural);
       break;
     case 102:
       ele = new (ba) LEIsoParamHexa(nnodes,n);
       ele->setCategory(Element::Structural);
       break;
     case 103:
       ele = new (ba) LEIsoParamTetra(nnodes,n);
       ele->setCategory(Element::Structural);
       break;
     case 105:
       ele = new (ba) HelmSpectralIsoParamHexa(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 106:
       ele = new (ba) RigidBeamWithMass106(n,1);
       ele->setCategory(Element::Structural);
       break;
     case 108:
       ele = new (ba) HelmSpectralIsoParamQuad(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 109:
       ele = new (ba) ThermIsoParamHexa(nnodes,n);
       ele->setCategory(Element::Thermal);
       break;
     case 111:
       ele = new (ba) TwoNodeTrussF(n);
       ele->setCategory(Element::Structural);
       break;
     case 113:
       ele = new (ba) RotationBlockerConstraint(n, 2, 1);
       ele->setCategory(Element::Structural);
       break;
     case 114:
       ele = new (ba) CommonPointConstraint(n);
       ele->setCategory(Element::Structural);
       break;
     case 115:
       ele = new (ba) ConstantDistanceConstraint(n);
       ele->setCategory(Element::Structural);
       break;
     case 116:
       ele = new (ba) ParallelAxesConstraint(n);
       ele->setCategory(Element::Structural);
       break;
     case 117:
       ele = new (ba) StraightLinePointFollowerConstraint(n);
       ele->setCategory(Element::Structural);
       break;
     case 118:
       ele = new (ba) PlanarJoint(n);
       ele->setCategory(Element::Structural);
       break;
     case 119:
       ele = new (ba) WeldedJoint(n);
       ele->setCategory(Element::Structural);
       break;
     case 120:
       ele = new (ba) SphericalJoint(n); 
       ele->setCategory(Element::Structural);
       break;
     case 121:
       ele = new (ba) TranslationalJoint(n);
       ele->setCategory(Element::Structural);
       break;
     case 122:
       ele = new (ba) UniversalJoint(n);
       ele->setCategory(Element::Structural);
       break;
     case 123:
       ele = new (ba) RevoluteJoint(n);
       ele->setCategory(Element::Structural);
       break;
     case 124:
       ele = new (ba) CylindricalJoint(n);
       ele->setCategory(Element::Structural);
       break;
     case 125:
       ele = new (ba) PrismaticJoint(n);
       ele->setCategory(Element::Structural);
       break;
     case 126:
       ele = new (ba) RevoluteDriver(n);
       ele->setCategory(Element::Structural);
       break;
     case 127:
       ele = new (ba) PinInSlotJoint(n);
       ele->setCategory(Element::Structural);
       break;
     case 131:
       ele = new (ba) DiscreteMass6Dof(n);
       ele->setCategory(Element::Structural);
       solInfo.inertiaLumping = 2;
       break;
     case 132:
       ele = new (ba) RigidBeam(n);
       ele->setCategory(Element::Structural);
       break;
     case 133:
       ele = new (ba) RigidBeam133(n,1);
       ele->setCategory(Element::Structural);
       break;
     case 134:
       ele = new (ba) PrismaticDriver(n);
       ele->setCategory(Element::Structural);
       break;
     case 140:
       ele = new (ba) FollowerForceElement(n);
       ele->setCategory(Element::Structural);
       break;
     case 143:
       ele = new (ba) FollowerMomentElement(n);
       ele->setCategory(Element::Structural);
       break;
     case 149:
       ele = new (ba) PseudoTangentialMomentElement(n);
       ele->setCategory(Element::Structural);
       break;
//     case 150:
//       ele = new (ba) RigidTetrahedral(n);
//       ele->setCategory(Element::Structural);
//       break;
     case 173:
       ele = new (ba) SegVariSegDistanceConstraintElement(n);
       ele->setCategory(Element::Structural);
       break;
     case 174:
       ele = new (ba) SegmentSegmentDistanceConstraintElement(n);
       ele->setCategory(Element::Structural);
       break;
     case 175:
       ele = new (ba) LineVariLineDistanceConstraintElement(n);
       ele->setCategory(Element::Structural);
       break;
     case 176:
       ele = new (ba) LineLineDistanceConstraintElement(n);
       ele->setCategory(Element::Structural);
       break;
     case 177:
       ele = new (ba) PointVariPointDistanceConstraintElement(n);
       ele->setCategory(Element::Structural);
       break;
     case 178:
       ele = new (ba) PointVariLineDistanceConstraintElement(n);
       ele->setCategory(Element::Structural);
       break;
     case 179:
       ele = new (ba) PointVariPlaneDistanceConstraintElement(n);
       ele->setCategory(Element::Structural);
       break;
#if ((__cplusplus >= 201103L) || defined(HACK_INTEL_COMPILER_ITS_CPP11)) && defined(HAS_CXX11_TEMPLATE_ALIAS)
     case 180:
       ele = new (ba) HexaQ1P0(n);
       ele->setCategory(Element::Structural);
       break;
     case 181:
       ele = new (ba) HexaQ2P0(n);
       ele->setCategory(Element::Structural);
       break;
     case 182:
       ele = new (ba) HexaQ2P1(n);
       ele->setCategory(Element::Structural);
       break;
     case 280:
       ele = new (ba) IncompressibleHexaQ1P0(n);
       ele->setCategory(Element::Structural);
       break;
#endif
     case 200:
       ele = new (ba) LinearTranslationalSpring(n);
       ele->setCategory(Element::Structural);
       break;
     case 201:
       ele = new (ba) NonlinearTranslationalSpring(n, 0);
       ele->setCategory(Element::Structural);
       break;
     case 202:
       ele = new (ba) NonlinearTorsionalSpring(n, 2, 1);
       ele->setCategory(Element::Structural);
       break;
     case 203:
       ele = new (ba) LinearTranslationalSpring203(n, 1);
       ele->setCategory(Element::Structural);
       domain->solInfo().freeplay = true;
       break;
     case 204:
       ele = new (ba) NonlinearTranslationalSpring204(n, 0, 0, 1);
       ele->setCategory(Element::Structural);
       domain->solInfo().freeplay = true;
       break;
     case 205:
       ele = new (ba) NonlinearTorsionalSpring(n, 2, 1, 0, 1);
       ele->setCategory(Element::Structural);
       domain->solInfo().freeplay = true;
       break;
     case 220:
       ele = new (ba) SphericalJointSpringCombo(n);
       ele->setCategory(Element::Structural);
       break;
     case 221:
       ele = new (ba) TranslationalJointSpringCombo(n);
       ele->setCategory(Element::Structural);
       break;
     case 222:
       ele = new (ba) UniversalJointSpringCombo(n);
       ele->setCategory(Element::Structural);
       break;
     case 223:
       ele = new (ba) RevoluteJointSpringCombo(n);
       ele->setCategory(Element::Structural);
       break;
     case 224:
       ele = new (ba) CylindricalJointSpringCombo(n);
       ele->setCategory(Element::Structural);
       break;
     case 225:
       ele = new (ba) PrismaticJointSpringCombo(n);
       ele->setCategory(Element::Structural);
       break;
     case 226:
       ele = new (ba) RevoluteActuator(n);
       ele->setCategory(Element::Structural);
       break;
     case 227:
       ele = new (ba) PinInSlotJointSpringCombo(n);
       ele->setCategory(Element::Structural);
       break;
     case 234:
       ele = new (ba) PrismaticActuator(n);
       ele->setCategory(Element::Structural);
       break;
     case 279:
       ele = new (ba) PointVariPlaneSegmentDistanceConstraintElement(n);
       ele->setCategory(Element::Structural);
       break;
     case 323:
       ele = new (ba) RevoluteJointSpringComboWithFreeplay(n);
       ele->setCategory(Element::Structural);
       domain->solInfo().freeplay = true;
       break;
     case 325:
       ele = new (ba) PrismaticJointSpringComboWithFreeplay(n);
       ele->setCategory(Element::Structural);
       domain->solInfo().freeplay = true;
       break;
     case 379:
       ele = new (ba) PointVariPlaneSegmentDistanceConstraintElement379(n);
       ele->setCategory(Element::Structural);
       break;
     case 128:
       ele = new (ba) NLMembrane4(n);
       ele->setCategory(Element::Structural);
       break;
     case 129:
       ele = new (ba) NLMembrane(n);
       ele->setCategory(Element::Structural);
       break;
     case 2020:
       ele = new (ba) Compo4NodeShell(n);
       ele->setCategory(Element::Structural);
       break;
     case 4646:
       ele = new (ba) Therm4NoShell(n);
       ele->setCategory(Element::Thermal);
       break; 

     case 1100:
       ele = new (ba) DGMHelm2d_4(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1110:
       ele = new (ba) DGMHelm2d_4t(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1101:
       ele = new (ba) DGMHelm2d_8(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1111:
       ele = new (ba) DGMHelm2d_8t(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1102:
       ele = new (ba) DGMHelm2d_16(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1103:
       ele = new (ba) DGMHelm2d_32(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1104:
       ele = new (ba) DGMHelm2d_Eva2_8(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1120:
       ele = new (ba) DEMHelm2d_4(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1130:
       ele = new (ba) DEMHelm2d_4t(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1121:
       ele = new (ba) DEMHelm2d_8(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1131:
       ele = new (ba) DEMHelm2d_8t(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1122:
       ele = new (ba) DEMHelm2d_16(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1123:
       ele = new (ba) DEMHelm2d_32(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;

     case 1200:
       ele = new (ba) DGMLE2d_4(nnodes,n);
       ele->setCategory(Element::Structural);
       break;
     case 1201:
       ele = new (ba) DGMLE2d_16(nnodes,n);
       ele->setCategory(Element::Structural);
       break;
     case 1220:
       ele = new (ba) DEMLE2d_4(nnodes,n);
       ele->setCategory(Element::Structural);
       break;

     case 1150:
       ele = new (ba) DGMHelm3d_6(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1160:
       ele = new (ba) DGMHelm3d_6t(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1151:
       ele = new (ba) DGMHelm3d_26(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1161:
       ele = new (ba) DGMHelm3d_26t(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1152:
       ele = new (ba) DGMHelm3d_56(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1162:
       ele = new (ba) DGMHelm3d_56t(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1153:
       ele = new (ba) DGMHelm3d_98(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1170:
       ele = new (ba) DEMHelm3d_6(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1171:
       ele = new (ba) DEMHelm3d_26(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1172:
       ele = new (ba) DEMHelm3d_56(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;
     case 1173:
       ele = new (ba) DEMHelm3d_98(nnodes,n);
       ele->setCategory(Element::Acoustic);
       break;

     case 1250:
       ele = new (ba) DGMLE3d_6(nnodes,n);
       ele->setCategory(Element::Structural);
       break;
     case 1251:
       ele = new (ba) DGMLE3d_26(nnodes,n);
       ele->setCategory(Element::Structural);
       break;
     case 1252:
       ele = new (ba) DGMLE3d_50(nnodes,n);
       ele->setCategory(Element::Structural);
       break;

     case 1270:
       ele = new (ba) DEMLE3d_6(nnodes,n);
       ele->setCategory(Element::Structural);
       break;
     case 1271:
       ele = new (ba) DEMLE3d_26(nnodes,n);
       ele->setCategory(Element::Structural);
       break;
     case 1272:
       ele = new (ba) DEMLE3d_50(nnodes,n);
       ele->setCategory(Element::Structural);
       break;

     case 301:
       ele = new (ba) SloshQuadGal(n);
       ele->setCategory(Element::Fluid);
       break;
     case 302:
       ele = new (ba) BarSloshFS(n);
       ele->setCategory(Element::Fluid);
       break;
     case 311:
       ele = new (ba) SloshTetra(n);
       ele->setCategory(Element::Fluid);
       break;
     case 312:
       ele = new (ba) SloshTriangleFS(n);
       ele->setCategory(Element::Fluid);
       break;
     case 321:
       ele = new (ba) HEVibQuadGal(n);
       ele->setCategory(Element::Fluid);
       break;
     case 331:
       ele = new (ba) HEVibTetra(n);
       ele->setCategory(Element::Fluid);
       break;
//-----------------------------------------
     default:
       std::cerr << "Element Type " << etype << " is not Supported." << std::endl;
       assert(0);
       break;
   }
   ele->setElementType(etype);
   ele->setGlNum(num);
   return ele;
}

void
Elemset::mpcelemadd(int num, LMPCons *mpc, bool nlflag)
{
   Element *ele;
   if(mpc->getSource() == mpc::ContactSurfaces)
     ele = new MpcElement(mpc, nlflag);
   else
     ele = new (ba) MpcElement(mpc, nlflag);
   ele->setElementType(1001);
   elemadd(num, ele);

}

void
Elemset::fsielemadd(int num, LMPCons *fsi)
{
   Element *ele;
   ele = new (ba) FsiElement(fsi);
   ele->setElementType(1002);
   elemadd(num, ele);
}


