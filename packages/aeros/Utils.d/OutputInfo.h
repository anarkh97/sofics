#ifndef _OUTPUTINFO_H_
#define _OUTPUTINFO_H_

// To add an output type
//
// 1. add the type to the enum Type
// 2. add this new type to the lexer.l file
// 3. add an output message in Driver.d/Header.h
// 

#include <Utils.d/BinFileHandler.h>
#include <cstdio>

struct OutputInfo {
   enum Type
        { Displacement=0, Temperature, StressXX,    StressYY,    StressZZ,
          StressXY,     StressYZ,    StressXZ,    StrainXX,    StrainYY,
          StrainZZ,     StrainXY,    StrainYZ,    StrainXZ,    HeatFlXX,
          HeatFlXY,     HeatFlXZ,    GrdTempX,    GrdTempY,    GrdTempZ,
          StressVM,     AggrStVM,
          StressPR1,   StressPR2,   StressPR3,   StrainPR1,
          StrainPR2,    StrainPR3,   InXForce,    InYForce,    InZForce,
          AXMoment,     AYMoment,    AZMoment,    Energies,    AeroForce,
          EigenPair,    StrainVM,    Helmholtz,   Disp6DOF,    EigenPair6, 
          AeroXForce,   AcousticPressure, 
          AeroYForce,   AeroZForce,  AeroXMom,    AeroYMom,    AeroZMom,
          Velocity,     Acceleration,YModulus,    MDensity,    Thicknes,
          ShapeAtt,     ShapeStc,    Composit,    DispX,       DispY,
          DispZ,        RotX,        RotY,        RotZ,        DispMod,
          RotMod,       TotMod,      Rigid,       ElemToNode,  NodeToElem,
          NodeToNode,   AeroFlux,    HeatFlX,     GrdTemp,     Velocity6,
          Accel6,       ModeAlpha,   ModeError,   Reactions,   ModalDsp,
          ModalExF,     
          ContactPressure, NodalMassMat, NodalStifMat, Farfield, Kirchhoff,
          EigenPressure, FreqRespModes, HelmholtzModes,
          StressPR1Direc, StressPR2Direc, StressPR3Direc,
          StrainPR1Direc, StrainPR2Direc, StrainPR3Direc,
          EigenSlosh, SloshDispX, SloshDispY, SloshDispZ, SloshDisplacement,
          TDEnforcement, Damage, EquivalentPlasticStrain, 
          TemperatureFirstTimeDerivative, PressureFirstTimeDerivative, PressureSecondTimeDerivative,
          HeatReactions, Reactions6, Statevector, Residual, Jacobian, 
          RobData, SampleMesh, Accelvector, Forcevector,
          RotationMatrix, ExternalXForce, ExternalYForce, ExternalZForce,
          ExternalXMom, ExternalYMom, ExternalZMom, Velocvector, InternalStateVar, Quaternion,
          PlasticStrainXX, PlasticStrainYY, PlasticStrainZZ, PlasticStrainXY,
          PlasticStrainYZ, PlasticStrainXZ, BackStressXX, BackStressYY,
          BackStressZZ, BackStressXY, BackStressYZ, BackStressXZ, 
          WeigThic, WeigShap, VMstThic, VMstShap, VMstMach, VMstAlpha, VMstBeta,
          DispThic, DispShap, DispMach, DispAlph, DispBeta,
          AGstShap, AGstThic, 
          DissipatedEnergy, DeletedElements, DualStateVar, MuStateVar,
          Constraintvector, Constraintviolation, RomResidual, RomResidual6,
          RomExtForce, RomExtForce6, ModalMass, ModalStiffness, ModalDamping, ModalDynamicMatrix, ModalMatrices, ErrorIndicator,
          RotationVector, AngularVelocity, AngularAcceleration };

   enum Group  { Nodal, Attribute, NodeGroup };
   Type  type;
   int   interval;
   char* filename;
   FILE* filptr;
   BinFileHandler *binfilptr;
   int   averageFlg;
   int   surface;
   int   str_therm_option;
   int   width;
   int   precision;
   int   nodeNumber;	// To output just one node's information to a file.
   int   groupNumber;   // To output the nodes associated w/this group.
   int   dim;           // dimension of output data
   int   dataType;      // 1 for nodal (nodal_full or nodal_partial), 2 for elemental
   int   loadcase;      // for loadcase =-1 file is assigned to all load cases
   double ylayer,zlayer;
   int timeSliceRank; // Global rank of the corresponding time-slice (used by PITA)
   int ndtype; // 0 = deterministic, 1 = mean, 2 = stddev, 3 = pdf
   int sentype;  // 0 = non-sensitivity, 1 = sensitivity
   enum { realimag, modulusphase, animate };
   int complexouttype;
   int ncomplexout;   
   enum { spatial=0, convected, total };
   int angularouttype;
   bool rescaling; // whether or not the rotation vector is rescaled such that -pi <= phi <= pi (default is true)
   enum { Euler=0, Complement, Linear, ReducedEulerRodrigues, CayleyGibbsRodrigues, WienerMilenkovic, BauchauTrainelli };
   int rotvecouttype; // rotation vector parameterization (default is Euler, not to be confused with "Euler angles")
                      // reference: "The vectorial parameterization of rotation" by Bauchau and Trainelli
   bool matlab;
   bool PodRomfile;
   bool isFirst;
   int tdenforc_var; // CONFACE=1, NORMAL_FORCE_MAG, NORMAL_TRACTION_MAG, TANGENTIAL_FORCE_MAG, TANGENTIAL_TRACTION_MAG,
                     // CDIRNORX, CDIRNORY, CDIRNORZ, CDIRTANX, CDIRTANY, CDIRTANZ, SLIP_MAG, NODAL_DISSIPATION,
                     // CONTACT_AREA, GAP_CUR, GAP_OLD
   int topFlag; // if this is set to 1 then the output file should use the compressed numbering (i.e. with gaps removed)
                // compatible with top files generated using -T command line argument
   enum FrameType { Global=0, Local, Material };
   FrameType oframe;

   void initialize() {
     width = 10; 
     precision = 4;
     interval = 0;
     filptr = 0;
     binfilptr = 0;
     filename = 0;
     averageFlg = 1;
     surface = 2;
     str_therm_option = 0;
     nodeNumber = -1;
     groupNumber = -1;
     ylayer = 0.0;
     zlayer = 0.0;
     timeSliceRank = 0;
     ndtype = 0;
     sentype = 0;
     complexouttype = OutputInfo::realimag;
     ncomplexout = 16;
     angularouttype = OutputInfo::convected;
     rescaling = true;
     rotvecouttype = OutputInfo::Euler;
     tdenforc_var = 3;
     matlab = false;
     isFirst = true;
     PodRomfile = false;
     topFlag = 0;
     oframe = OutputInfo::Global;
   }

   void finalize(int numColumns) {
     if(numColumns == 6 && (type == Displacement)) type = Disp6DOF;
     else if(numColumns == 6 && (type == EigenPair)) type = EigenPair6;
     else if(numColumns == 6 && (type == Velocity)) type = Velocity6;
     else if(numColumns == 6 && (type == Acceleration)) type = Accel6;
     else if(numColumns == 6 && (type == Reactions)) type = Reactions6;
     else if(numColumns == 6 && (type == RomResidual)) type = RomResidual6;
     else if(numColumns == 6 && (type == RomExtForce)) type = RomExtForce6;

     if (averageFlg == 1 || averageFlg == 2 || averageFlg == 3)
       dataType = 1;
     else if (averageFlg == 0)
       dataType = 2;
 
     if (type == InXForce || type == InYForce ||
         type == InZForce || type == AXMoment ||
         type == AYMoment || type == AZMoment)
       dataType = 2;

     switch(type) {

       case EigenPair:
       case Displacement:
       case Reactions:
       case RotationVector:
         dim = 3;
         break;
  
       case AeroForce:
         dim = 3;
         break;

       case Disp6DOF:
       case EigenPair6:
       case Reactions6:
         dim = 6;
         break;

       case Energies:
         dim = 5;
         break;

       case Velocity:
       case Acceleration:
       case AngularVelocity:
       case AngularAcceleration:
         dim = 3;
         break;
 
       case Velocity6:
       case Accel6:
         dim = 6;
         break;

       case StressPR1Direc:
       case StressPR2Direc:
       case StressPR3Direc:
       case StrainPR1Direc:
       case StrainPR2Direc:
       case StrainPR3Direc:
         dim = 3;
         break;

       default:
         dim = 1;
         break;
     }
   };

 bool isStressOrStrain() {
   switch(type) {
     case StressXX: 
     case StressYY:
     case StressZZ:
     case StressXY:
     case StressYZ:
     case StressXZ:
     case StrainXX:
     case StrainYY:
     case StrainZZ:
     case StrainXY:
     case StrainYZ:
     case StrainXZ:
     case StressVM:
     case StressPR1:
     case StressPR2:
     case StressPR3:
     case StrainPR1:
     case StrainPR2:
     case StrainPR3:
     case StrainVM:
     case StressPR1Direc:
     case StressPR2Direc:
     case StressPR3Direc:
     case StrainPR1Direc:
     case StrainPR2Direc:
     case StrainPR3Direc:
     case EquivalentPlasticStrain:
     case PlasticStrainXX:
     case PlasticStrainYY:
     case PlasticStrainZZ:
     case PlasticStrainXY:
     case PlasticStrainYZ:
     case PlasticStrainXZ:
     case BackStressXX:
     case BackStressYY:
     case BackStressZZ:
     case BackStressXY:
     case BackStressYZ:
     case BackStressXZ:
     case Energies:
     case DissipatedEnergy:
     case Reactions:
     case Reactions6:
     case Damage:
       return true; 
       break;
     default:
       return false;
   }
 }

 bool isRotation() {
   switch(type) {
     case Disp6DOF: 
     case RotationVector:
     case RotationMatrix:
     case Quaternion:
     case Velocity6:
     case AngularVelocity:
     case Accel6:
     case AngularAcceleration:
     case RotX:
     case RotY:
     case RotZ:
     case RotMod:
     case TotMod:
       return true;
       break;
     default:
       return false;
   }
 }

 void copyParam(const OutputInfo& oI) {
    *this = oI;  
    filptr = 0;
    binfilptr = 0;
 }
  
};


#endif
