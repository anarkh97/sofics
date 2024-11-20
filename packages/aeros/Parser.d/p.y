%{
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <limits>
#include <map>
#include <cstdlib>
#include <Parser.d/AuxDefs.h>
#include <Element.d/NonLinearity.d/BilinPlasKinHardMat.h>
#include <Element.d/NonLinearity.d/CrushableFoam.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/2DMat.h>
#include <Element.d/NonLinearity.d/ExpMat.h>
#include <Element.d/NonLinearity.d/MaterialWrapper.h>
#include <Element.d/NonLinearity.d/PronyViscoElastic.h>
#include <Element.d/NonLinearity.d/OgdenMat.h>
#include <Element.d/NonLinearity.d/SimoElasticMat.h>
#include <Element.d/NonLinearity.d/SimoPlasticMat.h>
#include <Element.d/NonLinearity.d/NeoHookeanMat.h>
#include <Element.d/NonLinearity.d/MooneyRivlinMat.h>
#include <Element.d/NonLinearity.d/BrittleFractureTB.h>
#include <Element.d/NonLinearity.d/PlaneStressMat.h>
#include <Driver.d/Domain.h>
#include <Sfem.d/Sfem.h>
#include <Utils.d/DistHelper.h>
#include <Driver.d/GeoSource.h>
#include <Utils.d/Conwep.d/BlastLoading.h>
#ifdef STRUCTOPT
#include <Structopt.d/Driver_opt.d/Domain_opt.h>
#endif
#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif

 int numColumns = 3;
 double amplitude = 1.0;
 int PitaTS = 1;         //CD: Pita
 extern std::string clusterData_;
 extern std::string subdomains_;
 extern std::string decomposition_;
 extern std::string connectivity_;
 extern bool randomShuffle;
 extern bool allowMechanisms;
 extern bool useScotch;
%}

%union
{
 int ival;
 double fval;
 char * strval;
 NumedNode nval;
 NumedElem eval;
 NumList nl;
 BCond bcval;
 BCList *bclist;
 LMPCons *lmpcons;
 MFTTData *ymtt;
 MFTTData *ctett;
 MFTTData *sdetaft;
#ifdef USE_EIGEN3
 GenMFTTData<Eigen::Vector4d> *rubdaft;
#endif
 SS2DTData *ss2dt;
 ComplexBCList *cxbclist;
 ComplexBCond cxbcval;
 FrameData frame;
 NodalFrameData nframe;
 TrivialPair<MFTTData*,int> mftval;
 TrivialPair<MFTTData*,int> hftval;
 LayerData ldata;
 LayInfo *linfo;
 CoefData coefdata;
 DoubleList dlist;
 StringList slist;
 SurfaceEntity* SurfObj;
 MortarHandler* MortarCondObj;
 LMPCTerm* mpcterm;
 GeoSource::Rprop rprop;
 OutputInfo oinfo;
 ConstraintOptions copt;
 BlastLoading::BlastData blastData;
 SolverCntl* scntl;
 FreeplayProps freeplayProps;
 ModalParams::Type mpt;
}

%expect 6

%token ACTUATORS ADJOINTBASIS AERO AEROH AEROTYPE AGRESSIVETOLERANCES ALPROC AMAT ANALYSIS ARCLENGTH ATTRIBUTES ANGULAROUTTYPE ARUBBERMAT 
%token AUGMENT AUGMENTTYPE AUTOTOL AUXILIARY AVERAGED ATDARB ACOU ATDDNB ATDROB ARPACK ATDDIR ATDNEU ALLOWMECHANISMS
%token AUXCOARSESOLVER ACMECNTL ADDEDMASS AEROEMBED ANDESCLR ANDESCQR ANDESBETAB ANDESALPHA ANDESBETAM AUGMENTED
%token BLOCKDIAG BOFFSET BUCKLE BGTL BMPC BINARYINPUT BINARYOUTPUT BLOCKSIZE
%token CHECKTOKEN COARSESOLVER COEF CFRAMES COLLOCATEDTYPE CONVECTION COMPOSITE CONDITION CONTACT
%token CONTROL CORNER CORNERTYPE CURVE CCTTOL CCTSOLVER CRHS COUPLEDSCALE CONTACTSURFACES CMPC CNORM
%token COMPLEXOUTTYPE CONSTRMAT CASES CONSTRAINEDSURFACES CSFRAMES CSTYPE
%token CONSTANT CONWEP
%token DAMPING DblConstant DEFAULTPENALTY DELETEELEMENTS DEM DIMASS DISP DIRECT DLAMBDA DP DYNAM DETER DECOMPOSE DECOMPFILE DMPC DEBUGCNTL DEBUGICNTL DOCLUSTERING DOROWCLUSTERING ANGLE DUALBASIS DUALRB KMEANS CRANDOM
%token CONSTRAINTS MULTIPLIERS PENALTY
%token ELLUMP EIGEN EFRAMES ELSCATTERER END ELHSOMMERFELD ENGINEERING ETEMP EXPLICIT EXTFOL EPSILON ELEMENTARYFUNCTIONTYPE
%token FABMAT FABRICMAP FABRICMAT FACE FACOUSTICS FETI FETI2TYPE FETIPREC FFP FFPDIR FITALG FNAME FLAGMN FLUX FORCE FRONTAL FETIH FIELDWEIGHTLIST FILTEREIG FLUID FREEPLAY
%token FREQSWEEP FREQSWEEP1 FREQSWEEP2 FREQSWEEPA FREQSWEEPAW FSGL FSINTERFACE FSISCALING FSIELEMENT NOLOCALFSISPLITING FSICORNER FFIDEBUG FAILSAFE FRAMETYPE
%token GEPS GLOBALSEARCHCULL GLOBALTOL GRAVITY GRBM GTGSOLVER GLOBALCRBMTOL GROUP GROUPTYPE GOLDFARBTOL
%token HDIRICHLET HEAT HFETI HNEUMAN HSOMMERFELD HFTT
%token HELMHOLTZ HNBO HELMMF HELMSO HSCBO HWIBO HZEM HZEMFILTER HLMPC 
%token HERMITIAN HESSIAN
%token IACC IDENTITY IDIS IDIS6 ILUDROPTOL IntConstant INTERFACELUMPED ITEMP ITERTYPE IVEL IMESH 
%token INCIDENCE IHDIRICHLET IHDSWEEP IHNEUMANN ISOLVERTYPE INPC INFINTY
%token JACOBI KEYLETTER KRYLOVTYPE KIRLOC
%token LAYC LAYN LAYD LAYO LAYMAT LFACTOR LISRBM LMPC LOAD LOADCASE LOBPCG LOCALREDUCEDORDERBASES LOCALSOLVER LINESEARCH LUMPED
%token KSPARAM KSMAX 
%token MASS MASSAUGMENTATION MATERIALS MATLAB MAXITR MAXELEM MAXORTHO MAXVEC MODAL MPCPRECNO MPCPRECNOID MPCTYPE MPCTYPEID MPCSCALING MPCELEMENT MPCBLOCKID 
%token MPCBLK_OVERLAP MFTT MRHS MPCCHECK MUMPSICNTL MUMPSCNTL MUMPSMINEQ MUMPSSTRIDE MECH MODDAMP MODEFILTER MOMENTTYPE MPROJECT MAXIMUM
%token NOWARPEDVOLUME NDTYPE NEIGPA NEWMARK NewLine NEWTON NL NLMAT NLPREC NLMEMPTYPE NOCOARSE NODETOKEN NOMULTIPLEINTERACTIONS NONINPC
%token NSBSPV NLTOL NUMCGM NONORMALSMOOTHING NOSECONDARY NOGHOSTING NFRAMES NORMALSMOOTHINGDISTANCE
%token SENSITIVITY SENSITIVITYMETHOD SKIPPHYSICALFACES OUTPUT OUTPUT6 OUTPUTFRAME OLDDYNAMICSEARCH
%token QSTATIC QLOAD QUASISTATIC
%token PITA PITADISP6 PITAVEL6 NOFORCE MDPITA GLOBALBASES LOCALBASES TIMEREVERSIBLE REMOTECOARSE ORTHOPROJTOL READINITSEED JUMPCVG JUMPOUTPUT
%token PRECNO PRECONDITIONER PRELOAD PRESSURE PRINTMATLAB printmatlab PRINTNUMBER PROJ PIVOT PRECTYPE PRECTYPEID PICKANYCORNER PADEPIVOT PROPORTIONING PLOAD PADEPOLES POINTSOURCE PLANEWAVE PTOL PLANTOL PMAXIT PIECEWISE PARAMETERS
%token RADIATION RBMFILTER RBMSET READMODE READSENSITIVITY REBUILD RESOLUTIONMETHOD REVERSEORDER REDFOL RENUM RENUMBERID REORTHO RESTART RECONS RECONSALG REBUILDCCT RANDOM RPROP RNORM REVERSENORMALS ROBTYPE ROTVECOUTTYPE RESCALING RUBDAFT
%token SCALING SCALINGTYPE SDETAFT SENSORS SHELLFABRICMAP SHELLFABRICMAT SHELLSIMPLELOFTING SOLVERCNTL SOLVERHANDLE SOLVERTYPE SHIFT SHARPNONSHARPANGLE
%token SPOOLESTAU SPOOLESSEED SPOOLESMAXSIZE SPOOLESMAXDOMAINSIZE SPOOLESMAXZEROS SPOOLESMSGLVL SPOOLESSCALE SPOOLESPIVOT SPOOLESRENUM SPARSEMAXSUP SPARSEDEFBLK SPARSERENUM
%token STATS STRESSID SUBCYCLE SUBSPACE SURFACE STR_THERM_OPTION SAVEMEMCOARSE SPACEDIMENSION SCATTERER STAGTOL SCALED SWITCH STABLE SUBTYPE STEP SOWER SHELLTHICKNESS SURF SPRINGMAT SENSITIVITYID
%token TABLE TANGENT TDENFORCE TEMP TIME TOLEIG TOLFETI TOLJAC TOLPCG TOLSEN TOPFILE TOPOLOGY TRBM TRBMlc THERMOE THERMOH RATIOTOLSEN
%token TETT TOLCGM TURKEL TIEDSURFACES THETA PROJSOL CENTER POSELEM HOTSTART HRC THIRDNODE THERMMAT TDENFORC TESTULRICH THRU TRIVIAL NUMROMCPUS THICKNESSGROUPLIST
%token USE USERDEFINEDISP USERDEFINEFORCE UPROJ UNSYMMETRIC USING USESCOTCH
%token VERBOSE VERSION WETCORNERS YMTT SS1DT SS2DT YMST YSST YSSRT
%token ZERO BINARY GEOMETRY DECOMPOSITION GLOBAL MATCHER CPUMAP
%token NODALCONTACT MODE FRIC GAP
%token OUTERLOOP EDGEWS WAVETYPE ORTHOTOL IMPE FREQ DPH WAVEMETHOD
%token MATSPEC MATUSAGE BILINEARPLASTIC FINITESTRAINPLASTIC CRUSHABLEFOAM LINEARELASTIC STVENANTKIRCHHOFF TULERBUTCHER LINPLSTRESS READ OPTCTV ISOTROPICLINEARELASTIC VISCOLINEARELASTIC VISCOSTVENANTKIRCHHOFF NEOHOOKEAN VISCONEOHOOKEAN ISOTROPICLINEARELASTICJ2PLASTIC ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS HYPERELASTIC MOONEYRIVLIN VISCOMOONEYRIVLIN HENCKY OGDEN SIMOELASTIC SIMOPLASTIC LOGSTRAINPLASTIC SVKPLSTRESS VISCOLINPLSTRESS VISCOSVKPLSTRESS VISCOFABRICMAP SHELLVISCOFABRICMAP VISCOFABRICMAT SHELLVISCOFABRICMAT
%token PARTITIONGAP PLANESTRESSLINEAR PLANESTRESSSTVENANTKIRCHHOFF PLANESTRESSNEOHOOKEAN PLANESTRESSMOONEYRIVLIN PLANESTRESSBILINEARPLASTIC PLANESTRESSFINITESTRAINPLASTIC PLANESTRESSVISCOLINEARELASTIC PLANESTRESSVISCOSTVENANTKIRCHHOFF PLANESTRESSVISCONEOHOOKEAN PLANESTRESSVISCOMOONEYRIVLIN
%token SURFACETOPOLOGY MORTARTIED MORTARSCALING MORTARINTEGRATIONRULE SEARCHTOL STDMORTAR DUALMORTAR WETINTERFACE
%token NSUBS EXITAFTERDEC SKIP ROBCSOLVE RANDOMSAMPLE OUTPUTMEMORY OUTPUTWEIGHT SOLVER SPNNLSSOLVERTYPE MAXSIZE CLUSTERSOLVER CLUSTERSOLVERTYPE
%token WEIGHTLIST GMRESRESIDUAL 
%token SLOSH SLGRAV SLZEM SLZEMFILTER
%token PDIR HEFSB HEFRS HEINTERFACE
%token SNAPFI VELSNAPFI ACCSNAPFI DSVSNAPFI MUVSNAPFI PODROB ROMENERGY TRNVCT OFFSET ORTHOG SVDTOKEN CONVERSIONTOKEN CONVFI ROMRES SAMPLING SNAPSHOTPROJECT PODSIZEMAX REFSUBTRACT TOLER NORMALIZETOKEN FNUMBER SNAPWEIGHT ROBFI STAVCT VELVCT ACCVCT CONWEPCFG SCALEPOSCOORDS NODEPOSCOORDS MESHSCALEFACTOR PSEUDOGNAT PSEUDOGNATELEM USENMF USENMFC USEGREEDY USEPQN FILTERROWS
%token VECTORNORM REBUILDFORCE REBUILDCONSTRAINT SAMPNODESLOT REDUCEDSTIFFNESS UDEIMBASIS FORCEROB CONSTRAINTROB DEIMINDICES UDEIMINDICES SVDFORCESNAP SVDCONSTRAINTSNAP
%token USEMASSNORMALIZEDBASIS USECONSTANTMASS ONLINEMASSNORMALIZEBASIS STACKED
%token NUMTHICKNESSGROUP STRESSNODELIST DISPNODELIST RELAXATIONSEN 
%token QRFACTORIZATION QMATRIX RMATRIX XMATRIX EIGENVALUE
%token NPMAX BSSPLH PGSPLH
%token LIB

%type <bclist>   BCDataList IDisp6 TBCDataList PBCDataList AtdDirScatterer AtdNeuScatterer IDisp6Pita IVel6Pita
%type <bclist>   DirichletBC NeumanBC TempDirichletBC TempNeumanBC TempConvection TempRadiation ModalValList
%type <bclist>   HEVDirichletBC HEVDBCDataList HEVFRSBCList HEVFRSBC HEVFRSBCElem 
%type <bclist>   UsddLocations
%type <bcval>    BC_Data TBC_Data ModalVal PBC_Data HEVDBC_Data
%type <coefdata> CoefList
%type <cxbcval>  ComplexBC_Data ComplexMPCHeader
%type <mpcterm>  MPCLine ComplexMPCLine
%type <cxbclist> ComplexBCDataList ComplexNeumanBC ComplexDirichletBC 
%type <frame>    Frame
%type <nframe>   NodalFrame
%type <fval>     Float DblConstant DeleteElementsList
%type <ival>     AEROTYPE ALPROC AlProc Attributes AUGMENTTYPE AVERAGED 
%type <ival>     COLLOCATEDTYPE CORNERTYPE COMPLEXOUTTYPE TDENFORC CSTYPE ANGULAROUTTYPE ROTVECOUTTYPE
%type <ival>     DAMPING ELEMENTARYFUNCTIONTYPE FETIPREC FETI2TYPE FRAMETYPE
%type <ival>     Integer IntConstant ITERTYPE LoadCase
%type <ival>     RBMSET RENUMBERID OPTCTV
%type <mpt>      ROBTYPE
%type <rprop>    RPROP
%type <ival>     WAVETYPE WAVEMETHOD
%type <ival>     SCALINGTYPE SOLVERTYPE STRESSID SURFACE STR_THERM_OPTION MOMENTTYPE SPNNLSSOLVERTYPE CLUSTERSOLVERTYPE SENSITIVITYID
%type <ldata>    LayData LayoData LayMatData
%type <linfo>    LaycInfo LaynInfo LaydInfo LayoInfo
%type <mftval>   MFTTInfo
%type <hftval>   HFTTInfo
%type <ival>     NDTYPE OUTPUTFRAME
%type <ival>     GROUPTYPE
%type <nl>       NodeNums SommNodeNums 
%type <nval>     Node
%type <lmpcons>  MPCList ComplexMPCList MPCHeader
%type <strval>   FNAME 
%type <ymtt>     YMTTList SS1DTList YMSTList YSSTList YSSRTList
%type <ss2dt>    SS2DTList
%type <ctett>    TETTList
%type <sdetaft>  SDETAFList
%type <rubdaft>  RUBDAFList
%type <dlist>    FloatList
%type <slist>    StringList
%type <SurfObj>  FaceSet
%type <MortarCondObj> MortarCondition TiedSurfaces ContactSurfacesInfo
%type <ival>     MPCTYPEID MPCPRECNOID MPCBLOCKID
%type <ival>     ISOLVERTYPE RECONSALG
%type <ival>     PRECTYPEID SWITCH KEYLETTER Pressure FLAGMN
%type <oinfo>    OutInfo
%type <copt>     ConstraintOptionsData
%type <blastData> ConwepData
%type <scntl>    Solver SolverMethod Solvercntl
%type <freeplayProps> FreeplayProps
%%
FinalizedData:
	All END
	 { 
          if(domain->solInfo().piecewise || domain->solInfo().freeplay) domain->solInfo().activatePiecewise();
          return 0;
         }
	;
All:
	Component
	| All Component
	;
Component:
	NodeSet 
        | DirichletBC
        { if(geoSource->setDirichlet($1->n,$1->d) < 0) return -1; delete $1; }
        | NeumanBC
        { if(geoSource->setNeuman($1->n,$1->d) < 0) return -1; }
        | ModalNeumanBC
        | ModalLMPConstrain
        | LMPConstrain 
        | ComplexLMPConstrain 
	| ElemSet
	| FrameDList
        | NodalFrameDList
        | ConstrainedSurfaceFrameDList
	| Attributes
	{ int j = geoSource->getLocalIndex();
          if(geoSource->elementLumpingWeightLocalSize(j)>0) geoSource->setLocalIndex(j+1); }
        | Ellump
        | LocalReducedOrderBases
        | SampNodeSlot
 	| ReducedStiffness
	| UDeimBasis
	| Materials
        | Statics
        | Solvercntl
	| Pressure
	| Lumped
        {}
        | MassAugmentation
        | Preload
        {}
	| Renumbering
	| IDisp
	| Mode
        | IDisp6Pita
        | IVel6Pita
	| IDisp6
	{}
	| IVel
	| ITemp
        | ETemp
	| SensorLocations
	| ActuatorLocations
	| UsddLocations
        { if(geoSource->setUsddLocation($1->n,$1->d) < 0) return -1;
          if(geoSource->setDirichlet($1->n,$1->d) < 0)    return -1; delete $1; }
	| UsdfLocations
        | DynInfo
        | PrintMat
        | DeleteElements
        | Conwep
	| SloshInfo 
	| HEVibInfo 
        | QstaticInfo
	| Output
        | NLInfo
	| DiscrMasses
	| Composites
	| Cframes
        | LayMat
	| MFTTInfo
	  { domain->setMFTT($1.first, $1.second); }
        | HFTTInfo
          { domain->setHFTT($1.first, $1.second); }
        | LoadCase
        | YMTTable
        | TETTable
        | SS1DTable
        | SS2DTable
        | YSSTable
        | YMSTable
        | YSSRTable
        | SDETAFTable
        | RUBDAFTable
	| RbmTolerance
        | ToleranceInfo
        | Gravity
	| RbmFilterInfo
	| HzemFilterInfo
	| SlzemFilterInfo
        | ModeFilterInfo
        | Restart
	| LoadCInfo
	| UseCInfo
        | Control
        | Sensitivity
        | AnalysisInfo
        | OrthoInfo
	| AeroInfo
        | AeroHeatInfo
        | ThermohInfo
        | ThermoeInfo
        | HzemInfo
        | SlzemInfo
	| CondInfo
	| MassInfo
	| TempDirichletBC
        { if(geoSource->setDirichlet($1->n,$1->d) < 0) return -1; }
        | TempNeumanBC
        { if(geoSource->setNeuman($1->n,$1->d) < 0) return -1; }
        | TempConvection
        { if(geoSource->setNeuman($1->n,$1->d) < 0) return -1; }
	| HEVDirichletBC
	{ if(geoSource->setDirichletFluid($1->n,$1->d) < 0) return -1; }
	| HEVFRSBC
	{ if(geoSource->setDirichletFluid($1->n,$1->d) < 0) return -1; }
	| TempRadiation
        { if(geoSource->setNeuman($1->n,$1->d) < 0) return -1; }
        | TopInfo
        | OldHelmInfo
        | DEMInfo
	| HelmInfo
        | HelmMFInfo
        | HelmSOInfo
        | Scatterer
        | FarFieldPattern
        | FarFieldPatternDirs
        | KirchhoffLocations
        | HelmHoltzBC
        | EleHelmHoltzBC
        | EleScatterer
        | NeumScatterer
        | WetScatterer
        | HelmScatterer
        | ComplexNeumanBC
        { if(domain->setComplexNeuman($1->n,$1->d) < 0) return -1; }
        | IComplexNeumannBC
        | ComplexDirichletBC
        { if(domain->setComplexDirichlet($1->n,$1->d) < 0) return -1; }
        | IComplexDirichletBC
        | IComplexDirichletBCSweep
        | BinarySpec
        | AtdDirScatterer
        { if(geoSource->setDirichlet($1->n,$1->d) < 0) return -1; }
        | AtdNeuScatterer
        { if(geoSource->setNeuman($1->n,$1->d) < 0) return -1; }
	| AtdArbScatterer
	| AtdNeumScatterer
	| AtdRobinScatterer
        | Decompose
	| WeightList
        | FieldWeightList
	| NodalContact
	{}
        | ModeInfo
        | Noninpc
        | Inpc
        | Group
        | Random
	//| Shift
        | Impe
	| MatSpec
	| MatUsage
	| FaceSet
	{}
	| MortarCondition
	{}
        | TiedSurfaces
        {}
        | WetInterface
        {}
        | FSInterface
        {}
        | HEInterface
        {}
        | ContactSurfaces
        {}
        | ConstrainedSurfaces
        {}
	| BoffsetList
	| ParallelInTimeInfo 
        | Parameters 
        | Constraints
	| SvdToken
        | DeimIndices
        | UDeimIndices
	| Sampling
        | MeshScaleFactor
        | SnapshotProject
        | ConversionToken
        ;
Noninpc:
        NONINPC NewLine Integer Integer NewLine
          { domain->solInfo().noninpc = true;
            sfem->setOrder($3); 
            domain->solInfo().nsample = $4;
          }
         ;
Inpc:
        INPC NewLine Integer NewLine
          { domain->solInfo().inpc = true;
            sfem->setOrder($3);
          }
        ;
Group:
        GROUP NewLine
        | Group GROUPTYPE Integer Integer NewLine
        { if ($2 == OutputInfo::Attribute)  geoSource->setAttributeGroup($3-1, $4-1);
          else if ($2 == OutputInfo::Nodal)  geoSource->setNodeGroup($3-1, $4);
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Group Type: %d\n", $2);  exit(-1); }
        }
        | Group GROUPTYPE Integer Integer Integer NewLine
        { int i;
          if ($2 == OutputInfo::Attribute)  {
            for(i=$3; i<$4+1; ++i)
              geoSource->setAttributeGroup(i-1,$5-1);
          }
          else if ($2 == OutputInfo::Nodal)  {
            for(i=$3; i<$4+1; ++i)
              geoSource->setNodeGroup(i-1, $5);
          }
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Group Type: %d\n", $2);  exit(-1); }
        }
        | Group GROUPTYPE SURF Integer Integer NewLine
        { if ($2 == OutputInfo::Nodal) geoSource->setSurfaceGroup($4-1, $5);
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Surface Group Type: %d\n", $2);  exit(-1); }
        }
        ;
Random:
        RANDOM NewLine
        | Random Integer RPROP Float Float NewLine
        { geoSource->setGroupRandomProperty($2-1,$3,$4,$5); }
        ; 
Impe:
        IMPE NewLine
          { domain->solInfo().curSweepParam = 0; }
        | IMPE Integer NewLine
          { domain->solInfo().curSweepParam = $2; }
        | Impe FREQ Float NewLine
          { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe($3); }
        | Impe FREQSWEEP1 Float Float Integer NewLine
          { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe($3);
            domain->setFrequencySet(domain->solInfo().curSweepParam);
            domain->addFrequencies1(2.0*PI*$3, 2.0*PI*$4, $5); }
        | Impe FREQSWEEP2 Float Float Integer NewLine
          { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe($3);
            domain->setFrequencySet(domain->solInfo().curSweepParam);
            domain->addFrequencies2(2.0*PI*$3, 2.0*PI*$4, $5); }
        | Impe FREQSWEEP Float Float Integer Integer NewLine
          { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe($3);
            domain->setFrequencySet(domain->solInfo().curSweepParam);
            domain->addFrequencies(2.0*PI*$3, 2.0*PI*$4, $5, $6); }
        | Impe FREQSWEEPA Float Float Integer RECONSALG Float Integer Integer Integer Integer NewLine
        {
          if(domain->solInfo().curSweepParam == 0) geoSource->setImpe($3);
          domain->setFrequencySet(domain->solInfo().curSweepParam);
          domain->addFrequencies(2.0*PI*$3, 2.0*PI*$4, 2, $5);
          domain->solInfo().getSweepParams()->isAdaptSweep = true;
          domain->solInfo().getSweepParams()->adaptSweep.maxP = $8;
          domain->solInfo().getSweepParams()->adaptSweep.numS = $5;
          if ($6 == SweepParams::KrylovGalProjection) 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 0; 
          else if ($6 == SweepParams::WCAWEGalProjection) 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 2;
          else
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 1;
          domain->solInfo().getSweepParams()->adaptSweep.w1 = 2.0*PI*$3;
          domain->solInfo().getSweepParams()->adaptSweep.w2 = 2.0*PI*$4;
          domain->solInfo().getSweepParams()->adaptSweep.atol = $7;
          domain->solInfo().getSweepParams()->adaptSweep.minRHS = $9;
          domain->solInfo().getSweepParams()->adaptSweep.maxRHS = $10;
          domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = $11;
          domain->solInfo().getSweepParams()->nFreqSweepRHS = $10;
        }
        | Impe FREQSWEEPAW Float Float Integer RECONSALG Float Integer Integer Float Float NewLine
        {
          if(domain->solInfo().curSweepParam == 0) geoSource->setImpe($3);
          domain->setFrequencySet(domain->solInfo().curSweepParam);
          domain->addFrequencies(2.0*PI*$3, 2.0*PI*$4, 2, $5);
          domain->solInfo().getSweepParams()->isAdaptSweep = true;
          domain->solInfo().getSweepParams()->adaptSweep.maxP = $8;
          domain->solInfo().getSweepParams()->adaptSweep.numS = $5;
          if ($6 == SweepParams::KrylovGalProjection) 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 0; 
          else if ($6 == SweepParams::WCAWEGalProjection) 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 2;
          else
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 1;
          domain->solInfo().getSweepParams()->adaptSweep.w1 = 2.0*PI*$3;
          domain->solInfo().getSweepParams()->adaptSweep.w2 = 2.0*PI*$4;
          domain->solInfo().getSweepParams()->adaptSweep.atol = $7;
          domain->solInfo().getSweepParams()->adaptSweep.minRHS = $9;
          domain->solInfo().getSweepParams()->adaptSweep.maxRHS = $9;
          domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = -1;
          domain->solInfo().getSweepParams()->nFreqSweepRHS = $9;
          domain->solInfo().getSweepParams()->adaptSweep.ctolf = $10;
          domain->solInfo().getSweepParams()->adaptSweep.tol1f = $11;
        }
        | Impe FREQSWEEPA Float Float Integer RECONSALG NewLine
        {
          if(domain->solInfo().curSweepParam == 0) geoSource->setImpe($3);
          domain->setFrequencySet(domain->solInfo().curSweepParam);
          domain->addFrequencies(2.0*PI*$3, 2.0*PI*$4, 2, $5);
          domain->solInfo().getSweepParams()->isAdaptSweep = true;
          domain->solInfo().getSweepParams()->adaptSweep.maxP = 6;
          domain->solInfo().getSweepParams()->adaptSweep.numS = $5;
          if ($6 == SweepParams::KrylovGalProjection) {
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 0; 
             domain->solInfo().getSweepParams()->adaptSweep.atol = 1e-2;
             domain->solInfo().getSweepParams()->adaptSweep.minRHS = 8;
             domain->solInfo().getSweepParams()->adaptSweep.maxRHS = 48;
             domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = 4;
          }
          else if ($6 == SweepParams::WCAWEGalProjection) {
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 2;
             domain->solInfo().getSweepParams()->adaptSweep.atol = 1e-2;
             domain->solInfo().getSweepParams()->adaptSweep.minRHS = 8;
             domain->solInfo().getSweepParams()->adaptSweep.maxRHS = 48;
             domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = 4;
          }
          else {
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 1;
             domain->solInfo().getSweepParams()->adaptSweep.atol = 1e-2;
             domain->solInfo().getSweepParams()->adaptSweep.minRHS = 8;
             domain->solInfo().getSweepParams()->adaptSweep.maxRHS = 16;
             domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = 4;
          }
          domain->solInfo().getSweepParams()->adaptSweep.w1 = 2.0*PI*$3;
          domain->solInfo().getSweepParams()->adaptSweep.w2 = 2.0*PI*$4;
          domain->solInfo().getSweepParams()->nFreqSweepRHS = domain->solInfo().getSweepParams()->adaptSweep.maxRHS;
        }
        | Impe FreqSweep 
        | Impe ReconsInfo
        | Impe ImpeDampInfo
        | Impe PadePivotInfo
        | Impe PadePolesInfo
        ;
ImpeDampInfo:
        DAMPING Float Float NewLine
        {
          if($1 == 1) {
            domain->solInfo().setDamping($2,$3);
            domain->solInfo().getSweepParams()->alphaD = $3;
            domain->solInfo().getSweepParams()->betaD = $2;
            domain->solInfo().setDamping($2,$3);
          }
          else return -1; // only RAYDAMP is allowed here
        }
        ;
PadePivotInfo:
        PADEPIVOT Float NewLine
        { domain->solInfo().getSweepParams()->pade_pivot = true; domain->solInfo().getSweepParams()->pade_tol = $2; }
        ;
PadePolesInfo:
        PADEPOLES NewLine
        { domain->solInfo().getSweepParams()->pade_poles = true; }
        | PADEPOLES Float Float NewLine
        { domain->solInfo().getSweepParams()->pade_poles = true; 
          domain->solInfo().getSweepParams()->pade_poles_sigmaL = $2; domain->solInfo().getSweepParams()->pade_poles_sigmaU = $3; }
FreqSweep:  
        /* FREQSWEEP f_0 */
        FREQSWEEP Float NewLine
        { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe($2);
          domain->setFrequencySet(domain->solInfo().curSweepParam);
          domain->addCoarseFrequency(2.0*PI*$2); }
        /* f_i num_fine_intervals (i-1 to i) */
        | FreqSweep Float Integer NewLine
        { domain->addFrequencies(2.0*PI*$2, $3); }
        ;
ReconsInfo:
        RECONS RECONSALG Integer NewLine 
        { domain->solInfo().getSweepParams()->freqSweepMethod = $2; 
          int &l = domain->solInfo().getSweepParams()->padeL,
              &m = domain->solInfo().getSweepParams()->padeM,
              &n = domain->solInfo().getSweepParams()->padeN;
          switch($2) {
            case SweepParams::Taylor:
              n = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = $3+1; // taylor
              break;
            case SweepParams::Pade1:
              n = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+m+1;
              break;
            case SweepParams::Pade:
            case SweepParams::Fourier:
              n = $3;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = (int) ceil(float(l+m+1)/float(n));
              break;
            case SweepParams::PadeLanczos:
              n = $3;
              if(m%n != 0) m = (m/n+1)*n; // round m up to the nearest multiple of n
              l = m-1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = m/n;
              break;
            case SweepParams::GalProjection:
              n = $3;
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::KrylovGalProjection:
              n = $3;
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::WCAWEGalProjection:
              n = $3;
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
          }
        }
        | RECONS RECONSALG Integer Integer Integer NewLine  
        { domain->solInfo().getSweepParams()->freqSweepMethod = $2;
          int &l = domain->solInfo().getSweepParams()->padeL,
              &m = domain->solInfo().getSweepParams()->padeM,
              &n = domain->solInfo().getSweepParams()->padeN;
          switch($2) {
            case SweepParams::Taylor:
              n = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = $3+1; // taylor
              break;
            case SweepParams::Pade1:
              n = 1;
              l = $4;
              m = $5;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+m+1;
              break;
            case SweepParams::Pade:
            case SweepParams::Fourier:
              n = $3;
              l = $4; 
              m = $5;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = (int) ceil(float(l+m+1)/float(n));
              break;
            case SweepParams::PadeLanczos:
              n = $3;
              m = $5;
              if(m%n != 0) m = (m/n+1)*n; // round m up to the nearest multiple of n
              l = m-1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = m/n;
              break;
            case SweepParams::GalProjection:
              n = $3;
              l = $4;
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::KrylovGalProjection:
              n = $3;
              l = $4;
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::WCAWEGalProjection:
              n = $3;
              l = $4;
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
          }
        }
        ;
BinarySpec:
        BINARY NewLine
        | BinarySpec BINARYINPUT SWITCH NewLine 
          { geoSource->binaryInput = bool($3); }
        | BinarySpec BINARYINPUT SWITCH FNAME NewLine
          { geoSource->binaryInput = bool($3);
            std::string prefix = $4;
            clusterData_ = prefix + ".msh";
            decomposition_ = prefix + ".dec";
            connectivity_ = prefix + ".con";
            subdomains_ = prefix + ".sub";
          }
        | BinarySpec BINARYOUTPUT SWITCH NewLine 
          { geoSource->binaryOutput = bool($3); }
        | BinarySpec BINARYOUTPUT SWITCH FNAME NewLine
          { geoSource->binaryOutput = bool($3);
            int len = strlen($4);
            char *file = new char[len+5];
            strcpy(file, $4);
            strcat(file,".con");
            geoSource->setGlob(file);
          }
        | BinarySpec GEOMETRY FNAME NewLine
          { geoSource->setGeo($3); }
        | BinarySpec DECOMPOSITION FNAME NewLine
          { geoSource->setDecomp($3); }
        | BinarySpec GLOBAL   FNAME NewLine
          { geoSource->setGlob($3); }
        | BinarySpec MATCHER  FNAME NewLine
          { geoSource->setMatch($3); }
        | BinarySpec CPUMAP  FNAME NewLine
          { geoSource->setCpuMap($3); }
        ;
AnalysisInfo:
        ANALYSIS Integer NewLine
        { 
#ifdef STRUCTOPT	  
	  dynamic_cast<Domain_opt*>(domain)->addAnalysis($2); 
#endif
	}
        ;
Decompose :
       DECOMPOSE NewLine
         {if(decInit==0) decInit = new DecInit(); }
       | Decompose DECOMPFILE FNAME NewLine
         {decInit->file = strdup($3);}
       | Decompose NSUBS Integer NewLine
         {decInit->nsubs = $3; }
       | Decompose OUTPUTWEIGHT NewLine
         {decInit->weight = true; }
       | Decompose OUTPUTMEMORY NewLine
         {decInit->memory = true; }
       | Decompose EXITAFTERDEC NewLine
         {decInit->exitAfterDec = true;}
       | Decompose SKIP NewLine
         {decInit->skip = true;}
       | Decompose DETER NewLine
         {decInit->nosa = true; }
       | Decompose TRIVIAL NewLine
         {decInit->trivial = true; }
       | Decompose TRIVIAL Integer NewLine
         {decInit->trivial = true; randomShuffle = bool($3); }
       | Decompose FSGL NewLine
         {decInit->fsgl = true; }
       | Decompose ALLOWMECHANISMS NewLine
         {allowMechanisms = true; }
       | Decompose USESCOTCH NewLine
         {useScotch = true; }
       ;
WeightList :
       WEIGHTLIST NewLine
        {}
       | WeightList Integer Float NewLine
        { weightList[$2] = $3; }
        ;
FieldWeightList :
       FIELDWEIGHTLIST NewLine
        {}
       | FieldWeightList ACOU Integer NewLine
        { fieldWeightList[(int)Element::Acoustic] = $3; }
       | FieldWeightList MECH Integer NewLine
        { fieldWeightList[(int)Element::Structural] = $3; }
       | FieldWeightList HEAT Integer NewLine
        { fieldWeightList[(int)Element::Thermal] = $3; }
       | FieldWeightList FLUID Integer NewLine
        { fieldWeightList[(int)Element::Fluid] = $3; }
MFTTInfo:
	MFTT NewLine
	{ $$.first = new MFTTData; $$.second = 0; }
        | MFTT Integer NewLine
        { $$.first = new MFTTData; $$.second = $2; }
	| MFTTInfo Float Float NewLine
	{ $$.first->add($2,$3); }
	;
HFTTInfo:
        HFTT NewLine
        { $$.first = new MFTTData; $$.second = 0; }
        | HFTT Integer NewLine
        { $$.first = new MFTTData; $$.second = $2; }
        | HFTTInfo Float Float NewLine
        { $$.first->add($2,$3); }
        ;
LoadCase:
        LOADCASE NewLine
        { $$ = 0; }
        | LOADCASE Integer NewLine
        { $$ = $2; }
        | LOADCASE Integer GRAVITY SWITCH NewLine 
        { $$ = $2; domain->setLoadFactorGrav($2, $4); }
        | LOADCASE Integer TEMP SWITCH NewLine 
        { $$ = $2; domain->setLoadFactorTemp($2, $4); }
        | LOADCASE Integer GRAVITY SWITCH TEMP SWITCH NewLine
        { $$ = $2; domain->setLoadFactorGrav($2, $4); domain->setLoadFactorTemp($2, $6); }
        | LoadCase Integer Float NewLine
        { domain->setLoadFactor($$, $2, $3); }
        | LoadCase Integer MFTT Integer NewLine
        { domain->setLoadFactorMFTT($$, $2, $4); }
        | LoadCase Integer HFTT Integer NewLine
        { domain->setLoadFactorHFTT($$, $2, $4); }
        ;
Composites:
	COMPOSITE NewLine
	| Composites CoefInfo
	| Composites LaycInfo
	| Composites LaynInfo
        | Composites LaydInfo
        | Composites LayoInfo
	;
Cframes:
	CFRAMES NewLine
	| Cframes Frame
	{ geoSource->addCFrame($2.num,$2.d); }
	;
CoefInfo:
	COEF Integer NewLine CoefList
	{ $4.coefFlag = false; geoSource->addCoefInfo($2-1,$4); }
        | COEF Integer Float Float Float Float Float Float NewLine CoefList
        { $10.c[6][0] = $3;
          $10.c[6][1] = $4;
          $10.c[6][2] = $5;
          $10.c[6][3] = $6;
          $10.c[6][4] = $7;
          $10.c[6][5] = $8;
          $10.coefFlag = false;
          geoSource->addCoefInfo($2-1,$10); }
        | COEF Integer SWITCH NewLine CoefList
        { $5.coefFlag = $3; geoSource->addCoefInfo($2-1,$5); }
        | COEF Integer Float Float Float Float Float Float SWITCH NewLine CoefList
        { $11.c[6][0] = $3;
          $11.c[6][1] = $4;
          $11.c[6][2] = $5;
          $11.c[6][3] = $6;
          $11.c[6][4] = $7;
          $11.c[6][5] = $8;
          $11.coefFlag = $9;
          geoSource->addCoefInfo($2-1,$11); }
	;
CoefList:
	Integer Integer Float NewLine
	{ $$.zero(); $$.setCoef($1-1,$2-1,$3); }
	| CoefList Integer Integer Float NewLine
	{ $$.setCoef($2-1,$3-1,$4); }
	;
LaycInfo:
	LAYC Integer NewLine
	{ $$ = new LayInfo(0); geoSource->addLay($2-1,$$); }
	| LaycInfo LayData
	{ $1->add($2.lnum,$2.d,$2.matid); }
	;
LaynInfo:
	LAYN Integer NewLine 
	{ $$ = new LayInfo(1); geoSource->addLay($2-1,$$); }
	| LaynInfo LayData
	{ $1->add($2.lnum,$2.d,$2.matid); }
	;
LaydInfo:
        LAYD Integer NewLine
        { $$ = new LayInfo(0); geoSource->addLay($2-1,$$); }
        | LaydInfo LayoData
        { $1->add($2.lnum,$2.d,$2.matid); }
        ;
LayoInfo:
        LAYO Integer NewLine
        { $$ = new LayInfo(1); geoSource->addLay($2-1,$$); }
        | LayoInfo LayoData
        { $1->add($2.lnum,$2.d,$2.matid); }
        ;
LayData:
	Integer Float Float Float Float Float Float Float Float Float NewLine
	{ $$.lnum = $1-1;
          $$.matid = -1; // this means elastic constants are defined
          $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
	  $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7;
	  $$.d[6] = $8; $$.d[7] = $9; $$.d[8] = $10;
          $$.d[9] = 0;  $$.d[10] = 0; $$.d[11] = 0; }
        | Integer Float Float Float Float Float Float Float Float Float Float Float NewLine
        { $$.lnum = $1-1;
          $$.matid = -1; // this means elastic constants are defined
          $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7;
          $$.d[6] = $8; $$.d[7] = $9; $$.d[8] = $10;
          $$.d[9] = $11;$$.d[10]= $12;$$.d[11] = 0; }
        | Integer Float Float Float Float Float Float Float Float Float Float Float Float NewLine
        { $$.lnum = $1-1;
          $$.matid = -1; // this means elastic constants are defined
          $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7;
          $$.d[6] = $8; $$.d[7] = $9; $$.d[8] = $10;
          $$.d[9] = $11;$$.d[10]= $12; $$.d[11] = $13; }
	;
LayoData:
        Integer Integer Float Float NewLine // elastic constants to be read later from LAYMAT
        { $$.lnum = $1-1;  $$.matid = $2-1; $$.d[7] = $3; $$.d[8] = $4; }
        ;
LayMat:
        LAYMAT NewLine
        | LayMat LayMatData
        { geoSource->addLayMat($2.matid, $2.d); }
        ;
LayMatData:
        Integer Float Float Float Float Float NewLine // E1 E2 nu12 G12 rho
                                                      // note: coefficients of mutual influence are zero
                                                      //       coefficients of thermal expansion are zero
        { $$.matid = $1-1; $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = 0.0; $$.d[5] = 0.0; $$.d[6] = $6; 
          $$.d[7] = 0; $$.d[8] = 0; $$.d[9] = 0; }
        | Integer Float Float Float Float Float Float Float NewLine // E1 E2 nu12 G12 mu1,12 mu2,12 rho
                                                                    // note: coefficients of thermal expansion are zero
        { $$.matid = $1-1; $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7; $$.d[6] = $8;
          $$.d[7] = 0; $$.d[8] = 0; $$.d[9] = 0; }
        | Integer Float Float Float Float Float Float Float Float Float NewLine // E1 E2 nu12 G12 mu1,12 mu2,12 rho, cte1, cte2
                                                                                // note: the shear coefficient of thermal expansion is zero
        { $$.matid = $1-1; $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7; $$.d[6] = $8;
          $$.d[7] = $9; $$.d[8] = $10; $$.d[9] = 0; }
        | Integer Float Float Float Float Float Float Float Float Float Float NewLine // E1 E2 nu12 G12 mu1,12 mu2,12 rho, cte1, cte2, cte12
        { $$.matid = $1-1; $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7; $$.d[6] = $8; 
          $$.d[7] = $9; $$.d[8] = $10; $$.d[9] = $11; }
        ;
DiscrMasses:
	DIMASS NewLine
	| DiscrMasses Integer Integer Float NewLine
	  { domain->addDMass($2-1,$3-1,$4); }
        | DiscrMasses Integer Integer Integer Float NewLine
          { domain->addDMass($2-1,$3-1,$5,$4-1); }
        | DiscrMasses MODAL NewLine FNAME NewLine
        { domain->solInfo().modalDIMASS = true;
          domain->solInfo().reducedMassFile = $4; }
        ;
Gravity:
	GRAVITY NewLine
	| Gravity Float Float Float NewLine
	  { domain->setGravity($2,$3,$4); }
	;
Restart:
        RESTART NewLine
        | Restart FNAME FNAME FNAME NewLine
        { geoSource->getCheckFileInfo()->lastRestartFile = $2;
          geoSource->getCheckFileInfo()->outputExt = $3;
          geoSource->getCheckFileInfo()->FlagRST = $4; }
        | Restart FNAME FNAME NewLine
        { geoSource->getCheckFileInfo()->lastRestartFile = $2;
          geoSource->getCheckFileInfo()->outputExt = $3;}
        | Restart FNAME Integer NewLine
        { geoSource->getCheckFileInfo()->currentRestartFile = $2;
          domain->solInfo().nRestart = $3; }
        ;
LoadCInfo:
	LOAD FNAME NewLine
       { geoSource->setControlFile($2);
         geoSource->setControlRoutine((char *) "controlObj");}
	;
UseCInfo:
	USE FNAME NewLine
	{ geoSource->setControlRoutine($2); }
	;
SensorLocations:
        SENSORS NewLine BCDataList
        { for(int i=0; i<$3->n; ++i) $3->d[i].type = BCond::Sensors;
          if(geoSource->setSensorLocations($3->n,$3->d) < 0) return -1; }
	;
ActuatorLocations:
        ACTUATORS NewLine BCDataList
        { for(int i=0; i<$3->n; ++i) { $3->d[i].type = BCond::Actuators; }
          if(geoSource->setActuatorLocations($3->n,$3->d) < 0) return -1; 
          if(geoSource->setNeuman($3->n,$3->d) < 0)            return -1; }
	;
UsdfLocations:
	USERDEFINEFORCE NewLine BCDataList
        { geoSource->binaryInputControlLeft = true;
          for(int i=0; i<$3->n; ++i) { $3->d[i].type = BCond::Usdf; }
          if(geoSource->setUsdfLocation($3->n,$3->d) < 0) return -1;
          if(geoSource->setNeuman($3->n,$3->d) < 0)       return -1; } 
	;
UsddLocations:
        USERDEFINEDISP NewLine
        { geoSource->binaryInputControlLeft = true;
          $$ = new BCList; }
        | UsddLocations BC_Data
        { $2.type = BCond::Usdd; $$->add($2); }
        | UsddLocations Integer THRU Integer Integer NewLine
        { for(int i=$2; i<=$4; ++i) { BCond bc; bc.setData(i-1, $5-1, 0., BCond::Usdd); $$->add(bc); } }
        | UsddLocations Integer THRU Integer STEP Integer Integer NewLine
        { for(int i=$2; i<=$4; i+=$6) { BCond bc; bc.setData(i-1, $7-1, 0., BCond::Usdd); $$->add(bc); } }
        | UsddLocations SURF BC_Data
        { BCond *surf_bc = new BCond[1];
          surf_bc[0] = $3;
          surf_bc[0].type = BCond::Usdd;
          geoSource->addSurfaceDirichlet(1,surf_bc);
          if(geoSource->getNumSurfaceDirichlet() > 1) delete [] surf_bc; }
	;
Output:
	OUTPUT NewLine
        { numColumns = 3; } // set number of output columns to 3 
	| OUTPUT6 NewLine
	{ numColumns = 6; } // set number of output columns to 6 
        | OUTPUT Integer NewLine
        { numColumns = 3; geoSource->setOutLimit($2); }
        | OUTPUT6 Integer NewLine
        { numColumns = 6; geoSource->setOutLimit($2); }
        | OUTPUT KEYLETTER NewLine
        { numColumns = 3; domain->outFlag = $2; }
        | OUTPUT6 KEYLETTER NewLine
        { numColumns = 6; domain->outFlag = $2; }
        | OUTPUT KEYLETTER Integer NewLine
        { numColumns = 3; domain->outFlag = $2; geoSource->setOutLimit($3); }
        | OUTPUT6 KEYLETTER Integer NewLine
        { numColumns = 6; domain->outFlag = $2; geoSource->setOutLimit($3); }
        | OUTPUT FNAME NewLine
        { numColumns = 3; geoSource->getCheckFileInfo()->outputExt = $2; }
        | OUTPUT6 FNAME NewLine
        { numColumns = 6; geoSource->getCheckFileInfo()->outputExt = $2; }
        | Output OutInfo NewLine
        { $2.finalize(numColumns); geoSource->addOutput($2); }
        ;
OutInfo:
          STRESSID FNAME Integer // unformatted output for all nodes
        { $$.initialize(); $$.type = (OutputInfo::Type) $1; $$.filename = $2; $$.interval = $3; }
        | STRESSID Integer Integer FNAME Integer // formatted output for all nodes
        { $$.initialize(); $$.type = (OutputInfo::Type) $1; $$.width = $2; $$.precision = $3; $$.filename = $4; $$.interval = $5; }
        | STRESSID FNAME Integer Integer // unformatted output for one node
        { $$.initialize(); $$.type = (OutputInfo::Type) $1; $$.filename = $2; $$.interval = $3; $$.nodeNumber = $4-1; }
        | STRESSID FNAME Integer GROUPTYPE Integer // TDL: unformatted output for attribute, node group or node
        { $$.initialize(); $$.type = (OutputInfo::Type) $1; $$.filename = $2; $$.interval = $3; 
          if ($4 == OutputInfo::Nodal) $$.nodeNumber = $5; else $$.groupNumber = $5-1;}
        | STRESSID Integer Integer FNAME Integer Integer // formatted output for one node
        { $$.initialize(); $$.type = (OutputInfo::Type) $1; $$.width = $2; $$.precision = $3; $$.filename = $4; $$.interval = $5; $$.nodeNumber = $6-1; }
        | STRESSID Integer Integer FNAME Integer GROUPTYPE Integer // TDL: formatted output for node group or node
        { $$.initialize(); $$.type = (OutputInfo::Type) $1; $$.width = $2; $$.precision = $3; $$.filename = $4; $$.interval = $5; if ($6 == OutputInfo::NodeGroup) $$.groupNumber = $7; else $$.nodeNumber = $7-1; }
        | TDENFORC FNAME Integer // unformatted output for all nodes (for explicit dynamics tied/contact surfaces)
        { $$.initialize(); $$.type = OutputInfo::TDEnforcement; $$.tdenforc_var = $1; $$.filename = $2; $$.interval = $3; }
        | TDENFORC Integer Integer FNAME Integer // formatted output for all nodes (for explicit dynamics tied/contact surfaces)
        { $$.initialize(); $$.type = OutputInfo::TDEnforcement; $$.tdenforc_var = $1; $$.width = $2; $$.precision = $3; $$.filename = $4; $$.interval = $5; }
        | TDENFORC FNAME Integer Integer // unformatted output for one node (for explicit dynamics tied/contact surfaces)
        { $$.initialize(); $$.type = OutputInfo::TDEnforcement; $$.tdenforc_var = $1; $$.filename = $2; $$.interval = $3; $$.nodeNumber = $4-1; }
        | TDENFORC FNAME Integer GROUPTYPE Integer // TDL: unformatted output for group or one node (for explicit dynamics tied/contact surfaces)
        { $$.initialize(); $$.type = OutputInfo::TDEnforcement; $$.tdenforc_var = $1; $$.filename = $2; $$.interval = $3; if ($4 == OutputInfo::NodeGroup) $$.groupNumber = $5; else $$.nodeNumber = $5-1; }
        | TDENFORC Integer Integer FNAME Integer Integer // formatted output for one node (for explicit dynamics tied/contact surfaces)
        { $$.initialize(); $$.type = OutputInfo::TDEnforcement; $$.tdenforc_var = $1; $$.width = $2; $$.precision = $3; $$.filename = $4; $$.interval = $5; $$.nodeNumber = $6-1; }
        | TDENFORC Integer Integer FNAME Integer GROUPTYPE Integer // TDL: formatted output for group or one node (for explicit dynamics tied/contact surfaces)
        { $$.initialize(); $$.type = OutputInfo::TDEnforcement; $$.tdenforc_var = $1; $$.width = $2; $$.precision = $3; $$.filename = $4; $$.interval = $5; if ($6 == OutputInfo::NodeGroup) $$.groupNumber = $7; else $$.nodeNumber = $7-1; }

        | OutInfo NODETOKEN Integer
        { $$.nodeNumber = $3-1; }
        | OutInfo SURFACE
        { $$.surface = $2; }
        | OutInfo STR_THERM_OPTION
        { $$.str_therm_option = $2; }
        | OutInfo Float Float
        { $$.ylayer = $2; $$.zlayer = $3; }
        | OutInfo AVERAGED
        { $$.averageFlg = $2; }
        | OutInfo COMPLEXOUTTYPE
        { $$.complexouttype = $2; }
        | OutInfo COMPLEXOUTTYPE Integer
        { $$.complexouttype = $2; $$.ncomplexout = $3; }
        | OutInfo ANGULAROUTTYPE
        { $$.angularouttype = $2; }
        | OutInfo ROTVECOUTTYPE
        { $$.rotvecouttype = $2; }
        | OutInfo LINEARELASTIC
        { $$.rotvecouttype = OutputInfo::Linear; }
        | OutInfo RESCALING SWITCH
        { $$.rescaling = bool($3); }
        | OutInfo NDTYPE
        { $$.ndtype = $2; }
        | OutInfo NDTYPE Integer
        { $$.ndtype = $2; sfem->setnsamp_out($3); }
        | OutInfo OUTPUTFRAME
        { $$.oframe = (OutputInfo::FrameType) $2; }
        | OutInfo MATLAB 
        { $$.matlab = true; }
        | XMATRIX FNAME
        { domain->solInfo().xmatrixname = $2; }
        | QMATRIX FNAME
        { domain->solInfo().qmatrixname = $2; }
        | RMATRIX FNAME
        { domain->solInfo().rmatrixname = $2; }
        | EIGENVALUE FNAME
        { domain->solInfo().eigenvaluename = $2; }
        ;
DynInfo:
        DynamInfo
	| EIGEN NewLine
	{ domain->solInfo().setProbType(SolverInfo::Modal);
          domain->solInfo().eigenSolverType = SolverInfo::SubSpace; }
	| EIGEN Integer NewLine
	{ domain->solInfo().setProbType(SolverInfo::Modal);
	  domain->solInfo().nEig = $2;}
	| QRFACTORIZATION SWITCH NewLine
	{ domain->solInfo().qrfactorization = $2;}
	| NEIGPA Integer NewLine
	{ domain->solInfo().nEig = $2; }
	| SUBSPACE NewLine
	{ domain->solInfo().eigenSolverType = SolverInfo::SubSpace;}
        | SUBSPACE Integer Float Float NewLine
        { domain->solInfo().setSubSpaceInfo($2,$3,$4); }
	| NSBSPV Integer NewLine
	{ domain->solInfo().subspaceSize = $2;}
	| TOLEIG Float NewLine
	{ domain->solInfo().tolEig = $2; }
	| TOLJAC Float NewLine
	{ domain->solInfo().tolJac = $2; }
        | EXPLICIT NewLine
        { domain->solInfo().explicitK = true; }
        | SHIFT Float NewLine
        { geoSource->setShift($2); }
        | ARPACK NewLine
        { domain->solInfo().eigenSolverType = SolverInfo::Arpack; }
        | ARPACK FNAME NewLine
        { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = $2; }
        | ARPACK FNAME Integer NewLine
        { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = $2; 
          domain->solInfo().arpack_mode = $3; }
        | ARPACK Float Float NewLine
        { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValue($2, int($3)); }
        | ARPACK Float Float Integer NewLine
        { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValues($2, $3, $4);}
        | FILTEREIG SWITCH NewLine
        { domain->solInfo().filtereig = bool($2); }
        | JACOBI Integer NewLine
        { domain->solInfo().eigenSolverSubType = $2; } 
        | LOBPCG NewLine
        { domain->solInfo().eigenSolverType = SolverInfo::LobPcg;
          domain->solInfo().explicitK = true;}
        | DynInfo MAXITR Integer NewLine
        { domain->solInfo().maxitEig = $3; }
        | TESTULRICH NewLine
        { domain->solInfo().test_ulrich = true; }
        | ADDEDMASS Integer NewLine
        { domain->solInfo().addedMass = $2; }
        ;
PrintMat:
        PRINTMATLAB FNAME NewLine
        { domain->solInfo().printMatLab = true;
          domain->solInfo().printMatLabFile = $2;
          domain->solInfo().printMatLabExit = true; }
        | PRINTMATLAB FNAME EXITAFTERDEC NewLine
        { domain->solInfo().printMatLab = true;
          domain->solInfo().printMatLabFile = $2; 
          domain->solInfo().printMatLabExit = true; }
	;
DeleteElements:
        DELETEELEMENTS NewLine DeleteElementsList NewLine
        { domain->solInfo().elementDeletion = true; }
        | DeleteElements DeleteElementsList NewLine
        ;
DeleteElementsList:
        Float Integer
        { $$ = $1; domain->solInfo().deleteElements.insert(std::pair<int,double>($2-1,$1)); }
        | DeleteElementsList Integer
        { domain->solInfo().deleteElements.insert(std::pair<int,double>($2-1,$1)); }
        ;
SloshInfo:
        SLOSH NewLine
        { domain->solInfo().sloshing = 1; }
	| SLGRAV Float NewLine
	{ domain->setGravitySloshing($2); }
	;
MassInfo:
	MASS NewLine
        { domain->solInfo().massFlag = 1; }
        | MASS FNAME NewLine
        { domain->solInfo().massFlag = 1;
          domain->solInfo().massFile = std::string($2); }
	;
CondInfo:
	CONDITION NewLine Float Integer NewLine
        { domain->solInfo().setProbType(SolverInfo::ConditionNumber); 
	  domain->solInfo().setCondNumTol($3, $4); }
	| CONDITION NewLine
        { domain->solInfo().setProbType(SolverInfo::ConditionNumber);}
	;
TopInfo:
	TOPFILE NewLine
	{ domain->solInfo().setProbType(SolverInfo::Top); }
	;
ModalInfo:
	MODAL Integer
	{ domain->solInfo().modal_id.push_back($2); }
	| ModalInfo Integer
	{ domain->solInfo().modal_id.push_back($2); }
	| ModalInfo CONTACT Integer
	{ domain->solInfo().contact_modal_id.push_back($3); }
DynamInfo:
        DYNAM NewLine 
        { domain->solInfo().setProbType(SolverInfo::Dynamic); }
	| DynamInfo TimeIntegration
        | DynamInfo TimeInfo
        | DynamInfo DampInfo
        | DynamInfo MODAL NewLine
        { domain->solInfo().modal = true; domain->solInfo().modal_id.push_back(0); }
        | DynamInfo ModalInfo NewLine
        { domain->solInfo().modal = true; }
        | DynamInfo STABLE Integer NewLine
        { domain->solInfo().stable = $3; }
        | DynamInfo STABLE SWITCH NewLine
        { domain->solInfo().stable = $3; }
        | DynamInfo STABLE Integer Float Float Integer Integer NewLine
        { domain->solInfo().stable = $3;
          domain->solInfo().stable_cfl = $4;
          domain->solInfo().stable_tol = $5;
          domain->solInfo().stable_maxit = $6;
          domain->solInfo().stable_freq = $7;
        }
        | DynamInfo IACC SWITCH NewLine
        { domain->solInfo().iacc_switch = bool($3); }
        | DynamInfo ZERO SWITCH NewLine
        { domain->solInfo().zeroRot = bool($3); }
        | DynamInfo NOSECONDARY NewLine
        { domain->solInfo().no_secondary = true; }
        | DynamInfo TDENFORCE SWITCH NewLine
        { domain->solInfo().tdenforceFlag = $3; }
        | DynamInfo TDENFORCE SWITCH Integer Float NewLine
        { domain->solInfo().tdenforceFlag = $3;
          domain->solInfo().tdenforceMaxItr = $4;
          domain->solInfo().tdenforceTolAbs = $5; }
        | DynamInfo TDENFORCE SWITCH Integer Float Float NewLine
        { domain->solInfo().tdenforceFlag = $3;
          domain->solInfo().tdenforceMaxItr = $4;
          domain->solInfo().tdenforceTolAbs = $5;
          domain->solInfo().tdenforceInitia = $6; }
        | DynamInfo TDENFORCE SWITCH Integer Float Float Float NewLine
        { domain->solInfo().tdenforceFlag = $3;
          domain->solInfo().tdenforceMaxItr = $4;
          domain->solInfo().tdenforceTolAbs = $5;
          domain->solInfo().tdenforceInitia = $6;
          domain->solInfo().tdenforceFinal = $7; }
        | DynamInfo CHECKTOKEN NewLine
        { domain->solInfo().check_energy_balance = true; }
        | DynamInfo CHECKTOKEN Float Float NewLine
        { domain->solInfo().check_energy_balance = true;
          domain->solInfo().epsilon1 = $3; 
          domain->solInfo().epsilon2 = $4; }
        | DynamInfo QUASISTATIC SWITCH NewLine
        { domain->solInfo().quasistatic = bool($3); }
        ;
Conwep:
        CONWEP NewLine ConwepData NewLine
        { domain->solInfo().ConwepOnOff = true;
          BlastLoading::InputFileData = $3; 
        };
ConwepData:
        Float Float Float Float Float Integer Integer
        { // Note: chargeWeight must be entered in the units of mass of the problem, not units of force.
          $$.ExplosivePosition[0] = $1;
          $$.ExplosivePosition[1] = $2;
          $$.ExplosivePosition[2] = $3;
          $$.ExplosiveWeight = $4;
          $$.ExplosiveDetonationTime = $5;
          ($6 == 0 ? BlastLoading::BlastData::SurfaceBurst : BlastLoading::BlastData::AirBurst);
          $$.UnitConversionId = $7;
        };
TimeIntegration:
        NEWMARK NewLine
        { domain->solInfo().timeIntegration = SolverInfo::Newmark; }
        | MECH NewmarkSecondOrder NewLine
        | ACOU NewmarkSecondOrder NewLine
        { domain->solInfo().acoustic = true; }
        | HEAT NewmarkFirstOrder NewLine
        ;
NewmarkSecondOrder:
        Float Float // beta  gamma
        { domain->solInfo().setNewmarkSecondOrderInfo($1,$2); }
        | Float Float Float Float // beta  gamma alphaf alpham
        { domain->solInfo().setNewmarkSecondOrderInfo($1,$2,$3,$4); }
        | Float // rhoinfty
        { domain->solInfo().setNewmarkSecondOrderInfo(0.0,0.0,10.0,10.0,$1); }
        | Float Float Float // beta  gamma coef
        { domain->solInfo().setNewmarkSecondOrderInfo($1,$2);
          domain->solInfo().modifiedWaveEquation = true;
          domain->solInfo().modifiedWaveEquationCoef = $3; }
        ;
NewmarkFirstOrder:
        Float // epsilon
        { 
          if(domain->solInfo().probType == SolverInfo::NonLinDynam) {
            domain->solInfo().order = 1;
          }
          else 
            domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setNewmarkFirstOrderInfo($1); 
        }
        ;
QstaticInfo:
        QSTATIC NewLine
        { domain->solInfo().setProbType(SolverInfo::Dynamic); 
          domain->solInfo().timeIntegration = SolverInfo::Qstatic; }
        | QstaticInfo QMechInfo
        | QstaticInfo QHeatInfo
        | QstaticInfo MODAL NewLine
        { domain->solInfo().modal = true; domain->solInfo().modal_id.push_back(0); }
        | QstaticInfo MODAL Integer NewLine
        { domain->solInfo().modal = true; domain->solInfo().modal_id.push_back($3); }
	;
QMechInfo:
        MECH Float Float Integer NewLine
        { domain->solInfo().setQuasistaticInfo($2, 0, $3, $4); }
        | MECH Float Float Integer Float NewLine
        { domain->solInfo().setQuasistaticInfo($2, 0, $3, $4, $5); }
	;
QHeatInfo:
        HEAT Float Float Float Integer NewLine
        { domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setQuasistaticInfo($2, $3, $4, $5); }
/* need to check this
        | MECH Float Float Integer Float Float NewLine
        { domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setQuasistaticInfo($2, $3, $4, $5, $6); }
*/
	;
AeroInfo:
        AERO AEROTYPE NewLine
        { domain->solInfo().setAero($2); 
          domain->solInfo().isCollocated = 0; }
	| AERO NewLine AEROTYPE NewLine
        { domain->solInfo().setAero($3); 
          domain->solInfo().isCollocated = 0;
          if($3 == 20 || $3 == 22) { // set default alphas for C0
            domain->solInfo().alphas[0] = 0.5+0.375;
            domain->solInfo().alphas[1] = -0.375;
            domain->solInfo().alphasv    = 0.0;
          }
        } 
	| AERO NewLine AEROTYPE Float Float NewLine
        { domain->solInfo().setAero($3);
          domain->solInfo().isCollocated = 0;
          if($3 == 8) {
            // MPP uses only the first of the two inputted alphas
            domain->solInfo().mppFactor = $4;
          }
          else {
            // These alphas are used in FlExchanger::sendDisplacements and DistFlExchanger::sendDisplacements
            // As of 4/14/2014 the following schemes can use the displacement predictor on the structure side:
            // A0, A4, A5, A6, A7 and C0. Furthermore, we now apply a separate "anti-predictor" for A6 and A7
            // to compensate for the legacy predictor on the fluid side (in MatchNodeSet::getDisplacement)
            domain->solInfo().alphas[0] = $4+$5;
            domain->solInfo().alphas[1] = -$5;
            domain->solInfo().alphasv    = 0.0;
          }
        }
	| AERO NewLine AEROTYPE Float Float Float NewLine
        { domain->solInfo().setAero($3);
          domain->solInfo().isCollocated = 0;
          if($3 == 8) {
            // MPP uses only the first of the two inputted alphas
            domain->solInfo().mppFactor = $4;
          }
          else {
            // These alphas are used in FlExchanger::sendDisplacements and DistFlExchanger::sendDisplacements
            // As of 4/14/2014 the following schemes can use the displacement predictor on the structure side:
            // A0, A4, A5, A6, A7 and C0. Furthermore, we now apply a separate "anti-predictor" for A6 and A7
            // to compensate for the legacy predictor on the fluid side (in MatchNodeSet::getDisplacement)
            domain->solInfo().alphas[0] = $4+$5;
            domain->solInfo().alphas[1] = -$5;

            // This option added by Alex Main.  The last number on this line is the value used for velocity prediciton,
            // which is useful for embedded simulations.
            domain->solInfo().alphasv    = $6;
          }
        }

        | AERO NewLine AEROTYPE Float NewLine
        { domain->solInfo().setAero($3);
          domain->solInfo().isCollocated = 0;
          domain->solInfo().mppFactor = $4;
        }
	| AeroInfo COLLOCATEDTYPE NewLine
	{ domain->solInfo().isCollocated = $2; }
        | AeroInfo SUBCYCLE Integer NewLine
        { domain->solInfo().subcycle = $3; }
        | AeroInfo MATCHER FNAME NewLine
        { geoSource->setMatch($3); }
        | AeroInfo AeroEmbeddedSurfaceInfo NewLine
        {}
	;
AeroEmbeddedSurfaceInfo:
        AEROEMBED
        {}
        | AeroEmbeddedSurfaceInfo Integer
        { domain->AddAeroEmbedSurfaceId($2); }
	;
AeroHeatInfo:
        AEROH NewLine AEROTYPE Float Float NewLine
        { domain->solInfo().setAeroHeat($3, $4, $5); }
	;
ThermohInfo:
        THERMOH NewLine 
        { domain->solInfo().setThermoh(1); }
	;
ThermoeInfo:
        THERMOE NewLine 
        { domain->solInfo().setThermoe(1); }
	;
ModeInfo:
        MODE NewLine 
        { domain->solInfo().setModeDecomp(1); }
        | MODE Integer NewLine
        { domain->solInfo().setModeDecomp(1, $2); }
	;
HzemInfo:
        HZEM NewLine
        { domain->solInfo().hzemFlag=1; }
	;
SlzemInfo:
        SLZEM NewLine
        { domain->solInfo().slzemFlag=1; }
	;
RbmTolerance:
        TRBM NewLine Float NewLine
        { domain->solInfo().setTrbm($3); }
	;
ToleranceInfo:
        GRBM NewLine Float Float NewLine
        { domain->solInfo().setGrbm($3,$4); }
        | GRBM NewLine Float NewLine
        { domain->solInfo().setGrbm($3); }
        | GRBM NewLine 
        { domain->solInfo().setGrbm(); }
        | GRBM NewLine Float Float Integer NewLine
        { domain->solInfo().setGrbm($3,$4);
          domain->solInfo().grbm_use_lmpc = bool($5); }
        | GRBM Float Float Float NewLine Float Float NewLine
        { domain->solInfo().setGrbm($6,$7);
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = $2; grbm_ref[1] = $3; grbm_ref[2] = $4;
        }
        | GRBM Float Float Float NewLine Float NewLine
        { domain->solInfo().setGrbm($6);
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = $2; grbm_ref[1] = $3; grbm_ref[2] = $4;
        }
        | GRBM Float Float Float NewLine 
        { domain->solInfo().setGrbm();
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = $2; grbm_ref[1] = $3; grbm_ref[2] = $4;
        }
        | GRBM Float Float Float NewLine Float Float Integer NewLine
        { domain->solInfo().setGrbm($6,$7);
          domain->solInfo().grbm_use_lmpc = bool($8);
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = $2; grbm_ref[1] = $3; grbm_ref[2] = $4;
        }
	;
ModeFilterInfo:
        MODEFILTER Integer NewLine
        { domain->solInfo().modeFilterFlag = $2; }
        | MODEFILTER NewLine
        { domain->solInfo().modeFilterFlag = 1; }
        ;
RbmFilterInfo:
	RBMFILTER Integer NewLine
	{ domain->solInfo().useRbmFilter($2); }
	| RBMFILTER NewLine
	{ domain->solInfo().useRbmFilter(1); }
        | RBMFILTER Integer Integer NewLine
        { domain->solInfo().useRbmFilter($2);
          domain->solInfo().filterQ = $3; }
        | RbmFilterInfo RbmList NewLine
	;
RbmList:
    LISRBM Integer
    { if($2 < 1) {
        fprintf(stderr, " *** ERROR: mode %d specified under RBMFILTER is invalid.\n", $2);
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters.insert($2-1);
    }
    | RbmList Integer
    { if($2 < 1) {
        fprintf(stderr, " *** ERROR: mode %d specified under RBMFILTER is invalid.\n", $2);
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters.insert($2-1);
    }
    ;
HzemFilterInfo:
        HZEMFILTER NewLine
        { domain->solInfo().hzemFilterFlag=1; }
	;
SlzemFilterInfo:
        SLZEMFILTER NewLine
        { domain->solInfo().slzemFilterFlag=1; }
	;
TimeInfo:
	TIME Float Float Float NewLine
	{ domain->solInfo().setTimes($4,$3,$2); }
	;
ParallelInTimeInfo:
        PITA NewLine Integer Integer ParallelInTimeOptions
        {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime($3,$4,1);
        }
        |
        PITA NewLine Integer Integer Integer ParallelInTimeOptions
        {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime($3,$4,$5);
        }
        | MDPITA NewLine Integer Integer Integer Integer ParallelInTimeOptions
        {
          domain->solInfo().activatePita = true;
          domain->solInfo().mdPita = true;
          domain->solInfo().setParallelInTime($3,$4,$5); 
          /*domain->solInfo().numSpaceMPIProc = $6;*/
        }
        ;
ParallelInTimeOptions:
        NewLine ParallelInTimeKeyWord ParallelInTimeOptions
        | NewLine
        ;
ParallelInTimeKeyWord:    
        NOFORCE
        { domain->solInfo().pitaNoForce = true; }
        | GLOBALBASES Integer
        { domain->solInfo().pitaGlobalBasisImprovement = $2; }
        | LOCALBASES
        { domain->solInfo().pitaLocalBasisImprovement = 1; }
        | TIMEREVERSIBLE
        { domain->solInfo().pitaTimeReversible = true; }
        | REMOTECOARSE
        { domain->solInfo().pitaRemoteCoarse = true; }
        | ORTHOPROJTOL Float
        { domain->solInfo().pitaProjTol = $2; }
        | READINITSEED
        { domain->solInfo().pitaReadInitSeed = true; }
        | JUMPCVG
        { domain->solInfo().pitaJumpCvgRatio = 0.0; }
        | JUMPCVG Float
        { domain->solInfo().pitaJumpCvgRatio = $2; }
        | JUMPOUTPUT
        { domain->solInfo().pitaJumpMagnOutput = true; }
        ;
DampInfo:
	DAMPING Float Float NewLine
	{ if($1 == 1) domain->solInfo().setDamping($2,$3);
          else return -1; // only RAYDAMP is allowed here
        }
        | DAMPING Float Float MOMENTTYPE NewLine
        { if($1 == 1) {
            domain->solInfo().setDamping($2,$3);
            domain->solInfo().mtypeDamp = (int)$4;
          }
          else return -1;
        } 
	| MODDAMP NewLine ModalValList
	{ if(geoSource->setModalDamping($3->n, $3->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true;
        }
        ;
ComplexDirichletBC:
	HDIRICHLET NewLine ComplexBCDataList
	{ $$ = $3; }
	;
IComplexDirichletBC:
        IHDIRICHLET NewLine Float Float Float NewLine
        {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, $3, $4, $5);
        }
        | IHDIRICHLET NewLine Float Float NewLine
        {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, $3, $4, 0.0);
        }
        ;
IComplexDirichletBCSweep:
        IHDSWEEP NewLine Integer Float Float NewLine
        {
           domain->implicitFlag = 1;
           domain->solInfo().setProbType(SolverInfo::HelmholtzDirSweep);
           domain->setWaveDirections($3, $4, $5);
        }
        | IHDSWEEP Integer NewLine
        {
           domain->implicitFlag = 1;
           domain->setWaveDirections($2,0.0,0.0,0.0);
        }
        | IComplexDirichletBCSweep DirectionVector
        ;
DirectionVector:
        Integer Float Float Float NewLine
        {
           domain->setWaveDirections($1, $2, $3, $4);
        }
        ;
IComplexNeumannBC:
        IHNEUMANN NewLine Float Float Float NewLine
        {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, $3, $4, $5);
        }
        ;
DirichletBC:
        DISP NewLine
        { $$ = new BCList; }
        | DirichletBC BC_Data
        { $2.type = BCond::Displacements; $$->add($2); }
        | DirichletBC Integer THRU Integer Integer Float NewLine
        { for(int i=$2; i<=$4; ++i) { BCond bc; bc.setData(i-1, $5-1, $6, BCond::Displacements); $$->add(bc); } }
        | DirichletBC Integer THRU Integer STEP Integer Integer Float NewLine
        { for(int i=$2; i<=$4; i+=$6) { BCond bc; bc.setData(i-1, $7-1, $8, BCond::Displacements); $$->add(bc); } }
        | DirichletBC SURF BC_Data
        { BCond *surf_bc = new BCond[1];
          surf_bc[0] = $3;
          surf_bc[0].type = BCond::Displacements;
          geoSource->addSurfaceDirichlet(1,surf_bc);
          if(geoSource->getNumSurfaceDirichlet() > 1) delete [] surf_bc; }
/* TODO | DirichletBC SURF Integer Integer Float USING ConstraintOptionsData NewLine
        { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = $3-1; surf_bc[0].dofnum = $4-1; surf_bc[0].val = $5;
          surf_bc[0].type = BCond::Lmpc;
          geoSource->addSurfaceDirichlet(1,surf_bc); }
*/
        ;
ConstrainedSurfaces:
        CONSTRAINEDSURFACES NewLine
        /*                    surface_id type   attribute_id csframe_id */
        | ConstrainedSurfaces Integer    CSTYPE Integer      Integer      NewLine
        { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = $2-1;
          surf_bc[0].type = (BCond::BCType) $3; //BCond::PointPlaneDistance;
          surf_bc[0].dofnum = $4-1;
          surf_bc[0].val = $5-1;
          geoSource->addSurfaceConstraint(1,surf_bc);
          if(geoSource->getNumSurfaceConstraint() > 1) delete [] surf_bc;
        }
HEVDirichletBC:
        PDIR NewLine HEVDBCDataList
        { for(int i=0; i<$3->n; ++i) $3->d[i].type = BCond::Pdir; $$ = $3; }
	;
HEVDBCDataList:
	HEVDBC_Data
	{ $$ = new BCList; $$->add($1); }
	| HEVDBCDataList HEVDBC_Data
	{ $$ = $1; $$->add($2); }
	;
HEVDBC_Data:
	Integer Float NewLine
	{ $$.nnum = $1-1; $$.dofnum = 10; $$.val = $2; }
	| Integer NewLine
	{ $$.nnum = $1-1; $$.dofnum = 10; $$.val = 0.0; }
	;
HEVFRSBC:  // Added for HEV problem, EC, 20080512
        HEFRS NewLine HEVFRSBCList
        { domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; 

          int* allPDirNodes = new int[$3->n];

          for (int ii=0; ii < $3->n; ii++)
            allPDirNodes[ii]=($3->d[ii]).nnum;
          std::sort(allPDirNodes,allPDirNodes + $3->n);

          int maxFSNodes = 32;
          int* allHEVFSNodes = new int[maxFSNodes];
          allHEVFSNodes[0] = allPDirNodes[0];

          int numHEVFSNodes = 1;
          for (int ii = 1; ii < $3->n; ++ii)  {
            if (numHEVFSNodes == maxFSNodes)  {
              int nMaxFSNodes = maxFSNodes*3/2;
              int* nAllHEVFSNodes = new int[nMaxFSNodes];
              for (int kk = 0; kk < maxFSNodes; kk++)
                nAllHEVFSNodes[kk] = allHEVFSNodes[kk];
              delete [] allHEVFSNodes;
              allHEVFSNodes = nAllHEVFSNodes;
              maxFSNodes = nMaxFSNodes;
            }

            if (allPDirNodes[ii]!=allPDirNodes[ii-1])
              allHEVFSNodes[numHEVFSNodes++] = allPDirNodes[ii];
          }
          delete [] allPDirNodes;
          
          $$ = new BCList;

          for (int ii=0; ii<numHEVFSNodes; ii++)
            $$->add(allHEVFSNodes[ii],10,0.0); 
          delete [] allHEVFSNodes;

          for(int i=0; i<$$->n; ++i) $$->d[i].type = BCond::Hefrs;
        }
        ;
HEVFRSBCList:
        HEVFRSBCElem
        { $$ = new BCList; 
          for (int ii = 0; ii < $1->n; ii++) 
           $$->add(($1->d)[ii]); }
        | HEVFRSBCList HEVFRSBCElem
        { $$ = $1; 
          for (int ii = 0; ii < $2->n; ii++) 
           $$->add(($2->d)[ii]); }
        ;
HEVFRSBCElem: // Added for HEV problem, EC, 20080512
        Integer Integer SommNodeNums NewLine
        { $$ = new BCList;
          for(int i=0; i<$3.num; ++i) 
          { $$->add($3.nd[i],10,0.0); } }
        ;
TempDirichletBC:
        TEMP NewLine
        { $$ = new BCList; if(domain->solInfo().soltyp != 2) domain->solInfo().thermalLoadFlag = 1;}
        | TempDirichletBC Integer Float NewLine
        { $$ = $1; BCond bc; bc.nnum = $2-1; bc.dofnum = 6;
          bc.val = $3; bc.type = BCond::Temperatures; $$->add(bc); }
        | TempDirichletBC Integer THRU Integer Float NewLine
        { for(int i=$2; i<=$4; ++i) { BCond bc; bc.setData(i-1, 6, $5, BCond::Temperatures); $$->add(bc); } }
        | TempDirichletBC Integer THRU Integer STEP Integer Float NewLine
        { for(int i=$2; i<=$4; i+=$6) { BCond bc; bc.setData(i-1, 6, $7, BCond::Temperatures); $$->add(bc); } }
        | TempDirichletBC SURF Integer Float NewLine
        { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = $3-1;
          surf_bc[0].val = $4;
          surf_bc[0].dofnum = 6;
          surf_bc[0].type = BCond::Temperatures;
          geoSource->addSurfaceDirichlet(1,surf_bc); }
	;
TempNeumanBC:
        FLUX NewLine
        { $$ = new BCList; }
        | FLUX Integer NewLine
        { $$ = new BCList($2); }
        | TempNeumanBC Integer Float NewLine
        { $$ = $1; BCond bc; bc.nnum = $2-1; bc.dofnum = 6;
          bc.val = $3; bc.type = BCond::Flux; bc.loadsetid = $$->loadsetid; $$->add(bc); }
        | TempNeumanBC SURF Integer Float NewLine
        { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = $3-1;
          surf_bc[0].dofnum = 6;
          surf_bc[0].val = $4;
          surf_bc[0].type = BCond::Flux;
          surf_bc[0].loadsetid = $$->loadsetid;
          geoSource->addSurfaceNeuman(1,surf_bc);
          if(geoSource->getNumSurfaceNeuman() > 1) delete [] surf_bc; }
	;
TempConvection:
        CONVECTION NewLine
        { $$ = new BCList; }
        | CONVECTION Integer NewLine
        { $$ = new BCList($2); }
        | TempConvection Integer Float Float Float NewLine
        { $$ = $1; BCond bc; bc.nnum = $2-1; bc.dofnum = 6;
          bc.val = $3*$4*$5; bc.type = BCond::Convection; bc.loadsetid = $$->loadsetid; $$->add(bc); }
        | TempConvection SURF Integer Float Float Float NewLine
        { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = $3-1;
          surf_bc[0].dofnum = 6;
          surf_bc[0].val = $4*$5*$6;
          surf_bc[0].type = BCond::Convection;
          surf_bc[0].loadsetid = $$->loadsetid;
          geoSource->addSurfaceNeuman(1,surf_bc);
          if(geoSource->getNumSurfaceNeuman() > 1) delete [] surf_bc; }
	;
TempRadiation:
        RADIATION NewLine
        { $$ = new BCList; }
        | RADIATION Integer NewLine
        { $$ = new BCList($2); }
        | TempRadiation Integer Float Float Float NewLine
        { $$ = $1; BCond bc; bc.nnum = $2-1; bc.dofnum = 6;
          bc.val = 5.670400E-8*$3*$4*$5*$5*$5*$5; bc.type = BCond::Radiation; $$->add(bc); }
        ;
HelmHoltzBC:
        HSOMMERFELD NewLine SommerfeldBCDataList
	;
SommerfeldBCDataList:
        SommerfeldBC_Data
        | SommerfeldBCDataList SommerfeldBC_Data
	;
SommerfeldBC_Data:
        Integer Integer Float NewLine
          { domain->addSommer(new LineSommerBC($1-1, $2-1)); }
        | Integer Integer Integer Float NewLine
          { domain->addSommer(new TriangleSommerBC($1-1,$2-1,$3-1)); }
        | Integer Integer Integer Integer Float NewLine
          { domain->addSommer(new QuadSommerBC($1-1,$2-1,$3-1, $4-1)); }
	;
EleHelmHoltzBC:
        ELHSOMMERFELD NewLine SommerElement
        | EleHelmHoltzBC SommerElement
        ;
SommerElement:
        Integer Integer SommNodeNums NewLine
        { domain->addSommerElem($1-1, $2, 1.0, $3.num, $3.nd); 
          /*geoSource->addElem($1-1, $2, $3.num, $3.nd);include Sommer nodes in PackedEset -JF*/
        }
        ;
SommNodeNums:
        Integer
        { $$.num = 1; $$.nd[0] = $1-1; }
        | SommNodeNums Integer
        { if($$.num == 64) return -1;
          $$.nd[$$.num] = $2-1; $$.num++; }
        ;
Scatterer:
        SCATTERER NewLine ScattererEleList
        ;
ScattererEleList:
        ScattererEle
        | ScattererEleList ScattererEle
        ;
ScattererEle:
        Integer Integer Float NewLine
        { domain->addScatter(new LineSommerBC($1-1, $2-1));
          domain->addNeum(new LineSommerBC($1-1, $2-1)); }
        | Integer Integer Integer Float NewLine
        { domain->addScatter(new TriangleSommerBC($1-1,$2-1,$3-1));
          domain->addNeum(new TriangleSommerBC($1-1,$2-1,$3-1)); }
        | Integer Integer Integer Integer Float NewLine
        { domain->addScatter(new QuadSommerBC($1-1,$2-1,$3-1, $4-1));
          domain->addNeum(new QuadSommerBC($1-1,$2-1,$3-1, $4-1)); }
        ;
EleScatterer:
        ELSCATTERER NewLine ScatterElement
        | EleScatterer ScatterElement
        ;
ScatterElement:
        Integer Integer SommNodeNums NewLine
        { domain->addScatterElem($1-1, $2, 1.0, $3.num, $3.nd);
          domain->addNeumElem($1-1, $2, 1.0, $3.num, $3.nd); }
        ;
NeumScatterer:
        HNBO NewLine NeumElement
        | NeumScatterer NeumElement
        ;
NeumElement:
        Integer Integer SommNodeNums NewLine
        { domain->addNeumElem($1-1, $2, 1.0, $3.num, $3.nd); }
        ;
WetScatterer:
        HWIBO NewLine WetInterfaceElement
        | WetScatterer WetInterfaceElement
        ;
WetInterfaceElement:
        Integer Integer SommNodeNums NewLine
        { domain->addWetElem($1-1, $2, 1.0, $3.num, $3.nd); 
          domain->solInfo().isCoupled = true; 
          domain->solInfo().isMatching = true; }
        ;
HelmScatterer:
        HSCBO NewLine HelmScattererElement
        | HelmScatterer HelmScattererElement
        ;
HelmScattererElement:
        Integer Integer SommNodeNums NewLine
        { domain->addScatterElem($1-1, $2, 1.0, $3.num, $3.nd);}
        ;
PBC_Data:
	Integer Float NewLine
	{ $$.nnum = $1-1; $$.dofnum = 7; $$.val = $2; }
	;
PBCDataList:
	PBC_Data
	{ $$ = new BCList; $$->add($1); }
	| PBCDataList PBC_Data
	{ $$ = $1; $$->add($2); }
	;
AtdDirScatterer:
        ATDDIR NewLine PBCDataList
        { for(int i=0; i<$3->n; ++i) $3->d[i].type = BCond::Atddir; $$ = $3; }
        ;
AtdNeuScatterer:
        ATDNEU NewLine PBCDataList
        { for(int i=0; i<$3->n; ++i) { $3->d[i].type = BCond::Atdneu; } $$ = $3; }
        | ATDNEU Integer NewLine PBCDataList
        { for(int i=0; i<$4->n; ++i) { $4->d[i].type = BCond::Atdneu; } $$ = $4; }
        ;
AtdArbScatterer:
	ATDARB Float NewLine
        { domain->solInfo().ATDARBFlag = $2;}
        | AtdArbScatterer SommerElement
		;
AtdNeumScatterer:
	ATDDNB Float NewLine
        { domain->solInfo().ATDDNBVal = $2;}
        | AtdNeumScatterer NeumElement
        ;
AtdRobinScatterer:
	ATDROB Float Float Float NewLine
        { domain->solInfo().ATDROBVal = $2;
          domain->solInfo().ATDROBalpha = $3;
          domain->solInfo().ATDROBbeta = $4;}
        | AtdRobinScatterer HelmScattererElement
        ;
FarFieldPattern:
        FFP Integer NewLine
        { domain->setFFP($2); }
        | FFP Integer Integer NewLine
        { domain->setFFP($2,$3); }
        ;
FarFieldPatternDirs:
        FFPDIR Integer NewLine FFPDirList
        {
           domain->setFFP($2);
        }
        ;
ReadModeInfo:
	Integer EIGEN FNAME Integer NewLine
        { domain->solInfo().readInModes[$1] = ModalParams(ModalParams::Eigen, $3, $4); }
        | Integer ROBTYPE FNAME Integer NewLine
        { domain->solInfo().readInModes[$1] = ModalParams($2, $3, $4); }
        | Integer EIGEN FNAME Integer Float NewLine
        { domain->solInfo().readInModes[$1] = ModalParams(ModalParams::Eigen, $3, $4, $5); }
        | Integer ROBTYPE FNAME Integer Float NewLine
        { domain->solInfo().readInModes[$1] = ModalParams($2, $3, $4, $5); }
	;
Mode:
        READMODE FNAME NewLine
        { domain->solInfo().readInModes[0] = ModalParams(ModalParams::Undefined, $2); }
        | READMODE FNAME Integer NewLine
        { domain->solInfo().readInModes[0] = ModalParams(ModalParams::Undefined, $2, $3); }
	| READMODE NewLine ReadModeInfo
	| Mode ReadModeInfo
        | Mode ADJOINTBASIS FNAME Integer STRESSID NewLine
        { domain->solInfo().adjointMap[(OutputInfo::Type)$5] = domain->solInfo().readInAdjointROB.size();
          domain->solInfo().readInAdjointROB.push_back($3);
          domain->solInfo().maxSizeAdjointBasis.push_back($4); }
	;
IDisp:
        IDIS NewLine
        { }
        | IDIS ZERO NewLine
        { domain->solInfo().zeroInitialDisp = 1; }
	| IDIS NewLine BCDataList
	{ for(int i=0; i<$3->n; ++i) $3->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDis($3->n,$3->d) < 0)  return -1; }
	| IDisp MODAL NewLine ModalValList
	{ for(int i=0; i<$4->n; ++i) $4->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDisModal($4->n, $4->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; }
        | IDisp MODAL Integer NewLine ModalValList
        { for(int i=0; i<$5->n; ++i) $5->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDisModal($5->n, $5->d) < 0) return -1;
          domain->solInfo().modalCalled = true;
          domain->solInfo().idis_modal_id = $3; }
	;
IDisp6:
	IDIS6 Float NewLine
	{ $$ = new BCList; amplitude = $2;  }
	| IDIS6 NewLine
	{ $$ = new BCList; amplitude = 1.0; }
	| IDisp6 Integer Float Float Float Float Float Float NewLine
        { BCond bc; /* add 6 boundary conditions */
          bc.type = BCond::Idisp6;
          bc.nnum = $2-1; bc.dofnum = 0; bc.val = amplitude*$3; $$->add(bc);
                          bc.dofnum = 1; bc.val = amplitude*$4; $$->add(bc);
                          bc.dofnum = 2; bc.val = amplitude*$5; $$->add(bc);
                          bc.dofnum = 3; bc.val = amplitude*$6; $$->add(bc);
                          bc.dofnum = 4; bc.val = amplitude*$7; $$->add(bc);
                          bc.dofnum = 5; bc.val = amplitude*$8; $$->add(bc);
          geoSource->setIDis6($$->n, $$->d);
        }
        | IDisp6 Integer Float Float Float NewLine
        { BCond bc; /* add 6 boundary conditions */
          bc.type = BCond::Idisp6;
          bc.nnum = $2-1; bc.dofnum = 0; bc.val = amplitude*$3; $$->add(bc);
                          bc.dofnum = 1; bc.val = amplitude*$4; $$->add(bc);
                          bc.dofnum = 2; bc.val = amplitude*$5; $$->add(bc);
                          bc.dofnum = 3; bc.val = 0.0         ; $$->add(bc);
                          bc.dofnum = 4; bc.val = 0.0         ; $$->add(bc);
                          bc.dofnum = 5; bc.val = 0.0         ; $$->add(bc);
          geoSource->setIDis6($$->n, $$->d);
        }
	| GEPS NewLine
	{ fprintf(stderr," ... Geometric Pre-Stress Effects   ... \n"); 
          domain->solInfo().setGEPS(); }
	| BUCKLE NewLine
	{ domain->solInfo().buckling = 1; }
	;
IDisp6Pita:
        /* PITA: Get initial seed displacement for each time-slice. */
        PITADISP6 Integer NewLine
        { $$ = new BCList; PitaTS = $2; }
        | IDisp6Pita Integer Float Float Float Float Float Float NewLine
        { BCond bc;                          /* add 6 boundary conditions */
          bc.nnum = $2-1; bc.dofnum = 0; bc.val = $3; $$->add(bc);
                          bc.dofnum = 1; bc.val = $4; $$->add(bc);
                          bc.dofnum = 2; bc.val = $5; $$->add(bc);
                          bc.dofnum = 3; bc.val = $6; $$->add(bc);
                          bc.dofnum = 4; bc.val = $7; $$->add(bc);
                          bc.dofnum = 5; bc.val = $8; $$->add(bc);
          geoSource->setPitaIDis6($$->n, $$->d, PitaTS);
        } 
        ;
IVel6Pita:
        /* PITA: Get initial seed velocity for each time-slice. */
        PITAVEL6 Integer NewLine
        { $$ = new BCList; PitaTS = $2; }
        | IVel6Pita Integer Float Float Float Float Float Float NewLine
        { BCond bc;                          /* add 6 boundary conditions */
          bc.nnum = $2-1; bc.dofnum = 0; bc.val = $3; $$->add(bc);
                          bc.dofnum = 1; bc.val = $4; $$->add(bc);
                          bc.dofnum = 2; bc.val = $5; $$->add(bc);
                          bc.dofnum = 3; bc.val = $6; $$->add(bc);
                          bc.dofnum = 4; bc.val = $7; $$->add(bc);
                          bc.dofnum = 5; bc.val = $8; $$->add(bc);
          geoSource->setPitaIVel6($$->n, $$->d, PitaTS);
        }
        ;
IVel:
        IVEL NewLine
        { }
        | IVEL NewLine BCDataList
        { for(int i=0; i<$3->n; ++i) $3->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVel($3->n,$3->d) < 0) return -1; }
        | IVel MODAL NewLine ModalValList
        { for(int i=0; i<$4->n; ++i) $4->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVelModal($4->n, $4->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; }
        | IVel MODAL Integer NewLine ModalValList
        { for(int i=0; i<$5->n; ++i) $5->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVelModal($5->n, $5->d) < 0) return -1;
          domain->solInfo().modalCalled = true;
          domain->solInfo().ivel_modal_id = $3; }
	;
ITemp:
        ITEMP NewLine TBCDataList
        { for(int i=0; i<$3->n; ++i) $3->d[i].type = BCond::Itemperatures;
          if(geoSource->setIDis($3->n,$3->d) < 0) return -1; }
        | ITemp MODAL NewLine ModalValList
        { for(int i=0; i<$4->n; ++i) $4->d[i].type = BCond::Itemperatures;
          if(geoSource->setIDisModal($4->n, $4->d) < 0) return -1;
          domain->solInfo().modalCalled = true; }
	;
ETemp:
        ETEMP NewLine TBCDataList
        { domain->solInfo().setGEPS();
          for(int i=0; i<$3->n; ++i) $3->d[i].type = BCond::Etemperatures;
          if(geoSource->setIDis6($3->n,$3->d) < 0) return -1; }
        ;
NeumanBC:
        FORCE NewLine
        { $$ = new BCList; }
        | FORCE Integer NewLine
        { $$ = new BCList($2); }
        | NeumanBC BC_Data
        { $2.type = BCond::Forces; $2.loadsetid = $$->loadsetid; $$->add($2); }
        | NeumanBC Integer THRU Integer Integer Float NewLine
        { for(int i=$2; i<=$4; ++i) { BCond bc; bc.setData(i-1, $5-1, $6, BCond::Forces, $$->loadsetid); $$->add(bc); } }
        | NeumanBC Integer THRU Integer STEP Integer Integer Float NewLine
        { for(int i=$2; i<=$4; i+=$6) { BCond bc; bc.setData(i-1, $7-1, $8, BCond::Forces, $$->loadsetid); $$->add(bc); } }
        | NeumanBC Integer THRU Integer Integer Float MOMENTTYPE NewLine
        { for(int i=$2; i<=$4; ++i) { BCond bc; bc.setData(i-1, $5-1, $6, BCond::Forces, $$->loadsetid, (BCond::MomentType) $7); $$->add(bc); } }
        | NeumanBC Integer THRU Integer STEP Integer Integer Float MOMENTTYPE NewLine
        { for(int i=$2; i<=$4; i+=$6) { BCond bc; bc.setData(i-1, $7-1, $8, BCond::Forces, $$->loadsetid, (BCond::MomentType) $7); $$->add(bc); } }
        | NeumanBC SURF BC_Data
        { BCond *surf_bc = new BCond[1];
          surf_bc[0] = $3;
          surf_bc[0].type = BCond::Forces;
          surf_bc[0].loadsetid = $$->loadsetid;
          geoSource->addSurfaceNeuman(1,surf_bc);
          if(geoSource->getNumSurfaceNeuman() > 1) delete [] surf_bc; }
        ;
ModalNeumanBC:
        FORCE NewLine MODAL NewLine ModalValList
        { for(int i=0; i<$5->n; ++i) $5->d[i].type = BCond::Forces;
          if(geoSource->setNeumanModal($5->n, $5->d) < 0) return -1; }
        | FORCE Integer NewLine MODAL NewLine ModalValList
        { for(int i=0; i<$6->n; ++i) { $6->d[i].type = BCond::Forces; $6->d[i].loadsetid = $2; }
          if(geoSource->setNeumanModal($6->n, $6->d) < 0) return -1; }
        | FORCE NewLine MODAL Integer NewLine ModalValList
        { for(int i=0; i<$6->n; ++i) $6->d[i].type = BCond::Forces;
          if(geoSource->setNeumanModal($6->n, $6->d) < 0) return -1; }
        | FORCE Integer NewLine MODAL Integer NewLine ModalValList
        { for(int i=0; i<$7->n; ++i) { $7->d[i].type = BCond::Forces; $7->d[i].loadsetid = $2; }
          if(geoSource->setNeumanModal($7->n, $7->d) < 0) return -1; }
        ;
BCDataList:
	BC_Data
	{ $$ = new BCList; $$->add($1); }
	| BCDataList BC_Data
	{ $$ = $1; $$->add($2); }
        | Integer THRU Integer Integer Float NewLine
        { $$ = new BCList; for(int i=$1; i<=$3; ++i) { BCond bc; bc.setData(i-1, $4-1, $5); $$->add(bc); }} 
        | BCDataList Integer THRU Integer Integer Float NewLine
        { $$ = $1; for(int i=$2; i<=$4; ++i) { BCond bc; bc.setData(i-1, $5-1, $6); $$->add(bc); }} 
        | Integer THRU Integer STEP Integer Integer Float NewLine
        { $$ = new BCList; for(int i=$1; i<=$3; i+=$5) { BCond bc; bc.setData(i-1, $6-1, $7); $$->add(bc); } }
        | BCDataList Integer THRU Integer STEP Integer Integer Float NewLine
        { $$ = $1; for(int i=$2; i<=$4; i+=$6) { BCond bc; bc.setData(i-1, $7-1, $8); $$->add(bc); } }
	;
ModalValList:
	ModalVal
	{ $$ = new BCList; $$->add($1); }
	| ModalValList ModalVal
	{ $$ = $1; $$->add($2); }
	;
TBCDataList:
	TBC_Data
	{ $$ = new BCList; $$->add($1); }
	| TBCDataList TBC_Data
	{ $$ = $1; $$->add($2); }
	;
YMTTable:
        YMTT NewLine
        | YMTT NewLine YMTTList
        ;
YMTTList:
        CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($2); $$->add($4, $5); domain->addYMTT($$);}
        | YMTTList Float Float NewLine
        { $$->add($2, $3); }
        | YMTTList CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($3); $$->add($5, $6); domain->addYMTT($$);}
        ;
TETTable:
        TETT NewLine
        | TETT NewLine TETTList
        ;
TETTList:
        CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($2); $$->add($4, $5); domain->addCTETT($$);}
        | TETTList Float Float NewLine
        { $$->add($2, $3); }
        | TETTList CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($3); $$->add($5, $6); domain->addCTETT($$);}
        ;
SS1DTable:
        SS1DT NewLine
        | SS1DT NewLine SS1DTList
        ;
SS1DTList:
        CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($2); $$->add($4, $5); domain->addSS1DT($$);}
        | SS1DTList Float Float NewLine
        { $$->add($2, $3); }
        | SS1DTList CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($3); $$->add($5, $6); domain->addSS1DT($$);}
        ;
SS2DTable:
        SS2DT NewLine
        | SS2DT NewLine SS2DTList
        ;
SS2DTList:
        SURF Integer NewLine FloatList NewLine Float FloatList NewLine
        { $$ = new SS2DTData($2, $4, true); $$->add($6, $7); domain->addSS2DT($$); }
        | SURF Integer ENGINEERING SWITCH NewLine FloatList NewLine Float FloatList NewLine
        { $$ = new SS2DTData($2, $6, bool($4)); $$->add($8, $9); domain->addSS2DT($$); }
        | SS2DTList Float FloatList NewLine
        { $$->add($2, $3); }
        | SS2DTList SURF Integer NewLine FloatList NewLine Float FloatList NewLine
        { $$ = new SS2DTData($3, $5, true); $$->add($7, $8); domain->addSS2DT($$); }
        | SS2DTList SURF Integer ENGINEERING SWITCH NewLine FloatList NewLine Float FloatList NewLine
        { $$ = new SS2DTData($3, $7, bool($5)); $$->add($9, $10); domain->addSS2DT($$); }
        ;
YSSTable:
        YSST NewLine
        | YSST NewLine YSSTList
        ;
YSSTList:
        CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($2); $$->add($4, $5); domain->addYSST($$);}
        | YSSTList Float Float NewLine
        { $$->add($2, $3); }
        | YSSTList CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($3); $$->add($5, $6); domain->addYSST($$);}
        ;
YSSRTable:
        YSSRT NewLine
        | YSSRT NewLine YSSRTList
        ;
YSSRTList:
        CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($2); $$->add($4, $5); domain->addYSSRT($$);}
        | YSSRTList Float Float NewLine
        { $$->add($2, $3); }
        | YSSRTList CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($3); $$->add($5, $6); domain->addYSSRT($$);}
        ;
YMSTable:
        YMST NewLine
        | YMST NewLine YMSTList
        ;
YMSTList:
        CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($2); $$->add($4, $5); domain->addYMST($$);}
        | YMSTList Float Float NewLine
        { $$->add($2, $3); }
        | YMSTList CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($3); $$->add($5, $6); domain->addYMST($$);}
        ;
SDETAFTable:
        SDETAFT NewLine
        | SDETAFT NewLine SDETAFList
        ;
SDETAFList:
        CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($2); $$->add($4, $5); domain->addSDETAFT($$);}
        | SDETAFList Float Float NewLine
        { $$->add($2, $3); }
        | SDETAFList CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($3); $$->add($5, $6); domain->addSDETAFT($$);}
        ;
RUBDAFTable:
        RUBDAFT NewLine
        | RUBDAFT NewLine RUBDAFList
        ;
RUBDAFList:
        TABLE Integer NewLine Float Float Float Float Float NewLine
        { 
#ifdef USE_EIGEN3
          $$ = new GenMFTTData<Eigen::Vector4d>($2); $$->add($4, Eigen::Vector4d($5,$6,$7,$8)); domain->addRUBDAFT($$);
#else
          std::cerr << " *** ERROR: RUBDAFT command requires AERO-S configured with Eigen library. Exiting...\n"; exit(-1);
#endif
        }
        | RUBDAFList Float Float Float Float Float NewLine
        {
#ifdef USE_EIGEN3
          $$->add($2, Eigen::Vector4d($3,$4,$5,$6));
#else
          std::cerr << " *** ERROR: RUBDAFT command requires AERO-S configured with Eigen library. Exiting...\n"; exit(-1);
#endif
        }
        | RUBDAFList TABLE Integer NewLine Float Float Float Float Float NewLine
        {
#ifdef USE_EIGEN3
          $$ = new GenMFTTData<Eigen::Vector4d>($3); $$->add($5, Eigen::Vector4d($6,$7,$8,$9)); domain->addRUBDAFT($$);
#else
          std::cerr << " *** ERROR: RUBDAFT command requires AERO-S configured with Eigen library. Exiting...\n"; exit(-1);
#endif
        }
        ;
LMPConstrain:
        LMPC NewLine
        | LMPC NewLine MPCList
        | LMPC Float Float Float NewLine
        { domain->solInfo().xLMPCFactor = $2;
          domain->solInfo().yLMPCFactor = $3;
          domain->solInfo().zLMPCFactor = $4; }
	;
ModalLMPConstrain:
        LMPC NewLine MODAL Integer NewLine
        { domain->solInfo().localDualBasisSize.push_back($4);
          domain->solInfo().modalLMPC = true; }
        | ModalLMPConstrain Float NewLine
        { geoSource->pushBackROMLMPCVec($2); }
        
MPCList:
        MPCHeader MPCLine
        { $$ = $1;
          $$->addterm($2);
          domain->addLMPC($$); }
        | MPCList MPCLine
        { $$->addterm($2); }
        | MPCList MPCHeader MPCLine
        { $$ = $2;
          $$->addterm($3);
          domain->addLMPC($$); }
	;
MPCHeader:
        Integer NewLine
        { $$ = new LMPCons($1, 0.0); 
          $$->setSource(mpc::Lmpc); }
        | Integer Float NewLine
        { $$ = new LMPCons($1, $2); 
          $$->setSource(mpc::Lmpc); }
        | Integer Float MODE Integer NewLine
        { $$ = new LMPCons($1, $2);
          $$->type = $4; 
          $$->setSource(mpc::Lmpc); }
        | Integer Float ConstraintOptionsData NewLine
        { $$ = new LMPCons($1, $2);
          $$->lagrangeMult = $3.lagrangeMult;
          $$->penalty = $3.penalty; 
          $$->setSource(mpc::Lmpc); }
        | Integer Float MODE Integer ConstraintOptionsData NewLine
        { $$ = new LMPCons($1, $2);
          $$->type = $4;
          $$->lagrangeMult = $5.lagrangeMult;
          $$->penalty = $5.penalty;
          $$->setSource(mpc::Lmpc); }
        ;
MPCLine:
        Integer Integer Float NewLine
        { if($3 == 0.0) {
            fprintf(stderr," *** WARNING: zero coefficient in LMPC\n");
            fprintf(stderr," ***          node %d dof %d\n",$1,$2);
          }
          $$ = new LMPCTerm();
          $$->nnum = $1-1;
          $$->dofnum = $2-1;
          $$->coef.r_value = $3;
          if($2 == 1){
            $$->coef.r_value /= domain->solInfo().xLMPCFactor;
          } else if ($2 == 2) {
            $$->coef.r_value /= domain->solInfo().yLMPCFactor;
          } else if ($2 == 3) {
            $$->coef.r_value /= domain->solInfo().zLMPCFactor;
          }
        }
	;
ComplexLMPConstrain:
        HLMPC NewLine
        | HLMPC NewLine ComplexMPCList
        ;
ComplexMPCList:
        ComplexMPCHeader ComplexMPCLine
        { $$ = new LMPCons($1.nnum,$1.reval,$1.imval,$2); domain->addLMPC($$); }
        | ComplexMPCList ComplexMPCLine
        { $$->addterm($2); }
        | ComplexMPCList ComplexMPCHeader ComplexMPCLine
        { $$ = new LMPCons($2.nnum,$2.reval,$2.imval,$3); domain->addLMPC($$); }
        ;
ComplexMPCHeader:
        Integer CRHS Float Float NewLine // PJSA: added CRHS keyword to eliminate conflict
        { $$.nnum=$1; $$.reval=$3; $$.imval=$4; }
        | Integer Float NewLine
        { $$.nnum=$1; $$.reval=$2; $$.imval=0.0; }
        | Integer NewLine
        { $$.nnum=$1; $$.reval=0.0; $$.imval=0.0; }
        ;
ComplexMPCLine:
        Integer Integer Float Float NewLine
        { if(($3==0.0) && ($4==0.0)) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",$1,$2);
          return -1;
          }
          else { $$ = new LMPCTerm(true); $$->nnum=($1-1); $$->dofnum=($2-1); $$->coef.c_value=DComplex($3, $4); }
        }
        | Integer Integer Float NewLine
        { if($3==0.0) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",$1,$2);
          return -1;
          }
          else { $$ = new LMPCTerm(true); $$->nnum=($1-1); $$->dofnum=($2-1); $$->coef.c_value=DComplex($3,0.0); }
        }
        ;
ComplexNeumanBC:
	HNEUMAN NewLine ComplexBCDataList
        { $$ = $3; }
        | HNEUMAN Integer NewLine ComplexBCDataList
        { for(int i=0; i<$4->n; ++i) $4->d[i].loadsetid = $2;
          $$ = $4; }
        ;
ComplexBCDataList:
	ComplexBC_Data
	{ $$ = new ComplexBCList; $$->add($1); }
        | ComplexBCDataList ComplexBC_Data
        { $$ = $1; $$->add($2); }
	;
Materials:
	MATERIALS NewLine MatData
	| Materials MatData
	;
MatData:
	Integer Float Float Float Float Float Float Float Float Float Float Float Float Float Float NewLine
	{ StructProp sp; 
	  sp.A = $2;  sp.E = $3;  sp.nu  = $4;  sp.rho = $5;
          sp.c = $6;  sp.k = $7;  sp.eh  = $8;  sp.P   = $9;  sp.Ta  = $10; 
          sp.Q = $11; sp.W = $12; sp.Ixx = $13; sp.Iyy = $14; sp.Izz = $15;
          geoSource->addMat( $1-1, sp );
        }
        | Integer Float Float Float Float Float Float Float Float Float Float Float Float Float Float DAMPING Float Float NewLine
        { StructProp sp;
          sp.A = $2;  sp.E = $3;  sp.nu  = $4;  sp.rho = $5;
          sp.c = $6;  sp.k = $7;  sp.eh  = $8;  sp.P   = $9;  sp.Ta  = $10;
          sp.Q = $11; sp.W = $12; sp.Ixx = $13; sp.Iyy = $14; sp.Izz = $15;
          switch($16) {
            case 1 : // RAYDAMP
              sp.betaDamp = $17; sp.alphaDamp = $18;
              break;
            case 2 : // STRDAMP
              sp.etaDamp = $17; sp.betaDamp = $18;
              break;
            case 3 : // RUBDAMP
              sp.eta_E = $17;
              if(sp.eta_E >= 0) sp.eta_mu = $18;
              sp.E0 = $3;
              sp.mu0 = $3/(2*(1+$4));
              break;
            default :
              return -1;
          }
          geoSource->addMat( $1-1, sp );
        }
	| Integer Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float NewLine
	{ StructProp sp; 
	  sp.A = $2;  sp.E = $3;  sp.nu  = $4;  sp.rho = $5;
          sp.c = $6;  sp.k = $7;  sp.eh  = $8;  sp.P   = $9;  sp.Ta  = $10; 
          sp.Q = $11; sp.W = $12; sp.Ixx = $13; sp.Iyy = $14; sp.Izz = $15;
	  sp.ymin = $16; sp.ymax = $17; sp.zmin = $18; sp.zmax = $19;
          geoSource->addMat( $1-1, sp );
        }
        | Integer Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float DAMPING Float Float NewLine
        { StructProp sp;
          sp.A = $2;  sp.E = $3;  sp.nu  = $4;  sp.rho = $5;
          sp.c = $6;  sp.k = $7;  sp.eh  = $8;  sp.P   = $9;  sp.Ta  = $10;
          sp.Q = $11; sp.W = $12; sp.Ixx = $13; sp.Iyy = $14; sp.Izz = $15;
          sp.ymin = $16; sp.ymax = $17; sp.zmin = $18; sp.zmax = $19;
          switch($20) {
            case 1 : // RAYDAMP
              sp.betaDamp = $21; sp.alphaDamp = $22;
              break;
            case 2 : // STRDAMP
              sp.etaDamp = $21; sp.betaDamp = $22;
              break;
            case 3 : // RUBDAMP
              sp.eta_E = $21; 
              if(sp.eta_E >= 0) sp.eta_mu = $22;
              sp.E0 = $3;
              sp.mu0 = $3/(2*(1+$4));
              break;
            default :
              return -1;
          }
          geoSource->addMat( $1-1, sp );
        }
        | Integer Float Float Float Float Float Float Float NewLine
        { StructProp sp;
          sp.A = $2; sp.E = $3; sp.nu = $4; sp.rho = $5;
          sp.c = $6; sp.k = $7; sp.eh = $8;
          geoSource->addMat( $1-1, sp ); 
        }
        | Integer Float Float Float Float Float Float Float Float Float Float Float Float NewLine
        { StructProp sp;  // this is for spring: GID Kx Ky Kz lx1 ...
          sp.A = $2;  sp.E = $3;  sp.nu  = $4;  sp.rho = $5;
          sp.c = $6;  sp.k = $7;  sp.eh  = $8;  sp.P   = $9;  sp.Ta  = $10;
          sp.Q = $11; sp.W = $12; sp.Ixx = $13;  
          geoSource->addMat( $1-1, sp );
        }
        | Integer Float Float Float Float Float Float Float Float Float Float Float Float DAMPING Float NewLine
        { StructProp sp;  // this is for spring with stiffness-proportional damping : GID Kx Ky Kz lx1 ...
          sp.A = $2;  sp.E = $3;  sp.nu  = $4;  sp.rho = $5;
          sp.c = $6;  sp.k = $7;  sp.eh  = $8;  sp.P   = $9;  sp.Ta  = $10;
          sp.Q = $11; sp.W = $12; sp.Ixx = $13;
          if($14 == 0) sp.betaDamp = $15; else return -1;
          geoSource->addMat( $1-1, sp );
        }
        | Integer Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float
          Float Float Integer Float Float Integer Integer Float Float Float NewLine
        { StructProp sp; // this is used for the reduced mesh file output in Rom.d/MeshOutput.C
                         // all properties relevant to structural nonlinear dynamics should be included
          sp.A = $2;  sp.E = $3;  sp.nu  = $4;  sp.rho = $5;
          sp.c = $6;  sp.k = $7;  sp.eh  = $8;  sp.P   = $9;  sp.Ta  = $10;
          sp.Q = $11; sp.W = $12; sp.Ixx = $13; sp.Iyy = $14; sp.Izz = $15;
          sp.ymin = $16; sp.ymax = $17; sp.zmin = $18; sp.zmax = $19;
          sp.betaDamp = $20; sp.alphaDamp = $21; 
          sp.lagrangeMult = bool($22); sp.penalty = $23; sp.initialPenalty = $24;
          sp.funtype = $25; sp.type = StructProp::PropType($26); sp.k1 = $27; sp.k2 = $28; sp.k3 = $29;
          geoSource->addMat( $1-1, sp );
        }
        | Integer AMAT Float Float Float Float Float Float Float Float Float NewLine
        { StructProp sp;
          sp.soundSpeed = complex<double>($3,0.0);
          sp.fp.PMLtype = int($4);
          sp.fp.gamma = $5;
          sp.fp.Rx = $6;
          sp.fp.Sx = $7;
          sp.fp.Ry = $8;
          sp.fp.Sy = $9;
          sp.fp.Rz = $10;
          sp.fp.Sz = $11;
          sp.type = StructProp::Fluid;
          geoSource->addMat( $1-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        }
        | Integer AMAT Float Float Float Float Float Float Float Float Float Float NewLine
        { StructProp sp;
          sp.soundSpeed = complex<double>($3,0.0);
          sp.rho = $4;
          sp.fp.PMLtype = int($5);
          sp.fp.gamma = $6;
          sp.fp.Rx = $7;
          sp.fp.Sx = $8;
          sp.fp.Ry = $9;
          sp.fp.Sy = $10;
          sp.fp.Rz = $11;
          sp.fp.Sz = $12;
          sp.type = StructProp::Fluid;
          geoSource->addMat( $1-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        }
        | Integer AMAT Float Float Float Float Float Float Float Float Float Float Float NewLine
        { StructProp sp;
          sp.soundSpeed = complex<double>($3,$4);
          sp.rho = $5;
          sp.fp.PMLtype = int($6);
          sp.fp.gamma = $7;
          sp.fp.Rx = $8;
          sp.fp.Sx = $9;
          sp.fp.Ry = $10;
          sp.fp.Sy = $11;
          sp.fp.Rz = $12;
          sp.fp.Sz = $13;
          sp.type = StructProp::Fluid;
          geoSource->addMat( $1-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        }
        | Integer AMAT Float Float NewLine
        { StructProp sp;
          sp.soundSpeed = complex<double>($3,0.0);
          sp.rho = $4;
          sp.type = StructProp::Fluid;
          geoSource->addMat( $1-1, sp );
          domain->solInfo().acoustic = true;
        }
        | Integer AMAT Float Float Float NewLine
        { StructProp sp;
          sp.soundSpeed = complex<double>($3,$4);
          sp.rho = $5;
          sp.type = StructProp::Fluid;
          geoSource->addMat( $1-1, sp );
          domain->solInfo().acoustic = true;
        }
	| Integer FABMAT Integer Float Float Float Float Float Float Float Float Float Integer Integer Integer NewLine
        { StructProp sp;
          sp.F_op = $3;
          sp.E = $4;
          sp.rho = $5;
	  sp.A = $6;
          sp.F_Uc = $7;
          sp.F_Uf = $8;
          sp.lambda = $9;
          sp.F_h = $10;
          sp.F_d = $11;
          sp.F_dlambda = $12;
          sp.F_np = $13;
          sp.F_Nf = $14;
	  sp.Seed = $15;
          sp.type = StructProp::Fabric;
          geoSource->addMat( $1-1, sp );
        }
        | Integer THERMMAT Float Float Float Float Float Float Float Float Float NewLine
        { StructProp sp; 
          sp.A = $3; sp.rho = $4; sp.Q = $5; sp.c = $6; 
          sp.sigma = $7; sp.k = $8; sp.eh = $9; sp.P = $10; sp.Ta = $11;
          sp.type = StructProp::Thermal;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT NewLine
        { // rigid element or joint with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ConstraintOptionsData NewLine
        { // rigid element or joint
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer MASS Float Float Float Float Float Float Float Float Float Float NewLine
        { //           m     Ixx   Iyy   Izz   Ixy   Iyz   Ixz   cx    cy    cz
          // discrete mass with offset
          StructProp sp;
          sp.rho = $3;
          sp.Ixx = $4;
          sp.Iyy = $5;
          sp.Izz = $6;
          sp.Ixy = $7;
          sp.Iyz = $8;
          sp.Ixz = $9;
          sp.cx  = $10;
          sp.cy  = $11;
          sp.cz  = $12;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT MASS Float NewLine
        { // rigid 8-node brick element with mass, and default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = $4;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ConstraintOptionsData MASS Float NewLine
        { // rigid 8-node brick element with mass
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.rho = $5;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT MASS Float Float NewLine
        { // rigid beam or shell element with mass, and default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = $4;
          sp.A = sp.eh = $5;
          geoSource->addMat( $1-1, sp );
        } 
        | Integer CONSTRMAT ConstraintOptionsData MASS Float Float NewLine
        { // rigid beam or shell element with mass
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.rho = $5;
          sp.A = sp.eh = $6;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT Integer Float Float Float Float Float Float Integer NewLine
        { // constraint function element XXX deprecated format
          StructProp sp;
          sp.lagrangeMult = bool($3);
          sp.initialPenalty = sp.penalty = $4;
          sp.amplitude = $5;
          sp.omega = $6;
          sp.phase = $7;
          sp.B = $8;
          sp.C = $9;
          sp.relop = $10;
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT MODE Integer NewLine
        { // constraint function element with default constraint options
          StructProp sp;
          sp.amplitude = 0;
          sp.omega = 0;
          sp.phase = 0;
          sp.B = 0;
          sp.C = 0;
          sp.relop = $4;
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ConstraintOptionsData MODE Integer NewLine
        { // constraint function element
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.amplitude = 0;
          sp.omega = 0;
          sp.phase = 0;
          sp.B = 0;
          sp.C = 0;
          sp.relop = $5;
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ELEMENTARYFUNCTIONTYPE Float Float Float Float NewLine
        { // joint-with-driver, 2-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = $3;
          sp.amplitude = $4;
          sp.offset = $5;
          sp.c1 = $6;
          sp.c2 = $7;
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ELEMENTARYFUNCTIONTYPE Float Float Float Float Float NewLine
        { // joint-with-driver, 3-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = $3;
          sp.amplitude = $4;
          sp.offset = $5;
          sp.c1 = $6;
          sp.c2 = $7;
          sp.c3 = $8;
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ELEMENTARYFUNCTIONTYPE Float Float Float Float Float Float NewLine
        { // joint-with-driver, 4-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = $3;
          sp.amplitude = $4;
          sp.offset = $5;
          sp.c1 = $6;
          sp.c2 = $7;
          sp.c3 = $8;
          sp.c4 = $9;
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ConstraintOptionsData ELEMENTARYFUNCTIONTYPE Float Float Float Float NewLine
        { // joint-with-driver, 2-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.funtype = $4;
          sp.amplitude = $5;
          sp.offset = $6;
          sp.c1 = $7;
          sp.c2 = $8;
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ConstraintOptionsData ELEMENTARYFUNCTIONTYPE Float Float Float Float Float NewLine
        { // joint-with-driver, 3-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.funtype = $4;
          sp.amplitude = $5;
          sp.offset = $6;
          sp.c1 = $7;
          sp.c2 = $8;
          sp.c3 = $9;
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ConstraintOptionsData ELEMENTARYFUNCTIONTYPE Float Float Float Float Float Float NewLine
        { // joint-with-driver, 4-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.funtype = $4;
          sp.amplitude = $5;
          sp.offset = $6;
          sp.c1 = $7;
          sp.c2 = $8;
          sp.c3 = $9;
          sp.c4 = $10;
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ELEMENTARYFUNCTIONTYPE Float Float Float Float SPRINGMAT Float NewLine
        { // actuated joint, 2-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = $3;
          sp.amplitude = $4;
          sp.offset = $5;
          sp.c1 = $6;
          sp.c2 = $7;
          sp.k1 = $9;
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ELEMENTARYFUNCTIONTYPE Float Float Float Float Float SPRINGMAT Float NewLine
        { // actuated joint, 3-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = $3;
          sp.amplitude = $4;
          sp.offset = $5;
          sp.c1 = $6;
          sp.c2 = $7;
          sp.c3 = $8;
          sp.k1 = $10;
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ELEMENTARYFUNCTIONTYPE Float Float Float Float Float Float SPRINGMAT Float NewLine
        { // actuated joints, 4-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = $3;
          sp.amplitude = $4;
          sp.offset = $5;
          sp.c1 = $6;
          sp.c2 = $7;
          sp.c3 = $8;
          sp.c4 = $9;
          sp.k1 = $11;
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ConstraintOptionsData ELEMENTARYFUNCTIONTYPE Float Float Float Float SPRINGMAT Float NewLine
        { // actuated joint, 2-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.funtype = $4;
          sp.amplitude = $5;
          sp.offset = $6;
          sp.c1 = $7;
          sp.c2 = $8;
          sp.k1 = $10;
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ConstraintOptionsData ELEMENTARYFUNCTIONTYPE Float Float Float Float Float SPRINGMAT Float NewLine
        { // actuated joint, 3-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.funtype = $4;
          sp.amplitude = $5;
          sp.offset = $6;
          sp.c1 = $7;
          sp.c2 = $8;
          sp.c3 = $9;
          sp.k1 = $11;
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ConstraintOptionsData ELEMENTARYFUNCTIONTYPE Float Float Float Float Float Float SPRINGMAT Float NewLine
        { // actuated joints, 4-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.funtype = $4;
          sp.amplitude = $5;
          sp.offset = $6;
          sp.c1 = $7;
          sp.c2 = $8;
          sp.c3 = $9;
          sp.c4 = $10;
          sp.k1 = $12;
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT SPRINGMAT Float NewLine
        { // RevoluteJointSpringCombo with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = $4;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT SPRINGMAT Float FREEPLAY FreeplayProps NewLine
        { // RevoluteJointSpringComboWithFreeplay with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = $4;
          sp.freeplay[0] = $6;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ConstraintOptionsData SPRINGMAT Float NewLine
        { // RevoluteJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = $5;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ConstraintOptionsData SPRINGMAT Float FREEPLAY FreeplayProps NewLine
        { // RevoluteJointSpringComboWithFreeplay
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = $5;
          sp.freeplay[0] = $7;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT SPRINGMAT Float Float NewLine
        { // UniversalJointSpringCombo with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = $4;
          sp.k2 = $5;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT SPRINGMAT Float Float FREEPLAY FreeplayProps FreeplayProps NewLine
        { // UniversalJointSpringComboWithFreeplay with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = $4;
          sp.k2 = $5;
          sp.freeplay[0] = $7;
          sp.freeplay[1] = $8;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ConstraintOptionsData SPRINGMAT Float Float NewLine
        { // UniversalJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = $5;
          sp.k2 = $6;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ConstraintOptionsData SPRINGMAT Float Float FREEPLAY FreeplayProps FreeplayProps NewLine
        { // UniversalJointSpringComboWithFreeplay
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = $5;
          sp.k2 = $6;
          sp.freeplay[0] = $8;
          sp.freeplay[1] = $9;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT SPRINGMAT Float Float Float NewLine
        { // SphericalJointSpringCombo with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = $4;
          sp.k2 = $5;
          sp.k3 = $6;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT SPRINGMAT Float Float Float FREEPLAY FreeplayProps FreeplayProps FreeplayProps NewLine
        { // SphericalJointSpringComboWithFreeplay with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = $4;
          sp.k2 = $5;
          sp.k3 = $6;
          sp.freeplay[0] = $8;
          sp.freeplay[1] = $9;
          sp.freeplay[2] = $10;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ConstraintOptionsData SPRINGMAT Float Float Float NewLine
        { // SphericalJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = $5;
          sp.k2 = $6;
          sp.k3 = $7;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT ConstraintOptionsData SPRINGMAT Float Float Float FREEPLAY FreeplayProps FreeplayProps FreeplayProps NewLine
        { // SphericalJointSpringComboWithFreeplay
          StructProp sp;
          sp.lagrangeMult = $3.lagrangeMult;
          sp.initialPenalty = sp.penalty = $3.penalty;
          sp.constraint_hess = $3.constraint_hess;
          sp.constraint_hess_eps = $3.constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = $5;
          sp.k2 = $6;
          sp.k3 = $7;
          sp.freeplay[0] = $9;
          sp.freeplay[1] = $10;
          sp.freeplay[2] = $11;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer SPRINGMAT Float NewLine
        { // TorsionalSpring or TranslationalSpring (types 200,201,202)
          StructProp sp;
          sp.k1 = $3;
          sp.rho = 0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer SPRINGMAT Float FREEPLAY Float NewLine
        { // 1-sided TorsionalSpring or TranslationalSpring with freeplay (types 203,204,205)
          StructProp sp;
          sp.k1 = $3;
          sp.rho = 0;
          sp.freeplay[0].ul = $5;
          sp.freeplay[0].dz = 0.0;
          sp.freeplay[0].uz = 1.0;
          geoSource->addMat( $1-1, sp );
        }
        | Integer ARUBBERMAT Float Float Float Float Float Float Float Float Float NewLine
        { // Acoustic rubber
          StructProp sp;
          sp.E0 = sp.E = $3; sp.dE = $4/(2*M_PI); sp.eta_E = $5; sp.deta_E = $6/(2*M_PI);
          sp.mu0 = $7; sp.dmu = $8/(2*M_PI); sp.eta_mu = $9; sp.deta_mu = $10/(2*M_PI);
          sp.rho = $11;
          geoSource->addMat( $1-1, sp );
        }
        | Integer ARUBBERMAT Float NewLine
        { // Acoustic rubber
          StructProp sp;
          sp.E0 = sp.E = $3;
          geoSource->addMat( $1-1, sp );
        }
	;
FreeplayProps:
      /*Float NewLine
        { // 1-parameter freeplay model
          $$.ll = -$1;
          $$.ul =  $1;
          $$.lz = 1.0;
          $$.dz = 0.0;
          $$.uz = 1.0;
        }*/
        Float Float Float Float Float
        { // 5-parameter freeplay model
          $$.ll = $1;
          $$.ul = $2;
          $$.lz = $3;
          $$.dz = $4;
          $$.uz = $5;
        }
ElemSet:
	TOPOLOGY NewLine Element
	| ElemSet Element
	;
FaceSet:
	SURFACETOPOLOGY Integer NewLine 
	{ if($2 == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          $$ = new SurfaceEntity($2);
          $$->SetReverseNormals(false);
          domain->AddSurfaceEntity($$);
        }
        | SURFACETOPOLOGY Integer REVERSENORMALS NewLine 
        { if($2 == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          $$ = new SurfaceEntity($2);
          $$->SetReverseNormals(true);
          domain->AddSurfaceEntity($$);
        }
        | SURFACETOPOLOGY Integer SHELLTHICKNESS Float NewLine 
        { if($2 == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          $$ = new SurfaceEntity($2);
          $$->SetIsShellFace(true);
          $$->SetShellThickness($4);
          domain->AddSurfaceEntity($$);
        }
        | SURFACETOPOLOGY Integer SHELLTHICKNESS Float REVERSENORMALS NewLine 
        { if($2 == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          $$ = new SurfaceEntity($2);
          $$->SetIsShellFace(true);
          $$->SetShellThickness($4);
          $$->SetReverseNormals(true);
          domain->AddSurfaceEntity($$);
        }
        | FaceSet Integer Integer NodeNums NewLine 
        { if($$->GetReverseNormals()) { // reverse the node numbering
            int *nodes = new int[$4.num];
            for(int i=0; i<$4.num; ++i) nodes[$4.num-1-i] = $4.nd[i];
            $$->AddFaceElement($2-1, $3, $4.num, nodes);
            delete [] nodes;
          }
          else $$->AddFaceElement($2-1, $3, $4.num, $4.nd);
        }
	; 
MortarCondition:
	MORTARTIED Integer Integer NewLine
	{ $$ = new MortarHandler($2, $3); domain->AddMortarCond($$); }
	| MORTARTIED Integer Integer STDMORTAR NewLine
	{ $$ = new MortarHandler($2, $3); domain->AddMortarCond($$); }
	| MORTARTIED Integer Integer DUALMORTAR NewLine
	{ $$ = new MortarHandler($2, $3); $$->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond($$); 
        }
	| MORTARTIED Integer Integer SEARCHTOL Float NewLine
	{ $$ = new MortarHandler($2, $3, $5); domain->AddMortarCond($$); }
	| MORTARTIED Integer Integer SEARCHTOL Float Float NewLine
	{ $$ = new MortarHandler($2, $3, $5, $6); domain->AddMortarCond($$); }
	| MORTARTIED Integer Integer STDMORTAR SEARCHTOL Float NewLine
	{ $$ = new MortarHandler($2, $3, $6); domain->AddMortarCond($$); }
	| MORTARTIED Integer Integer STDMORTAR SEARCHTOL Float Float NewLine
	{ $$ = new MortarHandler($2, $3, $6, $7); domain->AddMortarCond($$); }
	| MORTARTIED Integer Integer DUALMORTAR SEARCHTOL Float NewLine
	{ $$ = new MortarHandler($2, $3, $6); $$->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond($$); 
        }
	| MORTARTIED Integer Integer DUALMORTAR SEARCHTOL Float Float NewLine
	{ $$ = new MortarHandler($2, $3, $6, $7); $$->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond($$); 
        }
	;
WetInterface:
        WETINTERFACE Integer Integer NewLine
        { domain->addWetInterface($2, $3); domain->solInfo().isCoupled = true; }
        | WETINTERFACE Integer NewLine
        { domain->addWetInterface($2, $2); 
          domain->solInfo().isCoupled  = true; 
          domain->solInfo().isMatching = true; }
        ;
TiedSurfaces:
        // $2 = pair id, $3 = master surface, $4 = slave surface, $5 = mortar type, $6 = normal search tolerace, $7 = tangential search tolerance
        // $8 = number of iterations (TD enforcement), $9 = convergence tolerance (TD enforcement), $10 = friction coefficient
        TIEDSURFACES NewLine
        { }
        | TiedSurfaces Integer Integer Integer NewLine
        {
          $$ = new MortarHandler($3, $4);
          $$->SetInteractionType(MortarHandler::TIED);
          domain->AddMortarCond($$); 
        }
        | TiedSurfaces Integer Integer Integer Integer NewLine
        { 
          $$ = new MortarHandler($3, $4); 
          $$->SetInteractionType(MortarHandler::TIED);
          $$->SetMortarType($5);
          domain->AddMortarCond($$); 
        }
        | TiedSurfaces Integer Integer Integer Integer Float NewLine
        {
          $$ = new MortarHandler($3, $4, $6);
          $$->SetInteractionType(MortarHandler::TIED);
          $$->SetMortarType($5);
          domain->AddMortarCond($$);
        }
        | TiedSurfaces Integer Integer Integer Integer Float Float NewLine
        {
          $$ = new MortarHandler($3, $4, $6, $7);
          $$->SetInteractionType(MortarHandler::TIED);
          $$->SetMortarType($5);
          domain->AddMortarCond($$);
        }
        | TiedSurfaces Integer Integer Integer Integer Float Float Integer Float NewLine
        {
          $$ = new MortarHandler($3, $4, $6, $7);
          $$->SetInteractionType(MortarHandler::TIED);
          $$->SetMortarType($5);
          $$->SetTDEnfParams($8, $9);
          domain->AddMortarCond($$);
        }
        | TiedSurfaces Integer Integer Integer ConstraintOptionsData NewLine
        {
          $$ = new MortarHandler($3, $4);
          $$->SetInteractionType(MortarHandler::TIED);
          $$->SetConstraintOptions($5);
          domain->AddMortarCond($$);
        }
        | TiedSurfaces Integer Integer Integer Integer ConstraintOptionsData NewLine
        {
          $$ = new MortarHandler($3, $4);
          $$->SetInteractionType(MortarHandler::TIED);
          $$->SetMortarType($5);
          $$->SetConstraintOptions($6);
          domain->AddMortarCond($$);
        }
        | TiedSurfaces Integer Integer Integer Integer Float ConstraintOptionsData NewLine
        {
          $$ = new MortarHandler($3, $4, $6);
          $$->SetInteractionType(MortarHandler::TIED);
          $$->SetMortarType($5);
          $$->SetConstraintOptions($7);
          domain->AddMortarCond($$);
        }
        | TiedSurfaces Integer Integer Integer Integer Float Float ConstraintOptionsData NewLine
        {
          $$ = new MortarHandler($3, $4, $6, $7);
          $$->SetInteractionType(MortarHandler::TIED);
          $$->SetMortarType($5);
          $$->SetConstraintOptions($8);
          domain->AddMortarCond($$);
        }
        | TiedSurfaces Integer Integer Integer Integer Float Float Integer Float ConstraintOptionsData NewLine
        {
          $$ = new MortarHandler($3, $4, $6, $7);
          $$->SetInteractionType(MortarHandler::TIED);
          $$->SetMortarType($5);
          $$->SetTDEnfParams($8, $9);
          $$->SetConstraintOptions($10);
          domain->AddMortarCond($$);
        }
        ;
FSInterface:
        FSINTERFACE NewLine
        { }
        | FSInterface Integer Integer Integer NewLine
        { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface($3, $4); domain->solInfo().isCoupled = true; 
          if($3 == $4) domain->solInfo().isMatching = true;
        }
        | FSInterface Integer Integer Integer Float Float NewLine
        { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface($3, $4, $5, $6); domain->solInfo().isCoupled = true;
          if($3 == $4) domain->solInfo().isMatching = true;
        }
        ;
HEVibInfo:
        HEFSB NewLine HEVInterfaceElement
        { domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; }
        | HEVibInfo HEVInterfaceElement
        ;
HEVInterfaceElement:
        Integer Integer SommNodeNums NewLine
        { domain->addWetElem($1-1, $2, 1.0, $3.num, $3.nd);
          domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; }
        ;
HEInterface:
        HEINTERFACE NewLine
        { }
        | HEInterface Integer Integer Integer NewLine
        { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface($3, $4); domain->solInfo().HEV = 1;
          if($3 == $4) domain->solInfo().isMatching = true;
        }
        | HEInterface Integer Integer Integer Float Float NewLine
        { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface($3, $4, $5, $6); domain->solInfo().HEV = 1;
          if($3 == $4) domain->solInfo().isMatching = true;
        }
        ;
ContactSurfaces:
        // $2 = pair id, $3 = master surface, $4 = slave surface, $5 = mortar type, $6 = normal search tolerace, $7 = tangential search tolerance
        // $8 = number of iterations (TD enforcement), $9 = convergence tolerance (TD enforcement), $10 = friction coefficient
        CONTACTSURFACES NewLine
        { }
        | CONTACTSURFACES MODE Integer NewLine
        { domain->solInfo().contactsurface_mode = $3; }
        | ContactSurfaces ContactSurfacesInfo NewLine
        { domain->AddMortarCond($2); }
        | ContactSurfaces ContactSurfacesInfo ConstraintOptionsData NewLine
        { $2->SetConstraintOptions($3); domain->AddMortarCond($2); }
        | ContactSurfaces ContactSurfacesInfo ConstraintOptionsData MODE Integer NewLine
        { $2->SetConstraintOptions($3); $2->SetCtcMode($5); domain->AddMortarCond($2); }
        | ContactSurfaces ContactSurfacesInfo MODE Integer NewLine
        { $2->SetCtcMode($4); domain->AddMortarCond($2); }
        ;
ContactSurfacesInfo:
        // $2 = pair id, $3 = master surface, $4 = slave surface, $5 = mortar type, $6 = normal search tolerace, $7 = tangential search tolerance
        // $8 = number of iterations (TD enforcement), $9 = convergence tolerance (TD enforcement), $10 = friction coefficient
        TRIVIAL
        { domain->solInfo().trivial_detection = true;
          $$ = new MortarHandler(0, 0);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetMortarType(MortarHandler::STD);
          $$->SetCtcMode(domain->solInfo().contactsurface_mode); }
        | Integer Integer Integer
        {
          $$ = new MortarHandler($2, $3);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetMortarType(MortarHandler::STD);
          $$->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
        | Integer Integer Integer Integer
        {
          $$ = new MortarHandler($2, $3);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetMortarType($4);
          $$->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
        | Integer Integer Integer Integer Float
        {
          $$ = new MortarHandler($2, $3, $5);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetMortarType($4);
          $$->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
        | Integer Integer Integer Integer Float Float
        {
          $$ = new MortarHandler($2, $3, $5, $6);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetMortarType($4);
          $$->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
        | Integer Integer Integer Integer Float Float Integer Float
        { /* frictionless */
          $$ = new MortarHandler($2, $3, $5, $6);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetTDEnfParams($7, $8);
          $$->SetMortarType($4);
          $$->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
        | Integer Integer Integer Integer Float Float Integer Float Float
        { /* constant friction */
          $$ = new MortarHandler($2, $3, $5, $6);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetTDEnfParams($7, $8);
          $$->SetFrictionCoef($9);
          $$->SetMortarType($4);
          $$->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
        | Integer Integer Integer Integer Float Float Integer Float Float Float Float
        { /* velocity dependent friction */
          $$ = new MortarHandler($2, $3, $5, $6);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetTDEnfParams($7, $8);
          $$->SetFrictionCoef($9, $10, $11);
          $$->SetMortarType($4);
          $$->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
        | Integer Integer Integer Integer Float Float Integer Float Float Float Float Float
        { /* pressure dependent friction */
          $$ = new MortarHandler($2, $3, $5, $6);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetTDEnfParams($7, $8);
          $$->SetFrictionCoef($9, $10, $11, $12);
          $$->SetMortarType($4);
          $$->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
        ;
Parameters:
        PARAMETERS NewLine
        | Parameters ACMECNTL Integer NewLine
        { domain->solInfo().dist_acme = $3; }
        | Parameters NOSECONDARY NewLine
        { domain->solInfo().no_secondary = true; }
        | Parameters NOGHOSTING NewLine
        { domain->solInfo().no_ghosting = true; }
        | Parameters SHELLSIMPLELOFTING NewLine
        { domain->solInfo().shell_simple_lofting = true; }
        | Parameters NOMULTIPLEINTERACTIONS NewLine
        { domain->solInfo().no_multiple_interactions = true; }
        | Parameters SHARPNONSHARPANGLE Float NewLine
        { domain->solInfo().sharp_non_sharp_angle = $3; }
        | Parameters NONORMALSMOOTHING NewLine
        { domain->solInfo().normal_smoothing = false; }
        | Parameters NORMALSMOOTHINGDISTANCE Float NewLine
        { domain->solInfo().normal_smoothing_distance = $3; }
        | Parameters RESOLUTIONMETHOD Integer NewLine
        { domain->solInfo().resolution_method = $3; }
        | Parameters OLDDYNAMICSEARCH NewLine
        { domain->solInfo().old_dynamic_search = true; }
        | Parameters PARTITIONGAP NewLine
        { domain->solInfo().partition_gap = true; }
        | Parameters DEFAULTPENALTY NewLine
        { domain->solInfo().default_penalty = true; }
        | Parameters GLOBALSEARCHCULL NewLine
        { domain->solInfo().global_search_cull = true; }
        | Parameters NOWARPEDVOLUME NewLine
        { domain->solInfo().no_warped_volume = true; }
        | Parameters AUTOTOL NewLine
        { domain->solInfo().auto_tol = true; }
        | Parameters AGRESSIVETOLERANCES NewLine
        { domain->solInfo().auto_tol = true;
          domain->solInfo().agressive_tolerances = true; }
        | Parameters SKIPPHYSICALFACES NewLine
        { domain->solInfo().skip_physical_faces = true; }
        | Parameters FFIDEBUG Integer NewLine
        { domain->solInfo().ffi_debug = bool($3); }
        | Parameters MORTARSCALING Float NewLine
        { domain->solInfo().mortar_scaling = $3; }
        | Parameters MORTARINTEGRATIONRULE Integer NewLine
        { domain->solInfo().mortar_integration_rule = $3; }
        | Parameters ANDESCLR Float NewLine
        { domain->solInfo().andes_clr = $3; }
        | Parameters ANDESCQR Float NewLine
        { domain->solInfo().andes_cqr = $3; }
        | Parameters ANDESBETAB Float NewLine
        { domain->solInfo().andes_betab = $3; }
        | Parameters ANDESALPHA Float NewLine
        { domain->solInfo().andes_alpha = $3; }
        | Parameters ANDESBETAM Float NewLine
        { domain->solInfo().andes_betam = $3; }
        | Parameters NLMEMPTYPE Integer NewLine
        { domain->solInfo().nlmembrane_pressure_type = $3; }
NodeSet:
	NODETOKEN NewLine Node
	{ geoSource->addNode($3.num, $3.xyz, $3.cp, $3.cd); }
        | NODETOKEN Float Float Float NewLine Node
        { domain->solInfo().scalePosCoords = true;
          domain->solInfo().xScaleFactor = $2;
          domain->solInfo().yScaleFactor = $3;
          domain->solInfo().zScaleFactor = $4;
          geoSource->addNode($6.num, $6.xyz, $6.cp, $6.cd); }
	| NodeSet Node
	{ geoSource->addNode($2.num, $2.xyz, $2.cp, $2.cd); }
	;
Node:
	Integer Float Float Float NewLine
	{ $$.num = $1-1; $$.xyz[0] = $2; $$.xyz[1] = $3;  $$.xyz[2] = $4;  $$.cp = 0;  $$.cd = 0; }
	| Integer Float Float NewLine
	{ $$.num = $1-1; $$.xyz[0] = $2; $$.xyz[1] = $3;  $$.xyz[2] = 0.0; $$.cp = 0;  $$.cd = 0; }
	| Integer Float NewLine
	{ $$.num = $1-1; $$.xyz[0] = $2; $$.xyz[1] = 0.0; $$.xyz[2] = 0.0; $$.cp = 0;  $$.cd = 0; }
        | Integer Float Float Float Integer Integer NewLine
        { $$.num = $1-1; $$.xyz[0] = $2; $$.xyz[1] = $3;  $$.xyz[2] = $4;  $$.cp = $5; $$.cd = $6;
          if($5 != 0) domain->solInfo().basicPosCoords = false;
          if($6 != 0) domain->solInfo().basicDofCoords = false; }
        | Integer Float Float Float Integer NewLine
        { $$.num = $1-1; $$.xyz[0] = $2; $$.xyz[1] = $3;  $$.xyz[2] = $4;  $$.cp = $5; $$.cd = $5;
          if($5 != 0) { domain->solInfo().basicPosCoords = false; domain->solInfo().basicDofCoords = false; } }
	;
Element:
	Integer Integer NodeNums NewLine
	{ /* Define each Element */
          geoSource->addElem($1-1, $2, $3.num, $3.nd);}
	;
NodeNums:
	Integer
	{ $$.num = 1; $$.nd[0] = $1-1;}
	| NodeNums Integer
	{ if($$.num == 500) return -1; 
          $$.nd[$$.num] = $2-1; $$.num++;}
	;
BC_Data:
	Integer Integer Float NewLine
	{ $$.nnum = $1-1; $$.dofnum = $2-1; $$.val = $3; $$.mtype = BCond::Axial; }
	| Integer Integer NewLine
	{ $$.nnum = $1-1; $$.dofnum = $2-1; $$.val = 0.0; $$.mtype = BCond::Axial; }
        | Integer Integer Float MOMENTTYPE NewLine
        { $$.nnum = $1-1; $$.dofnum = $2-1; $$.val = $3; $$.mtype = (BCond::MomentType) $4; }
        | Integer Integer MOMENTTYPE NewLine
        { $$.nnum = $1-1; $$.dofnum = $2-1; $$.val = 0.0; $$.mtype = (BCond::MomentType) $3; }
	;
ModalVal:
	Integer Float NewLine
	{ $$.nnum = $1-1;  $$.dofnum = -1;  $$.val = $2; }
	;
TBC_Data:
	Integer Float NewLine
	{ $$.nnum = $1-1; $$.dofnum = 6; $$.val = $2; }
	;
ComplexBC_Data:
	Integer Integer Float Float NewLine
	{ $$.nnum = $1-1; $$.dofnum = $2-1; $$.reval = $3; $$.imval = $4;  }
	| Integer Integer Float NewLine
	{ $$.nnum = $1-1; $$.dofnum = $2-1; $$.reval = $3; $$.imval = 0.0; }
	;
ConstrainedSurfaceFrameDList:
        CSFRAMES NewLine
        | ConstrainedSurfaceFrameDList Frame
        { geoSource->setCSFrame($2.num,$2.d); }
        ;
FrameDList:
        EFRAMES NewLine
	| FrameDList Frame
	{ geoSource->setFrame($2.num,$2.d);  }
	;
Frame:
	Integer Float Float Float Float Float Float Float Float Float NewLine
	{ $$.num = $1-1; 
          $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7;
          $$.d[6] = $8; $$.d[7] = $9; $$.d[8] = $10; }
        | Integer THIRDNODE Integer NewLine
        { $$.num = $1-1;
          geoSource->makeEframe($1-1, $3, $$.d); }
	;
NodalFrameDList:
        NFRAMES NewLine
        | NodalFrameDList NodalFrame
        { geoSource->setNodalFrame($2.id,$2.o,$2.d,$2.type); }
        ;
NodalFrame:
        Integer Float Float Float Float Float Float Float Float Float NewLine
        { $$.id = $1;
          $$.type = NFrameData::Rectangular;
          $$.o[0] = 0;  $$.o[1] = 0;  $$.o[2] = 0;
          $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7;
          $$.d[6] = $8; $$.d[7] = $9; $$.d[8] = $10; }
        | Integer Float Float Float Float Float Float Float Float Float Float Float Float NewLine
        { $$.id = $1;
          $$.type = NFrameData::Rectangular;
          $$.o[0] = $2;  $$.o[1] = $3;  $$.o[2] = $4;
          $$.d[0] = $5;  $$.d[1] = $6;  $$.d[2] = $7;
          $$.d[3] = $8;  $$.d[4] = $9;  $$.d[5] = $10;
          $$.d[6] = $11; $$.d[7] = $12; $$.d[8] = $13; }
        | Integer FRAMETYPE Float Float Float Float Float Float Float Float Float NewLine
        { $$.id = $1;
          $$.type = $2;
          $$.o[0] = 0;  $$.o[1] = 0;   $$.o[2] = 0;
          $$.d[0] = $3; $$.d[1] = $4;  $$.d[2] = $5;
          $$.d[3] = $6; $$.d[4] = $7;  $$.d[5] = $8;
          $$.d[6] = $9; $$.d[7] = $10; $$.d[8] = $11; }
        | Integer FRAMETYPE Float Float Float Float Float Float Float Float Float Float Float Float NewLine
        { $$.id = $1;
          $$.type = $2;
          $$.o[0] = $3;  $$.o[1] = $4;  $$.o[2] = $5;
          $$.d[0] = $6;  $$.d[1] = $7;  $$.d[2] = $8;
          $$.d[3] = $9;  $$.d[4] = $10; $$.d[5] = $11;
          $$.d[6] = $12; $$.d[7] = $13; $$.d[8] = $14; }
        ;
BoffsetList:
	BOFFSET NewLine
	| BoffsetList Integer Integer Float Float Float NewLine
	{ OffsetData od;
	  od.first = $2-1; od.last = $3-1;
	  od.o[0] = $4; od.o[1] = $5; od.o[2] = $6; 
	  geoSource->addOffset(od); }
	;
Attributes:
	ATTRIBUTES NewLine 
        { $$ = 0; }
        | ATTRIBUTES Integer NewLine
        { geoSource->setLocalIndex($2-1); }
        /* define hyper reduction coefficient for element with no attribute */
        | Attributes Integer HRC Float NewLine
        { geoSource->setElementLumpingWeight($2-1,$4);
          domain->solInfo().elemLumpPodRom = true; }
        | Attributes Integer HRC Float EXTFOL NewLine
        { geoSource->setElementLumpingWeight($2-1,$4);
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
        /* define attribute and hyper reduction coefficient for non-composite element */
        | Attributes Integer Integer NewLine
        { geoSource->setAttrib($2-1,$3-1); }
        | Attributes Integer Integer HRC Float NewLine
        { geoSource->setAttrib($2-1,$3-1); 
          geoSource->setElementLumpingWeight($2-1,$5);
          domain->solInfo().elemLumpPodRom = true; }
        | Attributes Integer Integer HRC REDFOL Float NewLine /* deprecated syntax */
        { geoSource->setAttrib($2-1,$3-1);
          geoSource->setElementLumpingWeight($2-1,$6);
          domain->solInfo().elemLumpPodRom = true; 
          domain->solInfo().reduceFollower = true; }
        | Attributes Integer Integer HRC Float EXTFOL NewLine /* new syntax */
        { geoSource->setAttrib($2-1,$3-1);
          geoSource->setElementLumpingWeight($2-1,$5);
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
        /* define attribute and hyper reduction coefficient for composite element with CFRAME */
        | Attributes Integer Integer Integer Integer NewLine
        { geoSource->setAttrib($2-1,$3-1,$4-1,$5-1); }
        | Attributes Integer Integer Integer Integer HRC Float NewLine
        { geoSource->setAttrib($2-1,$3-1,$4-1,$5-1);
          geoSource->setElementLumpingWeight($2-1,$7); 
          domain->solInfo().elemLumpPodRom = true; }
        | Attributes Integer Integer Integer Integer HRC Float EXTFOL NewLine
        { geoSource->setAttrib($2-1,$3-1,$4-1,$5-1);
          geoSource->setElementLumpingWeight($2-1,$7);   
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
        /* define attribute and hyper reduction coefficient for composite element with reference angle  */
        | Attributes Integer Integer Integer THETA Float NewLine
        { geoSource->setAttrib($2-1,$3-1,$4-1,-2,$6); }
        | Attributes Integer Integer Integer THETA Float HRC Float NewLine
        { geoSource->setAttrib($2-1,$3-1,$4-1,-2,$6);
          geoSource->setElementLumpingWeight($2-1,$8);
          domain->solInfo().elemLumpPodRom = true; }
        | Attributes Integer Integer Integer THETA Float HRC Float EXTFOL NewLine
        { geoSource->setAttrib($2-1,$3-1,$4-1,-2,$6); 
          geoSource->setElementLumpingWeight($2-1,$8);
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
        /* define attributes for a range of elements */
        | Attributes Integer Integer IDENTITY NewLine
        { int i;
          for(i=$2; i<$3+1; ++i)
            geoSource->setAttrib(i-1,i-1);
        }
        | Attributes Integer Integer Integer NewLine
        { int i;
          for(i=$2; i<$3+1; ++i)
            geoSource->setAttrib(i-1,$4-1);
        }
        | Attributes Integer Integer Integer Integer Integer NewLine
        { int i;
          for(i=$2; i<$3+1; ++i)
            geoSource->setAttrib(i-1, $4-1, $5-1, $6-1);
        }
        | Attributes Integer Integer Integer Integer THETA Float NewLine
        { int i;
          for(i=$2; i<$3+1; ++i)
            geoSource->setAttrib(i-1, $4-1, $5-1, -2, $7);
        }
        ;
Ellump:
        ELLUMP NewLine
        { domain->solInfo().elemLumpPodRom = true;
          geoSource->setLocalIndex(0); }
        | ELLUMP Integer NewLine
        { domain->solInfo().elemLumpPodRom = true;
          geoSource->setLocalIndex($2-1); }
        | Ellump Integer Float NewLine
        { geoSource->setElementLumpingWeight($2 - 1, $3); }
        | Ellump REDFOL NewLine /* deprecated */
        { domain->solInfo().reduceFollower = true;}
        | Ellump EXTFOL NewLine
        { domain->solInfo().reduceFollower = true;}
        ;
LocalReducedOrderBases:
        LOCALREDUCEDORDERBASES NewLine
        | LocalReducedOrderBases LocalBasesAuxi
        | LocalReducedOrderBases LocalBasesCent
        ;
LocalBasesAuxi:
        AUXILIARY Integer Integer FNAME NewLine
        { domain->solInfo().readInLocalBasesAuxi[std::make_pair($2-1,$3-1)] = std::string($4); }
        ;
LocalBasesCent:
        CENTER FNAME NewLine
        { domain->solInfo().readInLocalBasesCent.push_back(std::string($2)); }
        ;
ReducedStiffness:
	REDUCEDSTIFFNESS NewLine
	{ domain->solInfo().ReducedStiffness = true;}
        | ReducedStiffness Float NewLine
        { geoSource->pushBackStiffVec($2);}

UDeimBasis:
	UDEIMBASIS NewLine
        | UDeimBasis Integer Integer NewLine
        { domain->solInfo().forcePodSize = $2;
          domain->solInfo().maxDeimBasisSize = $3;}
        | UDeimBasis Float NewLine
        { geoSource->pushBackUDEIMVec($2);}

SampNodeSlot:
        SAMPNODESLOT NewLine
        | SampNodeSlot Integer Integer NewLine
        { domain->solInfo().DEIMPodRom = true;
          geoSource->setSampleNodesAndSlots($2-1,$3);}
        | SampNodeSlot Integer Integer Integer NewLine
        { geoSource->setSampleElemsAndDOFs($3-1,$4);
          domain->solInfo().UDEIMPodRom = true;}
        ;

Pressure:
	PRESSURE NewLine
        { $$ = 0; }
        | PRESSURE Integer NewLine
        { $$ = $2; }
	| Pressure Integer Float NewLine
	{ PressureBCond pbc;
          pbc.setData($2-1, $3, $$, true);
          geoSource->setElementPressure(pbc); }
	| Pressure Integer Integer Float NewLine
	{ for(int i = $2; i < ($3+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, $4, $$, true);
            geoSource->setElementPressure(pbc);
          } }
        | Pressure SURF Integer Float NewLine
        { PressureBCond *pbc = new PressureBCond[1];
          pbc[0].setData($3-1, $4, $$, true);
          geoSource->addSurfacePressure(1, pbc);
          if(geoSource->getNumSurfacePressure() > 1) delete [] pbc; }
        | Pressure SURF Integer THRU Integer Float NewLine
        { for(int i = $3; i <= $5; ++i) {
            PressureBCond *pbc = new PressureBCond[1];
            pbc[0].setData(i-1, $6, $$, true);
            geoSource->addSurfacePressure(1, pbc);
            if(geoSource->getNumSurfacePressure() > 1) delete [] pbc;
          } }
        | Pressure Integer Float SWITCH NewLine
        { PressureBCond pbc;
          pbc.setData($2-1, $3, $$, $4);
          geoSource->setElementPressure(pbc); }
        | Pressure Integer Integer Float SWITCH NewLine
        { for(int i = $2; i < ($3+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, $4, $$, $5);
            geoSource->setElementPressure(pbc);
          } }
        | Pressure SURF Integer Float SWITCH NewLine
        { PressureBCond *pbc = new PressureBCond[1];
          pbc[0].setData($3-1, $4, $$, $5);
          geoSource->addSurfacePressure(1, pbc);
          if(geoSource->getNumSurfacePressure() > 1) delete [] pbc; }
        /* New option to set element pressure using "FACE" keyword */
        | Pressure Integer FACE Integer Float NewLine
        { PressureBCond pbc;
          pbc.setData($2-1, $5, $$, true);
          pbc.face = $4-1;
          geoSource->setElementPressure(pbc); }
        | Pressure Integer Integer FACE Integer Float NewLine
        { for(int i = $2; i < ($3+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, $6, $$, true);
            pbc.face = $5-1;
            geoSource->setElementPressure(pbc);
          } }
        | Pressure Integer FACE Integer Float SWITCH NewLine
        { PressureBCond pbc;
          pbc.setData($2-1, $5, $$, $6);
          pbc.face = $4-1;
          geoSource->setElementPressure(pbc); }
        | Pressure Integer Integer FACE Integer Float SWITCH NewLine
        { for(int i = $2; i < ($3+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, $6, $$, $7);
            pbc.face = $5-1;
            geoSource->setElementPressure(pbc);
          } }
	;
Lumped:
	LUMPED NewLine
	{ geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false); 
          geoSource->setConsistentPFlag(false); 
        }
        | LUMPED Integer NewLine
        { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false, $2);
          geoSource->setConsistentPFlag(false);
        }
/*      note: SolverInfo::inertialLumping is now set automatically
        | LUMPED Integer Integer NewLine
        { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false, $2);
          geoSource->setConsistentPFlag(false);
          domain->solInfo().inertiaLumping = $3;
        }*/
	;
MassAugmentation:
        MASSAUGMENTATION SWITCH NewLine
        { domain->solInfo().useMassAugmentation = $2; }
        ;
Preload:
        PRELOAD NewLine
        { }
        | Preload Integer Float NewLine
        { geoSource->setElementPreLoad( $2-1, $3 ); }
        | Preload Integer THRU Integer Float NewLine
        { int i;
          for(i=$2; i<($4+1); ++i)
            geoSource->setElementPreLoad( i-1, $5 );
        }
        | Preload Integer Float Float Float NewLine
        { double load[3] = { $3, $4, $5 };
          geoSource->setElementPreLoad( $2-1, load ); }
        | Preload Integer THRU Integer Float Float Float NewLine
        { double load[3] = { $5, $6, $7 };
          int i;
          for(i=$2; i<($4+1); ++i)
            geoSource->setElementPreLoad( i-1, load );
        }
	;
Sensitivity:
        SENSITIVITY NewLine
        { domain->solInfo().sensitivity = true; }
        | Sensitivity SENSITIVITYMETHOD SENSITIVITYID NewLine
        { domain->solInfo().sensitivityMethod = (SolverInfo::SensitivityMethod) $3; }
        | Sensitivity READSENSITIVITY FNAME NewLine
        { domain->solInfo().readShapeSen = true;
          domain->solInfo().readInShapeSen = $3; }  
        | Sensitivity TOLSEN Float NewLine
        { domain->solInfo().sensitivityTol = $3; }
        | Sensitivity RELAXATIONSEN Float NewLine
        { domain->solInfo().qsMaxvelSen = $3; }
        | Sensitivity RATIOTOLSEN Float NewLine
        { domain->solInfo().ratioSensitivityTol = $3; }    
        | Sensitivity KSPARAM Float NewLine
        { domain->solInfo().ksParameter = $3; } 
        | Sensitivity KSMAX Float NewLine
        { domain->solInfo().ksMax = $3; } 
        | Sensitivity THICKNESSGROUPLIST ThicknessGroup NewLine
        { } 
        | Sensitivity STRESSNODELIST StressNode NewLine
        { } 
        | Sensitivity DISPNODELIST NewLine DispNode
        { } 
  ;
DispNode:
        Integer Integer NewLine 
        { domain->setDispNode($1-1, $2-1); }
        | DispNode Integer Integer NewLine
        { domain->setDispNode($2-1, $3-1); }
        | Integer Integer Integer NewLine 
        { domain->setDispNode($1-1, $2-1, $3-1); }
        | DispNode Integer Integer Integer NewLine
        { domain->setDispNode($2-1, $3-1, $4-1); }
        | Integer Integer Integer Integer NewLine 
        { domain->setDispNode($1-1, $2-1, $3-1, $4-1); }
        | DispNode Integer Integer Integer Integer NewLine
        { domain->setDispNode($2-1, $3-1, $4-1, $5-1); }
        | Integer Integer Integer Integer Integer NewLine 
        { domain->setDispNode($1-1, $2-1, $3-1, $4-1, $5-1); }
        | DispNode Integer Integer Integer Integer Integer NewLine
        { domain->setDispNode($2-1, $3-1, $4-1, $5-1, $6-1); }
        | Integer Integer Integer Integer Integer Integer NewLine 
        { domain->setDispNode($1-1, $2-1, $3-1, $4-1, $5-1, $6-1); }
        | DispNode Integer Integer Integer Integer Integer Integer NewLine
        { domain->setDispNode($2-1, $3-1, $4-1, $5-1, $6-1, $7-1); }
        | Integer Integer Integer Integer Integer Integer Integer NewLine 
        { domain->setDispNode($1-1, $2-1, $3-1, $4-1, $5-1, $6-1, $7-1); }
        | DispNode Integer Integer Integer Integer Integer Integer Integer NewLine
        { domain->setDispNode($2-1, $3-1, $4-1, $5-1, $6-1, $7-1, $8-1); }
  ;
StressNode:
        Integer 
        { domain->setStressNodes($1); }
        | StressNode Integer 
        { domain->setStressNodes($2); }
  ;
ThicknessGroup:
        Integer 
        { domain->setThicknessGroup($1); }
        | ThicknessGroup Integer 
        { domain->setThicknessGroup($2); }
  ;
Statics:
        STATS NewLine
        { domain->solInfo().setProbType(SolverInfo::Static); }
        | Statics Solver
        { domain->solInfo().solvercntl = $2; }
        | Statics CASES CasesList NewLine
        | Statics PIECEWISE NewLine
        { // activate piecewise constant configuration dependent external forces for a linear dynamic analysis
          domain->solInfo().piecewise = true;
        }
        | Statics PIECEWISE Float Float NewLine
        { // activate piecewise constant configuration dependent external forces for a linear static analysis
          domain->solInfo().piecewise = true;
          domain->solInfo().piecewise_dlambda = $3;
          domain->solInfo().piecewise_maxLambda = $4;
        }
        | Statics COUPLEDSCALE Float NewLine
        { domain->solInfo().coupled_scale = $3; }
        | Statics TURKEL Integer NewLine
        { domain->sommerfeldType = $3;
          domain->curvatureFlag = 0; }
        | Statics TURKEL Integer Float NewLine
        { domain->sommerfeldType = $3;
          domain->curvatureConst1 = $4;
          domain->curvatureFlag = 1; }
        | Statics TURKEL Integer Float Float NewLine
        { domain->sommerfeldType = $3;
          domain->curvatureConst1 = $4;
          domain->curvatureConst2 = $5;
          domain->curvatureFlag = 2; }
        | Statics DMPC SWITCH NewLine
        { domain->solInfo().dmpc = bool($3); }
        | Statics MPCCHECK Integer NewLine
        { domain->solInfo().dbccheck = bool($3); }
        ;
CasesList:
        Integer
        { domain->solInfo().loadcases.push_back($1); }
        | CasesList Integer
        { domain->solInfo().loadcases.push_back($2); }
        ;
SolverMethod:
	DIRECT NewLine
        { $$ = new SolverCntl(default_cntl);
          $$->type = SolverSelection::Direct;
          $$->subtype = 0; }
	| DIRECT Integer NewLine
	{ $$ = new SolverCntl(default_cntl);
          $$->type = SolverSelection::Direct;
          $$->subtype = $2; }
	| SOLVERTYPE NewLine
	{ $$ = new SolverCntl(default_cntl);
          $$->type = SolverSelection::Direct;
          $$->subtype = $1; }
    | SOLVERTYPE PIVOT NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Direct;
      $$->subtype = $1;
      $$->pivot = true; }
    | SOLVERTYPE SCALED NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Direct;
      $$->subtype = $1;
      $$->scaled = true; }
    | SOLVERTYPE UNSYMMETRIC NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Direct;
      $$->subtype = $1;
      $$->unsymmetric = true; }
    | ITERTYPE NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Iterative;
      $$->iterType = $1; }
    | ITERTYPE Integer NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Iterative;
      $$->iterType = $1;
      $$->precond = $2; }
    | ITERTYPE Integer Float NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Iterative;
      $$->iterType = $1;
      $$->precond = $2;
      $$->tol=$3; }
    | ITERTYPE Integer Float Integer NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Iterative;
      $$->iterType = $1;
      $$->precond = $2;
      $$->tol = $3;
      $$->maxit = $4; }
    | ITERTYPE Integer Float Integer Integer NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Iterative;
      $$->iterType = $1;
      $$->precond = $2;
      $$->tol = $3;
      $$->maxit = $4;
      $$->iterSubtype = $5; }
    | ITERTYPE Integer Float Integer Integer Integer NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Iterative;
      $$->iterType = $1;
      $$->precond = $2;
      $$->tol = $3;
      $$->maxit = $4;
      $$->iterSubtype = $5;
      $$->maxvecsize = $6; }
    | FETI NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Feti; }
    | FETI Integer NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Feti;
      $$->fetiInfo.version = (FetiInfo::Version) ($2-1); }
    | FETI Integer Float NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Feti;
      $$->fetiInfo.maxit = $2;
      $$->fetiInfo.tol = $3;
      $$->fetiInfo.maxortho = $2; }
    | FETI Integer Float Integer NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Feti;
      $$->fetiInfo.maxit = $2;
      $$->fetiInfo.tol = $3;
      $$->fetiInfo.maxortho = $4; }
    | FETI Integer FETI2TYPE NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Feti;
      $$->fetiInfo.version = (FetiInfo::Version) ($2-1);
      $$->fetiInfo.feti2version = (FetiInfo::Feti2Version) $3; }
	| FETI DP NewLine
	{ $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Feti;
      $$->fetiInfo.corners = FetiInfo::allCorners6;
      $$->fetiInfo.version = FetiInfo::fetidp; }
	| FETI LIB NewLine
	{ $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::FetiLib;
      $$->fetiInfo.corners = FetiInfo::allCorners6;
      $$->fetiInfo.version = FetiInfo::fetidp; }
    | FETI DPH NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Feti;
      $$->fetiInfo.corners = FetiInfo::allCorners6;
      $$->fetiInfo.version = FetiInfo::fetidp;
      $$->fetiInfo.dph_flag = true; }
    | HFETI NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Feti;
      $$->fetiInfo.scaling = FetiInfo::tscaling;
      $$->fetiInfo.corners = FetiInfo::allCorners3;
      $$->fetiInfo.version = FetiInfo::fetidp;
      $$->fetiInfo.dph_flag = true;
      $$->fetiInfo.augment = FetiInfo::Edges;
      $$->fetiInfo.rbmType = FetiInfo::None;
      $$->fetiInfo.nGs = 0; }
    | HFETI Integer Float NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Feti;
      $$->fetiInfo.maxit = $2;
      $$->fetiInfo.tol = $3;
      $$->fetiInfo.krylovtype = 1;
      $$->fetiInfo.scaling = FetiInfo::tscaling;
      $$->fetiInfo.corners = FetiInfo::allCorners6;
      $$->fetiInfo.version = FetiInfo::fetidp;
      $$->fetiInfo.dph_flag = true;
      $$->fetiInfo.augment = FetiInfo::Edges;
      $$->fetiInfo.rbmType = FetiInfo::None;
      $$->fetiInfo.nGs = 0; }
    | HFETI Integer Float Integer NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Feti;
      $$->fetiInfo.maxit = $2;
      $$->fetiInfo.tol = $3;
      $$->fetiInfo.numcgm = $4;
      $$->fetiInfo.krylovtype = 1;
      $$->fetiInfo.scaling = FetiInfo::tscaling;
      $$->fetiInfo.corners = FetiInfo::allCorners6;
      $$->fetiInfo.version = FetiInfo::fetidp;
      $$->fetiInfo.dph_flag = true;
      $$->fetiInfo.augment = FetiInfo::Edges;
      $$->fetiInfo.rbmType = FetiInfo::None;
      $$->fetiInfo.nGs = 0; }
    | HFETI Integer Float Integer Float NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Feti;
      $$->fetiInfo.maxit = $2;
      $$->fetiInfo.tol = $3;
      $$->fetiInfo.numcgm = $4;
      $$->fetiInfo.tolcgm = $5;
      $$->fetiInfo.spaceDimension = 2;
      $$->fetiInfo.krylovtype = 1;
      $$->fetiInfo.scaling = FetiInfo::tscaling;
      $$->fetiInfo.corners = FetiInfo::allCorners6;
      $$->fetiInfo.version = FetiInfo::fetidp;
      $$->fetiInfo.dph_flag = true;
      $$->fetiInfo.augment = FetiInfo::Edges;
      $$->fetiInfo.rbmType = FetiInfo::None;
      $$->fetiInfo.nGs = 0; }
    | HFETI Integer Float Integer Float Integer NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::Feti;
      $$->fetiInfo.maxit = $2;
      $$->fetiInfo.tol = $3;
      $$->fetiInfo.numcgm = $4;
      $$->fetiInfo.tolcgm = $5;
      $$->fetiInfo.spaceDimension = $6;
      $$->fetiInfo.krylovtype = 1;
      $$->fetiInfo.scaling = FetiInfo::tscaling;
      $$->fetiInfo.corners = FetiInfo::allCorners6;
      $$->fetiInfo.version = FetiInfo::fetidp;
      $$->fetiInfo.dph_flag = true;
      $$->fetiInfo.augment = FetiInfo::Edges;
      $$->fetiInfo.rbmType = FetiInfo::None;
      $$->fetiInfo.nGs = 0; }
    | BLOCKDIAG NewLine
    { $$ = new SolverCntl(default_cntl);
      $$->type = SolverSelection::BlockDiag; }
    | SOLVERHANDLE Integer NewLine
    { $$ = &domain->solInfo().solvercntls[$2]; }
    ;
Solvercntl:
        SOLVERCNTL Integer NewLine Solver
        { $$ = &domain->solInfo().solvercntls[$2];
          *$$ = *$4; }
Solver:
        SolverMethod
        { $$ = $1; }
        | Solver VERBOSE NewLine
        { $$->verbose = 1; }
        | Solver VERBOSE Integer NewLine
        { $$->verbose = $3; }
        | Solver PRINTNUMBER Integer NewLine
        { $$->fetiInfo.printNumber = $3; }
        | Solver TRBMlc Float NewLine
        { $$->trbm = $3; }
        | Solver SPARSERENUM Integer NewLine
        { $$->sparse_renum = $3; }
        | Solver SPARSEMAXSUP Integer NewLine
        { $$->sparse_maxsup = $3; }
        | Solver SPARSEDEFBLK Integer NewLine
        { $$->sparse_defblk = $3; }
        | Solver SPOOLESTAU Float NewLine
        { $$->spooles_tau = $3; }
        | Solver SPOOLESMAXSIZE Integer NewLine
        { $$->spooles_maxsize = $3; }
        | Solver SPOOLESMAXDOMAINSIZE Integer NewLine
        { if($3 < 0) {
            $3 = 24;
            fprintf(stderr," *** WARNING: spooles_maxdomainsize must be > 0,"
                           " using 24\n");
          }
          $$->spooles_maxdomainsize = $3; }
        | Solver SPOOLESSEED Integer NewLine
        { $$->spooles_seed = $3; }
        | Solver SPOOLESMAXZEROS Float NewLine
        { if(($3 < 0.0) || ($3 > 1.0)) {
            $3 = 0.04;
            fprintf(stderr," *** WARNING: spooles_maxzeros outside acceptable limits (0..1),"
                           " using 0.04\n");
          }
          $$->spooles_maxzeros = $3; }
        | Solver SPOOLESMSGLVL Integer NewLine
        { $$->spooles_msglvl = $3; }
        | Solver SPOOLESSCALE Integer NewLine
        { $$->spooles_scale = $3; }
        | Solver SPOOLESRENUM Integer NewLine
        { $$->spooles_renum = $3; }
	| Solver MUMPSICNTL Integer Integer NewLine
	{ $$->mumps_icntl[$3] = $4; }
	| Solver MUMPSCNTL Integer Float NewLine
	{ $$->mumps_cntl[$3] = $4; }
	| Solver MUMPSMINEQ Integer NewLine
	{ $$->mumps_mineq = $3; }
	| Solver MUMPSSTRIDE Integer NewLine
	{ $$->mumps_stride = $3; }
        | Solver GOLDFARBTOL Float NewLine
        { $$->goldfarb_tol = $3; }
        | Solver TOLPCG Float NewLine
        { $$->tol = $3; }
        | Solver SUBTYPE Integer NewLine
        { $$->iterSubtype = $3; }
        | Solver ILUDROPTOL Float NewLine
        { $$->ilu_droptol = $3; }
        | Solver MAXITR Integer NewLine
        { $$->maxit = $3; 
          $$->fetiInfo.maxit = $3; }
        | Solver MAXORTHO Integer NewLine
        { $$->maxvecsize = $3;
          $$->fetiInfo.maxortho = $3; }
        | Solver FETIPREC NewLine
        { $$->fetiInfo.precno = (FetiInfo::Preconditioner) $2; }
	| Solver PRECNO FETIPREC NewLine
        { $$->fetiInfo.precno = (FetiInfo::Preconditioner) $3; }
        | Solver PRECNO LUMPED NewLine
        { $$->fetiInfo.precno = FetiInfo::lumped; }
        | Solver PRECNO Integer NewLine
        { $$->precond = $3;
          if(isFeti($$->type) && (($3 < 0) || ($3 > 3))) {
            $3 = 1;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner selected, using lumped\n");
          }
          $$->fetiInfo.precno = (FetiInfo::Preconditioner) $3; }
        | Solver PRECTYPE Integer NewLine
        { if(($3 < 0) || ($3 > 1)) {
            $3 = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          $$->fetiInfo.prectype = (FetiInfo::PreconditionerType) $3; }
        | Solver PRECTYPE PRECTYPEID NewLine
        { if(($3 < 0) || ($3 > 1)) {
            $3 = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          $$->fetiInfo.prectype = (FetiInfo::PreconditionerType) $3; }
        | Solver TOLFETI Float NewLine
  	{ $$->fetiInfo.tol = $3; }
        | Solver TOLFETI Float Float NewLine
        { $$->fetiInfo.tol = $3; 
          $$->fetiInfo.absolute_tol = $4; }
        | Solver STAGTOL Float NewLine
        { $$->fetiInfo.stagnation_tol = $3; }
        | Solver STAGTOL Float Float NewLine
        { $$->fetiInfo.stagnation_tol = $3;
          $$->fetiInfo.absolute_stagnation_tol = $4; }
        | Solver PTOL Float Float NewLine
        { $$->fetiInfo.primal_proj_tol = $3;
          $$->fetiInfo.dual_proj_tol = $4; }
        | Solver PMAXIT Integer Integer NewLine
        { $$->fetiInfo.primal_plan_maxit = $3;
          $$->fetiInfo.dual_plan_maxit = $4; }
        | Solver PLANTOL Float Float NewLine
        { $$->fetiInfo.primal_plan_tol = $3;
          $$->fetiInfo.dual_plan_tol = $4; }
	| Solver NOCOARSE NewLine
	{ $$->fetiInfo.noCoarse = 1; }
	| Solver PROJ Integer NewLine
	{ if($3 == 1) 
            $$->fetiInfo.nonLocalQ = 0;
          else if($3 == 2) {
            $$->fetiInfo.nonLocalQ = 1;
            if($$->fetiInfo.version == FetiInfo::feti2) {
              $$->fetiInfo.nonLocalQ = 0;
              fprintf(stderr," *** WARNING: Basic projector is used"
                             " with FETI 2\n");
            }
          } else if($3 == 3) {
            $$->fetiInfo.nonLocalQ = 1;
            $$->fetiInfo.nQ = 3;
          } else if($3 == 4) {
            $$->fetiInfo.nonLocalQ = 1;
            $$->fetiInfo.nQ = 4;
          } else
            fprintf(stderr," *** WARNING: This projector does not exist,"
                           " using basic projector\n"); }
	| Solver SCALING Integer NewLine
	{ if(($3 < 0) || ($3 > 2)) $3 = 1; 
          $$->fetiInfo.scaling = (FetiInfo::Scaling) $3; }
	| Solver SCALING SCALINGTYPE NewLine
	{ $$->fetiInfo.scaling = (FetiInfo::Scaling) $3; }
        | Solver MPCSCALING Integer NewLine
        { if(($3 < 0) || ($3 > 2)) $3 = 2;
          $$->fetiInfo.mpc_scaling = (FetiInfo::Scaling) $3; }
        | Solver MPCSCALING SCALINGTYPE NewLine
        { $$->fetiInfo.mpc_scaling = (FetiInfo::Scaling) $3; }
        | Solver FSISCALING Integer NewLine
        { if(($3 < 0) || ($3 > 2)) $3 = 2;
          $$->fetiInfo.fsi_scaling = (FetiInfo::Scaling) $3; }
        | Solver FSISCALING SCALINGTYPE NewLine
        { $$->fetiInfo.fsi_scaling = (FetiInfo::Scaling) $3; }
        | Solver MPCELEMENT NewLine
        { $$->fetiInfo.mpc_element = true; }
        | Solver FSIELEMENT NewLine
        { $$->fetiInfo.fsi_element = true; }
        | Solver FSICORNER Integer NewLine
        { $$->fetiInfo.fsi_corner = $3; }
        | Solver NOLOCALFSISPLITING NewLine
        { $$->fetiInfo.splitLocalFsi = false; } 
        | Solver WETCORNERS NewLine
        { $$->fetiInfo.wetcorners = true; }
	| Solver CORNER CORNERTYPE NewLine
        { $$->fetiInfo.corners = (FetiInfo::CornerType) $3; }
        | Solver CORNER CORNERTYPE Integer NewLine
        { $$->fetiInfo.corners = (FetiInfo::CornerType) $3; 
          $$->fetiInfo.pick_unsafe_corners = bool($4); }
        | Solver CORNER AUGMENTTYPE NewLine
        { if($3 == 0) {
            $$->fetiInfo.corners = FetiInfo::noCorners;
            $$->fetiInfo.pickAnyCorner = 0; 
            $$->fetiInfo.bmpc = true;
            $$->fetiInfo.pick_unsafe_corners = false;
            $$->fetiInfo.augment = FetiInfo::none;
          } }
        | Solver AUGMENT AUGMENTTYPE NewLine
        { if($$->fetiInfo.dph_flag && ($3 == 1)) {
            std::cerr << "WARNING: Selected augment type is unsupported for FETI-DPH, set to EdgeGs \n";
            $$->fetiInfo.augment = FetiInfo::Edges;
          }
          else $$->fetiInfo.augment = (FetiInfo::AugmentType) $3;

          if($$->fetiInfo.augment == FetiInfo::Edges) {
            $$->fetiInfo.rbmType = FetiInfo::translation;
            $$->fetiInfo.nGs = 3;
          } 
          else if($$->fetiInfo.augment == FetiInfo::WeightedEdges) {
            $$->fetiInfo.rbmType = FetiInfo::translation;
            $$->fetiInfo.nGs = 3;
          }
          else if($$->fetiInfo.augment == FetiInfo::Gs) {
            $$->fetiInfo.rbmType = FetiInfo::all;
            $$->fetiInfo.nGs = 6;
          } }
        | Solver AUGMENT AUGMENTTYPE RBMSET NewLine
        { if($$->fetiInfo.dph_flag && ($3 == 1)) {
            std::cerr << "WARNING: Selected augment type is unsupported for FETI-DPH, set to EdgeGs \n";
            $$->fetiInfo.augment = FetiInfo::Edges;
          }
          else $$->fetiInfo.augment = (FetiInfo::AugmentType) $3;
          if($$->fetiInfo.dph_flag && ($4 > 2) && ($4 < 6)) {
            std::cerr << "WARNING: Selected rbm type is unsupported for FETI-DPH, set to translation \n";
            $$->fetiInfo.rbmType = FetiInfo::translation;
          }
          else $$->fetiInfo.rbmType = (FetiInfo::RbmType) $4;

          if($$->fetiInfo.rbmType == FetiInfo::all)
            $$->fetiInfo.nGs = 6;
          else if($$->fetiInfo.rbmType != FetiInfo::None)
            $$->fetiInfo.nGs = 3; }
        | Solver AUGMENT EDGEWS Integer NewLine
        { $$->fetiInfo.numdir = $4; 
          if($$->fetiInfo.augment == FetiInfo::none)
            $$->fetiInfo.augment = FetiInfo::Edges; }
        | Solver AUGMENT EDGEWS WAVETYPE Integer NewLine
        { $$->fetiInfo.waveType = (FetiInfo::WaveType) $4;
          $$->fetiInfo.numdir = $5; 
          if($$->fetiInfo.augment == FetiInfo::none)
            $$->fetiInfo.augment = FetiInfo::Edges; }
        | Solver AUGMENT EDGEWS Integer WAVEMETHOD NewLine
        { $$->fetiInfo.numdir = $4;
          $$->fetiInfo.waveMethod = (FetiInfo::WaveMethod) $5;
          if($$->fetiInfo.augment == FetiInfo::none)
            $$->fetiInfo.augment = FetiInfo::Edges; }
        | Solver AUGMENT EDGEWS WAVETYPE Integer WAVEMETHOD NewLine
        { $$->fetiInfo.waveType = (FetiInfo::WaveType) $4;
          $$->fetiInfo.waveMethod = (FetiInfo::WaveMethod) $6;
          $$->fetiInfo.numdir = $5;
          if($$->fetiInfo.augment == FetiInfo::none)
            $$->fetiInfo.augment = FetiInfo::Edges; }
        | Solver AUGMENT MPCTYPEID NewLine
        { if ($3 == 2)
           $$->fetiInfo.augmentimpl = FetiInfo::Primal;
        }
	| Solver ORTHOTOL Float NewLine
	{ $$->fetiInfo.orthotol = $3; }
        | Solver ORTHOTOL Float Float NewLine
        { $$->fetiInfo.orthotol = $3; 
          $$->fetiInfo.orthotol2 = $4; }
	| Solver GLOBALTOL Float NewLine
	{ $$->fetiInfo.grbm_tol = $3; }
        | Solver GLOBALCRBMTOL Float NewLine
	{ $$->fetiInfo.crbm_tol = $3; }
        | Solver CCTTOL Float NewLine
        { $$->fetiInfo.cct_tol = $3; }
        | Solver REBUILDCCT SWITCH NewLine
        { $$->fetiInfo.rebuildcct = int($3); }
        | Solver UPROJ Integer NewLine
        { $$->fetiInfo.uproj = $3; }
        | Solver printmatlab NewLine
        { $$->fetiInfo.printMatLab = 1; }
        | Solver printmatlab FNAME NewLine
        { $$->printMatLab = 1;
          $$->printMatLabFile = $3; }
        | Solver LOCALSOLVER SolverMethod
        { $$->fetiInfo.local_cntl = $$->fetiInfo.kii_cntl = $3; }
        | Solver COARSESOLVER SolverMethod
        { $$->fetiInfo.coarse_cntl = $3; }
        | Solver AUXCOARSESOLVER SolverMethod
        { $$->fetiInfo.auxcoarse_cntl = $3; }
        | Solver CCTSOLVER SolverMethod
        { $$->fetiInfo.cct_cntl = $3; }
	| Solver VERSION Integer NewLine
	{ if($3 == 1)
            $$->fetiInfo.version = FetiInfo::feti1; 
          else if($3 == 2) {
            $$->fetiInfo.version = FetiInfo::feti2;
            if($$->fetiInfo.nonLocalQ == 1) {
              $$->fetiInfo.nonLocalQ = 0;
              $$->fetiInfo.nQ = 0;
              fprintf(stderr," *** WARNING: Basic projector is used "
                             "with FETI 2\n");
            }
          } else if($3 == 3) {
            $$->fetiInfo.version = FetiInfo::feti3;
            if($$->fetiInfo.nonLocalQ == 1) {
              $$->fetiInfo.nonLocalQ = 0;
              fprintf(stderr," *** WARNING: Basic projector is used "
                             "with FETI 2\n");
            }
          } 
          else {
            $$->fetiInfo.version = FetiInfo::feti1;
	    fprintf(stderr," *** WARNING: Version does not exist,"
                           " using FETI 1\n");
          } }
        | Solver GMRESRESIDUAL NewLine
        { $$->fetiInfo.gmresResidual = true; }
        | Solver GMRESRESIDUAL SWITCH NewLine
        { $$->fetiInfo.gmresResidual = bool($3); }
        | Solver PICKANYCORNER Integer NewLine
        { $$->fetiInfo.pickAnyCorner = $3; }
        | Solver NLPREC NewLine
        { $$->fetiInfo.type = FetiInfo::nonlinear;
          $$->fetiInfo.nlPrecFlg = 1; }
	| Solver NLPREC Integer NewLine
        { $$->fetiInfo.type = FetiInfo::nonlinear;
	  $$->fetiInfo.nlPrecFlg = $3; }
        | Solver NUMCGM Integer NewLine
        { $$->fetiInfo.numcgm = $3; }
        | Solver TOLCGM Float NewLine
        { $$->fetiInfo.tolcgm = $3; }
        | Solver SPACEDIMENSION Integer NewLine
        { $$->fetiInfo.spaceDimension = $3; }
        | Solver KRYLOVTYPE Integer NewLine
        { $$->fetiInfo.krylovtype = $3; }
        | Solver KRYLOVTYPE ISOLVERTYPE NewLine
        { $$->fetiInfo.krylovtype = $3; }
        | Solver INTERFACELUMPED NewLine
        { $$->fetiInfo.lumpedinterface = 1; }
        | Solver SAVEMEMCOARSE NewLine
        { $$->fetiInfo.saveMemCoarse = 1; }
        | Solver OUTERLOOP ITERTYPE NewLine
        { $$->fetiInfo.outerloop = (FetiInfo::OuterloopType) $3; }
        | Solver OUTERLOOP ITERTYPE HERMITIAN NewLine
        { $$->fetiInfo.outerloop = (FetiInfo::OuterloopType) $3;
          $$->fetiInfo.complex_hermitian = true; }
        | Solver MPCTYPE Integer NewLine
        { $$->fetiInfo.mpcflag = $3; }
        | Solver MPCTYPE MPCTYPEID NewLine
        { $$->fetiInfo.mpcflag = $3; }
        | Solver MPCPRECNO Integer NewLine
        { $$->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $3; }
        | Solver MPCPRECNO MPCPRECNOID NewLine
        { $$->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $3; }
        | Solver MPCPRECNO Integer MPCBLK_OVERLAP Integer NewLine
        { $$->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $3;
          $$->fetiInfo.mpcBlkOverlap = $5; }
        | Solver MPCPRECNO MPCPRECNOID MPCBLK_OVERLAP Integer NewLine
        { $$->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $3;
          $$->fetiInfo.mpcBlkOverlap = $5; }
        | Solver MPCPRECNO Integer Integer NewLine
        { $$->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $3; 
          $$->fetiInfo.mpc_block = (FetiInfo::MpcBlock) $4; }
        | Solver MPCPRECNO MPCPRECNOID MPCBLOCKID NewLine
        { $$->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $3;
          $$->fetiInfo.mpc_block = (FetiInfo::MpcBlock) $4; }
        | Solver MPCPRECNO Integer Integer MPCBLK_OVERLAP Integer NewLine
        { $$->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $3;
          $$->fetiInfo.mpc_block = (FetiInfo::MpcBlock) $4;
          $$->fetiInfo.mpcBlkOverlap = $6; }
        | Solver MPCPRECNO MPCPRECNOID MPCBLOCKID MPCBLK_OVERLAP Integer NewLine
        { $$->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $3;
          $$->fetiInfo.mpc_block = (FetiInfo::MpcBlock) $4; 
          $$->fetiInfo.mpcBlkOverlap = $6; }
        | Solver MRHS Integer NewLine
        { if($3 < 1) $$->fetiInfo.useMRHS = false; }
        | Solver PROPORTIONING Float NewLine
        { $$->fetiInfo.gamma = $3; }
        | Solver LINESEARCH Integer Float NewLine
        { $$->fetiInfo.linesearch_maxit = $3;
          $$->fetiInfo.linesearch_tau = $4; }
        | Solver BMPC SWITCH NewLine
        { $$->fetiInfo.bmpc = bool($3); }
        | Solver CMPC SWITCH NewLine
        { $$->fetiInfo.cmpc = bool($3); }
        | Solver CNORM SWITCH NewLine
        { $$->fetiInfo.c_normalize = bool($3); }
        ;
OldHelmInfo:
        FACOUSTICS NewLine FAcousticData
	;
FAcousticData:
        Float NewLine
        {
          geoSource->setOmega($1);
          StructProp sp; 
          sp.kappaHelm = $1;
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::Helmholtz);
        }
	;
Constraints:
        CONSTRAINTS ConstraintOptionsData NewLine
        { if(!$2.lagrangeMult && $2.penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = $2.lagrangeMult;
          domain->solInfo().penalty = $2.penalty;
          domain->solInfo().constraint_hess = $2.constraint_hess; 
          domain->solInfo().constraint_hess_eps = $2.constraint_hess_eps; }
        | CONSTRAINTS NewLine ConstraintOptionsData NewLine
        { if(!$3.lagrangeMult && $3.penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = $3.lagrangeMult;
          domain->solInfo().penalty = $3.penalty;
          domain->solInfo().constraint_hess = $3.constraint_hess;
          domain->solInfo().constraint_hess_eps = $3.constraint_hess_eps; }
        ;
ConstraintOptionsData:
        DIRECT
        { // Direct elimination of slave dofs
          $$.lagrangeMult = false;
          $$.penalty = 0.0;
          $$.constraint_hess = 0;
          $$.constraint_hess_eps = 0.0;
        }
        | DIRECT Float
        { $$.lagrangeMult = false; 
          $$.penalty = 0.0;
          $$.constraint_hess = 0;
          $$.constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = $2; }
        | DIRECT Float Float
        { $$.lagrangeMult = false; 
          $$.penalty = 0.0;
          $$.constraint_hess = 0;
          $$.constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = $2;
          domain->solInfo().coefFilterTol = $3; }
        | DIRECT Float Float Float
        { $$.lagrangeMult = false; 
          $$.penalty = 0.0;
          $$.constraint_hess = 0;
          $$.constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = $2; 
          domain->solInfo().coefFilterTol = $3;
          domain->solInfo().rhsZeroTol = $4; }
        | DIRECT Float Float Float Float
        { $$.lagrangeMult = false; 
          $$.penalty = 0.0;
          $$.constraint_hess = 0;
          $$.constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = $2;
          domain->solInfo().coefFilterTol = $3; 
          domain->solInfo().rhsZeroTol = $4;
          domain->solInfo().inconsistentTol = $5; }
        | MULTIPLIERS
        { // Treatment of constraints through Lagrange multipliers method
          $$.lagrangeMult = true; 
          $$.penalty = 0.0;
          $$.constraint_hess = 1;
          $$.constraint_hess_eps = 0.0; }
        | PENALTY Float
        { // Treatment of constraints through penalty method
          $$.lagrangeMult = false;
          $$.penalty = $2;
          $$.constraint_hess = 1;
          $$.constraint_hess_eps = 0.0; }
        | MULTIPLIERS PENALTY Float
        { // Treatment of constraints through augmented Lagrangian method
          $$.lagrangeMult = true;
          $$.penalty = $3;
          $$.constraint_hess = 1;
          $$.constraint_hess_eps = 0.0; }
        | AUGMENTED Float
        { // Alternative input syntax for treatment of constraints through augmented Lagrangian method
          $$.lagrangeMult = true;
          $$.penalty = $2;
          $$.constraint_hess = 1;
          $$.constraint_hess_eps = 0.0; }
        | ConstraintOptionsData HESSIAN Integer
        { $$.constraint_hess = $3;
          $$.constraint_hess_eps = 0; }
        | ConstraintOptionsData HESSIAN Integer Float
        { $$.constraint_hess = $3;
          $$.constraint_hess_eps = $4; }
        ;
HelmInfo:
        HELMHOLTZ NewLine
        { // hack??
	  domain->solInfo().acoustic = true;
          if(domain->solInfo().probType != SolverInfo::HelmholtzDirSweep) domain->solInfo().setProbType(SolverInfo::Helmholtz);
        }
/* shouldn't be here
        | FETIH NewLine
        { domain->solInfo().solvercntl->type = (2); 
          domain->solInfo().solvercntl->fetiInfo.scaling = FetiInfo::tscaling; 
          domain->solInfo().solvercntl->fetiInfo.corners = FetiInfo::allCorners3;
          domain->solInfo().solvercntl->fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().solvercntl->fetiInfo.dph_flag = true;
          domain->solInfo().solvercntl->fetiInfo.augment = FetiInfo::Edges;
          domain->solInfo().solvercntl->fetiInfo.rbmType = FetiInfo::None;
          domain->solInfo().solvercntl->fetiInfo.nGs = 0;
        }
*/
        | BGTL Integer NewLine
        {
          domain->sommerfeldType = $2;
          domain->curvatureFlag = 0;
        }
        | BGTL Integer Float NewLine
        {
          domain->sommerfeldType = $2;
          domain->curvatureConst1 = $3;
          domain->curvatureFlag = 1;
        }
        | BGTL Integer Float Float NewLine
        {
          domain->sommerfeldType = $2;
          domain->curvatureConst1 = $3;
          domain->curvatureConst2 = $4;
          domain->curvatureFlag = 2;
        }
        | POINTSOURCE Integer NewLine IncidenceList
        {
          domain->pointSourceFlag = 1;
          domain->implicitFlag = 1;
        }
        | PLANEWAVE Integer NewLine IncidenceList
        {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        }
        | INCIDENCE Integer NewLine IncidenceList
        {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        }
        | HelmInfo DampInfo
        | HelmInfo ReconsInfo 
        | HelmInfo PadePivotInfo
        ;
IncidenceList:
        IncidenceVector
        | IncidenceList IncidenceVector
        ;
IncidenceVector:
        Float Float Float NewLine
        { domain->setWaveDirections(0, $1, $2, $3); }
        ;
KirchhoffLocations:
        KIRLOC NewLine Float Float Float NewLine{
          domain->setKirchhoffLocations($3, $4, $5);
        } 
        | KirchhoffLocations Float Float Float NewLine {
          domain->setKirchhoffLocations($2, $3, $4);
        }
        ;
FFPDirList:
        FFPDirVector
        | FFPDirList FFPDirVector
        ;
FFPDirVector:
        Float Float Float NewLine
        { domain->setFFPDirections($1, $2, $3); }
        ;
HelmMFInfo:
        HELMMF NewLine FAcousticDataMF
        ;
FAcousticDataMF:
        Float NewLine
        {
          /*domain->omega = $1;*/ geoSource->setOmega($1);
          StructProp sp;
          sp.kappaHelm = $1;
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzMF);
        }
        ;
HelmSOInfo:
        HELMSO NewLine FAcousticDataSO
        ;
FAcousticDataSO:
        Float NewLine
        {
          /*domain->omega = $1;*/ geoSource->setOmega($1);
          StructProp sp;
          sp.kappaHelm = $1;
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzSO);
        }
        ;
DEMInfo:
        DEM NewLine
        {
          domain->solInfo().setProbType(SolverInfo::DisEnrM);
        }
        ;
AlProc:
        ALPROC
        { $$ = $1; }
        | Integer
        { $$ = $1; }
NLInfo:
        NL NewLine
        { if(domain->solInfo().probType == SolverInfo::Static || domain->solInfo().probType == SolverInfo::None)
            domain->solInfo().probType = SolverInfo::NonLinStatic;
          else if(domain->solInfo().probType == SolverInfo::Dynamic)
            domain->solInfo().probType = SolverInfo::NonLinDynam;
          else if(domain->solInfo().probType == SolverInfo::TempDynamic) {
            domain->solInfo().order = 1;
            domain->solInfo().probType = SolverInfo::NonLinDynam;
          }
          domain->solInfo().solvercntl->fetiInfo.type = FetiInfo::nonlinear;
          domain->solInfo().getNLInfo().setDefaults(); /* just in case PIECEWISE is used under statics */
          domain->solInfo().nlFlag = 1; /* can be used for decomposition when a different treatment is required, e.g. contact */ }
        | NLInfo NEWTON NewLine
        {}
        | NLInfo ARCLENGTH NewLine
        { if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::ArcLength; }
        | NLInfo NLMAT NewLine
        { if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::MatNonLinStatic;
          else if(domain->solInfo().probType == SolverInfo::NonLinDynam)
            domain->solInfo().probType = SolverInfo::MatNonLinDynam; }
        | NLInfo LINEARELASTIC NewLine
        { domain->solInfo().getNLInfo().linearelastic = 1; }
        | NLInfo LINEARELASTIC Integer NewLine
        { domain->solInfo().getNLInfo().linearelastic = $3; }
        | NLInfo MAXITR Integer NewLine
        { domain->solInfo().getNLInfo().maxiter = $3; }
        | NLInfo NLTOL Float NewLine
        { domain->solInfo().getNLInfo().tolRes = $3; }
        | NLInfo NLTOL Float Float NewLine
        { domain->solInfo().getNLInfo().tolRes = $3;
          domain->solInfo().getNLInfo().tolInc = $4; }
        | NLInfo NLTOL Float Float Float Float NewLine
        { domain->solInfo().getNLInfo().tolRes = $3;
          domain->solInfo().getNLInfo().tolInc = $4;
          domain->solInfo().getNLInfo().absTolRes = $5;
          domain->solInfo().getNLInfo().absTolInc = $6; }
        | NLInfo DLAMBDA Float Float NewLine
        { domain->solInfo().getNLInfo().dlambda = $3;
          domain->solInfo().getNLInfo().maxLambda = $4; }
        | NLInfo DLAMBDA Float Float Integer Integer NewLine
        { domain->solInfo().getNLInfo().dlambda = $3; 
          domain->solInfo().getNLInfo().maxLambda = $4;
          domain->solInfo().getNLInfo().extMin = $5;
          domain->solInfo().getNLInfo().extMax = $6; }
        | NLInfo FITALG Integer NewLine
        { domain->solInfo().getNLInfo().fitAlgShell = $3;
          domain->solInfo().getNLInfo().fitAlgBeam  = $3; }
        | NLInfo FITALG Integer Integer NewLine
        { domain->solInfo().getNLInfo().fitAlgShell = $3;
          domain->solInfo().getNLInfo().fitAlgBeam  = $4; }
        | NLInfo UNSYMMETRIC NewLine
        { domain->solInfo().getNLInfo().unsymmetric = true; }
        | NLInfo LFACTOR Float NewLine
        { domain->solInfo().getNLInfo().lfactor = $3; }
        | NLInfo LINESEARCH AlProc NewLine
        { domain->solInfo().getNLInfo().linesearch.type = $3; }
        | NLInfo LINESEARCH AlProc Integer NewLine
        { domain->solInfo().getNLInfo().linesearch.type = $3; 
          domain->solInfo().getNLInfo().linesearch.maxit = $4; }
        | NLInfo LINESEARCH AlProc Integer Float NewLine
        { domain->solInfo().getNLInfo().linesearch.type = $3;
          domain->solInfo().getNLInfo().linesearch.maxit = $4;
          // note: currently we use either c1 or c2, but never both
          domain->solInfo().getNLInfo().linesearch.c1 = $5;
          domain->solInfo().getNLInfo().linesearch.c2 = $5; }
        | NLInfo LINESEARCH AlProc Integer Float Float NewLine
        { domain->solInfo().getNLInfo().linesearch.type = $3;
          domain->solInfo().getNLInfo().linesearch.maxit = $4;
          domain->solInfo().getNLInfo().linesearch.c1 = $5; 
          domain->solInfo().getNLInfo().linesearch.c2 = $5;
          domain->solInfo().getNLInfo().linesearch.tau = $6; }
        | NLInfo LINESEARCH AlProc Integer Float Float Integer NewLine
        { domain->solInfo().getNLInfo().linesearch.type = $3;
          domain->solInfo().getNLInfo().linesearch.maxit = $4;
          domain->solInfo().getNLInfo().linesearch.c1 = $5; 
          domain->solInfo().getNLInfo().linesearch.c2 = $5;
          domain->solInfo().getNLInfo().linesearch.tau = $6;
          domain->solInfo().getNLInfo().linesearch.verbose = bool($7); }
        | NLInfo FAILSAFE NewLine
        { domain->solInfo().getNLInfo().failsafe = true; }
        | NLInfo FAILSAFE Float NewLine
        { domain->solInfo().getNLInfo().failsafe = true;
          domain->solInfo().getNLInfo().failsafe_tol = $3; }
        | NLInfo PENALTY Integer Float Float NewLine
        { domain->solInfo().num_penalty_its = $3; 
          domain->solInfo().penalty_tol = $4;
          domain->solInfo().penalty_beta = $5; }
        | NLInfo PENALTY Integer Float Float Integer Integer NewLine
        { domain->solInfo().num_penalty_its = $3;
          domain->solInfo().penalty_tol = $4;
          domain->solInfo().penalty_beta = $5;
          domain->solInfo().reinit_lm = bool($6);
          domain->solInfo().lm_update_flag = $7; }
        | NLInfo NUMROMCPUS Integer NewLine
        { domain->solInfo().numberOfRomCPUs = $3; }
        | NLInfo NewtonInfo
        ;
NewtonInfo:
        REBUILD Integer NewLine
        {
          domain->solInfo().setNewton($2); 
          domain->solInfo().solvercntl->fetiInfo.type  = FetiInfo::nonlinear; 
        }
        | REBUILD Integer Integer NewLine
        {
          domain->solInfo().setNewton($2);
          domain->solInfo().getNLInfo().stepUpdateK = $3;
          domain->solInfo().solvercntl->fetiInfo.type  = FetiInfo::nonlinear;
        }
        | REBUILD Integer Integer Integer NewLine
        {
          domain->solInfo().setNewton($2);
          domain->solInfo().getNLInfo().stepUpdateK = $3;
          domain->solInfo().piecewise_contact = bool($4);
          domain->solInfo().solvercntl->fetiInfo.type  = FetiInfo::nonlinear;
        }
	;
OrthoInfo:
	REORTHO NewLine
        { domain->solInfo().setReOrtho(); }
	;
Control:
	CONTROL NewLine
        | CONTROL NewLine FNAME NewLine Integer NewLine FNAME NewLine FNAME NewLine
         { geoSource->setControl($3,$7,$9); domain->solInfo().soltyp = $5; }
        | CONTROL NewLine FNAME NewLine Integer NewLine FNAME NewLine FNAME NewLine FNAME NewLine
         { geoSource->setControl($3,$7,$9,$11); domain->solInfo().soltyp = $5; }
/*
        | Control FNAME NewLine Integer NewLine FNAME NewLine FNAME NewLine
         { geoSource->setControl($2,$6,$8); }
        | Control FNAME NewLine Integer NewLine FNAME NewLine FNAME NewLine FNAME NewLine
         { geoSource->setControl($2,$6,$8,$10); }
*/
	;
NodalContact: 
        NODALCONTACT NewLine
        | NODALCONTACT MODE Integer NewLine
        { domain->solInfo().contact_mode = $3; }

        // glbNode01 glbNode02 nx ny nz
        | NodalContact Integer Integer Float Float Float NewLine
        { domain->addNodalCTC($2-1, $3-1, $4, $5, $6); }
        // glbNode01 glbNode02 nx ny nz GAP normalGap
        | NodalContact Integer Integer Float Float Float GAP Float NewLine
        { domain->addNodalCTC($2-1, $3-1, $4, $5, $6, $8);}
        // glbNode01 glbNode02 nx ny nz MODE mode
        | NodalContact Integer Integer Float Float Float MODE Integer NewLine
        { domain->addNodalCTC($2-1, $3-1, $4, $5, $6, 0.0, $8);}
        | NodalContact Integer Integer Float Float Float ConstraintOptionsData NewLine
        { domain->addNodalCTC($2-1, $3-1, $4, $5, $6, 0.0, -1, $7.lagrangeMult, $7.penalty);}
        // glbNode01 glbNode02 nx ny nz GAP normalGap MODE mode
        | NodalContact Integer Integer Float Float Float GAP Float MODE Integer NewLine
        { domain->addNodalCTC($2-1, $3-1, $4, $5, $6, $8, $10);}
        | NodalContact Integer Integer Float Float Float GAP Float MODE Integer ConstraintOptionsData NewLine
        { domain->addNodalCTC($2-1, $3-1, $4, $5, $6, $8, $10, $11.lagrangeMult, $11.penalty);}
	;
MatSpec:
	MATSPEC NewLine
	| MatSpec Integer BILINEARPLASTIC Float Float Float Float Float NewLine
	 { 
           geoSource->addMaterial($2-1, 
             new BilinPlasKinHardMat($4, $5, $6, $7, $8) );
         }
        | MatSpec Integer BILINEARPLASTIC Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new BilinPlasKinHardMat($4, $5, $6, $7, $8, $9) );
         }
        | MatSpec Integer BILINEARPLASTIC Float Float Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new BilinPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11) );
         }
        | MatSpec Integer BILINEARPLASTIC Float Float Float Float Float Float Float Float Float NewLine
         {
           if($12 > 0 && $12 < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial($2-1, new BilinPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11, $12) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial($2-1, new BilinPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11) );
           }
         }
        | MatSpec Integer BILINEARPLASTIC Float Float Float Float Float Float Float Float Float Float NewLine
         {
           if($12 > 0 && $12 < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial($2-1, new BilinPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11, $12, $13) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial($2-1, new BilinPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11,
                                    std::numeric_limits<double>::infinity(), $13) );
           }
         }
        | MatSpec Integer BILINEARPLASTIC Float Float Float Float Float Float Float Float Float Float Integer NewLine
         {
           if($12 > 0 && $12 < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial($2-1, new BilinPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial($2-1, new BilinPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11,
                                    std::numeric_limits<double>::infinity(), $13, $14) );
           }
         }
        | MatSpec Integer FINITESTRAINPLASTIC Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new FiniteStrainPlasKinHardMat($4, $5, $6, $7, $8) );
         }
        | MatSpec Integer FINITESTRAINPLASTIC Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new FiniteStrainPlasKinHardMat($4, $5, $6, $7, $8, $9) );
         }
        | MatSpec Integer FINITESTRAINPLASTIC Float Float Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new FiniteStrainPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11) );
         }
        | MatSpec Integer FINITESTRAINPLASTIC Float Float Float Float Float Float Float Float Float NewLine
         {
           if($12 > 0 && $12 < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial($2-1, new FiniteStrainPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11, $12) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial($2-1, new FiniteStrainPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11) );
           }
         }
        | MatSpec Integer FINITESTRAINPLASTIC Float Float Float Float Float Float Float Float Float Float NewLine
         {
           if($12 > 0 && $12 < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial($2-1, new FiniteStrainPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11, $12, $13) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial($2-1, new FiniteStrainPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11, 
                                    std::numeric_limits<double>::infinity(), $13) );
           }
         }
        | MatSpec Integer FINITESTRAINPLASTIC Float Float Float Float Float Float Float Float Float Float Integer NewLine
         {
           if($12 > 0 && $12 < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial($2-1, new FiniteStrainPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11, $12, $13, 414) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial($2-1, new FiniteStrainPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11,
                                    std::numeric_limits<double>::infinity(), $13, $14) );
           }
         }
        | MatSpec Integer CRUSHABLEFOAM Float Float Float Float Float Float Float Float NewLine
          {
            geoSource->addMaterial($2-1, new CrushableFoam($4, $5, $6, $7, $8, $9, $10, $11) );
          }
        | MatSpec Integer CRUSHABLEFOAM Float Float Float Float Float Float Float Float Float NewLine
          {
            geoSource->addMaterial($2-1, new CrushableFoam($4, $5, $6, $7, $8, $9, $10, $11, $12) );
          }
        | MatSpec Integer CRUSHABLEFOAM Float Float Float Float Float Float Float Float Float Float NewLine
          {
            geoSource->addMaterial($2-1, new CrushableFoam($4, $5, $6, $7, $8, $9, $10, $11, $12, $13) );
          }
        | MatSpec Integer CRUSHABLEFOAM Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            if($14 > 0 && $14 < std::numeric_limits<double>::infinity()) {
              geoSource->addMaterial($2-1, new CrushableFoam($4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14) );
              domain->solInfo().elementDeletion = true;
            }
            else {
              geoSource->addMaterial($2-1, new CrushableFoam($4, $5, $6, $7, $8, $9, $10, $11, $12, 
                                     $13, std::numeric_limits<double>::infinity()) );
            }
          }
        | MatSpec Integer LOGSTRAINPLASTIC Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new LogStrainPlasKinHardMat($4, $5, $6, $7, $8) );
         }
        | MatSpec Integer LOGSTRAINPLASTIC Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new LogStrainPlasKinHardMat($4, $5, $6, $7, $8, $9) );
         }
        | MatSpec Integer LOGSTRAINPLASTIC Float Float Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new LogStrainPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11) );
         }
        | MatSpec Integer LOGSTRAINPLASTIC Float Float Float Float Float Float Float Float Float NewLine
         {
           if($12 > 0 && $12 < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial($2-1, new LogStrainPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11, $12) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial($2-1, new LogStrainPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11) );
           }
         }
        | MatSpec Integer LOGSTRAINPLASTIC Float Float Float Float Float Float Float Float Float Float NewLine
         {
           if($12 > 0 && $12 < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial($2-1, new LogStrainPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11, $12, $13) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial($2-1, new LogStrainPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11, 
                                    std::numeric_limits<double>::infinity(), $13) );
           }
         }
        | MatSpec Integer LOGSTRAINPLASTIC Float Float Float Float Float Float Float Float Float Float Integer NewLine
         {
           if($12 > 0 && $12 < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial($2-1, new LogStrainPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial($2-1, new LogStrainPlasKinHardMat($4, $5, $6, $7, $8, $9, $10, $11,
                                    std::numeric_limits<double>::infinity(), $13, $14) );
           }
         }
        | MatSpec Integer LINEARELASTIC Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ElaLinIsoMat($4, $5, $6, $7, $8));
         }
	| MatSpec Integer LINEARELASTIC Float Float Float NewLine
	 { 
           geoSource->addMaterial($2-1, 
             new ElaLinIsoMat($4, $5, $6));
	 }
        | MatSpec Integer LINEARELASTIC Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ElaLinIsoMat($4));
         }
        | MatSpec Integer LINEARELASTIC Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new BrittleFractureTB<ElaLinIsoMat>($4, $5, $6, $7, $8, $10, $11, $12));
           domain->solInfo().elementDeletion = true;
         }
        | MatSpec Integer LINEARELASTIC Float Float Float TULERBUTCHER Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new BrittleFractureTB<ElaLinIsoMat>($4, $5, $6, $8, $9, $10));
           domain->solInfo().elementDeletion = true;
         }
        | MatSpec Integer LINEARELASTIC Float TULERBUTCHER Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new BrittleFractureTB<ElaLinIsoMat>($4, $6, $7, $8));
           domain->solInfo().elementDeletion = true;
         }
        | MatSpec Integer STVENANTKIRCHHOFF Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new StVenantKirchhoffMat($4, $5, $6, $7, $8));
         }
        | MatSpec Integer STVENANTKIRCHHOFF Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new StVenantKirchhoffMat($4, $5, $6));
         }
        | MatSpec Integer STVENANTKIRCHHOFF Float NewLine
         {
           geoSource->addMaterial($2-1,
             new StVenantKirchhoffMat($4));
         }
        | MatSpec Integer STVENANTKIRCHHOFF Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new BrittleFractureTB<StVenantKirchhoffMat>($4, $5, $6, $7, $8, $10, $11, $12));
           domain->solInfo().elementDeletion = true;
         }
        | MatSpec Integer STVENANTKIRCHHOFF Float Float Float TULERBUTCHER Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new BrittleFractureTB<StVenantKirchhoffMat>($4, $5, $6, $8, $9, $10));
           domain->solInfo().elementDeletion = true;
         }
        | MatSpec Integer STVENANTKIRCHHOFF Float TULERBUTCHER Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new BrittleFractureTB<StVenantKirchhoffMat>($4, $6, $7, $8));
           domain->solInfo().elementDeletion = true;
         }
        | MatSpec Integer HENCKY Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new HenckyMat($4, $5, $6, $7, $8));
         }
        | MatSpec Integer HENCKY Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new HenckyMat($4, $5, $6));
         }
        | MatSpec Integer HENCKY Float NewLine
         {
           geoSource->addMaterial($2-1,
             new HenckyMat($4));
         }
        | MatSpec Integer HENCKY Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new BrittleFractureTB<HenckyMat>($4, $5, $6, $7, $8, $10, $11, $12));
           domain->solInfo().elementDeletion = true;
         }
        | MatSpec Integer HENCKY Float Float Float TULERBUTCHER Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new BrittleFractureTB<HenckyMat>($4, $5, $6, $8, $9, $10));
           domain->solInfo().elementDeletion = true;
         }
        | MatSpec Integer HENCKY Float TULERBUTCHER Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new BrittleFractureTB<HenckyMat>($4, $6, $7, $8));
           domain->solInfo().elementDeletion = true;
         }
        | MatSpec Integer LINPLSTRESS Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ElaLinIsoMat2D($4, $5, $6, $7, 0, 0));
         }
        | MatSpec Integer LINPLSTRESS Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ElaLinIsoMat2D($4, $5, $6, $7, $8, $9));
         }
        | MatSpec Integer SVKPLSTRESS Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new StVenantKirchhoffMat2D($4, $5, $6, $7, 0, 0));
         }
        | MatSpec Integer SVKPLSTRESS Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new StVenantKirchhoffMat2D($4, $5, $6, $7, $8, $9));
         }
        | MatSpec Integer FABRICMAP Float Integer Integer Integer Float NewLine
         {
           geoSource->addMaterial($2-1,
             new FabricMap($4, $5, $6, $7, $8, 0, 0, FabricMap::GREEN_LAGRANGE));
         }
        | MatSpec Integer FABRICMAP Float Integer Integer Integer Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new FabricMap($4, $5, $6, $7, $8, $9, $10, FabricMap::GREEN_LAGRANGE));
         }
        | MatSpec Integer SHELLFABRICMAP Float Integer Integer Integer Float NewLine
         {
           geoSource->addMaterial($2-1,
             new FabricMap($4, $5, $6, $7, $8, 0, 0, FabricMap::INFINTESIMAL));
         }
        | MatSpec Integer SHELLFABRICMAP Float Integer Integer Integer Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new FabricMap($4, $5, $6, $7, $8, $9, $10, FabricMap::INFINTESIMAL));
         }
        | MatSpec Integer FABRICMAT Float Integer Integer Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new FabricMat($4, $5, $6, $7, $8, $9, $10, 0, 0, FabricMat::GREEN_LAGRANGE));
         }
        | MatSpec Integer FABRICMAT Float Integer Integer Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new FabricMat($4, $5, $6, $7, $8, $9, $10, $11, $12, FabricMat::GREEN_LAGRANGE));
         }
        | MatSpec Integer SHELLFABRICMAT Float Integer Integer Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new FabricMat($4, $5, $6, $7, $8, $9, $10, 0, 0, FabricMat::INFINTESIMAL));
         }
        | MatSpec Integer SHELLFABRICMAT Float Integer Integer Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new FabricMat($4, $5, $6, $7, $8, $9, $10, $11, $12, FabricMat::INFINTESIMAL));
         }
        | MatSpec Integer PLANESTRESSLINEAR Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<ElaLinIsoMat>($4, $5, $6, $7));
         }
        | MatSpec Integer PLANESTRESSLINEAR Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<ElaLinIsoMat>($4, $5, $6, $7, $8, $9));
         }
        | MatSpec Integer PLANESTRESSLINEAR Float Float Float Float TULERBUTCHER Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<BrittleFractureTB<ElaLinIsoMat> >($4, $5, $6, $9, $10, $11, $7));
           domain->solInfo().elementDeletion = true;
         }
        | MatSpec Integer PLANESTRESSLINEAR Float Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<BrittleFractureTB<ElaLinIsoMat> >($4, $5, $6, $7, $8, $11, $12, $13, $9));
           domain->solInfo().elementDeletion = true;
         }
        | MatSpec Integer PLANESTRESSSTVENANTKIRCHHOFF Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<StVenantKirchhoffMat>($4, $5, $6, $7));
         }
        | MatSpec Integer PLANESTRESSSTVENANTKIRCHHOFF Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<StVenantKirchhoffMat>($4, $5, $6, $7, $8, $9));
         }
        | MatSpec Integer PLANESTRESSSTVENANTKIRCHHOFF Float Float Float Float TULERBUTCHER Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<BrittleFractureTB<StVenantKirchhoffMat> >($4, $5, $6, $9, $10, $11, $7));
           domain->solInfo().elementDeletion = true;
         }
        | MatSpec Integer PLANESTRESSSTVENANTKIRCHHOFF Float Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<BrittleFractureTB<StVenantKirchhoffMat> >($4, $5, $6, $7, $8, $11, $12, $13, $9));
           domain->solInfo().elementDeletion = true;
         }
        | MatSpec Integer PLANESTRESSNEOHOOKEAN Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<NeoHookeanMat>($4, $5, $6, $7));
         }
        | MatSpec Integer PLANESTRESSNEOHOOKEAN Float Float Float Float TULERBUTCHER Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<BrittleFractureTB<NeoHookeanMat> >($4, $5, $6, $9, $10, $11, $7));
           domain->solInfo().elementDeletion = true;
         }
        | MatSpec Integer PLANESTRESSMOONEYRIVLIN Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<MooneyRivlinMat>($4, $5, $6, $7, $8));
         }
        | MatSpec Integer PLANESTRESSMOONEYRIVLIN Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<BrittleFractureTB<MooneyRivlinMat> >($4, $5, $6, $7, $10, $11, $12, $8));
           domain->solInfo().elementDeletion = true;
         }
        | MatSpec Integer PLANESTRESSBILINEARPLASTIC Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<BilinPlasKinHardMat>($4, $5, $6, $7, $8, $9) );
         }
        | MatSpec Integer PLANESTRESSBILINEARPLASTIC Float Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<BilinPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10) );
         }
        | MatSpec Integer PLANESTRESSBILINEARPLASTIC Float Float Float Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<BilinPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10, $11, $12) );
         }
        | MatSpec Integer PLANESTRESSBILINEARPLASTIC Float Float Float Float Float Float Float Float Float Float NewLine
         {
           if($13 > 0 && $13 < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial($2-1, new PlaneStressMat<BilinPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10, $11, $13, $12) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial($2-1, new PlaneStressMat<BilinPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10, $11, $12) );
           }
         }
        | MatSpec Integer PLANESTRESSBILINEARPLASTIC Float Float Float Float Float Float Float Float Float Float Float NewLine
         {
           if($13 > 0 && $13 < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial($2-1, new PlaneStressMat<BilinPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10, $11, $13, $14, $12) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial($2-1, new PlaneStressMat<BilinPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10, $11,
                                    std::numeric_limits<double>::infinity(), $14, $12) );
           }
         }
        | MatSpec Integer PLANESTRESSBILINEARPLASTIC Float Float Float Float Float Float Float Float Float Float Float Integer NewLine
         {
           if($13 > 0 && $13 < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial($2-1, new PlaneStressMat<BilinPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10, $11, $13, $14, $15, $12) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial($2-1, new PlaneStressMat<BilinPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10, $11,
                                    std::numeric_limits<double>::infinity(), $14, $15, $12) );
           }
         }
        | MatSpec Integer PLANESTRESSFINITESTRAINPLASTIC Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<FiniteStrainPlasKinHardMat>($4, $5, $6, $7, $8, $9) );
         }
        | MatSpec Integer PLANESTRESSFINITESTRAINPLASTIC Float Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<FiniteStrainPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10) );
         }
        | MatSpec Integer PLANESTRESSFINITESTRAINPLASTIC Float Float Float Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new PlaneStressMat<FiniteStrainPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10, $11, $12) );
         }
        | MatSpec Integer PLANESTRESSFINITESTRAINPLASTIC Float Float Float Float Float Float Float Float Float Float NewLine
         {
           if($13 > 0 && $13 < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial($2-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10, $11, $13, $12) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial($2-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10, $11, $12) );
           }
         }
        | MatSpec Integer PLANESTRESSFINITESTRAINPLASTIC Float Float Float Float Float Float Float Float Float Float Float NewLine
         {
           if($13 > 0 && $13 < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial($2-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10, $11, $13, $14, $12) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial($2-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10, $11,
                                    std::numeric_limits<double>::infinity(), $14, $12) );
           }
         }
        | MatSpec Integer PLANESTRESSFINITESTRAINPLASTIC Float Float Float Float Float Float Float Float Float Float Float Integer NewLine
         {
           if($13 > 0 && $13 < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial($2-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10, $11, $13, $14, $15, $12) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial($2-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>($4, $5, $6, $7, $8, $9, $10, $11,
                                    std::numeric_limits<double>::infinity(), $14, $15, $12) );
           }
         }
        | MatSpec Integer ISOTROPICLINEARELASTIC Float Float Float NewLine
          {
            double params[3] = { $4, $5, $6 };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<IsotropicLinearElastic>(params));
          }
        | MatSpec Integer VISCOLINEARELASTIC Float Float Float Float Float Float Float Float Float NewLine
          {
            double gthree = $11;
            double gtwo   = $9;
            double gone   = $7;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PronyViscoElastic<ElaLinIsoMat>($4, $5, $6, ginf, $7, $8, $9, $10, $11, $12));
          }
        | MatSpec Integer VISCOLINEARELASTIC Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double gthree = $13;
            double gtwo   = $11;
            double gone   = $9;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PronyViscoElastic<ElaLinIsoMat>($4, $5, $6, $7, $8, ginf, $9, $10, $11, $12, $13, $14));
          }
        | MatSpec Integer VISCOLINEARELASTIC Float Float Float Float Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
          {
            double gthree = $11;
            double gtwo   = $9;
            double gone   = $7;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> >($4, $5, $6, ginf, $7, $8, $9, $10, $11, $12, $14, $15, $16));
            domain->solInfo().elementDeletion = true;
          }
        | MatSpec Integer VISCOLINEARELASTIC Float Float Float Float Float Float Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
          {
            double gthree = $13;
            double gtwo   = $11;
            double gone   = $9;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> >($4, $5, $6, $7, $8, ginf, $9, $10, $11, $12, $13, $14, $16, $17, $18));
            domain->solInfo().elementDeletion = true;
          }
        | MatSpec Integer VISCOLINPLSTRESS Float Float Float Float Float Float Float Float Float Float NewLine
         {
           double gthree = $12;
           double gtwo   = $10;
           double gone   = $8;
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial($2-1,
             new PronyViscoElastic<ElaLinIsoMat2D>($4, $5, $6, $7, 0, 0, ginf, $8, $9, $10, $11, $12, $13));
         }
        | MatSpec Integer VISCOLINPLSTRESS Float Float Float Float Float Float Float Float Float Float Float Float NewLine
         {
           double gthree = $14;
           double gtwo   = $12;
           double gone   = $10;
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial($2-1,
             new PronyViscoElastic<ElaLinIsoMat2D>($4, $5, $6, $7, $8, $9, ginf, $10, $11, $12, $13, $14, $15));
         }
        | MatSpec Integer VISCOSVKPLSTRESS Float Float Float Float Float Float Float Float Float Float NewLine
         {
           double gthree = $12;
           double gtwo   = $10;
           double gone   = $8;
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial($2-1,
             new PronyViscoElastic<StVenantKirchhoffMat2D>($4, $5, $6, $7, 0, 0, ginf, $8, $9, $10, $11, $12, $13));
         }
        | MatSpec Integer VISCOSVKPLSTRESS Float Float Float Float Float Float Float Float Float Float Float Float NewLine
         {
           double gthree = $14;
           double gtwo   = $12;
           double gone   = $10;
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial($2-1,
             new PronyViscoElastic<StVenantKirchhoffMat2D>($4, $5, $6, $7, $8, $9, ginf, $10, $11, $12, $13, $14, $15));
         }
        | MatSpec Integer VISCOFABRICMAP Float Integer Integer Integer Float Float Float Float Float Float Float NewLine
         {
           double gthree = $13;
           double gtwo   = $11;
           double gone   = $9;
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial($2-1,
             new PronyViscoElastic<FabricMap>($4, $5, $6, $7, $8, 0, 0, FabricMap::GREEN_LAGRANGE, ginf, $9, $10, $11, $12, $13, $14));
         }
        | MatSpec Integer VISCOFABRICMAP Float Integer Integer Integer Float Float Float Float Float Float Float Float Float NewLine
         {
           double gthree = $15;
           double gtwo   = $13;
           double gone   = $11;
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial($2-1,
             new PronyViscoElastic<FabricMap>($4, $5, $6, $7, $8, $9, $10, FabricMap::GREEN_LAGRANGE, ginf, $11, $12, $13, $14, $15, $16));
         }
        | MatSpec Integer SHELLVISCOFABRICMAP Float Integer Integer Integer Float Float Float Float Float Float Float NewLine
         {
           double gthree = $13;
           double gtwo   = $11;
           double gone   = $9;
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial($2-1,
             new PronyViscoElastic<FabricMap>($4, $5, $6, $7, $8, 0, 0, FabricMap::INFINTESIMAL, ginf, $9, $10, $11, $12, $13, $14));
         }
        | MatSpec Integer SHELLVISCOFABRICMAP Float Integer Integer Integer Float Float Float Float Float Float Float Float Float NewLine
         {
           double gthree = $15;
           double gtwo   = $13;
           double gone   = $11;
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial($2-1,
             new PronyViscoElastic<FabricMap>($4, $5, $6, $7, $8, $9, $10, FabricMap::INFINTESIMAL, ginf, $11, $12, $13, $14, $15, $16));
         }
        | MatSpec Integer VISCOFABRICMAT Float Integer Integer Float Float Float Float Float Float Float Float Float Float NewLine
         {
           double gthree = $13;
           double gtwo   = $11;
           double gone   = $9;
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial($2-1,
             new PronyViscoElastic<FabricMat>($4, $5, $6, $7, $8, $9, $10, 0, 0, FabricMat::GREEN_LAGRANGE, ginf, $9, $10, $11, $12, $13, $14));
         }
        | MatSpec Integer VISCOFABRICMAT Float Integer Integer Float Float Float Float Float Float Float Float Float Float Float Float NewLine
         {
           double gthree = $15;
           double gtwo   = $13;
           double gone   = $11;
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial($2-1,
             new PronyViscoElastic<FabricMat>($4, $5, $6, $7, $8, $9, $10, $11, $12, FabricMat::GREEN_LAGRANGE, ginf, $11, $12, $13, $14, $15, $16));
         }
        | MatSpec Integer SHELLVISCOFABRICMAT Float Integer Integer Float Float Float Float Float Float Float Float Float Float NewLine
         {
           double gthree = $13;
           double gtwo   = $11;
           double gone   = $9;
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial($2-1,
             new PronyViscoElastic<FabricMat>($4, $5, $6, $7, $8, $9, $10, 0, 0, FabricMat::INFINTESIMAL, ginf, $9, $10, $11, $12, $13, $14));
         }
        | MatSpec Integer SHELLVISCOFABRICMAT Float Integer Integer Float Float Float Float Float Float Float Float Float Float Float Float NewLine
         {
           double gthree = $15;
           double gtwo   = $13;
           double gone   = $11;
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial($2-1,
             new PronyViscoElastic<FabricMat>($4, $5, $6, $7, $8, $9, $10, $11, $12, FabricMat::INFINTESIMAL, ginf, $11, $12, $13, $14, $15, $16));
         }
        | MatSpec Integer PLANESTRESSVISCOLINEARELASTIC Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double gthree = $12;
            double gtwo   = $10;
            double gone   = $8;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PlaneStressMat<PronyViscoElastic<ElaLinIsoMat> >($4, $5, $6, $7, ginf, $8, $9, $10, $11, $12, $13));
          }
        | MatSpec Integer PLANESTRESSVISCOLINEARELASTIC Float Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double gthree = $14;
            double gtwo   = $12;
            double gone   = $10;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PlaneStressMat<PronyViscoElastic<ElaLinIsoMat> >($4, $5, $6, $7, $8, $9, ginf, $10, $11, $12, $13, $14, $15));
          }
        | MatSpec Integer PLANESTRESSVISCOLINEARELASTIC Float Float Float Float Float Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
          {
            double gthree = $12;
            double gtwo   = $10;
            double gone   = $8; 
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> > >($4, $5, $6, $7, ginf, $8, $9, $10, $11, $12, $15, $16, $17, $13));
            domain->solInfo().elementDeletion = true;
          }
        | MatSpec Integer PLANESTRESSVISCOLINEARELASTIC Float Float Float Float Float Float Float Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
          {
            double gthree = $14;
            double gtwo   = $12;
            double gone   = $10;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> > >($4, $5, $6, $7, $8, $9, ginf, $10, $11, $12, $13, $14, $17, $18, $19, $15));
            domain->solInfo().elementDeletion = true;
          }
        | MatSpec Integer VISCOSTVENANTKIRCHHOFF Float Float Float Float Float Float Float Float Float NewLine
          {
            double gthree = $11;
            double gtwo   = $9;
            double gone   = $7;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PronyViscoElastic<StVenantKirchhoffMat>($4, $5, $6, ginf, $7, $8, $9, $10, $11, $12));
          }
        | MatSpec Integer VISCOSTVENANTKIRCHHOFF Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double gthree = $13;
            double gtwo   = $11;
            double gone   = $9;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PronyViscoElastic<StVenantKirchhoffMat>($4, $5, $6, $7, $8, ginf, $9, $10, $11, $12, $13, $14));
          }
        | MatSpec Integer VISCOSTVENANTKIRCHHOFF Float Float Float Float Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
          {
            double gthree = $11;
            double gtwo   = $9;
            double gone   = $7;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> >($4, $5, $6, ginf, $7, $8, $9, $10, $11, $12, $14, $15, $16));
            domain->solInfo().elementDeletion = true;
          }
        | MatSpec Integer VISCOSTVENANTKIRCHHOFF Float Float Float Float Float Float Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
          {
            double gthree = $13;
            double gtwo   = $11;
            double gone   = $9;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> >($4, $5, $6, $7, $8, ginf, $9, $10, $11, $12, $13, $14, $16, $17, $18));
            domain->solInfo().elementDeletion = true;
          }
        | MatSpec Integer PLANESTRESSVISCOSTVENANTKIRCHHOFF Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double gthree = $12;
            double gtwo   = $10;
            double gone   = $8;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PlaneStressMat<PronyViscoElastic<StVenantKirchhoffMat> >($4, $5, $6, $7, ginf, $8, $9, $10, $11, $12, $13));
          }
        | MatSpec Integer PLANESTRESSVISCOSTVENANTKIRCHHOFF Float Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double gthree = $14;
            double gtwo   = $12;
            double gone   = $10;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PlaneStressMat<PronyViscoElastic<StVenantKirchhoffMat> >($4, $5, $6, $7, $8, $9, ginf, $10, $11, $12, $13, $14, $15));
          }
        | MatSpec Integer PLANESTRESSVISCOSTVENANTKIRCHHOFF Float Float Float Float Float Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
          {
            double gthree = $12;
            double gtwo   = $10;
            double gone   = $8;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> > >($4, $5, $6, $7, ginf, $8, $9, $10, $11, $12, $15, $16, $17, $13));
            domain->solInfo().elementDeletion = true;
          }
        | MatSpec Integer PLANESTRESSVISCOSTVENANTKIRCHHOFF Float Float Float Float Float Float Float Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
          {
            double gthree = $14;
            double gtwo   = $12;
            double gone   = $10;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> > >($4, $5, $6, $7, $8, $9, ginf, $10, $11, $12, $13, $14, $17, $18, $19, $15));
            domain->solInfo().elementDeletion = true;
          }
        | MatSpec Integer NEOHOOKEAN Float Float Float NewLine
          {
            geoSource->addMaterial($2-1,
              new NeoHookeanMat($4, $5, $6));
          }
        | MatSpec Integer NEOHOOKEAN Float Float Float TULERBUTCHER Float Float Float NewLine
          {
            geoSource->addMaterial($2-1,
              new BrittleFractureTB<NeoHookeanMat>($4, $5, $6, $8, $9, $10));
            domain->solInfo().elementDeletion = true;
          }
        | MatSpec Integer VISCONEOHOOKEAN Float Float Float Float Float Float Float Float Float NewLine
          {
            double gthree = $11;
            double gtwo   = $9;
            double gone   = $7;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PronyViscoElastic<NeoHookeanMat>($4, $5, $6, ginf, $7, $8, $9, $10, $11, $12));
          }
        | MatSpec Integer VISCONEOHOOKEAN Float Float Float Float Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
          {
            double gthree = $11;
            double gtwo   = $9;
            double gone   = $7;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new BrittleFractureTB<PronyViscoElastic<NeoHookeanMat> >($4, $5, $6, ginf, $7, $8, $9, $10, $11, $12, $14, $15, $16));
            domain->solInfo().elementDeletion = true;
          }
        | MatSpec Integer PLANESTRESSVISCONEOHOOKEAN Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double gthree = $12;
            double gtwo   = $10;
            double gone   = $8;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PlaneStressMat<PronyViscoElastic<NeoHookeanMat> >($4, $5, $6, $7, ginf, $8, $9, $10, $11, $12, $13));
          }
        | MatSpec Integer PLANESTRESSVISCONEOHOOKEAN Float Float Float Float Float Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
          {
            double gthree = $12;
            double gtwo   = $10;
            double gone   = $8;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<NeoHookeanMat> > >($4, $5, $6, $7, ginf, $8, $9, $10, $11, $12, $15, $16, $17, $13));
            domain->solInfo().elementDeletion = true;
          }
        | MatSpec Integer MOONEYRIVLIN Float Float Float Float NewLine
          {
            geoSource->addMaterial($2-1,
              new MooneyRivlinMat($4, $5, $6, $7));
          }
        | MatSpec Integer MOONEYRIVLIN Float Float Float Float TULERBUTCHER Float Float Float NewLine
          {
            geoSource->addMaterial($2-1,
              new BrittleFractureTB<MooneyRivlinMat>($4, $5, $6, $7, $9, $10, $11));
            domain->solInfo().elementDeletion = true;
          }
        | MatSpec Integer VISCOMOONEYRIVLIN Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double gthree = $12;
            double gtwo   = $10;
            double gone   = $8;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PronyViscoElastic<MooneyRivlinMat>($4, $5, $6, $7, ginf, $8, $9, $10, $11, $12, $13));
          }
        | MatSpec Integer VISCOMOONEYRIVLIN Float Float Float Float Float Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
          {
            double gthree = $12;
            double gtwo   = $10;
            double gone   = $8;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new BrittleFractureTB<PronyViscoElastic<MooneyRivlinMat> >($4, $5, $6, $7, ginf, $8, $9, $10, $11, $12, $13, $15, $16, $17));
            domain->solInfo().elementDeletion = true;
          }
        | MatSpec Integer PLANESTRESSVISCOMOONEYRIVLIN Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double gthree = $13;
            double gtwo   = $11;
            double gone   = $9;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PlaneStressMat<PronyViscoElastic<MooneyRivlinMat> >($4, $5, $6, $7, $8, ginf, $9, $10, $11, $12, $13, $14));
          }
        | MatSpec Integer PLANESTRESSVISCOMOONEYRIVLIN Float Float Float Float Float Float Float Float Float Float Float TULERBUTCHER Float Float Float NewLine
          {
            double gthree = $13;
            double gtwo   = $11;
            double gone   = $9;
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial($2-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<MooneyRivlinMat> > >($4, $5, $6, $7, $8, ginf, $9, $10, $11, $12, $13, $16, $17, $18, $14));
            domain->solInfo().elementDeletion = true;
          }
        | MatSpec Integer OGDEN Float Float Float Float NewLine 
          {
            double mu[1] = { $5 };
            double alpha[1] = { $6 };
            double K[1] = { $7 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float NewLine
          {
            double mu[2] = { $5, $6 };
            double alpha[2] = { $7, $8 };
            double K[1] = { $9 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float Float NewLine
          {
            double mu[2] = { $5, $6 }; 
            double alpha[2] = { $7, $8 };
            double K[2] = { $9, $10 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float Float Float NewLine
          {
            double mu[3] = { $5, $6, $7 };
            double alpha[3] = { $8, $9, $10 };
            double K[1] = { $11 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float Float Float Float NewLine
          {
            double mu[3] = { $5, $6, $7 };
            double alpha[3] = { $8, $9, $10 };
            double K[2] = { $11, $12 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double mu[4] = { $5, $6, $7, $8 };
            double alpha[4] = { $9, $10, $11, $12 };
            double K[1] = { $13 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double mu[4] = { $5, $6, $7, $8 };
            double alpha[4] = { $9, $10, $11, $12 };
            double K[2] = { $13, $14 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double mu[5] = { $5, $6, $7, $8, $9 };
            double alpha[5] = { $10, $11, $12, $13, $14 };
            double K[1] = { $15 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double mu[5] = { $5, $6, $7, $8, $9 };
            double alpha[5] = { $10, $11, $12, $13, $14 };
            double K[2] = { $15, $16 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double mu[6] = { $5, $6, $7, $8, $9, $10 };
            double alpha[6] = { $11, $12, $13, $14, $15, $16 };
            double K[1] = { $17 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double mu[6] = { $5, $6, $7, $8, $9, $10 };
            double alpha[6] = { $11, $12, $13, $14, $15, $16 };
            double K[2] = { $17, $18 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double mu[7] = { $5, $6, $7, $8, $9, $10, $11 };
            double alpha[7] = { $12, $13, $14, $15, $16, $17, $18 };
            double K[1] = { $19 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double mu[7] = { $5, $6, $7, $8, $9, $10, $11 };
            double alpha[7] = { $12, $13, $14, $15, $16, $17, $18 };
            double K[2] = { $19, $20 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double mu[8] = { $5, $6, $7, $8, $9, $10, $11, $12 };
            double alpha[8] = { $13, $14, $15, $16, $17, $18, $19, $20 };
            double K[1] = { $21 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double mu[8] = { $5, $6, $7, $8, $9, $10, $11, $12 };
            double alpha[8] = { $13, $14, $15, $16, $17, $18, $19, $20 };
            double K[2] = { $21, $22 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double mu[9] = { $5, $6, $7, $8, $9, $10, $11, $12, $13 };
            double alpha[9] = { $14, $15, $16, $17, $18, $19, $20, $21, $22 };
            double K[1] = { $23 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer OGDEN Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float NewLine
          {
            double mu[9] = { $5, $6, $7, $8, $9, $10, $11, $12, $13 };
            double alpha[9] = { $14, $15, $16, $17, $18, $19, $20, $21, $22 };
            double K[2] = { $23, $24 };
            geoSource->addMaterial($2-1, new OgdenMat($4, mu, alpha, K));
          }
        | MatSpec Integer SIMOELASTIC Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new SimoElasticMat($4, $5, $6));
         }
        | MatSpec Integer SIMOPLASTIC Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new SimoPlasticMat($4, $5, $6, $7, $8) );
         }
        | MatSpec Integer SIMOPLASTIC Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new SimoPlasticMat($4, $5, $6, $7, $8, $9) );
         }
        | MatSpec Integer ISOTROPICLINEARELASTICJ2PLASTIC Float Float Float Float Float Float NewLine
          {
            double params[9] = { $4, $5, $6, $7, $8, $9, 1.0e-6, std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
          }
        | MatSpec Integer ISOTROPICLINEARELASTICJ2PLASTIC Float Float Float Float Float Float Float NewLine
          {
            double params[9] = { $4, $5, $6, $7, $8, $9, $10, std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
          }
        | MatSpec Integer ISOTROPICLINEARELASTICJ2PLASTIC Float Float Float Float Float Float Float Float NewLine
          {
            double params[9] = { $4, $5, $6, $7, $8, $9, $10, $11, 0. };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
            if($11 > 0 && $11 < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
        | MatSpec Integer ISOTROPICLINEARELASTICJ2PLASTIC Float Float Float Float Float Float Float Float Integer NewLine
          {
            double params[9] = { $4, $5, $6, $7, $8, $9, $10, $11, double($12) };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
            if($11 > 0 && $11 < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
        | MatSpec Integer ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS Float Float Float Float Float Float NewLine
          {
            double params[9] = { $4, $5, $6, $7, $8, $9, 1.0e-6, std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          }
        | MatSpec Integer ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS Float Float Float Float Float Float Float NewLine
          {
            double params[9] = { $4, $5, $6, $7, $8, $9, $10, std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          }
        | MatSpec Integer ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS Float Float Float Float Float Float Float Float NewLine
          {
            double params[9] = { $4, $5, $6, $7, $8, $9, $10, $11, 0. };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
            if($11 > 0 && $11 < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
        | MatSpec Integer ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS Float Float Float Float Float Float Float Float Integer NewLine
          {
            double params[9] = { $4, $5, $6, $7, $8, $9, $10, $11, double($12) };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
            if($11 > 0 && $11 < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
        | MatSpec Integer OPTCTV Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ExpMat($3, $4, $5, $6));
         }
        | MatSpec Integer OPTCTV Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ExpMat($3, $4, $5, $6, $7));
         }
        | MatSpec Integer OPTCTV Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ExpMat($3, $4, $5, $6, $7, $8));
         }
        | MatSpec Integer OPTCTV Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ExpMat($3, $4, $5, $6, $7, $8, $9));
         }
        | MatSpec Integer OPTCTV Float Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ExpMat($3, $4, $5, $6, $7, $8, $9, $10));
         }
        | MatSpec Integer OPTCTV Float Float Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ExpMat($3, $4, $5, $6, $7, $8, $9, $10, $11));
           if($11 > 0 && $11 < std::numeric_limits<double>::infinity()) {
             domain->solInfo().elementDeletion = true;
           }
         }
        | MatSpec Integer OPTCTV Float Float Float Float Float Float Float Float Integer NewLine
         {
           ExpMat *mat = new ExpMat($3, $4, $5, $6, $7, $8, $9, $10, $11);
           mat->yssrtid = $12;
           geoSource->addMaterial($2-1, mat);
           if($11 > 0 && $11 < std::numeric_limits<double>::infinity()) {
             domain->solInfo().elementDeletion = true;
           }
         }
        | MatSpec Integer OPTCTV Float Float Float Float Float Float Float Float Integer
                          Integer Integer Integer Integer Float Float Float Integer Float NewLine
         {
           ExpMat *mat = new ExpMat($3, $4, $5, $6, $7, $8, $9, $10, $11);
           mat->yssrtid = $12;
           mat->optcor0 = $13;
           mat->optcor1 = $14;
           mat->optprj = $15;
           mat->opthgc = $16;
           mat->prmhgc[0] = $17;
           mat->prmhgc[1] = $18;
           mat->prmhgc[2] = $19;
           mat->ngqpt2 = $20;
           mat->ematpro[18] = $21;
           geoSource->addMaterial($2-1, mat);
           if($11 > 0 && $11 < std::numeric_limits<double>::infinity()) {
             domain->solInfo().elementDeletion = true;
           }
         }
	| MatSpec READ FNAME FNAME NewLine
	 {
	   geoSource->loadMaterial($3, $4);
	 }
	| MatSpec Integer FNAME FloatList NewLine
	 {
	   geoSource->addMaterial($2-1, $3, $4);
	 }
	;
MatUsage:
	MATUSAGE NewLine
	| MatUsage Integer Integer NewLine
	  { geoSource->setMatUsage($2-1, $3-1); }
	| MatUsage Integer Integer Integer NewLine
	  {
            for(int i = $2-1; i < $3; ++i)
	      geoSource->setMatUsage(i, $4-1);
	  }
	;
FloatList:
	{ $$.nval = 0; }
	| FloatList Float
	{ 
          if($1.nval == 64) {
             fprintf(stderr, "You'd better invent another material model!\n");
	     exit(-1);
          }
          $$ = $1;
          $$.v[$$.nval++] = $2;
 	}
	;
StringList:
	{ $$.nval = 0; }
	| StringList FNAME
	{ 
          if($1.nval == 32) {
             fprintf(stderr, "Too many files!\n");
	     exit(-1);
          }
          $$ = $1;
          $$.v[$$.nval++] = $2;
 	}
	;
Renumbering:
	RENUM NewLine RENUMBERID NewLine
	{ domain->solInfo().setRenum($3);
          domain->solInfo().setSparseRenum($3); 
          domain->solInfo().setSpoolesRenum($3); }
        | RENUM NewLine RENUMBERID NewLine RENUMBERID NewLine
        { domain->solInfo().setRenum($3);
          domain->solInfo().setSparseRenum($5); }
        | RENUM NewLine RENUMBERID NewLine RENUMBERID NewLine RENUMBERID NewLine
        { domain->solInfo().setRenum($3);
          domain->solInfo().setSparseRenum($5); 
          domain->solInfo().setSpoolesRenum($7); }
	;

SvdToken:
    SVDTOKEN NewLine
  { domain->solInfo().activatePodRom = true; 
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().svdPodRom = true;}
  | SvdToken SvdOption NewLine
  ;

SvdOption:
  SNAPFI StringList
  { for(int i=0; i<$2.nval; ++i) domain->solInfo().snapfiPodRom.push_back(std::string($2.v[i])); }
  /*
  | SNAPFI FNAME Integer
  { domain->solInfo().snapfiPodRom.push_back($2);
    if ($3 == 1) domain->solInfo().statevectPodRom = true;
    if ($3 == 2) domain->solInfo().residvectPodRom = true;
    if ($3 == 3) domain->solInfo().jacobvectPodRom = true;
    if ($3 == 4) domain->solInfo().forcevectPodRom = true;
    if ($3 == 5) domain->solInfo().accelvectPodRom = true;}
    */
  | VELSNAPFI StringList
  { for(int i=0; i<$2.nval; ++i) domain->solInfo().velocPodRomFile.push_back(std::string($2.v[i])); }
  | ACCSNAPFI StringList
  { for(int i=0; i<$2.nval; ++i) domain->solInfo().accelPodRomFile.push_back(std::string($2.v[i])); }
  | DSVSNAPFI StringList
  { for(int i=0; i<$2.nval; ++i) domain->solInfo().dsvPodRomFile.push_back(std::string($2.v[i])); }
  | MUVSNAPFI StringList
  { for(int i=0; i<$2.nval; ++i) domain->solInfo().muvPodRomFile.push_back(std::string($2.v[i])); }
  | PODSIZEMAX Integer
  { domain->solInfo().maxSizePodRom = $2; }
  | NORMALIZETOKEN Integer
  { domain->solInfo().normalize = $2; }
  | NORMALIZETOKEN SWITCH
  { domain->solInfo().normalize = $2; }
  | NORMALIZETOKEN FLAGMN
  { domain->solInfo().normalize = $2; }
  | REFSUBTRACT FNAME
  { domain->solInfo().subtractRefPodRom = true;
    domain->solInfo().readInLocalBasesCent.push_back(std::string($2)); }
  | SNAPWEIGHT SWITCH
  { domain->solInfo().flagss = $2; }
  | SNAPWEIGHT SWITCH SWITCH
  { domain->solInfo().flagss = $2; domain->solInfo().flagrs = $3; }
  | SKIP Integer
  { domain->solInfo().skipPodRom = $2; } 
  | SKIP Integer Integer
  { domain->solInfo().skipPodRom = $2;
    domain->solInfo().skipOffSet = $3; }
  | ROBCSOLVE SWITCH
  { domain->solInfo().robcSolve = bool($2); }
  | ROBFI StringList
  { for(int i=0; i<$2.nval; ++i) domain->solInfo().robfi.push_back(std::string($2.v[i])); }
  | BLOCKSIZE Integer
  { domain->solInfo().svdBlockSize = $2; }
  | ROMENERGY Float
  { domain->solInfo().romEnergy = $2; }
  /*
  | USEPQN Integer Integer Integer Integer Float Integer Float
  { domain->solInfo().use_nmf = 3;
    domain->solInfo().nmfNumROBDim = $2;
    domain->solInfo().nmfDelROBDim = $3;
    domain->solInfo().nmfRandInit = $4;
    domain->solInfo().nmfMaxIter = $5;
    domain->solInfo().nmfTol = $6;
    domain->solInfo().nmfPqnNumInnerIter = $7;
    domain->solInfo().nmfPqnAlpha = $8; }
  */
  | USEPQN Integer Float Integer Float
  { domain->solInfo().use_nmf            = 3;
    domain->solInfo().nmfMaxIter         = $2;
    domain->solInfo().nmfTol             = $3;
    domain->solInfo().nmfPqnNumInnerIter = $4;
    domain->solInfo().nmfPqnAlpha        = $5; }
  | USENMF Integer Integer Integer Integer Float
  { domain->solInfo().use_nmf = 1;
    domain->solInfo().nmfNumROBDim = $2;
    domain->solInfo().nmfDelROBDim = $3;
    domain->solInfo().nmfRandInit  = $4;
    domain->solInfo().nmfMaxIter   = $5;
    domain->solInfo().nmfTol = $6; }
  | USENMF Integer Float
  { domain->solInfo().use_nmf    = 1;
    domain->solInfo().nmfMaxIter = $2;
    domain->solInfo().nmfTol     = $3; }
  | USENMFC Integer Float Float Float Float
  { domain->solInfo().use_nmf    = 4;   
    domain->solInfo().nmfMaxIter = $2; 
    domain->solInfo().nmfTol     = $3; 
    domain->solInfo().nmfcAlpha  = $4;
    domain->solInfo().nmfcBeta   = $5;
    domain->solInfo().nmfcGamma  = $6;}
  | USENMFC Integer Float
  { domain->solInfo().use_nmf    = 4;
    domain->solInfo().nmfMaxIter = $2;
    domain->solInfo().nmfTol     = $3; }
  | NSUBS Integer
  { domain->solInfo().nmfNumSub = $2; }
  | USEGREEDY
  { domain->solInfo().use_nmf = 2; }
  | DOCLUSTERING Integer
  { domain->solInfo().clustering = $2; }
  | DOCLUSTERING Integer ANGLE
  { domain->solInfo().clustering = $2; 
    domain->solInfo().clusterSubspaceAngle = true; }
  | CLUSTERSOLVER CLUSTERSOLVERTYPE
  { domain->solInfo().solverTypeCluster = $2; }
  | CLUSTERSOLVER CLUSTERSOLVERTYPE Float
  { domain->solInfo().solverTypeCluster = $2;
    domain->solInfo().tolPodRom = $3;}
  | CLUSTERSOLVER CLUSTERSOLVERTYPE Float SPNNLSSOLVERTYPE
  { domain->solInfo().solverTypeCluster = $2;
    domain->solInfo().tolPodRom = $3;
    domain->solInfo().solverTypeSpnnls = $4; }
  | HOTSTART SWITCH
  { domain->solInfo().hotstartSample = bool($2); }
  | DOROWCLUSTERING Integer
  { domain->solInfo().rowClustering = $2; }
  | DOROWCLUSTERING Integer ANGLE
  { domain->solInfo().rowClustering = $2;
    domain->solInfo().clusterSubspaceAngle = true; }
  | ConwepConfig
  ;

DeimIndices:
   DEIMINDICES NewLine
   { domain->solInfo().activatePodRom = true;
     domain->solInfo().setProbType(SolverInfo::PodRomOffline);
     domain->solInfo().DEIMBasisPod = true; } 
   | DeimIndices SamplingOption NewLine
   ;

UDeimIndices:
   UDEIMINDICES NewLine
   { domain->solInfo().activatePodRom = true;
     domain->solInfo().setProbType(SolverInfo::PodRomOffline);
     domain->solInfo().UDEIMBasisPod = true; }
   | UDeimIndices SamplingOption NewLine
   ;

Sampling:
    SAMPLING NewLine 
  { domain->solInfo().activatePodRom = true; 
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().samplingPodRom = true; }
  | Sampling SamplingOption NewLine
  | Sampling ConwepConfig
  | Sampling ScalePosCoords
  | Sampling NodePosCoords
  ;

SnapshotProject:
    SNAPSHOTPROJECT NewLine 
  { domain->solInfo().activatePodRom = true;
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().snapProjPodRom = true; }
  | SnapshotProject SamplingOption NewLine
  ;

SamplingOption:
    PODROB FNAME
  { domain->solInfo().readInROBorModes.push_back($2); }
  | PODROB FNAME Integer
  { domain->solInfo().readInROBorModes.push_back($2);
    domain->solInfo().localBasisSize.push_back($3); }
  | DUALBASIS FNAME Integer 
  { domain->solInfo().readInDualROB.push_back($2);
    domain->solInfo().localDualBasisSize.push_back($3); }
  | TRNVCT FNAME
  { domain->solInfo().statePodRomFile.push_back($2); }
  | TRNVCT FNAME FNAME
  { domain->solInfo().statePodRomFile.push_back($2);
    domain->solInfo().velocPodRomFile.push_back($3); }
  | TRNVCT FNAME FNAME FNAME
  { domain->solInfo().statePodRomFile.push_back($2);
    domain->solInfo().velocPodRomFile.push_back($3);
    domain->solInfo().accelPodRomFile.push_back($4); }
  | TOLER Float
  { domain->solInfo().tolPodRom = $2; }
  | SKIP Integer
  { domain->solInfo().skipPodRom = $2; }
  | RANDOMSAMPLE Integer
  { domain->solInfo().randomSampleSize = $2; 
    domain->solInfo().randomVecSampling = true; }
  | OFFSET Integer
  { domain->solInfo().skipOffSet = $2; }
  | PODSIZEMAX Integer
  { domain->solInfo().maxSizePodRom = $2; }
  | PODSIZEMAX Integer Integer
  { domain->solInfo().maxSizePodRom = $2; 
    domain->solInfo().forcePodSize = $3;}
  | PODSIZEMAX Integer Integer Integer
  { domain->solInfo().maxSizePodRom = $2; 
    domain->solInfo().forcePodSize = $3;
    domain->solInfo().maxDeimBasisSize = $4; }
  | USEMASSNORMALIZEDBASIS SWITCH
  { domain->solInfo().useMassNormalizedBasis = bool($2); }
  | USECONSTANTMASS SWITCH
  { domain->solInfo().useConstantMassForces = bool($2); }
  | STACKED SWITCH
  { domain->solInfo().stackedElementSampling = bool($2); }
  | MPROJECT SWITCH
  { domain->solInfo().useMassOrthogonalProjection = bool($2); }
  | REDFOL SWITCH /* deprecated */
  { domain->solInfo().reduceFollower = bool($2); }
  | EXTFOL SWITCH
  { domain->solInfo().reduceFollower = bool($2); }
  | VECTORNORM FNAME
  { domain->solInfo().PODerrornorm.push_back($2); }
  | VECTORNORM FNAME FNAME
  { domain->solInfo().PODerrornorm.push_back($2);
    domain->solInfo().PODerrornorm.push_back($3); }
  | VECTORNORM FNAME FNAME FNAME
  { domain->solInfo().PODerrornorm.push_back($2);
    domain->solInfo().PODerrornorm.push_back($3);
    domain->solInfo().PODerrornorm.push_back($4); }
  | SCALING SWITCH
  { domain->solInfo().useScalingSpnnls = bool($2); }
  | CENTER SWITCH
  { domain->solInfo().useCenterSpnnls = bool($2); }
  | PROJSOL SWITCH
  { domain->solInfo().projectSolution = bool($2); }
  | POSELEM SWITCH
  { domain->solInfo().positiveElements = bool($2); }
  | HOTSTART SWITCH
  { domain->solInfo().hotstartSample = bool($2); }
  | SOLVER SPNNLSSOLVERTYPE
  { domain->solInfo().solverTypeSpnnls = $2; }
  | MAXSIZE Float
  { domain->solInfo().maxSizeSpnnls = $2; }
  | MAXELEM Integer
  { domain->solInfo().maxElemSpnnls = $2; }
  | MAXITR Float
  { domain->solInfo().maxIterSpnnls = $2; }
  | FORCEROB FNAME
  { domain->solInfo().forcePodRomFile = $2; }
  | FORCEROB FNAME Integer
  { domain->solInfo().forcePodRomFile = $2;
    domain->solInfo().forcePodSize = $3; }
  | FORCEROB FNAME Integer Integer
  { domain->solInfo().forcePodRomFile = $2; 
    domain->solInfo().forcePodSize = $3; 
    domain->solInfo().maxDeimBasisSize = $4; }
  | CONSTRAINTROB FNAME
  { domain->solInfo().constraintPodRomFile = $2; 
    domain->solInfo().ConstraintBasisPod = true;}
  | CONSTRAINTROB FNAME Integer
  { domain->solInfo().constraintPodRomFile = $2;
    domain->solInfo().constraintPodSize = $3; 
    domain->solInfo().ConstraintBasisPod = true; }
  | CONSTRAINTROB FNAME Integer Integer
  { domain->solInfo().constraintPodRomFile = $2;
    domain->solInfo().constraintPodSize = $3;
    domain->solInfo().maxDeimBasisSize = $4; 
    domain->solInfo().ConstraintBasisPod = true; }
  | FILTERROWS SWITCH
  { domain->solInfo().filterSnapshotRows = bool($2); }
  | PSEUDOGNAT SWITCH
  { domain->solInfo().selectFullNode = bool($2); }
  | PSEUDOGNATELEM SWITCH
  { domain->solInfo().selectFullElem = bool($2); }
  | REBUILDFORCE SWITCH 
  { domain->solInfo().computeForceSnap = bool($2); }
  | REBUILDCONSTRAINT SWITCH
  { domain->solInfo().computeConstraintSnap = bool($2); }
  | SVDFORCESNAP SWITCH
  { domain->solInfo().orthogForceSnap = bool($2); }
  | SVDCONSTRAINTSNAP SWITCH
  { domain->solInfo().orthogConstraintSnap = bool($2); }
  | NPMAX Integer
  { domain->solInfo().npMax = $2; }
  | BSSPLH Integer Integer
  { domain->solInfo().scpkMB= $2;
    domain->solInfo().scpkNB= $3; }
  | PGSPLH Integer Integer
  { domain->solInfo().scpkMP= $2;
    domain->solInfo().scpkNP= $3; }
  | REVERSEORDER SWITCH
  { domain->solInfo().useReverseOrder = bool($2); }
  | USENMF Integer Integer Integer Integer Float
  { domain->solInfo().use_nmf = 1;
    domain->solInfo().nmfNumROBDim = $2;
    domain->solInfo().nmfDelROBDim = $3;
    domain->solInfo().nmfRandInit = $4;
    domain->solInfo().nmfMaxIter = $5;
    domain->solInfo().nmfTol = $6; }
  | USENMF Integer Float
  { domain->solInfo().use_nmf = 1;
    domain->solInfo().nmfMaxIter = $2;
    domain->solInfo().nmfTol = $3; }
  ;

ConwepConfig:
   CONWEPCFG NewLine
   | ConwepConfig ConwepData NewLine
   { domain->solInfo().conwepConfigurations.push_back($2); }
  ;

MeshScaleFactor:
     MESHSCALEFACTOR Float Float Float NewLine
   { domain->solInfo().scalePosCoords = true;
     domain->solInfo().xScaleFactor = $2;
     domain->solInfo().yScaleFactor = $3;
     domain->solInfo().zScaleFactor = $4;}
   ;


ScalePosCoords:
   SCALEPOSCOORDS NewLine
   { domain->solInfo().activatePOSCFG = true; } 
   | ScalePosCoords Float Float Float NewLine
   { domain->solInfo().xScaleFactors.push_back($2);
     domain->solInfo().yScaleFactors.push_back($3);
     domain->solInfo().zScaleFactors.push_back($4); }
   | ScalePosCoords FNAME NewLine
   { domain->solInfo().MassOrthogonalBasisFiles.push_back($2); }
   ;

NodePosCoords:
   NODEPOSCOORDS NewLine
   { domain->solInfo().activatePOSCFG = true; }
   | NodePosCoords StringList NewLine
   { domain->solInfo().NodeTrainingFiles.push_back(std::string($2.v[0]));
     for(int i=1; i<$2.nval; ++i) domain->solInfo().MassOrthogonalBasisFiles.push_back(std::string($2.v[i])); }
   ;
   
ConversionToken:
    CONVERSIONTOKEN NewLine
  { domain->solInfo().activatePodRom = true;
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().ROMPostProcess = true; }
  | ConversionToken ConversionOption NewLine
  ;

ConversionOption:
    CONVFI FNAME
  { domain->solInfo().RODConversionFiles.push_back($2); 
    domain->solInfo().numRODFile += 1; }
  | CONVFI FNAME Integer
  { domain->solInfo().RODConversionFiles.push_back($2);
    domain->solInfo().numRODFile += 1; 
    domain->solInfo().skipPodRom = $3;}
  | ROMRES Integer
  { domain->solInfo().romresidType = $2; }
  ;

Integer:
	IntConstant
	{ $$ = $1; }
        | MAXIMUM
        { $$ = std::numeric_limits<int>::max(); }
	;

Float:
	IntConstant
	{ $$ = $1; }
	| DblConstant 
	{ $$ = $1; }
        | INFINTY
        { $$ = std::numeric_limits<double>::infinity(); }
        | EPSILON
        { $$ = std::numeric_limits<double>::epsilon(); }
	;
%%
