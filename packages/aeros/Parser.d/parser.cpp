/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 1 "p.y" /* yacc.c:339  */

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

#line 112 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "parser.hpp".  */
#ifndef YY_YY_HOME_ANARKHEDE_TINKERCLIFFS_FEMWORKINGFOAM_PARSER_D_PARSER_HPP_INCLUDED
# define YY_YY_HOME_ANARKHEDE_TINKERCLIFFS_FEMWORKINGFOAM_PARSER_D_PARSER_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    ACTUATORS = 258,
    ADJOINTBASIS = 259,
    AERO = 260,
    AEROH = 261,
    AEROTYPE = 262,
    AGRESSIVETOLERANCES = 263,
    ALPROC = 264,
    AMAT = 265,
    ANALYSIS = 266,
    ARCLENGTH = 267,
    ATTRIBUTES = 268,
    ANGULAROUTTYPE = 269,
    ARUBBERMAT = 270,
    AUGMENT = 271,
    AUGMENTTYPE = 272,
    AUTOTOL = 273,
    AUXILIARY = 274,
    AVERAGED = 275,
    ATDARB = 276,
    ACOU = 277,
    ATDDNB = 278,
    ATDROB = 279,
    ARPACK = 280,
    ATDDIR = 281,
    ATDNEU = 282,
    ALLOWMECHANISMS = 283,
    AUXCOARSESOLVER = 284,
    ACMECNTL = 285,
    ADDEDMASS = 286,
    AEROEMBED = 287,
    ANDESCLR = 288,
    ANDESCQR = 289,
    ANDESBETAB = 290,
    ANDESALPHA = 291,
    ANDESBETAM = 292,
    AUGMENTED = 293,
    BLOCKDIAG = 294,
    BOFFSET = 295,
    BUCKLE = 296,
    BGTL = 297,
    BMPC = 298,
    BINARYINPUT = 299,
    BINARYOUTPUT = 300,
    BLOCKSIZE = 301,
    CHECKTOKEN = 302,
    COARSESOLVER = 303,
    COEF = 304,
    CFRAMES = 305,
    COLLOCATEDTYPE = 306,
    CONVECTION = 307,
    COMPOSITE = 308,
    CONDITION = 309,
    CONTACT = 310,
    CONTROL = 311,
    CORNER = 312,
    CORNERTYPE = 313,
    CURVE = 314,
    CCTTOL = 315,
    CCTSOLVER = 316,
    CRHS = 317,
    COUPLEDSCALE = 318,
    CONTACTSURFACES = 319,
    CMPC = 320,
    CNORM = 321,
    COMPLEXOUTTYPE = 322,
    CONSTRMAT = 323,
    CASES = 324,
    CONSTRAINEDSURFACES = 325,
    CSFRAMES = 326,
    CSTYPE = 327,
    CONSTANT = 328,
    CONWEP = 329,
    DAMPING = 330,
    DblConstant = 331,
    DEFAULTPENALTY = 332,
    DELETEELEMENTS = 333,
    DEM = 334,
    DIMASS = 335,
    DISP = 336,
    DIRECT = 337,
    DLAMBDA = 338,
    DP = 339,
    DYNAM = 340,
    DETER = 341,
    DECOMPOSE = 342,
    DECOMPFILE = 343,
    DMPC = 344,
    DEBUGCNTL = 345,
    DEBUGICNTL = 346,
    DOCLUSTERING = 347,
    DOROWCLUSTERING = 348,
    ANGLE = 349,
    DUALBASIS = 350,
    DUALRB = 351,
    KMEANS = 352,
    CRANDOM = 353,
    CONSTRAINTS = 354,
    MULTIPLIERS = 355,
    PENALTY = 356,
    ELLUMP = 357,
    EIGEN = 358,
    EFRAMES = 359,
    ELSCATTERER = 360,
    END = 361,
    ELHSOMMERFELD = 362,
    ENGINEERING = 363,
    ETEMP = 364,
    EXPLICIT = 365,
    EXTFOL = 366,
    EPSILON = 367,
    ELEMENTARYFUNCTIONTYPE = 368,
    FABMAT = 369,
    FABRICMAP = 370,
    FABRICMAT = 371,
    FACE = 372,
    FACOUSTICS = 373,
    FETI = 374,
    FETI2TYPE = 375,
    FETIPREC = 376,
    FFP = 377,
    FFPDIR = 378,
    FITALG = 379,
    FNAME = 380,
    FLAGMN = 381,
    FLUX = 382,
    FORCE = 383,
    FRONTAL = 384,
    FETIH = 385,
    FIELDWEIGHTLIST = 386,
    FILTEREIG = 387,
    FLUID = 388,
    FREEPLAY = 389,
    FREQSWEEP = 390,
    FREQSWEEP1 = 391,
    FREQSWEEP2 = 392,
    FREQSWEEPA = 393,
    FREQSWEEPAW = 394,
    FSGL = 395,
    FSINTERFACE = 396,
    FSISCALING = 397,
    FSIELEMENT = 398,
    NOLOCALFSISPLITING = 399,
    FSICORNER = 400,
    FFIDEBUG = 401,
    FAILSAFE = 402,
    FRAMETYPE = 403,
    GEPS = 404,
    GLOBALSEARCHCULL = 405,
    GLOBALTOL = 406,
    GRAVITY = 407,
    GRBM = 408,
    GTGSOLVER = 409,
    GLOBALCRBMTOL = 410,
    GROUP = 411,
    GROUPTYPE = 412,
    GOLDFARBTOL = 413,
    HDIRICHLET = 414,
    HEAT = 415,
    HFETI = 416,
    HNEUMAN = 417,
    HSOMMERFELD = 418,
    HFTT = 419,
    HELMHOLTZ = 420,
    HNBO = 421,
    HELMMF = 422,
    HELMSO = 423,
    HSCBO = 424,
    HWIBO = 425,
    HZEM = 426,
    HZEMFILTER = 427,
    HLMPC = 428,
    HERMITIAN = 429,
    HESSIAN = 430,
    IACC = 431,
    IDENTITY = 432,
    IDIS = 433,
    IDIS6 = 434,
    ILUDROPTOL = 435,
    IntConstant = 436,
    INTERFACELUMPED = 437,
    ITEMP = 438,
    ITERTYPE = 439,
    IVEL = 440,
    IMESH = 441,
    INCIDENCE = 442,
    IHDIRICHLET = 443,
    IHDSWEEP = 444,
    IHNEUMANN = 445,
    ISOLVERTYPE = 446,
    INPC = 447,
    INFINTY = 448,
    JACOBI = 449,
    KEYLETTER = 450,
    KRYLOVTYPE = 451,
    KIRLOC = 452,
    LAYC = 453,
    LAYN = 454,
    LAYD = 455,
    LAYO = 456,
    LAYMAT = 457,
    LFACTOR = 458,
    LISRBM = 459,
    LMPC = 460,
    LOAD = 461,
    LOADCASE = 462,
    LOBPCG = 463,
    LOCALREDUCEDORDERBASES = 464,
    LOCALSOLVER = 465,
    LINESEARCH = 466,
    LUMPED = 467,
    KSPARAM = 468,
    KSMAX = 469,
    MASS = 470,
    MASSAUGMENTATION = 471,
    MATERIALS = 472,
    MATLAB = 473,
    MAXITR = 474,
    MAXELEM = 475,
    MAXORTHO = 476,
    MAXVEC = 477,
    MODAL = 478,
    MPCPRECNO = 479,
    MPCPRECNOID = 480,
    MPCTYPE = 481,
    MPCTYPEID = 482,
    MPCSCALING = 483,
    MPCELEMENT = 484,
    MPCBLOCKID = 485,
    MPCBLK_OVERLAP = 486,
    MFTT = 487,
    MRHS = 488,
    MPCCHECK = 489,
    MUMPSICNTL = 490,
    MUMPSCNTL = 491,
    MUMPSMINEQ = 492,
    MUMPSSTRIDE = 493,
    MECH = 494,
    MODDAMP = 495,
    MODEFILTER = 496,
    MOMENTTYPE = 497,
    MPROJECT = 498,
    MAXIMUM = 499,
    NOWARPEDVOLUME = 500,
    NDTYPE = 501,
    NEIGPA = 502,
    NEWMARK = 503,
    NewLine = 504,
    NEWTON = 505,
    NL = 506,
    NLMAT = 507,
    NLPREC = 508,
    NLMEMPTYPE = 509,
    NOCOARSE = 510,
    NODETOKEN = 511,
    NOMULTIPLEINTERACTIONS = 512,
    NONINPC = 513,
    NSBSPV = 514,
    NLTOL = 515,
    NUMCGM = 516,
    NONORMALSMOOTHING = 517,
    NOSECONDARY = 518,
    NOGHOSTING = 519,
    NFRAMES = 520,
    NORMALSMOOTHINGDISTANCE = 521,
    SENSITIVITY = 522,
    SENSITIVITYMETHOD = 523,
    SKIPPHYSICALFACES = 524,
    OUTPUT = 525,
    OUTPUT6 = 526,
    OUTPUTFRAME = 527,
    OLDDYNAMICSEARCH = 528,
    QSTATIC = 529,
    QLOAD = 530,
    QUASISTATIC = 531,
    PITA = 532,
    PITADISP6 = 533,
    PITAVEL6 = 534,
    NOFORCE = 535,
    MDPITA = 536,
    GLOBALBASES = 537,
    LOCALBASES = 538,
    TIMEREVERSIBLE = 539,
    REMOTECOARSE = 540,
    ORTHOPROJTOL = 541,
    READINITSEED = 542,
    JUMPCVG = 543,
    JUMPOUTPUT = 544,
    PRECNO = 545,
    PRECONDITIONER = 546,
    PRELOAD = 547,
    PRESSURE = 548,
    PRINTMATLAB = 549,
    printmatlab = 550,
    PRINTNUMBER = 551,
    PROJ = 552,
    PIVOT = 553,
    PRECTYPE = 554,
    PRECTYPEID = 555,
    PICKANYCORNER = 556,
    PADEPIVOT = 557,
    PROPORTIONING = 558,
    PLOAD = 559,
    PADEPOLES = 560,
    POINTSOURCE = 561,
    PLANEWAVE = 562,
    PTOL = 563,
    PLANTOL = 564,
    PMAXIT = 565,
    PIECEWISE = 566,
    PARAMETERS = 567,
    RADIATION = 568,
    RBMFILTER = 569,
    RBMSET = 570,
    READMODE = 571,
    READSENSITIVITY = 572,
    REBUILD = 573,
    RESOLUTIONMETHOD = 574,
    REVERSEORDER = 575,
    REDFOL = 576,
    RENUM = 577,
    RENUMBERID = 578,
    REORTHO = 579,
    RESTART = 580,
    RECONS = 581,
    RECONSALG = 582,
    REBUILDCCT = 583,
    RANDOM = 584,
    RPROP = 585,
    RNORM = 586,
    REVERSENORMALS = 587,
    ROBTYPE = 588,
    ROTVECOUTTYPE = 589,
    RESCALING = 590,
    RUBDAFT = 591,
    SCALING = 592,
    SCALINGTYPE = 593,
    SDETAFT = 594,
    SENSORS = 595,
    SHELLFABRICMAP = 596,
    SHELLFABRICMAT = 597,
    SHELLSIMPLELOFTING = 598,
    SOLVERCNTL = 599,
    SOLVERHANDLE = 600,
    SOLVERTYPE = 601,
    SHIFT = 602,
    SHARPNONSHARPANGLE = 603,
    SPOOLESTAU = 604,
    SPOOLESSEED = 605,
    SPOOLESMAXSIZE = 606,
    SPOOLESMAXDOMAINSIZE = 607,
    SPOOLESMAXZEROS = 608,
    SPOOLESMSGLVL = 609,
    SPOOLESSCALE = 610,
    SPOOLESPIVOT = 611,
    SPOOLESRENUM = 612,
    SPARSEMAXSUP = 613,
    SPARSEDEFBLK = 614,
    SPARSERENUM = 615,
    STATS = 616,
    STRESSID = 617,
    SUBCYCLE = 618,
    SUBSPACE = 619,
    SURFACE = 620,
    STR_THERM_OPTION = 621,
    SAVEMEMCOARSE = 622,
    SPACEDIMENSION = 623,
    SCATTERER = 624,
    STAGTOL = 625,
    SCALED = 626,
    SWITCH = 627,
    STABLE = 628,
    SUBTYPE = 629,
    STEP = 630,
    SOWER = 631,
    SHELLTHICKNESS = 632,
    SURF = 633,
    SPRINGMAT = 634,
    SENSITIVITYID = 635,
    TABLE = 636,
    TANGENT = 637,
    TDENFORCE = 638,
    TEMP = 639,
    TIME = 640,
    TOLEIG = 641,
    TOLFETI = 642,
    TOLJAC = 643,
    TOLPCG = 644,
    TOLSEN = 645,
    TOPFILE = 646,
    TOPOLOGY = 647,
    TRBM = 648,
    TRBMlc = 649,
    THERMOE = 650,
    THERMOH = 651,
    RATIOTOLSEN = 652,
    TETT = 653,
    TOLCGM = 654,
    TURKEL = 655,
    TIEDSURFACES = 656,
    THETA = 657,
    PROJSOL = 658,
    CENTER = 659,
    POSELEM = 660,
    HOTSTART = 661,
    HRC = 662,
    THIRDNODE = 663,
    THERMMAT = 664,
    TDENFORC = 665,
    TESTULRICH = 666,
    THRU = 667,
    TRIVIAL = 668,
    NUMROMCPUS = 669,
    THICKNESSGROUPLIST = 670,
    USE = 671,
    USERDEFINEDISP = 672,
    USERDEFINEFORCE = 673,
    UPROJ = 674,
    UNSYMMETRIC = 675,
    USING = 676,
    USESCOTCH = 677,
    VERBOSE = 678,
    VERSION = 679,
    WETCORNERS = 680,
    YMTT = 681,
    SS1DT = 682,
    SS2DT = 683,
    YMST = 684,
    YSST = 685,
    YSSRT = 686,
    ZERO = 687,
    BINARY = 688,
    GEOMETRY = 689,
    DECOMPOSITION = 690,
    GLOBAL = 691,
    MATCHER = 692,
    CPUMAP = 693,
    NODALCONTACT = 694,
    MODE = 695,
    FRIC = 696,
    GAP = 697,
    OUTERLOOP = 698,
    EDGEWS = 699,
    WAVETYPE = 700,
    ORTHOTOL = 701,
    IMPE = 702,
    FREQ = 703,
    DPH = 704,
    WAVEMETHOD = 705,
    MATSPEC = 706,
    MATUSAGE = 707,
    BILINEARPLASTIC = 708,
    FINITESTRAINPLASTIC = 709,
    CRUSHABLEFOAM = 710,
    LINEARELASTIC = 711,
    STVENANTKIRCHHOFF = 712,
    TULERBUTCHER = 713,
    LINPLSTRESS = 714,
    READ = 715,
    OPTCTV = 716,
    ISOTROPICLINEARELASTIC = 717,
    VISCOLINEARELASTIC = 718,
    VISCOSTVENANTKIRCHHOFF = 719,
    NEOHOOKEAN = 720,
    VISCONEOHOOKEAN = 721,
    ISOTROPICLINEARELASTICJ2PLASTIC = 722,
    ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS = 723,
    HYPERELASTIC = 724,
    MOONEYRIVLIN = 725,
    VISCOMOONEYRIVLIN = 726,
    HENCKY = 727,
    OGDEN = 728,
    SIMOELASTIC = 729,
    SIMOPLASTIC = 730,
    LOGSTRAINPLASTIC = 731,
    SVKPLSTRESS = 732,
    VISCOLINPLSTRESS = 733,
    VISCOSVKPLSTRESS = 734,
    VISCOFABRICMAP = 735,
    SHELLVISCOFABRICMAP = 736,
    VISCOFABRICMAT = 737,
    SHELLVISCOFABRICMAT = 738,
    PARTITIONGAP = 739,
    PLANESTRESSLINEAR = 740,
    PLANESTRESSSTVENANTKIRCHHOFF = 741,
    PLANESTRESSNEOHOOKEAN = 742,
    PLANESTRESSMOONEYRIVLIN = 743,
    PLANESTRESSBILINEARPLASTIC = 744,
    PLANESTRESSFINITESTRAINPLASTIC = 745,
    PLANESTRESSVISCOLINEARELASTIC = 746,
    PLANESTRESSVISCOSTVENANTKIRCHHOFF = 747,
    PLANESTRESSVISCONEOHOOKEAN = 748,
    PLANESTRESSVISCOMOONEYRIVLIN = 749,
    SURFACETOPOLOGY = 750,
    MORTARTIED = 751,
    MORTARSCALING = 752,
    MORTARINTEGRATIONRULE = 753,
    SEARCHTOL = 754,
    STDMORTAR = 755,
    DUALMORTAR = 756,
    WETINTERFACE = 757,
    NSUBS = 758,
    EXITAFTERDEC = 759,
    SKIP = 760,
    ROBCSOLVE = 761,
    RANDOMSAMPLE = 762,
    OUTPUTMEMORY = 763,
    OUTPUTWEIGHT = 764,
    SOLVER = 765,
    SPNNLSSOLVERTYPE = 766,
    MAXSIZE = 767,
    CLUSTERSOLVER = 768,
    CLUSTERSOLVERTYPE = 769,
    WEIGHTLIST = 770,
    GMRESRESIDUAL = 771,
    SLOSH = 772,
    SLGRAV = 773,
    SLZEM = 774,
    SLZEMFILTER = 775,
    PDIR = 776,
    HEFSB = 777,
    HEFRS = 778,
    HEINTERFACE = 779,
    SNAPFI = 780,
    VELSNAPFI = 781,
    ACCSNAPFI = 782,
    DSVSNAPFI = 783,
    MUVSNAPFI = 784,
    PODROB = 785,
    ROMENERGY = 786,
    TRNVCT = 787,
    OFFSET = 788,
    ORTHOG = 789,
    SVDTOKEN = 790,
    CONVERSIONTOKEN = 791,
    CONVFI = 792,
    ROMRES = 793,
    SAMPLING = 794,
    SNAPSHOTPROJECT = 795,
    PODSIZEMAX = 796,
    REFSUBTRACT = 797,
    TOLER = 798,
    NORMALIZETOKEN = 799,
    FNUMBER = 800,
    SNAPWEIGHT = 801,
    ROBFI = 802,
    STAVCT = 803,
    VELVCT = 804,
    ACCVCT = 805,
    CONWEPCFG = 806,
    SCALEPOSCOORDS = 807,
    NODEPOSCOORDS = 808,
    MESHSCALEFACTOR = 809,
    PSEUDOGNAT = 810,
    PSEUDOGNATELEM = 811,
    USENMF = 812,
    USENMFC = 813,
    USEGREEDY = 814,
    USEPQN = 815,
    FILTERROWS = 816,
    VECTORNORM = 817,
    REBUILDFORCE = 818,
    REBUILDCONSTRAINT = 819,
    SAMPNODESLOT = 820,
    REDUCEDSTIFFNESS = 821,
    UDEIMBASIS = 822,
    FORCEROB = 823,
    CONSTRAINTROB = 824,
    DEIMINDICES = 825,
    UDEIMINDICES = 826,
    SVDFORCESNAP = 827,
    SVDCONSTRAINTSNAP = 828,
    USEMASSNORMALIZEDBASIS = 829,
    USECONSTANTMASS = 830,
    ONLINEMASSNORMALIZEBASIS = 831,
    STACKED = 832,
    NUMTHICKNESSGROUP = 833,
    STRESSNODELIST = 834,
    DISPNODELIST = 835,
    RELAXATIONSEN = 836,
    QRFACTORIZATION = 837,
    QMATRIX = 838,
    RMATRIX = 839,
    XMATRIX = 840,
    EIGENVALUE = 841,
    NPMAX = 842,
    BSSPLH = 843,
    PGSPLH = 844,
    LIB = 845
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 48 "p.y" /* yacc.c:355  */

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

#line 783 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_HOME_ANARKHEDE_TINKERCLIFFS_FEMWORKINGFOAM_PARSER_D_PARSER_HPP_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 800 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  542
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   8577

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  591
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  260
/* YYNRULES -- Number of rules.  */
#define YYNRULES  1430
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  3800

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   845

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint16 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   130,   131,   132,   133,   134,
     135,   136,   137,   138,   139,   140,   141,   142,   143,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160,   161,   162,   163,   164,
     165,   166,   167,   168,   169,   170,   171,   172,   173,   174,
     175,   176,   177,   178,   179,   180,   181,   182,   183,   184,
     185,   186,   187,   188,   189,   190,   191,   192,   193,   194,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,   207,   208,   209,   210,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,   239,   240,   241,   242,   243,   244,
     245,   246,   247,   248,   249,   250,   251,   252,   253,   254,
     255,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
     395,   396,   397,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,   412,   413,   414,
     415,   416,   417,   418,   419,   420,   421,   422,   423,   424,
     425,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     435,   436,   437,   438,   439,   440,   441,   442,   443,   444,
     445,   446,   447,   448,   449,   450,   451,   452,   453,   454,
     455,   456,   457,   458,   459,   460,   461,   462,   463,   464,
     465,   466,   467,   468,   469,   470,   471,   472,   473,   474,
     475,   476,   477,   478,   479,   480,   481,   482,   483,   484,
     485,   486,   487,   488,   489,   490,   491,   492,   493,   494,
     495,   496,   497,   498,   499,   500,   501,   502,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,   517,   518,   519,   520,   521,   522,   523,   524,
     525,   526,   527,   528,   529,   530,   531,   532,   533,   534,
     535,   536,   537,   538,   539,   540,   541,   542,   543,   544,
     545,   546,   547,   548,   549,   550,   551,   552,   553,   554,
     555,   556,   557,   558,   559,   560,   561,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,   578,   579,   580,   581,   582,   583,   584,
     585,   586,   587,   588,   589,   590
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   196,   196,   203,   204,   207,   208,   210,   212,   213,
     214,   215,   216,   217,   218,   219,   220,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   234,   235,   237,
     238,   239,   240,   241,   242,   244,   245,   246,   247,   248,
     249,   252,   253,   254,   255,   256,   257,   258,   259,   260,
     261,   262,   263,   264,   265,   266,   268,   270,   271,   272,
     273,   274,   275,   276,   277,   278,   279,   280,   281,   282,
     283,   284,   285,   286,   287,   288,   289,   290,   291,   292,
     293,   294,   295,   296,   297,   298,   299,   300,   301,   302,
     304,   306,   308,   310,   312,   314,   315,   316,   317,   318,
     319,   320,   321,   322,   323,   324,   325,   326,   327,   328,
     329,   330,   332,   333,   335,   336,   337,   338,   340,   342,
     343,   344,   345,   346,   347,   348,   350,   351,   352,   353,
     354,   356,   357,   358,   359,   361,   363,   365,   367,   369,
     371,   373,   375,   376,   377,   378,   379,   380,   381,   382,
     383,   384,   385,   388,   395,   401,   402,   407,   419,   425,
     426,   430,   432,   434,   436,   440,   444,   448,   470,   494,
     527,   528,   529,   530,   531,   534,   546,   550,   552,   557,
     562,   566,   608,   660,   661,   663,   671,   673,   681,   683,
     685,   687,   689,   693,   701,   703,   705,   707,   709,   711,
     713,   715,   717,   719,   721,   723,   725,   729,   731,   735,
     737,   739,   741,   743,   746,   748,   750,   754,   756,   758,
     762,   764,   766,   768,   770,   772,   774,   776,   780,   781,
     782,   783,   784,   785,   788,   789,   793,   795,   804,   806,
     817,   819,   823,   825,   829,   831,   835,   837,   841,   843,
     847,   854,   861,   870,   874,   875,   879,   885,   890,   895,
     901,   902,   904,   906,   911,   912,   916,   917,   921,   924,
     929,   934,   938,   943,   949,   956,   959,   961,   963,   965,
     973,   975,   977,   979,   981,   983,   985,   987,   989,   991,
     993,   997,   999,  1001,  1003,  1006,  1008,  1010,  1012,  1014,
    1016,  1018,  1020,  1023,  1025,  1027,  1029,  1031,  1033,  1035,
    1037,  1039,  1041,  1043,  1045,  1047,  1049,  1051,  1053,  1055,
    1057,  1059,  1063,  1064,  1067,  1070,  1072,  1074,  1076,  1078,
    1080,  1082,  1084,  1086,  1088,  1090,  1093,  1097,  1100,  1103,
    1105,  1107,  1110,  1112,  1114,  1118,  1122,  1128,  1130,  1133,
    1135,  1139,  1141,  1145,  1147,  1152,  1155,  1159,  1163,  1165,
    1167,  1170,  1172,  1173,  1174,  1175,  1177,  1179,  1181,  1183,
    1190,  1192,  1194,  1196,  1198,  1202,  1207,  1213,  1215,  1219,
    1223,  1228,  1239,  1241,  1242,  1244,  1247,  1249,  1251,  1253,
    1259,  1270,  1273,  1274,  1275,  1277,  1281,  1283,  1287,  1297,
    1300,  1309,  1326,  1347,  1352,  1354,  1356,  1358,  1362,  1364,
    1368,  1372,  1376,  1380,  1382,  1386,  1390,  1394,  1398,  1400,
    1402,  1404,  1407,  1412,  1417,  1422,  1430,  1432,  1436,  1438,
    1440,  1443,  1446,  1454,  1464,  1468,  1472,  1476,  1482,  1487,
    1496,  1497,  1500,  1502,  1504,  1506,  1508,  1510,  1512,  1514,
    1516,  1518,  1522,  1526,  1533,  1539,  1543,  1549,  1557,  1563,
    1568,  1571,  1577,  1585,  1587,  1589,  1591,  1593,  1607,  1609,
    1619,  1623,  1625,  1629,  1631,  1635,  1676,  1680,  1686,  1692,
    1694,  1697,  1699,  1701,  1710,  1712,  1714,  1717,  1728,  1730,
    1732,  1735,  1746,  1748,  1750,  1755,  1758,  1759,  1762,  1764,
    1766,  1770,  1771,  1774,  1780,  1782,  1787,  1790,  1791,  1794,
    1797,  1800,  1805,  1806,  1809,  1814,  1815,  1818,  1822,  1823,
    1826,  1832,  1833,  1836,  1840,  1844,  1846,  1850,  1854,  1856,
    1860,  1862,  1865,  1867,  1870,  1874,  1877,  1879,  1883,  1889,
    1891,  1893,  1895,  1899,  1901,  1903,  1904,  1905,  1911,  1913,
    1915,  1918,  1922,  1929,  1931,  1933,  1944,  1955,  1958,  1963,
    1965,  1978,  1980,  1992,  1994,  1997,  2001,  2008,  2011,  2017,
    2023,  2025,  2027,  2029,  2031,  2033,  2035,  2037,  2046,  2049,
    2052,  2055,  2060,  2062,  2064,  2066,  2068,  2070,  2074,  2076,
    2080,  2082,  2086,  2087,  2090,  2092,  2094,  2098,  2099,  2102,
    2104,  2106,  2110,  2111,  2114,  2116,  2118,  2122,  2123,  2126,
    2128,  2130,  2132,  2134,  2138,  2139,  2142,  2144,  2146,  2150,
    2151,  2154,  2156,  2158,  2162,  2163,  2166,  2168,  2170,  2174,
    2175,  2178,  2180,  2182,  2186,  2187,  2190,  2198,  2206,  2216,
    2217,  2218,  2224,  2227,  2231,  2235,  2237,  2243,  2246,  2249,
    2253,  2258,  2266,  2285,  2286,  2289,  2291,  2293,  2297,  2299,
    2301,  2305,  2313,  2323,  2325,  2330,  2332,  2336,  2337,  2340,
    2347,  2370,  2378,  2402,  2408,  2415,  2423,  2436,  2452,  2469,
    2486,  2494,  2502,  2520,  2527,  2534,  2545,  2561,  2568,  2579,
    2587,  2599,  2614,  2627,  2644,  2656,  2669,  2683,  2699,  2716,
    2734,  2747,  2761,  2776,  2793,  2811,  2830,  2838,  2847,  2859,
    2872,  2881,  2892,  2905,  2920,  2930,  2943,  2957,  2974,  2981,
    2991,  2999,  3015,  3024,  3025,  3028,  3034,  3040,  3047,  3055,
    3066,  3068,  3070,  3074,  3076,  3078,  3080,  3082,  3086,  3092,
    3094,  3102,  3104,  3110,  3117,  3124,  3131,  3139,  3146,  3154,
    3162,  3170,  3181,  3183,  3189,  3197,  3200,  3203,  3209,  3211,
    3217,  3227,  3229,  3231,  3233,  3235,  3237,  3243,  3249,  3256,
    3263,  3270,  3277,  3285,  3294,  3303,  3314,  3315,  3317,  3319,
    3321,  3323,  3325,  3327,  3329,  3331,  3333,  3335,  3337,  3339,
    3341,  3343,  3345,  3348,  3350,  3352,  3354,  3356,  3358,  3360,
    3362,  3364,  3366,  3369,  3371,  3377,  3381,  3383,  3385,  3387,
    3391,  3396,  3401,  3403,  3408,  3410,  3412,  3414,  3418,  3422,
    3426,  3428,  3432,  3433,  3437,  3438,  3442,  3447,  3452,  3453,
    3457,  3464,  3471,  3478,  3487,  3488,  3495,  3497,  3500,  3503,
    3508,  3510,  3514,  3519,  3525,  3527,  3531,  3537,  3539,  3543,
    3549,  3554,  3559,  3564,  3571,  3574,  3577,  3579,  3581,  3585,
    3586,  3587,  3590,  3594,  3598,  3600,  3604,  3605,  3608,  3612,
    3613,  3616,  3622,  3624,  3626,  3630,  3636,  3641,  3648,  3652,
    3658,  3664,  3669,  3676,  3681,  3690,  3695,  3709,  3713,  3715,
    3717,  3722,  3725,  3733,  3735,  3737,  3740,  3742,  3744,  3746,
    3748,  3750,  3752,  3754,  3758,  3760,  3762,  3764,  3766,  3768,
    3770,  3772,  3774,  3776,  3778,  3780,  3784,  3786,  3790,  3792,
    3796,  3798,  3800,  3801,  3805,  3811,  3813,  3816,  3820,  3825,
    3827,  3831,  3833,  3837,  3841,  3845,  3849,  3854,  3859,  3864,
    3868,  3873,  3879,  3886,  3894,  3903,  3906,  3910,  3916,  3922,
    3927,  3932,  3937,  3943,  3953,  3966,  3980,  3996,  4012,  4015,
    4019,  4023,  4025,  4027,  4029,  4031,  4033,  4035,  4037,  4039,
    4041,  4043,  4050,  4052,  4059,  4061,  4063,  4065,  4067,  4069,
    4071,  4073,  4075,  4077,  4079,  4081,  4084,  4087,  4089,  4091,
    4093,  4100,  4106,  4112,  4114,  4117,  4119,  4122,  4125,  4128,
    4131,  4133,  4152,  4155,  4157,  4160,  4162,  4165,  4167,  4169,
    4171,  4173,  4175,  4177,  4179,  4182,  4190,  4209,  4225,  4229,
    4234,  4239,  4245,  4249,  4251,  4254,  4256,  4258,  4260,  4262,
    4264,  4266,  4269,  4271,  4273,  4275,  4277,  4301,  4303,  4305,
    4307,  4310,  4313,  4315,  4317,  4319,  4321,  4323,  4325,  4327,
    4329,  4332,  4334,  4336,  4338,  4340,  4343,  4346,  4349,  4352,
    4356,  4360,  4362,  4364,  4367,  4369,  4371,  4375,  4378,  4389,
    4395,  4403,  4410,  4417,  4425,  4434,  4444,  4450,  4456,  4462,
    4468,  4471,  4476,  4493,  4498,  4504,  4511,  4516,  4521,  4526,
    4527,  4528,  4531,  4532,  4535,  4539,  4542,  4547,  4548,  4551,
    4555,  4558,  4569,  4572,  4583,  4589,  4591,  4594,  4606,  4608,
    4611,  4616,  4618,  4620,  4622,  4624,  4627,  4632,  4635,  4640,
    4643,  4646,  4648,  4650,  4652,  4655,  4661,  4667,  4674,  4676,
    4679,  4683,  4689,  4691,  4694,  4699,  4705,  4714,  4718,  4719,
    4721,  4731,  4732,  4736,  4739,  4742,  4744,  4747,  4749,  4753,
    4754,  4759,  4764,  4769,  4779,  4790,  4801,  4806,  4811,  4816,
    4826,  4837,  4848,  4852,  4856,  4860,  4871,  4876,  4881,  4886,
    4896,  4907,  4918,  4923,  4928,  4933,  4939,  4945,  4951,  4956,
    4961,  4966,  4972,  4978,  4984,  4989,  4994,  4999,  5005,  5011,
    5017,  5022,  5027,  5032,  5037,  5042,  5047,  5052,  5057,  5062,
    5067,  5072,  5077,  5082,  5087,  5093,  5099,  5104,  5109,  5115,
    5121,  5126,  5132,  5137,  5143,  5148,  5153,  5158,  5168,  5179,
    5190,  5195,  5200,  5205,  5215,  5226,  5237,  5243,  5252,  5261,
    5271,  5281,  5290,  5299,  5308,  5317,  5326,  5335,  5344,  5353,
    5362,  5371,  5380,  5389,  5398,  5407,  5417,  5427,  5436,  5445,
    5455,  5465,  5474,  5483,  5493,  5503,  5508,  5514,  5523,  5533,
    5542,  5552,  5557,  5563,  5572,  5582,  5591,  5601,  5608,  5615,
    5622,  5629,  5636,  5643,  5650,  5657,  5664,  5671,  5678,  5685,
    5692,  5699,  5706,  5713,  5720,  5725,  5730,  5735,  5741,  5747,
    5756,  5765,  5771,  5777,  5786,  5795,  5800,  5805,  5810,  5815,
    5820,  5828,  5837,  5856,  5860,  5866,  5867,  5869,  5876,  5877,
    5888,  5889,  5900,  5904,  5907,  5914,  5918,  5922,  5933,  5935,
    5937,  5939,  5941,  5943,  5945,  5947,  5949,  5952,  5954,  5956,
    5958,  5961,  5963,  5965,  5967,  5980,  5986,  5993,  5997,  6004,
    6008,  6010,  6012,  6014,  6017,  6019,  6022,  6026,  6028,  6030,
    6033,  6037,  6041,  6045,  6049,  6053,  6057,  6058,  6059,  6060,
    6064,  6068,  6072,  6074,  6077,  6080,  6082,  6085,  6089,  6091,
    6093,  6096,  6098,  6100,  6103,  6107,  6109,  6111,  6113,  6115,
    6117,  6119,  6121,  6124,  6128,  6130,  6132,  6134,  6136,  6138,
    6140,  6142,  6144,  6146,  6148,  6151,  6155,  6158,  6162,  6167,
    6169,  6171,  6173,  6175,  6177,  6179,  6181,  6183,  6186,  6189,
    6191,  6198,  6205,  6206,  6211,  6220,  6222,  6226,  6231,  6233,
    6239,  6243,  6247,  6250,  6254,  6259,  6261,  6266,  6268,  6270,
    6272
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "ACTUATORS", "ADJOINTBASIS", "AERO",
  "AEROH", "AEROTYPE", "AGRESSIVETOLERANCES", "ALPROC", "AMAT", "ANALYSIS",
  "ARCLENGTH", "ATTRIBUTES", "ANGULAROUTTYPE", "ARUBBERMAT", "AUGMENT",
  "AUGMENTTYPE", "AUTOTOL", "AUXILIARY", "AVERAGED", "ATDARB", "ACOU",
  "ATDDNB", "ATDROB", "ARPACK", "ATDDIR", "ATDNEU", "ALLOWMECHANISMS",
  "AUXCOARSESOLVER", "ACMECNTL", "ADDEDMASS", "AEROEMBED", "ANDESCLR",
  "ANDESCQR", "ANDESBETAB", "ANDESALPHA", "ANDESBETAM", "AUGMENTED",
  "BLOCKDIAG", "BOFFSET", "BUCKLE", "BGTL", "BMPC", "BINARYINPUT",
  "BINARYOUTPUT", "BLOCKSIZE", "CHECKTOKEN", "COARSESOLVER", "COEF",
  "CFRAMES", "COLLOCATEDTYPE", "CONVECTION", "COMPOSITE", "CONDITION",
  "CONTACT", "CONTROL", "CORNER", "CORNERTYPE", "CURVE", "CCTTOL",
  "CCTSOLVER", "CRHS", "COUPLEDSCALE", "CONTACTSURFACES", "CMPC", "CNORM",
  "COMPLEXOUTTYPE", "CONSTRMAT", "CASES", "CONSTRAINEDSURFACES",
  "CSFRAMES", "CSTYPE", "CONSTANT", "CONWEP", "DAMPING", "DblConstant",
  "DEFAULTPENALTY", "DELETEELEMENTS", "DEM", "DIMASS", "DISP", "DIRECT",
  "DLAMBDA", "DP", "DYNAM", "DETER", "DECOMPOSE", "DECOMPFILE", "DMPC",
  "DEBUGCNTL", "DEBUGICNTL", "DOCLUSTERING", "DOROWCLUSTERING", "ANGLE",
  "DUALBASIS", "DUALRB", "KMEANS", "CRANDOM", "CONSTRAINTS", "MULTIPLIERS",
  "PENALTY", "ELLUMP", "EIGEN", "EFRAMES", "ELSCATTERER", "END",
  "ELHSOMMERFELD", "ENGINEERING", "ETEMP", "EXPLICIT", "EXTFOL", "EPSILON",
  "ELEMENTARYFUNCTIONTYPE", "FABMAT", "FABRICMAP", "FABRICMAT", "FACE",
  "FACOUSTICS", "FETI", "FETI2TYPE", "FETIPREC", "FFP", "FFPDIR", "FITALG",
  "FNAME", "FLAGMN", "FLUX", "FORCE", "FRONTAL", "FETIH",
  "FIELDWEIGHTLIST", "FILTEREIG", "FLUID", "FREEPLAY", "FREQSWEEP",
  "FREQSWEEP1", "FREQSWEEP2", "FREQSWEEPA", "FREQSWEEPAW", "FSGL",
  "FSINTERFACE", "FSISCALING", "FSIELEMENT", "NOLOCALFSISPLITING",
  "FSICORNER", "FFIDEBUG", "FAILSAFE", "FRAMETYPE", "GEPS",
  "GLOBALSEARCHCULL", "GLOBALTOL", "GRAVITY", "GRBM", "GTGSOLVER",
  "GLOBALCRBMTOL", "GROUP", "GROUPTYPE", "GOLDFARBTOL", "HDIRICHLET",
  "HEAT", "HFETI", "HNEUMAN", "HSOMMERFELD", "HFTT", "HELMHOLTZ", "HNBO",
  "HELMMF", "HELMSO", "HSCBO", "HWIBO", "HZEM", "HZEMFILTER", "HLMPC",
  "HERMITIAN", "HESSIAN", "IACC", "IDENTITY", "IDIS", "IDIS6",
  "ILUDROPTOL", "IntConstant", "INTERFACELUMPED", "ITEMP", "ITERTYPE",
  "IVEL", "IMESH", "INCIDENCE", "IHDIRICHLET", "IHDSWEEP", "IHNEUMANN",
  "ISOLVERTYPE", "INPC", "INFINTY", "JACOBI", "KEYLETTER", "KRYLOVTYPE",
  "KIRLOC", "LAYC", "LAYN", "LAYD", "LAYO", "LAYMAT", "LFACTOR", "LISRBM",
  "LMPC", "LOAD", "LOADCASE", "LOBPCG", "LOCALREDUCEDORDERBASES",
  "LOCALSOLVER", "LINESEARCH", "LUMPED", "KSPARAM", "KSMAX", "MASS",
  "MASSAUGMENTATION", "MATERIALS", "MATLAB", "MAXITR", "MAXELEM",
  "MAXORTHO", "MAXVEC", "MODAL", "MPCPRECNO", "MPCPRECNOID", "MPCTYPE",
  "MPCTYPEID", "MPCSCALING", "MPCELEMENT", "MPCBLOCKID", "MPCBLK_OVERLAP",
  "MFTT", "MRHS", "MPCCHECK", "MUMPSICNTL", "MUMPSCNTL", "MUMPSMINEQ",
  "MUMPSSTRIDE", "MECH", "MODDAMP", "MODEFILTER", "MOMENTTYPE", "MPROJECT",
  "MAXIMUM", "NOWARPEDVOLUME", "NDTYPE", "NEIGPA", "NEWMARK", "NewLine",
  "NEWTON", "NL", "NLMAT", "NLPREC", "NLMEMPTYPE", "NOCOARSE", "NODETOKEN",
  "NOMULTIPLEINTERACTIONS", "NONINPC", "NSBSPV", "NLTOL", "NUMCGM",
  "NONORMALSMOOTHING", "NOSECONDARY", "NOGHOSTING", "NFRAMES",
  "NORMALSMOOTHINGDISTANCE", "SENSITIVITY", "SENSITIVITYMETHOD",
  "SKIPPHYSICALFACES", "OUTPUT", "OUTPUT6", "OUTPUTFRAME",
  "OLDDYNAMICSEARCH", "QSTATIC", "QLOAD", "QUASISTATIC", "PITA",
  "PITADISP6", "PITAVEL6", "NOFORCE", "MDPITA", "GLOBALBASES",
  "LOCALBASES", "TIMEREVERSIBLE", "REMOTECOARSE", "ORTHOPROJTOL",
  "READINITSEED", "JUMPCVG", "JUMPOUTPUT", "PRECNO", "PRECONDITIONER",
  "PRELOAD", "PRESSURE", "PRINTMATLAB", "printmatlab", "PRINTNUMBER",
  "PROJ", "PIVOT", "PRECTYPE", "PRECTYPEID", "PICKANYCORNER", "PADEPIVOT",
  "PROPORTIONING", "PLOAD", "PADEPOLES", "POINTSOURCE", "PLANEWAVE",
  "PTOL", "PLANTOL", "PMAXIT", "PIECEWISE", "PARAMETERS", "RADIATION",
  "RBMFILTER", "RBMSET", "READMODE", "READSENSITIVITY", "REBUILD",
  "RESOLUTIONMETHOD", "REVERSEORDER", "REDFOL", "RENUM", "RENUMBERID",
  "REORTHO", "RESTART", "RECONS", "RECONSALG", "REBUILDCCT", "RANDOM",
  "RPROP", "RNORM", "REVERSENORMALS", "ROBTYPE", "ROTVECOUTTYPE",
  "RESCALING", "RUBDAFT", "SCALING", "SCALINGTYPE", "SDETAFT", "SENSORS",
  "SHELLFABRICMAP", "SHELLFABRICMAT", "SHELLSIMPLELOFTING", "SOLVERCNTL",
  "SOLVERHANDLE", "SOLVERTYPE", "SHIFT", "SHARPNONSHARPANGLE",
  "SPOOLESTAU", "SPOOLESSEED", "SPOOLESMAXSIZE", "SPOOLESMAXDOMAINSIZE",
  "SPOOLESMAXZEROS", "SPOOLESMSGLVL", "SPOOLESSCALE", "SPOOLESPIVOT",
  "SPOOLESRENUM", "SPARSEMAXSUP", "SPARSEDEFBLK", "SPARSERENUM", "STATS",
  "STRESSID", "SUBCYCLE", "SUBSPACE", "SURFACE", "STR_THERM_OPTION",
  "SAVEMEMCOARSE", "SPACEDIMENSION", "SCATTERER", "STAGTOL", "SCALED",
  "SWITCH", "STABLE", "SUBTYPE", "STEP", "SOWER", "SHELLTHICKNESS", "SURF",
  "SPRINGMAT", "SENSITIVITYID", "TABLE", "TANGENT", "TDENFORCE", "TEMP",
  "TIME", "TOLEIG", "TOLFETI", "TOLJAC", "TOLPCG", "TOLSEN", "TOPFILE",
  "TOPOLOGY", "TRBM", "TRBMlc", "THERMOE", "THERMOH", "RATIOTOLSEN",
  "TETT", "TOLCGM", "TURKEL", "TIEDSURFACES", "THETA", "PROJSOL", "CENTER",
  "POSELEM", "HOTSTART", "HRC", "THIRDNODE", "THERMMAT", "TDENFORC",
  "TESTULRICH", "THRU", "TRIVIAL", "NUMROMCPUS", "THICKNESSGROUPLIST",
  "USE", "USERDEFINEDISP", "USERDEFINEFORCE", "UPROJ", "UNSYMMETRIC",
  "USING", "USESCOTCH", "VERBOSE", "VERSION", "WETCORNERS", "YMTT",
  "SS1DT", "SS2DT", "YMST", "YSST", "YSSRT", "ZERO", "BINARY", "GEOMETRY",
  "DECOMPOSITION", "GLOBAL", "MATCHER", "CPUMAP", "NODALCONTACT", "MODE",
  "FRIC", "GAP", "OUTERLOOP", "EDGEWS", "WAVETYPE", "ORTHOTOL", "IMPE",
  "FREQ", "DPH", "WAVEMETHOD", "MATSPEC", "MATUSAGE", "BILINEARPLASTIC",
  "FINITESTRAINPLASTIC", "CRUSHABLEFOAM", "LINEARELASTIC",
  "STVENANTKIRCHHOFF", "TULERBUTCHER", "LINPLSTRESS", "READ", "OPTCTV",
  "ISOTROPICLINEARELASTIC", "VISCOLINEARELASTIC", "VISCOSTVENANTKIRCHHOFF",
  "NEOHOOKEAN", "VISCONEOHOOKEAN", "ISOTROPICLINEARELASTICJ2PLASTIC",
  "ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS", "HYPERELASTIC",
  "MOONEYRIVLIN", "VISCOMOONEYRIVLIN", "HENCKY", "OGDEN", "SIMOELASTIC",
  "SIMOPLASTIC", "LOGSTRAINPLASTIC", "SVKPLSTRESS", "VISCOLINPLSTRESS",
  "VISCOSVKPLSTRESS", "VISCOFABRICMAP", "SHELLVISCOFABRICMAP",
  "VISCOFABRICMAT", "SHELLVISCOFABRICMAT", "PARTITIONGAP",
  "PLANESTRESSLINEAR", "PLANESTRESSSTVENANTKIRCHHOFF",
  "PLANESTRESSNEOHOOKEAN", "PLANESTRESSMOONEYRIVLIN",
  "PLANESTRESSBILINEARPLASTIC", "PLANESTRESSFINITESTRAINPLASTIC",
  "PLANESTRESSVISCOLINEARELASTIC", "PLANESTRESSVISCOSTVENANTKIRCHHOFF",
  "PLANESTRESSVISCONEOHOOKEAN", "PLANESTRESSVISCOMOONEYRIVLIN",
  "SURFACETOPOLOGY", "MORTARTIED", "MORTARSCALING",
  "MORTARINTEGRATIONRULE", "SEARCHTOL", "STDMORTAR", "DUALMORTAR",
  "WETINTERFACE", "NSUBS", "EXITAFTERDEC", "SKIP", "ROBCSOLVE",
  "RANDOMSAMPLE", "OUTPUTMEMORY", "OUTPUTWEIGHT", "SOLVER",
  "SPNNLSSOLVERTYPE", "MAXSIZE", "CLUSTERSOLVER", "CLUSTERSOLVERTYPE",
  "WEIGHTLIST", "GMRESRESIDUAL", "SLOSH", "SLGRAV", "SLZEM", "SLZEMFILTER",
  "PDIR", "HEFSB", "HEFRS", "HEINTERFACE", "SNAPFI", "VELSNAPFI",
  "ACCSNAPFI", "DSVSNAPFI", "MUVSNAPFI", "PODROB", "ROMENERGY", "TRNVCT",
  "OFFSET", "ORTHOG", "SVDTOKEN", "CONVERSIONTOKEN", "CONVFI", "ROMRES",
  "SAMPLING", "SNAPSHOTPROJECT", "PODSIZEMAX", "REFSUBTRACT", "TOLER",
  "NORMALIZETOKEN", "FNUMBER", "SNAPWEIGHT", "ROBFI", "STAVCT", "VELVCT",
  "ACCVCT", "CONWEPCFG", "SCALEPOSCOORDS", "NODEPOSCOORDS",
  "MESHSCALEFACTOR", "PSEUDOGNAT", "PSEUDOGNATELEM", "USENMF", "USENMFC",
  "USEGREEDY", "USEPQN", "FILTERROWS", "VECTORNORM", "REBUILDFORCE",
  "REBUILDCONSTRAINT", "SAMPNODESLOT", "REDUCEDSTIFFNESS", "UDEIMBASIS",
  "FORCEROB", "CONSTRAINTROB", "DEIMINDICES", "UDEIMINDICES",
  "SVDFORCESNAP", "SVDCONSTRAINTSNAP", "USEMASSNORMALIZEDBASIS",
  "USECONSTANTMASS", "ONLINEMASSNORMALIZEBASIS", "STACKED",
  "NUMTHICKNESSGROUP", "STRESSNODELIST", "DISPNODELIST", "RELAXATIONSEN",
  "QRFACTORIZATION", "QMATRIX", "RMATRIX", "XMATRIX", "EIGENVALUE",
  "NPMAX", "BSSPLH", "PGSPLH", "LIB", "$accept", "FinalizedData", "All",
  "Component", "Noninpc", "Inpc", "Group", "Random", "Impe",
  "ImpeDampInfo", "PadePivotInfo", "PadePolesInfo", "FreqSweep",
  "ReconsInfo", "BinarySpec", "AnalysisInfo", "Decompose", "WeightList",
  "FieldWeightList", "MFTTInfo", "HFTTInfo", "LoadCase", "Composites",
  "Cframes", "CoefInfo", "CoefList", "LaycInfo", "LaynInfo", "LaydInfo",
  "LayoInfo", "LayData", "LayoData", "LayMat", "LayMatData", "DiscrMasses",
  "Gravity", "Restart", "LoadCInfo", "UseCInfo", "SensorLocations",
  "ActuatorLocations", "UsdfLocations", "UsddLocations", "Output",
  "OutInfo", "DynInfo", "PrintMat", "DeleteElements", "DeleteElementsList",
  "SloshInfo", "MassInfo", "CondInfo", "TopInfo", "ModalInfo", "DynamInfo",
  "Conwep", "ConwepData", "TimeIntegration", "NewmarkSecondOrder",
  "NewmarkFirstOrder", "QstaticInfo", "QMechInfo", "QHeatInfo", "AeroInfo",
  "AeroEmbeddedSurfaceInfo", "AeroHeatInfo", "ThermohInfo", "ThermoeInfo",
  "ModeInfo", "HzemInfo", "SlzemInfo", "RbmTolerance", "ToleranceInfo",
  "ModeFilterInfo", "RbmFilterInfo", "RbmList", "HzemFilterInfo",
  "SlzemFilterInfo", "TimeInfo", "ParallelInTimeInfo",
  "ParallelInTimeOptions", "ParallelInTimeKeyWord", "DampInfo",
  "ComplexDirichletBC", "IComplexDirichletBC", "IComplexDirichletBCSweep",
  "DirectionVector", "IComplexNeumannBC", "DirichletBC",
  "ConstrainedSurfaces", "HEVDirichletBC", "HEVDBCDataList", "HEVDBC_Data",
  "HEVFRSBC", "HEVFRSBCList", "HEVFRSBCElem", "TempDirichletBC",
  "TempNeumanBC", "TempConvection", "TempRadiation", "HelmHoltzBC",
  "SommerfeldBCDataList", "SommerfeldBC_Data", "EleHelmHoltzBC",
  "SommerElement", "SommNodeNums", "Scatterer", "ScattererEleList",
  "ScattererEle", "EleScatterer", "ScatterElement", "NeumScatterer",
  "NeumElement", "WetScatterer", "WetInterfaceElement", "HelmScatterer",
  "HelmScattererElement", "PBC_Data", "PBCDataList", "AtdDirScatterer",
  "AtdNeuScatterer", "AtdArbScatterer", "AtdNeumScatterer",
  "AtdRobinScatterer", "FarFieldPattern", "FarFieldPatternDirs",
  "ReadModeInfo", "Mode", "IDisp", "IDisp6", "IDisp6Pita", "IVel6Pita",
  "IVel", "ITemp", "ETemp", "NeumanBC", "ModalNeumanBC", "BCDataList",
  "ModalValList", "TBCDataList", "YMTTable", "YMTTList", "TETTable",
  "TETTList", "SS1DTable", "SS1DTList", "SS2DTable", "SS2DTList",
  "YSSTable", "YSSTList", "YSSRTable", "YSSRTList", "YMSTable", "YMSTList",
  "SDETAFTable", "SDETAFList", "RUBDAFTable", "RUBDAFList", "LMPConstrain",
  "ModalLMPConstrain", "MPCList", "MPCHeader", "MPCLine",
  "ComplexLMPConstrain", "ComplexMPCList", "ComplexMPCHeader",
  "ComplexMPCLine", "ComplexNeumanBC", "ComplexBCDataList", "Materials",
  "MatData", "FreeplayProps", "ElemSet", "FaceSet", "MortarCondition",
  "WetInterface", "TiedSurfaces", "FSInterface", "HEVibInfo",
  "HEVInterfaceElement", "HEInterface", "ContactSurfaces",
  "ContactSurfacesInfo", "Parameters", "NodeSet", "Node", "Element",
  "NodeNums", "BC_Data", "ModalVal", "TBC_Data", "ComplexBC_Data",
  "ConstrainedSurfaceFrameDList", "FrameDList", "Frame", "NodalFrameDList",
  "NodalFrame", "BoffsetList", "Attributes", "Ellump",
  "LocalReducedOrderBases", "LocalBasesAuxi", "LocalBasesCent",
  "ReducedStiffness", "UDeimBasis", "SampNodeSlot", "Pressure", "Lumped",
  "MassAugmentation", "Preload", "Sensitivity", "DispNode", "StressNode",
  "ThicknessGroup", "Statics", "CasesList", "SolverMethod", "Solvercntl",
  "Solver", "OldHelmInfo", "FAcousticData", "Constraints",
  "ConstraintOptionsData", "HelmInfo", "IncidenceList", "IncidenceVector",
  "KirchhoffLocations", "FFPDirList", "FFPDirVector", "HelmMFInfo",
  "FAcousticDataMF", "HelmSOInfo", "FAcousticDataSO", "DEMInfo", "AlProc",
  "NLInfo", "NewtonInfo", "OrthoInfo", "Control", "NodalContact",
  "MatSpec", "MatUsage", "FloatList", "StringList", "Renumbering",
  "SvdToken", "SvdOption", "DeimIndices", "UDeimIndices", "Sampling",
  "SnapshotProject", "SamplingOption", "ConwepConfig", "MeshScaleFactor",
  "ScalePosCoords", "NodePosCoords", "ConversionToken", "ConversionOption",
  "Integer", "Float", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
     395,   396,   397,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,   412,   413,   414,
     415,   416,   417,   418,   419,   420,   421,   422,   423,   424,
     425,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     435,   436,   437,   438,   439,   440,   441,   442,   443,   444,
     445,   446,   447,   448,   449,   450,   451,   452,   453,   454,
     455,   456,   457,   458,   459,   460,   461,   462,   463,   464,
     465,   466,   467,   468,   469,   470,   471,   472,   473,   474,
     475,   476,   477,   478,   479,   480,   481,   482,   483,   484,
     485,   486,   487,   488,   489,   490,   491,   492,   493,   494,
     495,   496,   497,   498,   499,   500,   501,   502,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,   517,   518,   519,   520,   521,   522,   523,   524,
     525,   526,   527,   528,   529,   530,   531,   532,   533,   534,
     535,   536,   537,   538,   539,   540,   541,   542,   543,   544,
     545,   546,   547,   548,   549,   550,   551,   552,   553,   554,
     555,   556,   557,   558,   559,   560,   561,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,   578,   579,   580,   581,   582,   583,   584,
     585,   586,   587,   588,   589,   590,   591,   592,   593,   594,
     595,   596,   597,   598,   599,   600,   601,   602,   603,   604,
     605,   606,   607,   608,   609,   610,   611,   612,   613,   614,
     615,   616,   617,   618,   619,   620,   621,   622,   623,   624,
     625,   626,   627,   628,   629,   630,   631,   632,   633,   634,
     635,   636,   637,   638,   639,   640,   641,   642,   643,   644,
     645,   646,   647,   648,   649,   650,   651,   652,   653,   654,
     655,   656,   657,   658,   659,   660,   661,   662,   663,   664,
     665,   666,   667,   668,   669,   670,   671,   672,   673,   674,
     675,   676,   677,   678,   679,   680,   681,   682,   683,   684,
     685,   686,   687,   688,   689,   690,   691,   692,   693,   694,
     695,   696,   697,   698,   699,   700,   701,   702,   703,   704,
     705,   706,   707,   708,   709,   710,   711,   712,   713,   714,
     715,   716,   717,   718,   719,   720,   721,   722,   723,   724,
     725,   726,   727,   728,   729,   730,   731,   732,   733,   734,
     735,   736,   737,   738,   739,   740,   741,   742,   743,   744,
     745,   746,   747,   748,   749,   750,   751,   752,   753,   754,
     755,   756,   757,   758,   759,   760,   761,   762,   763,   764,
     765,   766,   767,   768,   769,   770,   771,   772,   773,   774,
     775,   776,   777,   778,   779,   780,   781,   782,   783,   784,
     785,   786,   787,   788,   789,   790,   791,   792,   793,   794,
     795,   796,   797,   798,   799,   800,   801,   802,   803,   804,
     805,   806,   807,   808,   809,   810,   811,   812,   813,   814,
     815,   816,   817,   818,   819,   820,   821,   822,   823,   824,
     825,   826,   827,   828,   829,   830,   831,   832,   833,   834,
     835,   836,   837,   838,   839,   840,   841,   842,   843,   844,
     845
};
# endif

#define YYPACT_NINF -2828

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-2828)))

#define YYTABLE_NINF -1426

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    7995,  -155,    31,  -113,   797,   890,  1128,  1128,  1128,   896,
    -109,   947,   797,  -101,   -87,   797,    99,  1062,   116,   188,
     209,  -107,   230,   235,   273,   277,   303,   319,   345,   410,
     429,   546,  1069,  1070,   432,   438,   462,   467,   476,   482,
     797,   797,  1071,  1100,   504,  -258,   514,   589,   598,   970,
     602,   608,  1112,   662,  1295,   684,   698,   722,   733,   756,
     780,   786,   798,   818,  -185,  1025,   829,   834,   797,   848,
    1322,   856,   864,   797,   886,   891,  1046,    44,  1395,   901,
     906,  1410,   284,     0,   924,  1429,  1437,   797,   928,  1053,
     941,   797,   945,   948,   745,   998,   963,   968,   797,   797,
     972,   974,  1446,   474,   797,   797,   977,  1466,  1606,   324,
     979,   981,   986,   989,   992,   999,  1017,   797,  1128,  1022,
    1621,  1026,  1029,  1128,  1128,  1033,  1034,  1036,  1041,  1045,
    1047,  1051,  1075,   518,  1082,  1084,  1087,  1088,  1090,  1094,
    1102,  1111,  1113,    -7,  1623,  1624,  1117,  1120,   797,   797,
     797,  1131,  1133,  1128,  1135,  1137,  1140,  1141,  1145,  1147,
    1201,  1216,  1218,  1219,  1128,  1222,  1223,  1230,  1231,  1234,
     404,  1486,  7425, -2828, -2828, -2828,   843,   797,   380,    40,
   -2828,     3,   797,    32,  1128,  1128,   797,   736,   797,   797,
     976,  1128,   767, -2828, -2828, -2828, -2828, -2828,   200,   151,
    1270, -2828,  1128, -2828, -2828, -2828, -2828,    58, -2828,   885,
       1, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
    1287, -2828, -2828, -2828, -2828, -2828,   797, -2828,   394,   797,
   -2828, -2828,   425,   446,   484,   797, -2828,   797, -2828,   797,
     797,   797,   797, -2828, -2828,   797,   797,   797, -2828, -2828,
      41,  1273,   797,   797,   797,  1274,  1281, -2828,   527, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
    1128, -2828, -2828,   797,   797,   797, -2828, -2828,   797,   797,
     797,   797,   -68,   138,   797,   797,   797,   797,   797,   797,
      28,    16,  1128,   898,   797,   536, -2828, -2828,   797,   350,
     720, -2828, -2828, -2828,    54,  1128, -2828, -2828, -2828,    10,
   -2828, -2828,   797,  -142,   797, -2828,   285,  5599,  5599,  5166,
    5599, -2828,  -385,   797,  1257,  1500,  1502, -2828, -2828,  1261,
   -2828,  1262, -2828, -2828, -2828, -2828,  1263,  1267,  1128,  1658,
   -2828,  1128,   797,   797,  1268,  1276, -2828, -2828,  1086, -2828,
   -2828,  1277, -2828,  1128,  1398, -2828,   797, -2828, -2828,  1128,
    1128, -2828, -2828, -2828, -2828, -2828,  1128,  1128,  1426,  1128,
    1011,   288, -2828,  1280, -2828,  1291, -2828,   797,   797,   797,
   -2828,  1128,  1660,  1292, -2828,  1293,  1320,  1296, -2828,  1301,
   -2828, -2828, -2828,  1128,  1128, -2828,   797,   797,  1304,   797,
   -2828,  1306, -2828,   797,  1128,  1128,   797,   797, -2828, -2828,
     797,   797,  1309, -2828,  1314,   797,   797,  1316,  1128,   797,
    1330,  1128,   797,  1331,  1128, -2828,  1124,  1128,  1338, -2828,
     -33, -2828, -2828, -2828,  1346,  1348, -2828,  1349,   797, -2828,
    1360, -2828,  1377,  1387, -2828,   797,  1128,   797,  1389, -2828,
   -2828,  1397,  1676, -2828,  1400,  1402,  1677, -2828,  1403, -2828,
     797,  1407,  1413,   797, -2828, -2828,  1415,  -209,  1417,  1420,
   -2828, -2828,  1425, -2828,  1700,  1703,   797,  1246, -2828, -2828,
   -2828,  1294,  1617,   797,  1430,  1435, -2828, -2828,  1128,   797,
   -2828,  1436,  1438, -2828,   797,  1128, -2828, -2828,  1630, -2828,
   -2828,  1443, -2828,   797,  1634,  1641,  1318,  1648,  1655,  1661,
   -2828, -2828,   797, -2828,  1455, -2828,  1468, -2828, -2828,   234,
     797,  1710, -2828, -2828,  1469, -2828, -2828,   797,   797,   797,
   -2828, -2828, -2828, -2828, -2828,  1128, -2828, -2828, -2828, -2828,
   -2828,  1473, -2828, -2828, -2828,   538,  1396,  1128,  1128,  1128,
    1128,  1128,  1128,  1128,  1091,  1404,  1128, -2828, -2828, -2828,
    1128, -2828,  1353,  1366,  1615,  1620,  1628,  1629,  1633,  1483,
    1493,  1638,  1497,  1717,  1506,   797,  1516,  1518,  1520,  1528,
    1128,   797,   797,   797,   797,  1128,  1128,   938,   797,   797,
     797,   797,   797, -2828,   797,   797,   797,   797, -2828,   -11,
   -2828,  1128,  1547,   797,  1128,   675,   797, -2828,   -43,   744,
     751,  1673,  1686,  1688,  1691,   224,   797,  1738,   797,  1128,
    1096,  1128,  1128,  1447,  1837,  1128,  1575,  1579,  1582,  1462,
     -46,  1471,  1128,  1477,   439, -2828, -2828, -2828,  1128,  1854,
    1128, -2828, -2828, -2828,  1583,   797,  1726,  1896,   797,  1899,
   -2828,  1128,   797, -2828,   -27,  1782,   797,   -26,   797,  1128,
     797,  1128,  1128, -2828,   797, -2828,   797, -2828,   797, -2828,
     797, -2828,   797, -2828, -2828, -2828,  1737, -2828,   -80,  1909,
    1128,  1128,  1128,  1911,  1610,   797, -2828,    76,  1631, -2828,
     112, -2828,   797,   797,   797,   797, -2828,   797,   797, -2828,
      -4,   797,  1639,  1643,   797,  1128,  1128,  1128,  1128,  1128,
    1644,   797,  1651,  1656,   797,  1662,  1663,  1664,  1667,  1128,
    1674,  1675,   797,  1682,  1128,  1692,  1128,   797, -2828,  1128,
   -2828, -2828, -2828,  1110,   797,   -78,  1699,  1704,  1128,   797,
    1744, -2828, -2828,  1709,  1941,   797,  1711,   797,   797,  2052,
      85,  1128,  1128,  1542,  1757,  1128,  1128,   797,   797,  1713,
    1128,  1739,  1128,   797,  1919,  1603,   -60,  1923,  1943,   797,
    1099,   797,   557,   797, -2828,  6462, -2828, -2828, -2828,  1128,
    1740,  1128,   797,   797,  1148,  1128,    68,   797,  1783,  1785,
    1128,   797,   797,  1788,  1944, -2828,   797,  1917,  1848,   797,
     797,   797,   797,  1672,   797,   797,  1687,  1532, -2828, -2828,
   -2828, -2828, -2828,  1128,   797,  1929,   -53,  1693, -2828,  1799,
     797,   797, -2828,   797,  1811,  1128,  1938,  1695,  1128,   797,
    1721,  1722,  1746,  1747,  1755,  1758,  1765,  1769,   797,   797,
    1553,  1128,  1947,  1951,   797,   797,  1128,  1772,  1774,   797,
    1778,  1984,  1793,  1803,  1987,  2046,  1822,  1825,  1827,  1831,
    1832,   797,   797,   797,  1817,  1819,  1921,  1934,  1957,  1128,
    1423,   343,  1961,  2091,   797,  1972,   797, -2828,   108, -2828,
    1149,  1128, -2828, -2828, -2828, -2828,  1128, -2828,  1985,  1971,
   -2828,   797,  1128,   797,   797, -2828, -2828,  1153, -2828,   797,
    2006,  2007,  2008,  1128,  1980, -2828,  1128,  1128, -2828,   650,
     797, -2828, -2828, -2828, -2828, -2828,   797, -2828,  1128, -2828,
    2009, -2828,  2010,  1128, -2828,  2042,  2038, -2828,  1157,  1128,
     797, -2828,   797,   797,   797,   797, -2828,   797, -2828, -2828,
   -2828,  2014, -2828,  2017, -2828, -2828,   797,   797,   881,   797,
   -2828, -2828,   797,   797,  1128,  1128,  1128, -2828,  1128,  2018,
   -2828,  1128,   797,   797,   797,  1210,  1128, -2828,  1902, -2828,
    1907, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,  1128,
     797, -2828, -2828, -2828,  2034, -2828, -2828, -2828,  2048, -2828,
     797, -2828, -2828,   797, -2828, -2828,  2049,  1128,  1128, -2828,
   -2828,  2094, -2828,  2095, -2828,  2100,   797,   157,   797,  1186,
     797,   488, -2828,  1128,   797, -2828,   797, -2828, -2828, -2828,
    2102,   797,  1200, -2828,   797,   797,  1838,   797,  1880,   797,
     198,   797,  2019,   797,  2288,   797,  2297,  2103, -2828, -2828,
   -2828,  2105,  1128,  -198, -2828,  2108, -2828,   797, -2828,  1211,
   -2828,   797, -2828,   797,  1128, -2828,   797,   797,  1128,  1128,
    1288,  1128,  1128,  1128,  1128,  2113, -2828,  1128,   797,  2114,
     797,   384,   428,  2116,  2117,  2118,  2119,  2121, -2828, -2828,
    2123, -2828, -2828,  2126, -2828,  2129, -2828, -2828, -2828, -2828,
    2130,  2132,  2133,  2134,  2135,  2136,  2138,   797,   797,  2139,
     448,  2140,  2141,  2142,  2143, -2828,  1128, -2828, -2828,   797,
   -2828,   797,  1128,  1128,  2276,   898,  1128,   430,  2154, -2828,
     797,   797,   911,   797,   797,   797,   797, -2828, -2828, -2828,
   -2828, -2828, -2828,   797, -2828,   797, -2828,   797, -2828, -2828,
    2032, -2828, -2828, -2828,  1128,  2157, -2828, -2828, -2828,  2158,
    1128, -2828,  1128,  1128,  2159, -2828,  2181, -2828, -2828,  2191,
     797, -2828, -2828,  2193,  2194,  1289,  2045,  1128,  2210,   797,
   -2828, -2828,  1128, -2828,  2211,  1128, -2828,  2213,  2218, -2828,
   -2828, -2828, -2828, -2828,  1128, -2828,   797,   797,  1128,   797,
    2221,  1128,  2222,  1128,  1128,  1128,   797,   797,   797,   797,
     797,   797,  2349,  2354,   797,  2235,  1128,  1128,  1128,   797,
    2238,   797, -2828,   797, -2828,  1128,  1128,    24,   797,  1128,
    1128,  1128,  1128,   797,   797,   797,   797,   797,   797, -2828,
     797,  -133,   797, -2828, -2828,  2239,  2240,  2243,  2247,  2249,
    2250, -2828,  2251, -2828, -2828,  2252, -2828, -2828, -2828, -2828,
    2255, -2828, -2828,  2258, -2828,  2260, -2828,  2287,  2290,  1297,
    1128,  1128,  1128,  1128,   -17, -2828, -2828,  2293,   797,  2294,
   -2828,  2295, -2828,  2101,   123,   797,  1749,   238,   797,  1375,
    2299,  2300,  2302,  2305,  2308,  2317,  2104, -2828,  2125, -2828,
     797,  2318, -2828,  2321,  2127, -2828, -2828,  2327,  2331,  2335,
   -2828,  2336,  2342,   922, -2828,  1128, -2828,  1380,  2346, -2828,
    1128,  2348, -2828,  2368,  2375,  2376,  1381,     4,   488,  2153,
     488,    87,  1128,   488,  2231,  2254,  2380,   449,  2381,  2383,
     797,  1128,  1128,  1128,  1128,  2387,  1126,   488,   797,   797,
     797,  1127,  1151,   523,  2388,   797,   797,   797,   797,   797,
    2212,  2391,   797,   993,   437,   797,   797,   712,   797,  1128,
    1128,  1128,   797,  2269,   552,  1128,   797,   797,   797,  1128,
     797,   797,   797,   797,   797,   797,  2393,   797,  1128,   797,
    1128,  1128,  1128,  1128,   797,  2314,   797,  2400,  2479,  1128,
     276,  1128, -2828,  1128,  1128,  2325, -2828,  2416,  2418, -2828,
    2337, -2828,  2419, -2828, -2828,  1401,  2338,  2420, -2828, -2828,
    2421,  1128,  2546,  1128,  1128, -2828,  1128,  1128,  1128,  1128,
    1128,  1128,  1128,  1128,  1128,  1128,  1128,  1128,  1128,  1128,
    1128,  1128,  1128,  1128,  1128,  1128,  1128,  1128,  1128,  1128,
    1128,  1128,  1128,  1128,  1128,  1128,  1128,  1128,  1128,  1128,
    1128,  1128,  1128,  1128,  1128,  1128,  2364, -2828,  2578,  2581,
   -2828, -2828,   797, -2828,  1128,  2551,  2551,  2551,  2551,  2551,
   -2828, -2828, -2828, -2828, -2828, -2828,  2307,  2551, -2828,   898,
    1128,  1128, -2828,  2428,   797, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
     797,  2568, -2828,   797, -2828, -2828, -2828,   898, -2828,  2569,
   -2828, -2828,   797,   797, -2828, -2828, -2828, -2828, -2828, -2828,
     797,   797, -2828, -2828, -2828, -2828, -2828,  2447,  1128,   461,
   -2828,   797, -2828, -2828, -2828,   225,   797, -2828,  1408,  1128,
    2476, -2828, -2828,  2477, -2828,  2494,   797, -2828,  2498,  2501,
     797, -2828, -2828,  1128, -2828,  1128, -2828, -2828,  1128, -2828,
    2509, -2828, -2828,  1128, -2828,  1128,   797,  2511,  2409, -2828,
    2411,  2512, -2828,  1128,   797, -2828,   898, -2828, -2828,   797,
   -2828,   832, -2828,   797,  1128, -2828,  2513,  1128, -2828,  1128,
    1409,  1128,  1128, -2828,  1128,  2517,   797, -2828,   905, -2828,
     797, -2828,    -1,  2529,  -173,  2532,  2534,  2536, -2828, -2828,
    2412,   797, -2828,  1128,  1128, -2828, -2828,  2464,  2540,   797,
    1128,  2544,   797,  1128,  6462,  2549, -2828,   898, -2828,  2550,
     797,  1128,  2553,   797,  1128,  2554,   797,  1128,    12,   797,
   -2828,  2555,   797,  1128,  2557,   797,  1128,  2560,   797,  1128,
   -2828, -2828,  -151, -2828,  1128,  -188,  -182, -2828, -2828, -2828,
    2565, -2828,   797,  2577,   797,  2413,  1128,  2593, -2828,   797,
     797,   797,   797,   797, -2828,  2594,  2415, -2828,  2596,  2597,
   -2828,  2599, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,  2601,  2608,
   -2828,   797,  2610,  1128, -2828, -2828, -2828, -2828,  1128,  1128,
    2612,  1128,  1128,  2613,  1128,  2619,  2621,  2622, -2828, -2828,
      78,  2624, -2828,   559,   758,  2675,   794,  2692, -2828, -2828,
   -2828, -2828, -2828, -2828, -2828,  1128,  2625,   647, -2828, -2828,
   -2828,   797, -2828,  1128, -2828, -2828, -2828,  1128, -2828,  1128,
    1128, -2828, -2828,  1128, -2828,   797, -2828, -2828,  1128,   321,
     797,  2627,    19, -2828,  2628, -2828,  1128,  1128,  1128,  2500,
   -2828,  2510,  2589,  2590,  2614,  2456,   797,   797,   797,   797,
    1128,  1128,  1128,   797,   797,   797,   406,  1128,  1412,  1128,
    1128, -2828,  1128,   797,   -66,  1128,  1128,  1128,   -59,  1128,
    1128,  2616, -2828,  2620,   669,  1418,  2631,  1419,  2636, -2828,
     797,   797, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828,  1422,  1128,  1128,  1128,
     -39,  2637, -2828,   178,   290, -2828,  2702, -2828, -2828, -2828,
    2639,   797,   322,  1128,   797,   323, -2828,  2640,  1128, -2828,
    1128, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828,   797,   797, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828,  2642, -2828,  2634,  2635, -2828,  2716, -2828,  2643,
   -2828, -2828, -2828, -2828, -2828,  1439,   390,  2644,  -137, -2828,
    2646, -2828,  2647,  2717,  2650, -2828,  2651,  2653, -2828,  2677,
    2682, -2828, -2828,  2685,  2689,  2715,  2719,  2721, -2828,  2726,
    2727, -2828,  1128,  2729,  2731,   623,  1000,  2732,  2733,  2742,
    2746, -2828,  2748,   797,  1128,  2750,  2751, -2828,  2757, -2828,
    2758,  2772,  2774,  2775,  2779, -2828,  2780,  2781,  2782,  2788,
    2789,  2790,  1128,  1128,   797,  2799,  2800,  2803,  2808,  2823,
    2825,  2828,  2829,  2842,  2845,  2857,  2863,  2874,  2875, -2828,
    2880,  1442,  2884,  1452,  2888,  2889,  2891,  2892, -2828,  2900,
    2902, -2828,   490,  1460, -2828,  2903,  2904,  2723,  1128, -2828,
    2911, -2828, -2828, -2828,  1484, -2828, -2828,  1501, -2828,  2740,
   -2828, -2828,  1128,  2914,   797,   797,  1530,   797,   797,  1128,
    1128,  1128,   -30,    14,  1128,  1128,  1128,  1128,  1128,  1128,
    1128,  1128,  1128,  1128,  1128,    47,  1128,  1128,  1128,  1128,
    1128,  1128,  1128,   797,   797,   797,   797,  1128,  1128,  1128,
    1128,  1128,  1128,  1128,  1128,  1128,  1128, -2828,  2915, -2828,
   -2828, -2828,  2466, -2828, -2828,   797, -2828,  1128,   797, -2828,
   -2828, -2828,  3040,   797,   797, -2828,  3048,   797,   797, -2828,
   -2828, -2828,  1128, -2828, -2828,   797,   471, -2828,  1543,  2929,
   -2828, -2828, -2828, -2828, -2828,  2930,  1128,  1128, -2828, -2828,
   -2828,  1128,   797,   797,   797,  2932, -2828,  2944,  1128,  1558,
     898,  2949, -2828,  1128,  1128, -2828, -2828,  1128, -2828,  2955,
    2957,  2968,  2973, -2828, -2828,  1128, -2828,   797,   731, -2828,
   -2828,  2862, -2828,   797, -2828,  2728, -2828,  2987,   797,  2999,
    1128,  3008,  1128,  1128,  3009,  3014, -2828,   898,  3017,  1128,
    3020,  3022,  1128,  3024,  3027,  1128,  3028,  3029,  2910, -2828,
     302,  1587,  1128,  3035,  3041,  1128,  3044,  3046,  1128,  3050,
    3060, -2828,  3064,  1604, -2828,  1128, -2828,  1128, -2828,  2798,
   -2828,  3066, -2828,  3077,  3087, -2828,   797,  3090,  3094,  2970,
    3018, -2828, -2828,   797, -2828, -2828, -2828, -2828, -2828,   797,
     797,   797,  1128,  1128,  1128, -2828,  1128,  1128, -2828,  3102,
   -2828, -2828, -2828,   797,  3103, -2828,  3104, -2828,   797, -2828,
     797,   797, -2828,   797,  1128, -2828,  3108, -2828, -2828,  3112,
     797,  1607,  3114,   797,  1622,  3118,   797,  1128,  3119, -2828,
     797,  3120, -2828,  1128,  3122,  3123, -2828, -2828, -2828, -2828,
   -2828, -2828,  3124,  1636,  1637,   797,  1645,  1128,  1128,   797,
     797,  1128,  1652, -2828,  1128,  1128,  1654,   903,  3129,  1128,
    1128, -2828,  1128,   797,  1128,  1128,  1128,  1128, -2828,  1128,
    1128, -2828, -2828, -2828, -2828,   735,   692, -2828,  1128, -2828,
   -2828,  1128, -2828,  3130,  1128, -2828,  2802,  1128,  1128,  1128,
    3148, -2828, -2828,  1128,   -19, -2828,  1128,   201,  3155, -2828,
    1128, -2828,  3157,   342,  1128, -2828,  3158, -2828,  1659,  3164,
     797,  2806, -2828, -2828,  3165, -2828,  1683, -2828,  2809, -2828,
   -2828,  3167, -2828,  3169, -2828,   797,  -200, -2828, -2828, -2828,
    3171, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828, -2828,  3180, -2828, -2828,   177,   797, -2828,   797,
   -2828,   629, -2828, -2828, -2828, -2828, -2828,  3184,  3186, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828, -2828,  3188,  3194,  3195, -2828, -2828, -2828, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828,  3200, -2828, -2828,  3204, -2828, -2828, -2828, -2828,
   -2828, -2828,  3206, -2828, -2828,  3207, -2828, -2828, -2828,   797,
    2812, -2828, -2828,  1694, -2828,  1128, -2828,  3208,  1128, -2828,
     797,   797, -2828, -2828,   797,   797,  1128,  1128,  1128, -2828,
    1128,  1128, -2828,  1128,  1128,  1128,  1128,  1128,  1128,  1128,
    1128,  1128,  1128,  1128,  1128,  1128, -2828,  1128,  1128,  1128,
    1128,  1128,  1128,  1128,  1128,  1128,   797,   797,   797,   797,
    1128,  1128,  1128,  1128,  1128,  1128,  1128,  1128,  1128,  1128,
   -2828, -2828,   797,  1128,  1128, -2828, -2828,   797, -2828, -2828,
   -2828,  3211,   571,   797,  1128, -2828,  3215, -2828,  3346,  1128,
   -2828,  3225,   797,   797,   797, -2828,  1696, -2828,  3227,  1128,
    3228, -2828,  1697,  3234,  3236, -2828, -2828, -2828, -2828,  3237,
     660, -2828,  3240, -2828, -2828,   797, -2828, -2828, -2828,  1128,
   -2828,  1128, -2828,  2987, -2828,  2987,  3174,  1128,  1128,  1128,
    1128,  1128, -2828,  1128,  3241, -2828,  1128,  1128, -2828,  1128,
    1128, -2828,  1128,  1128, -2828,  3262,  1702,  3140, -2828, -2828,
    1128,  1128, -2828,  1128,  1128, -2828,  1128,  1128, -2828, -2828,
   -2828,  3270,  1734,  1736, -2828, -2828, -2828, -2828,  3274, -2828,
   -2828,  1784,  1128,  3276,   797,  1128,   797,  1128,  1128,  3277,
    1128,  1128, -2828,   797, -2828, -2828, -2828,   828, -2828,   883,
   -2828, -2828, -2828,   797, -2828,  1787, -2828,  3285, -2828,  3289,
   -2828,   797,  3291, -2828,  1128, -2828,  3292, -2828, -2828, -2828,
   -2828,  3293, -2828,  3295, -2828,  1128,  1128,  1128,   797,   741,
   -2828,  1798,  1128,  1128, -2828,  3297,  1128, -2828,   983, -2828,
    1128,  1857,   987,  3299,  1128,  1128,  1128,  3301,  1128,  1128,
   -2828, -2828,   770,   783,  3307,  3309, -2828,  1128, -2828,  2815,
    1128,  1128,  3311, -2828,  3312,  3313, -2828,  -201, -2828,  1128,
    1128,  3324, -2828,  3325, -2828, -2828,  3332,   352, -2828, -2828,
    1128, -2828,  2844, -2828,  2881, -2828, -2828,  2899, -2828,  2917,
   -2828, -2828,  -141, -2828,  3335, -2828, -2828,   797, -2828,  3336,
    3337,   797, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828,  3338, -2828,   797, -2828,  2928,  3342, -2828,   -13,
     797,  1128,   797,  1128,  1128,  1128,  1128,  1128,    65,  1128,
      74,  1128,  1858,  3343,  1128,  1128,  -193,  1128,  1128,  1128,
    1128,  1128,  1128,   144,  1128,  3344,  1128,  1128,  1128,  1128,
    1128,   797,   797,  1128,  1128,  1128,  1128,  1128,  1128,  1128,
    1128,  1128,  1128,  1128,  1128,  1128,  1128, -2828,  1128, -2828,
     797,  1128,   797,  3345, -2828,  3347,   797, -2828,   797, -2828,
    2941, -2828,  3349, -2828, -2828,  3351, -2828, -2828, -2828, -2828,
     778, -2828, -2828, -2828, -2828, -2828, -2828,  3353,  1128,  1128,
    1128,  3355,  1128,  3357, -2828,  3359,  1128,  3361,  1128,  3369,
    1128, -2828,  1128,  3370,  1859,  3371,  1128,  3373,  1128,  3374,
    1128, -2828, -2828,  3378, -2828,  3383, -2828, -2828,   797,   797,
   -2828,  1128,  3384,  1128,  1128, -2828,  1128,  1862,  3385,   797,
   -2828,   797, -2828,  3386, -2828,  3389, -2828, -2828,  1128, -2828,
    3391, -2828, -2828, -2828,  1128,  1128,  1128,  1128,  3392, -2828,
   -2828,  1128,  1128,  1128, -2828,  3393,  1128,  1128, -2828,   -55,
    1128, -2828,  3394,  1128, -2828,  1014, -2828,  1128,  1128,  1128,
   -2828,  1128,  1128, -2828, -2828,   790,   717, -2828, -2828,   797,
   -2828,  3396,  1128,  1128, -2828, -2828, -2828, -2828,  1128,  3398,
      33, -2828, -2828, -2828, -2828,  3401,  3408, -2828,  2951, -2828,
    2958, -2828,  3409, -2828,  3410, -2828,  3418, -2828,  3419, -2828,
   -2828,  3425, -2828,  3427, -2828,  3428, -2828, -2828,   797,  1128,
     811,  1128,  1128,  1128,  1128,  1128,  1128,  1128,  1128, -2828,
    1128,  1128,  1128, -2828,  1128,  1128,  1864, -2828,  1908, -2828,
    1128,  1128, -2828,  1128,  1128,  1128,  1128,  -183,  1128,  1128,
   -2828,  1128,  1128,  1973, -2828,  1128,  1128,  1993,  1128,  1128,
    1128,  1128,  1128,  1128,   165,   203,  -175,  1128,  1128,  1128,
    1128,  1128,  1128,  1128, -2828, -2828, -2828,   797,  3429,  1128,
   -2828,  3470,   797, -2828,  3433, -2828, -2828, -2828, -2828,  1128,
    1128,  3435, -2828,  3437, -2828, -2828,  3439, -2828,  3442, -2828,
    3443,  1998, -2828, -2828,  1128, -2828,  3444, -2828,  3447, -2828,
    3450, -2828, -2828,   797,   797,  3453, -2828,  1128,  1128,  1128,
   -2828,  1128, -2828, -2828, -2828, -2828, -2828,  3458, -2828,  3460,
    3467,  3468,   769, -2828,  1128,  1128,   473, -2828,  1128,  1128,
    1128, -2828,  1128, -2828,  3469,  1128, -2828,   118,  1128,  1128,
    1128,  1128,  1128, -2828, -2828,   821,  1128,  1128, -2828,  1128,
    1128,   244, -2828,  3473, -2828, -2828, -2828, -2828,  2975, -2828,
    3002, -2828, -2828, -2828, -2828, -2828, -2828, -2828,  3477,   127,
   -2828,  2003,  1128,  2020,  1128,  2021,  2026,  1128,  3479,  1128,
    -134,  3481,  1128,  -119, -2828,  1128, -2828,  2044,  1128,  1128,
    1128,  1128,  1128,  1128, -2828,  1128,  1128,  3487,  1128,  -117,
   -2828,  1128,  2083,  2097, -2828,  1128,  1128,  1128,  1128,  1128,
    1128,  1128, -2828,  1128,  1128, -2828,  1128,  1128, -2828,  1128,
     -45,  1128,  1128,  1128,  1128,  1128,  1128,  1128, -2828,  3488,
    3489, -2828, -2828,  1128,  1128, -2828, -2828, -2828, -2828, -2828,
    1128,  2106,  2205, -2828, -2828, -2828, -2828,   797,  1128, -2828,
     376,  1128,  1128,  2282, -2828, -2828, -2828, -2828,  3493, -2828,
    1128,  1128, -2828,  1128,   528,  1128,  3495,  1128,   545, -2828,
    1128,  1128, -2828,  1128,  1128,  1128,  1128,  2283, -2828,   748,
    1128,  1128,  1128,  3497, -2828, -2828, -2828,  3003, -2828,  3045,
   -2828, -2828,   797, -2828,  1128,  1128, -2828,  1128,  1128, -2828,
    2284, -2828,  2301,  1128, -2828,  1128, -2828,  1128, -2828,  1128,
   -2828,  1128,  3502, -2828,  2304,  1128,  1128,  1128,  1128,  2329,
    2379,  1128,  1128, -2828,  1128, -2828,  1128,  2390, -2828,  3503,
   -2828,  2396,  3506,  1128,  1128,  1128,  1128,  1128,  1128,  1128,
     -12,  1128,     2,  1128, -2828,  1128,  2417,  2440,  1128,  1128,
    1128,  1128,  3509, -2828,  3643,  3514,  1128, -2828, -2828,  1128,
    2492,   797,  1128,   797,  3521,  1128,  1128, -2828,  1128, -2828,
    1128,  1128,  3525, -2828,  1128,  -186,  1128, -2828,  1128, -2828,
    1128,   548,  3531,  1128,   797,  1128,  1128,  1128, -2828,  1128,
   -2828,   869,  1128,  1128,  1128, -2828, -2828,  3057, -2828,  3533,
     820,  3539,  2497,  3545,  2499, -2828,  1128, -2828,  1128,  1128,
    3546,  1128,  3548,  1128, -2828, -2828,  2526,  1128,  1128,  3550,
    1128, -2828,  2570, -2828,  2571,  1128,  1128,  3553,  1128, -2828,
    2572, -2828, -2828,  1128, -2828,  1128,  1128,  1128,  1128,  1128,
    1128,  1128, -2828,  1128,  1128, -2828,  1128,  1128,  1128, -2828,
    2574, -2828,  2575,  1128,  1128,  1128,  1128, -2828,  3555, -2828,
    3556,  2576, -2828, -2828,  3557,  3558,   797,   797,  1128,  3559,
    2598,  1128,  1128, -2828,  3563, -2828,  1128, -2828,  3566,  3567,
   -2828,  1128,   218, -2828,  1128,  3568,  1128,  1128,  1128,  1128,
   -2828,  1128,  1128,  2701, -2828,  3569, -2828, -2828,   876, -2828,
   -2828,  1128, -2828, -2828,  1128,  2720,  2778,  2805, -2828,  1128,
   -2828,  1128, -2828,  3110,  1128,  1128, -2828,  1128, -2828,  3144,
   -2828,  3147,  3571,  1128, -2828,  1128, -2828,  2811,  2882,  1128,
    1128,  1128,  1128,  1128,  1128,  3581,  1128,  3583,  1128,  3590,
    1128, -2828,  1128, -2828,  1128,  1128,  1128,  1128,  1128, -2828,
   -2828, -2828,  2933, -2828, -2828,   797,  3068, -2828, -2828,  3591,
    3074,  3593, -2828,  3595, -2828, -2828,  3599, -2828,  1128,  3601,
   -2828,  1128,  1128,  3603,  1128,  1128,  3111, -2828,  1128, -2828,
   -2828,  3606,  3607, -2828,  3196, -2828,  3209, -2828,  3280,  3609,
    3611, -2828,  3203,   254,   347,    27, -2828,  3612, -2828,  3620,
   -2828,  1128,  3624, -2828,  3282, -2828,  3298,  1128,  1128,  1128,
    1128,  1128,  1128, -2828,  1128, -2828,  1128, -2828,  3625,  3300,
    3339,  1128,  1128,  1128,  1128, -2828, -2828,  1128, -2828, -2828,
    3372, -2828, -2828, -2828,  3628, -2828,   797,  3629, -2828,  1128,
   -2828, -2828,  1128,  1128, -2828, -2828, -2828,  3250, -2828,  3278,
   -2828,  3416, -2828, -2828, -2828,   797, -2828,  1128,  1128, -2828,
    1128,  1128, -2828,  1128, -2828, -2828,    43, -2828, -2828,  3431,
   -2828,  3286,  3463,  3554,  1128,  1128,  1128,  1128,  3632,  3635,
   -2828, -2828,  3560, -2828,  3578,   419,   424,   119,  1128,  3585,
   -2828,  3641, -2828,   797, -2828,   887,  1128,  3642, -2828,  3646,
   -2828,  3647, -2828,  3651,   797,  1128,   124,  1128,   173,  1128,
   -2828,  1128, -2828,  3594, -2828,  3653, -2828,  1128, -2828,  1128,
    3638,  3644,  1128,  1128, -2828, -2828, -2828,  3315, -2828,  3328,
   -2828,  1128,  1128, -2828,  1128,  1128, -2828,  1128,   205, -2828,
    3654, -2828,   797,  1128, -2828,  1128,  3655, -2828, -2828, -2828,
   -2828,   797,  1128, -2828,  1128,  1128, -2828,  1128,  1128,  1128,
   -2828,  3645, -2828,  3657,  3659, -2828,  1128, -2828,  1128,  1128,
    1128, -2828,  3660, -2828,  3661,  1128,   212,  1128,   213,  1128,
   -2828,  1128, -2828,  3662,  3663,   892, -2828,  1128,  3665,  1128,
    3669,  1128,  3670,  1128, -2828,  3648, -2828, -2828,  3671,  3673,
    3664,  3672, -2828, -2828,  1128, -2828,  1128,  1128, -2828,  1128,
    1128,  1128, -2828, -2828,  1128, -2828,  1128,  1128, -2828,  1128,
   -2828,  1128, -2828,  3674, -2828,  3689, -2828, -2828, -2828,  1128,
   -2828,  1128,  3675,  1128,  3683,  1128,  3685,  1128,  1128,  1128,
    1128,  3686,  3687, -2828, -2828,  3724,  3688,  3690, -2828,  1128,
   -2828,  1128, -2828,  3691,  3692,  1128,   797, -2828, -2828, -2828,
    3734, -2828, -2828,  3693,  3694, -2828, -2828,   894,  1128, -2828,
    3735, -2828, -2828,  1128, -2828,  1128,  3695, -2828,  3752,  1128,
     797, -2828, -2828,  3787,  3697,  1128, -2828,  3795, -2828,  1128,
   -2828,  3698,   797, -2828,   797,  1128,  1128,  1128,  3699, -2828
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint16 yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     3,   127,   128,   129,   130,   131,   116,
      79,   122,   123,   124,    55,    56,    57,    52,    53,    54,
      51,    69,    74,    75,    76,    38,    39,    41,    40,    49,
      42,    43,    44,    46,    88,    87,    95,   322,    45,    48,
      81,    82,    83,    84,   126,    85,    86,    67,    68,    73,
      70,    71,    72,   143,   113,   114,   115,   112,     6,   141,
      92,    93,    89,    90,    91,    94,   105,   106,   101,   107,
     108,   109,   110,   117,   118,   119,   120,   121,   102,   103,
      31,    30,    34,    32,    33,    35,    36,    37,     7,     8,
      58,    59,    60,    61,    62,    64,    63,    65,    66,    10,
       9,    11,   111,    22,    12,   134,   135,   137,   136,   138,
      47,   139,   140,   144,     5,    15,    13,    14,   142,    16,
      17,    18,    20,    21,    19,    25,    26,    27,    28,    78,
      23,    24,    96,   145,    98,   104,    99,   100,    97,    50,
      80,    77,   125,   132,   133,    29,   146,   147,   148,   149,
     151,   150,   152,     0,     0,     0,     0,  1425,  1426,     0,
     836,     0,  1428,  1430,  1427,  1429,     0,     0,     0,     0,
     334,     0,     0,     0,     0,     0,   834,   558,     0,   234,
     488,     0,   228,   356,  1138,   761,     0,   468,   822,     0,
       0,  1104,   260,   463,   361,   194,     0,  1071,  1076,     0,
       0,     0,   854,     0,   323,     0,   824,     0,     0,     0,
     332,     0,     0,     0,   484,     0,   570,     0,   209,     0,
     752,   557,   264,   420,     0,   155,     0,     0,     0,     0,
     217,     0,  1082,     0,     0,     0,     0,     0,   415,   434,
     653,   548,     0,   554,     0,     0,   563,     0,     0,     0,
       0,     0,     0,     0,     0,   254,   639,     0,     0,   220,
       0,   341,   859,   885,     0,     0,   353,     0,     0,   214,
       0,   427,     0,     0,  1107,     0,     0,     0,     0,   828,
     893,     0,     0,   280,     0,     0,     0,   281,     0,   391,
       0,     0,     0,     0,   888,   872,     0,     0,     0,     0,
     776,   492,     0,   429,     0,     0,     0,     0,  1137,   266,
     159,   634,   629,     0,     0,     0,   920,   327,     0,     0,
     479,     0,     0,   357,     0,     0,   412,   411,   597,   741,
     343,     0,   275,     0,   592,   602,   607,   624,   614,   619,
     183,  1141,     0,   413,     0,   161,     0,  1149,  1305,     0,
       0,     0,   207,   351,     0,   416,   435,     0,     0,     0,
     758,  1315,  1420,  1355,  1360,     0,   869,   864,   866,  1351,
    1353,     0,     1,     2,     4,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   172,   173,   174,
     170,   171,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   229,   230,   231,   232,   233,   235,     0,
     255,     0,     0,     0,     0,     0,     0,   276,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   362,   363,   364,     0,     0,
       0,   392,   393,   408,     0,     0,     0,     0,     0,     0,
     460,     0,     0,   464,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   502,     0,   513,     0,   516,     0,   519,
       0,   522,     0,   531,   533,   535,     0,   546,     0,     0,
       0,     0,     0,     0,     0,     0,   572,     0,     0,   668,
       0,   724,     0,     0,     0,     0,   756,     0,     0,   767,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   805,     0,
     823,   825,   829,     0,     0,     0,     0,     0,     0,     0,
       0,   860,   861,     0,  1427,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   961,   921,  1091,  1090,  1089,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1133,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1310,  1310,
    1310,  1310,  1310,     0,     0,     0,     0,     0,  1310,     0,
       0,     0,  1341,     0,     0,  1350,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1357,
    1358,  1359,     0,     0,     0,     0,   273,   582,     0,   399,
       0,     0,   193,   837,   530,   532,     0,   335,     0,     0,
     525,   527,     0,   528,     0,   344,  1083,     0,   489,     0,
       0,     0,     0,     0,     0,  1079,  1072,     0,  1077,     0,
       0,  1069,   855,   324,   512,   501,   569,   590,     0,  1067,
       0,   536,     0,     0,   485,     0,   571,   339,     0,     0,
     455,   665,     0,   663,     0,   495,   496,     0,   218,   515,
    1100,     0,  1102,     0,   521,   518,   654,     0,     0,   550,
     549,   553,   567,   564,     0,     0,     0,   459,     0,     0,
     340,     0,     0,   640,     0,     0,     0,   270,     0,   221,
       0,   886,   354,   887,   667,   215,   426,   326,   803,     0,
       0,   329,   288,   284,     0,   282,   289,   285,     0,   283,
       0,   559,   561,     0,   873,   345,     0,     0,     0,   493,
     428,     0,   543,     0,   545,     0,     0,   635,     0,   630,
     272,     0,   333,     0,   506,   507,     0,   330,   331,   723,
       0,     0,   598,   271,   274,     0,   593,     0,   603,     0,
     608,     0,   625,     0,   615,     0,   620,     0,   414,   162,
     725,     0,     0,     0,   740,     0,   352,   470,   471,     0,
     755,   475,   476,     0,     0,   325,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   177,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   205,   201,
       0,   204,   202,     0,   206,     0,   199,   200,   198,   197,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   243,     0,   245,   247,     0,
     249,     0,     0,     0,     0,     0,     0,     0,     0,   279,
       0,     0,     0,     0,     0,     0,     0,   319,   320,   318,
     321,   310,   307,   308,   317,   314,   290,     0,   316,   311,
       0,   304,   305,   312,     0,     0,   348,   350,   349,     0,
     388,   377,     0,     0,     0,   390,     0,   365,   358,     0,
       0,   382,   372,     0,     0,     0,     0,     0,     0,     0,
     366,   359,     0,   394,     0,     0,   404,     0,     0,   407,
     409,   432,   431,   433,     0,   467,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   577,     0,   643,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   763,
       0,     0,     0,   792,   791,     0,     0,     0,     0,     0,
       0,   788,     0,   789,   790,     0,   781,   783,   778,   779,
       0,   793,   786,     0,   780,     0,   787,     0,     0,     0,
       0,     0,     0,     0,     0,   858,   857,     0,     0,     0,
     865,     0,   868,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   918,     0,   916,
       0,     0,   958,     0,     0,   931,   933,     0,     0,     0,
     945,     0,     0,     0,   953,     0,   939,     0,     0,   923,
       0,     0,   935,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1109,     0,     0,     0,  1128,     0,     0,  1105,
       0,  1106,     0,  1108,  1110,     0,     0,     0,  1121,  1111,
       0,     0,     0,     0,     0,  1308,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1333,  1342,  1348,
    1347,  1340,  1329,  1331,  1344,  1317,  1318,  1319,  1320,  1321,
    1334,  1322,  1326,  1325,  1324,  1323,  1327,  1332,  1412,     0,
       0,     0,  1316,     0,     0,  1380,  1392,  1391,  1378,  1409,
    1379,  1384,  1386,  1385,  1387,  1388,  1369,  1370,  1389,  1390,
    1362,  1365,  1371,  1372,  1368,  1400,  1401,     0,  1399,  1381,
    1402,  1403,  1393,  1396,  1404,  1405,  1375,  1376,  1377,  1406,
       0,     0,  1352,  1354,  1415,  1418,  1356,     0,     0,     0,
    1361,  1422,  1424,  1421,   583,     0,     0,   400,     0,     0,
       0,   336,   337,     0,   526,     0,   529,  1084,     0,     0,
       0,   762,   380,     0,   347,  1073,  1078,  1070,  1080,   591,
       0,  1068,   537,   538,  1097,     0,     0,     0,     0,   419,
       0,     0,   666,     0,   664,   497,     0,  1101,  1103,     0,
     656,     0,   655,     0,     0,   660,     0,  1088,  1092,     0,
       0,     0,     0,   154,     0,     0,     0,   645,     0,   644,
       0,   647,     0,     0,     0,     0,     0,     0,   286,   287,
       0,     0,   346,  1086,  1087,   430,   544,  1312,     0,     0,
       0,     0,     0,     0,   960,     0,   508,     0,   417,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1308,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1142,   726,     0,   730,     0,     0,     0,   739,   472,   474,
       0,   477,     0,     0,     0,     0,     0,     0,   179,     0,
       0,     0,     0,     0,   176,     0,     0,   163,     0,     0,
     184,     0,   186,   188,   189,   190,   191,   192,   195,   203,
     196,   208,   210,   213,   212,   211,   216,   219,     0,     0,
     225,     0,     0,     0,   242,   244,   246,   248,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   268,   269,
       0,     0,   815,     0,   291,     0,   297,     0,   309,   315,
     303,   313,   306,   342,   384,   386,     0,     0,   385,   370,
     383,   454,   588,     0,   379,   368,   367,     0,   373,     0,
       0,   371,   360,     0,   395,     0,   405,   406,     0,     0,
       0,     0,     0,   480,     0,   486,     0,     0,     0,     0,
     504,     0,     0,     0,     0,     0,     0,     0,   551,     0,
       0,     0,     0,   565,     0,   568,     0,     0,     0,     0,
       0,   684,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   812,     0,     0,     0,     0,     0,     0,   764,
       0,   768,   777,   797,   798,   799,   800,   801,   794,   802,
     784,   785,   782,   795,   796,   808,     0,     0,     0,     0,
       0,     0,   840,     0,     0,   856,     0,   863,   867,   870,
       0,     0,     0,     0,     0,     0,   874,     0,     0,   889,
       0,   899,   900,   894,   895,   896,   898,   901,   919,   902,
     917,   903,     0,   897,   925,   922,   932,   934,   929,   950,
     952,   951,     0,   946,     0,     0,   940,     0,   930,     0,
     959,   936,   937,   938,   926,     0,     0,     0,     0,  1034,
       0,  1033,     0,     0,     0,  1035,     0,     0,   987,     0,
       0,  1009,  1011,     0,     0,     0,     0,     0,  1047,     0,
       0,  1032,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1008,     0,     0,     0,     0,     0,  1040,     0,  1000,
       0,     0,     0,     0,     0,  1030,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1048,
       0,     0,     0,     0,     0,     0,     0,     0,   962,     0,
       0,  1012,     0,     0,  1037,     0,     0,     0,     0,  1119,
       0,  1129,  1122,  1123,     0,  1113,  1114,     0,  1134,     0,
    1132,  1112,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1306,     0,  1343,
    1349,  1330,  1345,  1311,  1328,     0,  1337,  1339,     0,  1413,
    1364,  1363,  1366,  1373,     0,  1411,  1382,  1394,  1397,  1407,
    1408,  1417,     0,  1419,  1423,     0,     0,   403,     0,     0,
     534,   338,   524,  1085,   355,     0,     0,  1074,  1081,   819,
    1098,     0,   578,     0,     0,     0,   418,     0,   424,     0,
       0,     0,   657,     0,     0,   659,  1093,     0,   457,     0,
       0,     0,     0,   642,   646,     0,   648,     0,     0,   641,
     222,     0,   223,     0,   153,   441,   437,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   328,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1308,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   727,     0,     0,   731,     0,   732,     0,   473,     0,
    1414,     0,   156,     0,     0,   175,     0,     0,     0,     0,
       0,   178,   181,     0,   180,   185,   187,   227,   226,   236,
       0,     0,     0,     0,     0,   827,     0,     0,   263,     0,
     261,   265,   267,     0,     0,   817,     0,   814,     0,   293,
       0,     0,   299,     0,   389,   378,     0,   452,   589,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   483,
       0,     0,   487,     0,     0,     0,   503,   505,   514,   517,
     520,   523,     0,     0,     0,   552,     0,     0,     0,   566,
       0,     0,     0,   721,     0,     0,     0,     0,     0,     0,
       0,   685,     0,     0,     0,     0,     0,     0,   718,     0,
       0,   811,   813,   729,   742,     0,     0,   753,     0,   757,
     759,     0,   766,     0,   769,   807,     0,     0,     0,     0,
       0,   838,   850,     0,     0,   851,     0,     0,     0,   871,
       0,   876,     0,     0,     0,   875,     0,   878,     0,     0,
       0,     0,   949,   947,     0,   954,     0,   941,     0,   924,
     927,     0,  1016,     0,  1022,     0,     0,  1064,  1015,  1013,
       0,  1027,  1065,  1066,  1007,  1006,  1010,  1025,  1026,   981,
     984,  1046,  1045,     0,   985,   986,     0,     0,  1054,     0,
    1053,     0,  1052,  1051,  1005,  1004,  1061,     0,     0,   979,
     980,  1041,  1042,   988,   989,   990,  1031,   964,  1001,   992,
     991,  1039,  1062,     0,     0,     0,  1028,  1003,  1002,   969,
     972,   970,   971,   973,   974,   975,   976,   967,   968,   966,
    1044,   995,     0,   983,   993,     0,   982,   965,  1043,  1029,
     963,  1036,     0,  1049,  1023,     0,  1038,  1096,  1117,     0,
       0,  1120,  1124,     0,  1115,     0,  1135,     0,     0,  1303,
       0,     0,  1304,  1309,     0,     0,     0,     0,     0,  1174,
       0,     0,  1180,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1186,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1307,  1346,     0,     0,     0,  1367,  1374,     0,  1383,  1395,
    1398,     0,     0,     0,     0,   401,     0,   410,     0,     0,
    1075,     0,   580,   579,     0,   421,     0,   821,     0,     0,
       0,   498,     0,     0,     0,   456,   458,   462,  1095,     0,
       0,   650,     0,   804,   442,     0,   444,   445,   446,     0,
     448,   449,   451,     0,   438,     0,  1313,     0,     0,     0,
       0,     0,   632,     0,     0,   509,     0,     0,   600,     0,
       0,   595,     0,     0,   605,     0,     0,     0,  1308,   611,
       0,     0,   627,     0,     0,   617,     0,     0,   622,   728,
     733,     0,     0,     0,   478,   158,   157,   160,     0,   164,
     165,     0,     0,     0,     0,     0,   238,     0,     0,     0,
       0,     0,   262,     0,   277,   816,   294,   292,   300,   298,
     387,   453,   818,     0,   374,     0,   436,     0,   396,     0,
     461,     0,     0,   469,     0,   481,     0,   490,   494,   547,
     539,     0,   540,     0,   556,     0,     0,     0,     0,     0,
     680,     0,     0,     0,   687,     0,     0,   706,     0,   692,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     747,   743,     0,     0,     0,     0,   765,   770,   806,     0,
       0,     0,     0,   839,     0,     0,   841,     0,   844,     0,
       0,     0,   862,     0,   880,   881,     0,     0,   879,   890,
       0,   891,     0,   904,     0,   948,   955,     0,   942,     0,
     928,  1017,     0,  1018,     0,  1014,  1063,     0,  1058,     0,
       0,     0,  1057,   977,   978,   997,   999,   998,   996,   994,
    1050,  1024,     0,  1130,     0,  1125,     0,     0,  1136,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1335,     0,  1416,
       0,     0,     0,     0,   402,     0,     0,  1099,   581,   423,
       0,   820,     0,   499,   662,     0,   658,  1094,   652,   649,
       0,   224,   443,   447,   450,   440,   439,     0,     0,     0,
       0,     0,     0,     0,   510,     0,     0,     0,     0,     0,
       0,  1308,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   734,   735,     0,   737,     0,   166,   169,     0,     0,
     182,     0,     0,     0,     0,   253,     0,     0,     0,     0,
     295,     0,   301,     0,   375,     0,   398,   397,     0,   465,
       0,   491,   541,   542,     0,     0,     0,     0,     0,   573,
     681,     0,     0,     0,   689,     0,     0,     0,   710,     0,
       0,   688,     0,     0,   708,     0,   693,     0,     0,     0,
     719,     0,     0,   748,   744,     0,     0,   754,   760,   771,
     810,     0,     0,     0,   835,   842,   843,   847,     0,     0,
       0,   852,   877,   883,   882,     0,     0,   905,     0,   906,
       0,   956,     0,   943,     0,  1019,     0,  1020,     0,  1056,
    1055,     0,  1118,     0,  1126,     0,  1116,  1143,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1173,
       0,     0,     0,  1179,     0,     0,     0,  1295,     0,  1226,
       0,     0,  1255,     0,     0,     0,     0,     0,     0,     0,
    1185,     0,     0,     0,  1284,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1336,  1338,  1410,     0,     0,     0,
     584,     0,     0,   422,     0,   500,   661,   651,  1314,     0,
       0,     0,   631,     0,   511,   599,     0,   594,     0,   604,
       0,     0,  1308,  1308,     0,   626,     0,   616,     0,   621,
       0,   736,   738,     0,     0,     0,   240,     0,     0,     0,
     256,     0,   278,   296,   302,   369,   376,     0,   482,     0,
       0,     0,     0,   575,     0,     0,     0,   707,     0,     0,
       0,   714,     0,   690,     0,     0,   712,     0,     0,     0,
       0,     0,     0,   749,   745,     0,     0,     0,   809,     0,
       0,     0,   853,     0,   845,   884,   892,   907,     0,   908,
       0,   957,   944,  1021,  1060,  1059,  1131,  1127,     0,     0,
    1146,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1190,     0,  1296,     0,     0,     0,
       0,     0,     0,     0,  1261,     0,     0,     0,     0,     0,
    1267,     0,     0,     0,  1192,     0,     0,     0,     0,     0,
       0,     0,  1202,     0,     0,  1206,     0,     0,  1210,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   585,     0,
       0,   381,   425,     0,     0,   637,   633,   601,   596,   606,
       0,     0,     0,  1308,   628,   618,   623,     0,     0,   241,
       0,     0,     0,     0,   466,   555,   560,   562,     0,   574,
       0,     0,   694,     0,     0,     0,     0,     0,     0,   709,
       0,     0,   716,     0,     0,     0,     0,     0,   750,     0,
     772,     0,     0,     0,   848,   846,   909,     0,   910,     0,
    1145,  1144,     0,  1194,     0,     0,  1196,     0,     0,  1150,
       0,  1156,     0,     0,  1177,     0,  1172,     0,  1183,     0,
    1178,     0,     0,  1297,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1189,     0,  1184,     0,     0,  1285,     0,
    1166,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1212,     0,     0,     0,     0,     0,
       0,     0,     0,   586,  1139,     0,     0,  1308,   609,     0,
       0,     0,     0,     0,     0,     0,     0,   257,     0,   576,
       0,     0,     0,   695,     0,     0,     0,   711,     0,   697,
       0,     0,     0,     0,     0,     0,     0,     0,   673,     0,
     746,     0,   773,     0,     0,   849,   911,     0,   912,     0,
       0,     0,     0,     0,     0,  1151,     0,  1157,     0,     0,
       0,     0,     0,     0,  1191,  1298,     0,     0,     0,     0,
       0,  1287,     0,  1291,     0,     0,     0,     0,     0,  1268,
       0,  1286,  1167,     0,  1193,     0,     0,     0,     0,     0,
       0,     0,  1203,     0,     0,  1207,     0,     0,     0,  1214,
       0,  1220,     0,     0,     0,     0,     0,   587,     0,   636,
       0,     0,  1308,   612,     0,     0,   237,     0,     0,     0,
       0,     0,     0,   700,     0,   696,     0,   722,     0,     0,
     698,     0,     0,   713,     0,     0,     0,     0,     0,     0,
     751,     0,     0,     0,   913,     0,   914,  1147,     0,  1195,
    1198,     0,  1197,  1200,     0,     0,     0,     0,  1176,     0,
    1182,     0,  1299,     0,     0,     0,  1256,     0,  1288,     0,
    1292,     0,     0,     0,  1188,     0,  1269,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1215,     0,  1221,     0,     0,     0,     0,     0,  1140,
     638,   610,     0,   167,   168,   239,     0,   826,   258,     0,
       0,     0,   701,     0,   715,   703,     0,   699,     0,     0,
     691,     0,     0,     0,     0,   774,     0,   830,     0,   915,
    1148,     0,     0,  1152,     0,  1158,     0,  1162,     0,     0,
       0,  1300,     0,     0,     0,     0,  1289,     0,  1293,     0,
    1262,     0,     0,  1270,     0,  1168,     0,     0,     0,     0,
       0,     0,     0,  1204,     0,  1208,     0,  1211,     0,     0,
       0,     0,     0,     0,     0,   613,   250,     0,   259,   677,
       0,   720,   702,   704,     0,   717,     0,     0,   683,     0,
     775,   832,     0,     0,  1199,  1201,  1153,     0,  1159,     0,
    1163,     0,  1175,  1181,  1301,     0,  1227,     0,     0,  1247,
       0,     0,  1257,     0,  1290,  1294,     0,  1187,  1271,     0,
    1169,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1213,  1216,     0,  1222,     0,     0,     0,     0,     0,     0,
     678,     0,   705,     0,   686,     0,     0,     0,  1154,     0,
    1160,     0,  1164,     0,     0,     0,     0,     0,     0,     0,
    1263,     0,  1272,     0,  1170,     0,  1231,     0,  1233,     0,
       0,     0,     0,     0,  1205,  1209,  1217,     0,  1223,     0,
    1243,     0,     0,  1251,     0,     0,  1259,     0,     0,   251,
       0,   679,     0,     0,   674,     0,     0,   831,  1155,  1161,
    1165,     0,     0,  1228,     0,     0,  1248,     0,     0,     0,
    1273,     0,  1171,     0,     0,  1235,     0,  1237,     0,     0,
       0,  1218,     0,  1224,     0,     0,     0,     0,     0,     0,
    1265,     0,   252,     0,     0,     0,   833,     0,     0,     0,
       0,     0,     0,     0,  1274,     0,  1232,  1234,     0,     0,
       0,     0,  1219,  1225,     0,  1244,     0,     0,  1252,     0,
       0,     0,   682,   675,     0,   669,     0,     0,  1229,     0,
    1249,     0,  1258,     0,  1275,     0,  1236,  1238,  1239,     0,
    1241,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1264,  1276,     0,     0,     0,  1245,     0,
    1253,     0,  1260,     0,     0,     0,     0,  1230,  1250,  1277,
       0,  1240,  1242,     0,     0,  1266,   670,     0,     0,  1278,
       0,  1246,  1254,     0,   671,     0,     0,  1279,     0,     0,
       0,  1302,  1280,     0,     0,     0,  1281,     0,   672,     0,
    1282,     0,     0,  1283,     0,     0,     0,     0,     0,   676
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
   -2828, -2828, -2828,  3720, -2828, -2828, -2828, -2828, -2828, -2828,
    3649, -2828, -2828,  3650, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828, -2828, -2828, -2828, -2175, -2828, -2828, -2828, -2828,
    3354,  3358, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,  3592, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828,  3597, -2828,  3326, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2081, -2828,  3656, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828,  2912, -2828, -2828,  2906, -2828, -2828, -2828, -2828,
   -2828, -2828,  3015, -2828,  -135, -1170, -2828, -2828,  2947, -2828,
    3587, -2828,  -178, -2828,  3565, -2828,  -192,  -884,  -325, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828,  3486, -2828, -2828, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828, -2828,   455, -1196,  3551,
   -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828,  3011,  -954, -2828, -2828,  3019,  -935, -2828,  -378, -2828,
    3537, -2827, -2828, -2828, -2828, -2828, -2828, -2828, -2828,  3449,
   -2828, -2828, -2828, -2828, -2828,  -441,  3484,  2761,  -187, -1725,
    -729,  -916, -2828, -2828,   847, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828,   182, -2828,  2971, -2828,
   -2828, -2828,  -338, -2828,   186,  -746, -2828, -2828,  2423, -2828,
   -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828, -2828, -1635,   643, -2828, -2828, -2828, -2828, -2828, -2828,
   -2828,  1144,  3668, -2828, -2828, -2828, -2828, -2828,  3781,    -6
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,   171,   172,   173,   174,   175,   176,   177,   178,   557,
     558,   559,   560,   561,   179,   180,   181,   182,   183,   184,
     185,   186,   187,   188,   593,  2179,   594,   595,   596,   597,
    1105,  1108,   189,   600,   190,   191,   192,   193,   194,   195,
     196,   197,   198,   199,   615,   200,   201,   202,   617,   203,
     204,   205,   206,   634,   207,   208,  1483,   635,  1149,  1154,
     209,   641,   642,   210,   647,   211,   212,   213,   214,   215,
     216,   217,   218,   219,   220,   649,   221,   222,   636,   223,
    2116,  2513,   637,   224,   225,   226,   650,   227,   228,   229,
     230,  1047,  1048,   231,  1051,  1052,   232,   233,   234,   235,
     236,   935,   936,   237,   663,  1769,   238,  1014,  1015,   239,
     665,   240,   667,   241,   669,   242,   671,   890,   891,   243,
     244,   245,   246,   247,   248,   249,   677,   250,   251,   252,
     253,   254,   255,   256,   257,   258,   259,   876,  1741,   916,
     260,  1026,   261,  1022,   262,  1028,   263,  1030,   264,  1034,
     265,  1036,   266,  1032,   267,  1009,   268,  1007,   269,   270,
     963,   964,  1597,   271,   946,   947,  1580,   272,   930,   273,
     689,  2845,   274,   275,   276,   277,   278,   279,   280,   696,
     281,   282,   700,   283,   284,   728,   691,  1801,   877,  1742,
     917,   931,   285,   286,   598,   287,   732,   288,   289,   290,
     291,   741,   742,   292,   293,   294,   295,   296,   297,   298,
     299,  1861,  1288,  1286,   300,  1294,   774,   301,   775,   302,
     919,   303,   371,   304,  1587,  1588,   305,  1563,  1564,   306,
     940,   307,   942,   308,  1400,   309,   795,   310,   311,   312,
     313,   314,  1996,  1465,   315,   316,   824,   317,   318,   319,
     320,   864,   825,   321,   870,   871,   322,   875,  1743,  2846
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
     336,   337,   338,   341,   978,  2141,  2566,  1544,  1778,  1544,
    1599,   607,  1582,  1783,  1572,  1785,  2208,  1572,   893,   933,
    3029,  1886,   780,  1202,  1299,   366,  3034,  1771,  1772,  1773,
    1774,   569,   909,   643,   366,   739,  2514,   366,   324,   327,
     995,   653,   910,   394,   327,   676,   332,  2249,  2877,  2673,
     332,  1653,   644,  2208,   581,   675,  2932,  1806,  2208,   414,
    2208,  2154,   366,  3375,   411,   332,  3094,  2156,   674,   367,
     427,   686,  2280,  1473,  3118,  2257,  2110,  1399,   367,  3030,
     619,   367,   333,   446,   562,   563,   333,   368,   369,   570,
     332,   571,  2645,   781,   323,   332,   368,   369,  2151,   368,
     369,   333,   328,   327,  1892,   620,   367,   328,  2895,   910,
     673,   782,   485,   327,   389,  3206,  1809,   491,   492,   968,
    2138,   327,  1215,   332,   368,   369,   333,  1216,   327,   621,
    3210,   333,  3225,   621,   783,   327,   326,  1789,   327,   736,
     342,   332,   355,   572,  3053,  1893,   702,   524,   346,  2250,
     332,   334,   873,   874,   327,   334,   703,   784,   535,   333,
    1831,   332,   347,   335,   327,   582,   328,   335,   704,   428,
     334,   705,   706,   707,   708,   709,   328,   333,   585,   586,
    1217,  2152,   335,  2251,   328,   604,   333,  1559,   332,  1300,
    2258,   328,   583,  3376,  3031,   334,   618,   333,   328,   332,
     334,   328,  3166,  3167,  3244,   327,  2878,   335,  3170,   327,
    2281,  2111,   335,   785,   944,   710,   969,   328,   622,  2419,
     332,   786,   327,  1559,   333,   939,  1218,   328,   334,   787,
    2646,  1887,  1832,   332,   623,   333,  2907,  3342,  1131,  1790,
     335,   332,   511,   915,  1132,  1229,   334,   412,  2106,   327,
    2674,  3345,  3171,  1203,   332,   334,   333,   327,   335,   327,
     788,  2139,   789,  2422,   688,  2933,   334,   335,   328,   333,
     790,   584,   328,  1791,   332,  3095,  3562,   333,   335,   332,
     325,   624,  3054,  3119,   711,   328,   743,   746,   712,   327,
     333,  1133,  3610,   334,   626,   996,  2436,   625,   626,   779,
     332,  1654,  1655,  1656,   334,   335,   627,  1810,  2315,  2896,
     333,  2155,   328,  2252,  2919,   333,   335,  2157,   797,  1474,
     328,   628,   328,  2923,  3207,   334,  1164,  1219,   791,  1263,
     332,   800,   886,   356,   629,   889,   333,   335,   334,  3211,
    3278,  3226,   897,  3282,  3283,   699,   334,   899,   349,   737,
     335,   970,   328,   903,   618,  3183,   553,  2208,   335,   334,
     905,   906,  1231,   908,   645,   352,   333,  3172,  3636,  1121,
    2082,   335,   437,  3653,  2253,   920,  3191,   801,   802,   334,
     555,   327,   327,   713,   334,  1186,  1189,   928,   929,  1301,
    1833,   335,   714,  2940,  2220,   715,   335,  1111,   941,   943,
     716,   717,   718,  1792,   719,   334,   327,   720,  2677,   435,
    2537,   721,   955,  3245,  3112,   958,   573,   335,   961,  1119,
     740,   966,  3656,   332,   792,   574,  2678,  2908,  2420,  2909,
     793,   630,  2775,   512,  2776,   334,  1230,   353,   646,  2107,
     979,   631,  1134,   632,   328,   328,  3343,   335,  1888,   475,
    2648,  3378,  3115,  2193,  3680,   547,  3384,   722,   354,   333,
    3346,  3705,  3708,   910,  1793,  1185,   794,  3467, -1310,   328,
    1135,   327,  2423,  1136,   564,   565,   566,   567,   568,   357,
    1137,   723,  1013,  1040,   358,  3563,   724,  1846,  1213,  1020,
     633,  1220,  2159,  3184,  1169,   332,  1138,  1278,  1212,  2283,
     332,  3611,   327,  3556,  2536,  2437,   575,   576,   577,  1679,
    2208,   578,   579,   609,  2208,   548,   549,   550,   551,   552,
    1536,  1221,   359,  2920,   332,  1974,   360,   761,   334,  1054,
    1302,   333,  2924,   436,   328,  1841,   333,   911,  1619,  2285,
     335,  1059,  1060,  1061,  1062,  1063,  1064,  1065,  1067,   332,
    1069,  2538,   361,  1681,  1070,  1717,  1574,  3469,  1139,  1140,
     333,   610,  1934,   751,   752,   328,  1041,  3192,   362,  1546,
     764,  2291,  2295,   476,  1090,   327,  1639,  3637,   606,  1095,
    1096,  1099,  3654,  2235,   366,   333,  2043,   327,  2239,  1141,
    1142,  2655, -1310,  1112,   363,  1113,  3559,  3468,  1116,   467,
     334,  2884,  2941,  2649,   332,   334,   327,   766,  2650,  1144,
    1847,  1042,   335,  1150,  1152,  1153,  1155,   335,   753,  1150,
     327,   332,   725,  3113,   332,  3263,  1167,   327,   367,   334,
     327,  3657,  1172,  1680,  1175,   726,   727,  2065,   328,  2312,
     333,   335,  2104,   501,  2092,  1184,   368,   369,  1975,   767,
     328,  1190,   327,  1192,   334,  1194,  1195,   333,  1572,   364,
     333,  3116,  1544,  3681,  2392,   327,   335,   754,  3630,   328,
    3706,  3709,   768,  3633,  1206,  1207,  1208,  1682,   365,  1718,
    1143,   376,   553,   328,  1222,   554,  1935,   377,  1170,  1534,
     328,   803,  2286,   328,  2292,  2296,  2216,  1701,   366,  1236,
    1237,  1238,  1239,  1240,   327,  2313,   555,   366,   327,   334,
    2063,   378,  3557,  1250,  2656,   328,   379,   327,  1255,   327,
    1257,   335,  3162,  1259,  2885,   380,   334,  1261,   328,   334,
     366,   381,  1267,   327,   611,   612,   613,   614,   335,  2393,
     755,   335,   367,  1277,  1279,  1280,  1281,   756,  3264,  1284,
    1285,   367,   327,   388,  1291,   366,  1293,  2208,  2208,   761,
     368,   369,  1534,   390,  1310,   757,  1534,   328,   332,   368,
     369,   328,   652,  1391,   367,  1393,   541,  3273,  1397,  1398,
     328,  2240,   328,   762,  1405,   588,   366,  1899,   804,   763,
     805,   806,   368,   369,  3279,   370,   328,  3380,   807,   367,
    1117,  2196,   764,   656,   333,  3560,  1312,  1470,  2197,   765,
     808,   809,   810,   811,   812,   328,   813,   368,   369,   903,
    1702,   366,  1486,  1534,   658,   910,   814,   815,   556,   816,
     367,   817,   818,   771,   772,  1499,   819,  1534,   391,   766,
    1504,  2096,   820,   821,   822,   823,  2473,   392,   368,   369,
     327,   395,  3163,  2336,  2337,  1313,   327,   396,   366,   332,
    2681,  1919,   660,   903,  1528,   367,   949,  2096,  2096,  1123,
     451,   953,  2338,   334,  1538,  1539,  1125,  3631,  2682,  1794,
    1540,   767,  3634,   368,   369,   335,  1545,  2482,  2483,  2206,
    1946,  1548,   605,   327,  1584,   333,  2207,  1553,   327,  1557,
    1555,  1556,   367,  2794,   768,   685,   910,  3274,   332,  2769,
     910,   399,  1560,   328,   748,  2198,  1056,  1565,  2264,   328,
     368,   369,  1570,  1571,  3280,   327,   327,  3381,  1314,   758,
     759,   760,   327,   402,   589,   590,   591,   592,  1010,   327,
     452,  2631,  1586,  1584,   333,   910,  2750,   403,  1589,  1590,
    1591,  2201,  1592,   910,   769,  1594,   328,   332,  1024,  1602,
    1603,   328,  3643,   332,   334,   910,  3044,  3714,   332,  3773,
     332,   404,   332,  1606,   332,   327,   335,  1315,   327,   332,
    2501,   332,   405,  2838,  2630,  2819,   910,   332,   328,   328,
    2839,  1589,  1589,   333,   453,   328,   910,  3290,   332,   333,
     545,  1620,   328,  1623,   333,   406,   333,  1625,   333,   327,
     333,  3158,  1938,   744,   332,   333,  1631,   333,  3159,  2863,
    1634,   339,  1637,   333,  1640,   335,  1643,  2977,  1646,   407,
    1649,   770,  2864,  2208,   333,   408,  1652,  2616,   328,  3043,
    2821,   328,  1872,  1660,   910,   638,   332,   409,  1663,   366,
     333,   910,  1666,  1667,  1669,  1670,  1671,  1672,  1673,   332,
    3070,  1675,   334,   332,   327,   771,   772,   410,   334,  3397,
    3178,   327,   328,   334,   335,   334,   328,   334,   415,   744,
     335,  1585,   333,   416,   334,   335,   744,   335,  3366,   335,
     332,   335,   334,   367,  1703,   333,   335,   418,   335,   333,
    1708,   332,  1097,   334,   335,   421,  1711,  1712,   639,  1715,
    1716,   368,   369,   422,  1931,   335,  1723,  2847,  3390,   334,
     773,  2853,   332,   455,   640,  3480,   333,   328,   327,   332,
    1585,   335,   730,   731,   328,   424,  3644,   333,  1732,   330,
     425,  3715,   328,  3774,  1735,   340,  1736,  1737,  3035,   328,
     431,   334,  2617,  1721,  1601,   432,  2991,   327,   333,  1747,
    1722,  1750,   332,   335,   334,   333,  1753,   332,   334,  1755,
    1098,  1873,   332,   438,   327,   332,   335,   444,  1758,   327,
     335,   327,  1761,  1613,  1614,  1764,   332,  1766,  1767,  1768,
     447,   328,  3455,   456,   449,   334,   343,   450,   333,   602,
    1780,  1781,  1782,   333,   332,  1932,   334,   335,   333,  1787,
    1788,   333,   459,  1797,  1798,  1799,  1800,   460,   335,   393,
     328,   463,   333,   464,   332,   332,   470,   334,   477,   332,
     478,  2339,  2848,   332,   334,   479,  2854,   328,   480,   335,
     333,   481,   328,   327,   328,  1622,   335,   457,   482,  2340,
     327,   327,   327,  1826,  1827,  1828,  1829,  1830,  1260,  1630,
     333,   333,   332,  3036,  2108,   333,   483,   334,  1842,   333,
    1845,   486,   334,  1850,   413,   489,   332,   334,   490,   335,
     334,   327,   493,   494,   335,   495,   332,   332,  2758,   335,
     496,   334,   335,   327,   497,   426,   498,  1874,   333,  1875,
     499,  1877,   445,   335,  1879,   327,   328,   327,   327,   334,
    1885,   350,   333,   328,   328,   328,  1894,  1909,   372,   374,
     384,   335,   333,   333,   500,  1904,  1905,  1906,  1907,   334,
     334,   502,   327,   503,   334,   896,   504,   505,   334,   506,
    1066,   335,   335,   507,   328,  1151,   335,   962,  1309,   386,
     335,   508,  1915,  1941,  1942,  1943,   328,  3141,  3142,  1948,
     509,   397,   510,  1952,   332,   332,   517,   334,   328,   518,
     328,   328,  1961,   332,  1963,  1964,  1965,  1966,  1917,   335,
     522,   334,   523,  1973,   525,  1976,   526,  1977,  1978,   527,
     528,   334,   334,   335,   529,   328,   530,  1396,  1537,  1987,
     333,   333,  1547,   335,   335,  1992,  1569,  1994,  1995,   333,
    1997,  1998,  1999,  2000,  2001,  2002,  2003,  2004,  2005,  2006,
    2007,  2008,  2009,  2010,  2011,  2012,  2013,  2014,  2015,  2016,
    2017,  2018,  2019,  2020,  2021,  2022,  2023,  2024,  2025,  2026,
    2027,  2028,  2029,  2030,  2031,  2032,  2033,  2034,  2035,  2036,
     531,   332,  1466,  1467,  1468,  1469,   332,   332,  2042,  1601,
    1659,  1477,   865,   868,   872,   532,  2265,   533,   534,   334,
     334,   536,   537,  2046,  2047,  2048,   327,   332,   334,   538,
     539,   335,   335,   540,   332,   332,   542,   333,   332,   616,
     335,   648,   333,   333,   332,   332,   679,   683,   332,   332,
    1889,  2055,  1891,   327,   684,  1895,   879,   880,  3260,   881,
     882,   883,   884,   333,  1529,   332,   885,   894,   332,  1911,
     333,   333,  2062,   900,   333,   895,   898,   907,   332,   912,
     333,   333,  2068,  2069,   333,   333,   332,  1668,  1746,   328,
     913,   923,   924,   925,   400,   926,  1825,  2076,  1527,  2077,
     927,   333,  2078,   934,   333,   938,   334,  1565,   950,  2081,
     332,   334,   334,   951,   333,   954,   328,  2089,   335,  1005,
    2091,   419,   333,   335,   335,  1586,   327,   332,  2094,   957,
     960,  1589,   334,  2097,  2099,  2100,  2101,   967,  2102,   334,
     334,   327,  1602,   334,   335,   971,   333,   972,   973,   334,
     334,   335,   335,   334,   334,   335,   332,  1589,  1589,   975,
     327,   335,   335,   333,  2122,   335,   335,  2125,   327,   332,
     334,  2128,  3361,   334,  1849,  2131,   976,   327,  2134,  1876,
    1884,  2137,   335,   334,   332,   335,   977,  2144,   981,   328,
    2147,   334,   333,  2150,   429,   335,   982,   327,  2153,   985,
    1986,   986,   989,   335,   328,   333,   991,  2067,  2098,   433,
    2164,  2243,   992,   332,   994,   334,   997,  2267,  2270,   998,
     333,  2275,  2503,   328,   999,  1006,  1008,   335,   439,  1011,
     332,   328,   334,   332,  1012,  1017,   441,  1018,  2310,  1021,
     328,  2381,  1023,  1025,   335,   465,  1029,  2182,   332,   333,
    1027,  2384,  2183,  2184,  1038,  2186,  2187,  1031,  2189,  2394,
     328,   334,   332,   332,  1033,   471,   333,  1039,  1046,   333,
    1035,   332,  1055,   335,   334,  1071,  1058,  3452,   332,  2204,
     332,  1068,  1078,  2402,   333,   332,   335,  2209,  1072,   334,
    1073,  2210,  1079,  2211,  2212,  1074,  1081,  2213,   333,   333,
    2404,   335,  2215,  1075,  1076,  1084,  2221,   333,  1077,   332,
    2223,  2224,  2225,  1080,   333,  1086,   333,  1087,   334,  1088,
     332,   333,   332,   332,  2236,  2237,  2238,  1089,   332,  2412,
     335,  2242,  2244,  2245,  2246,   334,  2247,   327,   334,  2254,
    2255,  2256,  2475,  2259,  2260,   333,  1114,   335,  1127,  2268,
     335,  2271,   327,   334,   327,   327,   333,  2487,   333,   333,
     332,  1128,   332,  1129,   333,   335,  1130,   334,   334,  1156,
    2276,  2277,  2278,  2279,  1160,   332,   334,  2284,  1161,   335,
     335,  1162,  1176,   334,  1163,   334,  2539,  2293,   335,   327,
     334,   327,  2298,  1166,  2299,   335,   333,   335,   333,  1168,
     328,  1178,   335,  2550,  1187,   473,  2584,   327,   327,  1211,
     332,   333,  1201,   332,   334,   328,  1844,   328,   328,  1269,
     487,  2588,   513,   515,   332,   334,   335,   334,   334,  2311,
    1214,   327,  1283,   334,   327,  2600,  2602,   335,  1233,   335,
     335,   327,  1234,  1241,  2604,   335,   333,  1633,   327,   333,
    1243,  2610,   328,  2614,   328,  1244,  2333,   887,  2659,   921,
     333,  1246,  1247,  1248,   332,   334,  1249,   334,  2348,   327,
     328,   328,  1282,  1251,  1252,   983,   987,   335,  2632,   335,
     334,  1254,  2666,   332,   332,   332,  2363,  2364,   332,  1636,
     332,  1256,   335,  2695,   328,  2759,  2764,   328,  1265,  1000,
     333,  2792,  1002,  1266,   328,  2382,   332,  2385,  1270,  1044,
    1272,   328,  1290,  1413,  1414,   334,  1082,  2395,   334,   333,
     333,   333,  2400,  1415,   333,  1298,   333,   335,  2403,   334,
     335,  2405,   328,  2802,   332,  2804,  2408,  1146,  1292,  1392,
    2413,   335,   333,  2416,  2417,  2418,  2421,  2424,  2425,  2426,
    2427,  2428,  2429,  2430,  2431,  2432,  2433,  2434,  2435,  2438,
    2439,  2440,  2441,  2442,  2443,  2444,  2445, -1425,   327,   334,
     333,  2450,  2451,  2452,  2453,  2454,  2455,  2456,  2457,  2458,
    2459,   335,  1403,  2807,  1404,   327,  2824,  1408,   334,   334,
     334,  2463,  1412,   334,  1460,   334,  1464,  2840,  1478,   332,
     335,   335,   335, -1425,  1472,   335,  2471,   335, -1425,  1463,
    1482,   334,  2476,  1484,  1498,  1476,  1522,  1485,  1523,   332,
    2479,  2480,  1500,   335,   332,  2481,  1501,   327,  1642,   332,
     327,   328,  2486,  2488,  2490,   333,  1157,  2492,  2493,   334,
     327,  2494,   327,  1488,  1489,   332,   332,   332,   328,  2499,
     327,   335,   332,  1173,   327,   333,  2851,  2927,  2994,  1509,
     333,  3010,  1512,  3084,  2517,   333,  2519,  2520,  1490,  1491,
     332,  2524, -1425,  2526,   327,   327,  2529,  1492,   332,  2532,
    1493,   333,   333,   333, -1425,  2413,  2540,  1494,   333,  2543,
     328,  1495,  2546,   328,  1505,  1179,  1506,  2551,  1182,  2552,
    1508,  2553,   327,   328,   334,   328,   333,  3086,  1204,   332,
    1209,   327,  2770,   328,   333,  1510,   335,   328,  1296,  1275,
    1524,  1513,  1304,   332,   334,  1511,  2567,  2568,  2569,   334,
    2570,  2571,   332,  1525,   334, -1425,   335,   328,   328,  1416,
    1417,   335,  1306,  1409,  1514,   333,   335,  1515,  2580,  1516,
     334,   334,   334,  1517,  1518,  2585,  1526,   334,  2589,   333,
    1530,  2592,   335,   335,   335,   328,  1531,  2596,   333,   335,
    1542,  1533,  3100,   327,   328,   334,   327,  2601,  2603,  1554,
    2605,  2606,  2607,   744,  1541,  2609,  2611,   335,  2612,  2613,
    2615,  2618,  3104,  2620,  2621,   335,  2622,  3140,  2624,  2625,
    2626,  2627,  3193,  2628,  2629,  1550,  1551,  1552,  1561,  1562,
    2633,  1568,  2634,  1577,   334,  2635,  1578,  1593,  2637,  3196,
    3199,  2640,  2641,  2642,  1604,  3201,   335,  2644,   334,  1605,
    2647,   332,   327,  1608,  2653,   327,   328,   334,  2657,   328,
     335,  1566,  2660,  3213,  1748,  2865,   328,  1609,  1612,   335,
    2667,  1418,  1419,  1420,  1421,  1422,   327,  1423,   327,  1424,
    1425,  1426,  1427,  1428,  1429,  1430,  1431,   333,  1432,  1433,
    1434,  1435,  1436,  1437,  1438,  1439,  1440,  1441,  1442,  1443,
    1444,  1445,  3228,  1446,  1447,  1448,  1449,  1450,  1451,  1452,
    1453,  1454,  1455,  1615,  1616,   328,  3230,  1645,   328,  1617,
    1839,  1628,  1650,  1857,  1651,  3258,  1648,  1657,   332,   332,
     332,  2910,  1674,  1677,   332,  1683,  1684,  1685,  1686,   328,
    1687,   328,  1688,   332,  1859,  1689,  1865,   332,  1690,  1691,
     332,  1692,  1693,  1694,  1695,  1696,   334,  1697,  1700,  1704,
    1705,  1706,  1707,   327,   333,   333,   333,  2696,   335,  2697,
     333,  1713,  2699,  1719,  1731,   332,  1733,  1734,  1738,   333,
    2704,  2705,  2706,   333,  2707,  2708,   333,  2709,  2710,  2711,
    2712,  2713,  2714,  2715,  2716,  2717,  2718,  2719,  2720,  2721,
    1739,  2722,  2723,  2724,  2725,  2726,  2727,  2728,  2729,  2730,
    1740,   333,  1744,  1745,  2735,  2736,  2737,  2738,  2739,  2740,
    2741,  2742,  2743,  2744,  3259,   332,   328,  2746,  2747,  1751,
    1754,  1927,  1756,   334,   334,   334,   332,  1757,  2753,   334,
    1763,  1765,   332,  2756,  1776,   335,   335,   335,   334,  1777,
    2760,   335,   334,  2762,  1779,   334,  2765,  1784,  1812,  1813,
     335,   333,  1814,   332,   335,   327,  1815,   335,  1816,  1817,
    1818,  1819,   333,  2773,  1820,  2774,   327,  1821,   333,  1822,
     334,  2778,  2779,  2780,  2781,  2782,   332,  2783,   327,   327,
    2785,  2786,   335,  2787,  2788,  1890,  2789,  2790,  3045,   333,
    2413,  3267,  3288,  3305,  2795,  2796,  1823,  2797,  2798,  1824,
    2799,  2800,  1835,  1837,  1838,   327,  2803,  2805,  1851,  1852,
    3307,  1853,   333,  3315,  1854,  2808,  2809,  1855,   328,  2812,
     334,  2813,  2814,  1968,  2816,  2817,  1856,  1863,   332,   328,
    1864,   334,   335,   332,  1979,   332,  1867,   334,  3321,  2825,
    1868,   328,   328,   335,  1869,  1870,  1983,  1988,  2830,   335,
     327,  1871,   327,   327,   327,  1878,   327,  1880,   334,  2834,
    2835,  2836,   332,  1896,   333,  2841,  2842,  2843,   328,   333,
     335,   333,  2849,  2037,  2850,  2852,  2855,  1881,  2857,  2858,
    2859,   334,  2861,  2862,  1882,  1883,  1897,  2866,  3323,  1898,
    1901,  2869,  1902,   335,  2872,  2873,  1908,  1921,   333,  3329,
    1929,  1945,  1959,  2879,  2880,  3332,   332,   332,   332,  1971,
     332,   332,   332,   328,  2886,   328,   328,   328,  2084,   328,
    2086,  2115,  2162,  1972,  2172,  1981,  3349,  1982,  1985,  1990,
    1991,  1993,  2039,   334,   332,  2040,  2043,  2049,   334,  2044,
     334,   327,   333,   333,   333,   335,   333,   333,   333,  3351,
     335,   327,   335,  2052,  2056,  2912,  2061,  2914,  2915,  2916,
    2917,  2918,  2921,  2922,  2925,  2926,  2928,   334,  2930,  2931,
     333,  2934,  2935,  2936,  2937,  2938,  2939,  2942,  2943,   335,
    2945,  2946,  2947,  2948,  2949,  2070,  2071,  2952,  2953,  2954,
    2955,  2956,  2957,  2958,  2959,  2960,  2961,  2962,  2963,  2964,
    2965,  3363,  2966,  2072,   328,  2968,  3400,  2073,  3403,  2226,
    2074,   334,   334,   334,   328,   334,   334,   334,  2079,  2228,
    2083,  2088,  2095,   335,   335,   335,  2103,   335,   335,   335,
     327,   327,  2979,  2980,  2981,  3412,  2983,   332,  2109,   334,
    2986,  2112,  2988,  2113,  2990,  2114,  2992,  2119,  2413,  2120,
    2996,   335,  2998,  2123,  3000,   327,   332,   327,  2126,  2129,
    2200,   327,  2132,  2135,  2142,  3005,  2145,  3007,  3008,  2148,
    3009,  3011,   327,   333,  2158,   327,   327,  2203,  2232,  3418,
    3420,  3426,  3017,  3441,  3443,  3451,  2160,  2288,  3019,  3020,
    3021,  3022,   333,   328,   328,  3024,  3025,  3026,  2229,  2230,
    3028,  3291,  2165,  2171,  3032,  2174,  2175,  3458,  2176,  3037,
    2177,  3038,  3039,  3040,   332,  3041,  3042,  2178,   328,  2181,
     328,  2185,  2188,  2231,   328,  2261,  3049,  3050,  2190,  2263,
    2191,  2192,  3051,  2195,  2205,   328,  2219,  2222,   328,   328,
    2269,   332,   334,  2303,  2305,  2272,  2282,   332,  2289,  2297,
     333,  2302,  2309,  2314,   335,  2317,  2318,   327,   327,  2321,
    2322,   334,  2323,  3069,   327,  3071,  3072,  3073,  3074,  3075,
    3076,  3077,  3078,   335,  3079,  3080,  3081,   333,  3082,  3083,
    3085,   327,  3087,   333,  3088,  3089,  2324,  3090,  3091,  3092,
    3093,  2325,  3096,  3097,  2326,  3098,  3099,  3101,  2327,  3102,
    3103,  3105,  3106,  3107,  3108,  3109,  3110,  3111,  3114,  3117,
    3477,  3120,  3121,  3122,  3123,  3124,  3125,  3126,   332,   334,
     328,   328,  3398,  3129,  2328,  2307,  2319,   328,  2329,  3483,
    2330,   335,  2398,  3133,  3134,  2331,  2332,  2461,  2334,   327,
    2335,  2342,  2343,   327,   328,  2413,   334,   327,  3143,  2406,
     327,  2344,   334,   327,   333,  2345,   327,  2346,   335,  2349,
    2350,  3150,  3151,  3152,   335,  3153,  2351,  2352,  2504,   332,
    2505,  2506,  2507,  2508,  2509,  2510,  2511,  2512,  3160,  3161,
    3164,  2353,  3165,  2354,  2355,   327,  3168,  3485,  2356,  2357,
    2358,  2359,  3173,  3174,  3175,  3176,  3177,  2360,  2361,  2362,
    3179,  3180,   328,  3181,  3182,   333,   328,  2554,  2366,  2367,
     328,  2638,  2368,   328,  3487,  2663,   328,  2369,  2668,   328,
    3503,  2693,   327,   334,  2870,  3194,  3195,  3197,  3198,  3200,
    3202,  3203,  2370,  3205,  2371,   335,  3209,  2372,  2373,  3212,
     327,  3214,  3215,  3216,  3217,  3218,  3219,  3220,   328,  3221,
    3222,  2374,  3224,  2887,  2375,  3227,  3229,  3231,   327,  3232,
    3233,  3234,  3235,  3236,  3237,  3238,  2376,  3239,  3240,   327,
    3241,  3242,  2377,  3243,   334,  3246,  3247,  3248,  3249,  3250,
    3251,  3252,   327,  2378,  2379,   328,   335,  3255,  3256,  2380,
    2889,  3505,   327,  2383,  3257,  2413,  2413,  2386,  2387,   327,
    2388,  2389,  3262,   328,   332,  3265,  3266,  3268,  2891,  2390,
     332,  2391,  2396,  2397,  3270,  3271,   327,  3272,  3275,  3276,
    2401,   328,  3281,  2409,  2460,  2465,  2893,  3284,  3285,  3286,
    3287,  3289,   328,  2468,  3292,  3293,  3294,  2904,  2477,  2478,
     333,  2484,  3525,   327,   327,   328,   333,   332,  3301,  3302,
    2973,  3303,  3304,  2485,  3306,   328,  3308,  3309,  2491,  3310,
    3057,  3311,   328,  3312,  2495,  3313,  2496,  3059,  3316,  3317,
    3318,  3319,  3320,  3322,  3324,  3325,  3326,  2497,  3327,   328,
    3328,  3330,  2498,   333,  3186,  3333,   327,  3335,  3336,  3337,
    3338,  3339,  3340,  3341,  2502,  3344,  2115,  3347,   327,  3348,
    3350,  3352,  3353,  3354,  3355,  3356,   328,   328,  2516,   334,
    3360,  3188,  3296,  3362,  2413,   334,  3365,  2518,  2521,  3368,
    3369,   335,  3370,  2522,  3371,  3372,  2525,   335,  3374,  2527,
    3377,  2528,   332,  2530,  3379,  3382,  2531,  2533,  2534,  3386,
    3387,  3388,  2535,  3389,  2541,   332,  3391,  3392,  3393,   328,
    2542,   327,   334,  2544,  3298,  2545,  3401,  2561,  3404,  2547,
    3405,   328,  3406,  3407,   335,  3409,  3394,  3411,   333,  2548,
    3413,  3414,  3415,  2549,  3417,  2555,  3419,  3526,  3421,  3422,
    3423,   333,  3425,  3529,  3427,   327,  2556,  3428,   327,  3429,
    3430,  3431,  3432,  3433,  3434,  3435,  2557,  3436,  3437,  2559,
    3438,  3439,  3440,  2560,  3442,  2562,  3444,  3445,  3446,  3447,
    3448,  2572,  2574,  2575,   328,  2413,   332,  2581,   332,  3491,
    3541,  2582,  3456,  2586,  3459,  3460,  3461,  2590,  2593,  2595,
    3463,  2597,  2598,  2599,   332,  3466,   332,   334,  2619,  2636,
    3471,  3472,  3473,  3474,   327,  3475,  3476,  3478,   328,   335,
     334,   328,   333,  3496,   333,  3481,  3498,  2643,  3482,  3484,
    3486,  3488,   335,  3489,  2652,  3490,  2654,  2658,  3493,  3494,
     333,  3495,   333,  2661,  2665,   332,  2670,  3501,  2671,  3502,
    2675,  3504,  3506,  3507,  3508,  3509,  3510,  3511,  3512,  2676,
    3514,   327,  3516,  2683,  3518,  2684,  3519,  2685,  3520,  3521,
    3522,  3523,  3524,  2686,  2687,  3546,  2413,   328,   332,  2688,
    3527,   333,  3554,  2689,  3530,  2690,  2691,  2698,  3548,   327,
    2749,   334,  3534,   334,  2754,  3536,  3537,   327,  3539,  3540,
    3542,  2755,  3543,   335,  2757,   335,  2761,  2763,  3547,   334,
    3549,   334,  3551,  2766,   333,  2767,  2768,  3558,  3561,  2771,
    2784,   335,   332,   335,   328,  3566,   327,  2777,  3569,  3598,
    3571,  3572,  3573,  3574,  3575,  3576,  3577,   332,  3578,   327,
    3579,  2791,  2793,  3582,  3584,  3585,  3586,  3587,  3588,  2801,
     334,  3589,   328,  2806,  3591,  2810,  2815,  3600,   333,  3550,
     328,  3568,   335,  3595,  2826,  3614,  3596,  3597,  2827,   332,
    2829,  2831,  2832,   333,  2833,  3603,  2844,  3570,  2856,  3581,
    2860,  3605,  3606,   334,  3607,  3608,  2867,  3609,  2868,   328,
    2874,  2875,  2876,  3613,  3671,   335,  3617,  3619,  3620,  3621,
    3622,  3623,   328,  2881,  2882,   333,  3627,  3673,  3629,  3632,
    3635,  2883,  3638,  3640,  2897,  2899,  2900,  2902,  3583,  3645,
    3646,  2906,  2929,  2944,  2970,  3130,  2971,   334,  2975,  3652,
    2976,  3655,  2978,  3658,  2982,  3659,  2984,  3661,  2985,   335,
    2987,  3663,   334,  3664,  3666,  3668,  3669,  3670,  2989,  2993,
    2995,  3590,  2997,  2999,   335,  3675,  3676,  3001,  3677,  3678,
     332,  3679,  3002,  3006,  3012,  3015,   332,  3684,  3016,  3685,
    3018,  3023,  3027,  3033,   334,  3048,  3688,  3052,  3689,  3690,
    3055,  3691,  3692,  3693,   332,  3695,   335,  3056,  3061,  3062,
    3698,   332,  3699,  3700,  3701,  3602,   333,  3063,  3064,  3704,
     332,  3707,   333,  3710,  3065,  3711,  3066,  3067,  3128,  3716,
    3612,  3717,  3132,  3719,  3135,  3721,  3136,  3723,  3137,  3725,
     333,  3138,  3139,  3144,  3729,  3731,  3145,   333,  3732,  3146,
    3733,  3734,  3149,  3735,  3736,  3737,   333,  3154,  3738,  3155,
    3739,  3740,  3616,  3741,   332,  3742,  3156,  3157,  3169,  3745,
     332,   332,  3185,  3746,   332,  3747,  3190,  3749,  3204,  3751,
    3208,  3753,  3754,  3755,  3756,   334,  3223,  3253,  3254,  3760,
     332,   334,  3269,  3763,  3277,  3764,  3295,   335,   332,  3767,
     333,  3314,  3331,   335,  3770,  3334,   333,   333,  3357,   334,
     333,  3775,  3776,  3359,  3778,   332,   334,  3779,  3358,  3780,
    3367,   335,  3783,  3784,  3373,   334,   333,  3787,   335,  3789,
    3383,  3791,  3396,  3792,   333,   329,   331,   335,  3399,  3796,
    3797,  3798,   344,   345,  3402,  3408,   348,  3410,   351,  3416,
     332,   333,  3424,  3618,  3449,  3450,  3453,  3454,  3457,  3626,
     332,   332,  3462,   373,   375,  3464,  3465,  3470,  3479,   334,
    3500,   382,   383,   385,   387,   334,   334,  3628,   332,   334,
    3513,   335,  3515,   398,  3639,   401,   333,   335,   335,  3517,
    3528,   335,  3531,  3660,  3532,   334,   333,   333,  3533,   417,
    3535,   420,  3538,   334,   423,  3544,  3545,   335,  3552,   430,
    3553,  3564,   434,   332,   333,   335,   440,   442,   443,  3565,
     334,   332,   448,  3567,  3580,   454,   458,  3592,  3594,   461,
     462,  3624,   335,   466,  3625,   468,   469,  3665,   472,   474,
    3641,  3647,   544,  3667,  3694,  3648,  3649,  3724,   484,   333,
    3650,   488,  3662,  3682,  3686,   334,  3696,   333,  3697,  3702,
    3703,  3712,  3713,  3728,  3718,   334,   334,   335,  3720,  3722,
    3726,  3730,  3727,  3743,  3748,   514,   516,   335,   335,   519,
     520,   521,  3750,   334,  3752,  3757,  3758,  3761,  3744,  3762,
    3765,  3766,  3771,  3772,  3781,   335,  3788,  3793,  3799,  1107,
    1575,  1159,   904,   776,   777,  1110,   902,  1661,   546,  1658,
     778,  1626,  1004,   580,   914,  1579,   952,   587,   334,   599,
     601,   603,   945,  3759,  1596,   974,   334,  1050,  1019,   608,
     335,     0,  1624,  3769,  3777,  1803,  2080,   869,   335,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  3782,     0,     0,     0,     0,     0,   651,     0,   654,
     655,     0,     0,   657,   659,   661,   662,     0,   664,     0,
     666,   668,   670,   672,     0,     0,   664,   668,   672,     0,
       0,   678,     0,   680,   681,   682,  3786,     0,     0,   687,
       0,     0,     0,     0,  3790,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   690,   692,   693,     0,     0,   694,
     695,   697,   698,   701,     0,   729,   599,   599,   733,   734,
     735,   738,     0,     0,   745,   747,   749,     0,     0,   750,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   796,   798,   799,     0,     0,     0,     0,
       0,     0,     0,     0,   878,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     888,     0,     0,   892,   892,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   901,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   666,   664,
     918,     0,     0,   922,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   932,   932,     0,
     937,     0,     0,     0,   668,     0,     0,   672,   670,     0,
       0,   948,   878,     0,     0,     0,   918,   878,     0,     0,
     956,     0,     0,   959,     0,     0,     0,   965,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   690,
       0,     0,     0,     0,     0,     0,   729,     0,   980,     0,
       0,     0,     0,   984,     0,     0,     0,   988,     0,     0,
       0,   990,     0,     0,   993,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1001,  1003,   678,     0,     0,
       0,     0,     0,     0,   878,     0,     0,     0,     0,     0,
    1016,     0,     0,     0,     0,   692,     0,     0,     0,     0,
       0,     0,     0,     0,   878,     0,     0,     0,     0,     0,
       0,     0,     0,  1037,     0,     0,     0,     0,     0,     0,
       0,  1043,  1045,     0,     0,     0,     0,     0,  1049,   697,
    1053,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1057,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1083,     0,  1085,     0,     0,     0,
       0,     0,  1091,  1092,  1093,  1094,     0,     0,     0,  1100,
    1101,  1102,  1103,  1104,     0,  1106,  1106,  1109,  1109,     0,
       0,     0,     0,     0,  1115,     0,  1118,  1120,     0,  1122,
    1124,  1126,     0,     0,     0,     0,     0,  1145,  1147,  1148,
       0,     0,     0,     0,     0,  1158,     0,     0,     0,     0,
       0,  1165,     0,     0,     0,  1171,     0,     0,     0,     0,
    1174,     0,     0,     0,     0,     0,  1177,     0,  1180,  1181,
    1183,     0,     0,  1120,     0,  1122,     0,  1188,     0,  1191,
       0,  1193,     0,     0,     0,  1196,     0,  1197,     0,  1198,
       0,  1199,     0,  1200,     0,     0,     0,     0,     0,     0,
    1205,     0,     0,     0,  1210,     0,  1120,     0,  1122,     0,
       0,     0,     0,  1223,  1224,  1225,  1226,     0,  1227,  1228,
       0,     0,  1232,     0,     0,  1235,     0,     0,     0,     0,
       0,     0,  1242,     0,     0,  1245,     0,     0,     0,     0,
       0,     0,     0,  1253,     0,     0,     0,     0,  1258,     0,
       0,     0,     0,     0,     0,  1262,  1264,     0,     0,     0,
    1268,     0,     0,     0,     0,     0,  1271,     0,  1273,  1274,
    1276,     0,     0,     0,     0,     0,     0,     0,  1287,  1289,
       0,     0,     0,     0,  1295,  1297,     0,  1303,  1305,  1307,
    1308,     0,  1311,     0,  1316,     0,     0,     0,     0,     0,
       0,     0,     0,  1394,  1395,     0,     0,  1401,  1402,     0,
       0,     0,  1406,  1407,     0,  1410,     0,  1411,     0,     0,
    1456,  1457,  1458,  1459,     0,  1461,  1462,     0,     0,     0,
       0,     0,     0,     0,     0,  1471,     0,  1475,     0,     0,
       0,  1479,  1480,     0,  1481,     0,     0,     0,     0,     0,
    1487,     0,     0,     0,     0,     0,     0,     0,     0,  1496,
    1497,     0,     0,     0,     0,  1502,  1503,     0,     0,     0,
    1507,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1519,  1520,  1521,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1532,     0,  1535,     0,  1122,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1543,     0,   892,     0,   892,   892,     0,     0,     0,     0,
    1549,     0,     0,     0,     0,  1147,     0,     0,     0,     0,
       0,  1558,     0,     0,     0,     0,     0,   918,     0,     0,
       0,     0,     0,     0,     0,     0,  1567,     0,     0,     0,
       0,   932,     0,  1573,   932,   932,   937,     0,  1576,     0,
       0,     0,     0,     0,     0,     0,     0,  1581,  1583,     0,
    1535,     0,     0,   918,  1535,     0,     0,     0,     0,     0,
       0,     0,     0,  1595,  1598,  1600,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1607,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1610,     0,     0,  1611,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1618,     0,  1621,
       0,  1535,     0,     0,     0,  1016,     0,  1627,     0,     0,
       0,     0,  1629,     0,     0,  1535,  1632,     0,  1635,     0,
    1638,     0,  1641,     0,  1644,     0,  1647,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1049,     0,
       0,     0,  1053,     0,  1662,     0,     0,  1664,  1665,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1676,
       0,  1678,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1698,  1699,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1709,     0,  1710,     0,     0,     0,  1714,     0,     0,     0,
       0,  1122,  1720,     0,  1724,  1725,  1726,  1727,     0,     0,
       0,     0,     0,     0,  1728,     0,  1729,     0,  1730,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1749,     0,     0,
    1752,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1759,  1760,     0,
    1762,     0,     0,     0,     0,     0,     0,  1770,  1770,  1770,
    1770,  1770,  1775,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1786,     0,     0,     0,  1795,  1796,
       0,     0,     0,     0,  1802,  1802,  1804,  1805,  1770,  1807,
       0,  1808,     0,  1811,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1834,     0,     0,     0,  1836,
       0,     0,     0,     0,  1840,     0,  1843,     0,     0,  1848,
       0,     0,     0,     0,     0,     0,     0,  1858,     0,  1860,
       0,  1862,     0,     0,     0,  1866,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1900,     0,
       0,  1903,     0,     0,     0,     0,     0,  1910,     0,  1912,
    1913,  1914,  1916,  1918,  1920,     0,  1922,  1923,  1924,  1925,
    1926,  1928,     0,  1930,  1933,     0,  1936,  1937,  1939,  1940,
       0,     0,     0,  1944,     0,  1947,     0,  1949,  1950,  1951,
       0,  1953,  1954,  1955,  1956,  1957,  1958,     0,  1960,     0,
    1962,     0,     0,     0,     0,  1967,  1969,  1970,     0,     0,
       0,     0,     0,     0,     0,     0,  1980,     0,     0,     0,
       0,  1984,     0,     0,     0,     0,     0,  1989,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2038,     0,     0,
       0,     0,     0,  2041,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2045,   826,     0,     0,     0,  2050,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   827,     0,     0,
       0,  2051,     0,     0,  2053,     0,     0,     0,  2054,     0,
       0,     0,     0,  2057,  2058,     0,     0,     0,     0,     0,
       0,  2059,  2060,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  2064,     0,     0,     0,  1122,  2066,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   892,     0,     0,
       0,  2075,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  2085,
       0,  2087,     0,     0,     0,   932,     0,  2090,     0,     0,
    1583,     0,  2093,     0,  2093,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1600,     0,  2105,
       0,  2105,     0,     0,     0,   828,   829,     0,     0,     0,
       0,  2117,  2118,     0,     0,     0,     0,     0,     0,     0,
    2121,     0,     0,  2124,     0,     0,     0,     0,  2127,   830,
       0,  2130,     0,     0,  2133,     0,     0,  2136,     0,     0,
    2140,     0,     0,  2143,     0,     0,  2146,     0,     0,  2149,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1770,     0,  2161,  2163,     0,     0,     0,
    2166,  2167,  2168,  2169,  2170,     0,     0,  2173,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  2180,     0,     0,     0,   831,   832,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2194,     0,   833,     0,  2199,     0,  2202,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2214,     0,     0,     0,
    2217,  2218,     0,     0,     0,     0,     0,     0,     0,     0,
    2227,     0,  2227,  2227,  2227,  2227,     0,  2233,  2234,     0,
       0,     0,     0,     0,     0,     0,     0,  2241,     0,   834,
     835,   836,   837,     0,  2248,     0,     0,     0,     0,     0,
       0,     0,  2262,     0,  2262,  2266,     0,  2227,     0,     0,
       0,  2273,  2274,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2287,     0,     0,     0,     0,
       0,     0,  2290,     0,     0,  2294,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  2300,  2301,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2304,  2306,     0,  2308,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  2316,
       0,   838,     0,   839,  2320,     0,   840,     0,   841,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   826,     0,   842,  2341,   843,   844,
       0,     0,     0,     0,  2347,     0,     0,   845,     0,   846,
     827,     0,     0,     0,     0,     0,     0,   819,   866,   867,
       0,   847,   848,   849,     0,  2365,     0,   850,   851,   852,
     853,     0,     0,     0,   854,   855,     0,     0,   856,   857,
     858,   859,     0,   860,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   861,   862,   863,     0,     0,  2399,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2407,     0,     0,     0,     0,  2410,  2411,     0,  2414,  2415,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2446,  2447,  2448,  2449,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   828,   829,
       0,     0,     0,     0,     0,     0,  2462,     0,     0,  2464,
       0,     0,     0,     0,  2466,  2467,     0,     0,  2469,  2470,
       0,     0,   830,     0,     0,     0,  2472,  2474,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2489,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2500,     0,
       0,     0,     0,     0,   729,     0,     0,     0,     0,  2515,
       0,     0,     0,     0,     0,     0,     0,     0,  2523,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   831,
     832,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   833,     0,     0,     0,
    2227,     0,     0,     0,     0,     0,     0,  2558,     0,     0,
       0,     0,     0,     0,  2563,     0,     0,     0,     0,     0,
    2564,  2565,  2180,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2573,     0,     0,     0,     0,  2576,
       0,  2577,  2578,     0,  2579,     0,     0,     0,     0,     0,
       0,  2583,     0,     0,  2587,     0,     0,  2591,     0,     0,
       0,  2594,   834,   835,   836,   837,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2608,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2623,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2639,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2651,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2662,  2664,     0,     0,     0,     0,     0,     0,  2669,
       0,     0,     0,     0,     0,     0,  2672,     0,     0,     0,
       0,     0,     0,     0,   838,     0,   839,     0,     0,   840,
       0,   841,     0,     0,     0,     0,     0,     0,  2679,     0,
    2680,     0,     0,     0,     0,     0,     0,     0,     0,   842,
       0,   843,   844,     0,     0,     0,     0,     0,     0,     0,
     845,     0,   846,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   847,   848,   849,     0,     0,     0,
     850,   851,   852,   853,     0,     0,     0,   854,   855,     0,
       0,   856,   857,   858,   859,     0,   860,     0,     0,     0,
    2692,  2694,     0,     0,     0,     0,   861,   862,   863,     0,
       0,  2700,  2701,     0,     0,  2702,  2703,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2731,  2732,  2733,
    2734,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  2745,     0,     0,     0,     0,  2748,     0,
       0,     0,     0,  2751,  2752,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2772,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2811,     0,  2564,     0,     0,
       0,     0,     0,     0,  2818,     0,     0,     0,  2820,     0,
    2822,     0,     0,     0,  2823,     0,     0,     0,     0,     0,
       0,     0,  2828,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  2837,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2871,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  2888,     0,  2890,     0,     0,  2892,     0,
    2894,     0,     0,     0,     0,     0,     0,     0,  2898,     0,
       0,     0,  2901,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2903,     0,  2905,  1317,     0,
       0,  2911,     0,  2913,     0,     0,     0,     0,     0,     0,
       0,  1318,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1319,     0,     0,     0,     0,
    1320,     0,  2950,  2951,     0,     0,     0,     0,     0,  1321,
       0,     0,  1322,  1323,     0,     0,     0,  1324,  1325,     0,
       0,  2967,     0,  2969,     0,     0,     0,  2972,     0,     0,
       0,  2974,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1326,     0,     0,     0,     0,     0,  3003,
    3004,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    3013,     0,  3014,     0,  1327,  1328,  1329,  1330,     0,     0,
       0,     0,     0,  1331,     0,     0,     0,  1332,     0,     0,
    1333,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1334,     0,  1335,     0,     0,  3046,     0,     0,
    3047,     0,     0,     0,     0,     0,     0,     0,  1336,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  3058,
       0,  3060,  1337,  1338,     0,     0,     0,     0,     0,     0,
       0,  1339,     0,  1340,     0,     0,  1341,     0,  1342,  3068,
    1343,  1344,     0,     0,     0,  1345,     0,  1346,  1347,  1348,
    1349,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1350,     0,  1351,     0,     0,
       0,     0,     0,  1352,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  3127,     0,
       0,     0,  1353,  3131,     0,     0,     0,  1354,  1355,  1356,
       0,  1357,     0,  1358,     0,  1359,     0,     0,     0,     0,
    1360,  1361,  1362,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  3147,  3148,     0,     0,     0,     0,
    1363,     0,     0,     0,     0,     0,     0,     0,     0,  1364,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1365,  1366,  1367,  1368,  1369,  1370,  1371,     0,  1372,
    1373,  1374,  1375,     0,     0,     0,     0,     0,     0,  1376,
    1377,     0,  1378,     0,     0,     0,  1379,     0,     0,  3187,
       0,  3189,     0,     0,     0,     0,     0,     0,     0,  1380,
       0,  1381,     0,     0,     0,     0,  1382,     0,     0,     0,
       0,  1383,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1384,     0,     0,     0,  1385,  1386,  1387,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1388,     0,     0,  1389,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  3261,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  3297,     0,
    3299,     0,     0,  3300,     0,     0,     0,     0,  1390,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  3364,     0,  2180,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  3385,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  3395,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2564,  2180,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  3492,     0,     0,     0,     0,     0,
    3497,     0,  3499,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2564,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  3555,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  3593,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  3599,     0,
    3601,     0,     0,     0,     0,     0,  3604,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  3615,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  3642,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  3651,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  3672,     0,
    3674,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  3683,     0,     0,     0,     0,     1,     0,
       2,     3,  3687,     0,     0,     0,     4,     0,     5,     0,
       0,     0,     0,     0,     0,     0,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,     0,     0,
       0,     0,     0,     0,     0,    13,    14,    15,     0,     0,
       0,     0,     0,     0,     0,    16,     0,    17,    18,    19,
       0,    20,     0,     0,     0,     0,     0,     0,     0,    21,
       0,     0,     0,     0,     0,    22,    23,     0,     0,    24,
       0,     0,     0,    25,    26,    27,    28,     0,     0,     0,
      29,     0,    30,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    31,     0,     0,    32,    33,    34,
      35,   543,    36,     0,    37,    38,     0,  3768,     0,     0,
       0,     0,     0,    39,     0,     0,     0,    40,    41,     0,
       0,     0,    42,    43,     0,     0,    44,    45,     0,     0,
       0,  3785,     0,     0,     0,     0,    46,     0,     0,     0,
       0,     0,     0,  3794,    47,  3795,     0,    48,    49,     0,
       0,    50,     0,     0,    51,     0,     0,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,     0,
       0,     0,     0,    64,    65,     0,     0,     0,    66,     0,
      67,     0,    68,    69,    70,    71,     0,    72,     0,    73,
       0,     0,    74,     0,     0,     0,     0,    75,     0,     0,
      76,    77,    78,    79,    80,     0,     0,    81,     0,     0,
      82,    83,    84,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    85,     0,     0,
       0,     0,     0,     0,     0,     0,    86,     0,     0,     0,
       0,     0,    87,     0,     0,     0,    88,     0,     0,     0,
       0,    89,     0,    90,    91,     0,     0,     0,     0,     0,
      92,     0,    93,     0,     0,    94,    95,     0,     0,    96,
       0,     0,    97,    98,    99,     0,   100,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   101,   102,   103,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   104,   105,     0,     0,     0,     0,   106,   107,   108,
       0,   109,     0,     0,     0,     0,     0,   110,     0,   111,
     112,     0,     0,     0,   113,     0,     0,     0,     0,     0,
       0,   114,     0,     0,   115,   116,     0,     0,     0,   117,
       0,     0,   118,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   119,     0,     0,   120,
       0,     0,     0,     0,   121,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   122,
       0,   123,     0,   124,     0,     0,   125,   126,   127,     0,
     128,   129,     0,   130,     0,     0,   131,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   132,     0,     0,     0,
       0,   133,   134,   135,     0,     0,     0,     0,     0,     0,
       0,   136,   137,   138,   139,   140,   141,     0,   142,     0,
       0,     0,     0,     0,   143,   144,     0,     0,     0,     0,
       0,     0,   145,     0,     0,     0,   146,   147,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     148,   149,     0,     0,     0,     0,     0,   150,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     151,     0,   152,   153,   154,   155,   156,   157,   158,   159,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     160,   161,     0,     0,   162,   163,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   164,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     165,   166,   167,     0,     0,   168,   169,     0,     1,     0,
       2,     3,     0,     0,     0,     0,     4,   170,     5,     0,
       0,     0,     0,     0,     0,     0,     6,     0,     7,     8,
       9,    10,    11,     0,     0,     0,    12,     0,     0,     0,
       0,     0,     0,     0,     0,    13,    14,    15,     0,     0,
       0,     0,     0,     0,     0,    16,     0,    17,    18,    19,
       0,    20,     0,     0,     0,     0,     0,     0,     0,    21,
       0,     0,     0,     0,     0,    22,    23,     0,     0,    24,
       0,     0,     0,    25,    26,    27,    28,     0,     0,     0,
      29,     0,    30,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    31,     0,     0,    32,    33,    34,
      35,     0,    36,     0,    37,    38,     0,     0,     0,     0,
       0,     0,     0,    39,     0,     0,     0,    40,    41,     0,
       0,     0,    42,    43,     0,     0,    44,    45,     0,     0,
       0,     0,     0,     0,     0,     0,    46,     0,     0,     0,
       0,     0,     0,     0,    47,     0,     0,    48,    49,     0,
       0,    50,     0,     0,    51,     0,     0,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,     0,
       0,     0,     0,    64,    65,     0,     0,     0,    66,     0,
      67,     0,    68,    69,    70,    71,     0,    72,     0,    73,
       0,     0,    74,     0,     0,     0,     0,    75,     0,     0,
      76,    77,    78,    79,    80,     0,     0,    81,     0,     0,
      82,    83,    84,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    85,     0,     0,
       0,     0,     0,     0,     0,     0,    86,     0,     0,     0,
       0,     0,    87,     0,     0,     0,    88,     0,     0,     0,
       0,    89,     0,    90,    91,     0,     0,     0,     0,     0,
      92,     0,    93,     0,     0,    94,    95,     0,     0,    96,
       0,     0,    97,    98,    99,     0,   100,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   101,   102,   103,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   104,   105,     0,     0,     0,     0,   106,   107,   108,
       0,   109,     0,     0,     0,     0,     0,   110,     0,   111,
     112,     0,     0,     0,   113,     0,     0,     0,     0,     0,
       0,   114,     0,     0,   115,   116,     0,     0,     0,   117,
       0,     0,   118,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   119,     0,     0,   120,
       0,     0,     0,     0,   121,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   122,
       0,   123,     0,   124,     0,     0,   125,   126,   127,     0,
     128,   129,     0,   130,     0,     0,   131,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   132,     0,     0,     0,
       0,   133,   134,   135,     0,     0,     0,     0,     0,     0,
       0,   136,   137,   138,   139,   140,   141,     0,   142,     0,
       0,     0,     0,     0,   143,   144,     0,     0,     0,     0,
       0,     0,   145,     0,     0,     0,   146,   147,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     148,   149,     0,     0,     0,     0,     0,   150,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     151,     0,   152,   153,   154,   155,   156,   157,   158,   159,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     160,   161,     0,     0,   162,   163,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   164,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     165,   166,   167,     0,     0,   168,   169,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   170
};

static const yytype_int16 yycheck[] =
{
       6,     7,     8,     9,   445,  1640,  2181,   891,  1204,   893,
     964,   198,   947,  1209,   930,  1211,  1741,   933,   343,   397,
    2847,    17,    12,   103,    84,    38,  2853,  1197,  1198,  1199,
    1200,    28,   370,    32,    38,    19,  2117,    38,     7,   181,
     249,   228,   175,    49,   181,     4,    76,   113,   249,   249,
      76,   249,    51,  1778,    22,   247,   249,  1227,  1783,    65,
    1785,   249,    38,   249,   249,    76,   249,   249,   246,    82,
      76,   258,   111,   126,   249,   134,   249,     9,    82,   134,
      22,    82,   112,    89,    44,    45,   112,   100,   101,    86,
      76,    88,   111,    83,   249,    76,   100,   101,   249,   100,
     101,   112,   244,   181,    17,    47,    82,   244,   249,   175,
     245,   101,   118,   181,   372,   249,   249,   123,   124,   152,
     108,   181,    10,    76,   100,   101,   112,    15,   181,    75,
     249,   112,   249,    75,   124,   181,   249,   113,   181,   111,
     249,    76,   249,   140,   111,    58,     8,   153,   249,   215,
      76,   181,   537,   538,   181,   181,    18,   147,   164,   112,
     177,    76,   249,   193,   181,   133,   244,   193,    30,   125,
     181,    33,    34,    35,    36,    37,   244,   112,   184,   185,
      68,   332,   193,   249,   244,   191,   112,   916,    76,   249,
     249,   244,   160,   379,   249,   181,   202,   112,   244,    76,
     181,   244,  3029,  3030,   249,   181,   407,   193,  3035,   181,
     249,   384,   193,   203,   406,    77,   249,   244,   160,   249,
      76,   211,   181,   952,   112,   403,   114,   244,   181,   219,
     249,   227,   249,    76,   176,   112,   249,   249,    14,   215,
     193,    76,   249,   378,    20,   249,   181,   432,   249,   181,
     450,   249,   134,   333,    76,   181,   112,   181,   193,   181,
     250,   249,   252,   249,   270,   458,   181,   193,   244,   112,
     260,   239,   244,   249,    76,   458,   249,   112,   193,    76,
     249,   223,   249,   458,   146,   244,   292,   293,   150,   181,
     112,    67,   249,   181,   240,   504,   249,   239,   240,   305,
      76,   499,   500,   501,   181,   193,   248,   440,   445,   450,
     112,   499,   244,   379,   249,   112,   193,   499,   460,   372,
     244,   263,   244,   249,   458,   181,   372,   215,   318,   407,
      76,    46,   338,   440,   276,   341,   112,   193,   181,   458,
    3167,   458,   348,  3170,  3171,   413,   181,   353,   249,   321,
     193,   384,   244,   359,   360,   111,   302,  2082,   193,   181,
     366,   367,   700,   369,   363,   249,   112,   249,   249,   412,
    1566,   193,   372,   249,   440,   381,   249,    92,    93,   181,
     326,   181,   181,   245,   181,   412,   412,   393,   394,   449,
     407,   193,   254,   249,   375,   257,   193,   408,   404,   405,
     262,   263,   264,   379,   266,   181,   181,   269,   231,   125,
     108,   273,   418,   458,   249,   421,   413,   193,   424,   606,
     404,   427,   249,    76,   414,   422,   249,   440,   458,   442,
     420,   373,  2513,   440,  2515,   181,   440,   249,   437,   440,
     446,   383,   218,   385,   244,   244,   458,   193,   444,   125,
     249,  3278,   249,   375,   249,    75,  3283,   319,   249,   112,
     458,   249,   249,   175,   440,   652,   456,   249,   125,   244,
     246,   181,   458,   249,   434,   435,   436,   437,   438,   249,
     256,   343,   488,   249,   249,   458,   348,   249,   412,   495,
     432,   379,  1662,   249,    55,    76,   272,   412,   685,   321,
      76,   458,   181,   249,  2139,   458,   503,   504,   505,   125,
    2235,   508,   509,   362,  2239,   135,   136,   137,   138,   139,
     412,   409,   249,   458,    76,   249,   249,    39,   181,   535,
     590,   112,   458,   249,   244,   412,   112,   249,   381,   249,
     193,   547,   548,   549,   550,   551,   552,   553,   554,    76,
     556,   249,   249,   125,   560,   125,   934,  3384,   334,   335,
     112,   410,   125,   213,   214,   244,   332,   440,   249,   894,
      82,   249,   249,   249,   580,   181,   378,   458,   378,   585,
     586,   587,   458,  1779,    38,   112,   125,   181,  1784,   365,
     366,   249,   249,   599,   249,   601,   249,   379,   604,   125,
     181,   249,   458,   402,    76,   181,   181,   119,   407,   615,
     372,   377,   193,   619,   620,   621,   622,   193,   268,   625,
     181,    76,   484,   458,    76,   249,   632,   181,    82,   181,
     181,   458,   638,   249,   640,   497,   498,   412,   244,   249,
     112,   193,  1596,   125,  1579,   651,   100,   101,   372,   161,
     244,   657,   181,   659,   181,   661,   662,   112,  1574,   249,
     112,   458,  1546,   458,   174,   181,   193,   317,   249,   244,
     458,   458,   184,   249,   680,   681,   682,   249,   249,   249,
     456,   249,   302,   244,   690,   305,   249,   249,   249,   876,
     244,   406,   402,   244,   372,   372,   375,   249,    38,   705,
     706,   707,   708,   709,   181,   315,   326,    38,   181,   181,
     249,   249,   458,   719,   372,   244,   249,   181,   724,   181,
     726,   193,   249,   729,   372,   249,   181,   733,   244,   181,
      38,   249,   738,   181,   583,   584,   585,   586,   193,   249,
     390,   193,    82,   749,   750,   751,   752,   397,   372,   755,
     756,    82,   181,   249,   760,    38,   762,  2482,  2483,    39,
     100,   101,   949,   249,   770,   415,   953,   244,    76,   100,
     101,   244,   378,   779,    82,   781,   372,   249,   784,   785,
     244,   375,   244,    63,   790,    49,    38,   338,   503,    69,
     505,   506,   100,   101,   249,   249,   244,   249,   513,    82,
     125,   242,    82,   378,   112,   458,   249,   813,   249,    89,
     525,   526,   527,   528,   529,   244,   531,   100,   101,   825,
     372,    38,   828,  1010,   378,   175,   541,   542,   448,   544,
      82,   546,   547,   345,   346,   841,   551,  1024,   249,   119,
     846,  1587,   557,   558,   559,   560,   375,   249,   100,   101,
     181,   249,   379,   230,   231,   298,   181,   249,    38,    76,
     231,   338,   378,   869,   870,    82,   411,  1613,  1614,   125,
     125,   416,   249,   181,   880,   881,   125,   458,   249,  1217,
     886,   161,   458,   100,   101,   193,   892,  2083,  2084,   242,
     338,   897,   125,   181,    62,   112,   249,   903,   181,   249,
     906,   907,    82,  2538,   184,   378,   175,   379,    76,   249,
     175,   249,   918,   244,   378,   157,   378,   923,   249,   244,
     100,   101,   928,   929,   379,   181,   181,   379,   371,   579,
     580,   581,   181,   249,   198,   199,   200,   201,   483,   181,
     195,   249,   948,    62,   112,   175,   375,   249,   954,   955,
     956,   157,   958,   175,   234,   961,   244,    76,   503,   965,
     966,   244,    75,    76,   181,   175,   249,    75,    76,    75,
      76,   249,    76,   979,    76,   181,   193,   420,   181,    76,
     249,    76,   249,   242,   249,   157,   175,    76,   244,   244,
     249,   997,   998,   112,   249,   244,   175,   249,    76,   112,
     157,  1007,   244,  1009,   112,   249,   112,  1013,   112,   181,
     112,   242,   300,   181,    76,   112,  1022,   112,   249,   249,
    1026,   125,  1028,   112,  1030,   193,  1032,   249,  1034,   249,
    1036,   311,   249,  2758,   112,   249,  1042,   134,   244,   249,
     157,   244,   120,  1049,   175,   160,    76,   249,  1054,    38,
     112,   175,  1058,  1059,  1060,  1061,  1062,  1063,  1064,    76,
     249,  1067,   181,    76,   181,   345,   346,   249,   181,   249,
     249,   181,   244,   181,   193,   181,   244,   181,   249,   181,
     193,   249,   112,   249,   181,   193,   181,   193,  3263,   193,
      76,   193,   181,    82,  1100,   112,   193,   249,   193,   112,
    1106,    76,   164,   181,   193,   249,  1112,  1113,   223,  1115,
    1116,   100,   101,   249,   121,   193,  1122,   134,   249,   181,
     400,   134,    76,   125,   239,   249,   112,   244,   181,    76,
     249,   193,   285,   286,   244,   249,   249,   112,  1144,   249,
     249,   249,   244,   249,  1150,   249,  1152,  1153,   134,   244,
     249,   181,   249,   242,   249,   249,  2791,   181,   112,  1165,
     249,  1167,    76,   193,   181,   112,  1172,    76,   181,  1175,
     232,   249,    76,   249,   181,    76,   193,   249,  1184,   181,
     193,   181,  1188,   997,   998,  1191,    76,  1193,  1194,  1195,
     249,   244,  3367,   195,   249,   181,   249,   249,   112,   223,
    1206,  1207,  1208,   112,    76,   212,   181,   193,   112,  1215,
    1216,   112,   249,  1219,  1220,  1221,  1222,   249,   193,   249,
     244,   249,   112,   249,    76,    76,   249,   181,   249,    76,
     249,   231,   249,    76,   181,   249,   249,   244,   249,   193,
     112,   249,   244,   181,   244,    59,   193,   249,   249,   249,
     181,   181,   181,  1259,  1260,  1261,  1262,  1263,   148,    59,
     112,   112,    76,   249,  1602,   112,   249,   181,  1274,   112,
    1276,   249,   181,  1279,   249,   249,    76,   181,   249,   193,
     181,   181,   249,   249,   193,   249,    76,    76,  2484,   193,
     249,   181,   193,   181,   249,   249,   249,  1303,   112,  1305,
     249,  1307,   249,   193,  1310,   181,   244,   181,   181,   181,
    1316,   249,   112,   244,   244,   244,  1322,   191,   249,   249,
     249,   193,   112,   112,   249,  1331,  1332,  1333,  1334,   181,
     181,   249,   181,   249,   181,   249,   249,   249,   181,   249,
     249,   193,   193,   249,   244,   249,   193,   223,   249,   249,
     193,   249,   225,  1359,  1360,  1361,   244,  2992,  2993,  1365,
     249,   249,   249,  1369,    76,    76,   249,   181,   244,   249,
     244,   244,  1378,    76,  1380,  1381,  1382,  1383,   227,   193,
     249,   181,   249,  1389,   249,  1391,   249,  1393,  1394,   249,
     249,   181,   181,   193,   249,   244,   249,   249,   249,  1405,
     112,   112,   249,   193,   193,  1411,   249,  1413,  1414,   112,
    1416,  1417,  1418,  1419,  1420,  1421,  1422,  1423,  1424,  1425,
    1426,  1427,  1428,  1429,  1430,  1431,  1432,  1433,  1434,  1435,
    1436,  1437,  1438,  1439,  1440,  1441,  1442,  1443,  1444,  1445,
    1446,  1447,  1448,  1449,  1450,  1451,  1452,  1453,  1454,  1455,
     249,    76,   809,   810,   811,   812,    76,    76,  1464,   249,
     249,   818,   318,   319,   320,   249,  1804,   249,   249,   181,
     181,   249,   249,  1479,  1480,  1481,   181,    76,   181,   249,
     249,   193,   193,   249,    76,    76,     0,   112,    76,   219,
     193,   204,   112,   112,    76,    76,   223,   223,    76,    76,
    1318,  1507,  1320,   181,   223,  1323,   249,     7,  3143,     7,
     249,   249,   249,   112,   871,    76,   249,   249,    76,  1337,
     112,   112,  1528,   125,   112,   249,   249,   101,    76,   249,
     112,   112,  1538,  1539,   112,   112,    76,   249,   249,   244,
     249,   249,   249,   223,   249,   249,   249,  1553,   125,  1555,
     249,   112,  1558,   249,   112,   249,   181,  1563,   249,  1565,
      76,   181,   181,   249,   112,   249,   244,  1573,   193,   323,
    1576,   249,   112,   193,   193,  1581,   181,    76,  1584,   249,
     249,  1587,   181,  1589,  1590,  1591,  1592,   249,  1594,   181,
     181,   181,  1598,   181,   193,   249,   112,   249,   249,   181,
     181,   193,   193,   181,   181,   193,    76,  1613,  1614,   249,
     181,   193,   193,   112,  1620,   193,   193,  1623,   181,    76,
     181,  1627,  3257,   181,   249,  1631,   249,   181,  1634,   249,
     249,  1637,   193,   181,    76,   193,   249,  1643,   249,   244,
    1646,   181,   112,  1649,   249,   193,   249,   181,  1654,   249,
     249,   249,   249,   193,   244,   112,   249,   249,   249,   249,
    1666,   249,   249,    76,   249,   181,   249,   249,   249,   249,
     112,   249,  2113,   244,   249,   381,    59,   193,   249,   249,
      76,   244,   181,    76,   249,   249,   249,   249,   249,    59,
     244,   249,   249,    59,   193,   249,   378,  1703,    76,   112,
      59,   249,  1708,  1709,   249,  1711,  1712,    59,  1714,   249,
     244,   181,    76,    76,    59,   249,   112,   249,   249,   112,
      59,    76,   249,   193,   181,   372,   330,  3362,    76,  1735,
      76,   327,   249,   249,   112,    76,   193,  1743,   372,   181,
     125,  1747,   249,  1749,  1750,   125,   249,  1753,   112,   112,
     249,   193,  1758,   125,   125,   249,  1762,   112,   125,    76,
    1766,  1767,  1768,   125,   112,   249,   112,   249,   181,   249,
      76,   112,    76,    76,  1780,  1781,  1782,   249,    76,   249,
     193,  1787,  1788,  1789,  1790,   181,  1792,   181,   181,  1795,
    1796,  1797,   249,  1799,  1800,   112,   249,   193,   125,  1805,
     193,  1807,   181,   181,   181,   181,   112,   249,   112,   112,
      76,   125,    76,   125,   112,   193,   125,   181,   181,   372,
    1826,  1827,  1828,  1829,   249,    76,   181,  1833,   249,   193,
     193,   249,   249,   181,   372,   181,   249,  1843,   193,   181,
     181,   181,  1848,   372,  1850,   193,   112,   193,   112,   372,
     244,   125,   193,   249,    72,   249,   249,   181,   181,   249,
      76,   112,   125,    76,   181,   244,   117,   244,   244,   125,
     249,   249,   249,   249,    76,   181,   193,   181,   181,  1885,
     249,   181,   125,   181,   181,   249,   249,   193,   249,   193,
     193,   181,   249,   249,   249,   193,   112,    59,   181,   112,
     249,   249,   244,   249,   244,   249,  1912,   249,   249,   249,
     112,   249,   249,   249,    76,   181,   249,   181,  1924,   181,
     244,   244,   380,   249,   249,   249,   249,   193,  2266,   193,
     181,   249,   249,    76,    76,    76,  1942,  1943,    76,    59,
      76,   249,   193,   249,   244,   249,   249,   244,   249,   249,
     112,   249,   249,   249,   244,  1961,    76,  1963,   249,   249,
     249,   244,   249,   115,   116,   181,   249,  1973,   181,   112,
     112,   112,  1978,   125,   112,   372,   112,   193,  1984,   181,
     193,  1987,   244,   249,    76,   249,  1992,   249,   249,   249,
    1996,   193,   112,  1999,  2000,  2001,  2002,  2003,  2004,  2005,
    2006,  2007,  2008,  2009,  2010,  2011,  2012,  2013,  2014,  2015,
    2016,  2017,  2018,  2019,  2020,  2021,  2022,    76,   181,   181,
     112,  2027,  2028,  2029,  2030,  2031,  2032,  2033,  2034,  2035,
    2036,   193,   249,   249,   249,   181,   249,   249,   181,   181,
     181,  2047,   125,   181,   372,   181,   514,   249,   249,    76,
     193,   193,   193,   112,   125,   193,  2062,   193,   117,   372,
     249,   181,  2068,   125,   511,   372,   249,   372,   249,    76,
    2076,  2077,   125,   193,    76,  2081,   125,   181,    59,    76,
     181,   244,  2088,  2089,  2090,   112,   249,  2093,  2094,   181,
     181,  2097,   181,   372,   372,    76,    76,    76,   244,  2105,
     181,   193,    76,   249,   181,   112,   249,   249,   249,   125,
     112,   249,   125,   249,  2120,   112,  2122,  2123,   372,   372,
      76,  2127,   181,  2129,   181,   181,  2132,   372,    76,  2135,
     372,   112,   112,   112,   193,  2141,  2142,   372,   112,  2145,
     244,   372,  2148,   244,   372,   249,   372,  2153,   249,  2155,
     372,  2157,   181,   244,   181,   244,   112,   249,   249,    76,
     249,   181,  2500,   244,   112,   372,   193,   244,   249,   117,
     249,   125,   249,    76,   181,   372,  2182,  2183,  2184,   181,
    2186,  2187,    76,   249,   181,   244,   193,   244,   244,   341,
     342,   193,   249,   249,   372,   112,   193,   372,  2204,   372,
     181,   181,   181,   372,   372,  2211,   249,   181,  2214,   112,
     249,  2217,   193,   193,   193,   244,   125,  2223,   112,   193,
     249,   249,   249,   181,   244,   181,   181,  2233,  2234,   249,
    2236,  2237,  2238,   181,   249,  2241,  2242,   193,  2244,  2245,
    2246,  2247,   249,  2249,  2250,   193,  2252,   249,  2254,  2255,
    2256,  2257,   249,  2259,  2260,   249,   249,   249,   249,   249,
    2266,   223,  2268,   249,   181,  2271,   249,   249,  2274,   249,
     249,  2277,  2278,  2279,   372,   249,   193,  2283,   181,   372,
    2286,    76,   181,   249,  2290,   181,   244,   181,  2294,   244,
     193,   249,  2298,   249,   249,  2633,   244,   249,   249,   193,
    2306,   453,   454,   455,   456,   457,   181,   459,   181,   461,
     462,   463,   464,   465,   466,   467,   468,   112,   470,   471,
     472,   473,   474,   475,   476,   477,   478,   479,   480,   481,
     482,   483,   249,   485,   486,   487,   488,   489,   490,   491,
     492,   493,   494,   249,   249,   244,   249,    59,   244,   249,
     249,   249,   249,   249,   249,   249,    59,   249,    76,    76,
      76,  2699,   249,   249,    76,   249,   249,   249,   249,   244,
     249,   244,   249,    76,   249,   249,   249,    76,   249,   249,
      76,   249,   249,   249,   249,   249,   181,   249,   249,   249,
     249,   249,   249,   181,   112,   112,   112,  2403,   193,  2405,
     112,   125,  2408,   249,   372,    76,   249,   249,   249,   112,
    2416,  2417,  2418,   112,  2420,  2421,   112,  2423,  2424,  2425,
    2426,  2427,  2428,  2429,  2430,  2431,  2432,  2433,  2434,  2435,
     249,  2437,  2438,  2439,  2440,  2441,  2442,  2443,  2444,  2445,
     249,   112,   249,   249,  2450,  2451,  2452,  2453,  2454,  2455,
    2456,  2457,  2458,  2459,   249,    76,   244,  2463,  2464,   249,
     249,   249,   249,   181,   181,   181,    76,   249,  2474,   181,
     249,   249,    76,  2479,   125,   193,   193,   193,   181,   125,
    2486,   193,   181,  2489,   249,   181,  2492,   249,   249,   249,
     193,   112,   249,    76,   193,   181,   249,   193,   249,   249,
     249,   249,   112,  2509,   249,  2511,   181,   249,   112,   249,
     181,  2517,  2518,  2519,  2520,  2521,    76,  2523,   181,   181,
    2526,  2527,   193,  2529,  2530,   372,  2532,  2533,  2866,   112,
    2536,   249,   249,   249,  2540,  2541,   249,  2543,  2544,   249,
    2546,  2547,   249,   249,   249,   181,  2552,  2553,   249,   249,
     249,   249,   112,   249,   249,  2561,  2562,   249,   244,  2565,
     181,  2567,  2568,   249,  2570,  2571,   249,   249,    76,   244,
     249,   181,   193,    76,   249,    76,   249,   181,   249,  2585,
     249,   244,   244,   193,   249,   249,   249,   249,  2594,   193,
     181,   249,   181,   181,   181,   249,   181,   249,   181,  2605,
    2606,  2607,    76,   372,   112,  2611,  2612,  2613,   244,   112,
     193,   112,  2618,   249,  2620,  2621,  2622,   249,  2624,  2625,
    2626,   181,  2628,  2629,   249,   249,   372,  2633,   249,   249,
     249,  2637,   249,   193,  2640,  2641,   249,   249,   112,   249,
     249,   372,   249,  2649,  2650,   249,    76,    76,    76,   249,
      76,    76,    76,   244,  2660,   244,   244,   244,   249,   244,
     249,   249,   249,   184,   249,   249,   249,   249,   249,   249,
     249,   125,    94,   181,    76,    94,   125,   249,   181,   372,
     181,   181,   112,   112,   112,   193,   112,   112,   112,   249,
     193,   181,   193,   125,   125,  2701,   249,  2703,  2704,  2705,
    2706,  2707,  2708,  2709,  2710,  2711,  2712,   181,  2714,  2715,
     112,  2717,  2718,  2719,  2720,  2721,  2722,  2723,  2724,   193,
    2726,  2727,  2728,  2729,  2730,   249,   249,  2733,  2734,  2735,
    2736,  2737,  2738,  2739,  2740,  2741,  2742,  2743,  2744,  2745,
    2746,   249,  2748,   249,   244,  2751,   249,   249,   249,   249,
     249,   181,   181,   181,   244,   181,   181,   181,   249,   249,
     249,   249,   249,   193,   193,   193,   249,   193,   193,   193,
     181,   181,  2778,  2779,  2780,   249,  2782,    76,   249,   181,
    2786,   249,  2788,   249,  2790,   249,  2792,   323,  2794,   249,
    2796,   193,  2798,   249,  2800,   181,    76,   181,   249,   249,
     125,   181,   249,   249,   249,  2811,   249,  2813,  2814,   249,
    2816,  2817,   181,   112,   249,   181,   181,   125,   362,   249,
     249,   249,  2828,   249,   249,   249,   249,   125,  2834,  2835,
    2836,  2837,   112,   244,   244,  2841,  2842,  2843,   249,   249,
    2846,  3179,   249,   249,  2850,   249,   249,   249,   249,  2855,
     249,  2857,  2858,  2859,    76,  2861,  2862,   249,   244,   249,
     244,   249,   249,   249,   244,   249,  2872,  2873,   249,   249,
     249,   249,  2878,   249,   249,   244,   249,   249,   244,   244,
     249,    76,   181,   249,   249,   249,   249,    76,   249,   249,
     112,   249,   249,   249,   193,   249,   249,   181,   181,   249,
     249,   181,   249,  2909,   181,  2911,  2912,  2913,  2914,  2915,
    2916,  2917,  2918,   193,  2920,  2921,  2922,   112,  2924,  2925,
    2926,   181,  2928,   112,  2930,  2931,   249,  2933,  2934,  2935,
    2936,   249,  2938,  2939,   249,  2941,  2942,  2943,   249,  2945,
    2946,  2947,  2948,  2949,  2950,  2951,  2952,  2953,  2954,  2955,
     249,  2957,  2958,  2959,  2960,  2961,  2962,  2963,    76,   181,
     244,   244,  3300,  2969,   249,   249,   249,   244,   249,   249,
     249,   193,   249,  2979,  2980,   249,   249,   511,   249,   181,
     249,   249,   249,   181,   244,  2991,   181,   181,  2994,   249,
     181,   249,   181,   181,   112,   249,   181,   249,   193,   249,
     249,  3007,  3008,  3009,   193,  3011,   249,   249,   280,    76,
     282,   283,   284,   285,   286,   287,   288,   289,  3024,  3025,
    3026,   249,  3028,   249,   249,   181,  3032,   249,   249,   249,
     249,   249,  3038,  3039,  3040,  3041,  3042,   249,   249,   249,
    3046,  3047,   244,  3049,  3050,   112,   244,   249,   249,   249,
     244,   249,   249,   244,   249,   249,   244,   249,   249,   244,
     249,   249,   181,   181,   249,  3071,  3072,  3073,  3074,  3075,
    3076,  3077,   249,  3079,   249,   193,  3082,   249,   249,  3085,
     181,  3087,  3088,  3089,  3090,  3091,  3092,  3093,   244,  3095,
    3096,   249,  3098,   249,   249,  3101,  3102,  3103,   181,  3105,
    3106,  3107,  3108,  3109,  3110,  3111,   249,  3113,  3114,   181,
    3116,  3117,   249,  3119,   181,  3121,  3122,  3123,  3124,  3125,
    3126,  3127,   181,   249,   249,   244,   193,  3133,  3134,   249,
     249,   249,   181,   249,  3140,  3141,  3142,   249,   249,   181,
     249,   249,  3148,   244,    76,  3151,  3152,  3153,   249,   249,
      76,   249,   249,   249,  3160,  3161,   181,  3163,  3164,  3165,
     249,   244,  3168,   249,   249,   125,   249,  3173,  3174,  3175,
    3176,  3177,   244,   125,  3180,  3181,  3182,   249,   249,   249,
     112,   249,   249,   181,   181,   244,   112,    76,  3194,  3195,
     249,  3197,  3198,   249,  3200,   244,  3202,  3203,   249,  3205,
     249,  3207,   244,  3209,   249,  3211,   249,   249,  3214,  3215,
    3216,  3217,  3218,  3219,  3220,  3221,  3222,   249,  3224,   244,
    3226,  3227,   249,   112,   249,  3231,   181,  3233,  3234,  3235,
    3236,  3237,  3238,  3239,   372,  3241,   249,  3243,   181,  3245,
    3246,  3247,  3248,  3249,  3250,  3251,   244,   244,   249,   181,
    3256,   249,   249,  3259,  3260,   181,  3262,   249,   249,  3265,
    3266,   193,  3268,   249,  3270,  3271,   249,   193,  3274,   249,
    3276,   249,    76,   249,  3280,  3281,   249,   249,   249,  3285,
    3286,  3287,   372,  3289,   249,    76,  3292,  3293,  3294,   244,
     249,   181,   181,   249,   249,   249,  3302,   327,  3304,   249,
    3306,   244,  3308,  3309,   193,  3311,   249,  3313,   112,   249,
    3316,  3317,  3318,   249,  3320,   249,  3322,   249,  3324,  3325,
    3326,   112,  3328,   249,  3330,   181,   249,  3333,   181,  3335,
    3336,  3337,  3338,  3339,  3340,  3341,   249,  3343,  3344,   249,
    3346,  3347,  3348,   249,  3350,   327,  3352,  3353,  3354,  3355,
    3356,   249,   249,   249,   244,  3361,    76,   249,    76,   249,
     249,   249,  3368,   249,  3370,  3371,  3372,   249,   249,   249,
    3376,   249,   249,   249,    76,  3381,    76,   181,   249,   249,
    3386,  3387,  3388,  3389,   181,  3391,  3392,  3393,   244,   193,
     181,   244,   112,   249,   112,  3401,   249,   249,  3404,  3405,
    3406,  3407,   193,  3409,   249,  3411,   249,   249,  3414,  3415,
     112,  3417,   112,   249,   249,    76,   249,  3423,   249,  3425,
     249,  3427,  3428,  3429,  3430,  3431,  3432,  3433,  3434,   249,
    3436,   181,  3438,   249,  3440,   249,  3442,   249,  3444,  3445,
    3446,  3447,  3448,   249,   249,   249,  3452,   244,    76,   249,
    3456,   112,   249,   249,  3460,   249,   249,   249,   249,   181,
     249,   181,  3468,   181,   249,  3471,  3472,   181,  3474,  3475,
    3476,   125,  3478,   193,   249,   193,   249,   249,  3484,   181,
    3486,   181,  3488,   249,   112,   249,   249,  3493,  3494,   249,
     249,   193,    76,   193,   244,  3501,   181,   323,  3504,   249,
    3506,  3507,  3508,  3509,  3510,  3511,  3512,    76,  3514,   181,
    3516,   249,   372,  3519,  3520,  3521,  3522,  3523,  3524,   249,
     181,  3527,   244,   249,  3530,   249,   249,   249,   112,   249,
     244,   249,   193,  3539,   249,   249,  3542,  3543,   249,    76,
     249,   249,   249,   112,   249,  3551,   249,   249,   249,   249,
     249,  3557,  3558,   181,  3560,  3561,   249,  3563,   249,   244,
     249,   249,   249,  3569,   249,   193,  3572,  3573,  3574,  3575,
    3576,  3577,   244,   249,   249,   112,  3582,   249,  3584,  3585,
    3586,   249,  3588,  3589,   249,   249,   249,   249,   249,  3595,
    3596,   249,   249,   249,   249,   125,   249,   181,   249,  3605,
     249,  3607,   249,  3609,   249,  3611,   249,  3613,   249,   193,
     249,  3617,   181,  3619,  3620,  3621,  3622,  3623,   249,   249,
     249,   249,   249,   249,   193,  3631,  3632,   249,  3634,  3635,
      76,  3637,   249,   249,   249,   249,    76,  3643,   249,  3645,
     249,   249,   249,   249,   181,   249,  3652,   249,  3654,  3655,
     249,  3657,  3658,  3659,    76,  3661,   193,   249,   249,   249,
    3666,    76,  3668,  3669,  3670,   249,   112,   249,   249,  3675,
      76,  3677,   112,  3679,   249,  3681,   249,   249,   249,  3685,
     249,  3687,   249,  3689,   249,  3691,   249,  3693,   249,  3695,
     112,   249,   249,   249,  3700,  3701,   249,   112,  3704,   249,
    3706,  3707,   249,  3709,  3710,  3711,   112,   249,  3714,   249,
    3716,  3717,   249,  3719,    76,  3721,   249,   249,   249,  3725,
      76,    76,   249,  3729,    76,  3731,   249,  3733,   249,  3735,
     249,  3737,  3738,  3739,  3740,   181,   249,   249,   249,  3745,
      76,   181,   249,  3749,   249,  3751,   249,   193,    76,  3755,
     112,   249,   249,   193,  3760,   249,   112,   112,   249,   181,
     112,  3767,  3768,   249,  3770,    76,   181,  3773,   125,  3775,
     249,   193,  3778,  3779,   249,   181,   112,  3783,   193,  3785,
     249,  3787,   249,  3789,   112,     4,     5,   193,   249,  3795,
    3796,  3797,    11,    12,   249,   249,    15,   249,    17,   249,
      76,   112,   249,   249,   249,   249,   249,   249,   249,   249,
      76,    76,   249,    32,    33,   249,   249,   249,   249,   181,
     249,    40,    41,    42,    43,   181,   181,   249,    76,   181,
     249,   193,   249,    52,   249,    54,   112,   193,   193,   249,
     249,   193,   249,   249,   249,   181,   112,   112,   249,    68,
     249,    70,   249,   181,    73,   249,   249,   193,   249,    78,
     249,   249,    81,    76,   112,   193,    85,    86,    87,   249,
     181,    76,    91,   249,   249,    94,    95,   249,   249,    98,
      99,   249,   193,   102,   249,   104,   105,   249,   107,   108,
     249,   249,   172,   249,   249,   249,   249,   249,   117,   112,
     249,   120,   249,   249,   249,   181,   249,   112,   249,   249,
     249,   249,   249,   249,   249,   181,   181,   193,   249,   249,
     249,   249,   249,   249,   249,   144,   145,   193,   193,   148,
     149,   150,   249,   181,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   193,   249,   249,   249,   595,
     935,   625,   360,   304,   304,   597,   359,  1051,   177,  1047,
     304,  1014,   476,   182,   377,   946,   415,   186,   181,   188,
     189,   190,   407,   249,   963,   438,   181,   528,   494,   198,
     193,    -1,  1011,   249,   249,  1224,  1563,   319,   193,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   249,    -1,    -1,    -1,    -1,    -1,   226,    -1,   228,
     229,    -1,    -1,   232,   233,   234,   235,    -1,   237,    -1,
     239,   240,   241,   242,    -1,    -1,   245,   246,   247,    -1,
      -1,   250,    -1,   252,   253,   254,   249,    -1,    -1,   258,
      -1,    -1,    -1,    -1,   249,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   273,   274,   275,    -1,    -1,   278,
     279,   280,   281,   282,    -1,   284,   285,   286,   287,   288,
     289,   290,    -1,    -1,   293,   294,   295,    -1,    -1,   298,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   312,   313,   314,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   323,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     339,    -1,    -1,   342,   343,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   356,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   377,   378,
     379,    -1,    -1,   382,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   396,   397,    -1,
     399,    -1,    -1,    -1,   403,    -1,    -1,   406,   407,    -1,
      -1,   410,   411,    -1,    -1,    -1,   415,   416,    -1,    -1,
     419,    -1,    -1,   422,    -1,    -1,    -1,   426,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   438,
      -1,    -1,    -1,    -1,    -1,    -1,   445,    -1,   447,    -1,
      -1,    -1,    -1,   452,    -1,    -1,    -1,   456,    -1,    -1,
      -1,   460,    -1,    -1,   463,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   474,   475,   476,    -1,    -1,
      -1,    -1,    -1,    -1,   483,    -1,    -1,    -1,    -1,    -1,
     489,    -1,    -1,    -1,    -1,   494,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   503,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   512,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   520,   521,    -1,    -1,    -1,    -1,    -1,   527,   528,
     529,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   545,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   573,    -1,   575,    -1,    -1,    -1,
      -1,    -1,   581,   582,   583,   584,    -1,    -1,    -1,   588,
     589,   590,   591,   592,    -1,   594,   595,   596,   597,    -1,
      -1,    -1,    -1,    -1,   603,    -1,   605,   606,    -1,   608,
     609,   610,    -1,    -1,    -1,    -1,    -1,   616,   617,   618,
      -1,    -1,    -1,    -1,    -1,   624,    -1,    -1,    -1,    -1,
      -1,   630,    -1,    -1,    -1,   634,    -1,    -1,    -1,    -1,
     639,    -1,    -1,    -1,    -1,    -1,   645,    -1,   647,   648,
     649,    -1,    -1,   652,    -1,   654,    -1,   656,    -1,   658,
      -1,   660,    -1,    -1,    -1,   664,    -1,   666,    -1,   668,
      -1,   670,    -1,   672,    -1,    -1,    -1,    -1,    -1,    -1,
     679,    -1,    -1,    -1,   683,    -1,   685,    -1,   687,    -1,
      -1,    -1,    -1,   692,   693,   694,   695,    -1,   697,   698,
      -1,    -1,   701,    -1,    -1,   704,    -1,    -1,    -1,    -1,
      -1,    -1,   711,    -1,    -1,   714,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   722,    -1,    -1,    -1,    -1,   727,    -1,
      -1,    -1,    -1,    -1,    -1,   734,   735,    -1,    -1,    -1,
     739,    -1,    -1,    -1,    -1,    -1,   745,    -1,   747,   748,
     749,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   757,   758,
      -1,    -1,    -1,    -1,   763,   764,    -1,   766,   767,   768,
     769,    -1,   771,    -1,   773,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   782,   783,    -1,    -1,   786,   787,    -1,
      -1,    -1,   791,   792,    -1,   794,    -1,   796,    -1,    -1,
     799,   800,   801,   802,    -1,   804,   805,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   814,    -1,   816,    -1,    -1,
      -1,   820,   821,    -1,   823,    -1,    -1,    -1,    -1,    -1,
     829,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   838,
     839,    -1,    -1,    -1,    -1,   844,   845,    -1,    -1,    -1,
     849,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   861,   862,   863,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   874,    -1,   876,    -1,   878,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     889,    -1,   891,    -1,   893,   894,    -1,    -1,    -1,    -1,
     899,    -1,    -1,    -1,    -1,   904,    -1,    -1,    -1,    -1,
      -1,   910,    -1,    -1,    -1,    -1,    -1,   916,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   925,    -1,    -1,    -1,
      -1,   930,    -1,   932,   933,   934,   935,    -1,   937,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   946,   947,    -1,
     949,    -1,    -1,   952,   953,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   962,   963,   964,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   980,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   990,    -1,    -1,   993,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1006,    -1,  1008,
      -1,  1010,    -1,    -1,    -1,  1014,    -1,  1016,    -1,    -1,
      -1,    -1,  1021,    -1,    -1,  1024,  1025,    -1,  1027,    -1,
    1029,    -1,  1031,    -1,  1033,    -1,  1035,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1047,    -1,
      -1,    -1,  1051,    -1,  1053,    -1,    -1,  1056,  1057,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1068,
      -1,  1070,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1097,  1098,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1109,    -1,  1111,    -1,    -1,    -1,  1115,    -1,    -1,    -1,
      -1,  1120,  1121,    -1,  1123,  1124,  1125,  1126,    -1,    -1,
      -1,    -1,    -1,    -1,  1133,    -1,  1135,    -1,  1137,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1166,    -1,    -1,
    1169,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1186,  1187,    -1,
    1189,    -1,    -1,    -1,    -1,    -1,    -1,  1196,  1197,  1198,
    1199,  1200,  1201,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1213,    -1,    -1,    -1,  1217,  1218,
      -1,    -1,    -1,    -1,  1223,  1224,  1225,  1226,  1227,  1228,
      -1,  1230,    -1,  1232,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1264,    -1,    -1,    -1,  1268,
      -1,    -1,    -1,    -1,  1273,    -1,  1275,    -1,    -1,  1278,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1286,    -1,  1288,
      -1,  1290,    -1,    -1,    -1,  1294,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1327,    -1,
      -1,  1330,    -1,    -1,    -1,    -1,    -1,  1336,    -1,  1338,
    1339,  1340,  1341,  1342,  1343,    -1,  1345,  1346,  1347,  1348,
    1349,  1350,    -1,  1352,  1353,    -1,  1355,  1356,  1357,  1358,
      -1,    -1,    -1,  1362,    -1,  1364,    -1,  1366,  1367,  1368,
      -1,  1370,  1371,  1372,  1373,  1374,  1375,    -1,  1377,    -1,
    1379,    -1,    -1,    -1,    -1,  1384,  1385,  1386,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1395,    -1,    -1,    -1,
      -1,  1400,    -1,    -1,    -1,    -1,    -1,  1406,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1456,    -1,    -1,
      -1,    -1,    -1,  1462,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1479,    95,    -1,    -1,    -1,  1484,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   111,    -1,    -1,
      -1,  1500,    -1,    -1,  1503,    -1,    -1,    -1,  1507,    -1,
      -1,    -1,    -1,  1512,  1513,    -1,    -1,    -1,    -1,    -1,
      -1,  1520,  1521,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1531,    -1,    -1,    -1,  1535,  1536,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1546,    -1,    -1,
      -1,  1550,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1568,
      -1,  1570,    -1,    -1,    -1,  1574,    -1,  1576,    -1,    -1,
    1579,    -1,  1581,    -1,  1583,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1596,    -1,  1598,
      -1,  1600,    -1,    -1,    -1,   219,   220,    -1,    -1,    -1,
      -1,  1610,  1611,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1619,    -1,    -1,  1622,    -1,    -1,    -1,    -1,  1627,   243,
      -1,  1630,    -1,    -1,  1633,    -1,    -1,  1636,    -1,    -1,
    1639,    -1,    -1,  1642,    -1,    -1,  1645,    -1,    -1,  1648,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1662,    -1,  1664,  1665,    -1,    -1,    -1,
    1669,  1670,  1671,  1672,  1673,    -1,    -1,  1676,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1701,    -1,    -1,    -1,   320,   321,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1720,    -1,   337,    -1,  1724,    -1,  1726,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1755,    -1,    -1,    -1,
    1759,  1760,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1769,    -1,  1771,  1772,  1773,  1774,    -1,  1776,  1777,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1786,    -1,   403,
     404,   405,   406,    -1,  1793,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1801,    -1,  1803,  1804,    -1,  1806,    -1,    -1,
      -1,  1810,  1811,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1834,    -1,    -1,    -1,    -1,
      -1,    -1,  1841,    -1,    -1,  1844,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1861,  1862,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1874,  1875,    -1,  1877,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1888,
      -1,   505,    -1,   507,  1893,    -1,   510,    -1,   512,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    95,    -1,   530,  1916,   532,   533,
      -1,    -1,    -1,    -1,  1923,    -1,    -1,   541,    -1,   543,
     111,    -1,    -1,    -1,    -1,    -1,    -1,   551,   552,   553,
      -1,   555,   556,   557,    -1,  1944,    -1,   561,   562,   563,
     564,    -1,    -1,    -1,   568,   569,    -1,    -1,   572,   573,
     574,   575,    -1,   577,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   587,   588,   589,    -1,    -1,  1977,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1989,    -1,    -1,    -1,    -1,  1994,  1995,    -1,  1997,  1998,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2023,  2024,  2025,  2026,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   219,   220,
      -1,    -1,    -1,    -1,    -1,    -1,  2045,    -1,    -1,  2048,
      -1,    -1,    -1,    -1,  2053,  2054,    -1,    -1,  2057,  2058,
      -1,    -1,   243,    -1,    -1,    -1,  2065,  2066,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2090,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2107,    -1,
      -1,    -1,    -1,    -1,  2113,    -1,    -1,    -1,    -1,  2118,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2127,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   320,
     321,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   337,    -1,    -1,    -1,
    2159,    -1,    -1,    -1,    -1,    -1,    -1,  2166,    -1,    -1,
      -1,    -1,    -1,    -1,  2173,    -1,    -1,    -1,    -1,    -1,
    2179,  2180,  2181,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2193,    -1,    -1,    -1,    -1,  2198,
      -1,  2200,  2201,    -1,  2203,    -1,    -1,    -1,    -1,    -1,
      -1,  2210,    -1,    -1,  2213,    -1,    -1,  2216,    -1,    -1,
      -1,  2220,   403,   404,   405,   406,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2240,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2253,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  2276,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2287,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2300,  2301,    -1,    -1,    -1,    -1,    -1,    -1,  2308,
      -1,    -1,    -1,    -1,    -1,    -1,  2315,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   505,    -1,   507,    -1,    -1,   510,
      -1,   512,    -1,    -1,    -1,    -1,    -1,    -1,  2337,    -1,
    2339,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   530,
      -1,   532,   533,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     541,    -1,   543,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   555,   556,   557,    -1,    -1,    -1,
     561,   562,   563,   564,    -1,    -1,    -1,   568,   569,    -1,
      -1,   572,   573,   574,   575,    -1,   577,    -1,    -1,    -1,
    2399,  2400,    -1,    -1,    -1,    -1,   587,   588,   589,    -1,
      -1,  2410,  2411,    -1,    -1,  2414,  2415,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  2446,  2447,  2448,
    2449,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  2462,    -1,    -1,    -1,    -1,  2467,    -1,
      -1,    -1,    -1,  2472,  2473,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  2505,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2564,    -1,  2566,    -1,    -1,
      -1,    -1,    -1,    -1,  2573,    -1,    -1,    -1,  2577,    -1,
    2579,    -1,    -1,    -1,  2583,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  2591,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2608,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2639,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  2662,    -1,  2664,    -1,    -1,  2667,    -1,
    2669,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2677,    -1,
      -1,    -1,  2681,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2694,    -1,  2696,    16,    -1,
      -1,  2700,    -1,  2702,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    43,    -1,    -1,    -1,    -1,
      48,    -1,  2731,  2732,    -1,    -1,    -1,    -1,    -1,    57,
      -1,    -1,    60,    61,    -1,    -1,    -1,    65,    66,    -1,
      -1,  2750,    -1,  2752,    -1,    -1,    -1,  2756,    -1,    -1,
      -1,  2760,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   121,    -1,    -1,    -1,    -1,    -1,  2808,
    2809,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2819,    -1,  2821,    -1,   142,   143,   144,   145,    -1,    -1,
      -1,    -1,    -1,   151,    -1,    -1,    -1,   155,    -1,    -1,
     158,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   180,    -1,   182,    -1,    -1,  2866,    -1,    -1,
    2869,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   196,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2888,
      -1,  2890,   210,   211,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   219,    -1,   221,    -1,    -1,   224,    -1,   226,  2908,
     228,   229,    -1,    -1,    -1,   233,    -1,   235,   236,   237,
     238,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   253,    -1,   255,    -1,    -1,
      -1,    -1,    -1,   261,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2967,    -1,
      -1,    -1,   290,  2972,    -1,    -1,    -1,   295,   296,   297,
      -1,   299,    -1,   301,    -1,   303,    -1,    -1,    -1,    -1,
     308,   309,   310,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  3003,  3004,    -1,    -1,    -1,    -1,
     328,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   337,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   349,   350,   351,   352,   353,   354,   355,    -1,   357,
     358,   359,   360,    -1,    -1,    -1,    -1,    -1,    -1,   367,
     368,    -1,   370,    -1,    -1,    -1,   374,    -1,    -1,  3058,
      -1,  3060,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   387,
      -1,   389,    -1,    -1,    -1,    -1,   394,    -1,    -1,    -1,
      -1,   399,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   419,    -1,    -1,    -1,   423,   424,   425,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   443,    -1,    -1,   446,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  3147,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  3187,    -1,
    3189,    -1,    -1,  3192,    -1,    -1,    -1,    -1,   516,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  3261,    -1,  3263,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  3284,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  3297,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  3366,  3367,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  3413,    -1,    -1,    -1,    -1,    -1,
    3419,    -1,  3421,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  3455,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  3492,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  3536,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  3547,    -1,
    3549,    -1,    -1,    -1,    -1,    -1,  3555,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  3571,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  3593,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  3604,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  3627,    -1,
    3629,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  3642,    -1,    -1,    -1,    -1,     3,    -1,
       5,     6,  3651,    -1,    -1,    -1,    11,    -1,    13,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    21,    -1,    23,    24,
      25,    26,    27,    -1,    -1,    -1,    31,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    40,    41,    42,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    50,    -1,    52,    53,    54,
      -1,    56,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    64,
      -1,    -1,    -1,    -1,    -1,    70,    71,    -1,    -1,    74,
      -1,    -1,    -1,    78,    79,    80,    81,    -1,    -1,    -1,
      85,    -1,    87,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    99,    -1,    -1,   102,   103,   104,
     105,   106,   107,    -1,   109,   110,    -1,  3756,    -1,    -1,
      -1,    -1,    -1,   118,    -1,    -1,    -1,   122,   123,    -1,
      -1,    -1,   127,   128,    -1,    -1,   131,   132,    -1,    -1,
      -1,  3780,    -1,    -1,    -1,    -1,   141,    -1,    -1,    -1,
      -1,    -1,    -1,  3792,   149,  3794,    -1,   152,   153,    -1,
      -1,   156,    -1,    -1,   159,    -1,    -1,   162,   163,   164,
     165,   166,   167,   168,   169,   170,   171,   172,   173,    -1,
      -1,    -1,    -1,   178,   179,    -1,    -1,    -1,   183,    -1,
     185,    -1,   187,   188,   189,   190,    -1,   192,    -1,   194,
      -1,    -1,   197,    -1,    -1,    -1,    -1,   202,    -1,    -1,
     205,   206,   207,   208,   209,    -1,    -1,   212,    -1,    -1,
     215,   216,   217,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   232,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   241,    -1,    -1,    -1,
      -1,    -1,   247,    -1,    -1,    -1,   251,    -1,    -1,    -1,
      -1,   256,    -1,   258,   259,    -1,    -1,    -1,    -1,    -1,
     265,    -1,   267,    -1,    -1,   270,   271,    -1,    -1,   274,
      -1,    -1,   277,   278,   279,    -1,   281,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   292,   293,   294,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   306,   307,    -1,    -1,    -1,    -1,   312,   313,   314,
      -1,   316,    -1,    -1,    -1,    -1,    -1,   322,    -1,   324,
     325,    -1,    -1,    -1,   329,    -1,    -1,    -1,    -1,    -1,
      -1,   336,    -1,    -1,   339,   340,    -1,    -1,    -1,   344,
      -1,    -1,   347,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   361,    -1,    -1,   364,
      -1,    -1,    -1,    -1,   369,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   384,
      -1,   386,    -1,   388,    -1,    -1,   391,   392,   393,    -1,
     395,   396,    -1,   398,    -1,    -1,   401,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   411,    -1,    -1,    -1,
      -1,   416,   417,   418,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   426,   427,   428,   429,   430,   431,    -1,   433,    -1,
      -1,    -1,    -1,    -1,   439,   440,    -1,    -1,    -1,    -1,
      -1,    -1,   447,    -1,    -1,    -1,   451,   452,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     495,   496,    -1,    -1,    -1,    -1,    -1,   502,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     515,    -1,   517,   518,   519,   520,   521,   522,   523,   524,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     535,   536,    -1,    -1,   539,   540,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   554,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     565,   566,   567,    -1,    -1,   570,   571,    -1,     3,    -1,
       5,     6,    -1,    -1,    -1,    -1,    11,   582,    13,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    21,    -1,    23,    24,
      25,    26,    27,    -1,    -1,    -1,    31,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    40,    41,    42,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    50,    -1,    52,    53,    54,
      -1,    56,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    64,
      -1,    -1,    -1,    -1,    -1,    70,    71,    -1,    -1,    74,
      -1,    -1,    -1,    78,    79,    80,    81,    -1,    -1,    -1,
      85,    -1,    87,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    99,    -1,    -1,   102,   103,   104,
     105,    -1,   107,    -1,   109,   110,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   118,    -1,    -1,    -1,   122,   123,    -1,
      -1,    -1,   127,   128,    -1,    -1,   131,   132,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   141,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   149,    -1,    -1,   152,   153,    -1,
      -1,   156,    -1,    -1,   159,    -1,    -1,   162,   163,   164,
     165,   166,   167,   168,   169,   170,   171,   172,   173,    -1,
      -1,    -1,    -1,   178,   179,    -1,    -1,    -1,   183,    -1,
     185,    -1,   187,   188,   189,   190,    -1,   192,    -1,   194,
      -1,    -1,   197,    -1,    -1,    -1,    -1,   202,    -1,    -1,
     205,   206,   207,   208,   209,    -1,    -1,   212,    -1,    -1,
     215,   216,   217,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   232,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   241,    -1,    -1,    -1,
      -1,    -1,   247,    -1,    -1,    -1,   251,    -1,    -1,    -1,
      -1,   256,    -1,   258,   259,    -1,    -1,    -1,    -1,    -1,
     265,    -1,   267,    -1,    -1,   270,   271,    -1,    -1,   274,
      -1,    -1,   277,   278,   279,    -1,   281,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   292,   293,   294,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   306,   307,    -1,    -1,    -1,    -1,   312,   313,   314,
      -1,   316,    -1,    -1,    -1,    -1,    -1,   322,    -1,   324,
     325,    -1,    -1,    -1,   329,    -1,    -1,    -1,    -1,    -1,
      -1,   336,    -1,    -1,   339,   340,    -1,    -1,    -1,   344,
      -1,    -1,   347,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   361,    -1,    -1,   364,
      -1,    -1,    -1,    -1,   369,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   384,
      -1,   386,    -1,   388,    -1,    -1,   391,   392,   393,    -1,
     395,   396,    -1,   398,    -1,    -1,   401,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   411,    -1,    -1,    -1,
      -1,   416,   417,   418,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   426,   427,   428,   429,   430,   431,    -1,   433,    -1,
      -1,    -1,    -1,    -1,   439,   440,    -1,    -1,    -1,    -1,
      -1,    -1,   447,    -1,    -1,    -1,   451,   452,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     495,   496,    -1,    -1,    -1,    -1,    -1,   502,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     515,    -1,   517,   518,   519,   520,   521,   522,   523,   524,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     535,   536,    -1,    -1,   539,   540,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   554,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     565,   566,   567,    -1,    -1,   570,   571,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   582
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint16 yystos[] =
{
       0,     3,     5,     6,    11,    13,    21,    23,    24,    25,
      26,    27,    31,    40,    41,    42,    50,    52,    53,    54,
      56,    64,    70,    71,    74,    78,    79,    80,    81,    85,
      87,    99,   102,   103,   104,   105,   107,   109,   110,   118,
     122,   123,   127,   128,   131,   132,   141,   149,   152,   153,
     156,   159,   162,   163,   164,   165,   166,   167,   168,   169,
     170,   171,   172,   173,   178,   179,   183,   185,   187,   188,
     189,   190,   192,   194,   197,   202,   205,   206,   207,   208,
     209,   212,   215,   216,   217,   232,   241,   247,   251,   256,
     258,   259,   265,   267,   270,   271,   274,   277,   278,   279,
     281,   292,   293,   294,   306,   307,   312,   313,   314,   316,
     322,   324,   325,   329,   336,   339,   340,   344,   347,   361,
     364,   369,   384,   386,   388,   391,   392,   393,   395,   396,
     398,   401,   411,   416,   417,   418,   426,   427,   428,   429,
     430,   431,   433,   439,   440,   447,   451,   452,   495,   496,
     502,   515,   517,   518,   519,   520,   521,   522,   523,   524,
     535,   536,   539,   540,   554,   565,   566,   567,   570,   571,
     582,   592,   593,   594,   595,   596,   597,   598,   599,   605,
     606,   607,   608,   609,   610,   611,   612,   613,   614,   623,
     625,   626,   627,   628,   629,   630,   631,   632,   633,   634,
     636,   637,   638,   640,   641,   642,   643,   645,   646,   651,
     654,   656,   657,   658,   659,   660,   661,   662,   663,   664,
     665,   667,   668,   670,   674,   675,   676,   678,   679,   680,
     681,   684,   687,   688,   689,   690,   691,   694,   697,   700,
     702,   704,   706,   710,   711,   712,   713,   714,   715,   716,
     718,   719,   720,   721,   722,   723,   724,   725,   726,   727,
     731,   733,   735,   737,   739,   741,   743,   745,   747,   749,
     750,   754,   758,   760,   763,   764,   765,   766,   767,   768,
     769,   771,   772,   774,   775,   783,   784,   786,   788,   789,
     790,   791,   794,   795,   796,   797,   798,   799,   800,   801,
     805,   808,   810,   812,   814,   817,   820,   822,   824,   826,
     828,   829,   830,   831,   832,   835,   836,   838,   839,   840,
     841,   844,   847,   249,     7,   249,   249,   181,   244,   849,
     249,   849,    76,   112,   181,   193,   850,   850,   850,   125,
     249,   850,   249,   249,   849,   849,   249,   249,   849,   249,
     249,   849,   249,   249,   249,   249,   440,   249,   249,   249,
     249,   249,   249,   249,   249,   249,    38,    82,   100,   101,
     249,   813,   249,   849,   249,   849,   249,   249,   249,   249,
     249,   249,   849,   849,   249,   849,   249,   849,   249,   372,
     249,   249,   249,   249,   850,   249,   249,   249,   849,   249,
     249,   849,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   432,   249,   850,   249,   249,   849,   249,   249,
     849,   249,   249,   849,   249,   249,   249,   850,   125,   249,
     849,   249,   249,   249,   849,   125,   249,   372,   249,   249,
     849,   249,   849,   849,   249,   249,   850,   249,   849,   249,
     249,   125,   195,   249,   849,   125,   195,   249,   849,   249,
     249,   849,   849,   249,   249,   249,   849,   125,   849,   849,
     249,   249,   849,   249,   849,   125,   249,   249,   249,   249,
     249,   249,   249,   249,   849,   850,   249,   249,   849,   249,
     249,   850,   850,   249,   249,   249,   249,   249,   249,   249,
     249,   125,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   440,   249,   849,   249,   849,   249,   249,   849,
     849,   849,   249,   249,   850,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   850,   249,   249,   249,   249,
     249,   372,     0,   106,   594,   157,   849,    75,   135,   136,
     137,   138,   139,   302,   305,   326,   448,   600,   601,   602,
     603,   604,    44,    45,   434,   435,   436,   437,   438,    28,
      86,    88,   140,   413,   422,   503,   504,   505,   508,   509,
     849,    22,   133,   160,   239,   850,   850,   849,    49,   198,
     199,   200,   201,   615,   617,   618,   619,   620,   785,   849,
     624,   849,   223,   849,   850,   125,   378,   779,   849,   362,
     410,   583,   584,   585,   586,   635,   219,   639,   850,    22,
      47,    75,   160,   176,   223,   239,   240,   248,   263,   276,
     373,   383,   385,   432,   644,   648,   669,   673,   160,   223,
     239,   652,   653,    32,    51,   363,   437,   655,   204,   666,
     677,   849,   378,   779,   849,   849,   378,   849,   378,   849,
     378,   849,   849,   695,   849,   701,   849,   703,   849,   705,
     849,   707,   849,   695,   703,   707,     4,   717,   849,   223,
     849,   849,   849,   223,   223,   378,   779,   849,   850,   761,
     849,   777,   849,   849,   849,   849,   770,   849,   849,   413,
     773,   849,     8,    18,    30,    33,    34,    35,    36,    37,
      77,   146,   150,   245,   254,   257,   262,   263,   264,   266,
     269,   273,   319,   343,   348,   484,   497,   498,   776,   849,
     785,   785,   787,   849,   849,   849,   111,   321,   849,    19,
     404,   792,   793,   850,   181,   849,   850,   849,   378,   849,
     849,   213,   214,   268,   317,   390,   397,   415,   579,   580,
     581,    39,    63,    69,    82,    89,   119,   161,   184,   234,
     311,   345,   346,   400,   807,   809,   601,   604,   673,   850,
      12,    83,   101,   124,   147,   203,   211,   219,   250,   252,
     260,   318,   414,   420,   456,   827,   849,   460,   849,   849,
      46,    92,    93,   406,   503,   505,   506,   513,   525,   526,
     527,   528,   529,   531,   541,   542,   544,   546,   547,   551,
     557,   558,   559,   560,   837,   843,    95,   111,   219,   220,
     243,   320,   321,   337,   403,   404,   405,   406,   505,   507,
     510,   512,   530,   532,   533,   541,   543,   555,   556,   557,
     561,   562,   563,   564,   568,   569,   572,   573,   574,   575,
     577,   587,   588,   589,   842,   842,   552,   553,   842,   843,
     845,   846,   842,   537,   538,   848,   728,   779,   849,   249,
       7,     7,   249,   249,   249,   249,   850,   249,   849,   850,
     708,   709,   849,   709,   249,   249,   249,   850,   249,   850,
     125,   849,   647,   850,   639,   850,   850,   101,   850,   813,
     175,   249,   249,   249,   701,   695,   730,   781,   849,   811,
     850,   249,   849,   249,   249,   223,   249,   249,   850,   850,
     759,   782,   849,   759,   249,   692,   693,   849,   249,   703,
     821,   850,   823,   850,   707,   705,   755,   756,   849,   728,
     249,   249,   730,   728,   249,   850,   849,   249,   850,   849,
     249,   850,   223,   751,   752,   849,   850,   249,   152,   249,
     384,   249,   249,   249,   761,   249,   249,   249,   776,   850,
     849,   249,   249,   249,   849,   249,   249,   249,   849,   249,
     849,   249,   249,   849,   249,   249,   504,   249,   249,   249,
     249,   849,   249,   849,   717,   323,   381,   748,    59,   746,
     728,   249,   249,   850,   698,   699,   849,   249,   249,   777,
     850,    59,   734,   249,   728,    59,   732,    59,   736,   378,
     738,    59,   744,    59,   740,    59,   742,   849,   249,   249,
     249,   332,   377,   849,   249,   849,   249,   682,   683,   849,
     770,   685,   686,   849,   850,   249,   378,   849,   330,   850,
     850,   850,   850,   850,   850,   850,   249,   850,   327,   850,
     850,   372,   372,   125,   125,   125,   125,   125,   249,   249,
     125,   249,   249,   849,   249,   849,   249,   249,   249,   249,
     850,   849,   849,   849,   849,   850,   850,   164,   232,   850,
     849,   849,   849,   849,   849,   621,   849,   621,   622,   849,
     622,   408,   850,   850,   249,   849,   850,   125,   849,   779,
     849,   412,   849,   125,   849,   125,   849,   125,   125,   125,
     125,    14,    20,    67,   218,   246,   249,   256,   272,   334,
     335,   365,   366,   456,   850,   849,   249,   849,   849,   649,
     850,   249,   850,   850,   650,   850,   372,   249,   849,   649,
     249,   249,   249,   372,   372,   849,   372,   850,   372,    55,
     249,   849,   850,   249,   849,   850,   249,   849,   125,   249,
     849,   849,   249,   849,   850,   779,   412,    72,   849,   412,
     850,   849,   850,   849,   850,   850,   849,   849,   849,   849,
     849,   125,   103,   333,   249,   849,   850,   850,   850,   249,
     849,   249,   779,   412,   249,    10,    15,    68,   114,   215,
     379,   409,   850,   849,   849,   849,   849,   849,   849,   249,
     440,   813,   849,   249,   249,   849,   850,   850,   850,   850,
     850,   249,   849,   249,   249,   849,   249,   249,   249,   249,
     850,   249,   249,   849,   249,   850,   249,   850,   849,   850,
     148,   850,   849,   407,   849,   249,   249,   850,   849,   125,
     249,   849,   249,   849,   849,   117,   849,   850,   412,   850,
     850,   850,   380,   125,   850,   850,   804,   849,   803,   849,
     249,   850,   249,   850,   806,   849,   249,   849,   372,    84,
     249,   449,   590,   849,   249,   849,   249,   849,   849,   249,
     850,   849,   249,   298,   371,   420,   849,    16,    29,    43,
      48,    57,    60,    61,    65,    66,   121,   142,   143,   144,
     145,   151,   155,   158,   180,   182,   196,   210,   211,   219,
     221,   224,   226,   228,   229,   233,   235,   236,   237,   238,
     253,   255,   261,   290,   295,   296,   297,   299,   301,   303,
     308,   309,   310,   328,   337,   349,   350,   351,   352,   353,
     354,   355,   357,   358,   359,   360,   367,   368,   370,   374,
     387,   389,   394,   399,   419,   423,   424,   425,   443,   446,
     516,   850,   249,   850,   849,   849,   249,   850,   850,     9,
     825,   849,   849,   249,   249,   850,   849,   849,   249,   249,
     849,   849,   125,   115,   116,   125,   341,   342,   453,   454,
     455,   456,   457,   459,   461,   462,   463,   464,   465,   466,
     467,   468,   470,   471,   472,   473,   474,   475,   476,   477,
     478,   479,   480,   481,   482,   483,   485,   486,   487,   488,
     489,   490,   491,   492,   493,   494,   849,   849,   849,   849,
     372,   849,   849,   372,   514,   834,   834,   834,   834,   834,
     850,   849,   125,   126,   372,   849,   372,   834,   249,   849,
     849,   849,   249,   647,   125,   372,   850,   849,   372,   372,
     372,   372,   372,   372,   372,   372,   849,   849,   511,   850,
     125,   125,   849,   849,   850,   372,   372,   849,   372,   125,
     372,   372,   125,   125,   372,   372,   372,   372,   372,   849,
     849,   849,   249,   249,   249,   249,   249,   125,   850,   834,
     249,   125,   849,   249,   779,   849,   412,   249,   850,   850,
     850,   249,   249,   849,   708,   850,   709,   249,   850,   849,
     249,   249,   249,   850,   249,   850,   850,   249,   849,   781,
     850,   249,   249,   818,   819,   850,   249,   849,   223,   249,
     850,   850,   782,   849,   759,   693,   849,   249,   249,   756,
     757,   849,   757,   849,    62,   249,   850,   815,   816,   850,
     850,   850,   850,   249,   850,   849,   752,   753,   849,   753,
     849,   249,   850,   850,   372,   372,   850,   849,   249,   249,
     849,   849,   249,   815,   815,   249,   249,   249,   849,   381,
     850,   849,    59,   850,   809,   850,   699,   849,   249,   849,
      59,   850,   849,    59,   850,   849,    59,   850,   849,   378,
     850,   849,    59,   850,   849,    59,   850,   849,    59,   850,
     249,   249,   850,   249,   499,   500,   501,   249,   683,   249,
     850,   686,   849,   850,   849,   849,   850,   850,   249,   850,
     850,   850,   850,   850,   249,   850,   849,   249,   849,   125,
     249,   125,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   849,   849,
     249,   249,   372,   850,   249,   249,   249,   249,   850,   849,
     849,   850,   850,   125,   849,   850,   850,   125,   249,   249,
     849,   242,   249,   850,   849,   849,   849,   849,   849,   849,
     849,   372,   850,   249,   249,   850,   850,   850,   249,   249,
     249,   729,   780,   849,   249,   249,   249,   850,   249,   849,
     850,   249,   849,   850,   249,   850,   249,   249,   850,   849,
     849,   850,   849,   249,   850,   249,   850,   850,   850,   696,
     849,   696,   696,   696,   696,   849,   125,   125,   729,   249,
     850,   850,   850,   729,   249,   729,   849,   850,   850,   113,
     215,   249,   379,   440,   813,   849,   849,   850,   850,   850,
     850,   778,   849,   778,   849,   849,   696,   849,   849,   249,
     440,   849,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   850,   850,   850,   850,
     850,   177,   249,   407,   849,   249,   849,   249,   249,   249,
     849,   412,   850,   849,   117,   850,   249,   372,   849,   249,
     850,   249,   249,   249,   249,   249,   249,   249,   849,   249,
     849,   802,   849,   249,   249,   249,   849,   249,   249,   249,
     249,   249,   120,   249,   850,   850,   249,   850,   249,   850,
     249,   249,   249,   249,   249,   850,    17,   227,   444,   807,
     372,   807,    17,    58,   850,   807,   372,   372,   249,   338,
     849,   249,   249,   849,   850,   850,   850,   850,   249,   191,
     849,   807,   849,   849,   849,   225,   849,   227,   849,   338,
     849,   249,   849,   849,   849,   849,   849,   249,   849,   249,
     849,   121,   212,   849,   125,   249,   849,   849,   300,   849,
     849,   850,   850,   850,   849,   372,   338,   849,   850,   849,
     849,   849,   850,   849,   849,   849,   849,   849,   849,   249,
     849,   850,   849,   850,   850,   850,   850,   849,   249,   849,
     849,   249,   184,   850,   249,   372,   850,   850,   850,   249,
     849,   249,   249,   249,   849,   249,   249,   850,   249,   849,
     249,   249,   850,   125,   850,   850,   833,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   249,   849,    94,
      94,   849,   850,   125,   372,   849,   850,   850,   850,   249,
     849,   849,   125,   849,   849,   850,   125,   849,   849,   849,
     849,   249,   850,   249,   849,   412,   849,   249,   850,   850,
     249,   249,   249,   249,   249,   849,   850,   850,   850,   249,
     819,   850,   729,   249,   249,   849,   249,   849,   249,   850,
     849,   850,   757,   849,   850,   249,   816,   850,   249,   850,
     850,   850,   850,   249,   753,   849,   249,   440,   813,   249,
     249,   384,   249,   249,   249,   249,   671,   849,   849,   323,
     249,   849,   850,   249,   849,   850,   249,   849,   850,   249,
     849,   850,   249,   849,   850,   249,   849,   850,   108,   249,
     849,   833,   249,   849,   850,   249,   849,   850,   249,   849,
     850,   249,   332,   850,   249,   499,   249,   499,   249,   696,
     249,   849,   249,   849,   850,   249,   849,   849,   849,   849,
     849,   249,   249,   849,   249,   249,   249,   249,   249,   616,
     849,   249,   850,   850,   850,   249,   850,   850,   249,   850,
     249,   249,   249,   375,   849,   249,   242,   249,   157,   849,
     125,   157,   849,   125,   850,   249,   242,   249,   780,   850,
     850,   850,   850,   850,   849,   850,   375,   849,   849,   249,
     375,   850,   249,   850,   850,   850,   249,   849,   249,   249,
     249,   249,   362,   849,   849,   729,   850,   850,   850,   729,
     375,   849,   850,   249,   850,   850,   850,   850,   849,   113,
     215,   249,   379,   440,   850,   850,   850,   134,   249,   850,
     850,   249,   849,   249,   249,   813,   849,   249,   850,   249,
     249,   850,   249,   849,   849,   249,   850,   850,   850,   850,
     111,   249,   249,   321,   850,   249,   402,   849,   125,   249,
     849,   249,   372,   850,   849,   249,   372,   249,   850,   850,
     849,   849,   249,   249,   849,   249,   849,   249,   849,   249,
     249,   850,   249,   315,   249,   445,   849,   249,   249,   249,
     849,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   850,   249,   249,   230,   231,   249,   231,
     249,   849,   249,   249,   249,   249,   249,   849,   850,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   850,   850,   849,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   850,   249,   249,   850,   249,   249,   249,   249,
     249,   249,   174,   249,   249,   850,   249,   249,   249,   849,
     850,   249,   249,   850,   249,   850,   249,   849,   850,   249,
     849,   849,   249,   850,   849,   849,   850,   850,   850,   249,
     458,   850,   249,   458,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   249,   458,   850,   850,
     850,   850,   850,   850,   850,   850,   849,   849,   849,   849,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     249,   511,   849,   850,   849,   125,   849,   849,   125,   849,
     849,   850,   849,   375,   849,   249,   850,   249,   249,   850,
     850,   850,   729,   729,   249,   249,   850,   249,   850,   849,
     850,   249,   850,   850,   850,   249,   249,   249,   249,   850,
     849,   249,   372,   776,   280,   282,   283,   284,   285,   286,
     287,   288,   289,   672,   671,   849,   249,   850,   249,   850,
     850,   249,   249,   849,   850,   249,   850,   249,   249,   850,
     249,   249,   850,   249,   249,   372,   833,   108,   249,   249,
     850,   249,   249,   850,   249,   249,   850,   249,   249,   249,
     249,   850,   850,   850,   249,   249,   249,   249,   849,   249,
     249,   327,   327,   849,   849,   849,   616,   850,   850,   850,
     850,   850,   249,   849,   249,   249,   849,   849,   849,   849,
     850,   249,   249,   849,   249,   850,   249,   849,   249,   850,
     249,   849,   850,   249,   849,   249,   850,   249,   249,   249,
     249,   850,   249,   850,   249,   850,   850,   850,   849,   850,
     249,   850,   850,   850,   249,   850,   134,   249,   850,   249,
     850,   850,   850,   849,   850,   850,   850,   850,   850,   850,
     249,   249,   813,   850,   850,   850,   249,   850,   249,   849,
     850,   850,   850,   249,   850,   111,   249,   850,   249,   402,
     407,   849,   249,   850,   249,   249,   372,   850,   249,   249,
     850,   249,   849,   249,   849,   249,   249,   850,   249,   849,
     249,   249,   849,   249,   450,   249,   249,   231,   249,   849,
     849,   231,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   849,   249,   849,   249,   850,   850,   249,   850,
     849,   849,   849,   849,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   849,   849,   849,   849,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   849,   850,   850,   849,   249,
     375,   849,   849,   850,   249,   125,   850,   249,   729,   249,
     850,   249,   850,   249,   249,   850,   249,   249,   249,   249,
     813,   249,   849,   850,   850,   671,   671,   323,   850,   850,
     850,   850,   850,   850,   249,   850,   850,   850,   850,   850,
     850,   249,   249,   372,   833,   850,   850,   850,   850,   850,
     850,   249,   249,   850,   249,   850,   249,   249,   850,   850,
     249,   849,   850,   850,   850,   249,   850,   850,   849,   157,
     849,   157,   849,   849,   249,   850,   249,   249,   849,   249,
     850,   249,   249,   249,   850,   850,   850,   849,   242,   249,
     249,   850,   850,   850,   249,   762,   850,   134,   249,   850,
     850,   249,   850,   134,   249,   850,   249,   850,   850,   850,
     249,   850,   850,   249,   249,   813,   850,   249,   249,   850,
     249,   849,   850,   850,   249,   249,   249,   249,   407,   850,
     850,   249,   249,   249,   249,   372,   850,   249,   849,   249,
     849,   249,   849,   249,   849,   249,   450,   249,   849,   249,
     249,   849,   249,   849,   249,   849,   249,   249,   440,   442,
     813,   849,   850,   849,   850,   850,   850,   850,   850,   249,
     458,   850,   850,   249,   458,   850,   850,   249,   850,   249,
     850,   850,   249,   458,   850,   850,   850,   850,   850,   850,
     249,   458,   850,   850,   249,   850,   850,   850,   850,   850,
     849,   849,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   849,   850,   849,
     249,   249,   849,   249,   849,   249,   249,   249,   249,   850,
     850,   850,   249,   850,   249,   249,   850,   249,   850,   249,
     850,   833,   850,   249,   249,   249,   850,   249,   850,   249,
     850,   249,   249,   849,   849,   850,   249,   850,   850,   850,
     249,   850,   249,   849,   849,   249,   249,   850,   249,   850,
     850,   850,   850,   249,   850,   850,   850,   249,   850,   762,
     134,   249,   850,   249,   762,   134,   249,   850,   850,   850,
     850,   850,   850,   249,   249,   813,   849,   849,   249,   850,
     850,   850,   249,   111,   249,   249,   249,   249,   849,   249,
     849,   249,   249,   249,   249,   249,   249,   249,   849,   850,
     249,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   249,   850,   249,   850,   850,   850,
     850,   850,   850,   850,   249,   458,   850,   850,   850,   850,
     249,   850,   850,   850,   249,   850,   850,   850,   850,   850,
     850,   850,   249,   458,   850,   249,   458,   850,   249,   458,
     850,   850,   850,   850,   850,   850,   850,   849,   249,   850,
     125,   849,   249,   850,   850,   249,   249,   249,   249,   249,
     249,   833,   833,   850,   249,   249,   249,   849,   849,   249,
     850,   850,   850,   850,   249,   249,   249,   249,   242,   249,
     850,   850,   249,   379,   850,   850,   762,   762,   850,   249,
     762,   134,   249,   850,   850,   850,   850,   850,   249,   850,
     850,   850,   850,   111,   249,   249,   249,   849,   249,   849,
     249,   249,   440,   249,   850,   850,   249,   850,   850,   249,
     850,   249,   850,   850,   249,   850,   249,   458,   249,   850,
     249,   458,   850,   249,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   249,   850,   249,   458,   850,   249,   850,
     249,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   249,   458,   850,   850,   850,   850,
     850,   850,   850,   249,   249,   850,   850,   850,   249,   249,
     833,   849,   850,   249,   372,   850,   850,   249,   850,   249,
     850,   850,   850,   249,   379,   850,   850,   249,   762,   249,
     379,   850,   762,   762,   850,   850,   850,   850,   249,   850,
     249,   813,   850,   850,   850,   249,   249,   849,   249,   849,
     849,   850,   850,   850,   850,   249,   850,   249,   850,   850,
     850,   850,   850,   850,   249,   249,   850,   850,   850,   850,
     850,   249,   850,   249,   850,   850,   850,   850,   850,   249,
     850,   249,   249,   850,   249,   850,   850,   850,   850,   850,
     850,   850,   249,   458,   850,   249,   458,   850,   850,   249,
     850,   249,   850,   850,   850,   850,   850,   249,   125,   249,
     850,   833,   850,   249,   849,   850,   616,   249,   850,   850,
     850,   850,   850,   249,   850,   249,   379,   850,   762,   850,
     249,   379,   850,   249,   762,   849,   850,   850,   850,   850,
     249,   850,   850,   850,   249,   849,   249,   249,   813,   249,
     249,   850,   249,   249,   850,   850,   850,   850,   249,   850,
     249,   850,   249,   850,   850,   850,   249,   850,   249,   850,
     249,   850,   850,   850,   249,   850,   249,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   249,   850,   249,   850,   850,   850,   850,   850,   249,
     249,   249,   833,   249,   249,   616,   850,   249,   249,   850,
     850,   850,   249,   850,   249,   249,   850,   249,   379,   762,
     249,   850,   850,   850,   850,   850,   850,   249,   850,   249,
     249,   850,   850,   249,   850,   249,   850,   249,   850,   850,
     850,   249,   849,   850,   850,   850,   249,   849,   249,   849,
     249,   850,   850,   249,   850,   249,   850,   850,   850,   850,
     850,   850,   850,   249,   850,   249,   850,   249,   850,   850,
     850,   850,   850,   850,   850,   249,   249,   850,   249,   249,
     850,   249,   249,   249,   850,   249,   850,   850,   249,   850,
     850,   249,   850,   850,   249,   249,   249,   850,   249,   850,
     249,   850,   249,   249,   249,   849,   249,   458,   850,   249,
     458,   850,   249,   458,   249,   249,   850,   249,   249,   850,
     249,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     249,   249,   850,   249,   850,   850,   850,   850,   850,   850,
     249,   850,   249,   849,   249,   850,   850,   850,   249,   849,
     249,   849,   249,   850,   849,   850,   850,   850,   850,   850,
     249,   458,   249,   850,   249,   849,   249,   850,   249,   850,
     850,   850,   850,   850,   249,   249,   249,   850,   249,   850,
     249,   458,   850,   249,   458,   850,   249,   458,   850,   249,
     850,   249,   849,    75,   249,   850,   850,   249,   249,   249,
     249,   849,   850,   249,   458,   850,   249,   458,   850,   850,
     249,   850,   249,   850,   850,   249,   850,   249,   850,   850,
     850,   249,   849,   249,   849,   850,   850,   850,   850,   850,
     249,   458,   249,   849,   850,   850,   249,   849,   850,   850,
     850,   850,   850,   850,   249,   850,   249,   249,   850,   850,
     850,   850,   249,   249,   850,   249,   458,   850,   249,   458,
     850,   850,   249,   249,    75,   249,   850,   850,   249,   850,
     249,   850,   249,   850,   249,   850,   249,   249,   249,   850,
     249,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   249,   249,   850,   850,   850,   249,   850,
     249,   850,   249,   850,   850,   850,   850,   249,   249,   249,
     850,   249,   249,   850,   850,   249,   249,   850,   849,   249,
     850,   249,   249,    75,   249,   850,   850,   249,   850,   850,
     850,   249,   249,   850,   850,   849,   249,   850,   249,   850,
     249,   850,   850,   249,   849,   849,   850,   850,   850,   249
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   591,   592,   593,   593,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   595,   596,   597,   597,   597,   597,   598,
     598,   599,   599,   599,   599,   599,   599,   599,   599,   599,
     599,   599,   599,   599,   599,   600,   601,   602,   602,   603,
     603,   604,   604,   605,   605,   605,   605,   605,   605,   605,
     605,   605,   605,   606,   607,   607,   607,   607,   607,   607,
     607,   607,   607,   607,   607,   607,   607,   608,   608,   609,
     609,   609,   609,   609,   610,   610,   610,   611,   611,   611,
     612,   612,   612,   612,   612,   612,   612,   612,   613,   613,
     613,   613,   613,   613,   614,   614,   615,   615,   615,   615,
     616,   616,   617,   617,   618,   618,   619,   619,   620,   620,
     621,   621,   621,   622,   623,   623,   624,   624,   624,   624,
     625,   625,   625,   625,   626,   626,   627,   627,   627,   627,
     628,   629,   630,   631,   632,   633,   633,   633,   633,   633,
     634,   634,   634,   634,   634,   634,   634,   634,   634,   634,
     634,   635,   635,   635,   635,   635,   635,   635,   635,   635,
     635,   635,   635,   635,   635,   635,   635,   635,   635,   635,
     635,   635,   635,   635,   635,   635,   635,   635,   635,   635,
     635,   635,   636,   636,   636,   636,   636,   636,   636,   636,
     636,   636,   636,   636,   636,   636,   636,   636,   636,   636,
     636,   636,   636,   636,   636,   637,   637,   638,   638,   639,
     639,   640,   640,   641,   641,   642,   642,   643,   644,   644,
     644,   645,   645,   645,   645,   645,   645,   645,   645,   645,
     645,   645,   645,   645,   645,   645,   645,   645,   645,   645,
     646,   647,   648,   648,   648,   648,   649,   649,   649,   649,
     650,   651,   651,   651,   651,   651,   652,   652,   653,   654,
     654,   654,   654,   654,   654,   654,   654,   654,   655,   655,
     656,   657,   658,   659,   659,   660,   661,   662,   663,   663,
     663,   663,   663,   663,   663,   663,   664,   664,   665,   665,
     665,   665,   666,   666,   667,   668,   669,   670,   670,   670,
     671,   671,   672,   672,   672,   672,   672,   672,   672,   672,
     672,   672,   673,   673,   673,   674,   675,   675,   676,   676,
     676,   677,   678,   679,   679,   679,   679,   679,   680,   680,
     681,   682,   682,   683,   683,   684,   685,   685,   686,   687,
     687,   687,   687,   687,   688,   688,   688,   688,   689,   689,
     689,   689,   690,   690,   690,   691,   692,   692,   693,   693,
     693,   694,   694,   695,   696,   696,   697,   698,   698,   699,
     699,   699,   700,   700,   701,   702,   702,   703,   704,   704,
     705,   706,   706,   707,   708,   709,   709,   710,   711,   711,
     712,   712,   713,   713,   714,   714,   715,   715,   716,   717,
     717,   717,   717,   718,   718,   718,   718,   718,   719,   719,
     719,   719,   719,   720,   720,   720,   720,   720,   720,   721,
     721,   722,   722,   723,   723,   723,   723,   724,   724,   725,
     726,   726,   726,   726,   726,   726,   726,   726,   727,   727,
     727,   727,   728,   728,   728,   728,   728,   728,   729,   729,
     730,   730,   731,   731,   732,   732,   732,   733,   733,   734,
     734,   734,   735,   735,   736,   736,   736,   737,   737,   738,
     738,   738,   738,   738,   739,   739,   740,   740,   740,   741,
     741,   742,   742,   742,   743,   743,   744,   744,   744,   745,
     745,   746,   746,   746,   747,   747,   748,   748,   748,   749,
     749,   749,   750,   750,   751,   751,   751,   752,   752,   752,
     752,   752,   753,   754,   754,   755,   755,   755,   756,   756,
     756,   757,   757,   758,   758,   759,   759,   760,   760,   761,
     761,   761,   761,   761,   761,   761,   761,   761,   761,   761,
     761,   761,   761,   761,   761,   761,   761,   761,   761,   761,
     761,   761,   761,   761,   761,   761,   761,   761,   761,   761,
     761,   761,   761,   761,   761,   761,   761,   761,   761,   761,
     761,   761,   761,   761,   761,   761,   761,   761,   761,   761,
     761,   761,   762,   763,   763,   764,   764,   764,   764,   764,
     765,   765,   765,   765,   765,   765,   765,   765,   765,   766,
     766,   767,   767,   767,   767,   767,   767,   767,   767,   767,
     767,   767,   768,   768,   768,   769,   769,   770,   771,   771,
     771,   772,   772,   772,   772,   772,   772,   773,   773,   773,
     773,   773,   773,   773,   773,   773,   774,   774,   774,   774,
     774,   774,   774,   774,   774,   774,   774,   774,   774,   774,
     774,   774,   774,   774,   774,   774,   774,   774,   774,   774,
     774,   774,   774,   775,   775,   775,   776,   776,   776,   776,
     776,   777,   778,   778,   779,   779,   779,   779,   780,   781,
     782,   782,   783,   783,   784,   784,   785,   785,   786,   786,
     787,   787,   787,   787,   788,   788,   789,   789,   789,   789,
     789,   789,   789,   789,   789,   789,   789,   789,   789,   789,
     789,   789,   789,   789,   790,   790,   790,   790,   790,   791,
     791,   791,   792,   793,   794,   794,   795,   795,   795,   796,
     796,   796,   797,   797,   797,   797,   797,   797,   797,   797,
     797,   797,   797,   797,   797,   798,   798,   799,   800,   800,
     800,   800,   800,   801,   801,   801,   801,   801,   801,   801,
     801,   801,   801,   801,   802,   802,   802,   802,   802,   802,
     802,   802,   802,   802,   802,   802,   803,   803,   804,   804,
     805,   805,   805,   805,   805,   805,   805,   805,   805,   805,
     805,   806,   806,   807,   807,   807,   807,   807,   807,   807,
     807,   807,   807,   807,   807,   807,   807,   807,   807,   807,
     807,   807,   807,   807,   807,   807,   807,   807,   807,   807,
     808,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   810,   811,   812,
     812,   813,   813,   813,   813,   813,   813,   813,   813,   813,
     813,   813,   814,   814,   814,   814,   814,   814,   814,   814,
     814,   814,   815,   815,   816,   817,   817,   818,   818,   819,
     820,   821,   822,   823,   824,   825,   825,   826,   826,   826,
     826,   826,   826,   826,   826,   826,   826,   826,   826,   826,
     826,   826,   826,   826,   826,   826,   826,   826,   826,   826,
     826,   826,   826,   826,   827,   827,   827,   828,   829,   829,
     829,   830,   830,   830,   830,   830,   830,   830,   830,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   832,   832,   832,   833,   833,
     834,   834,   835,   835,   835,   836,   836,   837,   837,   837,
     837,   837,   837,   837,   837,   837,   837,   837,   837,   837,
     837,   837,   837,   837,   837,   837,   837,   837,   837,   837,
     837,   837,   837,   837,   837,   837,   837,   837,   837,   837,
     837,   838,   838,   839,   839,   840,   840,   840,   840,   840,
     841,   841,   842,   842,   842,   842,   842,   842,   842,   842,
     842,   842,   842,   842,   842,   842,   842,   842,   842,   842,
     842,   842,   842,   842,   842,   842,   842,   842,   842,   842,
     842,   842,   842,   842,   842,   842,   842,   842,   842,   842,
     842,   842,   842,   842,   842,   842,   842,   842,   842,   842,
     842,   842,   843,   843,   844,   845,   845,   845,   846,   846,
     847,   847,   848,   848,   848,   849,   849,   850,   850,   850,
     850
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     2,     1,     2,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     5,     4,     2,     5,     6,     6,     2,
       6,     2,     3,     4,     6,     6,     7,    12,    12,     7,
       2,     2,     2,     2,     2,     4,     3,     2,     4,     3,
       4,     4,     6,     2,     4,     5,     4,     5,     4,     4,
       4,     4,     4,     3,     2,     4,     4,     3,     3,     3,
       3,     3,     3,     4,     3,     3,     3,     2,     4,     2,
       4,     4,     4,     4,     2,     3,     4,     2,     3,     4,
       2,     3,     5,     5,     7,     4,     5,     5,     2,     2,
       2,     2,     2,     2,     2,     2,     4,    10,     5,    11,
       4,     5,     3,     2,     3,     2,     3,     2,     3,     2,
      11,    13,    14,     5,     2,     2,     7,     9,    11,    12,
       2,     5,     6,     5,     2,     5,     2,     5,     4,     4,
       3,     3,     3,     3,     3,     2,     2,     6,     8,     3,
       2,     2,     3,     3,     3,     3,     4,     4,     3,     3,
       3,     3,     5,     4,     5,     6,     7,     3,     5,     4,
       5,     6,     7,     3,     2,     2,     3,     2,     2,     3,
       2,     2,     2,     3,     2,     3,     2,     2,     2,     2,
       2,     2,     1,     2,     3,     3,     3,     2,     5,     3,
       3,     3,     2,     3,     2,     3,     4,     4,     5,     3,
       3,     2,     4,     2,     3,     3,     4,     4,     3,     2,
       2,     2,     3,     2,     3,     5,     2,     2,     2,     2,
       3,     2,     2,     2,     2,     3,     3,     4,     4,     8,
       4,     4,     3,     4,     6,     7,     8,     3,     5,     4,
       4,     7,     2,     3,     3,     3,     2,     4,     1,     3,
       1,     2,     2,     2,     3,     4,     5,     6,     6,     3,
       4,     6,     7,     5,     3,     4,     4,     3,     1,     2,
       6,     2,     2,     2,     3,     2,     2,     4,     5,     4,
       2,     6,     8,     7,     5,     9,     3,     2,     3,     2,
       4,     3,     2,     2,     2,     2,     5,     5,     6,     7,
       3,     1,     1,     2,     1,     1,     1,     2,     1,     1,
       2,     1,     4,     5,     3,     3,     6,     5,     6,     3,
       2,     5,     6,     2,     2,     7,     9,     3,     2,     6,
       3,     1,     2,     3,     2,     3,     1,     2,     4,     2,
       4,     6,     8,     5,     2,     3,     4,     5,     2,     3,
       6,     7,     2,     3,     6,     3,     1,     2,     4,     5,
       6,     3,     2,     4,     1,     2,     3,     1,     2,     4,
       5,     6,     3,     2,     4,     3,     2,     4,     3,     2,
       4,     3,     2,     4,     3,     1,     2,     3,     3,     4,
       3,     2,     3,     2,     5,     2,     3,     4,     4,     5,
       5,     6,     6,     3,     4,     3,     2,     6,     2,     3,
       3,     4,     5,     3,     2,     9,     6,     2,     2,     3,
       9,     3,     9,     2,     3,     4,     5,     3,     4,     3,
       2,     3,     2,     7,     9,     8,    10,     3,     5,     6,
       6,     7,     1,     2,     6,     7,     8,     9,     1,     2,
       1,     2,     2,     3,     6,     4,     7,     2,     3,     6,
       4,     7,     2,     3,     6,     4,     7,     2,     3,     8,
      10,     4,     9,    11,     2,     3,     6,     4,     7,     2,
       3,     6,     4,     7,     2,     3,     6,     4,     7,     2,
       3,     6,     4,     7,     2,     3,     9,     7,    10,     2,
       3,     5,     5,     3,     2,     2,     3,     2,     3,     5,
       4,     6,     4,     2,     3,     2,     2,     3,     5,     3,
       2,     5,     4,     3,     4,     1,     2,     3,     2,    16,
      19,    20,    23,     9,    14,    16,    30,    12,    13,    14,
       5,     6,    16,    12,     3,     4,    13,     5,     6,     6,
       7,    11,     5,     6,     8,     9,    10,     9,    10,    11,
      10,    11,    12,    11,    12,    13,     5,     7,     6,     8,
       6,     9,     7,    10,     7,    11,     8,    12,     4,     6,
      12,     4,     5,     3,     2,     3,     4,     5,     6,     5,
       4,     5,     5,     6,     7,     7,     8,     7,     8,     4,
       3,     2,     5,     6,     7,     8,    10,     6,     7,     8,
       9,    11,     2,     5,     7,     3,     2,     4,     2,     5,
       7,     2,     4,     3,     4,     6,     5,     1,     3,     4,
       5,     6,     8,     9,    11,    12,     2,     4,     3,     3,
       3,     3,     4,     3,     4,     4,     3,     3,     3,     3,
       3,     3,     3,     3,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     3,     6,     2,     5,     4,     3,     7,
       6,     4,     1,     2,     4,     3,     5,     4,     3,     3,
       5,     4,     2,     2,     2,     2,    11,     4,     2,     2,
      11,    14,    12,    15,     2,     7,     2,     3,     5,     6,
       4,     6,     7,     7,     6,     8,     9,     7,     9,    10,
       5,     5,     7,     8,     2,     3,     4,     3,     3,     2,
       2,     2,     5,     3,     2,     3,     2,     4,     3,     2,
       4,     5,     2,     3,     4,     5,     5,     7,     5,     6,
       6,     6,     7,     7,     8,     2,     3,     3,     2,     4,
       6,     6,     8,     2,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     3,     4,     4,     5,     5,     6,
       6,     7,     7,     8,     8,     9,     1,     2,     1,     2,
       2,     2,     4,     3,     5,     4,     4,     5,     6,     4,
       4,     1,     2,     2,     3,     2,     3,     3,     3,     2,
       3,     4,     5,     6,     7,     2,     3,     4,     5,     4,
       3,     3,     3,     2,     4,     5,     6,     7,     2,     3,
       4,     1,     3,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     5,     5,     4,
       4,     4,     4,     4,     4,     4,     4,     3,     4,     4,
       4,     4,     4,     4,     5,     4,     5,     5,     5,     5,
       3,     4,     4,     4,     4,     4,     4,     4,     3,     3,
       4,     3,     3,     4,     5,     4,     4,     5,     5,     6,
       6,     7,     4,     4,     5,     4,     4,     4,     4,     4,
       3,     4,     3,     3,     3,     3,     4,     3,     4,     4,
       3,     4,     4,     4,     4,     4,     4,     3,     3,     4,
       5,     4,     4,     4,     4,     6,     6,     5,     5,     7,
       7,     4,     4,     5,     4,     4,     4,     3,     2,     3,
       4,     1,     2,     3,     4,     5,     1,     2,     3,     2,
       3,     4,     2,     3,     4,     5,     4,     4,     4,     2,
       2,     2,     1,     2,     4,     6,     5,     1,     2,     4,
       3,     2,     3,     2,     2,     1,     1,     2,     3,     3,
       3,     3,     4,     4,     4,     5,     7,     5,     7,     4,
       5,     3,     4,     4,     5,     6,     7,     8,     3,     4,
       6,     8,     4,     2,     3,     4,     5,     2,     2,    10,
      12,     2,     4,     7,     9,     9,     8,    11,    12,     2,
       9,    10,    12,    13,    14,    15,     9,    10,    12,    13,
      14,    15,    12,    13,    14,    15,     9,    10,    12,    13,
      14,    15,     9,     7,     5,    13,    11,     9,     9,     7,
       5,    13,    11,     9,     9,     7,     5,    13,    11,     9,
       8,    10,     8,    10,     9,    11,     9,    11,    11,    13,
      11,    13,     8,    10,    12,    14,     8,    10,    12,    14,
       8,    12,     9,    13,    10,    11,    13,    14,    15,    16,
      10,    11,    13,    14,    15,    16,     7,    13,    15,    17,
      19,    14,    16,    14,    16,    15,    17,    15,    17,    17,
      19,    17,    19,    14,    16,    18,    20,    13,    15,    17,
      19,    14,    16,    18,    20,     7,    11,    13,    17,    14,
      18,     8,    12,    14,    18,    15,    19,     8,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,     7,     9,    10,    10,    11,    12,
      13,    10,    11,    12,    13,     7,     8,     9,    10,    11,
      12,    13,    22,     5,     5,     2,     4,     5,     0,     2,
       0,     2,     4,     6,     8,     2,     3,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     3,     2,
       3,     2,     2,     2,     2,     5,     6,     3,     6,     3,
       2,     1,     2,     3,     2,     3,     4,     2,     2,     3,
       1,     2,     3,     2,     3,     2,     3,     2,     2,     2,
       2,     3,     2,     3,     3,     2,     3,     4,     2,     2,
       2,     2,     2,     3,     4,     2,     2,     2,     2,     2,
       2,     2,     3,     4,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     3,     4,     2,     3,     4,     2,
       2,     2,     2,     2,     2,     2,     2,     3,     3,     2,
       6,     3,     2,     3,     5,     2,     5,     3,     2,     3,
       2,     3,     2,     3,     2,     1,     1,     1,     1,     1,
       1
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yystacksize);

        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 197 "p.y" /* yacc.c:1646  */
    { 
          if(domain->solInfo().piecewise || domain->solInfo().freeplay) domain->solInfo().activatePiecewise();
          return 0;
         }
#line 5490 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 6:
#line 209 "p.y" /* yacc.c:1646  */
    { if(geoSource->setDirichlet((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; delete (yyvsp[0].bclist); }
#line 5496 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 7:
#line 211 "p.y" /* yacc.c:1646  */
    { if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5502 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 16:
#line 221 "p.y" /* yacc.c:1646  */
    { int j = geoSource->getLocalIndex();
          if(geoSource->elementLumpingWeightLocalSize(j)>0) geoSource->setLocalIndex(j+1); }
#line 5509 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 26:
#line 233 "p.y" /* yacc.c:1646  */
    {}
#line 5515 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 28:
#line 236 "p.y" /* yacc.c:1646  */
    {}
#line 5521 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 34:
#line 243 "p.y" /* yacc.c:1646  */
    {}
#line 5527 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 40:
#line 250 "p.y" /* yacc.c:1646  */
    { if(geoSource->setUsddLocation((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1;
          if(geoSource->setDirichlet((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0)    return -1; delete (yyvsp[0].bclist); }
#line 5534 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 55:
#line 267 "p.y" /* yacc.c:1646  */
    { domain->setMFTT((yyvsp[0].mftval).first, (yyvsp[0].mftval).second); }
#line 5540 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 56:
#line 269 "p.y" /* yacc.c:1646  */
    { domain->setHFTT((yyvsp[0].hftval).first, (yyvsp[0].hftval).second); }
#line 5546 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 89:
#line 303 "p.y" /* yacc.c:1646  */
    { if(geoSource->setDirichlet((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5552 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 90:
#line 305 "p.y" /* yacc.c:1646  */
    { if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5558 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 91:
#line 307 "p.y" /* yacc.c:1646  */
    { if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5564 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 92:
#line 309 "p.y" /* yacc.c:1646  */
    { if(geoSource->setDirichletFluid((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5570 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 93:
#line 311 "p.y" /* yacc.c:1646  */
    { if(geoSource->setDirichletFluid((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5576 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 94:
#line 313 "p.y" /* yacc.c:1646  */
    { if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5582 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 111:
#line 331 "p.y" /* yacc.c:1646  */
    { if(domain->setComplexNeuman((yyvsp[0].cxbclist)->n,(yyvsp[0].cxbclist)->d) < 0) return -1; }
#line 5588 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 113:
#line 334 "p.y" /* yacc.c:1646  */
    { if(domain->setComplexDirichlet((yyvsp[0].cxbclist)->n,(yyvsp[0].cxbclist)->d) < 0) return -1; }
#line 5594 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 117:
#line 339 "p.y" /* yacc.c:1646  */
    { if(geoSource->setDirichlet((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5600 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 118:
#line 341 "p.y" /* yacc.c:1646  */
    { if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5606 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 125:
#line 349 "p.y" /* yacc.c:1646  */
    {}
#line 5612 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 134:
#line 360 "p.y" /* yacc.c:1646  */
    {}
#line 5618 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 135:
#line 362 "p.y" /* yacc.c:1646  */
    {}
#line 5624 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 136:
#line 364 "p.y" /* yacc.c:1646  */
    {}
#line 5630 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 137:
#line 366 "p.y" /* yacc.c:1646  */
    {}
#line 5636 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 138:
#line 368 "p.y" /* yacc.c:1646  */
    {}
#line 5642 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 139:
#line 370 "p.y" /* yacc.c:1646  */
    {}
#line 5648 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 140:
#line 372 "p.y" /* yacc.c:1646  */
    {}
#line 5654 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 141:
#line 374 "p.y" /* yacc.c:1646  */
    {}
#line 5660 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 153:
#line 389 "p.y" /* yacc.c:1646  */
    { domain->solInfo().noninpc = true;
            sfem->setOrder((yyvsp[-2].ival)); 
            domain->solInfo().nsample = (yyvsp[-1].ival);
          }
#line 5669 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 154:
#line 396 "p.y" /* yacc.c:1646  */
    { domain->solInfo().inpc = true;
            sfem->setOrder((yyvsp[-1].ival));
          }
#line 5677 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 156:
#line 403 "p.y" /* yacc.c:1646  */
    { if ((yyvsp[-3].ival) == OutputInfo::Attribute)  geoSource->setAttributeGroup((yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1);
          else if ((yyvsp[-3].ival) == OutputInfo::Nodal)  geoSource->setNodeGroup((yyvsp[-2].ival)-1, (yyvsp[-1].ival));
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Group Type: %d\n", (yyvsp[-3].ival));  exit(-1); }
        }
#line 5686 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 157:
#line 408 "p.y" /* yacc.c:1646  */
    { int i;
          if ((yyvsp[-4].ival) == OutputInfo::Attribute)  {
            for(i=(yyvsp[-3].ival); i<(yyvsp[-2].ival)+1; ++i)
              geoSource->setAttributeGroup(i-1,(yyvsp[-1].ival)-1);
          }
          else if ((yyvsp[-4].ival) == OutputInfo::Nodal)  {
            for(i=(yyvsp[-3].ival); i<(yyvsp[-2].ival)+1; ++i)
              geoSource->setNodeGroup(i-1, (yyvsp[-1].ival));
          }
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Group Type: %d\n", (yyvsp[-4].ival));  exit(-1); }
        }
#line 5702 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 158:
#line 420 "p.y" /* yacc.c:1646  */
    { if ((yyvsp[-4].ival) == OutputInfo::Nodal) geoSource->setSurfaceGroup((yyvsp[-2].ival)-1, (yyvsp[-1].ival));
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Surface Group Type: %d\n", (yyvsp[-4].ival));  exit(-1); }
        }
#line 5710 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 160:
#line 427 "p.y" /* yacc.c:1646  */
    { geoSource->setGroupRandomProperty((yyvsp[-4].ival)-1,(yyvsp[-3].rprop),(yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 5716 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 161:
#line 431 "p.y" /* yacc.c:1646  */
    { domain->solInfo().curSweepParam = 0; }
#line 5722 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 162:
#line 433 "p.y" /* yacc.c:1646  */
    { domain->solInfo().curSweepParam = (yyvsp[-1].ival); }
#line 5728 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 163:
#line 435 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-1].fval)); }
#line 5734 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 164:
#line 437 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-3].fval));
            domain->setFrequencySet(domain->solInfo().curSweepParam);
            domain->addFrequencies1(2.0*PI*(yyvsp[-3].fval), 2.0*PI*(yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 5742 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 165:
#line 441 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-3].fval));
            domain->setFrequencySet(domain->solInfo().curSweepParam);
            domain->addFrequencies2(2.0*PI*(yyvsp[-3].fval), 2.0*PI*(yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 5750 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 166:
#line 445 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-4].fval));
            domain->setFrequencySet(domain->solInfo().curSweepParam);
            domain->addFrequencies(2.0*PI*(yyvsp[-4].fval), 2.0*PI*(yyvsp[-3].fval), (yyvsp[-2].ival), (yyvsp[-1].ival)); }
#line 5758 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 167:
#line 449 "p.y" /* yacc.c:1646  */
    {
          if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-9].fval));
          domain->setFrequencySet(domain->solInfo().curSweepParam);
          domain->addFrequencies(2.0*PI*(yyvsp[-9].fval), 2.0*PI*(yyvsp[-8].fval), 2, (yyvsp[-7].ival));
          domain->solInfo().getSweepParams()->isAdaptSweep = true;
          domain->solInfo().getSweepParams()->adaptSweep.maxP = (yyvsp[-4].ival);
          domain->solInfo().getSweepParams()->adaptSweep.numS = (yyvsp[-7].ival);
          if ((yyvsp[-6].ival) == SweepParams::KrylovGalProjection) 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 0; 
          else if ((yyvsp[-6].ival) == SweepParams::WCAWEGalProjection) 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 2;
          else
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 1;
          domain->solInfo().getSweepParams()->adaptSweep.w1 = 2.0*PI*(yyvsp[-9].fval);
          domain->solInfo().getSweepParams()->adaptSweep.w2 = 2.0*PI*(yyvsp[-8].fval);
          domain->solInfo().getSweepParams()->adaptSweep.atol = (yyvsp[-5].fval);
          domain->solInfo().getSweepParams()->adaptSweep.minRHS = (yyvsp[-3].ival);
          domain->solInfo().getSweepParams()->adaptSweep.maxRHS = (yyvsp[-2].ival);
          domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = (yyvsp[-1].ival);
          domain->solInfo().getSweepParams()->nFreqSweepRHS = (yyvsp[-2].ival);
        }
#line 5784 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 168:
#line 471 "p.y" /* yacc.c:1646  */
    {
          if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-9].fval));
          domain->setFrequencySet(domain->solInfo().curSweepParam);
          domain->addFrequencies(2.0*PI*(yyvsp[-9].fval), 2.0*PI*(yyvsp[-8].fval), 2, (yyvsp[-7].ival));
          domain->solInfo().getSweepParams()->isAdaptSweep = true;
          domain->solInfo().getSweepParams()->adaptSweep.maxP = (yyvsp[-4].ival);
          domain->solInfo().getSweepParams()->adaptSweep.numS = (yyvsp[-7].ival);
          if ((yyvsp[-6].ival) == SweepParams::KrylovGalProjection) 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 0; 
          else if ((yyvsp[-6].ival) == SweepParams::WCAWEGalProjection) 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 2;
          else
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 1;
          domain->solInfo().getSweepParams()->adaptSweep.w1 = 2.0*PI*(yyvsp[-9].fval);
          domain->solInfo().getSweepParams()->adaptSweep.w2 = 2.0*PI*(yyvsp[-8].fval);
          domain->solInfo().getSweepParams()->adaptSweep.atol = (yyvsp[-5].fval);
          domain->solInfo().getSweepParams()->adaptSweep.minRHS = (yyvsp[-3].ival);
          domain->solInfo().getSweepParams()->adaptSweep.maxRHS = (yyvsp[-3].ival);
          domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = -1;
          domain->solInfo().getSweepParams()->nFreqSweepRHS = (yyvsp[-3].ival);
          domain->solInfo().getSweepParams()->adaptSweep.ctolf = (yyvsp[-2].fval);
          domain->solInfo().getSweepParams()->adaptSweep.tol1f = (yyvsp[-1].fval);
        }
#line 5812 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 169:
#line 495 "p.y" /* yacc.c:1646  */
    {
          if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-4].fval));
          domain->setFrequencySet(domain->solInfo().curSweepParam);
          domain->addFrequencies(2.0*PI*(yyvsp[-4].fval), 2.0*PI*(yyvsp[-3].fval), 2, (yyvsp[-2].ival));
          domain->solInfo().getSweepParams()->isAdaptSweep = true;
          domain->solInfo().getSweepParams()->adaptSweep.maxP = 6;
          domain->solInfo().getSweepParams()->adaptSweep.numS = (yyvsp[-2].ival);
          if ((yyvsp[-1].ival) == SweepParams::KrylovGalProjection) {
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 0; 
             domain->solInfo().getSweepParams()->adaptSweep.atol = 1e-2;
             domain->solInfo().getSweepParams()->adaptSweep.minRHS = 8;
             domain->solInfo().getSweepParams()->adaptSweep.maxRHS = 48;
             domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = 4;
          }
          else if ((yyvsp[-1].ival) == SweepParams::WCAWEGalProjection) {
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
          domain->solInfo().getSweepParams()->adaptSweep.w1 = 2.0*PI*(yyvsp[-4].fval);
          domain->solInfo().getSweepParams()->adaptSweep.w2 = 2.0*PI*(yyvsp[-3].fval);
          domain->solInfo().getSweepParams()->nFreqSweepRHS = domain->solInfo().getSweepParams()->adaptSweep.maxRHS;
        }
#line 5849 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 175:
#line 535 "p.y" /* yacc.c:1646  */
    {
          if((yyvsp[-3].ival) == 1) {
            domain->solInfo().setDamping((yyvsp[-2].fval),(yyvsp[-1].fval));
            domain->solInfo().getSweepParams()->alphaD = (yyvsp[-1].fval);
            domain->solInfo().getSweepParams()->betaD = (yyvsp[-2].fval);
            domain->solInfo().setDamping((yyvsp[-2].fval),(yyvsp[-1].fval));
          }
          else return -1; // only RAYDAMP is allowed here
        }
#line 5863 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 176:
#line 547 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getSweepParams()->pade_pivot = true; domain->solInfo().getSweepParams()->pade_tol = (yyvsp[-1].fval); }
#line 5869 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 177:
#line 551 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getSweepParams()->pade_poles = true; }
#line 5875 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 178:
#line 553 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getSweepParams()->pade_poles = true; 
          domain->solInfo().getSweepParams()->pade_poles_sigmaL = (yyvsp[-2].fval); domain->solInfo().getSweepParams()->pade_poles_sigmaU = (yyvsp[-1].fval); }
#line 5882 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 179:
#line 558 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-1].fval));
          domain->setFrequencySet(domain->solInfo().curSweepParam);
          domain->addCoarseFrequency(2.0*PI*(yyvsp[-1].fval)); }
#line 5890 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 180:
#line 563 "p.y" /* yacc.c:1646  */
    { domain->addFrequencies(2.0*PI*(yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 5896 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 181:
#line 567 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getSweepParams()->freqSweepMethod = (yyvsp[-2].ival); 
          int &l = domain->solInfo().getSweepParams()->padeL,
              &m = domain->solInfo().getSweepParams()->padeM,
              &n = domain->solInfo().getSweepParams()->padeN;
          switch((yyvsp[-2].ival)) {
            case SweepParams::Taylor:
              n = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = (yyvsp[-1].ival)+1; // taylor
              break;
            case SweepParams::Pade1:
              n = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+m+1;
              break;
            case SweepParams::Pade:
            case SweepParams::Fourier:
              n = (yyvsp[-1].ival);
              domain->solInfo().getSweepParams()->nFreqSweepRHS = (int) ceil(float(l+m+1)/float(n));
              break;
            case SweepParams::PadeLanczos:
              n = (yyvsp[-1].ival);
              if(m%n != 0) m = (m/n+1)*n; // round m up to the nearest multiple of n
              l = m-1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = m/n;
              break;
            case SweepParams::GalProjection:
              n = (yyvsp[-1].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::KrylovGalProjection:
              n = (yyvsp[-1].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::WCAWEGalProjection:
              n = (yyvsp[-1].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
          }
        }
#line 5942 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 182:
#line 609 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getSweepParams()->freqSweepMethod = (yyvsp[-4].ival);
          int &l = domain->solInfo().getSweepParams()->padeL,
              &m = domain->solInfo().getSweepParams()->padeM,
              &n = domain->solInfo().getSweepParams()->padeN;
          switch((yyvsp[-4].ival)) {
            case SweepParams::Taylor:
              n = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = (yyvsp[-3].ival)+1; // taylor
              break;
            case SweepParams::Pade1:
              n = 1;
              l = (yyvsp[-2].ival);
              m = (yyvsp[-1].ival);
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+m+1;
              break;
            case SweepParams::Pade:
            case SweepParams::Fourier:
              n = (yyvsp[-3].ival);
              l = (yyvsp[-2].ival); 
              m = (yyvsp[-1].ival);
              domain->solInfo().getSweepParams()->nFreqSweepRHS = (int) ceil(float(l+m+1)/float(n));
              break;
            case SweepParams::PadeLanczos:
              n = (yyvsp[-3].ival);
              m = (yyvsp[-1].ival);
              if(m%n != 0) m = (m/n+1)*n; // round m up to the nearest multiple of n
              l = m-1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = m/n;
              break;
            case SweepParams::GalProjection:
              n = (yyvsp[-3].ival);
              l = (yyvsp[-2].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::KrylovGalProjection:
              n = (yyvsp[-3].ival);
              l = (yyvsp[-2].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::WCAWEGalProjection:
              n = (yyvsp[-3].ival);
              l = (yyvsp[-2].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
          }
        }
#line 5996 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 184:
#line 662 "p.y" /* yacc.c:1646  */
    { geoSource->binaryInput = bool((yyvsp[-1].ival)); }
#line 6002 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 185:
#line 664 "p.y" /* yacc.c:1646  */
    { geoSource->binaryInput = bool((yyvsp[-2].ival));
            std::string prefix = (yyvsp[-1].strval);
            clusterData_ = prefix + ".msh";
            decomposition_ = prefix + ".dec";
            connectivity_ = prefix + ".con";
            subdomains_ = prefix + ".sub";
          }
#line 6014 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 186:
#line 672 "p.y" /* yacc.c:1646  */
    { geoSource->binaryOutput = bool((yyvsp[-1].ival)); }
#line 6020 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 187:
#line 674 "p.y" /* yacc.c:1646  */
    { geoSource->binaryOutput = bool((yyvsp[-2].ival));
            int len = strlen((yyvsp[-1].strval));
            char *file = new char[len+5];
            strcpy(file, (yyvsp[-1].strval));
            strcat(file,".con");
            geoSource->setGlob(file);
          }
#line 6032 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 188:
#line 682 "p.y" /* yacc.c:1646  */
    { geoSource->setGeo((yyvsp[-1].strval)); }
#line 6038 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 189:
#line 684 "p.y" /* yacc.c:1646  */
    { geoSource->setDecomp((yyvsp[-1].strval)); }
#line 6044 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 190:
#line 686 "p.y" /* yacc.c:1646  */
    { geoSource->setGlob((yyvsp[-1].strval)); }
#line 6050 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 191:
#line 688 "p.y" /* yacc.c:1646  */
    { geoSource->setMatch((yyvsp[-1].strval)); }
#line 6056 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 192:
#line 690 "p.y" /* yacc.c:1646  */
    { geoSource->setCpuMap((yyvsp[-1].strval)); }
#line 6062 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 193:
#line 694 "p.y" /* yacc.c:1646  */
    { 
#ifdef STRUCTOPT	  
	  dynamic_cast<Domain_opt*>(domain)->addAnalysis((yyvsp[-1].ival)); 
#endif
	}
#line 6072 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 194:
#line 702 "p.y" /* yacc.c:1646  */
    {if(decInit==0) decInit = new DecInit(); }
#line 6078 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 195:
#line 704 "p.y" /* yacc.c:1646  */
    {decInit->file = strdup((yyvsp[-1].strval));}
#line 6084 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 196:
#line 706 "p.y" /* yacc.c:1646  */
    {decInit->nsubs = (yyvsp[-1].ival); }
#line 6090 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 197:
#line 708 "p.y" /* yacc.c:1646  */
    {decInit->weight = true; }
#line 6096 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 198:
#line 710 "p.y" /* yacc.c:1646  */
    {decInit->memory = true; }
#line 6102 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 199:
#line 712 "p.y" /* yacc.c:1646  */
    {decInit->exitAfterDec = true;}
#line 6108 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 200:
#line 714 "p.y" /* yacc.c:1646  */
    {decInit->skip = true;}
#line 6114 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 201:
#line 716 "p.y" /* yacc.c:1646  */
    {decInit->nosa = true; }
#line 6120 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 202:
#line 718 "p.y" /* yacc.c:1646  */
    {decInit->trivial = true; }
#line 6126 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 203:
#line 720 "p.y" /* yacc.c:1646  */
    {decInit->trivial = true; randomShuffle = bool((yyvsp[-1].ival)); }
#line 6132 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 204:
#line 722 "p.y" /* yacc.c:1646  */
    {decInit->fsgl = true; }
#line 6138 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 205:
#line 724 "p.y" /* yacc.c:1646  */
    {allowMechanisms = true; }
#line 6144 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 206:
#line 726 "p.y" /* yacc.c:1646  */
    {useScotch = true; }
#line 6150 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 207:
#line 730 "p.y" /* yacc.c:1646  */
    {}
#line 6156 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 208:
#line 732 "p.y" /* yacc.c:1646  */
    { weightList[(yyvsp[-2].ival)] = (yyvsp[-1].fval); }
#line 6162 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 209:
#line 736 "p.y" /* yacc.c:1646  */
    {}
#line 6168 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 210:
#line 738 "p.y" /* yacc.c:1646  */
    { fieldWeightList[(int)Element::Acoustic] = (yyvsp[-1].ival); }
#line 6174 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 211:
#line 740 "p.y" /* yacc.c:1646  */
    { fieldWeightList[(int)Element::Structural] = (yyvsp[-1].ival); }
#line 6180 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 212:
#line 742 "p.y" /* yacc.c:1646  */
    { fieldWeightList[(int)Element::Thermal] = (yyvsp[-1].ival); }
#line 6186 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 213:
#line 744 "p.y" /* yacc.c:1646  */
    { fieldWeightList[(int)Element::Fluid] = (yyvsp[-1].ival); }
#line 6192 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 214:
#line 747 "p.y" /* yacc.c:1646  */
    { (yyval.mftval).first = new MFTTData; (yyval.mftval).second = 0; }
#line 6198 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 215:
#line 749 "p.y" /* yacc.c:1646  */
    { (yyval.mftval).first = new MFTTData; (yyval.mftval).second = (yyvsp[-1].ival); }
#line 6204 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 216:
#line 751 "p.y" /* yacc.c:1646  */
    { (yyval.mftval).first->add((yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 6210 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 217:
#line 755 "p.y" /* yacc.c:1646  */
    { (yyval.hftval).first = new MFTTData; (yyval.hftval).second = 0; }
#line 6216 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 218:
#line 757 "p.y" /* yacc.c:1646  */
    { (yyval.hftval).first = new MFTTData; (yyval.hftval).second = (yyvsp[-1].ival); }
#line 6222 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 219:
#line 759 "p.y" /* yacc.c:1646  */
    { (yyval.hftval).first->add((yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 6228 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 220:
#line 763 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 6234 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 221:
#line 765 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-1].ival); }
#line 6240 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 222:
#line 767 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-3].ival); domain->setLoadFactorGrav((yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 6246 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 223:
#line 769 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-3].ival); domain->setLoadFactorTemp((yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 6252 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 224:
#line 771 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-5].ival); domain->setLoadFactorGrav((yyvsp[-5].ival), (yyvsp[-3].ival)); domain->setLoadFactorTemp((yyvsp[-5].ival), (yyvsp[-1].ival)); }
#line 6258 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 225:
#line 773 "p.y" /* yacc.c:1646  */
    { domain->setLoadFactor((yyval.ival), (yyvsp[-2].ival), (yyvsp[-1].fval)); }
#line 6264 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 226:
#line 775 "p.y" /* yacc.c:1646  */
    { domain->setLoadFactorMFTT((yyval.ival), (yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 6270 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 227:
#line 777 "p.y" /* yacc.c:1646  */
    { domain->setLoadFactorHFTT((yyval.ival), (yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 6276 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 235:
#line 790 "p.y" /* yacc.c:1646  */
    { geoSource->addCFrame((yyvsp[0].frame).num,(yyvsp[0].frame).d); }
#line 6282 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 236:
#line 794 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].coefdata).coefFlag = false; geoSource->addCoefInfo((yyvsp[-2].ival)-1,(yyvsp[0].coefdata)); }
#line 6288 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 237:
#line 796 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].coefdata).c[6][0] = (yyvsp[-7].fval);
          (yyvsp[0].coefdata).c[6][1] = (yyvsp[-6].fval);
          (yyvsp[0].coefdata).c[6][2] = (yyvsp[-5].fval);
          (yyvsp[0].coefdata).c[6][3] = (yyvsp[-4].fval);
          (yyvsp[0].coefdata).c[6][4] = (yyvsp[-3].fval);
          (yyvsp[0].coefdata).c[6][5] = (yyvsp[-2].fval);
          (yyvsp[0].coefdata).coefFlag = false;
          geoSource->addCoefInfo((yyvsp[-8].ival)-1,(yyvsp[0].coefdata)); }
#line 6301 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 238:
#line 805 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].coefdata).coefFlag = (yyvsp[-2].ival); geoSource->addCoefInfo((yyvsp[-3].ival)-1,(yyvsp[0].coefdata)); }
#line 6307 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 239:
#line 807 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].coefdata).c[6][0] = (yyvsp[-8].fval);
          (yyvsp[0].coefdata).c[6][1] = (yyvsp[-7].fval);
          (yyvsp[0].coefdata).c[6][2] = (yyvsp[-6].fval);
          (yyvsp[0].coefdata).c[6][3] = (yyvsp[-5].fval);
          (yyvsp[0].coefdata).c[6][4] = (yyvsp[-4].fval);
          (yyvsp[0].coefdata).c[6][5] = (yyvsp[-3].fval);
          (yyvsp[0].coefdata).coefFlag = (yyvsp[-2].ival);
          geoSource->addCoefInfo((yyvsp[-9].ival)-1,(yyvsp[0].coefdata)); }
#line 6320 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 240:
#line 818 "p.y" /* yacc.c:1646  */
    { (yyval.coefdata).zero(); (yyval.coefdata).setCoef((yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1,(yyvsp[-1].fval)); }
#line 6326 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 241:
#line 820 "p.y" /* yacc.c:1646  */
    { (yyval.coefdata).setCoef((yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1,(yyvsp[-1].fval)); }
#line 6332 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 242:
#line 824 "p.y" /* yacc.c:1646  */
    { (yyval.linfo) = new LayInfo(0); geoSource->addLay((yyvsp[-1].ival)-1,(yyval.linfo)); }
#line 6338 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 243:
#line 826 "p.y" /* yacc.c:1646  */
    { (yyvsp[-1].linfo)->add((yyvsp[0].ldata).lnum,(yyvsp[0].ldata).d,(yyvsp[0].ldata).matid); }
#line 6344 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 244:
#line 830 "p.y" /* yacc.c:1646  */
    { (yyval.linfo) = new LayInfo(1); geoSource->addLay((yyvsp[-1].ival)-1,(yyval.linfo)); }
#line 6350 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 245:
#line 832 "p.y" /* yacc.c:1646  */
    { (yyvsp[-1].linfo)->add((yyvsp[0].ldata).lnum,(yyvsp[0].ldata).d,(yyvsp[0].ldata).matid); }
#line 6356 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 246:
#line 836 "p.y" /* yacc.c:1646  */
    { (yyval.linfo) = new LayInfo(0); geoSource->addLay((yyvsp[-1].ival)-1,(yyval.linfo)); }
#line 6362 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 247:
#line 838 "p.y" /* yacc.c:1646  */
    { (yyvsp[-1].linfo)->add((yyvsp[0].ldata).lnum,(yyvsp[0].ldata).d,(yyvsp[0].ldata).matid); }
#line 6368 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 248:
#line 842 "p.y" /* yacc.c:1646  */
    { (yyval.linfo) = new LayInfo(1); geoSource->addLay((yyvsp[-1].ival)-1,(yyval.linfo)); }
#line 6374 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 249:
#line 844 "p.y" /* yacc.c:1646  */
    { (yyvsp[-1].linfo)->add((yyvsp[0].ldata).lnum,(yyvsp[0].ldata).d,(yyvsp[0].ldata).matid); }
#line 6380 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 250:
#line 848 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).lnum = (yyvsp[-10].ival)-1;
          (yyval.ldata).matid = -1; // this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[-9].fval); (yyval.ldata).d[1] = (yyvsp[-8].fval); (yyval.ldata).d[2] = (yyvsp[-7].fval);
	  (yyval.ldata).d[3] = (yyvsp[-6].fval); (yyval.ldata).d[4] = (yyvsp[-5].fval); (yyval.ldata).d[5] = (yyvsp[-4].fval);
	  (yyval.ldata).d[6] = (yyvsp[-3].fval); (yyval.ldata).d[7] = (yyvsp[-2].fval); (yyval.ldata).d[8] = (yyvsp[-1].fval);
          (yyval.ldata).d[9] = 0;  (yyval.ldata).d[10] = 0; (yyval.ldata).d[11] = 0; }
#line 6391 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 251:
#line 855 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).lnum = (yyvsp[-12].ival)-1;
          (yyval.ldata).matid = -1; // this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[-11].fval); (yyval.ldata).d[1] = (yyvsp[-10].fval); (yyval.ldata).d[2] = (yyvsp[-9].fval);
          (yyval.ldata).d[3] = (yyvsp[-8].fval); (yyval.ldata).d[4] = (yyvsp[-7].fval); (yyval.ldata).d[5] = (yyvsp[-6].fval);
          (yyval.ldata).d[6] = (yyvsp[-5].fval); (yyval.ldata).d[7] = (yyvsp[-4].fval); (yyval.ldata).d[8] = (yyvsp[-3].fval);
          (yyval.ldata).d[9] = (yyvsp[-2].fval);(yyval.ldata).d[10]= (yyvsp[-1].fval);(yyval.ldata).d[11] = 0; }
#line 6402 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 252:
#line 862 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).lnum = (yyvsp[-13].ival)-1;
          (yyval.ldata).matid = -1; // this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[-12].fval); (yyval.ldata).d[1] = (yyvsp[-11].fval); (yyval.ldata).d[2] = (yyvsp[-10].fval);
          (yyval.ldata).d[3] = (yyvsp[-9].fval); (yyval.ldata).d[4] = (yyvsp[-8].fval); (yyval.ldata).d[5] = (yyvsp[-7].fval);
          (yyval.ldata).d[6] = (yyvsp[-6].fval); (yyval.ldata).d[7] = (yyvsp[-5].fval); (yyval.ldata).d[8] = (yyvsp[-4].fval);
          (yyval.ldata).d[9] = (yyvsp[-3].fval);(yyval.ldata).d[10]= (yyvsp[-2].fval); (yyval.ldata).d[11] = (yyvsp[-1].fval); }
#line 6413 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 253:
#line 871 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).lnum = (yyvsp[-4].ival)-1;  (yyval.ldata).matid = (yyvsp[-3].ival)-1; (yyval.ldata).d[7] = (yyvsp[-2].fval); (yyval.ldata).d[8] = (yyvsp[-1].fval); }
#line 6419 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 255:
#line 876 "p.y" /* yacc.c:1646  */
    { geoSource->addLayMat((yyvsp[0].ldata).matid, (yyvsp[0].ldata).d); }
#line 6425 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 256:
#line 882 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).matid = (yyvsp[-6].ival)-1; (yyval.ldata).d[0] = (yyvsp[-5].fval); (yyval.ldata).d[1] = (yyvsp[-4].fval); (yyval.ldata).d[2] = (yyvsp[-3].fval);
          (yyval.ldata).d[3] = (yyvsp[-2].fval); (yyval.ldata).d[4] = 0.0; (yyval.ldata).d[5] = 0.0; (yyval.ldata).d[6] = (yyvsp[-1].fval); 
          (yyval.ldata).d[7] = 0; (yyval.ldata).d[8] = 0; (yyval.ldata).d[9] = 0; }
#line 6433 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 257:
#line 887 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).matid = (yyvsp[-8].ival)-1; (yyval.ldata).d[0] = (yyvsp[-7].fval); (yyval.ldata).d[1] = (yyvsp[-6].fval); (yyval.ldata).d[2] = (yyvsp[-5].fval);
          (yyval.ldata).d[3] = (yyvsp[-4].fval); (yyval.ldata).d[4] = (yyvsp[-3].fval); (yyval.ldata).d[5] = (yyvsp[-2].fval); (yyval.ldata).d[6] = (yyvsp[-1].fval);
          (yyval.ldata).d[7] = 0; (yyval.ldata).d[8] = 0; (yyval.ldata).d[9] = 0; }
#line 6441 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 258:
#line 892 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).matid = (yyvsp[-10].ival)-1; (yyval.ldata).d[0] = (yyvsp[-9].fval); (yyval.ldata).d[1] = (yyvsp[-8].fval); (yyval.ldata).d[2] = (yyvsp[-7].fval);
          (yyval.ldata).d[3] = (yyvsp[-6].fval); (yyval.ldata).d[4] = (yyvsp[-5].fval); (yyval.ldata).d[5] = (yyvsp[-4].fval); (yyval.ldata).d[6] = (yyvsp[-3].fval);
          (yyval.ldata).d[7] = (yyvsp[-2].fval); (yyval.ldata).d[8] = (yyvsp[-1].fval); (yyval.ldata).d[9] = 0; }
#line 6449 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 259:
#line 896 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).matid = (yyvsp[-11].ival)-1; (yyval.ldata).d[0] = (yyvsp[-10].fval); (yyval.ldata).d[1] = (yyvsp[-9].fval); (yyval.ldata).d[2] = (yyvsp[-8].fval);
          (yyval.ldata).d[3] = (yyvsp[-7].fval); (yyval.ldata).d[4] = (yyvsp[-6].fval); (yyval.ldata).d[5] = (yyvsp[-5].fval); (yyval.ldata).d[6] = (yyvsp[-4].fval); 
          (yyval.ldata).d[7] = (yyvsp[-3].fval); (yyval.ldata).d[8] = (yyvsp[-2].fval); (yyval.ldata).d[9] = (yyvsp[-1].fval); }
#line 6457 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 261:
#line 903 "p.y" /* yacc.c:1646  */
    { domain->addDMass((yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1,(yyvsp[-1].fval)); }
#line 6463 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 262:
#line 905 "p.y" /* yacc.c:1646  */
    { domain->addDMass((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,(yyvsp[-1].fval),(yyvsp[-2].ival)-1); }
#line 6469 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 263:
#line 907 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modalDIMASS = true;
          domain->solInfo().reducedMassFile = (yyvsp[-1].strval); }
#line 6476 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 265:
#line 913 "p.y" /* yacc.c:1646  */
    { domain->setGravity((yyvsp[-3].fval),(yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 6482 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 267:
#line 918 "p.y" /* yacc.c:1646  */
    { geoSource->getCheckFileInfo()->lastRestartFile = (yyvsp[-3].strval);
          geoSource->getCheckFileInfo()->outputExt = (yyvsp[-2].strval);
          geoSource->getCheckFileInfo()->FlagRST = (yyvsp[-1].strval); }
#line 6490 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 268:
#line 922 "p.y" /* yacc.c:1646  */
    { geoSource->getCheckFileInfo()->lastRestartFile = (yyvsp[-2].strval);
          geoSource->getCheckFileInfo()->outputExt = (yyvsp[-1].strval);}
#line 6497 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 269:
#line 925 "p.y" /* yacc.c:1646  */
    { geoSource->getCheckFileInfo()->currentRestartFile = (yyvsp[-2].strval);
          domain->solInfo().nRestart = (yyvsp[-1].ival); }
#line 6504 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 270:
#line 930 "p.y" /* yacc.c:1646  */
    { geoSource->setControlFile((yyvsp[-1].strval));
         geoSource->setControlRoutine((char *) "controlObj");}
#line 6511 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 271:
#line 935 "p.y" /* yacc.c:1646  */
    { geoSource->setControlRoutine((yyvsp[-1].strval)); }
#line 6517 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 272:
#line 939 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Sensors;
          if(geoSource->setSensorLocations((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 6524 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 273:
#line 944 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Actuators; }
          if(geoSource->setActuatorLocations((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; 
          if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0)            return -1; }
#line 6532 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 274:
#line 950 "p.y" /* yacc.c:1646  */
    { geoSource->binaryInputControlLeft = true;
          for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Usdf; }
          if(geoSource->setUsdfLocation((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1;
          if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0)       return -1; }
#line 6541 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 275:
#line 957 "p.y" /* yacc.c:1646  */
    { geoSource->binaryInputControlLeft = true;
          (yyval.bclist) = new BCList; }
#line 6548 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 276:
#line 960 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].bcval).type = BCond::Usdd; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 6554 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 277:
#line 962 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-4].ival); i<=(yyvsp[-2].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-1].ival)-1, 0., BCond::Usdd); (yyval.bclist)->add(bc); } }
#line 6560 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 278:
#line 964 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-6].ival); i<=(yyvsp[-4].ival); i+=(yyvsp[-2].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-1].ival)-1, 0., BCond::Usdd); (yyval.bclist)->add(bc); } }
#line 6566 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 279:
#line 966 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[0].bcval);
          surf_bc[0].type = BCond::Usdd;
          geoSource->addSurfaceDirichlet(1,surf_bc);
          if(geoSource->getNumSurfaceDirichlet() > 1) delete [] surf_bc; }
#line 6576 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 280:
#line 974 "p.y" /* yacc.c:1646  */
    { numColumns = 3; }
#line 6582 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 281:
#line 976 "p.y" /* yacc.c:1646  */
    { numColumns = 6; }
#line 6588 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 282:
#line 978 "p.y" /* yacc.c:1646  */
    { numColumns = 3; geoSource->setOutLimit((yyvsp[-1].ival)); }
#line 6594 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 283:
#line 980 "p.y" /* yacc.c:1646  */
    { numColumns = 6; geoSource->setOutLimit((yyvsp[-1].ival)); }
#line 6600 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 284:
#line 982 "p.y" /* yacc.c:1646  */
    { numColumns = 3; domain->outFlag = (yyvsp[-1].ival); }
#line 6606 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 285:
#line 984 "p.y" /* yacc.c:1646  */
    { numColumns = 6; domain->outFlag = (yyvsp[-1].ival); }
#line 6612 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 286:
#line 986 "p.y" /* yacc.c:1646  */
    { numColumns = 3; domain->outFlag = (yyvsp[-2].ival); geoSource->setOutLimit((yyvsp[-1].ival)); }
#line 6618 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 287:
#line 988 "p.y" /* yacc.c:1646  */
    { numColumns = 6; domain->outFlag = (yyvsp[-2].ival); geoSource->setOutLimit((yyvsp[-1].ival)); }
#line 6624 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 288:
#line 990 "p.y" /* yacc.c:1646  */
    { numColumns = 3; geoSource->getCheckFileInfo()->outputExt = (yyvsp[-1].strval); }
#line 6630 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 289:
#line 992 "p.y" /* yacc.c:1646  */
    { numColumns = 6; geoSource->getCheckFileInfo()->outputExt = (yyvsp[-1].strval); }
#line 6636 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 290:
#line 994 "p.y" /* yacc.c:1646  */
    { (yyvsp[-1].oinfo).finalize(numColumns); geoSource->addOutput((yyvsp[-1].oinfo)); }
#line 6642 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 291:
#line 998 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-2].ival); (yyval.oinfo).filename = (yyvsp[-1].strval); (yyval.oinfo).interval = (yyvsp[0].ival); }
#line 6648 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 292:
#line 1000 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-4].ival); (yyval.oinfo).width = (yyvsp[-3].ival); (yyval.oinfo).precision = (yyvsp[-2].ival); (yyval.oinfo).filename = (yyvsp[-1].strval); (yyval.oinfo).interval = (yyvsp[0].ival); }
#line 6654 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 293:
#line 1002 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-3].ival); (yyval.oinfo).filename = (yyvsp[-2].strval); (yyval.oinfo).interval = (yyvsp[-1].ival); (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6660 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 294:
#line 1004 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-4].ival); (yyval.oinfo).filename = (yyvsp[-3].strval); (yyval.oinfo).interval = (yyvsp[-2].ival); 
          if ((yyvsp[-1].ival) == OutputInfo::Nodal) (yyval.oinfo).nodeNumber = (yyvsp[0].ival); else (yyval.oinfo).groupNumber = (yyvsp[0].ival)-1;}
#line 6667 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 295:
#line 1007 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-5].ival); (yyval.oinfo).width = (yyvsp[-4].ival); (yyval.oinfo).precision = (yyvsp[-3].ival); (yyval.oinfo).filename = (yyvsp[-2].strval); (yyval.oinfo).interval = (yyvsp[-1].ival); (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6673 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 296:
#line 1009 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-6].ival); (yyval.oinfo).width = (yyvsp[-5].ival); (yyval.oinfo).precision = (yyvsp[-4].ival); (yyval.oinfo).filename = (yyvsp[-3].strval); (yyval.oinfo).interval = (yyvsp[-2].ival); if ((yyvsp[-1].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[0].ival); else (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6679 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 297:
#line 1011 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-2].ival); (yyval.oinfo).filename = (yyvsp[-1].strval); (yyval.oinfo).interval = (yyvsp[0].ival); }
#line 6685 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 298:
#line 1013 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-4].ival); (yyval.oinfo).width = (yyvsp[-3].ival); (yyval.oinfo).precision = (yyvsp[-2].ival); (yyval.oinfo).filename = (yyvsp[-1].strval); (yyval.oinfo).interval = (yyvsp[0].ival); }
#line 6691 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 299:
#line 1015 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-3].ival); (yyval.oinfo).filename = (yyvsp[-2].strval); (yyval.oinfo).interval = (yyvsp[-1].ival); (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6697 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 300:
#line 1017 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-4].ival); (yyval.oinfo).filename = (yyvsp[-3].strval); (yyval.oinfo).interval = (yyvsp[-2].ival); if ((yyvsp[-1].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[0].ival); else (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6703 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 301:
#line 1019 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-5].ival); (yyval.oinfo).width = (yyvsp[-4].ival); (yyval.oinfo).precision = (yyvsp[-3].ival); (yyval.oinfo).filename = (yyvsp[-2].strval); (yyval.oinfo).interval = (yyvsp[-1].ival); (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6709 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 302:
#line 1021 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-6].ival); (yyval.oinfo).width = (yyvsp[-5].ival); (yyval.oinfo).precision = (yyvsp[-4].ival); (yyval.oinfo).filename = (yyvsp[-3].strval); (yyval.oinfo).interval = (yyvsp[-2].ival); if ((yyvsp[-1].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[0].ival); else (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6715 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 303:
#line 1024 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6721 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 304:
#line 1026 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).surface = (yyvsp[0].ival); }
#line 6727 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 305:
#line 1028 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).str_therm_option = (yyvsp[0].ival); }
#line 6733 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 306:
#line 1030 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).ylayer = (yyvsp[-1].fval); (yyval.oinfo).zlayer = (yyvsp[0].fval); }
#line 6739 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 307:
#line 1032 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).averageFlg = (yyvsp[0].ival); }
#line 6745 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 308:
#line 1034 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).complexouttype = (yyvsp[0].ival); }
#line 6751 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 309:
#line 1036 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).complexouttype = (yyvsp[-1].ival); (yyval.oinfo).ncomplexout = (yyvsp[0].ival); }
#line 6757 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 310:
#line 1038 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).angularouttype = (yyvsp[0].ival); }
#line 6763 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 311:
#line 1040 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).rotvecouttype = (yyvsp[0].ival); }
#line 6769 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 312:
#line 1042 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).rotvecouttype = OutputInfo::Linear; }
#line 6775 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 313:
#line 1044 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).rescaling = bool((yyvsp[0].ival)); }
#line 6781 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 314:
#line 1046 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).ndtype = (yyvsp[0].ival); }
#line 6787 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 315:
#line 1048 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).ndtype = (yyvsp[-1].ival); sfem->setnsamp_out((yyvsp[0].ival)); }
#line 6793 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 316:
#line 1050 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).oframe = (OutputInfo::FrameType) (yyvsp[0].ival); }
#line 6799 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 317:
#line 1052 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).matlab = true; }
#line 6805 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 318:
#line 1054 "p.y" /* yacc.c:1646  */
    { domain->solInfo().xmatrixname = (yyvsp[0].strval); }
#line 6811 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 319:
#line 1056 "p.y" /* yacc.c:1646  */
    { domain->solInfo().qmatrixname = (yyvsp[0].strval); }
#line 6817 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 320:
#line 1058 "p.y" /* yacc.c:1646  */
    { domain->solInfo().rmatrixname = (yyvsp[0].strval); }
#line 6823 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 321:
#line 1060 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenvaluename = (yyvsp[0].strval); }
#line 6829 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 323:
#line 1065 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Modal);
          domain->solInfo().eigenSolverType = SolverInfo::SubSpace; }
#line 6836 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 324:
#line 1068 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Modal);
	  domain->solInfo().nEig = (yyvsp[-1].ival);}
#line 6843 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 325:
#line 1071 "p.y" /* yacc.c:1646  */
    { domain->solInfo().qrfactorization = (yyvsp[-1].ival);}
#line 6849 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 326:
#line 1073 "p.y" /* yacc.c:1646  */
    { domain->solInfo().nEig = (yyvsp[-1].ival); }
#line 6855 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 327:
#line 1075 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::SubSpace;}
#line 6861 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 328:
#line 1077 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setSubSpaceInfo((yyvsp[-3].ival),(yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 6867 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 329:
#line 1079 "p.y" /* yacc.c:1646  */
    { domain->solInfo().subspaceSize = (yyvsp[-1].ival);}
#line 6873 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 330:
#line 1081 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tolEig = (yyvsp[-1].fval); }
#line 6879 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 331:
#line 1083 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tolJac = (yyvsp[-1].fval); }
#line 6885 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 332:
#line 1085 "p.y" /* yacc.c:1646  */
    { domain->solInfo().explicitK = true; }
#line 6891 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 333:
#line 1087 "p.y" /* yacc.c:1646  */
    { geoSource->setShift((yyvsp[-1].fval)); }
#line 6897 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 334:
#line 1089 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack; }
#line 6903 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 335:
#line 1091 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = (yyvsp[-1].strval); }
#line 6910 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 336:
#line 1094 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = (yyvsp[-2].strval); 
          domain->solInfo().arpack_mode = (yyvsp[-1].ival); }
#line 6918 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 337:
#line 1098 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValue((yyvsp[-2].fval), int((yyvsp[-1].fval))); }
#line 6925 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 338:
#line 1101 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValues((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].ival));}
#line 6932 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 339:
#line 1104 "p.y" /* yacc.c:1646  */
    { domain->solInfo().filtereig = bool((yyvsp[-1].ival)); }
#line 6938 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 340:
#line 1106 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverSubType = (yyvsp[-1].ival); }
#line 6944 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 341:
#line 1108 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::LobPcg;
          domain->solInfo().explicitK = true;}
#line 6951 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 342:
#line 1111 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxitEig = (yyvsp[-1].ival); }
#line 6957 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 343:
#line 1113 "p.y" /* yacc.c:1646  */
    { domain->solInfo().test_ulrich = true; }
#line 6963 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 344:
#line 1115 "p.y" /* yacc.c:1646  */
    { domain->solInfo().addedMass = (yyvsp[-1].ival); }
#line 6969 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 345:
#line 1119 "p.y" /* yacc.c:1646  */
    { domain->solInfo().printMatLab = true;
          domain->solInfo().printMatLabFile = (yyvsp[-1].strval);
          domain->solInfo().printMatLabExit = true; }
#line 6977 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 346:
#line 1123 "p.y" /* yacc.c:1646  */
    { domain->solInfo().printMatLab = true;
          domain->solInfo().printMatLabFile = (yyvsp[-2].strval); 
          domain->solInfo().printMatLabExit = true; }
#line 6985 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 347:
#line 1129 "p.y" /* yacc.c:1646  */
    { domain->solInfo().elementDeletion = true; }
#line 6991 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 349:
#line 1134 "p.y" /* yacc.c:1646  */
    { (yyval.fval) = (yyvsp[-1].fval); domain->solInfo().deleteElements.insert(std::pair<int,double>((yyvsp[0].ival)-1,(yyvsp[-1].fval))); }
#line 6997 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 350:
#line 1136 "p.y" /* yacc.c:1646  */
    { domain->solInfo().deleteElements.insert(std::pair<int,double>((yyvsp[0].ival)-1,(yyvsp[-1].fval))); }
#line 7003 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 351:
#line 1140 "p.y" /* yacc.c:1646  */
    { domain->solInfo().sloshing = 1; }
#line 7009 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 352:
#line 1142 "p.y" /* yacc.c:1646  */
    { domain->setGravitySloshing((yyvsp[-1].fval)); }
#line 7015 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 353:
#line 1146 "p.y" /* yacc.c:1646  */
    { domain->solInfo().massFlag = 1; }
#line 7021 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 354:
#line 1148 "p.y" /* yacc.c:1646  */
    { domain->solInfo().massFlag = 1;
          domain->solInfo().massFile = std::string((yyvsp[-1].strval)); }
#line 7028 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 355:
#line 1153 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::ConditionNumber); 
	  domain->solInfo().setCondNumTol((yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 7035 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 356:
#line 1156 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::ConditionNumber);}
#line 7041 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 357:
#line 1160 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Top); }
#line 7047 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 358:
#line 1164 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal_id.push_back((yyvsp[0].ival)); }
#line 7053 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 359:
#line 1166 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal_id.push_back((yyvsp[0].ival)); }
#line 7059 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 360:
#line 1168 "p.y" /* yacc.c:1646  */
    { domain->solInfo().contact_modal_id.push_back((yyvsp[0].ival)); }
#line 7065 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 361:
#line 1171 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Dynamic); }
#line 7071 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 365:
#line 1176 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal = true; domain->solInfo().modal_id.push_back(0); }
#line 7077 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 366:
#line 1178 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal = true; }
#line 7083 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 367:
#line 1180 "p.y" /* yacc.c:1646  */
    { domain->solInfo().stable = (yyvsp[-1].ival); }
#line 7089 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 368:
#line 1182 "p.y" /* yacc.c:1646  */
    { domain->solInfo().stable = (yyvsp[-1].ival); }
#line 7095 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 369:
#line 1184 "p.y" /* yacc.c:1646  */
    { domain->solInfo().stable = (yyvsp[-5].ival);
          domain->solInfo().stable_cfl = (yyvsp[-4].fval);
          domain->solInfo().stable_tol = (yyvsp[-3].fval);
          domain->solInfo().stable_maxit = (yyvsp[-2].ival);
          domain->solInfo().stable_freq = (yyvsp[-1].ival);
        }
#line 7106 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 370:
#line 1191 "p.y" /* yacc.c:1646  */
    { domain->solInfo().iacc_switch = bool((yyvsp[-1].ival)); }
#line 7112 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 371:
#line 1193 "p.y" /* yacc.c:1646  */
    { domain->solInfo().zeroRot = bool((yyvsp[-1].ival)); }
#line 7118 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 372:
#line 1195 "p.y" /* yacc.c:1646  */
    { domain->solInfo().no_secondary = true; }
#line 7124 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 373:
#line 1197 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tdenforceFlag = (yyvsp[-1].ival); }
#line 7130 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 374:
#line 1199 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tdenforceFlag = (yyvsp[-3].ival);
          domain->solInfo().tdenforceMaxItr = (yyvsp[-2].ival);
          domain->solInfo().tdenforceTolAbs = (yyvsp[-1].fval); }
#line 7138 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 375:
#line 1203 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tdenforceFlag = (yyvsp[-4].ival);
          domain->solInfo().tdenforceMaxItr = (yyvsp[-3].ival);
          domain->solInfo().tdenforceTolAbs = (yyvsp[-2].fval);
          domain->solInfo().tdenforceInitia = (yyvsp[-1].fval); }
#line 7147 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 376:
#line 1208 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tdenforceFlag = (yyvsp[-5].ival);
          domain->solInfo().tdenforceMaxItr = (yyvsp[-4].ival);
          domain->solInfo().tdenforceTolAbs = (yyvsp[-3].fval);
          domain->solInfo().tdenforceInitia = (yyvsp[-2].fval);
          domain->solInfo().tdenforceFinal = (yyvsp[-1].fval); }
#line 7157 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 377:
#line 1214 "p.y" /* yacc.c:1646  */
    { domain->solInfo().check_energy_balance = true; }
#line 7163 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 378:
#line 1216 "p.y" /* yacc.c:1646  */
    { domain->solInfo().check_energy_balance = true;
          domain->solInfo().epsilon1 = (yyvsp[-2].fval); 
          domain->solInfo().epsilon2 = (yyvsp[-1].fval); }
#line 7171 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 379:
#line 1220 "p.y" /* yacc.c:1646  */
    { domain->solInfo().quasistatic = bool((yyvsp[-1].ival)); }
#line 7177 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 380:
#line 1224 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ConwepOnOff = true;
          BlastLoading::InputFileData = (yyvsp[-1].blastData); 
        }
#line 7185 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 381:
#line 1229 "p.y" /* yacc.c:1646  */
    { // Note: chargeWeight must be entered in the units of mass of the problem, not units of force.
          (yyval.blastData).ExplosivePosition[0] = (yyvsp[-6].fval);
          (yyval.blastData).ExplosivePosition[1] = (yyvsp[-5].fval);
          (yyval.blastData).ExplosivePosition[2] = (yyvsp[-4].fval);
          (yyval.blastData).ExplosiveWeight = (yyvsp[-3].fval);
          (yyval.blastData).ExplosiveDetonationTime = (yyvsp[-2].fval);
          ((yyvsp[-1].ival) == 0 ? BlastLoading::BlastData::SurfaceBurst : BlastLoading::BlastData::AirBurst);
          (yyval.blastData).UnitConversionId = (yyvsp[0].ival);
        }
#line 7199 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 382:
#line 1240 "p.y" /* yacc.c:1646  */
    { domain->solInfo().timeIntegration = SolverInfo::Newmark; }
#line 7205 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 384:
#line 1243 "p.y" /* yacc.c:1646  */
    { domain->solInfo().acoustic = true; }
#line 7211 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 386:
#line 1248 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[-1].fval),(yyvsp[0].fval)); }
#line 7217 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 387:
#line 1250 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[-3].fval),(yyvsp[-2].fval),(yyvsp[-1].fval),(yyvsp[0].fval)); }
#line 7223 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 388:
#line 1252 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setNewmarkSecondOrderInfo(0.0,0.0,10.0,10.0,(yyvsp[0].fval)); }
#line 7229 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 389:
#line 1254 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[-2].fval),(yyvsp[-1].fval));
          domain->solInfo().modifiedWaveEquation = true;
          domain->solInfo().modifiedWaveEquationCoef = (yyvsp[0].fval); }
#line 7237 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 390:
#line 1260 "p.y" /* yacc.c:1646  */
    { 
          if(domain->solInfo().probType == SolverInfo::NonLinDynam) {
            domain->solInfo().order = 1;
          }
          else 
            domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setNewmarkFirstOrderInfo((yyvsp[0].fval)); 
        }
#line 7250 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 391:
#line 1271 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Dynamic); 
          domain->solInfo().timeIntegration = SolverInfo::Qstatic; }
#line 7257 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 394:
#line 1276 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal = true; domain->solInfo().modal_id.push_back(0); }
#line 7263 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 395:
#line 1278 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal = true; domain->solInfo().modal_id.push_back((yyvsp[-1].ival)); }
#line 7269 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 396:
#line 1282 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setQuasistaticInfo((yyvsp[-3].fval), 0, (yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 7275 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 397:
#line 1284 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setQuasistaticInfo((yyvsp[-4].fval), 0, (yyvsp[-3].fval), (yyvsp[-2].ival), (yyvsp[-1].fval)); }
#line 7281 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 398:
#line 1288 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setQuasistaticInfo((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 7288 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 399:
#line 1298 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setAero((yyvsp[-1].ival)); 
          domain->solInfo().isCollocated = 0; }
#line 7295 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 400:
#line 1301 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setAero((yyvsp[-1].ival)); 
          domain->solInfo().isCollocated = 0;
          if((yyvsp[-1].ival) == 20 || (yyvsp[-1].ival) == 22) { // set default alphas for C0
            domain->solInfo().alphas[0] = 0.5+0.375;
            domain->solInfo().alphas[1] = -0.375;
            domain->solInfo().alphasv    = 0.0;
          }
        }
#line 7308 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 401:
#line 1310 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setAero((yyvsp[-3].ival));
          domain->solInfo().isCollocated = 0;
          if((yyvsp[-3].ival) == 8) {
            // MPP uses only the first of the two inputted alphas
            domain->solInfo().mppFactor = (yyvsp[-2].fval);
          }
          else {
            // These alphas are used in FlExchanger::sendDisplacements and DistFlExchanger::sendDisplacements
            // As of 4/14/2014 the following schemes can use the displacement predictor on the structure side:
            // A0, A4, A5, A6, A7 and C0. Furthermore, we now apply a separate "anti-predictor" for A6 and A7
            // to compensate for the legacy predictor on the fluid side (in MatchNodeSet::getDisplacement)
            domain->solInfo().alphas[0] = (yyvsp[-2].fval)+(yyvsp[-1].fval);
            domain->solInfo().alphas[1] = -(yyvsp[-1].fval);
            domain->solInfo().alphasv    = 0.0;
          }
        }
#line 7329 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 402:
#line 1327 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setAero((yyvsp[-4].ival));
          domain->solInfo().isCollocated = 0;
          if((yyvsp[-4].ival) == 8) {
            // MPP uses only the first of the two inputted alphas
            domain->solInfo().mppFactor = (yyvsp[-3].fval);
          }
          else {
            // These alphas are used in FlExchanger::sendDisplacements and DistFlExchanger::sendDisplacements
            // As of 4/14/2014 the following schemes can use the displacement predictor on the structure side:
            // A0, A4, A5, A6, A7 and C0. Furthermore, we now apply a separate "anti-predictor" for A6 and A7
            // to compensate for the legacy predictor on the fluid side (in MatchNodeSet::getDisplacement)
            domain->solInfo().alphas[0] = (yyvsp[-3].fval)+(yyvsp[-2].fval);
            domain->solInfo().alphas[1] = -(yyvsp[-2].fval);

            // This option added by Alex Main.  The last number on this line is the value used for velocity prediciton,
            // which is useful for embedded simulations.
            domain->solInfo().alphasv    = (yyvsp[-1].fval);
          }
        }
#line 7353 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 403:
#line 1348 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setAero((yyvsp[-2].ival));
          domain->solInfo().isCollocated = 0;
          domain->solInfo().mppFactor = (yyvsp[-1].fval);
        }
#line 7362 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 404:
#line 1353 "p.y" /* yacc.c:1646  */
    { domain->solInfo().isCollocated = (yyvsp[-1].ival); }
#line 7368 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 405:
#line 1355 "p.y" /* yacc.c:1646  */
    { domain->solInfo().subcycle = (yyvsp[-1].ival); }
#line 7374 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 406:
#line 1357 "p.y" /* yacc.c:1646  */
    { geoSource->setMatch((yyvsp[-1].strval)); }
#line 7380 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 407:
#line 1359 "p.y" /* yacc.c:1646  */
    {}
#line 7386 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 408:
#line 1363 "p.y" /* yacc.c:1646  */
    {}
#line 7392 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 409:
#line 1365 "p.y" /* yacc.c:1646  */
    { domain->AddAeroEmbedSurfaceId((yyvsp[0].ival)); }
#line 7398 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 410:
#line 1369 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setAeroHeat((yyvsp[-3].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 7404 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 411:
#line 1373 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setThermoh(1); }
#line 7410 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 412:
#line 1377 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setThermoe(1); }
#line 7416 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 413:
#line 1381 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setModeDecomp(1); }
#line 7422 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 414:
#line 1383 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setModeDecomp(1, (yyvsp[-1].ival)); }
#line 7428 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 415:
#line 1387 "p.y" /* yacc.c:1646  */
    { domain->solInfo().hzemFlag=1; }
#line 7434 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 416:
#line 1391 "p.y" /* yacc.c:1646  */
    { domain->solInfo().slzemFlag=1; }
#line 7440 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 417:
#line 1395 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setTrbm((yyvsp[-1].fval)); }
#line 7446 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 418:
#line 1399 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 7452 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 419:
#line 1401 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-1].fval)); }
#line 7458 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 420:
#line 1403 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm(); }
#line 7464 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 421:
#line 1405 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-3].fval),(yyvsp[-2].fval));
          domain->solInfo().grbm_use_lmpc = bool((yyvsp[-1].ival)); }
#line 7471 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 422:
#line 1408 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-2].fval),(yyvsp[-1].fval));
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = (yyvsp[-6].fval); grbm_ref[1] = (yyvsp[-5].fval); grbm_ref[2] = (yyvsp[-4].fval);
        }
#line 7480 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 423:
#line 1413 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-1].fval));
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = (yyvsp[-5].fval); grbm_ref[1] = (yyvsp[-4].fval); grbm_ref[2] = (yyvsp[-3].fval);
        }
#line 7489 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 424:
#line 1418 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm();
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = (yyvsp[-3].fval); grbm_ref[1] = (yyvsp[-2].fval); grbm_ref[2] = (yyvsp[-1].fval);
        }
#line 7498 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 425:
#line 1423 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-3].fval),(yyvsp[-2].fval));
          domain->solInfo().grbm_use_lmpc = bool((yyvsp[-1].ival));
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = (yyvsp[-7].fval); grbm_ref[1] = (yyvsp[-6].fval); grbm_ref[2] = (yyvsp[-5].fval);
        }
#line 7508 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 426:
#line 1431 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modeFilterFlag = (yyvsp[-1].ival); }
#line 7514 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 427:
#line 1433 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modeFilterFlag = 1; }
#line 7520 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 428:
#line 1437 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useRbmFilter((yyvsp[-1].ival)); }
#line 7526 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 429:
#line 1439 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useRbmFilter(1); }
#line 7532 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 430:
#line 1441 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useRbmFilter((yyvsp[-2].ival));
          domain->solInfo().filterQ = (yyvsp[-1].ival); }
#line 7539 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 432:
#line 1447 "p.y" /* yacc.c:1646  */
    { if((yyvsp[0].ival) < 1) {
        fprintf(stderr, " *** ERROR: mode %d specified under RBMFILTER is invalid.\n", (yyvsp[0].ival));
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters.insert((yyvsp[0].ival)-1);
    }
#line 7551 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 433:
#line 1455 "p.y" /* yacc.c:1646  */
    { if((yyvsp[0].ival) < 1) {
        fprintf(stderr, " *** ERROR: mode %d specified under RBMFILTER is invalid.\n", (yyvsp[0].ival));
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters.insert((yyvsp[0].ival)-1);
    }
#line 7563 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 434:
#line 1465 "p.y" /* yacc.c:1646  */
    { domain->solInfo().hzemFilterFlag=1; }
#line 7569 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 435:
#line 1469 "p.y" /* yacc.c:1646  */
    { domain->solInfo().slzemFilterFlag=1; }
#line 7575 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 436:
#line 1473 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setTimes((yyvsp[-1].fval),(yyvsp[-2].fval),(yyvsp[-3].fval)); }
#line 7581 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 437:
#line 1477 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime((yyvsp[-2].ival),(yyvsp[-1].ival),1);
        }
#line 7590 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 438:
#line 1483 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime((yyvsp[-3].ival),(yyvsp[-2].ival),(yyvsp[-1].ival));
        }
#line 7599 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 439:
#line 1488 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().mdPita = true;
          domain->solInfo().setParallelInTime((yyvsp[-4].ival),(yyvsp[-3].ival),(yyvsp[-2].ival)); 
          /*domain->solInfo().numSpaceMPIProc = $6;*/
        }
#line 7610 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 442:
#line 1501 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaNoForce = true; }
#line 7616 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 443:
#line 1503 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaGlobalBasisImprovement = (yyvsp[0].ival); }
#line 7622 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 444:
#line 1505 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaLocalBasisImprovement = 1; }
#line 7628 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 445:
#line 1507 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaTimeReversible = true; }
#line 7634 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 446:
#line 1509 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaRemoteCoarse = true; }
#line 7640 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 447:
#line 1511 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaProjTol = (yyvsp[0].fval); }
#line 7646 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 448:
#line 1513 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaReadInitSeed = true; }
#line 7652 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 449:
#line 1515 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaJumpCvgRatio = 0.0; }
#line 7658 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 450:
#line 1517 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaJumpCvgRatio = (yyvsp[0].fval); }
#line 7664 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 451:
#line 1519 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaJumpMagnOutput = true; }
#line 7670 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 452:
#line 1523 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-3].ival) == 1) domain->solInfo().setDamping((yyvsp[-2].fval),(yyvsp[-1].fval));
          else return -1; // only RAYDAMP is allowed here
        }
#line 7678 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 453:
#line 1527 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-4].ival) == 1) {
            domain->solInfo().setDamping((yyvsp[-3].fval),(yyvsp[-2].fval));
            domain->solInfo().mtypeDamp = (int)(yyvsp[-1].ival);
          }
          else return -1;
        }
#line 7689 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 454:
#line 1534 "p.y" /* yacc.c:1646  */
    { if(geoSource->setModalDamping((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true;
        }
#line 7697 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 455:
#line 1540 "p.y" /* yacc.c:1646  */
    { (yyval.cxbclist) = (yyvsp[0].cxbclist); }
#line 7703 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 456:
#line 1544 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 7713 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 457:
#line 1550 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[-2].fval), (yyvsp[-1].fval), 0.0);
        }
#line 7723 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 458:
#line 1558 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->solInfo().setProbType(SolverInfo::HelmholtzDirSweep);
           domain->setWaveDirections((yyvsp[-3].ival), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 7733 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 459:
#line 1564 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections((yyvsp[-1].ival),0.0,0.0,0.0);
        }
#line 7742 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 461:
#line 1572 "p.y" /* yacc.c:1646  */
    {
           domain->setWaveDirections((yyvsp[-4].ival), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 7750 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 462:
#line 1578 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 7760 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 463:
#line 1586 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; }
#line 7766 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 464:
#line 1588 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].bcval).type = BCond::Displacements; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 7772 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 465:
#line 1590 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-5].ival); i<=(yyvsp[-3].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval), BCond::Displacements); (yyval.bclist)->add(bc); } }
#line 7778 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 466:
#line 1592 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-7].ival); i<=(yyvsp[-5].ival); i+=(yyvsp[-3].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval), BCond::Displacements); (yyval.bclist)->add(bc); } }
#line 7784 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 467:
#line 1594 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[0].bcval);
          surf_bc[0].type = BCond::Displacements;
          geoSource->addSurfaceDirichlet(1,surf_bc);
          if(geoSource->getNumSurfaceDirichlet() > 1) delete [] surf_bc; }
#line 7794 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 469:
#line 1610 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[-4].ival)-1;
          surf_bc[0].type = (BCond::BCType) (yyvsp[-3].ival); //BCond::PointPlaneDistance;
          surf_bc[0].dofnum = (yyvsp[-2].ival)-1;
          surf_bc[0].val = (yyvsp[-1].ival)-1;
          geoSource->addSurfaceConstraint(1,surf_bc);
          if(geoSource->getNumSurfaceConstraint() > 1) delete [] surf_bc;
        }
#line 7807 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 470:
#line 1620 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Pdir; (yyval.bclist) = (yyvsp[0].bclist); }
#line 7813 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 471:
#line 1624 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 7819 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 472:
#line 1626 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 7825 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 473:
#line 1630 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-2].ival)-1; (yyval.bcval).dofnum = 10; (yyval.bcval).val = (yyvsp[-1].fval); }
#line 7831 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 474:
#line 1632 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-1].ival)-1; (yyval.bcval).dofnum = 10; (yyval.bcval).val = 0.0; }
#line 7837 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 475:
#line 1636 "p.y" /* yacc.c:1646  */
    { domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; 

          int* allPDirNodes = new int[(yyvsp[0].bclist)->n];

          for (int ii=0; ii < (yyvsp[0].bclist)->n; ii++)
            allPDirNodes[ii]=((yyvsp[0].bclist)->d[ii]).nnum;
          std::sort(allPDirNodes,allPDirNodes + (yyvsp[0].bclist)->n);

          int maxFSNodes = 32;
          int* allHEVFSNodes = new int[maxFSNodes];
          allHEVFSNodes[0] = allPDirNodes[0];

          int numHEVFSNodes = 1;
          for (int ii = 1; ii < (yyvsp[0].bclist)->n; ++ii)  {
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
          
          (yyval.bclist) = new BCList;

          for (int ii=0; ii<numHEVFSNodes; ii++)
            (yyval.bclist)->add(allHEVFSNodes[ii],10,0.0); 
          delete [] allHEVFSNodes;

          for(int i=0; i<(yyval.bclist)->n; ++i) (yyval.bclist)->d[i].type = BCond::Hefrs;
        }
#line 7880 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 476:
#line 1677 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; 
          for (int ii = 0; ii < (yyvsp[0].bclist)->n; ii++) 
           (yyval.bclist)->add(((yyvsp[0].bclist)->d)[ii]); }
#line 7888 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 477:
#line 1681 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); 
          for (int ii = 0; ii < (yyvsp[0].bclist)->n; ii++) 
           (yyval.bclist)->add(((yyvsp[0].bclist)->d)[ii]); }
#line 7896 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 478:
#line 1687 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList;
          for(int i=0; i<(yyvsp[-1].nl).num; ++i) 
          { (yyval.bclist)->add((yyvsp[-1].nl).nd[i],10,0.0); } }
#line 7904 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 479:
#line 1693 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; if(domain->solInfo().soltyp != 2) domain->solInfo().thermalLoadFlag = 1;}
#line 7910 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 480:
#line 1695 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-3].bclist); BCond bc; bc.nnum = (yyvsp[-2].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[-1].fval); bc.type = BCond::Temperatures; (yyval.bclist)->add(bc); }
#line 7917 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 481:
#line 1698 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-4].ival); i<=(yyvsp[-2].ival); ++i) { BCond bc; bc.setData(i-1, 6, (yyvsp[-1].fval), BCond::Temperatures); (yyval.bclist)->add(bc); } }
#line 7923 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 482:
#line 1700 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-6].ival); i<=(yyvsp[-4].ival); i+=(yyvsp[-2].ival)) { BCond bc; bc.setData(i-1, 6, (yyvsp[-1].fval), BCond::Temperatures); (yyval.bclist)->add(bc); } }
#line 7929 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 483:
#line 1702 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[-2].ival)-1;
          surf_bc[0].val = (yyvsp[-1].fval);
          surf_bc[0].dofnum = 6;
          surf_bc[0].type = BCond::Temperatures;
          geoSource->addSurfaceDirichlet(1,surf_bc); }
#line 7940 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 484:
#line 1711 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; }
#line 7946 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 485:
#line 1713 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList((yyvsp[-1].ival)); }
#line 7952 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 486:
#line 1715 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-3].bclist); BCond bc; bc.nnum = (yyvsp[-2].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[-1].fval); bc.type = BCond::Flux; bc.loadsetid = (yyval.bclist)->loadsetid; (yyval.bclist)->add(bc); }
#line 7959 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 487:
#line 1718 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[-2].ival)-1;
          surf_bc[0].dofnum = 6;
          surf_bc[0].val = (yyvsp[-1].fval);
          surf_bc[0].type = BCond::Flux;
          surf_bc[0].loadsetid = (yyval.bclist)->loadsetid;
          geoSource->addSurfaceNeuman(1,surf_bc);
          if(geoSource->getNumSurfaceNeuman() > 1) delete [] surf_bc; }
#line 7972 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 488:
#line 1729 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; }
#line 7978 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 489:
#line 1731 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList((yyvsp[-1].ival)); }
#line 7984 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 490:
#line 1733 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-5].bclist); BCond bc; bc.nnum = (yyvsp[-4].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[-3].fval)*(yyvsp[-2].fval)*(yyvsp[-1].fval); bc.type = BCond::Convection; bc.loadsetid = (yyval.bclist)->loadsetid; (yyval.bclist)->add(bc); }
#line 7991 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 491:
#line 1736 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[-4].ival)-1;
          surf_bc[0].dofnum = 6;
          surf_bc[0].val = (yyvsp[-3].fval)*(yyvsp[-2].fval)*(yyvsp[-1].fval);
          surf_bc[0].type = BCond::Convection;
          surf_bc[0].loadsetid = (yyval.bclist)->loadsetid;
          geoSource->addSurfaceNeuman(1,surf_bc);
          if(geoSource->getNumSurfaceNeuman() > 1) delete [] surf_bc; }
#line 8004 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 492:
#line 1747 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; }
#line 8010 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 493:
#line 1749 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList((yyvsp[-1].ival)); }
#line 8016 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 494:
#line 1751 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-5].bclist); BCond bc; bc.nnum = (yyvsp[-4].ival)-1; bc.dofnum = 6;
          bc.val = 5.670400E-8*(yyvsp[-3].fval)*(yyvsp[-2].fval)*(yyvsp[-1].fval)*(yyvsp[-1].fval)*(yyvsp[-1].fval)*(yyvsp[-1].fval); bc.type = BCond::Radiation; (yyval.bclist)->add(bc); }
#line 8023 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 498:
#line 1763 "p.y" /* yacc.c:1646  */
    { domain->addSommer(new LineSommerBC((yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1)); }
#line 8029 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 499:
#line 1765 "p.y" /* yacc.c:1646  */
    { domain->addSommer(new TriangleSommerBC((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1)); }
#line 8035 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 500:
#line 1767 "p.y" /* yacc.c:1646  */
    { domain->addSommer(new QuadSommerBC((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1)); }
#line 8041 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 503:
#line 1775 "p.y" /* yacc.c:1646  */
    { domain->addSommerElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd); 
          /*geoSource->addElem($1-1, $2, $3.num, $3.nd);include Sommer nodes in PackedEset -JF*/
        }
#line 8049 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 504:
#line 1781 "p.y" /* yacc.c:1646  */
    { (yyval.nl).num = 1; (yyval.nl).nd[0] = (yyvsp[0].ival)-1; }
#line 8055 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 505:
#line 1783 "p.y" /* yacc.c:1646  */
    { if((yyval.nl).num == 64) return -1;
          (yyval.nl).nd[(yyval.nl).num] = (yyvsp[0].ival)-1; (yyval.nl).num++; }
#line 8062 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 509:
#line 1795 "p.y" /* yacc.c:1646  */
    { domain->addScatter(new LineSommerBC((yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1));
          domain->addNeum(new LineSommerBC((yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1)); }
#line 8069 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 510:
#line 1798 "p.y" /* yacc.c:1646  */
    { domain->addScatter(new TriangleSommerBC((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1));
          domain->addNeum(new TriangleSommerBC((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1)); }
#line 8076 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 511:
#line 1801 "p.y" /* yacc.c:1646  */
    { domain->addScatter(new QuadSommerBC((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1));
          domain->addNeum(new QuadSommerBC((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1)); }
#line 8083 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 514:
#line 1810 "p.y" /* yacc.c:1646  */
    { domain->addScatterElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd);
          domain->addNeumElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd); }
#line 8090 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 517:
#line 1819 "p.y" /* yacc.c:1646  */
    { domain->addNeumElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd); }
#line 8096 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 520:
#line 1827 "p.y" /* yacc.c:1646  */
    { domain->addWetElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd); 
          domain->solInfo().isCoupled = true; 
          domain->solInfo().isMatching = true; }
#line 8104 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 523:
#line 1837 "p.y" /* yacc.c:1646  */
    { domain->addScatterElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd);}
#line 8110 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 524:
#line 1841 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-2].ival)-1; (yyval.bcval).dofnum = 7; (yyval.bcval).val = (yyvsp[-1].fval); }
#line 8116 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 525:
#line 1845 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8122 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 526:
#line 1847 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8128 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 527:
#line 1851 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Atddir; (yyval.bclist) = (yyvsp[0].bclist); }
#line 8134 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 528:
#line 1855 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Atdneu; } (yyval.bclist) = (yyvsp[0].bclist); }
#line 8140 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 529:
#line 1857 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Atdneu; } (yyval.bclist) = (yyvsp[0].bclist); }
#line 8146 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 530:
#line 1861 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ATDARBFlag = (yyvsp[-1].fval);}
#line 8152 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 532:
#line 1866 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ATDDNBVal = (yyvsp[-1].fval);}
#line 8158 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 534:
#line 1871 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ATDROBVal = (yyvsp[-3].fval);
          domain->solInfo().ATDROBalpha = (yyvsp[-2].fval);
          domain->solInfo().ATDROBbeta = (yyvsp[-1].fval);}
#line 8166 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 536:
#line 1878 "p.y" /* yacc.c:1646  */
    { domain->setFFP((yyvsp[-1].ival)); }
#line 8172 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 537:
#line 1880 "p.y" /* yacc.c:1646  */
    { domain->setFFP((yyvsp[-2].ival),(yyvsp[-1].ival)); }
#line 8178 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 538:
#line 1884 "p.y" /* yacc.c:1646  */
    {
           domain->setFFP((yyvsp[-2].ival));
        }
#line 8186 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 539:
#line 1890 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[(yyvsp[-4].ival)] = ModalParams(ModalParams::Eigen, (yyvsp[-2].strval), (yyvsp[-1].ival)); }
#line 8192 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 540:
#line 1892 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[(yyvsp[-4].ival)] = ModalParams((yyvsp[-3].mpt), (yyvsp[-2].strval), (yyvsp[-1].ival)); }
#line 8198 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 541:
#line 1894 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[(yyvsp[-5].ival)] = ModalParams(ModalParams::Eigen, (yyvsp[-3].strval), (yyvsp[-2].ival), (yyvsp[-1].fval)); }
#line 8204 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 542:
#line 1896 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[(yyvsp[-5].ival)] = ModalParams((yyvsp[-4].mpt), (yyvsp[-3].strval), (yyvsp[-2].ival), (yyvsp[-1].fval)); }
#line 8210 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 543:
#line 1900 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[0] = ModalParams(ModalParams::Undefined, (yyvsp[-1].strval)); }
#line 8216 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 544:
#line 1902 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[0] = ModalParams(ModalParams::Undefined, (yyvsp[-2].strval), (yyvsp[-1].ival)); }
#line 8222 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 547:
#line 1906 "p.y" /* yacc.c:1646  */
    { domain->solInfo().adjointMap[(OutputInfo::Type)(yyvsp[-1].ival)] = domain->solInfo().readInAdjointROB.size();
          domain->solInfo().readInAdjointROB.push_back((yyvsp[-3].strval));
          domain->solInfo().maxSizeAdjointBasis.push_back((yyvsp[-2].ival)); }
#line 8230 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 548:
#line 1912 "p.y" /* yacc.c:1646  */
    { }
#line 8236 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 549:
#line 1914 "p.y" /* yacc.c:1646  */
    { domain->solInfo().zeroInitialDisp = 1; }
#line 8242 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 550:
#line 1916 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDis((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0)  return -1; }
#line 8249 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 551:
#line 1919 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDisModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; }
#line 8257 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 552:
#line 1923 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDisModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1;
          domain->solInfo().modalCalled = true;
          domain->solInfo().idis_modal_id = (yyvsp[-2].ival); }
#line 8266 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 553:
#line 1930 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; amplitude = (yyvsp[-1].fval);  }
#line 8272 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 554:
#line 1932 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; amplitude = 1.0; }
#line 8278 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 555:
#line 1934 "p.y" /* yacc.c:1646  */
    { BCond bc; /* add 6 boundary conditions */
          bc.type = BCond::Idisp6;
          bc.nnum = (yyvsp[-7].ival)-1; bc.dofnum = 0; bc.val = amplitude*(yyvsp[-6].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = amplitude*(yyvsp[-5].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = amplitude*(yyvsp[-4].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = amplitude*(yyvsp[-3].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = amplitude*(yyvsp[-2].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = amplitude*(yyvsp[-1].fval); (yyval.bclist)->add(bc);
          geoSource->setIDis6((yyval.bclist)->n, (yyval.bclist)->d);
        }
#line 8293 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 556:
#line 1945 "p.y" /* yacc.c:1646  */
    { BCond bc; /* add 6 boundary conditions */
          bc.type = BCond::Idisp6;
          bc.nnum = (yyvsp[-4].ival)-1; bc.dofnum = 0; bc.val = amplitude*(yyvsp[-3].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = amplitude*(yyvsp[-2].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = amplitude*(yyvsp[-1].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = 0.0         ; (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = 0.0         ; (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = 0.0         ; (yyval.bclist)->add(bc);
          geoSource->setIDis6((yyval.bclist)->n, (yyval.bclist)->d);
        }
#line 8308 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 557:
#line 1956 "p.y" /* yacc.c:1646  */
    { fprintf(stderr," ... Geometric Pre-Stress Effects   ... \n"); 
          domain->solInfo().setGEPS(); }
#line 8315 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 558:
#line 1959 "p.y" /* yacc.c:1646  */
    { domain->solInfo().buckling = 1; }
#line 8321 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 559:
#line 1964 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; PitaTS = (yyvsp[-1].ival); }
#line 8327 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 560:
#line 1966 "p.y" /* yacc.c:1646  */
    { BCond bc;                          /* add 6 boundary conditions */
          bc.nnum = (yyvsp[-7].ival)-1; bc.dofnum = 0; bc.val = (yyvsp[-6].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = (yyvsp[-5].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = (yyvsp[-4].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = (yyvsp[-3].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = (yyvsp[-2].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = (yyvsp[-1].fval); (yyval.bclist)->add(bc);
          geoSource->setPitaIDis6((yyval.bclist)->n, (yyval.bclist)->d, PitaTS);
        }
#line 8341 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 561:
#line 1979 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; PitaTS = (yyvsp[-1].ival); }
#line 8347 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 562:
#line 1981 "p.y" /* yacc.c:1646  */
    { BCond bc;                          /* add 6 boundary conditions */
          bc.nnum = (yyvsp[-7].ival)-1; bc.dofnum = 0; bc.val = (yyvsp[-6].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = (yyvsp[-5].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = (yyvsp[-4].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = (yyvsp[-3].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = (yyvsp[-2].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = (yyvsp[-1].fval); (yyval.bclist)->add(bc);
          geoSource->setPitaIVel6((yyval.bclist)->n, (yyval.bclist)->d, PitaTS);
        }
#line 8361 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 563:
#line 1993 "p.y" /* yacc.c:1646  */
    { }
#line 8367 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 564:
#line 1995 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVel((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 8374 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 565:
#line 1998 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVelModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; }
#line 8382 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 566:
#line 2002 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVelModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1;
          domain->solInfo().modalCalled = true;
          domain->solInfo().ivel_modal_id = (yyvsp[-2].ival); }
#line 8391 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 567:
#line 2009 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Itemperatures;
          if(geoSource->setIDis((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 8398 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 568:
#line 2012 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Itemperatures;
          if(geoSource->setIDisModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1;
          domain->solInfo().modalCalled = true; }
#line 8406 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 569:
#line 2018 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGEPS();
          for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Etemperatures;
          if(geoSource->setIDis6((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 8414 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 570:
#line 2024 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; }
#line 8420 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 571:
#line 2026 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList((yyvsp[-1].ival)); }
#line 8426 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 572:
#line 2028 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].bcval).type = BCond::Forces; (yyvsp[0].bcval).loadsetid = (yyval.bclist)->loadsetid; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8432 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 573:
#line 2030 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-5].ival); i<=(yyvsp[-3].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval), BCond::Forces, (yyval.bclist)->loadsetid); (yyval.bclist)->add(bc); } }
#line 8438 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 574:
#line 2032 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-7].ival); i<=(yyvsp[-5].ival); i+=(yyvsp[-3].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval), BCond::Forces, (yyval.bclist)->loadsetid); (yyval.bclist)->add(bc); } }
#line 8444 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 575:
#line 2034 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-6].ival); i<=(yyvsp[-4].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-3].ival)-1, (yyvsp[-2].fval), BCond::Forces, (yyval.bclist)->loadsetid, (BCond::MomentType) (yyvsp[-1].ival)); (yyval.bclist)->add(bc); } }
#line 8450 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 576:
#line 2036 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-8].ival); i<=(yyvsp[-6].ival); i+=(yyvsp[-4].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-3].ival)-1, (yyvsp[-2].fval), BCond::Forces, (yyval.bclist)->loadsetid, (BCond::MomentType) (yyvsp[-3].ival)); (yyval.bclist)->add(bc); } }
#line 8456 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 577:
#line 2038 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[0].bcval);
          surf_bc[0].type = BCond::Forces;
          surf_bc[0].loadsetid = (yyval.bclist)->loadsetid;
          geoSource->addSurfaceNeuman(1,surf_bc);
          if(geoSource->getNumSurfaceNeuman() > 1) delete [] surf_bc; }
#line 8467 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 578:
#line 2047 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Forces;
          if(geoSource->setNeumanModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; }
#line 8474 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 579:
#line 2050 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Forces; (yyvsp[0].bclist)->d[i].loadsetid = (yyvsp[-4].ival); }
          if(geoSource->setNeumanModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; }
#line 8481 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 580:
#line 2053 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Forces;
          if(geoSource->setNeumanModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; }
#line 8488 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 581:
#line 2056 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Forces; (yyvsp[0].bclist)->d[i].loadsetid = (yyvsp[-5].ival); }
          if(geoSource->setNeumanModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; }
#line 8495 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 582:
#line 2061 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8501 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 583:
#line 2063 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8507 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 584:
#line 2065 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; for(int i=(yyvsp[-5].ival); i<=(yyvsp[-3].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval)); (yyval.bclist)->add(bc); }}
#line 8513 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 585:
#line 2067 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-6].bclist); for(int i=(yyvsp[-5].ival); i<=(yyvsp[-3].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval)); (yyval.bclist)->add(bc); }}
#line 8519 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 586:
#line 2069 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; for(int i=(yyvsp[-7].ival); i<=(yyvsp[-5].ival); i+=(yyvsp[-3].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval)); (yyval.bclist)->add(bc); } }
#line 8525 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 587:
#line 2071 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-8].bclist); for(int i=(yyvsp[-7].ival); i<=(yyvsp[-5].ival); i+=(yyvsp[-3].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval)); (yyval.bclist)->add(bc); } }
#line 8531 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 588:
#line 2075 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8537 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 589:
#line 2077 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8543 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 590:
#line 2081 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8549 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 591:
#line 2083 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8555 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 594:
#line 2091 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYMTT((yyval.ymtt));}
#line 8561 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 595:
#line 2093 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8567 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 596:
#line 2095 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYMTT((yyval.ymtt));}
#line 8573 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 599:
#line 2103 "p.y" /* yacc.c:1646  */
    { (yyval.ctett) = new MFTTData((yyvsp[-4].ival)); (yyval.ctett)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addCTETT((yyval.ctett));}
#line 8579 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 600:
#line 2105 "p.y" /* yacc.c:1646  */
    { (yyval.ctett)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8585 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 601:
#line 2107 "p.y" /* yacc.c:1646  */
    { (yyval.ctett) = new MFTTData((yyvsp[-4].ival)); (yyval.ctett)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addCTETT((yyval.ctett));}
#line 8591 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 604:
#line 2115 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addSS1DT((yyval.ymtt));}
#line 8597 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 605:
#line 2117 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8603 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 606:
#line 2119 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addSS1DT((yyval.ymtt));}
#line 8609 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 609:
#line 2127 "p.y" /* yacc.c:1646  */
    { (yyval.ss2dt) = new SS2DTData((yyvsp[-6].ival), (yyvsp[-4].dlist), true); (yyval.ss2dt)->add((yyvsp[-2].fval), (yyvsp[-1].dlist)); domain->addSS2DT((yyval.ss2dt)); }
#line 8615 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 610:
#line 2129 "p.y" /* yacc.c:1646  */
    { (yyval.ss2dt) = new SS2DTData((yyvsp[-8].ival), (yyvsp[-4].dlist), bool((yyvsp[-6].ival))); (yyval.ss2dt)->add((yyvsp[-2].fval), (yyvsp[-1].dlist)); domain->addSS2DT((yyval.ss2dt)); }
#line 8621 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 611:
#line 2131 "p.y" /* yacc.c:1646  */
    { (yyval.ss2dt)->add((yyvsp[-2].fval), (yyvsp[-1].dlist)); }
#line 8627 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 612:
#line 2133 "p.y" /* yacc.c:1646  */
    { (yyval.ss2dt) = new SS2DTData((yyvsp[-6].ival), (yyvsp[-4].dlist), true); (yyval.ss2dt)->add((yyvsp[-2].fval), (yyvsp[-1].dlist)); domain->addSS2DT((yyval.ss2dt)); }
#line 8633 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 613:
#line 2135 "p.y" /* yacc.c:1646  */
    { (yyval.ss2dt) = new SS2DTData((yyvsp[-8].ival), (yyvsp[-4].dlist), bool((yyvsp[-6].ival))); (yyval.ss2dt)->add((yyvsp[-2].fval), (yyvsp[-1].dlist)); domain->addSS2DT((yyval.ss2dt)); }
#line 8639 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 616:
#line 2143 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYSST((yyval.ymtt));}
#line 8645 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 617:
#line 2145 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8651 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 618:
#line 2147 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYSST((yyval.ymtt));}
#line 8657 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 621:
#line 2155 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYSSRT((yyval.ymtt));}
#line 8663 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 622:
#line 2157 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8669 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 623:
#line 2159 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYSSRT((yyval.ymtt));}
#line 8675 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 626:
#line 2167 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYMST((yyval.ymtt));}
#line 8681 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 627:
#line 2169 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8687 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 628:
#line 2171 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYMST((yyval.ymtt));}
#line 8693 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 631:
#line 2179 "p.y" /* yacc.c:1646  */
    { (yyval.sdetaft) = new MFTTData((yyvsp[-4].ival)); (yyval.sdetaft)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addSDETAFT((yyval.sdetaft));}
#line 8699 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 632:
#line 2181 "p.y" /* yacc.c:1646  */
    { (yyval.sdetaft)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8705 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 633:
#line 2183 "p.y" /* yacc.c:1646  */
    { (yyval.sdetaft) = new MFTTData((yyvsp[-4].ival)); (yyval.sdetaft)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addSDETAFT((yyval.sdetaft));}
#line 8711 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 636:
#line 2191 "p.y" /* yacc.c:1646  */
    { 
#ifdef USE_EIGEN3
          (yyval.rubdaft) = new GenMFTTData<Eigen::Vector4d>((yyvsp[-7].ival)); (yyval.rubdaft)->add((yyvsp[-5].fval), Eigen::Vector4d((yyvsp[-4].fval),(yyvsp[-3].fval),(yyvsp[-2].fval),(yyvsp[-1].fval))); domain->addRUBDAFT((yyval.rubdaft));
#else
          std::cerr << " *** ERROR: RUBDAFT command requires AERO-S configured with Eigen library. Exiting...\n"; exit(-1);
#endif
        }
#line 8723 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 637:
#line 2199 "p.y" /* yacc.c:1646  */
    {
#ifdef USE_EIGEN3
          (yyval.rubdaft)->add((yyvsp[-5].fval), Eigen::Vector4d((yyvsp[-4].fval),(yyvsp[-3].fval),(yyvsp[-2].fval),(yyvsp[-1].fval)));
#else
          std::cerr << " *** ERROR: RUBDAFT command requires AERO-S configured with Eigen library. Exiting...\n"; exit(-1);
#endif
        }
#line 8735 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 638:
#line 2207 "p.y" /* yacc.c:1646  */
    {
#ifdef USE_EIGEN3
          (yyval.rubdaft) = new GenMFTTData<Eigen::Vector4d>((yyvsp[-7].ival)); (yyval.rubdaft)->add((yyvsp[-5].fval), Eigen::Vector4d((yyvsp[-4].fval),(yyvsp[-3].fval),(yyvsp[-2].fval),(yyvsp[-1].fval))); domain->addRUBDAFT((yyval.rubdaft));
#else
          std::cerr << " *** ERROR: RUBDAFT command requires AERO-S configured with Eigen library. Exiting...\n"; exit(-1);
#endif
        }
#line 8747 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 641:
#line 2219 "p.y" /* yacc.c:1646  */
    { domain->solInfo().xLMPCFactor = (yyvsp[-3].fval);
          domain->solInfo().yLMPCFactor = (yyvsp[-2].fval);
          domain->solInfo().zLMPCFactor = (yyvsp[-1].fval); }
#line 8755 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 642:
#line 2225 "p.y" /* yacc.c:1646  */
    { domain->solInfo().localDualBasisSize.push_back((yyvsp[-1].ival));
          domain->solInfo().modalLMPC = true; }
#line 8762 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 643:
#line 2228 "p.y" /* yacc.c:1646  */
    { geoSource->pushBackROMLMPCVec((yyvsp[-1].fval)); }
#line 8768 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 644:
#line 2232 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = (yyvsp[-1].lmpcons);
          (yyval.lmpcons)->addterm((yyvsp[0].mpcterm));
          domain->addLMPC((yyval.lmpcons)); }
#line 8776 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 645:
#line 2236 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons)->addterm((yyvsp[0].mpcterm)); }
#line 8782 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 646:
#line 2238 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = (yyvsp[-1].lmpcons);
          (yyval.lmpcons)->addterm((yyvsp[0].mpcterm));
          domain->addLMPC((yyval.lmpcons)); }
#line 8790 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 647:
#line 2244 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-1].ival), 0.0); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
#line 8797 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 648:
#line 2247 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-2].ival), (yyvsp[-1].fval)); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
#line 8804 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 649:
#line 2250 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-4].ival), (yyvsp[-3].fval));
          (yyval.lmpcons)->type = (yyvsp[-1].ival); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
#line 8812 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 650:
#line 2254 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-3].ival), (yyvsp[-2].fval));
          (yyval.lmpcons)->lagrangeMult = (yyvsp[-1].copt).lagrangeMult;
          (yyval.lmpcons)->penalty = (yyvsp[-1].copt).penalty; 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
#line 8821 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 651:
#line 2259 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-5].ival), (yyvsp[-4].fval));
          (yyval.lmpcons)->type = (yyvsp[-2].ival);
          (yyval.lmpcons)->lagrangeMult = (yyvsp[-1].copt).lagrangeMult;
          (yyval.lmpcons)->penalty = (yyvsp[-1].copt).penalty;
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
#line 8831 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 652:
#line 2267 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].fval) == 0.0) {
            fprintf(stderr," *** WARNING: zero coefficient in LMPC\n");
            fprintf(stderr," ***          node %d dof %d\n",(yyvsp[-3].ival),(yyvsp[-2].ival));
          }
          (yyval.mpcterm) = new LMPCTerm();
          (yyval.mpcterm)->nnum = (yyvsp[-3].ival)-1;
          (yyval.mpcterm)->dofnum = (yyvsp[-2].ival)-1;
          (yyval.mpcterm)->coef.r_value = (yyvsp[-1].fval);
          if((yyvsp[-2].ival) == 1){
            (yyval.mpcterm)->coef.r_value /= domain->solInfo().xLMPCFactor;
          } else if ((yyvsp[-2].ival) == 2) {
            (yyval.mpcterm)->coef.r_value /= domain->solInfo().yLMPCFactor;
          } else if ((yyvsp[-2].ival) == 3) {
            (yyval.mpcterm)->coef.r_value /= domain->solInfo().zLMPCFactor;
          }
        }
#line 8852 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 655:
#line 2290 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-1].cxbcval).nnum,(yyvsp[-1].cxbcval).reval,(yyvsp[-1].cxbcval).imval,(yyvsp[0].mpcterm)); domain->addLMPC((yyval.lmpcons)); }
#line 8858 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 656:
#line 2292 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons)->addterm((yyvsp[0].mpcterm)); }
#line 8864 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 657:
#line 2294 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-1].cxbcval).nnum,(yyvsp[-1].cxbcval).reval,(yyvsp[-1].cxbcval).imval,(yyvsp[0].mpcterm)); domain->addLMPC((yyval.lmpcons)); }
#line 8870 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 658:
#line 2298 "p.y" /* yacc.c:1646  */
    { (yyval.cxbcval).nnum=(yyvsp[-4].ival); (yyval.cxbcval).reval=(yyvsp[-2].fval); (yyval.cxbcval).imval=(yyvsp[-1].fval); }
#line 8876 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 659:
#line 2300 "p.y" /* yacc.c:1646  */
    { (yyval.cxbcval).nnum=(yyvsp[-2].ival); (yyval.cxbcval).reval=(yyvsp[-1].fval); (yyval.cxbcval).imval=0.0; }
#line 8882 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 660:
#line 2302 "p.y" /* yacc.c:1646  */
    { (yyval.cxbcval).nnum=(yyvsp[-1].ival); (yyval.cxbcval).reval=0.0; (yyval.cxbcval).imval=0.0; }
#line 8888 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 661:
#line 2306 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-2].fval)==0.0) && ((yyvsp[-1].fval)==0.0)) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",(yyvsp[-4].ival),(yyvsp[-3].ival));
          return -1;
          }
          else { (yyval.mpcterm) = new LMPCTerm(true); (yyval.mpcterm)->nnum=((yyvsp[-4].ival)-1); (yyval.mpcterm)->dofnum=((yyvsp[-3].ival)-1); (yyval.mpcterm)->coef.c_value=DComplex((yyvsp[-2].fval), (yyvsp[-1].fval)); }
        }
#line 8900 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 662:
#line 2314 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].fval)==0.0) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",(yyvsp[-3].ival),(yyvsp[-2].ival));
          return -1;
          }
          else { (yyval.mpcterm) = new LMPCTerm(true); (yyval.mpcterm)->nnum=((yyvsp[-3].ival)-1); (yyval.mpcterm)->dofnum=((yyvsp[-2].ival)-1); (yyval.mpcterm)->coef.c_value=DComplex((yyvsp[-1].fval),0.0); }
        }
#line 8912 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 663:
#line 2324 "p.y" /* yacc.c:1646  */
    { (yyval.cxbclist) = (yyvsp[0].cxbclist); }
#line 8918 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 664:
#line 2326 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].cxbclist)->n; ++i) (yyvsp[0].cxbclist)->d[i].loadsetid = (yyvsp[-2].ival);
          (yyval.cxbclist) = (yyvsp[0].cxbclist); }
#line 8925 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 665:
#line 2331 "p.y" /* yacc.c:1646  */
    { (yyval.cxbclist) = new ComplexBCList; (yyval.cxbclist)->add((yyvsp[0].cxbcval)); }
#line 8931 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 666:
#line 2333 "p.y" /* yacc.c:1646  */
    { (yyval.cxbclist) = (yyvsp[-1].cxbclist); (yyval.cxbclist)->add((yyvsp[0].cxbcval)); }
#line 8937 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 669:
#line 2341 "p.y" /* yacc.c:1646  */
    { StructProp sp; 
	  sp.A = (yyvsp[-14].fval);  sp.E = (yyvsp[-13].fval);  sp.nu  = (yyvsp[-12].fval);  sp.rho = (yyvsp[-11].fval);
          sp.c = (yyvsp[-10].fval);  sp.k = (yyvsp[-9].fval);  sp.eh  = (yyvsp[-8].fval);  sp.P   = (yyvsp[-7].fval);  sp.Ta  = (yyvsp[-6].fval); 
          sp.Q = (yyvsp[-5].fval); sp.W = (yyvsp[-4].fval); sp.Ixx = (yyvsp[-3].fval); sp.Iyy = (yyvsp[-2].fval); sp.Izz = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-15].ival)-1, sp );
        }
#line 8948 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 670:
#line 2348 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.A = (yyvsp[-17].fval);  sp.E = (yyvsp[-16].fval);  sp.nu  = (yyvsp[-15].fval);  sp.rho = (yyvsp[-14].fval);
          sp.c = (yyvsp[-13].fval);  sp.k = (yyvsp[-12].fval);  sp.eh  = (yyvsp[-11].fval);  sp.P   = (yyvsp[-10].fval);  sp.Ta  = (yyvsp[-9].fval);
          sp.Q = (yyvsp[-8].fval); sp.W = (yyvsp[-7].fval); sp.Ixx = (yyvsp[-6].fval); sp.Iyy = (yyvsp[-5].fval); sp.Izz = (yyvsp[-4].fval);
          switch((yyvsp[-3].ival)) {
            case 1 : // RAYDAMP
              sp.betaDamp = (yyvsp[-2].fval); sp.alphaDamp = (yyvsp[-1].fval);
              break;
            case 2 : // STRDAMP
              sp.etaDamp = (yyvsp[-2].fval); sp.betaDamp = (yyvsp[-1].fval);
              break;
            case 3 : // RUBDAMP
              sp.eta_E = (yyvsp[-2].fval);
              if(sp.eta_E >= 0) sp.eta_mu = (yyvsp[-1].fval);
              sp.E0 = (yyvsp[-16].fval);
              sp.mu0 = (yyvsp[-16].fval)/(2*(1+(yyvsp[-15].fval)));
              break;
            default :
              return -1;
          }
          geoSource->addMat( (yyvsp[-18].ival)-1, sp );
        }
#line 8975 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 671:
#line 2371 "p.y" /* yacc.c:1646  */
    { StructProp sp; 
	  sp.A = (yyvsp[-18].fval);  sp.E = (yyvsp[-17].fval);  sp.nu  = (yyvsp[-16].fval);  sp.rho = (yyvsp[-15].fval);
          sp.c = (yyvsp[-14].fval);  sp.k = (yyvsp[-13].fval);  sp.eh  = (yyvsp[-12].fval);  sp.P   = (yyvsp[-11].fval);  sp.Ta  = (yyvsp[-10].fval); 
          sp.Q = (yyvsp[-9].fval); sp.W = (yyvsp[-8].fval); sp.Ixx = (yyvsp[-7].fval); sp.Iyy = (yyvsp[-6].fval); sp.Izz = (yyvsp[-5].fval);
	  sp.ymin = (yyvsp[-4].fval); sp.ymax = (yyvsp[-3].fval); sp.zmin = (yyvsp[-2].fval); sp.zmax = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-19].ival)-1, sp );
        }
#line 8987 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 672:
#line 2379 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.A = (yyvsp[-21].fval);  sp.E = (yyvsp[-20].fval);  sp.nu  = (yyvsp[-19].fval);  sp.rho = (yyvsp[-18].fval);
          sp.c = (yyvsp[-17].fval);  sp.k = (yyvsp[-16].fval);  sp.eh  = (yyvsp[-15].fval);  sp.P   = (yyvsp[-14].fval);  sp.Ta  = (yyvsp[-13].fval);
          sp.Q = (yyvsp[-12].fval); sp.W = (yyvsp[-11].fval); sp.Ixx = (yyvsp[-10].fval); sp.Iyy = (yyvsp[-9].fval); sp.Izz = (yyvsp[-8].fval);
          sp.ymin = (yyvsp[-7].fval); sp.ymax = (yyvsp[-6].fval); sp.zmin = (yyvsp[-5].fval); sp.zmax = (yyvsp[-4].fval);
          switch((yyvsp[-3].ival)) {
            case 1 : // RAYDAMP
              sp.betaDamp = (yyvsp[-2].fval); sp.alphaDamp = (yyvsp[-1].fval);
              break;
            case 2 : // STRDAMP
              sp.etaDamp = (yyvsp[-2].fval); sp.betaDamp = (yyvsp[-1].fval);
              break;
            case 3 : // RUBDAMP
              sp.eta_E = (yyvsp[-2].fval); 
              if(sp.eta_E >= 0) sp.eta_mu = (yyvsp[-1].fval);
              sp.E0 = (yyvsp[-20].fval);
              sp.mu0 = (yyvsp[-20].fval)/(2*(1+(yyvsp[-19].fval)));
              break;
            default :
              return -1;
          }
          geoSource->addMat( (yyvsp[-22].ival)-1, sp );
        }
#line 9015 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 673:
#line 2403 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.A = (yyvsp[-7].fval); sp.E = (yyvsp[-6].fval); sp.nu = (yyvsp[-5].fval); sp.rho = (yyvsp[-4].fval);
          sp.c = (yyvsp[-3].fval); sp.k = (yyvsp[-2].fval); sp.eh = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-8].ival)-1, sp ); 
        }
#line 9025 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 674:
#line 2409 "p.y" /* yacc.c:1646  */
    { StructProp sp;  // this is for spring: GID Kx Ky Kz lx1 ...
          sp.A = (yyvsp[-12].fval);  sp.E = (yyvsp[-11].fval);  sp.nu  = (yyvsp[-10].fval);  sp.rho = (yyvsp[-9].fval);
          sp.c = (yyvsp[-8].fval);  sp.k = (yyvsp[-7].fval);  sp.eh  = (yyvsp[-6].fval);  sp.P   = (yyvsp[-5].fval);  sp.Ta  = (yyvsp[-4].fval);
          sp.Q = (yyvsp[-3].fval); sp.W = (yyvsp[-2].fval); sp.Ixx = (yyvsp[-1].fval);  
          geoSource->addMat( (yyvsp[-13].ival)-1, sp );
        }
#line 9036 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 675:
#line 2416 "p.y" /* yacc.c:1646  */
    { StructProp sp;  // this is for spring with stiffness-proportional damping : GID Kx Ky Kz lx1 ...
          sp.A = (yyvsp[-14].fval);  sp.E = (yyvsp[-13].fval);  sp.nu  = (yyvsp[-12].fval);  sp.rho = (yyvsp[-11].fval);
          sp.c = (yyvsp[-10].fval);  sp.k = (yyvsp[-9].fval);  sp.eh  = (yyvsp[-8].fval);  sp.P   = (yyvsp[-7].fval);  sp.Ta  = (yyvsp[-6].fval);
          sp.Q = (yyvsp[-5].fval); sp.W = (yyvsp[-4].fval); sp.Ixx = (yyvsp[-3].fval);
          if((yyvsp[-2].ival) == 0) sp.betaDamp = (yyvsp[-1].fval); else return -1;
          geoSource->addMat( (yyvsp[-15].ival)-1, sp );
        }
#line 9048 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 676:
#line 2425 "p.y" /* yacc.c:1646  */
    { StructProp sp; // this is used for the reduced mesh file output in Rom.d/MeshOutput.C
                         // all properties relevant to structural nonlinear dynamics should be included
          sp.A = (yyvsp[-28].fval);  sp.E = (yyvsp[-27].fval);  sp.nu  = (yyvsp[-26].fval);  sp.rho = (yyvsp[-25].fval);
          sp.c = (yyvsp[-24].fval);  sp.k = (yyvsp[-23].fval);  sp.eh  = (yyvsp[-22].fval);  sp.P   = (yyvsp[-21].fval);  sp.Ta  = (yyvsp[-20].fval);
          sp.Q = (yyvsp[-19].fval); sp.W = (yyvsp[-18].fval); sp.Ixx = (yyvsp[-17].fval); sp.Iyy = (yyvsp[-16].fval); sp.Izz = (yyvsp[-15].fval);
          sp.ymin = (yyvsp[-14].fval); sp.ymax = (yyvsp[-13].fval); sp.zmin = (yyvsp[-12].fval); sp.zmax = (yyvsp[-11].fval);
          sp.betaDamp = (yyvsp[-10].fval); sp.alphaDamp = (yyvsp[-9].fval); 
          sp.lagrangeMult = bool((yyvsp[-8].ival)); sp.penalty = (yyvsp[-7].fval); sp.initialPenalty = (yyvsp[-6].fval);
          sp.funtype = (yyvsp[-5].ival); sp.type = StructProp::PropType((yyvsp[-4].ival)); sp.k1 = (yyvsp[-3].fval); sp.k2 = (yyvsp[-2].fval); sp.k3 = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-29].ival)-1, sp );
        }
#line 9064 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 677:
#line 2437 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[-9].fval),0.0);
          sp.fp.PMLtype = int((yyvsp[-8].fval));
          sp.fp.gamma = (yyvsp[-7].fval);
          sp.fp.Rx = (yyvsp[-6].fval);
          sp.fp.Sx = (yyvsp[-5].fval);
          sp.fp.Ry = (yyvsp[-4].fval);
          sp.fp.Sy = (yyvsp[-3].fval);
          sp.fp.Rz = (yyvsp[-2].fval);
          sp.fp.Sz = (yyvsp[-1].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[-11].ival)-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        }
#line 9084 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 678:
#line 2453 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[-10].fval),0.0);
          sp.rho = (yyvsp[-9].fval);
          sp.fp.PMLtype = int((yyvsp[-8].fval));
          sp.fp.gamma = (yyvsp[-7].fval);
          sp.fp.Rx = (yyvsp[-6].fval);
          sp.fp.Sx = (yyvsp[-5].fval);
          sp.fp.Ry = (yyvsp[-4].fval);
          sp.fp.Sy = (yyvsp[-3].fval);
          sp.fp.Rz = (yyvsp[-2].fval);
          sp.fp.Sz = (yyvsp[-1].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[-12].ival)-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        }
#line 9105 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 679:
#line 2470 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[-11].fval),(yyvsp[-10].fval));
          sp.rho = (yyvsp[-9].fval);
          sp.fp.PMLtype = int((yyvsp[-8].fval));
          sp.fp.gamma = (yyvsp[-7].fval);
          sp.fp.Rx = (yyvsp[-6].fval);
          sp.fp.Sx = (yyvsp[-5].fval);
          sp.fp.Ry = (yyvsp[-4].fval);
          sp.fp.Sy = (yyvsp[-3].fval);
          sp.fp.Rz = (yyvsp[-2].fval);
          sp.fp.Sz = (yyvsp[-1].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[-13].ival)-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        }
#line 9126 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 680:
#line 2487 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[-2].fval),0.0);
          sp.rho = (yyvsp[-1].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[-4].ival)-1, sp );
          domain->solInfo().acoustic = true;
        }
#line 9138 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 681:
#line 2495 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[-3].fval),(yyvsp[-2].fval));
          sp.rho = (yyvsp[-1].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[-5].ival)-1, sp );
          domain->solInfo().acoustic = true;
        }
#line 9150 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 682:
#line 2503 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.F_op = (yyvsp[-13].ival);
          sp.E = (yyvsp[-12].fval);
          sp.rho = (yyvsp[-11].fval);
	  sp.A = (yyvsp[-10].fval);
          sp.F_Uc = (yyvsp[-9].fval);
          sp.F_Uf = (yyvsp[-8].fval);
          sp.lambda = (yyvsp[-7].fval);
          sp.F_h = (yyvsp[-6].fval);
          sp.F_d = (yyvsp[-5].fval);
          sp.F_dlambda = (yyvsp[-4].fval);
          sp.F_np = (yyvsp[-3].ival);
          sp.F_Nf = (yyvsp[-2].ival);
	  sp.Seed = (yyvsp[-1].ival);
          sp.type = StructProp::Fabric;
          geoSource->addMat( (yyvsp[-15].ival)-1, sp );
        }
#line 9172 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 683:
#line 2521 "p.y" /* yacc.c:1646  */
    { StructProp sp; 
          sp.A = (yyvsp[-9].fval); sp.rho = (yyvsp[-8].fval); sp.Q = (yyvsp[-7].fval); sp.c = (yyvsp[-6].fval); 
          sp.sigma = (yyvsp[-5].fval); sp.k = (yyvsp[-4].fval); sp.eh = (yyvsp[-3].fval); sp.P = (yyvsp[-2].fval); sp.Ta = (yyvsp[-1].fval);
          sp.type = StructProp::Thermal;
          geoSource->addMat( (yyvsp[-11].ival)-1, sp );
        }
#line 9183 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 684:
#line 2528 "p.y" /* yacc.c:1646  */
    { // rigid element or joint with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-2].ival)-1, sp );
        }
#line 9194 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 685:
#line 2535 "p.y" /* yacc.c:1646  */
    { // rigid element or joint
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-1].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-1].copt).penalty;
          sp.constraint_hess = (yyvsp[-1].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-1].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-3].ival)-1, sp );
        }
#line 9209 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 686:
#line 2546 "p.y" /* yacc.c:1646  */
    { //           m     Ixx   Iyy   Izz   Ixy   Iyz   Ixz   cx    cy    cz
          // discrete mass with offset
          StructProp sp;
          sp.rho = (yyvsp[-10].fval);
          sp.Ixx = (yyvsp[-9].fval);
          sp.Iyy = (yyvsp[-8].fval);
          sp.Izz = (yyvsp[-7].fval);
          sp.Ixy = (yyvsp[-6].fval);
          sp.Iyz = (yyvsp[-5].fval);
          sp.Ixz = (yyvsp[-4].fval);
          sp.cx  = (yyvsp[-3].fval);
          sp.cy  = (yyvsp[-2].fval);
          sp.cz  = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-12].ival)-1, sp );
        }
#line 9229 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 687:
#line 2562 "p.y" /* yacc.c:1646  */
    { // rigid 8-node brick element with mass, and default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-4].ival)-1, sp );
        }
#line 9240 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 688:
#line 2569 "p.y" /* yacc.c:1646  */
    { // rigid 8-node brick element with mass
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-3].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-3].copt).penalty;
          sp.constraint_hess = (yyvsp[-3].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-3].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.rho = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-5].ival)-1, sp );
        }
#line 9255 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 689:
#line 2580 "p.y" /* yacc.c:1646  */
    { // rigid beam or shell element with mass, and default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = (yyvsp[-2].fval);
          sp.A = sp.eh = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-5].ival)-1, sp );
        }
#line 9267 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 690:
#line 2588 "p.y" /* yacc.c:1646  */
    { // rigid beam or shell element with mass
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-4].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-4].copt).penalty;
          sp.constraint_hess = (yyvsp[-4].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-4].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.rho = (yyvsp[-2].fval);
          sp.A = sp.eh = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-6].ival)-1, sp );
        }
#line 9283 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 691:
#line 2600 "p.y" /* yacc.c:1646  */
    { // constraint function element XXX deprecated format
          StructProp sp;
          sp.lagrangeMult = bool((yyvsp[-8].ival));
          sp.initialPenalty = sp.penalty = (yyvsp[-7].fval);
          sp.amplitude = (yyvsp[-6].fval);
          sp.omega = (yyvsp[-5].fval);
          sp.phase = (yyvsp[-4].fval);
          sp.B = (yyvsp[-3].fval);
          sp.C = (yyvsp[-2].fval);
          sp.relop = (yyvsp[-1].ival);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-10].ival)-1, sp );
        }
#line 9302 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 692:
#line 2615 "p.y" /* yacc.c:1646  */
    { // constraint function element with default constraint options
          StructProp sp;
          sp.amplitude = 0;
          sp.omega = 0;
          sp.phase = 0;
          sp.B = 0;
          sp.C = 0;
          sp.relop = (yyvsp[-1].ival);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-4].ival)-1, sp );
        }
#line 9319 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 693:
#line 2628 "p.y" /* yacc.c:1646  */
    { // constraint function element
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-3].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-3].copt).penalty;
          sp.constraint_hess = (yyvsp[-3].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-3].copt).constraint_hess_eps;
          sp.amplitude = 0;
          sp.omega = 0;
          sp.phase = 0;
          sp.B = 0;
          sp.C = 0;
          sp.relop = (yyvsp[-1].ival);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-5].ival)-1, sp );
        }
#line 9340 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 694:
#line 2645 "p.y" /* yacc.c:1646  */
    { // joint-with-driver, 2-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[-5].ival);
          sp.amplitude = (yyvsp[-4].fval);
          sp.offset = (yyvsp[-3].fval);
          sp.c1 = (yyvsp[-2].fval);
          sp.c2 = (yyvsp[-1].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-7].ival)-1, sp );
        }
#line 9356 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 695:
#line 2657 "p.y" /* yacc.c:1646  */
    { // joint-with-driver, 3-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[-6].ival);
          sp.amplitude = (yyvsp[-5].fval);
          sp.offset = (yyvsp[-4].fval);
          sp.c1 = (yyvsp[-3].fval);
          sp.c2 = (yyvsp[-2].fval);
          sp.c3 = (yyvsp[-1].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-8].ival)-1, sp );
        }
#line 9373 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 696:
#line 2670 "p.y" /* yacc.c:1646  */
    { // joint-with-driver, 4-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[-7].ival);
          sp.amplitude = (yyvsp[-6].fval);
          sp.offset = (yyvsp[-5].fval);
          sp.c1 = (yyvsp[-4].fval);
          sp.c2 = (yyvsp[-3].fval);
          sp.c3 = (yyvsp[-2].fval);
          sp.c4 = (yyvsp[-1].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-9].ival)-1, sp );
        }
#line 9391 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 697:
#line 2684 "p.y" /* yacc.c:1646  */
    { // joint-with-driver, 2-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-6].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-6].copt).penalty;
          sp.constraint_hess = (yyvsp[-6].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-6].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[-5].ival);
          sp.amplitude = (yyvsp[-4].fval);
          sp.offset = (yyvsp[-3].fval);
          sp.c1 = (yyvsp[-2].fval);
          sp.c2 = (yyvsp[-1].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-8].ival)-1, sp );
        }
#line 9411 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 698:
#line 2700 "p.y" /* yacc.c:1646  */
    { // joint-with-driver, 3-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-7].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-7].copt).penalty;
          sp.constraint_hess = (yyvsp[-7].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-7].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[-6].ival);
          sp.amplitude = (yyvsp[-5].fval);
          sp.offset = (yyvsp[-4].fval);
          sp.c1 = (yyvsp[-3].fval);
          sp.c2 = (yyvsp[-2].fval);
          sp.c3 = (yyvsp[-1].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-9].ival)-1, sp );
        }
#line 9432 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 699:
#line 2717 "p.y" /* yacc.c:1646  */
    { // joint-with-driver, 4-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-8].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-8].copt).penalty;
          sp.constraint_hess = (yyvsp[-8].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-8].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[-7].ival);
          sp.amplitude = (yyvsp[-6].fval);
          sp.offset = (yyvsp[-5].fval);
          sp.c1 = (yyvsp[-4].fval);
          sp.c2 = (yyvsp[-3].fval);
          sp.c3 = (yyvsp[-2].fval);
          sp.c4 = (yyvsp[-1].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-10].ival)-1, sp );
        }
#line 9454 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 700:
#line 2735 "p.y" /* yacc.c:1646  */
    { // actuated joint, 2-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[-7].ival);
          sp.amplitude = (yyvsp[-6].fval);
          sp.offset = (yyvsp[-5].fval);
          sp.c1 = (yyvsp[-4].fval);
          sp.c2 = (yyvsp[-3].fval);
          sp.k1 = (yyvsp[-1].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-9].ival)-1, sp );
        }
#line 9471 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 701:
#line 2748 "p.y" /* yacc.c:1646  */
    { // actuated joint, 3-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[-8].ival);
          sp.amplitude = (yyvsp[-7].fval);
          sp.offset = (yyvsp[-6].fval);
          sp.c1 = (yyvsp[-5].fval);
          sp.c2 = (yyvsp[-4].fval);
          sp.c3 = (yyvsp[-3].fval);
          sp.k1 = (yyvsp[-1].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-10].ival)-1, sp );
        }
#line 9489 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 702:
#line 2762 "p.y" /* yacc.c:1646  */
    { // actuated joints, 4-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[-9].ival);
          sp.amplitude = (yyvsp[-8].fval);
          sp.offset = (yyvsp[-7].fval);
          sp.c1 = (yyvsp[-6].fval);
          sp.c2 = (yyvsp[-5].fval);
          sp.c3 = (yyvsp[-4].fval);
          sp.c4 = (yyvsp[-3].fval);
          sp.k1 = (yyvsp[-1].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-11].ival)-1, sp );
        }
#line 9508 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 703:
#line 2777 "p.y" /* yacc.c:1646  */
    { // actuated joint, 2-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-8].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-8].copt).penalty;
          sp.constraint_hess = (yyvsp[-8].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-8].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[-7].ival);
          sp.amplitude = (yyvsp[-6].fval);
          sp.offset = (yyvsp[-5].fval);
          sp.c1 = (yyvsp[-4].fval);
          sp.c2 = (yyvsp[-3].fval);
          sp.k1 = (yyvsp[-1].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-10].ival)-1, sp );
        }
#line 9529 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 704:
#line 2794 "p.y" /* yacc.c:1646  */
    { // actuated joint, 3-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-9].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-9].copt).penalty;
          sp.constraint_hess = (yyvsp[-9].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-9].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[-8].ival);
          sp.amplitude = (yyvsp[-7].fval);
          sp.offset = (yyvsp[-6].fval);
          sp.c1 = (yyvsp[-5].fval);
          sp.c2 = (yyvsp[-4].fval);
          sp.c3 = (yyvsp[-3].fval);
          sp.k1 = (yyvsp[-1].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-11].ival)-1, sp );
        }
#line 9551 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 705:
#line 2812 "p.y" /* yacc.c:1646  */
    { // actuated joints, 4-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-10].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-10].copt).penalty;
          sp.constraint_hess = (yyvsp[-10].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-10].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[-9].ival);
          sp.amplitude = (yyvsp[-8].fval);
          sp.offset = (yyvsp[-7].fval);
          sp.c1 = (yyvsp[-6].fval);
          sp.c2 = (yyvsp[-5].fval);
          sp.c3 = (yyvsp[-4].fval);
          sp.c4 = (yyvsp[-3].fval);
          sp.k1 = (yyvsp[-1].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-12].ival)-1, sp );
        }
#line 9574 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 706:
#line 2831 "p.y" /* yacc.c:1646  */
    { // RevoluteJointSpringCombo with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[-1].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-4].ival)-1, sp );
        }
#line 9586 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 707:
#line 2839 "p.y" /* yacc.c:1646  */
    { // RevoluteJointSpringComboWithFreeplay with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[-3].fval);
          sp.freeplay[0] = (yyvsp[-1].freeplayProps);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-6].ival)-1, sp );
        }
#line 9599 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 708:
#line 2848 "p.y" /* yacc.c:1646  */
    { // RevoluteJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-3].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-3].copt).penalty;
          sp.constraint_hess = (yyvsp[-3].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-3].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[-1].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-5].ival)-1, sp );
        }
#line 9615 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 709:
#line 2860 "p.y" /* yacc.c:1646  */
    { // RevoluteJointSpringComboWithFreeplay
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-5].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-5].copt).penalty;
          sp.constraint_hess = (yyvsp[-5].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-5].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[-3].fval);
          sp.freeplay[0] = (yyvsp[-1].freeplayProps);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-7].ival)-1, sp );
        }
#line 9632 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 710:
#line 2873 "p.y" /* yacc.c:1646  */
    { // UniversalJointSpringCombo with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[-2].fval);
          sp.k2 = (yyvsp[-1].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-5].ival)-1, sp );
        }
#line 9645 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 711:
#line 2882 "p.y" /* yacc.c:1646  */
    { // UniversalJointSpringComboWithFreeplay with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[-5].fval);
          sp.k2 = (yyvsp[-4].fval);
          sp.freeplay[0] = (yyvsp[-2].freeplayProps);
          sp.freeplay[1] = (yyvsp[-1].freeplayProps);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-8].ival)-1, sp );
        }
#line 9660 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 712:
#line 2893 "p.y" /* yacc.c:1646  */
    { // UniversalJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-4].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-4].copt).penalty;
          sp.constraint_hess = (yyvsp[-4].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-4].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[-2].fval);
          sp.k2 = (yyvsp[-1].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-6].ival)-1, sp );
        }
#line 9677 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 713:
#line 2906 "p.y" /* yacc.c:1646  */
    { // UniversalJointSpringComboWithFreeplay
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-7].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-7].copt).penalty;
          sp.constraint_hess = (yyvsp[-7].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-7].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[-5].fval);
          sp.k2 = (yyvsp[-4].fval);
          sp.freeplay[0] = (yyvsp[-2].freeplayProps);
          sp.freeplay[1] = (yyvsp[-1].freeplayProps);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-9].ival)-1, sp );
        }
#line 9696 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 714:
#line 2921 "p.y" /* yacc.c:1646  */
    { // SphericalJointSpringCombo with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[-3].fval);
          sp.k2 = (yyvsp[-2].fval);
          sp.k3 = (yyvsp[-1].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-6].ival)-1, sp );
        }
#line 9710 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 715:
#line 2931 "p.y" /* yacc.c:1646  */
    { // SphericalJointSpringComboWithFreeplay with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[-7].fval);
          sp.k2 = (yyvsp[-6].fval);
          sp.k3 = (yyvsp[-5].fval);
          sp.freeplay[0] = (yyvsp[-3].freeplayProps);
          sp.freeplay[1] = (yyvsp[-2].freeplayProps);
          sp.freeplay[2] = (yyvsp[-1].freeplayProps);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-10].ival)-1, sp );
        }
#line 9727 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 716:
#line 2944 "p.y" /* yacc.c:1646  */
    { // SphericalJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-5].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-5].copt).penalty;
          sp.constraint_hess = (yyvsp[-5].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-5].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[-3].fval);
          sp.k2 = (yyvsp[-2].fval);
          sp.k3 = (yyvsp[-1].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-7].ival)-1, sp );
        }
#line 9745 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 717:
#line 2958 "p.y" /* yacc.c:1646  */
    { // SphericalJointSpringComboWithFreeplay
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-9].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-9].copt).penalty;
          sp.constraint_hess = (yyvsp[-9].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-9].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[-7].fval);
          sp.k2 = (yyvsp[-6].fval);
          sp.k3 = (yyvsp[-5].fval);
          sp.freeplay[0] = (yyvsp[-3].freeplayProps);
          sp.freeplay[1] = (yyvsp[-2].freeplayProps);
          sp.freeplay[2] = (yyvsp[-1].freeplayProps);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-11].ival)-1, sp );
        }
#line 9766 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 718:
#line 2975 "p.y" /* yacc.c:1646  */
    { // TorsionalSpring or TranslationalSpring (types 200,201,202)
          StructProp sp;
          sp.k1 = (yyvsp[-1].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-3].ival)-1, sp );
        }
#line 9777 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 719:
#line 2982 "p.y" /* yacc.c:1646  */
    { // 1-sided TorsionalSpring or TranslationalSpring with freeplay (types 203,204,205)
          StructProp sp;
          sp.k1 = (yyvsp[-3].fval);
          sp.rho = 0;
          sp.freeplay[0].ul = (yyvsp[-1].fval);
          sp.freeplay[0].dz = 0.0;
          sp.freeplay[0].uz = 1.0;
          geoSource->addMat( (yyvsp[-5].ival)-1, sp );
        }
#line 9791 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 720:
#line 2992 "p.y" /* yacc.c:1646  */
    { // Acoustic rubber
          StructProp sp;
          sp.E0 = sp.E = (yyvsp[-9].fval); sp.dE = (yyvsp[-8].fval)/(2*M_PI); sp.eta_E = (yyvsp[-7].fval); sp.deta_E = (yyvsp[-6].fval)/(2*M_PI);
          sp.mu0 = (yyvsp[-5].fval); sp.dmu = (yyvsp[-4].fval)/(2*M_PI); sp.eta_mu = (yyvsp[-3].fval); sp.deta_mu = (yyvsp[-2].fval)/(2*M_PI);
          sp.rho = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-11].ival)-1, sp );
        }
#line 9803 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 721:
#line 3000 "p.y" /* yacc.c:1646  */
    { // Acoustic rubber
          StructProp sp;
          sp.E0 = sp.E = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-3].ival)-1, sp );
        }
#line 9813 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 722:
#line 3016 "p.y" /* yacc.c:1646  */
    { // 5-parameter freeplay model
          (yyval.freeplayProps).ll = (yyvsp[-4].fval);
          (yyval.freeplayProps).ul = (yyvsp[-3].fval);
          (yyval.freeplayProps).lz = (yyvsp[-2].fval);
          (yyval.freeplayProps).dz = (yyvsp[-1].fval);
          (yyval.freeplayProps).uz = (yyvsp[0].fval);
        }
#line 9825 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 725:
#line 3029 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].ival) == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[-1].ival));
          (yyval.SurfObj)->SetReverseNormals(false);
          domain->AddSurfaceEntity((yyval.SurfObj));
        }
#line 9835 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 726:
#line 3035 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-2].ival) == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[-2].ival));
          (yyval.SurfObj)->SetReverseNormals(true);
          domain->AddSurfaceEntity((yyval.SurfObj));
        }
#line 9845 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 727:
#line 3041 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-3].ival) == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[-3].ival));
          (yyval.SurfObj)->SetIsShellFace(true);
          (yyval.SurfObj)->SetShellThickness((yyvsp[-1].fval));
          domain->AddSurfaceEntity((yyval.SurfObj));
        }
#line 9856 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 728:
#line 3048 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-4].ival) == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[-4].ival));
          (yyval.SurfObj)->SetIsShellFace(true);
          (yyval.SurfObj)->SetShellThickness((yyvsp[-2].fval));
          (yyval.SurfObj)->SetReverseNormals(true);
          domain->AddSurfaceEntity((yyval.SurfObj));
        }
#line 9868 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 729:
#line 3056 "p.y" /* yacc.c:1646  */
    { if((yyval.SurfObj)->GetReverseNormals()) { // reverse the node numbering
            int *nodes = new int[(yyvsp[-1].nl).num];
            for(int i=0; i<(yyvsp[-1].nl).num; ++i) nodes[(yyvsp[-1].nl).num-1-i] = (yyvsp[-1].nl).nd[i];
            (yyval.SurfObj)->AddFaceElement((yyvsp[-3].ival)-1, (yyvsp[-2].ival), (yyvsp[-1].nl).num, nodes);
            delete [] nodes;
          }
          else (yyval.SurfObj)->AddFaceElement((yyvsp[-3].ival)-1, (yyvsp[-2].ival), (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd);
        }
#line 9881 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 730:
#line 3067 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-2].ival), (yyvsp[-1].ival)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9887 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 731:
#line 3069 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-3].ival), (yyvsp[-2].ival)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9893 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 732:
#line 3071 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-3].ival), (yyvsp[-2].ival)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
#line 9901 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 733:
#line 3075 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-1].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9907 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 734:
#line 3077 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9913 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 735:
#line 3079 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-1].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9919 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 736:
#line 3081 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9925 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 737:
#line 3083 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-1].fval)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
#line 9933 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 738:
#line 3087 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
#line 9941 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 739:
#line 3093 "p.y" /* yacc.c:1646  */
    { domain->addWetInterface((yyvsp[-2].ival), (yyvsp[-1].ival)); domain->solInfo().isCoupled = true; }
#line 9947 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 740:
#line 3095 "p.y" /* yacc.c:1646  */
    { domain->addWetInterface((yyvsp[-1].ival), (yyvsp[-1].ival)); 
          domain->solInfo().isCoupled  = true; 
          domain->solInfo().isMatching = true; }
#line 9955 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 741:
#line 3103 "p.y" /* yacc.c:1646  */
    { }
#line 9961 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 742:
#line 3105 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-2].ival), (yyvsp[-1].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
#line 9971 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 743:
#line 3111 "p.y" /* yacc.c:1646  */
    { 
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-3].ival), (yyvsp[-2].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-1].ival));
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
#line 9982 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 744:
#line 3118 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-1].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-2].ival));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 9993 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 745:
#line 3125 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-2].fval), (yyvsp[-1].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-3].ival));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 10004 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 746:
#line 3132 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-7].ival), (yyvsp[-6].ival), (yyvsp[-4].fval), (yyvsp[-3].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-5].ival));
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[-2].ival), (yyvsp[-1].fval));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 10016 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 747:
#line 3140 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-3].ival), (yyvsp[-2].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[-1].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 10027 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 748:
#line 3147 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-4].ival), (yyvsp[-3].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-2].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[-1].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 10039 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 749:
#line 3155 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-2].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-3].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[-1].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 10051 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 750:
#line 3163 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-3].fval), (yyvsp[-2].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-4].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[-1].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 10063 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 751:
#line 3171 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-8].ival), (yyvsp[-7].ival), (yyvsp[-5].fval), (yyvsp[-4].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-6].ival));
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[-3].ival), (yyvsp[-2].fval));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[-1].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 10076 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 752:
#line 3182 "p.y" /* yacc.c:1646  */
    { }
#line 10082 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 753:
#line 3184 "p.y" /* yacc.c:1646  */
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[-2].ival), (yyvsp[-1].ival)); domain->solInfo().isCoupled = true; 
          if((yyvsp[-2].ival) == (yyvsp[-1].ival)) domain->solInfo().isMatching = true;
        }
#line 10092 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 754:
#line 3190 "p.y" /* yacc.c:1646  */
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); domain->solInfo().isCoupled = true;
          if((yyvsp[-4].ival) == (yyvsp[-3].ival)) domain->solInfo().isMatching = true;
        }
#line 10102 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 755:
#line 3198 "p.y" /* yacc.c:1646  */
    { domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; }
#line 10109 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 757:
#line 3204 "p.y" /* yacc.c:1646  */
    { domain->addWetElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd);
          domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; }
#line 10117 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 758:
#line 3210 "p.y" /* yacc.c:1646  */
    { }
#line 10123 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 759:
#line 3212 "p.y" /* yacc.c:1646  */
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[-2].ival), (yyvsp[-1].ival)); domain->solInfo().HEV = 1;
          if((yyvsp[-2].ival) == (yyvsp[-1].ival)) domain->solInfo().isMatching = true;
        }
#line 10133 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 760:
#line 3218 "p.y" /* yacc.c:1646  */
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); domain->solInfo().HEV = 1;
          if((yyvsp[-4].ival) == (yyvsp[-3].ival)) domain->solInfo().isMatching = true;
        }
#line 10143 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 761:
#line 3228 "p.y" /* yacc.c:1646  */
    { }
#line 10149 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 762:
#line 3230 "p.y" /* yacc.c:1646  */
    { domain->solInfo().contactsurface_mode = (yyvsp[-1].ival); }
#line 10155 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 763:
#line 3232 "p.y" /* yacc.c:1646  */
    { domain->AddMortarCond((yyvsp[-1].MortarCondObj)); }
#line 10161 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 764:
#line 3234 "p.y" /* yacc.c:1646  */
    { (yyvsp[-2].MortarCondObj)->SetConstraintOptions((yyvsp[-1].copt)); domain->AddMortarCond((yyvsp[-2].MortarCondObj)); }
#line 10167 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 765:
#line 3236 "p.y" /* yacc.c:1646  */
    { (yyvsp[-4].MortarCondObj)->SetConstraintOptions((yyvsp[-3].copt)); (yyvsp[-4].MortarCondObj)->SetCtcMode((yyvsp[-1].ival)); domain->AddMortarCond((yyvsp[-4].MortarCondObj)); }
#line 10173 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 766:
#line 3238 "p.y" /* yacc.c:1646  */
    { (yyvsp[-3].MortarCondObj)->SetCtcMode((yyvsp[-1].ival)); domain->AddMortarCond((yyvsp[-3].MortarCondObj)); }
#line 10179 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 767:
#line 3244 "p.y" /* yacc.c:1646  */
    { domain->solInfo().trivial_detection = true;
          (yyval.MortarCondObj) = new MortarHandler(0, 0);
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType(MortarHandler::STD);
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode); }
#line 10189 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 768:
#line 3250 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-1].ival), (yyvsp[0].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType(MortarHandler::STD);
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10200 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 769:
#line 3257 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-2].ival), (yyvsp[-1].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[0].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10211 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 770:
#line 3264 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-3].ival), (yyvsp[-2].ival), (yyvsp[0].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-1].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10222 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 771:
#line 3271 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-1].fval), (yyvsp[0].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-2].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10233 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 772:
#line 3278 "p.y" /* yacc.c:1646  */
    { /* frictionless */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-3].fval), (yyvsp[-2].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[-1].ival), (yyvsp[0].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-4].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10245 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 773:
#line 3286 "p.y" /* yacc.c:1646  */
    { /* constant friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-7].ival), (yyvsp[-6].ival), (yyvsp[-4].fval), (yyvsp[-3].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[-2].ival), (yyvsp[-1].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[0].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-5].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10258 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 774:
#line 3295 "p.y" /* yacc.c:1646  */
    { /* velocity dependent friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-9].ival), (yyvsp[-8].ival), (yyvsp[-6].fval), (yyvsp[-5].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[-4].ival), (yyvsp[-3].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[0].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-7].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10271 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 775:
#line 3304 "p.y" /* yacc.c:1646  */
    { /* pressure dependent friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-10].ival), (yyvsp[-9].ival), (yyvsp[-7].fval), (yyvsp[-6].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[-5].ival), (yyvsp[-4].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[0].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-8].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10284 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 777:
#line 3316 "p.y" /* yacc.c:1646  */
    { domain->solInfo().dist_acme = (yyvsp[-1].ival); }
#line 10290 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 778:
#line 3318 "p.y" /* yacc.c:1646  */
    { domain->solInfo().no_secondary = true; }
#line 10296 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 779:
#line 3320 "p.y" /* yacc.c:1646  */
    { domain->solInfo().no_ghosting = true; }
#line 10302 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 780:
#line 3322 "p.y" /* yacc.c:1646  */
    { domain->solInfo().shell_simple_lofting = true; }
#line 10308 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 781:
#line 3324 "p.y" /* yacc.c:1646  */
    { domain->solInfo().no_multiple_interactions = true; }
#line 10314 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 782:
#line 3326 "p.y" /* yacc.c:1646  */
    { domain->solInfo().sharp_non_sharp_angle = (yyvsp[-1].fval); }
#line 10320 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 783:
#line 3328 "p.y" /* yacc.c:1646  */
    { domain->solInfo().normal_smoothing = false; }
#line 10326 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 784:
#line 3330 "p.y" /* yacc.c:1646  */
    { domain->solInfo().normal_smoothing_distance = (yyvsp[-1].fval); }
#line 10332 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 785:
#line 3332 "p.y" /* yacc.c:1646  */
    { domain->solInfo().resolution_method = (yyvsp[-1].ival); }
#line 10338 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 786:
#line 3334 "p.y" /* yacc.c:1646  */
    { domain->solInfo().old_dynamic_search = true; }
#line 10344 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 787:
#line 3336 "p.y" /* yacc.c:1646  */
    { domain->solInfo().partition_gap = true; }
#line 10350 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 788:
#line 3338 "p.y" /* yacc.c:1646  */
    { domain->solInfo().default_penalty = true; }
#line 10356 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 789:
#line 3340 "p.y" /* yacc.c:1646  */
    { domain->solInfo().global_search_cull = true; }
#line 10362 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 790:
#line 3342 "p.y" /* yacc.c:1646  */
    { domain->solInfo().no_warped_volume = true; }
#line 10368 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 791:
#line 3344 "p.y" /* yacc.c:1646  */
    { domain->solInfo().auto_tol = true; }
#line 10374 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 792:
#line 3346 "p.y" /* yacc.c:1646  */
    { domain->solInfo().auto_tol = true;
          domain->solInfo().agressive_tolerances = true; }
#line 10381 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 793:
#line 3349 "p.y" /* yacc.c:1646  */
    { domain->solInfo().skip_physical_faces = true; }
#line 10387 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 794:
#line 3351 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ffi_debug = bool((yyvsp[-1].ival)); }
#line 10393 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 795:
#line 3353 "p.y" /* yacc.c:1646  */
    { domain->solInfo().mortar_scaling = (yyvsp[-1].fval); }
#line 10399 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 796:
#line 3355 "p.y" /* yacc.c:1646  */
    { domain->solInfo().mortar_integration_rule = (yyvsp[-1].ival); }
#line 10405 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 797:
#line 3357 "p.y" /* yacc.c:1646  */
    { domain->solInfo().andes_clr = (yyvsp[-1].fval); }
#line 10411 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 798:
#line 3359 "p.y" /* yacc.c:1646  */
    { domain->solInfo().andes_cqr = (yyvsp[-1].fval); }
#line 10417 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 799:
#line 3361 "p.y" /* yacc.c:1646  */
    { domain->solInfo().andes_betab = (yyvsp[-1].fval); }
#line 10423 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 800:
#line 3363 "p.y" /* yacc.c:1646  */
    { domain->solInfo().andes_alpha = (yyvsp[-1].fval); }
#line 10429 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 801:
#line 3365 "p.y" /* yacc.c:1646  */
    { domain->solInfo().andes_betam = (yyvsp[-1].fval); }
#line 10435 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 802:
#line 3367 "p.y" /* yacc.c:1646  */
    { domain->solInfo().nlmembrane_pressure_type = (yyvsp[-1].ival); }
#line 10441 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 803:
#line 3370 "p.y" /* yacc.c:1646  */
    { geoSource->addNode((yyvsp[0].nval).num, (yyvsp[0].nval).xyz, (yyvsp[0].nval).cp, (yyvsp[0].nval).cd); }
#line 10447 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 804:
#line 3372 "p.y" /* yacc.c:1646  */
    { domain->solInfo().scalePosCoords = true;
          domain->solInfo().xScaleFactor = (yyvsp[-4].fval);
          domain->solInfo().yScaleFactor = (yyvsp[-3].fval);
          domain->solInfo().zScaleFactor = (yyvsp[-2].fval);
          geoSource->addNode((yyvsp[0].nval).num, (yyvsp[0].nval).xyz, (yyvsp[0].nval).cp, (yyvsp[0].nval).cd); }
#line 10457 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 805:
#line 3378 "p.y" /* yacc.c:1646  */
    { geoSource->addNode((yyvsp[0].nval).num, (yyvsp[0].nval).xyz, (yyvsp[0].nval).cp, (yyvsp[0].nval).cd); }
#line 10463 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 806:
#line 3382 "p.y" /* yacc.c:1646  */
    { (yyval.nval).num = (yyvsp[-4].ival)-1; (yyval.nval).xyz[0] = (yyvsp[-3].fval); (yyval.nval).xyz[1] = (yyvsp[-2].fval);  (yyval.nval).xyz[2] = (yyvsp[-1].fval);  (yyval.nval).cp = 0;  (yyval.nval).cd = 0; }
#line 10469 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 807:
#line 3384 "p.y" /* yacc.c:1646  */
    { (yyval.nval).num = (yyvsp[-3].ival)-1; (yyval.nval).xyz[0] = (yyvsp[-2].fval); (yyval.nval).xyz[1] = (yyvsp[-1].fval);  (yyval.nval).xyz[2] = 0.0; (yyval.nval).cp = 0;  (yyval.nval).cd = 0; }
#line 10475 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 808:
#line 3386 "p.y" /* yacc.c:1646  */
    { (yyval.nval).num = (yyvsp[-2].ival)-1; (yyval.nval).xyz[0] = (yyvsp[-1].fval); (yyval.nval).xyz[1] = 0.0; (yyval.nval).xyz[2] = 0.0; (yyval.nval).cp = 0;  (yyval.nval).cd = 0; }
#line 10481 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 809:
#line 3388 "p.y" /* yacc.c:1646  */
    { (yyval.nval).num = (yyvsp[-6].ival)-1; (yyval.nval).xyz[0] = (yyvsp[-5].fval); (yyval.nval).xyz[1] = (yyvsp[-4].fval);  (yyval.nval).xyz[2] = (yyvsp[-3].fval);  (yyval.nval).cp = (yyvsp[-2].ival); (yyval.nval).cd = (yyvsp[-1].ival);
          if((yyvsp[-2].ival) != 0) domain->solInfo().basicPosCoords = false;
          if((yyvsp[-1].ival) != 0) domain->solInfo().basicDofCoords = false; }
#line 10489 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 810:
#line 3392 "p.y" /* yacc.c:1646  */
    { (yyval.nval).num = (yyvsp[-5].ival)-1; (yyval.nval).xyz[0] = (yyvsp[-4].fval); (yyval.nval).xyz[1] = (yyvsp[-3].fval);  (yyval.nval).xyz[2] = (yyvsp[-2].fval);  (yyval.nval).cp = (yyvsp[-1].ival); (yyval.nval).cd = (yyvsp[-1].ival);
          if((yyvsp[-1].ival) != 0) { domain->solInfo().basicPosCoords = false; domain->solInfo().basicDofCoords = false; } }
#line 10496 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 811:
#line 3397 "p.y" /* yacc.c:1646  */
    { /* Define each Element */
          geoSource->addElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd);}
#line 10503 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 812:
#line 3402 "p.y" /* yacc.c:1646  */
    { (yyval.nl).num = 1; (yyval.nl).nd[0] = (yyvsp[0].ival)-1;}
#line 10509 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 813:
#line 3404 "p.y" /* yacc.c:1646  */
    { if((yyval.nl).num == 500) return -1; 
          (yyval.nl).nd[(yyval.nl).num] = (yyvsp[0].ival)-1; (yyval.nl).num++;}
#line 10516 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 814:
#line 3409 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-3].ival)-1; (yyval.bcval).dofnum = (yyvsp[-2].ival)-1; (yyval.bcval).val = (yyvsp[-1].fval); (yyval.bcval).mtype = BCond::Axial; }
#line 10522 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 815:
#line 3411 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-2].ival)-1; (yyval.bcval).dofnum = (yyvsp[-1].ival)-1; (yyval.bcval).val = 0.0; (yyval.bcval).mtype = BCond::Axial; }
#line 10528 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 816:
#line 3413 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-4].ival)-1; (yyval.bcval).dofnum = (yyvsp[-3].ival)-1; (yyval.bcval).val = (yyvsp[-2].fval); (yyval.bcval).mtype = (BCond::MomentType) (yyvsp[-1].ival); }
#line 10534 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 817:
#line 3415 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-3].ival)-1; (yyval.bcval).dofnum = (yyvsp[-2].ival)-1; (yyval.bcval).val = 0.0; (yyval.bcval).mtype = (BCond::MomentType) (yyvsp[-1].ival); }
#line 10540 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 818:
#line 3419 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-2].ival)-1;  (yyval.bcval).dofnum = -1;  (yyval.bcval).val = (yyvsp[-1].fval); }
#line 10546 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 819:
#line 3423 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-2].ival)-1; (yyval.bcval).dofnum = 6; (yyval.bcval).val = (yyvsp[-1].fval); }
#line 10552 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 820:
#line 3427 "p.y" /* yacc.c:1646  */
    { (yyval.cxbcval).nnum = (yyvsp[-4].ival)-1; (yyval.cxbcval).dofnum = (yyvsp[-3].ival)-1; (yyval.cxbcval).reval = (yyvsp[-2].fval); (yyval.cxbcval).imval = (yyvsp[-1].fval);  }
#line 10558 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 821:
#line 3429 "p.y" /* yacc.c:1646  */
    { (yyval.cxbcval).nnum = (yyvsp[-3].ival)-1; (yyval.cxbcval).dofnum = (yyvsp[-2].ival)-1; (yyval.cxbcval).reval = (yyvsp[-1].fval); (yyval.cxbcval).imval = 0.0; }
#line 10564 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 823:
#line 3434 "p.y" /* yacc.c:1646  */
    { geoSource->setCSFrame((yyvsp[0].frame).num,(yyvsp[0].frame).d); }
#line 10570 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 825:
#line 3439 "p.y" /* yacc.c:1646  */
    { geoSource->setFrame((yyvsp[0].frame).num,(yyvsp[0].frame).d);  }
#line 10576 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 826:
#line 3443 "p.y" /* yacc.c:1646  */
    { (yyval.frame).num = (yyvsp[-10].ival)-1; 
          (yyval.frame).d[0] = (yyvsp[-9].fval); (yyval.frame).d[1] = (yyvsp[-8].fval); (yyval.frame).d[2] = (yyvsp[-7].fval);
          (yyval.frame).d[3] = (yyvsp[-6].fval); (yyval.frame).d[4] = (yyvsp[-5].fval); (yyval.frame).d[5] = (yyvsp[-4].fval);
          (yyval.frame).d[6] = (yyvsp[-3].fval); (yyval.frame).d[7] = (yyvsp[-2].fval); (yyval.frame).d[8] = (yyvsp[-1].fval); }
#line 10585 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 827:
#line 3448 "p.y" /* yacc.c:1646  */
    { (yyval.frame).num = (yyvsp[-3].ival)-1;
          geoSource->makeEframe((yyvsp[-3].ival)-1, (yyvsp[-1].ival), (yyval.frame).d); }
#line 10592 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 829:
#line 3454 "p.y" /* yacc.c:1646  */
    { geoSource->setNodalFrame((yyvsp[0].nframe).id,(yyvsp[0].nframe).o,(yyvsp[0].nframe).d,(yyvsp[0].nframe).type); }
#line 10598 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 830:
#line 3458 "p.y" /* yacc.c:1646  */
    { (yyval.nframe).id = (yyvsp[-10].ival);
          (yyval.nframe).type = NFrameData::Rectangular;
          (yyval.nframe).o[0] = 0;  (yyval.nframe).o[1] = 0;  (yyval.nframe).o[2] = 0;
          (yyval.nframe).d[0] = (yyvsp[-9].fval); (yyval.nframe).d[1] = (yyvsp[-8].fval); (yyval.nframe).d[2] = (yyvsp[-7].fval);
          (yyval.nframe).d[3] = (yyvsp[-6].fval); (yyval.nframe).d[4] = (yyvsp[-5].fval); (yyval.nframe).d[5] = (yyvsp[-4].fval);
          (yyval.nframe).d[6] = (yyvsp[-3].fval); (yyval.nframe).d[7] = (yyvsp[-2].fval); (yyval.nframe).d[8] = (yyvsp[-1].fval); }
#line 10609 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 831:
#line 3465 "p.y" /* yacc.c:1646  */
    { (yyval.nframe).id = (yyvsp[-13].ival);
          (yyval.nframe).type = NFrameData::Rectangular;
          (yyval.nframe).o[0] = (yyvsp[-12].fval);  (yyval.nframe).o[1] = (yyvsp[-11].fval);  (yyval.nframe).o[2] = (yyvsp[-10].fval);
          (yyval.nframe).d[0] = (yyvsp[-9].fval);  (yyval.nframe).d[1] = (yyvsp[-8].fval);  (yyval.nframe).d[2] = (yyvsp[-7].fval);
          (yyval.nframe).d[3] = (yyvsp[-6].fval);  (yyval.nframe).d[4] = (yyvsp[-5].fval);  (yyval.nframe).d[5] = (yyvsp[-4].fval);
          (yyval.nframe).d[6] = (yyvsp[-3].fval); (yyval.nframe).d[7] = (yyvsp[-2].fval); (yyval.nframe).d[8] = (yyvsp[-1].fval); }
#line 10620 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 832:
#line 3472 "p.y" /* yacc.c:1646  */
    { (yyval.nframe).id = (yyvsp[-11].ival);
          (yyval.nframe).type = (yyvsp[-10].ival);
          (yyval.nframe).o[0] = 0;  (yyval.nframe).o[1] = 0;   (yyval.nframe).o[2] = 0;
          (yyval.nframe).d[0] = (yyvsp[-9].fval); (yyval.nframe).d[1] = (yyvsp[-8].fval);  (yyval.nframe).d[2] = (yyvsp[-7].fval);
          (yyval.nframe).d[3] = (yyvsp[-6].fval); (yyval.nframe).d[4] = (yyvsp[-5].fval);  (yyval.nframe).d[5] = (yyvsp[-4].fval);
          (yyval.nframe).d[6] = (yyvsp[-3].fval); (yyval.nframe).d[7] = (yyvsp[-2].fval); (yyval.nframe).d[8] = (yyvsp[-1].fval); }
#line 10631 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 833:
#line 3479 "p.y" /* yacc.c:1646  */
    { (yyval.nframe).id = (yyvsp[-14].ival);
          (yyval.nframe).type = (yyvsp[-13].ival);
          (yyval.nframe).o[0] = (yyvsp[-12].fval);  (yyval.nframe).o[1] = (yyvsp[-11].fval);  (yyval.nframe).o[2] = (yyvsp[-10].fval);
          (yyval.nframe).d[0] = (yyvsp[-9].fval);  (yyval.nframe).d[1] = (yyvsp[-8].fval);  (yyval.nframe).d[2] = (yyvsp[-7].fval);
          (yyval.nframe).d[3] = (yyvsp[-6].fval);  (yyval.nframe).d[4] = (yyvsp[-5].fval); (yyval.nframe).d[5] = (yyvsp[-4].fval);
          (yyval.nframe).d[6] = (yyvsp[-3].fval); (yyval.nframe).d[7] = (yyvsp[-2].fval); (yyval.nframe).d[8] = (yyvsp[-1].fval); }
#line 10642 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 835:
#line 3489 "p.y" /* yacc.c:1646  */
    { OffsetData od;
	  od.first = (yyvsp[-5].ival)-1; od.last = (yyvsp[-4].ival)-1;
	  od.o[0] = (yyvsp[-3].fval); od.o[1] = (yyvsp[-2].fval); od.o[2] = (yyvsp[-1].fval); 
	  geoSource->addOffset(od); }
#line 10651 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 836:
#line 3496 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 10657 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 837:
#line 3498 "p.y" /* yacc.c:1646  */
    { geoSource->setLocalIndex((yyvsp[-1].ival)-1); }
#line 10663 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 838:
#line 3501 "p.y" /* yacc.c:1646  */
    { geoSource->setElementLumpingWeight((yyvsp[-3].ival)-1,(yyvsp[-1].fval));
          domain->solInfo().elemLumpPodRom = true; }
#line 10670 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 839:
#line 3504 "p.y" /* yacc.c:1646  */
    { geoSource->setElementLumpingWeight((yyvsp[-4].ival)-1,(yyvsp[-2].fval));
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
#line 10678 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 840:
#line 3509 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-2].ival)-1,(yyvsp[-1].ival)-1); }
#line 10684 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 841:
#line 3511 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1); 
          geoSource->setElementLumpingWeight((yyvsp[-4].ival)-1,(yyvsp[-1].fval));
          domain->solInfo().elemLumpPodRom = true; }
#line 10692 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 842:
#line 3515 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[-5].ival)-1,(yyvsp[-1].fval));
          domain->solInfo().elemLumpPodRom = true; 
          domain->solInfo().reduceFollower = true; }
#line 10701 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 843:
#line 3520 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[-5].ival)-1,(yyvsp[-2].fval));
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
#line 10710 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 844:
#line 3526 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1,(yyvsp[-1].ival)-1); }
#line 10716 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 845:
#line 3528 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-6].ival)-1,(yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[-6].ival)-1,(yyvsp[-1].fval)); 
          domain->solInfo().elemLumpPodRom = true; }
#line 10724 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 846:
#line 3532 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-7].ival)-1,(yyvsp[-6].ival)-1,(yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[-7].ival)-1,(yyvsp[-2].fval));   
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
#line 10733 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 847:
#line 3538 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,-2,(yyvsp[-1].fval)); }
#line 10739 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 848:
#line 3540 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-7].ival)-1,(yyvsp[-6].ival)-1,(yyvsp[-5].ival)-1,-2,(yyvsp[-3].fval));
          geoSource->setElementLumpingWeight((yyvsp[-7].ival)-1,(yyvsp[-1].fval));
          domain->solInfo().elemLumpPodRom = true; }
#line 10747 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 849:
#line 3544 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-8].ival)-1,(yyvsp[-7].ival)-1,(yyvsp[-6].ival)-1,-2,(yyvsp[-4].fval)); 
          geoSource->setElementLumpingWeight((yyvsp[-8].ival)-1,(yyvsp[-2].fval));
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
#line 10756 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 850:
#line 3550 "p.y" /* yacc.c:1646  */
    { int i;
          for(i=(yyvsp[-3].ival); i<(yyvsp[-2].ival)+1; ++i)
            geoSource->setAttrib(i-1,i-1);
        }
#line 10765 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 851:
#line 3555 "p.y" /* yacc.c:1646  */
    { int i;
          for(i=(yyvsp[-3].ival); i<(yyvsp[-2].ival)+1; ++i)
            geoSource->setAttrib(i-1,(yyvsp[-1].ival)-1);
        }
#line 10774 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 852:
#line 3560 "p.y" /* yacc.c:1646  */
    { int i;
          for(i=(yyvsp[-5].ival); i<(yyvsp[-4].ival)+1; ++i)
            geoSource->setAttrib(i-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1);
        }
#line 10783 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 853:
#line 3565 "p.y" /* yacc.c:1646  */
    { int i;
          for(i=(yyvsp[-6].ival); i<(yyvsp[-5].ival)+1; ++i)
            geoSource->setAttrib(i-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, -2, (yyvsp[-1].fval));
        }
#line 10792 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 854:
#line 3572 "p.y" /* yacc.c:1646  */
    { domain->solInfo().elemLumpPodRom = true;
          geoSource->setLocalIndex(0); }
#line 10799 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 855:
#line 3575 "p.y" /* yacc.c:1646  */
    { domain->solInfo().elemLumpPodRom = true;
          geoSource->setLocalIndex((yyvsp[-1].ival)-1); }
#line 10806 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 856:
#line 3578 "p.y" /* yacc.c:1646  */
    { geoSource->setElementLumpingWeight((yyvsp[-2].ival) - 1, (yyvsp[-1].fval)); }
#line 10812 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 857:
#line 3580 "p.y" /* yacc.c:1646  */
    { domain->solInfo().reduceFollower = true;}
#line 10818 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 858:
#line 3582 "p.y" /* yacc.c:1646  */
    { domain->solInfo().reduceFollower = true;}
#line 10824 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 862:
#line 3591 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInLocalBasesAuxi[std::make_pair((yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1)] = std::string((yyvsp[-1].strval)); }
#line 10830 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 863:
#line 3595 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInLocalBasesCent.push_back(std::string((yyvsp[-1].strval))); }
#line 10836 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 864:
#line 3599 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ReducedStiffness = true;}
#line 10842 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 865:
#line 3601 "p.y" /* yacc.c:1646  */
    { geoSource->pushBackStiffVec((yyvsp[-1].fval));}
#line 10848 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 867:
#line 3606 "p.y" /* yacc.c:1646  */
    { domain->solInfo().forcePodSize = (yyvsp[-2].ival);
          domain->solInfo().maxDeimBasisSize = (yyvsp[-1].ival);}
#line 10855 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 868:
#line 3609 "p.y" /* yacc.c:1646  */
    { geoSource->pushBackUDEIMVec((yyvsp[-1].fval));}
#line 10861 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 870:
#line 3614 "p.y" /* yacc.c:1646  */
    { domain->solInfo().DEIMPodRom = true;
          geoSource->setSampleNodesAndSlots((yyvsp[-2].ival)-1,(yyvsp[-1].ival));}
#line 10868 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 871:
#line 3617 "p.y" /* yacc.c:1646  */
    { geoSource->setSampleElemsAndDOFs((yyvsp[-2].ival)-1,(yyvsp[-1].ival));
          domain->solInfo().UDEIMPodRom = true;}
#line 10875 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 872:
#line 3623 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 10881 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 873:
#line 3625 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-1].ival); }
#line 10887 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 874:
#line 3627 "p.y" /* yacc.c:1646  */
    { PressureBCond pbc;
          pbc.setData((yyvsp[-2].ival)-1, (yyvsp[-1].fval), (yyval.ival), true);
          geoSource->setElementPressure(pbc); }
#line 10895 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 875:
#line 3631 "p.y" /* yacc.c:1646  */
    { for(int i = (yyvsp[-3].ival); i < ((yyvsp[-2].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[-1].fval), (yyval.ival), true);
            geoSource->setElementPressure(pbc);
          } }
#line 10905 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 876:
#line 3637 "p.y" /* yacc.c:1646  */
    { PressureBCond *pbc = new PressureBCond[1];
          pbc[0].setData((yyvsp[-2].ival)-1, (yyvsp[-1].fval), (yyval.ival), true);
          geoSource->addSurfacePressure(1, pbc);
          if(geoSource->getNumSurfacePressure() > 1) delete [] pbc; }
#line 10914 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 877:
#line 3642 "p.y" /* yacc.c:1646  */
    { for(int i = (yyvsp[-4].ival); i <= (yyvsp[-2].ival); ++i) {
            PressureBCond *pbc = new PressureBCond[1];
            pbc[0].setData(i-1, (yyvsp[-1].fval), (yyval.ival), true);
            geoSource->addSurfacePressure(1, pbc);
            if(geoSource->getNumSurfacePressure() > 1) delete [] pbc;
          } }
#line 10925 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 878:
#line 3649 "p.y" /* yacc.c:1646  */
    { PressureBCond pbc;
          pbc.setData((yyvsp[-3].ival)-1, (yyvsp[-2].fval), (yyval.ival), (yyvsp[-1].ival));
          geoSource->setElementPressure(pbc); }
#line 10933 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 879:
#line 3653 "p.y" /* yacc.c:1646  */
    { for(int i = (yyvsp[-4].ival); i < ((yyvsp[-3].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[-2].fval), (yyval.ival), (yyvsp[-1].ival));
            geoSource->setElementPressure(pbc);
          } }
#line 10943 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 880:
#line 3659 "p.y" /* yacc.c:1646  */
    { PressureBCond *pbc = new PressureBCond[1];
          pbc[0].setData((yyvsp[-3].ival)-1, (yyvsp[-2].fval), (yyval.ival), (yyvsp[-1].ival));
          geoSource->addSurfacePressure(1, pbc);
          if(geoSource->getNumSurfacePressure() > 1) delete [] pbc; }
#line 10952 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 881:
#line 3665 "p.y" /* yacc.c:1646  */
    { PressureBCond pbc;
          pbc.setData((yyvsp[-4].ival)-1, (yyvsp[-1].fval), (yyval.ival), true);
          pbc.face = (yyvsp[-2].ival)-1;
          geoSource->setElementPressure(pbc); }
#line 10961 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 882:
#line 3670 "p.y" /* yacc.c:1646  */
    { for(int i = (yyvsp[-5].ival); i < ((yyvsp[-4].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[-1].fval), (yyval.ival), true);
            pbc.face = (yyvsp[-2].ival)-1;
            geoSource->setElementPressure(pbc);
          } }
#line 10972 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 883:
#line 3677 "p.y" /* yacc.c:1646  */
    { PressureBCond pbc;
          pbc.setData((yyvsp[-5].ival)-1, (yyvsp[-2].fval), (yyval.ival), (yyvsp[-1].ival));
          pbc.face = (yyvsp[-3].ival)-1;
          geoSource->setElementPressure(pbc); }
#line 10981 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 884:
#line 3682 "p.y" /* yacc.c:1646  */
    { for(int i = (yyvsp[-6].ival); i < ((yyvsp[-5].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[-2].fval), (yyval.ival), (yyvsp[-1].ival));
            pbc.face = (yyvsp[-3].ival)-1;
            geoSource->setElementPressure(pbc);
          } }
#line 10992 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 885:
#line 3691 "p.y" /* yacc.c:1646  */
    { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false); 
          geoSource->setConsistentPFlag(false); 
        }
#line 11001 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 886:
#line 3696 "p.y" /* yacc.c:1646  */
    { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false, (yyvsp[-1].ival));
          geoSource->setConsistentPFlag(false);
        }
#line 11010 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 887:
#line 3710 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useMassAugmentation = (yyvsp[-1].ival); }
#line 11016 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 888:
#line 3714 "p.y" /* yacc.c:1646  */
    { }
#line 11022 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 889:
#line 3716 "p.y" /* yacc.c:1646  */
    { geoSource->setElementPreLoad( (yyvsp[-2].ival)-1, (yyvsp[-1].fval) ); }
#line 11028 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 890:
#line 3718 "p.y" /* yacc.c:1646  */
    { int i;
          for(i=(yyvsp[-4].ival); i<((yyvsp[-2].ival)+1); ++i)
            geoSource->setElementPreLoad( i-1, (yyvsp[-1].fval) );
        }
#line 11037 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 891:
#line 3723 "p.y" /* yacc.c:1646  */
    { double load[3] = { (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval) };
          geoSource->setElementPreLoad( (yyvsp[-4].ival)-1, load ); }
#line 11044 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 892:
#line 3726 "p.y" /* yacc.c:1646  */
    { double load[3] = { (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval) };
          int i;
          for(i=(yyvsp[-6].ival); i<((yyvsp[-4].ival)+1); ++i)
            geoSource->setElementPreLoad( i-1, load );
        }
#line 11054 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 893:
#line 3734 "p.y" /* yacc.c:1646  */
    { domain->solInfo().sensitivity = true; }
#line 11060 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 894:
#line 3736 "p.y" /* yacc.c:1646  */
    { domain->solInfo().sensitivityMethod = (SolverInfo::SensitivityMethod) (yyvsp[-1].ival); }
#line 11066 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 895:
#line 3738 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readShapeSen = true;
          domain->solInfo().readInShapeSen = (yyvsp[-1].strval); }
#line 11073 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 896:
#line 3741 "p.y" /* yacc.c:1646  */
    { domain->solInfo().sensitivityTol = (yyvsp[-1].fval); }
#line 11079 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 897:
#line 3743 "p.y" /* yacc.c:1646  */
    { domain->solInfo().qsMaxvelSen = (yyvsp[-1].fval); }
#line 11085 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 898:
#line 3745 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ratioSensitivityTol = (yyvsp[-1].fval); }
#line 11091 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 899:
#line 3747 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ksParameter = (yyvsp[-1].fval); }
#line 11097 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 900:
#line 3749 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ksMax = (yyvsp[-1].fval); }
#line 11103 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 901:
#line 3751 "p.y" /* yacc.c:1646  */
    { }
#line 11109 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 902:
#line 3753 "p.y" /* yacc.c:1646  */
    { }
#line 11115 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 903:
#line 3755 "p.y" /* yacc.c:1646  */
    { }
#line 11121 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 904:
#line 3759 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11127 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 905:
#line 3761 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11133 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 906:
#line 3763 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11139 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 907:
#line 3765 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11145 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 908:
#line 3767 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11151 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 909:
#line 3769 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11157 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 910:
#line 3771 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11163 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 911:
#line 3773 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11169 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 912:
#line 3775 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-6].ival)-1, (yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11175 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 913:
#line 3777 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-6].ival)-1, (yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11181 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 914:
#line 3779 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-7].ival)-1, (yyvsp[-6].ival)-1, (yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11187 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 915:
#line 3781 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-7].ival)-1, (yyvsp[-6].ival)-1, (yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11193 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 916:
#line 3785 "p.y" /* yacc.c:1646  */
    { domain->setStressNodes((yyvsp[0].ival)); }
#line 11199 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 917:
#line 3787 "p.y" /* yacc.c:1646  */
    { domain->setStressNodes((yyvsp[0].ival)); }
#line 11205 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 918:
#line 3791 "p.y" /* yacc.c:1646  */
    { domain->setThicknessGroup((yyvsp[0].ival)); }
#line 11211 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 919:
#line 3793 "p.y" /* yacc.c:1646  */
    { domain->setThicknessGroup((yyvsp[0].ival)); }
#line 11217 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 920:
#line 3797 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Static); }
#line 11223 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 921:
#line 3799 "p.y" /* yacc.c:1646  */
    { domain->solInfo().solvercntl = (yyvsp[0].scntl); }
#line 11229 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 923:
#line 3802 "p.y" /* yacc.c:1646  */
    { // activate piecewise constant configuration dependent external forces for a linear dynamic analysis
          domain->solInfo().piecewise = true;
        }
#line 11237 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 924:
#line 3806 "p.y" /* yacc.c:1646  */
    { // activate piecewise constant configuration dependent external forces for a linear static analysis
          domain->solInfo().piecewise = true;
          domain->solInfo().piecewise_dlambda = (yyvsp[-2].fval);
          domain->solInfo().piecewise_maxLambda = (yyvsp[-1].fval);
        }
#line 11247 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 925:
#line 3812 "p.y" /* yacc.c:1646  */
    { domain->solInfo().coupled_scale = (yyvsp[-1].fval); }
#line 11253 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 926:
#line 3814 "p.y" /* yacc.c:1646  */
    { domain->sommerfeldType = (yyvsp[-1].ival);
          domain->curvatureFlag = 0; }
#line 11260 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 927:
#line 3817 "p.y" /* yacc.c:1646  */
    { domain->sommerfeldType = (yyvsp[-2].ival);
          domain->curvatureConst1 = (yyvsp[-1].fval);
          domain->curvatureFlag = 1; }
#line 11268 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 928:
#line 3821 "p.y" /* yacc.c:1646  */
    { domain->sommerfeldType = (yyvsp[-3].ival);
          domain->curvatureConst1 = (yyvsp[-2].fval);
          domain->curvatureConst2 = (yyvsp[-1].fval);
          domain->curvatureFlag = 2; }
#line 11277 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 929:
#line 3826 "p.y" /* yacc.c:1646  */
    { domain->solInfo().dmpc = bool((yyvsp[-1].ival)); }
#line 11283 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 930:
#line 3828 "p.y" /* yacc.c:1646  */
    { domain->solInfo().dbccheck = bool((yyvsp[-1].ival)); }
#line 11289 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 931:
#line 3832 "p.y" /* yacc.c:1646  */
    { domain->solInfo().loadcases.push_back((yyvsp[0].ival)); }
#line 11295 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 932:
#line 3834 "p.y" /* yacc.c:1646  */
    { domain->solInfo().loadcases.push_back((yyvsp[0].ival)); }
#line 11301 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 933:
#line 3838 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
          (yyval.scntl)->type = SolverSelection::Direct;
          (yyval.scntl)->subtype = 0; }
#line 11309 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 934:
#line 3842 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
          (yyval.scntl)->type = SolverSelection::Direct;
          (yyval.scntl)->subtype = (yyvsp[-1].ival); }
#line 11317 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 935:
#line 3846 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
          (yyval.scntl)->type = SolverSelection::Direct;
          (yyval.scntl)->subtype = (yyvsp[-1].ival); }
#line 11325 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 936:
#line 3850 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Direct;
      (yyval.scntl)->subtype = (yyvsp[-2].ival);
      (yyval.scntl)->pivot = true; }
#line 11334 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 937:
#line 3855 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Direct;
      (yyval.scntl)->subtype = (yyvsp[-2].ival);
      (yyval.scntl)->scaled = true; }
#line 11343 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 938:
#line 3860 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Direct;
      (yyval.scntl)->subtype = (yyvsp[-2].ival);
      (yyval.scntl)->unsymmetric = true; }
#line 11352 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 939:
#line 3865 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-1].ival); }
#line 11360 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 940:
#line 3869 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-2].ival);
      (yyval.scntl)->precond = (yyvsp[-1].ival); }
#line 11369 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 941:
#line 3874 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-3].ival);
      (yyval.scntl)->precond = (yyvsp[-2].ival);
      (yyval.scntl)->tol=(yyvsp[-1].fval); }
#line 11379 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 942:
#line 3880 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-4].ival);
      (yyval.scntl)->precond = (yyvsp[-3].ival);
      (yyval.scntl)->tol = (yyvsp[-2].fval);
      (yyval.scntl)->maxit = (yyvsp[-1].ival); }
#line 11390 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 943:
#line 3887 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-5].ival);
      (yyval.scntl)->precond = (yyvsp[-4].ival);
      (yyval.scntl)->tol = (yyvsp[-3].fval);
      (yyval.scntl)->maxit = (yyvsp[-2].ival);
      (yyval.scntl)->iterSubtype = (yyvsp[-1].ival); }
#line 11402 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 944:
#line 3895 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-6].ival);
      (yyval.scntl)->precond = (yyvsp[-5].ival);
      (yyval.scntl)->tol = (yyvsp[-4].fval);
      (yyval.scntl)->maxit = (yyvsp[-3].ival);
      (yyval.scntl)->iterSubtype = (yyvsp[-2].ival);
      (yyval.scntl)->maxvecsize = (yyvsp[-1].ival); }
#line 11415 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 945:
#line 3904 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti; }
#line 11422 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 946:
#line 3907 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.version = (FetiInfo::Version) ((yyvsp[-1].ival)-1); }
#line 11430 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 947:
#line 3911 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.maxit = (yyvsp[-2].ival);
      (yyval.scntl)->fetiInfo.tol = (yyvsp[-1].fval);
      (yyval.scntl)->fetiInfo.maxortho = (yyvsp[-2].ival); }
#line 11440 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 948:
#line 3917 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.maxit = (yyvsp[-3].ival);
      (yyval.scntl)->fetiInfo.tol = (yyvsp[-2].fval);
      (yyval.scntl)->fetiInfo.maxortho = (yyvsp[-1].ival); }
#line 11450 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 949:
#line 3923 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.version = (FetiInfo::Version) ((yyvsp[-2].ival)-1);
      (yyval.scntl)->fetiInfo.feti2version = (FetiInfo::Feti2Version) (yyvsp[-1].ival); }
#line 11459 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 950:
#line 3928 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp; }
#line 11468 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 951:
#line 3933 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::FetiLib;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp; }
#line 11477 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 952:
#line 3938 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp;
      (yyval.scntl)->fetiInfo.dph_flag = true; }
#line 11487 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 953:
#line 3944 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.scaling = FetiInfo::tscaling;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners3;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp;
      (yyval.scntl)->fetiInfo.dph_flag = true;
      (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges;
      (yyval.scntl)->fetiInfo.rbmType = FetiInfo::None;
      (yyval.scntl)->fetiInfo.nGs = 0; }
#line 11501 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 954:
#line 3954 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.maxit = (yyvsp[-2].ival);
      (yyval.scntl)->fetiInfo.tol = (yyvsp[-1].fval);
      (yyval.scntl)->fetiInfo.krylovtype = 1;
      (yyval.scntl)->fetiInfo.scaling = FetiInfo::tscaling;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp;
      (yyval.scntl)->fetiInfo.dph_flag = true;
      (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges;
      (yyval.scntl)->fetiInfo.rbmType = FetiInfo::None;
      (yyval.scntl)->fetiInfo.nGs = 0; }
#line 11518 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 955:
#line 3967 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.maxit = (yyvsp[-3].ival);
      (yyval.scntl)->fetiInfo.tol = (yyvsp[-2].fval);
      (yyval.scntl)->fetiInfo.numcgm = (yyvsp[-1].ival);
      (yyval.scntl)->fetiInfo.krylovtype = 1;
      (yyval.scntl)->fetiInfo.scaling = FetiInfo::tscaling;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp;
      (yyval.scntl)->fetiInfo.dph_flag = true;
      (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges;
      (yyval.scntl)->fetiInfo.rbmType = FetiInfo::None;
      (yyval.scntl)->fetiInfo.nGs = 0; }
#line 11536 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 956:
#line 3981 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.maxit = (yyvsp[-4].ival);
      (yyval.scntl)->fetiInfo.tol = (yyvsp[-3].fval);
      (yyval.scntl)->fetiInfo.numcgm = (yyvsp[-2].ival);
      (yyval.scntl)->fetiInfo.tolcgm = (yyvsp[-1].fval);
      (yyval.scntl)->fetiInfo.spaceDimension = 2;
      (yyval.scntl)->fetiInfo.krylovtype = 1;
      (yyval.scntl)->fetiInfo.scaling = FetiInfo::tscaling;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp;
      (yyval.scntl)->fetiInfo.dph_flag = true;
      (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges;
      (yyval.scntl)->fetiInfo.rbmType = FetiInfo::None;
      (yyval.scntl)->fetiInfo.nGs = 0; }
#line 11556 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 957:
#line 3997 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.maxit = (yyvsp[-5].ival);
      (yyval.scntl)->fetiInfo.tol = (yyvsp[-4].fval);
      (yyval.scntl)->fetiInfo.numcgm = (yyvsp[-3].ival);
      (yyval.scntl)->fetiInfo.tolcgm = (yyvsp[-2].fval);
      (yyval.scntl)->fetiInfo.spaceDimension = (yyvsp[-1].ival);
      (yyval.scntl)->fetiInfo.krylovtype = 1;
      (yyval.scntl)->fetiInfo.scaling = FetiInfo::tscaling;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp;
      (yyval.scntl)->fetiInfo.dph_flag = true;
      (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges;
      (yyval.scntl)->fetiInfo.rbmType = FetiInfo::None;
      (yyval.scntl)->fetiInfo.nGs = 0; }
#line 11576 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 958:
#line 4013 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::BlockDiag; }
#line 11583 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 959:
#line 4016 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = &domain->solInfo().solvercntls[(yyvsp[-1].ival)]; }
#line 11589 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 960:
#line 4020 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = &domain->solInfo().solvercntls[(yyvsp[-2].ival)];
          *(yyval.scntl) = *(yyvsp[0].scntl); }
#line 11596 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 961:
#line 4024 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = (yyvsp[0].scntl); }
#line 11602 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 962:
#line 4026 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->verbose = 1; }
#line 11608 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 963:
#line 4028 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->verbose = (yyvsp[-1].ival); }
#line 11614 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 964:
#line 4030 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.printNumber = (yyvsp[-1].ival); }
#line 11620 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 965:
#line 4032 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->trbm = (yyvsp[-1].fval); }
#line 11626 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 966:
#line 4034 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->sparse_renum = (yyvsp[-1].ival); }
#line 11632 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 967:
#line 4036 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->sparse_maxsup = (yyvsp[-1].ival); }
#line 11638 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 968:
#line 4038 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->sparse_defblk = (yyvsp[-1].ival); }
#line 11644 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 969:
#line 4040 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_tau = (yyvsp[-1].fval); }
#line 11650 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 970:
#line 4042 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_maxsize = (yyvsp[-1].ival); }
#line 11656 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 971:
#line 4044 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].ival) < 0) {
            (yyvsp[-1].ival) = 24;
            fprintf(stderr," *** WARNING: spooles_maxdomainsize must be > 0,"
                           " using 24\n");
          }
          (yyval.scntl)->spooles_maxdomainsize = (yyvsp[-1].ival); }
#line 11667 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 972:
#line 4051 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_seed = (yyvsp[-1].ival); }
#line 11673 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 973:
#line 4053 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].fval) < 0.0) || ((yyvsp[-1].fval) > 1.0)) {
            (yyvsp[-1].fval) = 0.04;
            fprintf(stderr," *** WARNING: spooles_maxzeros outside acceptable limits (0..1),"
                           " using 0.04\n");
          }
          (yyval.scntl)->spooles_maxzeros = (yyvsp[-1].fval); }
#line 11684 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 974:
#line 4060 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_msglvl = (yyvsp[-1].ival); }
#line 11690 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 975:
#line 4062 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_scale = (yyvsp[-1].ival); }
#line 11696 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 976:
#line 4064 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_renum = (yyvsp[-1].ival); }
#line 11702 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 977:
#line 4066 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->mumps_icntl[(yyvsp[-2].ival)] = (yyvsp[-1].ival); }
#line 11708 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 978:
#line 4068 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->mumps_cntl[(yyvsp[-2].ival)] = (yyvsp[-1].fval); }
#line 11714 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 979:
#line 4070 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->mumps_mineq = (yyvsp[-1].ival); }
#line 11720 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 980:
#line 4072 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->mumps_stride = (yyvsp[-1].ival); }
#line 11726 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 981:
#line 4074 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->goldfarb_tol = (yyvsp[-1].fval); }
#line 11732 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 982:
#line 4076 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->tol = (yyvsp[-1].fval); }
#line 11738 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 983:
#line 4078 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->iterSubtype = (yyvsp[-1].ival); }
#line 11744 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 984:
#line 4080 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->ilu_droptol = (yyvsp[-1].fval); }
#line 11750 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 985:
#line 4082 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->maxit = (yyvsp[-1].ival); 
          (yyval.scntl)->fetiInfo.maxit = (yyvsp[-1].ival); }
#line 11757 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 986:
#line 4085 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->maxvecsize = (yyvsp[-1].ival);
          (yyval.scntl)->fetiInfo.maxortho = (yyvsp[-1].ival); }
#line 11764 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 987:
#line 4088 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[-1].ival); }
#line 11770 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 988:
#line 4090 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[-1].ival); }
#line 11776 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 989:
#line 4092 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.precno = FetiInfo::lumped; }
#line 11782 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 990:
#line 4094 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->precond = (yyvsp[-1].ival);
          if(isFeti((yyval.scntl)->type) && (((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 3))) {
            (yyvsp[-1].ival) = 1;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner selected, using lumped\n");
          }
          (yyval.scntl)->fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[-1].ival); }
#line 11793 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 991:
#line 4101 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 1)) {
            (yyvsp[-1].ival) = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          (yyval.scntl)->fetiInfo.prectype = (FetiInfo::PreconditionerType) (yyvsp[-1].ival); }
#line 11803 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 992:
#line 4107 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 1)) {
            (yyvsp[-1].ival) = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          (yyval.scntl)->fetiInfo.prectype = (FetiInfo::PreconditionerType) (yyvsp[-1].ival); }
#line 11813 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 993:
#line 4113 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.tol = (yyvsp[-1].fval); }
#line 11819 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 994:
#line 4115 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.tol = (yyvsp[-2].fval); 
          (yyval.scntl)->fetiInfo.absolute_tol = (yyvsp[-1].fval); }
#line 11826 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 995:
#line 4118 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.stagnation_tol = (yyvsp[-1].fval); }
#line 11832 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 996:
#line 4120 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.stagnation_tol = (yyvsp[-2].fval);
          (yyval.scntl)->fetiInfo.absolute_stagnation_tol = (yyvsp[-1].fval); }
#line 11839 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 997:
#line 4123 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.primal_proj_tol = (yyvsp[-2].fval);
          (yyval.scntl)->fetiInfo.dual_proj_tol = (yyvsp[-1].fval); }
#line 11846 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 998:
#line 4126 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.primal_plan_maxit = (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.dual_plan_maxit = (yyvsp[-1].ival); }
#line 11853 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 999:
#line 4129 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.primal_plan_tol = (yyvsp[-2].fval);
          (yyval.scntl)->fetiInfo.dual_plan_tol = (yyvsp[-1].fval); }
#line 11860 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1000:
#line 4132 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.noCoarse = 1; }
#line 11866 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1001:
#line 4134 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].ival) == 1) 
            (yyval.scntl)->fetiInfo.nonLocalQ = 0;
          else if((yyvsp[-1].ival) == 2) {
            (yyval.scntl)->fetiInfo.nonLocalQ = 1;
            if((yyval.scntl)->fetiInfo.version == FetiInfo::feti2) {
              (yyval.scntl)->fetiInfo.nonLocalQ = 0;
              fprintf(stderr," *** WARNING: Basic projector is used"
                             " with FETI 2\n");
            }
          } else if((yyvsp[-1].ival) == 3) {
            (yyval.scntl)->fetiInfo.nonLocalQ = 1;
            (yyval.scntl)->fetiInfo.nQ = 3;
          } else if((yyvsp[-1].ival) == 4) {
            (yyval.scntl)->fetiInfo.nonLocalQ = 1;
            (yyval.scntl)->fetiInfo.nQ = 4;
          } else
            fprintf(stderr," *** WARNING: This projector does not exist,"
                           " using basic projector\n"); }
#line 11889 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1002:
#line 4153 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 2)) (yyvsp[-1].ival) = 1; 
          (yyval.scntl)->fetiInfo.scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11896 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1003:
#line 4156 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11902 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1004:
#line 4158 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 2)) (yyvsp[-1].ival) = 2;
          (yyval.scntl)->fetiInfo.mpc_scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11909 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1005:
#line 4161 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11915 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1006:
#line 4163 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 2)) (yyvsp[-1].ival) = 2;
          (yyval.scntl)->fetiInfo.fsi_scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11922 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1007:
#line 4166 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.fsi_scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11928 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1008:
#line 4168 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_element = true; }
#line 11934 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1009:
#line 4170 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.fsi_element = true; }
#line 11940 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1010:
#line 4172 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.fsi_corner = (yyvsp[-1].ival); }
#line 11946 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1011:
#line 4174 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.splitLocalFsi = false; }
#line 11952 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1012:
#line 4176 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.wetcorners = true; }
#line 11958 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1013:
#line 4178 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.corners = (FetiInfo::CornerType) (yyvsp[-1].ival); }
#line 11964 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1014:
#line 4180 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.corners = (FetiInfo::CornerType) (yyvsp[-2].ival); 
          (yyval.scntl)->fetiInfo.pick_unsafe_corners = bool((yyvsp[-1].ival)); }
#line 11971 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1015:
#line 4183 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].ival) == 0) {
            (yyval.scntl)->fetiInfo.corners = FetiInfo::noCorners;
            (yyval.scntl)->fetiInfo.pickAnyCorner = 0; 
            (yyval.scntl)->fetiInfo.bmpc = true;
            (yyval.scntl)->fetiInfo.pick_unsafe_corners = false;
            (yyval.scntl)->fetiInfo.augment = FetiInfo::none;
          } }
#line 11983 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1016:
#line 4191 "p.y" /* yacc.c:1646  */
    { if((yyval.scntl)->fetiInfo.dph_flag && ((yyvsp[-1].ival) == 1)) {
            std::cerr << "WARNING: Selected augment type is unsupported for FETI-DPH, set to EdgeGs \n";
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges;
          }
          else (yyval.scntl)->fetiInfo.augment = (FetiInfo::AugmentType) (yyvsp[-1].ival);

          if((yyval.scntl)->fetiInfo.augment == FetiInfo::Edges) {
            (yyval.scntl)->fetiInfo.rbmType = FetiInfo::translation;
            (yyval.scntl)->fetiInfo.nGs = 3;
          } 
          else if((yyval.scntl)->fetiInfo.augment == FetiInfo::WeightedEdges) {
            (yyval.scntl)->fetiInfo.rbmType = FetiInfo::translation;
            (yyval.scntl)->fetiInfo.nGs = 3;
          }
          else if((yyval.scntl)->fetiInfo.augment == FetiInfo::Gs) {
            (yyval.scntl)->fetiInfo.rbmType = FetiInfo::all;
            (yyval.scntl)->fetiInfo.nGs = 6;
          } }
#line 12006 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1017:
#line 4210 "p.y" /* yacc.c:1646  */
    { if((yyval.scntl)->fetiInfo.dph_flag && ((yyvsp[-2].ival) == 1)) {
            std::cerr << "WARNING: Selected augment type is unsupported for FETI-DPH, set to EdgeGs \n";
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges;
          }
          else (yyval.scntl)->fetiInfo.augment = (FetiInfo::AugmentType) (yyvsp[-2].ival);
          if((yyval.scntl)->fetiInfo.dph_flag && ((yyvsp[-1].ival) > 2) && ((yyvsp[-1].ival) < 6)) {
            std::cerr << "WARNING: Selected rbm type is unsupported for FETI-DPH, set to translation \n";
            (yyval.scntl)->fetiInfo.rbmType = FetiInfo::translation;
          }
          else (yyval.scntl)->fetiInfo.rbmType = (FetiInfo::RbmType) (yyvsp[-1].ival);

          if((yyval.scntl)->fetiInfo.rbmType == FetiInfo::all)
            (yyval.scntl)->fetiInfo.nGs = 6;
          else if((yyval.scntl)->fetiInfo.rbmType != FetiInfo::None)
            (yyval.scntl)->fetiInfo.nGs = 3; }
#line 12026 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1018:
#line 4226 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.numdir = (yyvsp[-1].ival); 
          if((yyval.scntl)->fetiInfo.augment == FetiInfo::none)
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges; }
#line 12034 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1019:
#line 4230 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.waveType = (FetiInfo::WaveType) (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.numdir = (yyvsp[-1].ival); 
          if((yyval.scntl)->fetiInfo.augment == FetiInfo::none)
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges; }
#line 12043 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1020:
#line 4235 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.numdir = (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.waveMethod = (FetiInfo::WaveMethod) (yyvsp[-1].ival);
          if((yyval.scntl)->fetiInfo.augment == FetiInfo::none)
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges; }
#line 12052 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1021:
#line 4240 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.waveType = (FetiInfo::WaveType) (yyvsp[-3].ival);
          (yyval.scntl)->fetiInfo.waveMethod = (FetiInfo::WaveMethod) (yyvsp[-1].ival);
          (yyval.scntl)->fetiInfo.numdir = (yyvsp[-2].ival);
          if((yyval.scntl)->fetiInfo.augment == FetiInfo::none)
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges; }
#line 12062 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1022:
#line 4246 "p.y" /* yacc.c:1646  */
    { if ((yyvsp[-1].ival) == 2)
           (yyval.scntl)->fetiInfo.augmentimpl = FetiInfo::Primal;
        }
#line 12070 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1023:
#line 4250 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.orthotol = (yyvsp[-1].fval); }
#line 12076 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1024:
#line 4252 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.orthotol = (yyvsp[-2].fval); 
          (yyval.scntl)->fetiInfo.orthotol2 = (yyvsp[-1].fval); }
#line 12083 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1025:
#line 4255 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.grbm_tol = (yyvsp[-1].fval); }
#line 12089 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1026:
#line 4257 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.crbm_tol = (yyvsp[-1].fval); }
#line 12095 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1027:
#line 4259 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.cct_tol = (yyvsp[-1].fval); }
#line 12101 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1028:
#line 4261 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.rebuildcct = int((yyvsp[-1].ival)); }
#line 12107 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1029:
#line 4263 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.uproj = (yyvsp[-1].ival); }
#line 12113 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1030:
#line 4265 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.printMatLab = 1; }
#line 12119 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1031:
#line 4267 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->printMatLab = 1;
          (yyval.scntl)->printMatLabFile = (yyvsp[-1].strval); }
#line 12126 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1032:
#line 4270 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.local_cntl = (yyval.scntl)->fetiInfo.kii_cntl = (yyvsp[0].scntl); }
#line 12132 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1033:
#line 4272 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.coarse_cntl = (yyvsp[0].scntl); }
#line 12138 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1034:
#line 4274 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.auxcoarse_cntl = (yyvsp[0].scntl); }
#line 12144 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1035:
#line 4276 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.cct_cntl = (yyvsp[0].scntl); }
#line 12150 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1036:
#line 4278 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].ival) == 1)
            (yyval.scntl)->fetiInfo.version = FetiInfo::feti1; 
          else if((yyvsp[-1].ival) == 2) {
            (yyval.scntl)->fetiInfo.version = FetiInfo::feti2;
            if((yyval.scntl)->fetiInfo.nonLocalQ == 1) {
              (yyval.scntl)->fetiInfo.nonLocalQ = 0;
              (yyval.scntl)->fetiInfo.nQ = 0;
              fprintf(stderr," *** WARNING: Basic projector is used "
                             "with FETI 2\n");
            }
          } else if((yyvsp[-1].ival) == 3) {
            (yyval.scntl)->fetiInfo.version = FetiInfo::feti3;
            if((yyval.scntl)->fetiInfo.nonLocalQ == 1) {
              (yyval.scntl)->fetiInfo.nonLocalQ = 0;
              fprintf(stderr," *** WARNING: Basic projector is used "
                             "with FETI 2\n");
            }
          } 
          else {
            (yyval.scntl)->fetiInfo.version = FetiInfo::feti1;
	    fprintf(stderr," *** WARNING: Version does not exist,"
                           " using FETI 1\n");
          } }
#line 12178 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1037:
#line 4302 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.gmresResidual = true; }
#line 12184 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1038:
#line 4304 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.gmresResidual = bool((yyvsp[-1].ival)); }
#line 12190 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1039:
#line 4306 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.pickAnyCorner = (yyvsp[-1].ival); }
#line 12196 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1040:
#line 4308 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.type = FetiInfo::nonlinear;
          (yyval.scntl)->fetiInfo.nlPrecFlg = 1; }
#line 12203 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1041:
#line 4311 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.type = FetiInfo::nonlinear;
	  (yyval.scntl)->fetiInfo.nlPrecFlg = (yyvsp[-1].ival); }
#line 12210 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1042:
#line 4314 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.numcgm = (yyvsp[-1].ival); }
#line 12216 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1043:
#line 4316 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.tolcgm = (yyvsp[-1].fval); }
#line 12222 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1044:
#line 4318 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.spaceDimension = (yyvsp[-1].ival); }
#line 12228 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1045:
#line 4320 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.krylovtype = (yyvsp[-1].ival); }
#line 12234 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1046:
#line 4322 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.krylovtype = (yyvsp[-1].ival); }
#line 12240 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1047:
#line 4324 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.lumpedinterface = 1; }
#line 12246 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1048:
#line 4326 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.saveMemCoarse = 1; }
#line 12252 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1049:
#line 4328 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.outerloop = (FetiInfo::OuterloopType) (yyvsp[-1].ival); }
#line 12258 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1050:
#line 4330 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.outerloop = (FetiInfo::OuterloopType) (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.complex_hermitian = true; }
#line 12265 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1051:
#line 4333 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpcflag = (yyvsp[-1].ival); }
#line 12271 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1052:
#line 4335 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpcflag = (yyvsp[-1].ival); }
#line 12277 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1053:
#line 4337 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-1].ival); }
#line 12283 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1054:
#line 4339 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-1].ival); }
#line 12289 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1055:
#line 4341 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-3].ival);
          (yyval.scntl)->fetiInfo.mpcBlkOverlap = (yyvsp[-1].ival); }
#line 12296 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1056:
#line 4344 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-3].ival);
          (yyval.scntl)->fetiInfo.mpcBlkOverlap = (yyvsp[-1].ival); }
#line 12303 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1057:
#line 4347 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-2].ival); 
          (yyval.scntl)->fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[-1].ival); }
#line 12310 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1058:
#line 4350 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[-1].ival); }
#line 12317 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1059:
#line 4353 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-4].ival);
          (yyval.scntl)->fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[-3].ival);
          (yyval.scntl)->fetiInfo.mpcBlkOverlap = (yyvsp[-1].ival); }
#line 12325 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1060:
#line 4357 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-4].ival);
          (yyval.scntl)->fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[-3].ival); 
          (yyval.scntl)->fetiInfo.mpcBlkOverlap = (yyvsp[-1].ival); }
#line 12333 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1061:
#line 4361 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].ival) < 1) (yyval.scntl)->fetiInfo.useMRHS = false; }
#line 12339 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1062:
#line 4363 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.gamma = (yyvsp[-1].fval); }
#line 12345 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1063:
#line 4365 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.linesearch_maxit = (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.linesearch_tau = (yyvsp[-1].fval); }
#line 12352 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1064:
#line 4368 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.bmpc = bool((yyvsp[-1].ival)); }
#line 12358 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1065:
#line 4370 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.cmpc = bool((yyvsp[-1].ival)); }
#line 12364 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1066:
#line 4372 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.c_normalize = bool((yyvsp[-1].ival)); }
#line 12370 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1068:
#line 4379 "p.y" /* yacc.c:1646  */
    {
          geoSource->setOmega((yyvsp[-1].fval));
          StructProp sp; 
          sp.kappaHelm = (yyvsp[-1].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::Helmholtz);
        }
#line 12383 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1069:
#line 4390 "p.y" /* yacc.c:1646  */
    { if(!(yyvsp[-1].copt).lagrangeMult && (yyvsp[-1].copt).penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = (yyvsp[-1].copt).lagrangeMult;
          domain->solInfo().penalty = (yyvsp[-1].copt).penalty;
          domain->solInfo().constraint_hess = (yyvsp[-1].copt).constraint_hess; 
          domain->solInfo().constraint_hess_eps = (yyvsp[-1].copt).constraint_hess_eps; }
#line 12393 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1070:
#line 4396 "p.y" /* yacc.c:1646  */
    { if(!(yyvsp[-1].copt).lagrangeMult && (yyvsp[-1].copt).penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = (yyvsp[-1].copt).lagrangeMult;
          domain->solInfo().penalty = (yyvsp[-1].copt).penalty;
          domain->solInfo().constraint_hess = (yyvsp[-1].copt).constraint_hess;
          domain->solInfo().constraint_hess_eps = (yyvsp[-1].copt).constraint_hess_eps; }
#line 12403 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1071:
#line 4404 "p.y" /* yacc.c:1646  */
    { // Direct elimination of slave dofs
          (yyval.copt).lagrangeMult = false;
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
        }
#line 12414 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1072:
#line 4411 "p.y" /* yacc.c:1646  */
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[0].fval); }
#line 12425 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1073:
#line 4418 "p.y" /* yacc.c:1646  */
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[-1].fval);
          domain->solInfo().coefFilterTol = (yyvsp[0].fval); }
#line 12437 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1074:
#line 4426 "p.y" /* yacc.c:1646  */
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[-2].fval); 
          domain->solInfo().coefFilterTol = (yyvsp[-1].fval);
          domain->solInfo().rhsZeroTol = (yyvsp[0].fval); }
#line 12450 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1075:
#line 4435 "p.y" /* yacc.c:1646  */
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[-3].fval);
          domain->solInfo().coefFilterTol = (yyvsp[-2].fval); 
          domain->solInfo().rhsZeroTol = (yyvsp[-1].fval);
          domain->solInfo().inconsistentTol = (yyvsp[0].fval); }
#line 12464 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1076:
#line 4445 "p.y" /* yacc.c:1646  */
    { // Treatment of constraints through Lagrange multipliers method
          (yyval.copt).lagrangeMult = true; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; }
#line 12474 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1077:
#line 4451 "p.y" /* yacc.c:1646  */
    { // Treatment of constraints through penalty method
          (yyval.copt).lagrangeMult = false;
          (yyval.copt).penalty = (yyvsp[0].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; }
#line 12484 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1078:
#line 4457 "p.y" /* yacc.c:1646  */
    { // Treatment of constraints through augmented Lagrangian method
          (yyval.copt).lagrangeMult = true;
          (yyval.copt).penalty = (yyvsp[0].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; }
#line 12494 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1079:
#line 4463 "p.y" /* yacc.c:1646  */
    { // Alternative input syntax for treatment of constraints through augmented Lagrangian method
          (yyval.copt).lagrangeMult = true;
          (yyval.copt).penalty = (yyvsp[0].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; }
#line 12504 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1080:
#line 4469 "p.y" /* yacc.c:1646  */
    { (yyval.copt).constraint_hess = (yyvsp[0].ival);
          (yyval.copt).constraint_hess_eps = 0; }
#line 12511 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1081:
#line 4472 "p.y" /* yacc.c:1646  */
    { (yyval.copt).constraint_hess = (yyvsp[-1].ival);
          (yyval.copt).constraint_hess_eps = (yyvsp[0].fval); }
#line 12518 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1082:
#line 4477 "p.y" /* yacc.c:1646  */
    { // hack??
	  domain->solInfo().acoustic = true;
          if(domain->solInfo().probType != SolverInfo::HelmholtzDirSweep) domain->solInfo().setProbType(SolverInfo::Helmholtz);
        }
#line 12527 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1083:
#line 4494 "p.y" /* yacc.c:1646  */
    {
          domain->sommerfeldType = (yyvsp[-1].ival);
          domain->curvatureFlag = 0;
        }
#line 12536 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1084:
#line 4499 "p.y" /* yacc.c:1646  */
    {
          domain->sommerfeldType = (yyvsp[-2].ival);
          domain->curvatureConst1 = (yyvsp[-1].fval);
          domain->curvatureFlag = 1;
        }
#line 12546 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1085:
#line 4505 "p.y" /* yacc.c:1646  */
    {
          domain->sommerfeldType = (yyvsp[-3].ival);
          domain->curvatureConst1 = (yyvsp[-2].fval);
          domain->curvatureConst2 = (yyvsp[-1].fval);
          domain->curvatureFlag = 2;
        }
#line 12557 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1086:
#line 4512 "p.y" /* yacc.c:1646  */
    {
          domain->pointSourceFlag = 1;
          domain->implicitFlag = 1;
        }
#line 12566 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1087:
#line 4517 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        }
#line 12575 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1088:
#line 4522 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        }
#line 12584 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1094:
#line 4536 "p.y" /* yacc.c:1646  */
    { domain->setWaveDirections(0, (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 12590 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1095:
#line 4539 "p.y" /* yacc.c:1646  */
    {
          domain->setKirchhoffLocations((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 12598 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1096:
#line 4542 "p.y" /* yacc.c:1646  */
    {
          domain->setKirchhoffLocations((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 12606 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1099:
#line 4552 "p.y" /* yacc.c:1646  */
    { domain->setFFPDirections((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 12612 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1101:
#line 4559 "p.y" /* yacc.c:1646  */
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[-1].fval));
          StructProp sp;
          sp.kappaHelm = (yyvsp[-1].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzMF);
        }
#line 12625 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1103:
#line 4573 "p.y" /* yacc.c:1646  */
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[-1].fval));
          StructProp sp;
          sp.kappaHelm = (yyvsp[-1].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzSO);
        }
#line 12638 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1104:
#line 4584 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().setProbType(SolverInfo::DisEnrM);
        }
#line 12646 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1105:
#line 4590 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[0].ival); }
#line 12652 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1106:
#line 4592 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[0].ival); }
#line 12658 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1107:
#line 4595 "p.y" /* yacc.c:1646  */
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
#line 12674 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1108:
#line 4607 "p.y" /* yacc.c:1646  */
    {}
#line 12680 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1109:
#line 4609 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::ArcLength; }
#line 12687 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1110:
#line 4612 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::MatNonLinStatic;
          else if(domain->solInfo().probType == SolverInfo::NonLinDynam)
            domain->solInfo().probType = SolverInfo::MatNonLinDynam; }
#line 12696 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1111:
#line 4617 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linearelastic = 1; }
#line 12702 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1112:
#line 4619 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linearelastic = (yyvsp[-1].ival); }
#line 12708 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1113:
#line 4621 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().maxiter = (yyvsp[-1].ival); }
#line 12714 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1114:
#line 4623 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[-1].fval); }
#line 12720 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1115:
#line 4625 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[-2].fval);
          domain->solInfo().getNLInfo().tolInc = (yyvsp[-1].fval); }
#line 12727 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1116:
#line 4628 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[-4].fval);
          domain->solInfo().getNLInfo().tolInc = (yyvsp[-3].fval);
          domain->solInfo().getNLInfo().absTolRes = (yyvsp[-2].fval);
          domain->solInfo().getNLInfo().absTolInc = (yyvsp[-1].fval); }
#line 12736 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1117:
#line 4633 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().dlambda = (yyvsp[-2].fval);
          domain->solInfo().getNLInfo().maxLambda = (yyvsp[-1].fval); }
#line 12743 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1118:
#line 4636 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().dlambda = (yyvsp[-4].fval); 
          domain->solInfo().getNLInfo().maxLambda = (yyvsp[-3].fval);
          domain->solInfo().getNLInfo().extMin = (yyvsp[-2].ival);
          domain->solInfo().getNLInfo().extMax = (yyvsp[-1].ival); }
#line 12752 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1119:
#line 4641 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().fitAlgShell = (yyvsp[-1].ival);
          domain->solInfo().getNLInfo().fitAlgBeam  = (yyvsp[-1].ival); }
#line 12759 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1120:
#line 4644 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().fitAlgShell = (yyvsp[-2].ival);
          domain->solInfo().getNLInfo().fitAlgBeam  = (yyvsp[-1].ival); }
#line 12766 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1121:
#line 4647 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().unsymmetric = true; }
#line 12772 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1122:
#line 4649 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().lfactor = (yyvsp[-1].fval); }
#line 12778 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1123:
#line 4651 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[-1].ival); }
#line 12784 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1124:
#line 4653 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[-2].ival); 
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[-1].ival); }
#line 12791 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1125:
#line 4656 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[-3].ival);
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[-2].ival);
          // note: currently we use either c1 or c2, but never both
          domain->solInfo().getNLInfo().linesearch.c1 = (yyvsp[-1].fval);
          domain->solInfo().getNLInfo().linesearch.c2 = (yyvsp[-1].fval); }
#line 12801 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1126:
#line 4662 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[-4].ival);
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[-3].ival);
          domain->solInfo().getNLInfo().linesearch.c1 = (yyvsp[-2].fval); 
          domain->solInfo().getNLInfo().linesearch.c2 = (yyvsp[-2].fval);
          domain->solInfo().getNLInfo().linesearch.tau = (yyvsp[-1].fval); }
#line 12811 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1127:
#line 4668 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[-5].ival);
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[-4].ival);
          domain->solInfo().getNLInfo().linesearch.c1 = (yyvsp[-3].fval); 
          domain->solInfo().getNLInfo().linesearch.c2 = (yyvsp[-3].fval);
          domain->solInfo().getNLInfo().linesearch.tau = (yyvsp[-2].fval);
          domain->solInfo().getNLInfo().linesearch.verbose = bool((yyvsp[-1].ival)); }
#line 12822 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1128:
#line 4675 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().failsafe = true; }
#line 12828 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1129:
#line 4677 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().failsafe = true;
          domain->solInfo().getNLInfo().failsafe_tol = (yyvsp[-1].fval); }
#line 12835 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1130:
#line 4680 "p.y" /* yacc.c:1646  */
    { domain->solInfo().num_penalty_its = (yyvsp[-3].ival); 
          domain->solInfo().penalty_tol = (yyvsp[-2].fval);
          domain->solInfo().penalty_beta = (yyvsp[-1].fval); }
#line 12843 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1131:
#line 4684 "p.y" /* yacc.c:1646  */
    { domain->solInfo().num_penalty_its = (yyvsp[-5].ival);
          domain->solInfo().penalty_tol = (yyvsp[-4].fval);
          domain->solInfo().penalty_beta = (yyvsp[-3].fval);
          domain->solInfo().reinit_lm = bool((yyvsp[-2].ival));
          domain->solInfo().lm_update_flag = (yyvsp[-1].ival); }
#line 12853 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1132:
#line 4690 "p.y" /* yacc.c:1646  */
    { domain->solInfo().numberOfRomCPUs = (yyvsp[-1].ival); }
#line 12859 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1134:
#line 4695 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().setNewton((yyvsp[-1].ival)); 
          domain->solInfo().solvercntl->fetiInfo.type  = FetiInfo::nonlinear; 
        }
#line 12868 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1135:
#line 4700 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().setNewton((yyvsp[-2].ival));
          domain->solInfo().getNLInfo().stepUpdateK = (yyvsp[-1].ival);
          domain->solInfo().solvercntl->fetiInfo.type  = FetiInfo::nonlinear;
        }
#line 12878 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1136:
#line 4706 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().setNewton((yyvsp[-3].ival));
          domain->solInfo().getNLInfo().stepUpdateK = (yyvsp[-2].ival);
          domain->solInfo().piecewise_contact = bool((yyvsp[-1].ival));
          domain->solInfo().solvercntl->fetiInfo.type  = FetiInfo::nonlinear;
        }
#line 12889 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1137:
#line 4715 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setReOrtho(); }
#line 12895 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1139:
#line 4720 "p.y" /* yacc.c:1646  */
    { geoSource->setControl((yyvsp[-7].strval),(yyvsp[-3].strval),(yyvsp[-1].strval)); domain->solInfo().soltyp = (yyvsp[-5].ival); }
#line 12901 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1140:
#line 4722 "p.y" /* yacc.c:1646  */
    { geoSource->setControl((yyvsp[-9].strval),(yyvsp[-5].strval),(yyvsp[-3].strval),(yyvsp[-1].strval)); domain->solInfo().soltyp = (yyvsp[-7].ival); }
#line 12907 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1142:
#line 4733 "p.y" /* yacc.c:1646  */
    { domain->solInfo().contact_mode = (yyvsp[-1].ival); }
#line 12913 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1143:
#line 4737 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 12919 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1144:
#line 4740 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-7].ival)-1, (yyvsp[-6].ival)-1, (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-1].fval));}
#line 12925 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1145:
#line 4743 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-7].ival)-1, (yyvsp[-6].ival)-1, (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), 0.0, (yyvsp[-1].ival));}
#line 12931 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1146:
#line 4745 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-6].ival)-1, (yyvsp[-5].ival)-1, (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), 0.0, -1, (yyvsp[-1].copt).lagrangeMult, (yyvsp[-1].copt).penalty);}
#line 12937 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1147:
#line 4748 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-9].ival)-1, (yyvsp[-8].ival)-1, (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-1].ival));}
#line 12943 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1148:
#line 4750 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-10].ival)-1, (yyvsp[-9].ival)-1, (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-4].fval), (yyvsp[-2].ival), (yyvsp[-1].copt).lagrangeMult, (yyvsp[-1].copt).penalty);}
#line 12949 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1150:
#line 4755 "p.y" /* yacc.c:1646  */
    { 
           geoSource->addMaterial((yyvsp[-7].ival)-1, 
             new BilinPlasKinHardMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 12958 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1151:
#line 4760 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new BilinPlasKinHardMat((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 12967 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1152:
#line 4765 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new BilinPlasKinHardMat((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 12976 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1153:
#line 4770 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-11].ival)-1, new BilinPlasKinHardMat((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-11].ival)-1, new BilinPlasKinHardMat((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval)) );
           }
         }
#line 12990 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1154:
#line 4780 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new BilinPlasKinHardMat((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new BilinPlasKinHardMat((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-1].fval)) );
           }
         }
#line 13005 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1155:
#line 4791 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-3].fval) > 0 && (yyvsp[-3].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new BilinPlasKinHardMat((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].ival)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new BilinPlasKinHardMat((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-2].fval), (yyvsp[-1].ival)) );
           }
         }
#line 13020 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1156:
#line 4802 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13029 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1157:
#line 4807 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13038 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1158:
#line 4812 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13047 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1159:
#line 4817 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-11].ival)-1, new FiniteStrainPlasKinHardMat((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-11].ival)-1, new FiniteStrainPlasKinHardMat((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval)) );
           }
         }
#line 13061 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1160:
#line 4827 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new FiniteStrainPlasKinHardMat((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new FiniteStrainPlasKinHardMat((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), 
                                    std::numeric_limits<double>::infinity(), (yyvsp[-1].fval)) );
           }
         }
#line 13076 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1161:
#line 4838 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-3].fval) > 0 && (yyvsp[-3].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new FiniteStrainPlasKinHardMat((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), 414) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new FiniteStrainPlasKinHardMat((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-2].fval), (yyvsp[-1].ival)) );
           }
         }
#line 13091 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1162:
#line 4849 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-10].ival)-1, new CrushableFoam((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
          }
#line 13099 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1163:
#line 4853 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-11].ival)-1, new CrushableFoam((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
          }
#line 13107 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1164:
#line 4857 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-12].ival)-1, new CrushableFoam((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
          }
#line 13115 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1165:
#line 4861 "p.y" /* yacc.c:1646  */
    {
            if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
              geoSource->addMaterial((yyvsp[-13].ival)-1, new CrushableFoam((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
              domain->solInfo().elementDeletion = true;
            }
            else {
              geoSource->addMaterial((yyvsp[-13].ival)-1, new CrushableFoam((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), 
                                     (yyvsp[-2].fval), std::numeric_limits<double>::infinity()) );
            }
          }
#line 13130 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1166:
#line 4872 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13139 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1167:
#line 4877 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13148 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1168:
#line 4882 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13157 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1169:
#line 4887 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-11].ival)-1, new LogStrainPlasKinHardMat((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-11].ival)-1, new LogStrainPlasKinHardMat((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval)) );
           }
         }
#line 13171 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1170:
#line 4897 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new LogStrainPlasKinHardMat((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new LogStrainPlasKinHardMat((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), 
                                    std::numeric_limits<double>::infinity(), (yyvsp[-1].fval)) );
           }
         }
#line 13186 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1171:
#line 4908 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-3].fval) > 0 && (yyvsp[-3].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new LogStrainPlasKinHardMat((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].ival)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new LogStrainPlasKinHardMat((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-2].fval), (yyvsp[-1].ival)) );
           }
         }
#line 13201 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1172:
#line 4919 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new ElaLinIsoMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13210 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1173:
#line 4924 "p.y" /* yacc.c:1646  */
    { 
           geoSource->addMaterial((yyvsp[-5].ival)-1, 
             new ElaLinIsoMat((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
	 }
#line 13219 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1174:
#line 4929 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-3].ival)-1,
             new ElaLinIsoMat((yyvsp[-1].fval)));
         }
#line 13228 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1175:
#line 4934 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new BrittleFractureTB<ElaLinIsoMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13238 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1176:
#line 4940 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new BrittleFractureTB<ElaLinIsoMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13248 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1177:
#line 4946 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new BrittleFractureTB<ElaLinIsoMat>((yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13258 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1178:
#line 4952 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new StVenantKirchhoffMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13267 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1179:
#line 4957 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-5].ival)-1,
             new StVenantKirchhoffMat((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13276 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1180:
#line 4962 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-3].ival)-1,
             new StVenantKirchhoffMat((yyvsp[-1].fval)));
         }
#line 13285 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1181:
#line 4967 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new BrittleFractureTB<StVenantKirchhoffMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13295 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1182:
#line 4973 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new BrittleFractureTB<StVenantKirchhoffMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13305 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1183:
#line 4979 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new BrittleFractureTB<StVenantKirchhoffMat>((yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13315 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1184:
#line 4985 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new HenckyMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13324 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1185:
#line 4990 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-5].ival)-1,
             new HenckyMat((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13333 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1186:
#line 4995 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-3].ival)-1,
             new HenckyMat((yyvsp[-1].fval)));
         }
#line 13342 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1187:
#line 5000 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new BrittleFractureTB<HenckyMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13352 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1188:
#line 5006 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new BrittleFractureTB<HenckyMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13362 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1189:
#line 5012 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new BrittleFractureTB<HenckyMat>((yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13372 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1190:
#line 5018 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new ElaLinIsoMat2D((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0, 0));
         }
#line 13381 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1191:
#line 5023 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new ElaLinIsoMat2D((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13390 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1192:
#line 5028 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new StVenantKirchhoffMat2D((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0, 0));
         }
#line 13399 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1193:
#line 5033 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new StVenantKirchhoffMat2D((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13408 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1194:
#line 5038 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new FabricMap((yyvsp[-5].fval), (yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-2].ival), (yyvsp[-1].fval), 0, 0, FabricMap::GREEN_LAGRANGE));
         }
#line 13417 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1195:
#line 5043 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new FabricMap((yyvsp[-7].fval), (yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), FabricMap::GREEN_LAGRANGE));
         }
#line 13426 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1196:
#line 5048 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new FabricMap((yyvsp[-5].fval), (yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-2].ival), (yyvsp[-1].fval), 0, 0, FabricMap::INFINTESIMAL));
         }
#line 13435 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1197:
#line 5053 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new FabricMap((yyvsp[-7].fval), (yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), FabricMap::INFINTESIMAL));
         }
#line 13444 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1198:
#line 5058 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new FabricMat((yyvsp[-7].fval), (yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0, 0, FabricMat::GREEN_LAGRANGE));
         }
#line 13453 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1199:
#line 5063 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new FabricMat((yyvsp[-9].fval), (yyvsp[-8].ival), (yyvsp[-7].ival), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), FabricMat::GREEN_LAGRANGE));
         }
#line 13462 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1200:
#line 5068 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new FabricMat((yyvsp[-7].fval), (yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0, 0, FabricMat::INFINTESIMAL));
         }
#line 13471 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1201:
#line 5073 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new FabricMat((yyvsp[-9].fval), (yyvsp[-8].ival), (yyvsp[-7].ival), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), FabricMat::INFINTESIMAL));
         }
#line 13480 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1202:
#line 5078 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new PlaneStressMat<ElaLinIsoMat>((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13489 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1203:
#line 5083 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new PlaneStressMat<ElaLinIsoMat>((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13498 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1204:
#line 5088 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new PlaneStressMat<BrittleFractureTB<ElaLinIsoMat> >((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13508 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1205:
#line 5094 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-12].ival)-1,
             new PlaneStressMat<BrittleFractureTB<ElaLinIsoMat> >((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13518 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1206:
#line 5100 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new PlaneStressMat<StVenantKirchhoffMat>((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13527 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1207:
#line 5105 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new PlaneStressMat<StVenantKirchhoffMat>((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13536 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1208:
#line 5110 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new PlaneStressMat<BrittleFractureTB<StVenantKirchhoffMat> >((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13546 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1209:
#line 5116 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-12].ival)-1,
             new PlaneStressMat<BrittleFractureTB<StVenantKirchhoffMat> >((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13556 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1210:
#line 5122 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new PlaneStressMat<NeoHookeanMat>((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13565 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1211:
#line 5127 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new PlaneStressMat<BrittleFractureTB<NeoHookeanMat> >((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13575 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1212:
#line 5133 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new PlaneStressMat<MooneyRivlinMat>((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13584 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1213:
#line 5138 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new PlaneStressMat<BrittleFractureTB<MooneyRivlinMat> >((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13594 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1214:
#line 5144 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13603 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1215:
#line 5149 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13612 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1216:
#line 5154 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13621 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1217:
#line 5159 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-1].fval), (yyvsp[-2].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval)) );
           }
         }
#line 13635 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1218:
#line 5169 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-3].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-1].fval), (yyvsp[-3].fval)) );
           }
         }
#line 13650 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1219:
#line 5180 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-3].fval) > 0 && (yyvsp[-3].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-14].ival)-1, new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].ival), (yyvsp[-4].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-14].ival)-1, new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-2].fval), (yyvsp[-1].ival), (yyvsp[-4].fval)) );
           }
         }
#line 13665 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1220:
#line 5191 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13674 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1221:
#line 5196 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13683 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1222:
#line 5201 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13692 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1223:
#line 5206 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-1].fval), (yyvsp[-2].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval)) );
           }
         }
#line 13706 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1224:
#line 5216 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-3].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-1].fval), (yyvsp[-3].fval)) );
           }
         }
#line 13721 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1225:
#line 5227 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-3].fval) > 0 && (yyvsp[-3].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-14].ival)-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].ival), (yyvsp[-4].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-14].ival)-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-2].fval), (yyvsp[-1].ival), (yyvsp[-4].fval)) );
           }
         }
#line 13736 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1226:
#line 5238 "p.y" /* yacc.c:1646  */
    {
            double params[3] = { (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-5].ival)-1,
              new MaterialWrapper<IsotropicLinearElastic>(params));
          }
#line 13746 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1227:
#line 5244 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-11].ival)-1,
              new PronyViscoElastic<ElaLinIsoMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 13759 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1228:
#line 5253 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-13].ival)-1,
              new PronyViscoElastic<ElaLinIsoMat>((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 13772 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1229:
#line 5262 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-15].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> >((yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 13786 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1230:
#line 5272 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-17].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> >((yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 13800 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1231:
#line 5282 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-12].ival)-1,
             new PronyViscoElastic<ElaLinIsoMat2D>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), 0, 0, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13813 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1232:
#line 5291 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-14].ival)-1,
             new PronyViscoElastic<ElaLinIsoMat2D>((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13826 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1233:
#line 5300 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-12].ival)-1,
             new PronyViscoElastic<StVenantKirchhoffMat2D>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), 0, 0, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13839 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1234:
#line 5309 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-14].ival)-1,
             new PronyViscoElastic<StVenantKirchhoffMat2D>((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13852 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1235:
#line 5318 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-13].ival)-1,
             new PronyViscoElastic<FabricMap>((yyvsp[-11].fval), (yyvsp[-10].ival), (yyvsp[-9].ival), (yyvsp[-8].ival), (yyvsp[-7].fval), 0, 0, FabricMap::GREEN_LAGRANGE, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13865 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1236:
#line 5327 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-15].ival)-1,
             new PronyViscoElastic<FabricMap>((yyvsp[-13].fval), (yyvsp[-12].ival), (yyvsp[-11].ival), (yyvsp[-10].ival), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), FabricMap::GREEN_LAGRANGE, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13878 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1237:
#line 5336 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-13].ival)-1,
             new PronyViscoElastic<FabricMap>((yyvsp[-11].fval), (yyvsp[-10].ival), (yyvsp[-9].ival), (yyvsp[-8].ival), (yyvsp[-7].fval), 0, 0, FabricMap::INFINTESIMAL, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13891 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1238:
#line 5345 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-15].ival)-1,
             new PronyViscoElastic<FabricMap>((yyvsp[-13].fval), (yyvsp[-12].ival), (yyvsp[-11].ival), (yyvsp[-10].ival), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), FabricMap::INFINTESIMAL, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13904 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1239:
#line 5354 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-4].fval);
           double gtwo   = (yyvsp[-6].fval);
           double gone   = (yyvsp[-8].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-15].ival)-1,
             new PronyViscoElastic<FabricMat>((yyvsp[-13].fval), (yyvsp[-12].ival), (yyvsp[-11].ival), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), 0, 0, FabricMat::GREEN_LAGRANGE, ginf, (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval)));
         }
#line 13917 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1240:
#line 5363 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-4].fval);
           double gtwo   = (yyvsp[-6].fval);
           double gone   = (yyvsp[-8].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-17].ival)-1,
             new PronyViscoElastic<FabricMat>((yyvsp[-15].fval), (yyvsp[-14].ival), (yyvsp[-13].ival), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), FabricMat::GREEN_LAGRANGE, ginf, (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval)));
         }
#line 13930 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1241:
#line 5372 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-4].fval);
           double gtwo   = (yyvsp[-6].fval);
           double gone   = (yyvsp[-8].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-15].ival)-1,
             new PronyViscoElastic<FabricMat>((yyvsp[-13].fval), (yyvsp[-12].ival), (yyvsp[-11].ival), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), 0, 0, FabricMat::INFINTESIMAL, ginf, (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval)));
         }
#line 13943 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1242:
#line 5381 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-4].fval);
           double gtwo   = (yyvsp[-6].fval);
           double gone   = (yyvsp[-8].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-17].ival)-1,
             new PronyViscoElastic<FabricMat>((yyvsp[-15].fval), (yyvsp[-14].ival), (yyvsp[-13].ival), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), FabricMat::INFINTESIMAL, ginf, (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval)));
         }
#line 13956 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1243:
#line 5390 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-12].ival)-1,
              new PlaneStressMat<PronyViscoElastic<ElaLinIsoMat> >((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 13969 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1244:
#line 5399 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-14].ival)-1,
              new PlaneStressMat<PronyViscoElastic<ElaLinIsoMat> >((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 13982 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1245:
#line 5408 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval); 
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-16].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> > >((yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 13996 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1246:
#line 5418 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-18].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> > >((yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14010 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1247:
#line 5428 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-11].ival)-1,
              new PronyViscoElastic<StVenantKirchhoffMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14023 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1248:
#line 5437 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-13].ival)-1,
              new PronyViscoElastic<StVenantKirchhoffMat>((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14036 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1249:
#line 5446 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-15].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> >((yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14050 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1250:
#line 5456 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-17].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> >((yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14064 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1251:
#line 5466 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-12].ival)-1,
              new PlaneStressMat<PronyViscoElastic<StVenantKirchhoffMat> >((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14077 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1252:
#line 5475 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-14].ival)-1,
              new PlaneStressMat<PronyViscoElastic<StVenantKirchhoffMat> >((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14090 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1253:
#line 5484 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-16].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> > >((yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14104 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1254:
#line 5494 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-18].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> > >((yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14118 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1255:
#line 5504 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-5].ival)-1,
              new NeoHookeanMat((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14127 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1256:
#line 5509 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-9].ival)-1,
              new BrittleFractureTB<NeoHookeanMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14137 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1257:
#line 5515 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-11].ival)-1,
              new PronyViscoElastic<NeoHookeanMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14150 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1258:
#line 5524 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-15].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<NeoHookeanMat> >((yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14164 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1259:
#line 5534 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-12].ival)-1,
              new PlaneStressMat<PronyViscoElastic<NeoHookeanMat> >((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14177 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1260:
#line 5543 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-16].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<NeoHookeanMat> > >((yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14191 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1261:
#line 5553 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-6].ival)-1,
              new MooneyRivlinMat((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14200 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1262:
#line 5558 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-10].ival)-1,
              new BrittleFractureTB<MooneyRivlinMat>((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14210 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1263:
#line 5564 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-12].ival)-1,
              new PronyViscoElastic<MooneyRivlinMat>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14223 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1264:
#line 5573 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-16].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<MooneyRivlinMat> >((yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14237 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1265:
#line 5583 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-13].ival)-1,
              new PlaneStressMat<PronyViscoElastic<MooneyRivlinMat> >((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14250 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1266:
#line 5592 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-17].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<MooneyRivlinMat> > >((yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14264 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1267:
#line 5602 "p.y" /* yacc.c:1646  */
    {
            double mu[1] = { (yyvsp[-3].fval) };
            double alpha[1] = { (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-6].ival)-1, new OgdenMat((yyvsp[-4].fval), mu, alpha, K));
          }
#line 14275 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1268:
#line 5609 "p.y" /* yacc.c:1646  */
    {
            double mu[2] = { (yyvsp[-5].fval), (yyvsp[-4].fval) };
            double alpha[2] = { (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-8].ival)-1, new OgdenMat((yyvsp[-6].fval), mu, alpha, K));
          }
#line 14286 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1269:
#line 5616 "p.y" /* yacc.c:1646  */
    {
            double mu[2] = { (yyvsp[-6].fval), (yyvsp[-5].fval) }; 
            double alpha[2] = { (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-9].ival)-1, new OgdenMat((yyvsp[-7].fval), mu, alpha, K));
          }
#line 14297 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1270:
#line 5623 "p.y" /* yacc.c:1646  */
    {
            double mu[3] = { (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval) };
            double alpha[3] = { (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-10].ival)-1, new OgdenMat((yyvsp[-8].fval), mu, alpha, K));
          }
#line 14308 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1271:
#line 5630 "p.y" /* yacc.c:1646  */
    {
            double mu[3] = { (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval) };
            double alpha[3] = { (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-11].ival)-1, new OgdenMat((yyvsp[-9].fval), mu, alpha, K));
          }
#line 14319 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1272:
#line 5637 "p.y" /* yacc.c:1646  */
    {
            double mu[4] = { (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval) };
            double alpha[4] = { (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-12].ival)-1, new OgdenMat((yyvsp[-10].fval), mu, alpha, K));
          }
#line 14330 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1273:
#line 5644 "p.y" /* yacc.c:1646  */
    {
            double mu[4] = { (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval) };
            double alpha[4] = { (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-13].ival)-1, new OgdenMat((yyvsp[-11].fval), mu, alpha, K));
          }
#line 14341 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1274:
#line 5651 "p.y" /* yacc.c:1646  */
    {
            double mu[5] = { (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval) };
            double alpha[5] = { (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-14].ival)-1, new OgdenMat((yyvsp[-12].fval), mu, alpha, K));
          }
#line 14352 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1275:
#line 5658 "p.y" /* yacc.c:1646  */
    {
            double mu[5] = { (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval) };
            double alpha[5] = { (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-15].ival)-1, new OgdenMat((yyvsp[-13].fval), mu, alpha, K));
          }
#line 14363 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1276:
#line 5665 "p.y" /* yacc.c:1646  */
    {
            double mu[6] = { (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval) };
            double alpha[6] = { (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-16].ival)-1, new OgdenMat((yyvsp[-14].fval), mu, alpha, K));
          }
#line 14374 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1277:
#line 5672 "p.y" /* yacc.c:1646  */
    {
            double mu[6] = { (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval) };
            double alpha[6] = { (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-17].ival)-1, new OgdenMat((yyvsp[-15].fval), mu, alpha, K));
          }
#line 14385 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1278:
#line 5679 "p.y" /* yacc.c:1646  */
    {
            double mu[7] = { (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval) };
            double alpha[7] = { (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-18].ival)-1, new OgdenMat((yyvsp[-16].fval), mu, alpha, K));
          }
#line 14396 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1279:
#line 5686 "p.y" /* yacc.c:1646  */
    {
            double mu[7] = { (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval) };
            double alpha[7] = { (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-19].ival)-1, new OgdenMat((yyvsp[-17].fval), mu, alpha, K));
          }
#line 14407 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1280:
#line 5693 "p.y" /* yacc.c:1646  */
    {
            double mu[8] = { (yyvsp[-17].fval), (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval) };
            double alpha[8] = { (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-20].ival)-1, new OgdenMat((yyvsp[-18].fval), mu, alpha, K));
          }
#line 14418 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1281:
#line 5700 "p.y" /* yacc.c:1646  */
    {
            double mu[8] = { (yyvsp[-18].fval), (yyvsp[-17].fval), (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval) };
            double alpha[8] = { (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-21].ival)-1, new OgdenMat((yyvsp[-19].fval), mu, alpha, K));
          }
#line 14429 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1282:
#line 5707 "p.y" /* yacc.c:1646  */
    {
            double mu[9] = { (yyvsp[-19].fval), (yyvsp[-18].fval), (yyvsp[-17].fval), (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval) };
            double alpha[9] = { (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-22].ival)-1, new OgdenMat((yyvsp[-20].fval), mu, alpha, K));
          }
#line 14440 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1283:
#line 5714 "p.y" /* yacc.c:1646  */
    {
            double mu[9] = { (yyvsp[-20].fval), (yyvsp[-19].fval), (yyvsp[-18].fval), (yyvsp[-17].fval), (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval) };
            double alpha[9] = { (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-23].ival)-1, new OgdenMat((yyvsp[-21].fval), mu, alpha, K));
          }
#line 14451 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1284:
#line 5721 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-5].ival)-1,
             new SimoElasticMat((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14460 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1285:
#line 5726 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new SimoPlasticMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 14469 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1286:
#line 5731 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new SimoPlasticMat((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 14478 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1287:
#line 5736 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 1.0e-6, std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial((yyvsp[-8].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
          }
#line 14488 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1288:
#line 5742 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial((yyvsp[-9].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
          }
#line 14498 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1289:
#line 5748 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0. };
            geoSource->addMaterial((yyvsp[-10].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
            if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
#line 14511 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1290:
#line 5757 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), double((yyvsp[-1].ival)) };
            geoSource->addMaterial((yyvsp[-11].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
            if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
#line 14524 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1291:
#line 5766 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 1.0e-6, std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial((yyvsp[-8].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          }
#line 14534 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1292:
#line 5772 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial((yyvsp[-9].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          }
#line 14544 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1293:
#line 5778 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0. };
            geoSource->addMaterial((yyvsp[-10].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
            if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
#line 14557 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1294:
#line 5787 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), double((yyvsp[-1].ival)) };
            geoSource->addMaterial((yyvsp[-11].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
            if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
#line 14570 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1295:
#line 5796 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-5].ival)-1,
             new ExpMat((yyvsp[-4].ival), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14579 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1296:
#line 5801 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new ExpMat((yyvsp[-5].ival), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14588 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1297:
#line 5806 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new ExpMat((yyvsp[-6].ival), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14597 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1298:
#line 5811 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new ExpMat((yyvsp[-7].ival), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14606 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1299:
#line 5816 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new ExpMat((yyvsp[-8].ival), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14615 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1300:
#line 5821 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new ExpMat((yyvsp[-9].ival), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             domain->solInfo().elementDeletion = true;
           }
         }
#line 14627 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1301:
#line 5829 "p.y" /* yacc.c:1646  */
    {
           ExpMat *mat = new ExpMat((yyvsp[-10].ival), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval));
           mat->yssrtid = (yyvsp[-1].ival);
           geoSource->addMaterial((yyvsp[-11].ival)-1, mat);
           if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
             domain->solInfo().elementDeletion = true;
           }
         }
#line 14640 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1302:
#line 5839 "p.y" /* yacc.c:1646  */
    {
           ExpMat *mat = new ExpMat((yyvsp[-19].ival), (yyvsp[-18].fval), (yyvsp[-17].fval), (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval));
           mat->yssrtid = (yyvsp[-10].ival);
           mat->optcor0 = (yyvsp[-9].ival);
           mat->optcor1 = (yyvsp[-8].ival);
           mat->optprj = (yyvsp[-7].ival);
           mat->opthgc = (yyvsp[-6].ival);
           mat->prmhgc[0] = (yyvsp[-5].fval);
           mat->prmhgc[1] = (yyvsp[-4].fval);
           mat->prmhgc[2] = (yyvsp[-3].fval);
           mat->ngqpt2 = (yyvsp[-2].ival);
           mat->ematpro[18] = (yyvsp[-1].fval);
           geoSource->addMaterial((yyvsp[-20].ival)-1, mat);
           if((yyvsp[-11].fval) > 0 && (yyvsp[-11].fval) < std::numeric_limits<double>::infinity()) {
             domain->solInfo().elementDeletion = true;
           }
         }
#line 14662 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1303:
#line 5857 "p.y" /* yacc.c:1646  */
    {
	   geoSource->loadMaterial((yyvsp[-2].strval), (yyvsp[-1].strval));
	 }
#line 14670 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1304:
#line 5861 "p.y" /* yacc.c:1646  */
    {
	   geoSource->addMaterial((yyvsp[-3].ival)-1, (yyvsp[-2].strval), (yyvsp[-1].dlist));
	 }
#line 14678 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1306:
#line 5868 "p.y" /* yacc.c:1646  */
    { geoSource->setMatUsage((yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 14684 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1307:
#line 5870 "p.y" /* yacc.c:1646  */
    {
            for(int i = (yyvsp[-3].ival)-1; i < (yyvsp[-2].ival); ++i)
	      geoSource->setMatUsage(i, (yyvsp[-1].ival)-1);
	  }
#line 14693 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1308:
#line 5876 "p.y" /* yacc.c:1646  */
    { (yyval.dlist).nval = 0; }
#line 14699 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1309:
#line 5878 "p.y" /* yacc.c:1646  */
    { 
          if((yyvsp[-1].dlist).nval == 64) {
             fprintf(stderr, "You'd better invent another material model!\n");
	     exit(-1);
          }
          (yyval.dlist) = (yyvsp[-1].dlist);
          (yyval.dlist).v[(yyval.dlist).nval++] = (yyvsp[0].fval);
 	}
#line 14712 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1310:
#line 5888 "p.y" /* yacc.c:1646  */
    { (yyval.slist).nval = 0; }
#line 14718 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1311:
#line 5890 "p.y" /* yacc.c:1646  */
    { 
          if((yyvsp[-1].slist).nval == 32) {
             fprintf(stderr, "Too many files!\n");
	     exit(-1);
          }
          (yyval.slist) = (yyvsp[-1].slist);
          (yyval.slist).v[(yyval.slist).nval++] = (yyvsp[0].strval);
 	}
#line 14731 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1312:
#line 5901 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setRenum((yyvsp[-1].ival));
          domain->solInfo().setSparseRenum((yyvsp[-1].ival)); 
          domain->solInfo().setSpoolesRenum((yyvsp[-1].ival)); }
#line 14739 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1313:
#line 5905 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setRenum((yyvsp[-3].ival));
          domain->solInfo().setSparseRenum((yyvsp[-1].ival)); }
#line 14746 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1314:
#line 5908 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setRenum((yyvsp[-5].ival));
          domain->solInfo().setSparseRenum((yyvsp[-3].ival)); 
          domain->solInfo().setSpoolesRenum((yyvsp[-1].ival)); }
#line 14754 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1315:
#line 5915 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true; 
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().svdPodRom = true;}
#line 14762 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1317:
#line 5923 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().snapfiPodRom.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14768 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1318:
#line 5934 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().velocPodRomFile.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14774 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1319:
#line 5936 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().accelPodRomFile.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14780 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1320:
#line 5938 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().dsvPodRomFile.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14786 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1321:
#line 5940 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().muvPodRomFile.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14792 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1322:
#line 5942 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxSizePodRom = (yyvsp[0].ival); }
#line 14798 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1323:
#line 5944 "p.y" /* yacc.c:1646  */
    { domain->solInfo().normalize = (yyvsp[0].ival); }
#line 14804 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1324:
#line 5946 "p.y" /* yacc.c:1646  */
    { domain->solInfo().normalize = (yyvsp[0].ival); }
#line 14810 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1325:
#line 5948 "p.y" /* yacc.c:1646  */
    { domain->solInfo().normalize = (yyvsp[0].ival); }
#line 14816 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1326:
#line 5950 "p.y" /* yacc.c:1646  */
    { domain->solInfo().subtractRefPodRom = true;
    domain->solInfo().readInLocalBasesCent.push_back(std::string((yyvsp[0].strval))); }
#line 14823 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1327:
#line 5953 "p.y" /* yacc.c:1646  */
    { domain->solInfo().flagss = (yyvsp[0].ival); }
#line 14829 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1328:
#line 5955 "p.y" /* yacc.c:1646  */
    { domain->solInfo().flagss = (yyvsp[-1].ival); domain->solInfo().flagrs = (yyvsp[0].ival); }
#line 14835 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1329:
#line 5957 "p.y" /* yacc.c:1646  */
    { domain->solInfo().skipPodRom = (yyvsp[0].ival); }
#line 14841 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1330:
#line 5959 "p.y" /* yacc.c:1646  */
    { domain->solInfo().skipPodRom = (yyvsp[-1].ival);
    domain->solInfo().skipOffSet = (yyvsp[0].ival); }
#line 14848 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1331:
#line 5962 "p.y" /* yacc.c:1646  */
    { domain->solInfo().robcSolve = bool((yyvsp[0].ival)); }
#line 14854 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1332:
#line 5964 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().robfi.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14860 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1333:
#line 5966 "p.y" /* yacc.c:1646  */
    { domain->solInfo().svdBlockSize = (yyvsp[0].ival); }
#line 14866 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1334:
#line 5968 "p.y" /* yacc.c:1646  */
    { domain->solInfo().romEnergy = (yyvsp[0].fval); }
#line 14872 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1335:
#line 5981 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf            = 3;
    domain->solInfo().nmfMaxIter         = (yyvsp[-3].ival);
    domain->solInfo().nmfTol             = (yyvsp[-2].fval);
    domain->solInfo().nmfPqnNumInnerIter = (yyvsp[-1].ival);
    domain->solInfo().nmfPqnAlpha        = (yyvsp[0].fval); }
#line 14882 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1336:
#line 5987 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf = 1;
    domain->solInfo().nmfNumROBDim = (yyvsp[-4].ival);
    domain->solInfo().nmfDelROBDim = (yyvsp[-3].ival);
    domain->solInfo().nmfRandInit  = (yyvsp[-2].ival);
    domain->solInfo().nmfMaxIter   = (yyvsp[-1].ival);
    domain->solInfo().nmfTol = (yyvsp[0].fval); }
#line 14893 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1337:
#line 5994 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf    = 1;
    domain->solInfo().nmfMaxIter = (yyvsp[-1].ival);
    domain->solInfo().nmfTol     = (yyvsp[0].fval); }
#line 14901 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1338:
#line 5998 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf    = 4;   
    domain->solInfo().nmfMaxIter = (yyvsp[-4].ival); 
    domain->solInfo().nmfTol     = (yyvsp[-3].fval); 
    domain->solInfo().nmfcAlpha  = (yyvsp[-2].fval);
    domain->solInfo().nmfcBeta   = (yyvsp[-1].fval);
    domain->solInfo().nmfcGamma  = (yyvsp[0].fval);}
#line 14912 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1339:
#line 6005 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf    = 4;
    domain->solInfo().nmfMaxIter = (yyvsp[-1].ival);
    domain->solInfo().nmfTol     = (yyvsp[0].fval); }
#line 14920 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1340:
#line 6009 "p.y" /* yacc.c:1646  */
    { domain->solInfo().nmfNumSub = (yyvsp[0].ival); }
#line 14926 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1341:
#line 6011 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf = 2; }
#line 14932 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1342:
#line 6013 "p.y" /* yacc.c:1646  */
    { domain->solInfo().clustering = (yyvsp[0].ival); }
#line 14938 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1343:
#line 6015 "p.y" /* yacc.c:1646  */
    { domain->solInfo().clustering = (yyvsp[-1].ival); 
    domain->solInfo().clusterSubspaceAngle = true; }
#line 14945 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1344:
#line 6018 "p.y" /* yacc.c:1646  */
    { domain->solInfo().solverTypeCluster = (yyvsp[0].ival); }
#line 14951 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1345:
#line 6020 "p.y" /* yacc.c:1646  */
    { domain->solInfo().solverTypeCluster = (yyvsp[-1].ival);
    domain->solInfo().tolPodRom = (yyvsp[0].fval);}
#line 14958 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1346:
#line 6023 "p.y" /* yacc.c:1646  */
    { domain->solInfo().solverTypeCluster = (yyvsp[-2].ival);
    domain->solInfo().tolPodRom = (yyvsp[-1].fval);
    domain->solInfo().solverTypeSpnnls = (yyvsp[0].ival); }
#line 14966 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1347:
#line 6027 "p.y" /* yacc.c:1646  */
    { domain->solInfo().hotstartSample = bool((yyvsp[0].ival)); }
#line 14972 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1348:
#line 6029 "p.y" /* yacc.c:1646  */
    { domain->solInfo().rowClustering = (yyvsp[0].ival); }
#line 14978 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1349:
#line 6031 "p.y" /* yacc.c:1646  */
    { domain->solInfo().rowClustering = (yyvsp[-1].ival);
    domain->solInfo().clusterSubspaceAngle = true; }
#line 14985 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1351:
#line 6038 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true;
     domain->solInfo().setProbType(SolverInfo::PodRomOffline);
     domain->solInfo().DEIMBasisPod = true; }
#line 14993 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1353:
#line 6046 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true;
     domain->solInfo().setProbType(SolverInfo::PodRomOffline);
     domain->solInfo().UDEIMBasisPod = true; }
#line 15001 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1355:
#line 6054 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true; 
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().samplingPodRom = true; }
#line 15009 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1360:
#line 6065 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true;
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().snapProjPodRom = true; }
#line 15017 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1362:
#line 6073 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInROBorModes.push_back((yyvsp[0].strval)); }
#line 15023 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1363:
#line 6075 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInROBorModes.push_back((yyvsp[-1].strval));
    domain->solInfo().localBasisSize.push_back((yyvsp[0].ival)); }
#line 15030 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1364:
#line 6078 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInDualROB.push_back((yyvsp[-1].strval));
    domain->solInfo().localDualBasisSize.push_back((yyvsp[0].ival)); }
#line 15037 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1365:
#line 6081 "p.y" /* yacc.c:1646  */
    { domain->solInfo().statePodRomFile.push_back((yyvsp[0].strval)); }
#line 15043 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1366:
#line 6083 "p.y" /* yacc.c:1646  */
    { domain->solInfo().statePodRomFile.push_back((yyvsp[-1].strval));
    domain->solInfo().velocPodRomFile.push_back((yyvsp[0].strval)); }
#line 15050 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1367:
#line 6086 "p.y" /* yacc.c:1646  */
    { domain->solInfo().statePodRomFile.push_back((yyvsp[-2].strval));
    domain->solInfo().velocPodRomFile.push_back((yyvsp[-1].strval));
    domain->solInfo().accelPodRomFile.push_back((yyvsp[0].strval)); }
#line 15058 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1368:
#line 6090 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tolPodRom = (yyvsp[0].fval); }
#line 15064 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1369:
#line 6092 "p.y" /* yacc.c:1646  */
    { domain->solInfo().skipPodRom = (yyvsp[0].ival); }
#line 15070 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1370:
#line 6094 "p.y" /* yacc.c:1646  */
    { domain->solInfo().randomSampleSize = (yyvsp[0].ival); 
    domain->solInfo().randomVecSampling = true; }
#line 15077 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1371:
#line 6097 "p.y" /* yacc.c:1646  */
    { domain->solInfo().skipOffSet = (yyvsp[0].ival); }
#line 15083 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1372:
#line 6099 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxSizePodRom = (yyvsp[0].ival); }
#line 15089 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1373:
#line 6101 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxSizePodRom = (yyvsp[-1].ival); 
    domain->solInfo().forcePodSize = (yyvsp[0].ival);}
#line 15096 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1374:
#line 6104 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxSizePodRom = (yyvsp[-2].ival); 
    domain->solInfo().forcePodSize = (yyvsp[-1].ival);
    domain->solInfo().maxDeimBasisSize = (yyvsp[0].ival); }
#line 15104 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1375:
#line 6108 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useMassNormalizedBasis = bool((yyvsp[0].ival)); }
#line 15110 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1376:
#line 6110 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useConstantMassForces = bool((yyvsp[0].ival)); }
#line 15116 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1377:
#line 6112 "p.y" /* yacc.c:1646  */
    { domain->solInfo().stackedElementSampling = bool((yyvsp[0].ival)); }
#line 15122 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1378:
#line 6114 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useMassOrthogonalProjection = bool((yyvsp[0].ival)); }
#line 15128 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1379:
#line 6116 "p.y" /* yacc.c:1646  */
    { domain->solInfo().reduceFollower = bool((yyvsp[0].ival)); }
#line 15134 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1380:
#line 6118 "p.y" /* yacc.c:1646  */
    { domain->solInfo().reduceFollower = bool((yyvsp[0].ival)); }
#line 15140 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1381:
#line 6120 "p.y" /* yacc.c:1646  */
    { domain->solInfo().PODerrornorm.push_back((yyvsp[0].strval)); }
#line 15146 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1382:
#line 6122 "p.y" /* yacc.c:1646  */
    { domain->solInfo().PODerrornorm.push_back((yyvsp[-1].strval));
    domain->solInfo().PODerrornorm.push_back((yyvsp[0].strval)); }
#line 15153 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1383:
#line 6125 "p.y" /* yacc.c:1646  */
    { domain->solInfo().PODerrornorm.push_back((yyvsp[-2].strval));
    domain->solInfo().PODerrornorm.push_back((yyvsp[-1].strval));
    domain->solInfo().PODerrornorm.push_back((yyvsp[0].strval)); }
#line 15161 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1384:
#line 6129 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useScalingSpnnls = bool((yyvsp[0].ival)); }
#line 15167 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1385:
#line 6131 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useCenterSpnnls = bool((yyvsp[0].ival)); }
#line 15173 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1386:
#line 6133 "p.y" /* yacc.c:1646  */
    { domain->solInfo().projectSolution = bool((yyvsp[0].ival)); }
#line 15179 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1387:
#line 6135 "p.y" /* yacc.c:1646  */
    { domain->solInfo().positiveElements = bool((yyvsp[0].ival)); }
#line 15185 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1388:
#line 6137 "p.y" /* yacc.c:1646  */
    { domain->solInfo().hotstartSample = bool((yyvsp[0].ival)); }
#line 15191 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1389:
#line 6139 "p.y" /* yacc.c:1646  */
    { domain->solInfo().solverTypeSpnnls = (yyvsp[0].ival); }
#line 15197 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1390:
#line 6141 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxSizeSpnnls = (yyvsp[0].fval); }
#line 15203 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1391:
#line 6143 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxElemSpnnls = (yyvsp[0].ival); }
#line 15209 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1392:
#line 6145 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxIterSpnnls = (yyvsp[0].fval); }
#line 15215 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1393:
#line 6147 "p.y" /* yacc.c:1646  */
    { domain->solInfo().forcePodRomFile = (yyvsp[0].strval); }
#line 15221 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1394:
#line 6149 "p.y" /* yacc.c:1646  */
    { domain->solInfo().forcePodRomFile = (yyvsp[-1].strval);
    domain->solInfo().forcePodSize = (yyvsp[0].ival); }
#line 15228 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1395:
#line 6152 "p.y" /* yacc.c:1646  */
    { domain->solInfo().forcePodRomFile = (yyvsp[-2].strval); 
    domain->solInfo().forcePodSize = (yyvsp[-1].ival); 
    domain->solInfo().maxDeimBasisSize = (yyvsp[0].ival); }
#line 15236 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1396:
#line 6156 "p.y" /* yacc.c:1646  */
    { domain->solInfo().constraintPodRomFile = (yyvsp[0].strval); 
    domain->solInfo().ConstraintBasisPod = true;}
#line 15243 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1397:
#line 6159 "p.y" /* yacc.c:1646  */
    { domain->solInfo().constraintPodRomFile = (yyvsp[-1].strval);
    domain->solInfo().constraintPodSize = (yyvsp[0].ival); 
    domain->solInfo().ConstraintBasisPod = true; }
#line 15251 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1398:
#line 6163 "p.y" /* yacc.c:1646  */
    { domain->solInfo().constraintPodRomFile = (yyvsp[-2].strval);
    domain->solInfo().constraintPodSize = (yyvsp[-1].ival);
    domain->solInfo().maxDeimBasisSize = (yyvsp[0].ival); 
    domain->solInfo().ConstraintBasisPod = true; }
#line 15260 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1399:
#line 6168 "p.y" /* yacc.c:1646  */
    { domain->solInfo().filterSnapshotRows = bool((yyvsp[0].ival)); }
#line 15266 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1400:
#line 6170 "p.y" /* yacc.c:1646  */
    { domain->solInfo().selectFullNode = bool((yyvsp[0].ival)); }
#line 15272 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1401:
#line 6172 "p.y" /* yacc.c:1646  */
    { domain->solInfo().selectFullElem = bool((yyvsp[0].ival)); }
#line 15278 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1402:
#line 6174 "p.y" /* yacc.c:1646  */
    { domain->solInfo().computeForceSnap = bool((yyvsp[0].ival)); }
#line 15284 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1403:
#line 6176 "p.y" /* yacc.c:1646  */
    { domain->solInfo().computeConstraintSnap = bool((yyvsp[0].ival)); }
#line 15290 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1404:
#line 6178 "p.y" /* yacc.c:1646  */
    { domain->solInfo().orthogForceSnap = bool((yyvsp[0].ival)); }
#line 15296 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1405:
#line 6180 "p.y" /* yacc.c:1646  */
    { domain->solInfo().orthogConstraintSnap = bool((yyvsp[0].ival)); }
#line 15302 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1406:
#line 6182 "p.y" /* yacc.c:1646  */
    { domain->solInfo().npMax = (yyvsp[0].ival); }
#line 15308 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1407:
#line 6184 "p.y" /* yacc.c:1646  */
    { domain->solInfo().scpkMB= (yyvsp[-1].ival);
    domain->solInfo().scpkNB= (yyvsp[0].ival); }
#line 15315 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1408:
#line 6187 "p.y" /* yacc.c:1646  */
    { domain->solInfo().scpkMP= (yyvsp[-1].ival);
    domain->solInfo().scpkNP= (yyvsp[0].ival); }
#line 15322 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1409:
#line 6190 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useReverseOrder = bool((yyvsp[0].ival)); }
#line 15328 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1410:
#line 6192 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf = 1;
    domain->solInfo().nmfNumROBDim = (yyvsp[-4].ival);
    domain->solInfo().nmfDelROBDim = (yyvsp[-3].ival);
    domain->solInfo().nmfRandInit = (yyvsp[-2].ival);
    domain->solInfo().nmfMaxIter = (yyvsp[-1].ival);
    domain->solInfo().nmfTol = (yyvsp[0].fval); }
#line 15339 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1411:
#line 6199 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf = 1;
    domain->solInfo().nmfMaxIter = (yyvsp[-1].ival);
    domain->solInfo().nmfTol = (yyvsp[0].fval); }
#line 15347 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1413:
#line 6207 "p.y" /* yacc.c:1646  */
    { domain->solInfo().conwepConfigurations.push_back((yyvsp[-1].blastData)); }
#line 15353 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1414:
#line 6212 "p.y" /* yacc.c:1646  */
    { domain->solInfo().scalePosCoords = true;
     domain->solInfo().xScaleFactor = (yyvsp[-3].fval);
     domain->solInfo().yScaleFactor = (yyvsp[-2].fval);
     domain->solInfo().zScaleFactor = (yyvsp[-1].fval);}
#line 15362 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1415:
#line 6221 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePOSCFG = true; }
#line 15368 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1416:
#line 6223 "p.y" /* yacc.c:1646  */
    { domain->solInfo().xScaleFactors.push_back((yyvsp[-3].fval));
     domain->solInfo().yScaleFactors.push_back((yyvsp[-2].fval));
     domain->solInfo().zScaleFactors.push_back((yyvsp[-1].fval)); }
#line 15376 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1417:
#line 6227 "p.y" /* yacc.c:1646  */
    { domain->solInfo().MassOrthogonalBasisFiles.push_back((yyvsp[-1].strval)); }
#line 15382 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1418:
#line 6232 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePOSCFG = true; }
#line 15388 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1419:
#line 6234 "p.y" /* yacc.c:1646  */
    { domain->solInfo().NodeTrainingFiles.push_back(std::string((yyvsp[-1].slist).v[0]));
     for(int i=1; i<(yyvsp[-1].slist).nval; ++i) domain->solInfo().MassOrthogonalBasisFiles.push_back(std::string((yyvsp[-1].slist).v[i])); }
#line 15395 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1420:
#line 6240 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true;
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().ROMPostProcess = true; }
#line 15403 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1422:
#line 6248 "p.y" /* yacc.c:1646  */
    { domain->solInfo().RODConversionFiles.push_back((yyvsp[0].strval)); 
    domain->solInfo().numRODFile += 1; }
#line 15410 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1423:
#line 6251 "p.y" /* yacc.c:1646  */
    { domain->solInfo().RODConversionFiles.push_back((yyvsp[-1].strval));
    domain->solInfo().numRODFile += 1; 
    domain->solInfo().skipPodRom = (yyvsp[0].ival);}
#line 15418 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1424:
#line 6255 "p.y" /* yacc.c:1646  */
    { domain->solInfo().romresidType = (yyvsp[0].ival); }
#line 15424 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1425:
#line 6260 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[0].ival); }
#line 15430 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1426:
#line 6262 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = std::numeric_limits<int>::max(); }
#line 15436 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1427:
#line 6267 "p.y" /* yacc.c:1646  */
    { (yyval.fval) = (yyvsp[0].ival); }
#line 15442 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1428:
#line 6269 "p.y" /* yacc.c:1646  */
    { (yyval.fval) = (yyvsp[0].fval); }
#line 15448 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1429:
#line 6271 "p.y" /* yacc.c:1646  */
    { (yyval.fval) = std::numeric_limits<double>::infinity(); }
#line 15454 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1430:
#line 6273 "p.y" /* yacc.c:1646  */
    { (yyval.fval) = std::numeric_limits<double>::epsilon(); }
#line 15460 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;


#line 15464 "/home/anarkhede/tinkercliffs/FEMWorkingFoam/Parser.d/parser.cpp" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 6275 "p.y" /* yacc.c:1906  */

