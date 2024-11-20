%x incl

%{
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stack>
#include <string>
#include <Parser.d/AuxDefs.h>
#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>
#include "parser.hpp"
int line_num = 1;

// used for include
std::stack<YY_BUFFER_STATE> includedFiles;
std::stack<int> lineNumber; // for fancy error message

%}

%e 40000
%p 20000
%a 40000
%n 20000
%k 40000
%o 30000


DIGIT [0-9]
HEXDIGIT [0-9A-Fa-f]
SPACES [ \t]+
NAME [a-zA-Z_][a-zA-Z_0-9]*
TRAIL [a-zA-Z]*
%%
"INCLUDE"|"include" { BEGIN(incl);   } // go to include state 
<incl>[ \t]* 
<incl>[^ \t\n\r]+ {
        int c = yyinput();
        while(c == '\t' || c == ' ' || c == '\r')
          c = yyinput();
        if(c != '\n')
          yyunput(c,yytext_ptr);
        std::string fname = std::string(yytext);
        if(fname.length() > 0)
        {
          if(fname[0]=='"' && fname[ fname.length()-1 ]=='"' && fname.length()>2) {
            fname = fname.substr(1,fname.length()-2);
          }

          FILE * f = fopen(fname.c_str(),"r");
          if(f!=NULL)
          { // open the file and substitute it as the current buffer
            includedFiles.push(YY_CURRENT_BUFFER);
            lineNumber.push(line_num);
            line_num = 1;
            yy_switch_to_buffer( yy_create_buffer(f,YY_BUF_SIZE) );
          }
          else
          {
            std::cout << " **************************************************** " << std::endl 
                      << " ** ERROR :  Couldn't open file " << fname << std::endl
                      << " **************************************************** " << std::endl; 
             exit(-1);
          }
          BEGIN(INITIAL); // go back to parsing
        }  
      }
<<EOF>> {
      if(includedFiles.size()==0)
      { // end of everything
        yyterminate();
      }
      else
      { // switch to previous buffer
        yy_delete_buffer( YY_CURRENT_BUFFER );
        yy_switch_to_buffer( includedFiles.top());
        line_num = lineNumber.top()+1;
        lineNumber.pop();
        includedFiles.pop();
      }
    }

"NODE"[a-zA-Z]*|"node"[a-zA-Z]*    { return NODETOKEN;    }

"ZERO"[a-zA-Z]*|"zero"[a-zA-Z]*    { return ZERO;    }

"READMODE"[a-zA-Z]*|"readmode"[a-zA-Z]* { return READMODE; }

"MNORM"|"mnorm"                    { yylval.ival = ModalParams::Mnorm; return ROBTYPE; }
"INORM"|"inorm"                    { yylval.ival = ModalParams::Inorm; return ROBTYPE; }
"NONEG"|"noneg"                    { yylval.ival = ModalParams::Noneg; return ROBTYPE; }

"DUALBASIS"|"dualbasis"            { return DUALBASIS; }

"ADJOINTBASIS"|"adjointbasis"      { return ADJOINTBASIS; }

"USE_MASS_NORMALIZED_BASIS"|"use_mass_normalized_basis" { return USEMASSNORMALIZEDBASIS; }

"USE_CONSTANT_MASS_FORCES"|"use_constant_mass_forces" { return USECONSTANTMASS; }

"ONLINE_MASS_NORMALIZE_BASIS"|"online_mass_normalize_basis" { return ONLINEMASSNORMALIZEBASIS;}

"STACKED_TRAINING"|"stacked_training" { return STACKED; }

"PRES"[a-zA-Z]*|"pres"[a-zA-Z]*    { return PRESSURE; }

"FACE"|"face"                       { return FACE; }

"CONW"[a-zA-Z]*|"conw"[a-zA-Z]*    { return CONWEP; }

"CONSTANT"|"constant"              { return CONSTANT; }

"PREL"[a-zA-Z]*|"prel"[a-zA-Z]*    { return PRELOAD; }

"BUCK"[a-zA-Z]*|"buck"[a-zA-Z]*    { return BUCKLE; }

"CONTROL"[a-zA-Z]*|"control"[a-zA-Z]*    { return CONTROL;  }

"TOPO"[a-zA-Z]*|"topo"|"topology"  { return TOPOLOGY; }

"SOWE"[a-zA-Z]*|"sowe"[a-zA-Z]*    { return SOWER; }

"BINARY"|"binary"                  { return BINARY; }

"BINARYINPUT"|"binaryinput"        { return BINARYINPUT; }

"BINARYOUTPUT"|"binaryoutput"      { return BINARYOUTPUT; }

"e"                                { yylval.ival = 0; return KEYLETTER; }
"c"                                { yylval.ival = 1; return KEYLETTER; }

"LOCF"[a-zA-Z]*|"locf"[a-zA-Z]*    { yylval.ival = OutputInfo::Local;  return OUTPUTFRAME; }
"NODF"[a-zA-Z]*|"nodf"[a-zA-Z]*    { yylval.ival = OutputInfo::Local;  return OUTPUTFRAME; }
"GLOF"[a-zA-Z]*|"glof"[a-zA-Z]*    { yylval.ival = OutputInfo::Global; return OUTPUTFRAME; }
"MATF"[a-zA-Z]*|"matf"[a-zA-Z]*    { yylval.ival = OutputInfo::Material; return OUTPUTFRAME; }

"ON"|"on"|"On"                     { yylval.ival = 1; return SWITCH; }
"OFF"|"off"|"Off"                  { yylval.ival = 0; return SWITCH; }
"YES"|"yes"|"Yes"                  { yylval.ival = 1; return SWITCH; }
"NO"|"no"|"No"                     { yylval.ival = 0; return SWITCH; }
"TRUE"|"true"|"True"               { yylval.ival = 1; return SWITCH; }
"FALSE"|"false"|"False"            { yylval.ival = 0; return SWITCH; }

"INF"[a-zA-Z]*|"inf"[a-zA-Z]*      { return INFINTY; }
"EPS"[a-zA-Z]*|"eps"[a-zA-Z]*      { return EPSILON; }
"MAX"|"max"                        { return MAXIMUM; }

"DECO"[a-zA-Z]{0,7}|"deco"[a-zA-Z]{0,7} { return DECOMPOSE; }

"nsubs"|"NSUBS"                    { return NSUBS; }

"outfile"|"OUTFILE"                { return DECOMPFILE; }

"exit"|"EXIT"                      { return EXITAFTERDEC; }

"skip"|"SKIP"                      { return SKIP; }

"outload"|"OUTLOAD"                { return OUTPUTWEIGHT; }

"outmem"|"OUTMEM"                  { return OUTPUTMEMORY; }

"deter"|"DETER"                    { return DETER; }

"WEIGHTS"|"weights"                { return WEIGHTLIST; }

"FWEI"[a-zA-Z]*|"fwei"[a-zA-Z]*    { return FIELDWEIGHTLIST; }

"trivial"|"TRIVIAL"                { return TRIVIAL; }

"allow_mechanisms"|"ALLOW_MECHANISMS" { return ALLOWMECHANISMS; }

"use_scotch"|"USE_SCOTCH"          { return USESCOTCH; }

"decomposition"|"DECOMPOSITION"    { return DECOMPOSITION; }

"global"|"GLOBAL"                  { return GLOBAL; }

"matcher"|"MATCHER"|"Matcher"      { return MATCHER; }

"geometry"|"GEOMETRY"              { return GEOMETRY; }

"cpumap"|"CPUMAP"                  { return CPUMAP; }

"TOPF"[a-zA-Z]*|"topf"[a-zA-Z]*    { return TOPFILE;  }

"ATTR"[a-zA-Z]*|"attr"[a-zA-Z]*    { return ATTRIBUTES; }

"THETA"|"theta"                    { return THETA; }

"HRC"|"hrc"			   { return HRC; }

"PROJECT"|"project"		   { return PROJSOL; }
 
"POSITIVE"|"positive"		   { return POSELEM; }

"HOTSTART"|"hotstart"		   { return HOTSTART; }

"CENT"[a-zA-Z]*|"cent"[a-zA-Z]*    { return CENTER; }

"REVORD"|"revord"                  { return REVERSEORDER; } 

"REDFOL"|"redfol"                  { return REDFOL; }

"EXTFOL"|"extfol"                  { return EXTFOL; }

"FULLND"|"fullnd"		   { return PSEUDOGNAT; }

"FULLEL"|"fullel"		   { return PSEUDOGNATELEM; }

"ELLUMP"                           { return ELLUMP; }

"LOCROB"|"LocRob"|"locrob"         { return LOCALREDUCEDORDERBASES; }

"AUXI"[a-zA-Z]*|"auxi"[a-zA-Z]*    { return AUXILIARY; }

"IDENTITY"|"identity"              { return IDENTITY; }

"EFRA"[a-zA-Z]*|"efra"[a-zA-Z]*    { return EFRAMES; }

"NFRA"[a-zA-Z]*|"nfra"[a-zA-Z]*    { return NFRAMES; }

"CSFR"[a-zA-Z]*|"csfr"[a-zA-Z]*    { return CSFRAMES; }

"RECT"[a-zA-Z]*|"rect"[a-zA-Z]*    { yylval.ival = NFrameData::Rectangular; return FRAMETYPE; }
"CYLI"[a-zA-Z]*|"cyli"[a-zA-Z]*    { yylval.ival = NFrameData::Cylindrical; return FRAMETYPE; }
"SPHE"[a-zA-Z]*|"sphe"[a-zA-Z]*    { yylval.ival = NFrameData::Spherical;   return FRAMETYPE; }

"THIRDNODE"|"thirdnode"            { return THIRDNODE; }

"BOFF"[a-zA-Z]*|"boff"[a-zA-Z]*    { return BOFFSET; }

"idisp6pita"|"IDISP6PITA"          { return PITADISP6; }

"ivel6pita"|"IVEL6PITA"            { return PITAVEL6; }        

"IDISP"[a-zA-Z]*|"idisp"[a-zA-Z]*  { return IDIS;  }

"IDISP6"[a-zA-Z]*|"idisp6"[a-zA-Z]* { return IDIS6; }

"GEPS"[a-zA-Z]*|"geps"[a-zA-Z]*    { return GEPS; }

"IVEL"[a-zA-Z]*|"ivel"[a-zA-Z]*    { return IVEL; }

"IACC"[a-zA-Z]*|"iacc"[a-zA-Z]*    { return IACC; }

"ITEM"[a-zA-Z]*|"item"[a-zA-Z]*    { return ITEMP; }

"ETEM"[a-zA-Z]*|"etem"[a-zA-Z]*    { return ETEMP; }

"TEMP"[a-zA-Z]*|"temp"[a-zA-Z]*    { return TEMP; }

"PDIR"[a-zA-Z]*|"pdir"[a-zA-Z]*    { return PDIR; }

"HEFR"[a-zA-Z]*|"hefr"[a-zA-Z]*    { return HEFRS; }

"FLUX"[a-zA-Z]*|"flux"[a-zA-Z]*    { return FLUX; } 

"END"|"end"                        { return END;  }

"LOAD"|"load"                      { return LOAD; }

"LOADC"[a-zA-Z]*|"loadc"[a-zA-Z]*  { return LOADCASE; }

"USE"|"use"                        { return USE;  }

"SENS"[a-hj-zA-HJ-Z]*|"sens"[a-hj-zA-HJ-Z]*    { return SENSORS; }

"ACTU"[a-zA-Z]*|"actu"[a-zA-Z]*    { return ACTUATORS; }

"MATE"[a-zA-Z]*|"mate"[a-zA-Z]*    { return MATERIALS; }

"DISP"|"disp"|"DISPL"|"displ"|"DISPLA"|"displa"|"DISPLAC"|"displac"|"DISPLACE"|"displace"|"DISPLACEM"[a-zA-Z]*|"displacem"[a-zA-Z]* { return DISP; }

"THRU"[a-zA-Z]*|"thru"[a-zA-Z]*    { return THRU; }

"STEP"[a-zA-Z]*|"step"[a-zA-Z]*    { return STEP; }

"SURFACE"|"surface"                { return SURF; }

"LMPC"[a-zA-Z]*|"lmpc"[a-zA-Z]*    { return LMPC; }

"HLMPC"[a-zA-Z]*|"hlmpc"[a-zA-Z]*  { return HLMPC; }

"USING"|"using"                    { return USING; }

"CRHS"|"crhs"                      { return CRHS; }

"FORC"[a-zA-Z]*|"forc"[a-zA-Z]*    { return FORCE; }

"axial"|"AXIAL"                    { yylval.ival = BCond::Axial;      return MOMENTTYPE; }

"rotational"|"ROTATIONAL"          { yylval.ival = BCond::Rotational; return MOMENTTYPE; }

"follower"|"FOLLOWER"              { yylval.ival = BCond::Follower;   return MOMENTTYPE; }

"COND"[a-zA-Z]*|"cond"[a-zA-Z]*    { return CONDITION; }

"MASS"[a-zA-Z]*|"mass"[a-zA-Z]*    { return MASS; }

"QLOA"[a-zA-Z]*|"qloa"[a-zA-Z]*    { return QLOAD; }

"PLOA"[a-zA-Z]*|"ploa"[a-zA-Z]*    { return PLOAD; }

"STATI"[a-zA-Z]*|"stati"[a-zA-Z]*|"STAT"|"stat"    { return STATS; }

"CASE"[a-zA-Z]*|"case"[a-zA-Z]*    { return CASES; }

"PIEC"[a-zA-Z]*|"piec"[a-zA-Z]*    { return PIECEWISE; }

"ELSCATT"[a-zA-Z]*|"elscatt"[a-zA-Z]* { return ELSCATTERER; }

"FFP"|"ffp"                        {return FFP; }

"FFPDIR"|"ffpdir"                  {return FFPDIR; }

"KIRLOC"|"KIRLOC"                  {return KIRLOC; }

"HSOM"[a-zA-Z]*|"hsom"[a-zA-Z]*    { return HSOMMERFELD; }

"SCATT"[a-zA-Z]*|"scatt"[a-zA-Z]*  { return SCATTERER; }

"IHDI"[a-zA-Z]*|"ihdi"[a-zA-Z]*    { return IHDIRICHLET; }

"IHNE"[a-zA-Z]*|"ihne"[a-zA-Z]*    { return IHNEUMANN; }

"IHDS"[a-zA-Z]*|"ihds"[a-zA-Z]*    { return IHDSWEEP; }

"FACO"[a-zA-Z]*|"faco"[a-zA-Z]*    { return FACOUSTICS; }

"ELHSOM"[a-zA-Z]*|"elhsom"[a-zA-Z]* { return ELHSOMMERFELD; }

"HPNEU"[a-zA-Z]*|"hpneu"[a-zA-Z]*  { return HNEUMAN; }

"TURKEL"|"turkel"                  { return TURKEL; }

"SHIFT"|"shift"                    { return SHIFT; }

"IMPE"[a-zA-Z]*|"impe"[a-zA-Z]*    { return IMPE; }

"freq"|"FREQ"|"frequency"|"FREQUENCY" { return FREQ; }

"freqsweep"|"FREQSWEEP"            { return FREQSWEEP; }

"freqsweepa"|"FREQSWEEPA"            { return FREQSWEEPA; }

"freqsweepaw"|"FREQSWEEPAW"            { return FREQSWEEPAW; }

"freqsweep1"|"FREQSWEEP1"          { return FREQSWEEP1; }

"freqsweep2"|"FREQSWEEP2"          { return FREQSWEEP2; }

"RECONS"[a-zA-Z]*|"recons"[a-zA-Z]* { return RECONS; }

"Taylor"|"taylor"|"TAYLOR"         { yylval.ival = 0; return RECONSALG; }

"Pade1"|"pade1"|"PADE1"            { yylval.ival = 1; return RECONSALG; }

"Pade"|"pade"|"PADE"               { yylval.ival = 2; return RECONSALG; }

"padepivot"|"PADEPIVOT"            { return PADEPIVOT; }

"Fourier"|"fourier"|"FOURIER"      { yylval.ival = 3; return RECONSALG; }

"PadeLanczos"|"padelanczos"|"PADELANCZOS" { yylval.ival = 4; return RECONSALG; }

"GalProjection"|"galprojection"|"GALPROJECTION" { yylval.ival = 5; return RECONSALG; }

"KrylovGalProjection"|"krylovgalprojection"|"KRYLOVGALPROJECTION" { yylval.ival= 6; return RECONSALG; }

"WCAWEGalProjection"|"wcawegalprojection"|"WCAWEGALPROJECTION" { yylval.ival = 7; return RECONSALG; }

"padepoles"|"PADEPOLES"            { return PADEPOLES; }

"modefilter"|"MODEFILTER"          { return MODEFILTER; }

"test_ulrich"                      { return TESTULRICH; }

"numcgm"|"NUMCGM"                  { return NUMCGM; }

"interfacelumped"|"INTERFACELUMPED" { return INTERFACELUMPED; }

"MODE"|"mode"                      { return MODE; }

"CG"|"cg"|"PCG"|"pcg"              { yylval.ival = 0; return ITERTYPE; }

"GMRES"|"gmres"                    { yylval.ival = 1; return ITERTYPE; }

"GCR"|"gcr"                        { yylval.ival = 2; return ITERTYPE; }

"BCG"|"bcg"                        { yylval.ival = 4; return ITERTYPE; }

"CR"|"cr"                          { yylval.ival = 5; return ITERTYPE; }

"FETIH"|"fetih"                    { return FETIH; }

"HFETI"|"hfeti"                    { return HFETI; }

"tolcgm"|"TOLCGM"                  { return TOLCGM; }

"spacedim"|"SPACEDIM"              { return SPACEDIMENSION; }

"krylovtype"|"KRYLOVTYPE"          { return KRYLOVTYPE; }

"savememcoarse"|"SAVEMEMCOARSE"    { return SAVEMEMCOARSE; }

"HELMHOLTZ"|"helmholtz"            { return HELMHOLTZ; }

"DEM"                              { return DEM; }

"bgtl"|"BGTL"                      { return BGTL; }

"incidence"|"INCIDENCE"            { return INCIDENCE; }

"pointsource"|"POINTSOURCE"        { return POINTSOURCE; }

"planewave"|"PLANEWAVE"            { return PLANEWAVE; }

"HDIR"[a-zA-Z]*|"hdir"[a-zA-Z]*    { return HDIRICHLET; }

"HNEU"[a-zA-Z]*|"hneu"[a-zA-Z]*    { return HNEUMAN; }

"HNBO"[a-zA-Z]*|"hnbo"[a-zA-Z]*    { return HNBO; }

"HDNB"[a-zA-Z]*|"hdnb"[a-zA-Z]*    { return HNBO; }

"HARB"[a-zA-Z]*|"harb"[a-zA-Z]*    { return ELHSOMMERFELD; }

"HWIB"[a-zA-Z]*|"hwib"[a-zA-Z]*    { return HWIBO; }

"HELMMF"|"helmmf"                  { return HELMMF; }

"HELMSO"[a-zA-Z]*|"helmso"[a-zA-Z]* { return HELMSO; }

"HSCB"[a-zA-Z]*|"hscb"[a-zA-Z]*    { return HSCBO; }

"ATDARB"[a-zA-Z]*|"atdarb"[a-zA-Z]* { return ATDARB; }

"ATDNEU"[a-zA-Z]*|"atdneu"[a-zA-Z]* { return ATDNEU; }

"ATDDNB"[a-zA-Z]*|"atddnb"[a-zA-Z]* { return ATDDNB; }

"ATDDIR"[a-zA-Z]*|"atddir"[a-zA-Z]* { return ATDDIR; }

"ATDROB"[a-zA-Z]*|"atdrob"[a-zA-Z]* { return ATDROB; }

"BLOCKDIAG" 			   { return BLOCKDIAG; }

"FETI"|"feti"                      { return FETI; }

"OLD"|"old"                        { yylval.ival = 0; return FETI2TYPE; }

"NEW"|"new"                        { yylval.ival = 1; return FETI2TYPE; }

"DP"|"dp"                          { return DP; }

"DPH"|"dph"                        { return DPH; }

"LIB"|"lib"                        { return LIB; }

"noco"[a-zA-Z]*|"NOCO"[a-zA-Z]*    { return NOCOARSE; }

"proj"[a-zA-Z]*|"PROJ"[a-zA-Z]*    { return PROJ; }

"interf_solver"|"INTERF_SOLVER"    { return OUTERLOOP; }

"hermitian"                        { return HERMITIAN; }

"local_solver"|"LOCAL_SOLVER"      { return LOCALSOLVER; }

"coarse_solver"|"COARSE_SOLVER"    { return COARSESOLVER; }

"aux_coarse_solver"                { return AUXCOARSESOLVER; }

"cct_solver"|"CCT_SOLVER"          { return CCTSOLVER; }

"rebuildcct"|"REBUILDCCT"          { return REBUILDCCT; }

"solverc"[a-zA-Z]*|"SOLVERC"[a-zA-Z]* { return SOLVERCNTL; }

"solverh"[a-zA-Z]*|"SOLVERH"[a-zA-Z]* { return SOLVERHANDLE; }

"skyl"[a-zA-Z]*|"SKYL"[a-zA-Z]*    { yylval.ival = 0; return SOLVERTYPE; }

"spar"[a-zA-Z]*|"SPAR"[a-zA-Z]*    { yylval.ival = 1; return SOLVERTYPE; }

"bloc"[a-zA-Z]*|"BLOC"[a-zA-Z]*    { yylval.ival = 2; return SOLVERTYPE; }

"llt"[a-zA-Z]*|"LLT"[a-zA-Z]*      { yylval.ival = 3; return SOLVERTYPE; }

"ldlt"[a-zA-Z]*|"LDLT"[a-zA-Z]*    { yylval.ival = 4; return SOLVERTYPE; }

"splu"[a-zA-Z]*|"SPLU"[a-zA-Z]*    { yylval.ival = 15; return SOLVERTYPE; }

"ssqr"[a-zA-Z]*|"SSQR"[a-zA-Z]*    { yylval.ival = 16; return SOLVERTYPE; }

"spqr"[a-zA-Z]*|"SPQR"[a-zA-Z]*    { yylval.ival = 17; return SOLVERTYPE; }

"chol"[a-zA-Z]*|"CHOL"[a-zA-Z]*    { yylval.ival = 5; return SOLVERTYPE; }

"umfp"[a-zA-Z]*|"UMFP"[a-zA-Z]*    { yylval.ival = 6; return SOLVERTYPE; }

"supe"[a-zA-Z]*|"SUPE"[a-zA-Z]*    { yylval.ival = 7; return SOLVERTYPE; }

"spoo"[a-zA-Z]*|"SPOO"[a-zA-Z]*    { yylval.ival = 8; return SOLVERTYPE; }

"mump"[a-zA-Z]*|"MUMP"[a-zA-Z]*    { yylval.ival = 9; return SOLVERTYPE; }

"native_cholmod"                   { yylval.ival = 18; return SOLVERTYPE; }
"pardiso"|"PARDISO"                { yylval.ival = 19; return SOLVERTYPE; }

"purelydiag"[a-zA-Z]*|"PURELYDIAG"[a-zA-Z]* { yylval.ival = 10; return SOLVERTYPE; }

"dbsgal"[a-zA-Z]*|"DBSGAL"[a-zA-Z]* { yylval.ival = 12; return SOLVERTYPE; }

"eisgal"[a-zA-Z]*|"EISGAL"[a-zA-Z]* { yylval.ival = 13; return SOLVERTYPE; }

"gold"[a-zA-Z]*|"GOLD"[a-zA-Z]*    { yylval.ival = 14; return SOLVERTYPE; }

"verb"[a-zA-Z]*|"VERB"[a-zA-Z]*    { return VERBOSE; }

"print_number"                     { return PRINTNUMBER; }

"mumps_icntl"|"MUMPS_ICNTL"        { return MUMPSICNTL; }

"mumps_cntl"|"MUMPS_CNTL"          { return MUMPSCNTL; }

"mumps_mineq"|"MUMPS_MINEQ"        { return MUMPSMINEQ; }

"mumps_stride"|"MUMPS_STRIDE"      { return MUMPSSTRIDE; }

"pivot"[a-zA-Z]*|"PIVOT"[a-zA-Z]*  { return PIVOT; }

"scaled"[a-zA-Z]*|"SCALED"[a-zA-Z]* { return SCALED; }

"spooles_tau"|"SPOOLES_TAU"        { return SPOOLESTAU; }

"spooles_maxsize"|"SPOOLES_MAXSIZE" { return SPOOLESMAXSIZE; }

"spooles_maxdomainsize"|"SPOOLES_MAXDOMAINSIZE" { return SPOOLESMAXDOMAINSIZE; }

"spooles_seed"|"SPOOLES_SEED"      { return SPOOLESSEED; }

"spooles_maxzeros"|"SPOOLES_MAXZEROS" { return SPOOLESMAXZEROS; }

"spooles_msglvl"|"SPOOLES_MSGLVL"  { return SPOOLESMSGLVL; }

"spooles_scale"|"SPOOLES_SCALE"    { return SPOOLESSCALE; }

"spooles_renum"|"SPOOLES_RENUM"    { return SPOOLESRENUM; }

"sparse_maxsup"|"SPARSE_MAXSUP"    { return SPARSEMAXSUP; }

"sparse_defblk"|"SPARSE_DEFBLK"    { return SPARSEDEFBLK; }

"sparse_renum"|"SPARSE_RENUM"      { return SPARSERENUM; }

"goldfarb_tol"|"GOLDFARB_TOL"      { return GOLDFARBTOL; }

"ilu_"[a-zA-Z]*|"ILU_"[a-zA-Z]*    { return ILUDROPTOL; }

"kryp"[a-zA-Z]*|"KRYP"[a-zA-Z]*    { return NLPREC; }

"tolf"[a-zA-Z]*|"TOLF"[a-zA-Z]     { return TOLFETI; }

"stag"[a-zA-Z_]*|"STAG"[a-zA-Z_]*  { return STAGTOL; }

"global_rbm_tol"|"GLOBAL_RBM_TOL"  { return GLOBALTOL; }

"global_cor_rbm_tol"|"GLOBAL_COR_RBM_TOL" { return GLOBALCRBMTOL; }

"cct_tol"|"CCT_TOL"                { return CCTTOL; }

"maxor"[a-zA-Z]*|"MAXOR"[a-zA-Z]*  { return MAXORTHO; }

"vers"[a-zA-Z]*|"VERS"[a-zA-Z]*    { return VERSION; }

"USDF"[a-zA-Z]*|"usdf"[a-zA-Z]*    { return USERDEFINEFORCE; }

"USDD"[a-zA-Z]*|"usdd"[a-zA-Z]*    { return USERDEFINEDISP; }

"PRINTMAT"[a-zA-Z]*                { return PRINTMATLAB; }

"printmat"[a-zA-Z]*                { return printmatlab; }

"u_proj"[a-zA-Z]*|"U_PROJ"[a-zA-Z]* { return UPROJ; }

"mpc_type"|"MPC_TYPE"              { return MPCTYPE; }

"dual"                             { yylval.ival = 1; return MPCTYPEID; }

"primal"|"PRIMAL"                  { yylval.ival = 2; return MPCTYPEID; }

"mpc_precn"[a-zA-Z]*|"MPC__PRECN"[a-zA-Z]* { return MPCPRECNO; }

"diag"[a-zA-Z]*|"DIAG"[a-zA-Z]*    { yylval.ival = 1; return MPCPRECNOID; }

"full"[a-zA-Z]*|"FULL"[a-zA-Z]*    { yylval.ival = 2; return MPCPRECNOID; }

"tblock"[a-zA-Z]*|"TBLOCK"[a-zA-Z]* { yylval.ival = 3; return MPCPRECNOID; }

"sblock"[a-zA-Z]*|"SBLOCK"[a-zA-Z]* { yylval.ival = 4; return MPCPRECNOID; }

"mblock"[a-zA-Z]*|"MBLOCK"[a-zA-Z]* { yylval.ival = 5; return MPCPRECNOID; }

"blk_overlap"[a-zA-Z]*|"BLK_OVERLAP"[a-zA-Z]* { return MPCBLK_OVERLAP; }

"c_normalize" { return CNORM; }

"mpc_check" { return MPCCHECK; }

"gmres_residual"|"GMRES_RESIDUAL"  { return GMRESRESIDUAL; }

"pick_any_corner" { return PICKANYCORNER; }

"DIRE"[a-zA-Z]*|"dire"[a-zA-Z]*    { return DIRECT; } // deprecated
"ELIM"[a-zA-Z]*|"elim"[a-zA-Z]*    { return DIRECT; }

"MULT"[a-zA-Z]*|"mult"[a-zA-Z]*    { return MULTIPLIERS; }

"CONSTRAINTS"|"constraints"        { return CONSTRAINTS; }

"HESSI"[a-zA-Z]*|"hessi"[a-zA-Z]*  { return HESSIAN; }

"PENA"[a-zA-Z]*|"pena"[a-zA-Z]*    { return PENALTY; }

"AUGMENTE"[a-zA-Z]*|"augmente"[a-zA-Z]* { return AUGMENTED; }

"RENU"[a-zA-Z]*|"renu"[a-zA-Z]*    { return RENUM; }

"DYNA"[a-zA-Z]*|"dyna"[a-zA-Z]*    { return DYNAM; }

"STAB"[a-zA-Z]*|"stab"[a-zA-Z]*    { return STABLE; }

"DELE"[a-zA-Z_]*|"dele"[a-zA-Z_]*  { return DELETEELEMENTS; }

"CHEC"[a-zA-Z]*|"chec"[a-zA-Z]*    { return CHECKTOKEN; }

"rbmf"[a-zA-Z]*|"RBMF"[a-zA-Z]*    { return RBMFILTER; }

"lisr"[a-zA-Z]*|"LISR"[a-zA-Z]*    { return LISRBM; }

"RESTART"|"restart"                { return RESTART; }

"CONVECTI"[a-zA-Z]*|"convecti"[a-zA-Z]* { return CONVECTION; }

"RADI"[a-zA-Z]*|"radi"[a-zA-Z]*    { return RADIATION; }

"newm"[a-zA-Z]*|"NEWM"[a-zA-Z]*    { return NEWMARK; }

"QSTA"[a-zA-Z]*|"qsta"[a-zA-Z]*    { return QSTATIC; }

"QUASISTATIC"|"quasistatic"        { return QUASISTATIC; }

"modal"|"MODAL"|"srom"|"SROM"      { return MODAL; }

"contact"|"CONTACT"                { return CONTACT; }

"precn"o[a-zA-Z]*|"PRECN"O[a-zA-Z]* { return PRECNO; }

"prectype"|"PRECTYPE"              { return PRECTYPE; }

"nonshifted"|"NONSHIFTED"          { yylval.ival = 0; return PRECTYPEID; }

"shifted"|"SHIFTED"                { yylval.ival = 1; return PRECTYPEID; }

"tolp"[a-zA-Z]*|"TOLP"[a-zA-Z]*    { return TOLPCG; }

"tolb"[a-zA-Z]*|"TOLB"[a-zA-Z]*    { return TOLPCG; }

"tolcr"|"TOLCR"                    { return TOLPCG; }

"tolg"[a-zA-Z]*|"TOLG"[a-zA-Z]*    { return TOLPCG; }

"subt"[a-zA-Z]*|"SUBT"[a-zA-Z]*    { return SUBTYPE; }

"maxel"[a-zA-Z]*|"MAXEL"[a-zA-Z]*  { return MAXELEM; }

"maxit"[a-zA-Z]*|"MAXIT"[a-zA-Z]*  { return MAXITR; }

"propo"[a-zA-Z]*|"PROPO"[a-zA-Z]*  { return PROPORTIONING; }

"lines"[a-zA-Z]*|"LINES"[a-zA-Z]*  { return LINESEARCH; }

"fails"[a-zA-Z]*|"FAILS"[a-zA-Z]*  { return FAILSAFE; }

"debug_icntl"                      { return DEBUGICNTL; }

"debug_cntl"                       { return DEBUGCNTL; }

"ptol"[a-zA-Z]*|"PTOL"[a-zA-Z]*    { return PTOL; }

"pmax"[a-zA-Z]*|"PMAX"[a-zA-Z]*    { return PMAXIT; }

"plan"[a-zA-Z]*|"PLAN"[a-zA-Z]*    { return PLANTOL; }

"nsbs"[a-zA-Z]*|"NSBS"[a-zA-Z]*    { return NSBSPV; }

"neig"[a-zA-Z]*|"NEIG"[a-zA-Z]*    { return NEIGPA; }

"toleig"[a-zA-Z]*|"TOLEIG"[a-zA-Z]* { return TOLEIG; }

"tolj"[a-zA-Z]*|"TOLJ"[a-zA-Z]*    { return TOLJAC; }

"tolsen"[a-zA-Z]*|"TOLSEN"[a-zA-Z]* { return TOLSEN; }

"relaxationsen"|"RELAXATIONSEN"    { return RELAXATIONSEN; }

"ratiotolsen"|"RATIOTOLSEN"        { return RATIOTOLSEN; }

"ksparam"|"KSPARAM"                { return KSPARAM; }

"ksmax"|"KSMAX"                    { return KSMAX; }

"qrfactorization"|"QRFACTORIZATION" { return QRFACTORIZATION; }

"expl"[a-zA-Z]*|"EXPL"[a-zA-Z]*    { return EXPLICIT; }

"lobpcg"[a-zA-Z]*|"LOBPCG"[a-zA-Z]* { return LOBPCG; }

"arpack"|"ARPACK"                  { return ARPACK; }

"filtereig"|"FILTEREIG"            { return FILTEREIG; }

"jacobi"|"JACOBI"                  { return JACOBI; }

"damp"[a-zA-Z]*|"DAMP"[a-zA-Z]*    { yylval.ival = 0; return DAMPING; }

"rayd"[a-zA-Z]*|"RAYD"[a-zA-Z]*    { yylval.ival = 1; return DAMPING; }

"strd"[a-zA-Z]*|"STRD"[a-zA-Z]*    { yylval.ival = 2; return DAMPING; }

"rubdam"[a-zA-Z]*|"RUBDAM"[a-zA-Z]* { yylval.ival = 3; return DAMPING; }

"modd"[a-zA-Z]*|"MODD"[a-zA-Z]*    { return MODDAMP; }

"ARUBBERMAT"|"arubbermat"            { return ARUBBERMAT; }

"time"|"TIME"                      { return TIME; }

"pita"|"PITA"|"Pita"               { return PITA; }

"mdpita"|"MDPITA"|"MdPita"         { return MDPITA; }

"NoForce"|"noforce"|"NOFORCE"      { return NOFORCE; }

"globalbases"|"GlobalBases"|"GLOBALBASES" { return GLOBALBASES; }

"localbases"|"LocalBases"|"LOCALBASES" { return LOCALBASES; }

"remotecoarse"|"REMOTECOARSE"      { return REMOTECOARSE; }

"ORTHOPROJTOL"|"orthoprojtol"      { return ORTHOPROJTOL; }

"timereversible"|"TIMEREVERSIBLE"  { return TIMEREVERSIBLE; }

"readinitseed"|"READINITSEED"      { return READINITSEED; }

"jumpcvg"|"JUMPCVG"                { return JUMPCVG; }

"jumpoutput"|"JUMPOUTPUT"          { return JUMPOUTPUT; }

"heat"|"HEAT"                      { return HEAT; }

"mech"|"MECH"                      { return MECH; }

"acou"[a-zA-Z]*|"ACOU"[a-zA-Z]*    { return ACOU; }

"flui"|"FLUI"                      { return FLUID; }

"OUTPUT"|"output"                  { return OUTPUT; }

"OUTPUT6"|"output6"                { return OUTPUT6; }

"EIGE"[a-zA-Z]*|"eige"[a-zA-Z]*    { return EIGEN; }
 
"SENSITIVITY"|"sensitivity"        { return SENSITIVITY; }

"READSE"[a-zA-Z]*|"readse"[a-zA-Z]*  { return READSENSITIVITY; }

"THGRLI"|"thgrli"                  { return THICKNESSGROUPLIST; }

"STNOLI"|"stnoli"                  { return STRESSNODELIST; }

"DINOLI"|"dinoli"                  { return DISPNODELIST; }

"addedmass"                        { return ADDEDMASS; }

"SLOS"[a-zA-Z]*|"slos"[a-zA-Z]*    { return SLOSH; }

"SLGR"[a-zA-Z]*|"slgr"[a-zA-Z]*    { return SLGRAV; }

"HEFS"[a-zA-Z]*|"hefs"[a-zA-Z]*    { return HEFSB; }

"HEIN"[a-zA-Z]*|"hein"[a-zA-Z]*    { return HEINTERFACE; }

"subs"[a-zA-Z]*|"subs"[a-zA-Z]*    { return SUBSPACE;  }

"NONL"[a-zA-Z]*|"nonl"[a-zA-Z]*    { return NL; }

"ROMCPUS"|"romcpus"|"ROMcpus"      { return NUMROMCPUS; }

"nlmat"[a-zA-Z]*|"nlmat"[a-zA-Z]*  { return NLMAT; }

"FITA"[a-zA-Z]*|"fita"[a-zA-Z]*    { return FITALG; }

"nlto"[a-zA-Z]*|"NLTO"[a-zA-Z]*    { return NLTOL; }

"maxv"[a-zA-Z]*|"MAXV"[a-zA-Z]*    { return MAXVEC; }

"arcl"[a-zA-Z]*|"ARCL"[a-zA-Z]*    { return ARCLENGTH; }

"newt"[a-zA-Z]*|"NEWT"[a-zA-Z]*|"Newt"[a-zA-Z]* { return NEWTON; }

"dlam"[a-zA-Z]*|"DLAM"[a-zA-Z]*    { return DLAMBDA; }

"lfac"[a-zA-Z]*|"LFAC"[a-zA-Z]*    { return LFACTOR; }

"unsy"[a-zA-Z]*|"UNSY"[a-zA-Z]*    { return UNSYMMETRIC; }

"rebuild"|"REBUILD"                { return REBUILD; }

"back"[a-zA-Z]*|"BACK"[a-zA-Z]*    { yylval.ival = 1; return ALPROC; }

"bise"[a-zA-Z]*|"BISE"[a-zA-Z]*    { yylval.ival = 2; return ALPROC; }

"regu"[a-zA-Z]*|"REGU"[a-zA-Z]*    { yylval.ival = 3; return ALPROC; }

"seca"[a-zA-Z]*|"SECA"[a-zA-Z]*    { yylval.ival = 4; return ALPROC; }

"tang"[a-zA-Z]*|"TANG"[a-zA-Z]*    { return TANGENT; }

"preco"[a-zA-Z]*|"PRECO"[a-zA-Z]*  { return PRECONDITIONER; }

"DIMA"[a-zA-Z]*|"dima"[a-zA-Z]*    { return DIMASS; }

"COMPO"[a-zA-Z]*|"composite"       { return COMPOSITE; }

"COEF"|"coef"                      { return COEF; }

"CFRA"[a-zA-Z]*|"cfra"[a-zA-Z]*    { return CFRAMES; }

"LAYC"|"layc"                      { return LAYC; }

"LAYN"|"layn"                      { return LAYN; }

"LAYD"|"layd"                      { return LAYD; }

"LAYO"|"layo"                      { return LAYO; }

"LAYMAT"|"laymat"                  { return LAYMAT; }

"TRBM"[a-zA-Z]*                    { return TRBM; }

"trbm"[a-zA-Z]*                    { return TRBMlc; }

"GRBM"[a-zA-Z]*|"grbm"[a-zA-Z]*    { return GRBM; }

"HZEM"|"hzem"                      { return HZEM; }

"HZEMF"[a-zA-Z]*|"hzemf"[a-zA-Z]*  { return HZEMFILTER; }

"SZEM"[a-zA-Z]*|"szem"[a-zA-Z]*    { return SLZEM; }

"SZFI"[a-zA-Z]*|"szfi"[a-zA-Z]*    { return SLZEMFILTER; }

"GRAV"[a-zA-Z]*|"grav"[a-zA-Z]*    { return GRAVITY; }

"MFTT"|"mftt"                      { return MFTT; }

"HFTT"|"hftt"                      { return HFTT; }

"YMTT"|"ymtt"                      { return YMTT; }

"SDETAFT"|"sdetaft"                { return SDETAFT; }

"RUBDAFT"|"rubdaft"                { return RUBDAFT; }

"MATLAB"|"matlab"                  { return MATLAB; }

"curve"|"CURVE"                    { return CURVE; }

"table"|"TABLE"                    { return TABLE; }

"TETT"|"tett"                      { return TETT; }

"SS1DT"|"ss1dt"                    { return SS1DT; }

"SS2DT"|"ss2dt"                    { return SS2DT; }

"ENGI"[a-zA-Z]*|"engi"[a-zA-Z]*    { return ENGINEERING; }

"YSST"|"ysst"                      { return YSST; }

"YSSRT"|"yssrt"                    { return YSSRT; }

"YMST"|"ymst"                      { return YMST; }

"REOR"[a-zA-Z]*|"reor"[a-zA-Z]*    { return REORTHO; }

"scal"[a-zA-Z]*|"SCAL"[a-zA-Z]*    { return SCALING; }

"mpc_scal"[a-zA-Z]*|"MPC_SCAL"[a-zA-Z]* { return MPCSCALING; }

"fsi_scal"[a-zA-Z]*|"FSI_SCAL"[a-zA-Z]* { return FSISCALING; }

"mpc_ele"[a-zA-Z]*|"MPC_ELE"[a-zA-Z]* { return MPCELEMENT; }

"fsi_ele"[a-zA-Z]*|"FSI_ELE"[a-zA-Z]* { return FSIELEMENT; }

"fsi_cor"[a-zA-Z]*|"FSI_COR"[a-zA-Z]* { return FSICORNER; }


"FSGL"|"fsgl" { return FSGL; }

"coupled_scale"      { return COUPLEDSCALE; }

"wet_corners"        { return WETCORNERS; }

"noLocalFsiSpliting" { return NOLOCALFSISPLITING; }

"stif"[a-zA-Z]*|"STIF"[a-zA-Z]*    { yylval.ival = 1; return SCALINGTYPE; }

"topological"                      { yylval.ival = 2; return SCALINGTYPE; }

"sloan"|"SLOAN"                    { yylval.ival = 1; return RENUMBERID; }

"rcm"|"RCM"                        { yylval.ival = 2; return RENUMBERID; }

"mmd"|"MMD"                        { yylval.ival = 3; return RENUMBERID; }

"nd"|"ND"                          { yylval.ival = 4; return RENUMBERID; }

"ms"|"MS"                          { yylval.ival = 5; return RENUMBERID; }

"ndms"|"NDMS"                      { yylval.ival = 6; return RENUMBERID; }

"esmo"[a-zA-Z]*|"ESMO"[a-zA-Z]*    { yylval.ival = 3; return RENUMBERID; }

"meti"[a-zA-Z]*|"METI"[a-zA-Z]*    { yylval.ival = 4; return RENUMBERID; }

"corn"[a-zA-Z]*|"CORN"[a-zA-Z]*    { return CORNER; }

"sixdof"|"SIXDOF"                  { yylval.ival = 0; return CORNERTYPE; }

"threedof"|"THREEDOF"              { yylval.ival = 1; return CORNERTYPE; }

"noEnd6"|"noend6"|"NOEND6"         { yylval.ival = 2; return CORNERTYPE; }

"noEnd3"|"noend3"|"NOEND3"         { yylval.ival = 3; return CORNERTYPE; }

"be6"|"BE6"                        { yylval.ival = 0; return CORNERTYPE; }

"be3"|"BE3"                        { yylval.ival = 1; return CORNERTYPE; }

"cp6"|"CP6"                        { yylval.ival = 2; return CORNERTYPE; }

"cp3"|"CP3"                        { yylval.ival = 3; return CORNERTYPE; }

"augment"|"AUGMENT"                { return AUGMENT; }

"none"|"NONE"                      { yylval.ival = 0; return AUGMENTTYPE;}

"Gs"|"gs"|"GS"                     { yylval.ival = 1; return AUGMENTTYPE;}

"EdgeGs"|"edgegs"|"EDGEGS"         { yylval.ival = 2; return AUGMENTTYPE;}

"WeightedEdgeGs"|"weightededgegs"|"WEIGHTEDEDGEGS" { yylval.ival = 3; return AUGMENTTYPE;}

"EdgeWs"|"edgews"|"EDGEWS"         { return EDGEWS; }

"solid"|"SOLID"                    { yylval.ival = 0; return WAVETYPE;}

"shell"|"SHELL"                    { yylval.ival = 1; return WAVETYPE;}

"fluid"|"FLUID"                    { yylval.ival = 2; return WAVETYPE;}

"AMAT"|"amat"                      { return AMAT;}

"FABMAT"|"fabmat"                  { return FABMAT;}

"THRMAT"|"thrmat"                  { return THERMMAT; }

"CONMAT"|"conmat"                  { return CONSTRMAT; }

"SPRINGMAT"|"springmat"            { return SPRINGMAT; }

"FREEPLAY"|"freeplay"|"FreePlay"   { return FREEPLAY; }

"SINE"[a-zA-Z]*|"sine"[a-zA-Z]*|"Sine"[a-zA-Z]*  { yylval.ival = 0; return ELEMENTARYFUNCTIONTYPE; }
"RAMP"[a-zA-Z]*|"ramp"[a-zA-Z]*|"Ramp"[a-zA-Z]*  { yylval.ival = 1; return ELEMENTARYFUNCTIONTYPE; }
"TRIA"[a-zA-Z]*|"tria"[a-zA-Z]*|"Tria"[a-zA-Z]*  { yylval.ival = 2; return ELEMENTARYFUNCTIONTYPE; }
"TRAP"[a-zA-Z]*|"trap"[a-zA-Z]*|"Trap"[a-zA-Z]*  { yylval.ival = 3; return ELEMENTARYFUNCTIONTYPE; }
"SSHA"[a-zA-Z]*|"ssha"[a-zA-Z]*|"Ssha"[a-zA-Z]*  { yylval.ival = 4; return ELEMENTARYFUNCTIONTYPE; }

"any"|"ANY"                        { yylval.ival = 3; return WAVETYPE; }

"averageK"|"averagek"|"AVERAGEK"   { yylval.ival = 0; return WAVEMETHOD;}

"averageMat"|"averagemat"|"AVERAGEMAT" { yylval.ival = 1; return WAVEMETHOD;}

"uniform"|"UNIFORM"                { yylval.ival = 2; return WAVEMETHOD;}

"orthotol"|"ORTHOTOL"              { return ORTHOTOL; } 

"trans"[a-zA-Z]*|"TRANS"[a-zA-Z]*  { yylval.ival = FetiInfo::translation; return RBMSET; }

"rot"|"ROT"|"rota"|"ROTA"|"rotat"|"ROTAT"|"rotati"|"ROTATI"|"rotatio"|"ROTATIO"|"rotation"|"ROTATION" { yylval.ival = FetiInfo::rotation; return RBMSET; }

"all"[a-zA-Z]*|"ALL"[a-zA-Z]*      { yylval.ival = FetiInfo::all;         return RBMSET; }

"averageTrans"[a-zA-Z]*|"averagetrans"[a-zA-Z]*|"AVERAGETRANS"[a-zA-Z]* { yylval.ival = FetiInfo::averageTran; return RBMSET; }

"averageRot"[a-zA-Z]*|"averagerot"[a-zA-Z]*|"AVERAGEROT"[a-zA-Z]* { yylval.ival = FetiInfo::averageRot;  return RBMSET; }

"averageAll"[a-zA-Z]*|"averageall"[a-zA-Z]*|"AVERAGEALL"[a-zA-Z]* { yylval.ival = FetiInfo::averageAll;  return RBMSET; }

"noprec"[a-zA-Z]*|"NOPREC"[a-zA-Z]* { yylval.ival = FetiInfo::noPrec;    return FETIPREC; }

"lump"[a-zA-Z]*|"LUMP"[a-zA-Z]*    { return LUMPED; }

"maug"[a-zA-Z]*|"MAUG"[a-zA-Z]*    { return MASSAUGMENTATION; }

"diri"[a-zA-Z]*|"DIRI"[a-zA-Z]*    { yylval.ival = FetiInfo::dirichlet; return FETIPREC; }

"idprec"[a-zA-Z]*|"IDPREC"[a-zA-Z]*  { yylval.ival = FetiInfo::identity; return FETIPREC; }

"mrhs"|"MRHS"                      { return MRHS; }

"AERO"|"aero"                      { return AERO; }

"AEROH"|"aeroh"                    { return AEROH; }

"THERMOH"|"thermoh"                { return THERMOH; }

"THERMOE"|"thermoe"                { return THERMOE; }

"COLL"[a-zA-Z]*|"coll"[a-zA-Z]*    { yylval.ival = 1; return COLLOCATEDTYPE; }

"NONC"[a-zA-Z]*|"nonc"[a-zA-Z]*    { yylval.ival = 0; return COLLOCATEDTYPE; }

"A0"|"a0"                          { yylval.ival = 0; return AEROTYPE; }

"PP"|"pp"                          { yylval.ival = 1; return AEROTYPE; }

"A4"|"a4"                          { yylval.ival = 4; return AEROTYPE; }

"A5"|"a5"                          { yylval.ival = 5; return AEROTYPE; }

"A6"|"a6"                          { yylval.ival = 6; return AEROTYPE; }

"A7"|"a7"                          { yylval.ival = 7; return AEROTYPE; }

"MPP"|"mpp"                        { yylval.ival = 8; return AEROTYPE; }

"B0"|"b0"                          { yylval.ival = 10; return AEROTYPE; }

"C0_legacy"|"c0_legacy"            { yylval.ival = 20; return AEROTYPE; }

"C0"|"c0"                          { yylval.ival = 22; return AEROTYPE; }

"EMBE"[a-zA-Z]*|"embe"[a-zA-z]*    { return AEROEMBED; }

"SUBC"[a-zA-Z]*|"subc"[a-zA-z]*    { return SUBCYCLE; }

"method"|"METHOD"                  { return SENSITIVITYMETHOD; }

"Direct"                           { yylval.ival = SolverInfo::Direct;  return SENSITIVITYID; }

"Adjoint"                          { yylval.ival = SolverInfo::Adjoint; return SENSITIVITYID; }

"displacx"|"DISPLACX"              { yylval.ival = OutputInfo::DispX;  return STRESSID; }

"displacy"|"DISPLACY"              { yylval.ival = OutputInfo::DispY;  return STRESSID; }

"displacz"|"DISPLACZ"              { yylval.ival = OutputInfo::DispZ;  return STRESSID; }

"rotatiox"|"ROTATIOX"              { yylval.ival = OutputInfo::RotX;   return STRESSID; }

"rotatioy"|"ROTATIOY"              { yylval.ival = OutputInfo::RotY;   return STRESSID; }

"rotatioz"|"ROTATIOZ"              { yylval.ival = OutputInfo::RotZ;   return STRESSID; }

"displmod"|"DISPLMOD"              { yylval.ival = OutputInfo::DispMod;return STRESSID; }

"rotatmod"|"ROTATMOD"              { yylval.ival = OutputInfo::RotMod; return STRESSID; }

"gdispmod"|"GDISPMOD"              { yylval.ival = OutputInfo::TotMod; return STRESSID; }

"weigthic"|"WEIGTHIC"              { yylval.ival = OutputInfo::WeigThic; return STRESSID; }

"weigshap"|"WEIGSHAP"              { yylval.ival = OutputInfo::WeigShap; return STRESSID; }

"agstthic"|"AGSTTHIC"              { yylval.ival = OutputInfo::AGstThic; return STRESSID; }

"agstshap"|"AGSTSHAP"              { yylval.ival = OutputInfo::AGstShap; return STRESSID; }

"vmstthic"|"VMSTTHIC"              { yylval.ival = OutputInfo::VMstThic; return STRESSID; }

"vmstshap"|"VMSTSHAP"              { yylval.ival = OutputInfo::VMstShap; return STRESSID; }

"vmstmach"|"VMSTMACH"              { yylval.ival = OutputInfo::VMstMach; return STRESSID; }

"vmstalph"|"VMSTALPH"              { yylval.ival = OutputInfo::VMstAlpha; return STRESSID; }

"vmstbeta"|"VMSTBETA"              { yylval.ival = OutputInfo::VMstBeta; return STRESSID; }

"dispthic"|"DISPTHIC"              { yylval.ival = OutputInfo::DispThic; return STRESSID; }

"dispshap"|"DISPSHAP"              { yylval.ival = OutputInfo::DispShap; return STRESSID; }

"dispmach"|"DISPMACH"              { yylval.ival = OutputInfo::DispMach; return STRESSID; }

"dispalph"|"DISPALPH"              { yylval.ival = OutputInfo::DispAlph; return STRESSID; }

"dispbeta"|"DISPBTEA"              { yylval.ival = OutputInfo::DispBeta; return STRESSID; }

"elemental"[a-zA-Z]*|"ELEMENTAL"[a-zA-Z]* { yylval.ival =  0; return AVERAGED; }

"gauss"[a-zA-Z]*|"GAUSS"[a-zA-Z]*         { yylval.ival = -1; return AVERAGED; }

"nodal"|"NODAL"|"nodalfull"|"NODALFULL"   { yylval.ival =  1; return AVERAGED; }

"nodalpartial"|"NODALPARTIAL"             { yylval.ival =  2; return AVERAGED; }

"nodalpartialgroup"|"NODALPARTIALGROUP"   { yylval.ival =  3; return AVERAGED; }

"upper"[a-zA-Z]*|"UPPER"[a-zA-Z]*  { yylval.ival = 1; return SURFACE; }

"median"[a-zA-Z]*|"MEDIAN"[a-zA-Z]* { yylval.ival = 2; return SURFACE; }

"lower"[a-zA-Z]*|"LOWER"[a-zA-Z]*  { yylval.ival = 3; return SURFACE; }

"thermomechanical"[a-zA-Z]*|"THERMOMECHANICAL"[a-zA-Z]*  { yylval.ival = 0; return STR_THERM_OPTION; }

"thermal"[a-zA-Z]*|"THERMAL"[a-zA-Z]*  { yylval.ival = 1; return STR_THERM_OPTION; }

"mechanical"[a-zA-Z]*|"MECHANICAL"[a-zA-Z]* { yylval.ival = 2; return STR_THERM_OPTION; }

"gdisplac"|"GDISPLAC"              { yylval.ival = OutputInfo::Displacement; return STRESSID; }

"romatrix"|"ROMATRIX"              { yylval.ival = OutputInfo::RotationMatrix; return STRESSID; }

"rovector"|"ROVECTOR"              { yylval.ival = OutputInfo::RotationVector; return STRESSID; }

"angveloc"|"ANGVELOC"              { yylval.ival = OutputInfo::AngularVelocity; return STRESSID; }

"angaccel"|"ANGACCEL"              { yylval.ival = OutputInfo::AngularAcceleration; return STRESSID; }

"quaterni"|"QUATERNI"              { yylval.ival = OutputInfo::Quaternion;   return STRESSID; }

"modaldsp"|"MODALDSP"|"sromdisp"|"SROMDISP" { yylval.ival = OutputInfo::ModalDsp;     return STRESSID; }

"modalexf"|"MODALEXF"|"sromextf"|"SROMEXTF" { yylval.ival = OutputInfo::ModalExF;     return STRESSID; }

"modalmas"|"MODALMAS"|"srommass"|"SROMMASS" { yylval.ival = OutputInfo::ModalMass;    return STRESSID; }

"modalstf"|"MODALSTF"|"sromstif"|"SROMSTIF" { yylval.ival = OutputInfo::ModalStiffness; return STRESSID; }

"modaldmp"|"MODALDMP"|"sromdamp"|"SROMDAMP" { yylval.ival = OutputInfo::ModalDamping; return STRESSID; }

"modaldyn"|"MODALDYN"|"sromdyna"|"SROMDYNA" { yylval.ival = OutputInfo::ModalDynamicMatrix; return STRESSID; }

"modalmat"|"MODALMAT"|"srommatr"|"SROMMATR" { yylval.ival = OutputInfo::ModalMatrices; return STRESSID; }

"errindic"|"ERRINDIC"              { yylval.ival = OutputInfo::ErrorIndicator; return STRESSID; }

"gtempera"|"GTEMPERA"              { yylval.ival = OutputInfo::Temperature;  return STRESSID; }

"statevct"|"STATEVCT"		   { yylval.ival = OutputInfo::Statevector;  return STRESSID; }

"velocvct"|"VELOCVCT"              { yylval.ival = OutputInfo::Velocvector;  return STRESSID; }

"accelvct"|"ACCELVCT"		   { yylval.ival = OutputInfo::Accelvector;  return STRESSID; }

"frcgvect"|"FRCGVECT"		   { yylval.ival = OutputInfo::Forcevector;  return STRESSID; } /* deprecated */
"frcngvct"|"FRCNGVCT"              { yylval.ival = OutputInfo::Forcevector;  return STRESSID; }

"constvct"|"CONSTVCT"		   { yylval.ival = OutputInfo::Constraintvector;  return STRESSID; }

"cviolate"|"CVIOLATE"              { yylval.ival = OutputInfo::Constraintviolation;  return STRESSID; }

"residvct"|"RESIDVCT"		   { yylval.ival = OutputInfo::Residual;     return STRESSID; } /* deprecated */
"fdiffvct"|"FDIFFVCT"              { yylval.ival = OutputInfo::Residual;     return STRESSID; }

"jacobvct"|"JACOBVCT"		   { yylval.ival = OutputInfo::Jacobian;     return STRESSID; }

"isvarvct"|"ISVARVCT"              { yylval.ival = OutputInfo::InternalStateVar; return STRESSID; }

"dsvarvct"|"DSVARVCT"              { yylval.ival = OutputInfo::DualStateVar; return STRESSID; }

"muvarvct"|"MUVARVCT"              { yylval.ival = OutputInfo::MuStateVar; return STRESSID; }

"robdataf"|"ROBDATAF"		   { yylval.ival = OutputInfo::RobData;	     return STRESSID; }
"clusterf"|"CLUSTERF"              { yylval.ival = OutputInfo::RobData;      return STRESSID; }

"samplmsh"|"SAMPLMSH" 		   { yylval.ival = OutputInfo::SampleMesh;   return STRESSID; }

"stressxx"|"STRESSXX"              { yylval.ival = OutputInfo::StressXX;     return STRESSID; }

"stressyy"|"STRESSYY"              { yylval.ival = OutputInfo::StressYY;     return STRESSID; }

"stresszz"|"STRESSZZ"              { yylval.ival = OutputInfo::StressZZ;     return STRESSID; }

"stressxy"|"STRESSXY"              { yylval.ival = OutputInfo::StressXY;     return STRESSID; }

"stressxz"|"STRESSXZ"              { yylval.ival = OutputInfo::StressXZ;     return STRESSID; }

"stresszx"|"STRESSZX"              { yylval.ival = OutputInfo::StressXZ;     return STRESSID; }

"stressyz"|"STRESSYZ"              { yylval.ival = OutputInfo::StressYZ;     return STRESSID; }

"strainxx"|"STRAINXX"              { yylval.ival = OutputInfo::StrainXX;     return STRESSID; }

"strainyy"|"STRAINYY"              { yylval.ival = OutputInfo::StrainYY;     return STRESSID; }

"strainzz"|"STRAINZZ"              { yylval.ival = OutputInfo::StrainZZ;     return STRESSID; }

"strainxy"|"STRAINXY"              { yylval.ival = OutputInfo::StrainXY;     return STRESSID; }

"strainzx"|"STRAINZX"              { yylval.ival = OutputInfo::StrainXZ;     return STRESSID; }

"strainxz"|"STRAINXZ"              { yylval.ival = OutputInfo::StrainXZ;     return STRESSID; }

"strainyz"|"STRAINYZ"              { yylval.ival = OutputInfo::StrainYZ;     return STRESSID; }

"heatflxx"|"HEATFLXX"              { yylval.ival = OutputInfo::HeatFlXX;     return STRESSID; }

"heatflxy"|"HEATFLXY"              { yylval.ival = OutputInfo::HeatFlXY;     return STRESSID; }

"heatflxz"|"HEATFLXZ"              { yylval.ival = OutputInfo::HeatFlXZ;     return STRESSID; }

"grdtempx"|"GRDTEMPX"              { yylval.ival = OutputInfo::GrdTempX;     return STRESSID; }

"grdtempy"|"GRDTEMPY"              { yylval.ival = OutputInfo::GrdTempY;     return STRESSID; }

"grdtempz"|"GRDTEMPZ"              { yylval.ival = OutputInfo::GrdTempZ;     return STRESSID; }

"slshdisp"|"SLSHDISP"              { yylval.ival = OutputInfo::SloshDisplacement;     return STRESSID; }

"slshdspx"|"SLSHDSPX"              { yylval.ival = OutputInfo::SloshDispX;     return STRESSID; }

"slshdspy"|"SLSHDSPY"              { yylval.ival = OutputInfo::SloshDispY;     return STRESSID; }

"slshdspz"|"SLSHDSPZ"              { yylval.ival = OutputInfo::SloshDispZ;     return STRESSID; }

"stressvm"|"STRESSVM"              { yylval.ival = OutputInfo::StressVM;     return STRESSID; }

"aggrstvm"|"AGGRSTVM"              { yylval.ival = OutputInfo::AggrStVM;     return STRESSID; }

"stressp1"|"STRESSP1"              { yylval.ival = OutputInfo::StressPR1;    return STRESSID; }

"stressp2"|"STRESSP2"              { yylval.ival = OutputInfo::StressPR2;    return STRESSID; }

"stressp3"|"STRESSP3"              { yylval.ival = OutputInfo::StressPR3;    return STRESSID; }

"strainp1"|"STRAINP1"              { yylval.ival = OutputInfo::StrainPR1;    return STRESSID; }

"strainp2"|"STRAINP2"              { yylval.ival = OutputInfo::StrainPR2;    return STRESSID; }

"strainp3"|"STRAINP3"              { yylval.ival = OutputInfo::StrainPR3;    return STRESSID; }

"inxforce"|"INXFORCE"              { yylval.ival = OutputInfo::InXForce;     return STRESSID; }

"inyforce"|"INYFORCE"              { yylval.ival = OutputInfo::InYForce;     return STRESSID; }

"inzforce"|"INZFORCE"              { yylval.ival = OutputInfo::InZForce;     return STRESSID; }

"axmoment"|"AXMOMENT"              { yylval.ival = OutputInfo::AXMoment;     return STRESSID; }

"aymoment"|"AYMOMENT"              { yylval.ival = OutputInfo::AYMoment;     return STRESSID; }

"azmoment"|"AZMOMENT"              { yylval.ival = OutputInfo::AZMoment;     return STRESSID; }

"energies"|"ENERGIES"              { yylval.ival = OutputInfo::Energies;     return STRESSID; }

"raerofor"|"RAEROFOR"              { yylval.ival = OutputInfo::AeroForce;    return STRESSID; }

"romresid"|"ROMRESID"              { yylval.ival = OutputInfo::RomResidual;  return STRESSID; }

"romextfo"|"ROMEXTFO"              { yylval.ival = OutputInfo::RomExtForce;  return STRESSID; }

"geigenpa"|"GEIGENPA"              { yylval.ival = OutputInfo::EigenPair;    return STRESSID; }

"strainvm"|"STRAINVM"              { yylval.ival = OutputInfo::StrainVM;     return STRESSID; }

"ghelmhol"|"gacoustd"|"pressure"   { yylval.ival = OutputInfo::Helmholtz;    return STRESSID; }

"GHELMHOL"|"GACOUSTD"|"PRESSURE"   { yylval.ival = OutputInfo::Helmholtz;    return STRESSID; }

"GACOUPRE"|"gacoupre"              { yylval.ival = OutputInfo::AcousticPressure; return STRESSID; }

"aeroforx"|"AEROFORX"              { yylval.ival = OutputInfo::AeroXForce;   return STRESSID; }

"aerofory"|"AEROFORY"              { yylval.ival = OutputInfo::AeroYForce;   return STRESSID; }

"aeroforz"|"AEROFORZ"              { yylval.ival = OutputInfo::AeroZForce;   return STRESSID; }

"aeromomx"|"AEROMOMX"              { yylval.ival = OutputInfo::AeroXMom;     return STRESSID; }

"aeromomy"|"AEROMOMY"              { yylval.ival = OutputInfo::AeroYMom;     return STRESSID; }

"aeromomz"|"AEROMOMZ"              { yylval.ival = OutputInfo::AeroZMom;     return STRESSID; }

"exteforx"|"EXTEFORX"              { yylval.ival = OutputInfo::ExternalXForce; return STRESSID; }

"extefory"|"EXTEFORY"              { yylval.ival = OutputInfo::ExternalYForce; return STRESSID; }

"exteforz"|"EXTEFORZ"              { yylval.ival = OutputInfo::ExternalZForce; return STRESSID; }

"extemomx"|"EXTEMOMX"              { yylval.ival = OutputInfo::ExternalXMom;  return STRESSID; }

"extemomy"|"EXTEMOMY"              { yylval.ival = OutputInfo::ExternalYMom;  return STRESSID; }

"extemomz"|"EXTEMOMZ"              { yylval.ival = OutputInfo::ExternalZMom;  return STRESSID; }

"gvelocit"|"GVELOCIT"              { yylval.ival = OutputInfo::Velocity;     return STRESSID; }

"gacceler"|"GACCELER"              { yylval.ival = OutputInfo::Acceleration; return STRESSID; }

"gtempvel"|"GTEMPVEL"              { yylval.ival = OutputInfo::TemperatureFirstTimeDerivative; return STRESSID; }

"gacouvel"|"GACOUVEL"              { yylval.ival = OutputInfo::PressureFirstTimeDerivative; return STRESSID; }

"gacouacc"|"GACOUACC"              { yylval.ival = OutputInfo::PressureSecondTimeDerivative; return STRESSID; }

"ymodulus"|"YMODULUS"              { yylval.ival = OutputInfo::YModulus;     return STRESSID; }

"mdensity"|"MDENSITY"              { yylval.ival = OutputInfo::MDensity;     return STRESSID; }

"thicknes"|"THICKNES"              { yylval.ival = OutputInfo::Thicknes;     return STRESSID; }

"shapeatt"|"SHAPEATT"              { yylval.ival = OutputInfo::ShapeAtt;     return STRESSID; }

"shapestc"|"SHAPESTC"              { yylval.ival = OutputInfo::ShapeStc;     return STRESSID; }

"rigidbod"|"RIGIDBOD"              { yylval.ival = OutputInfo::Rigid;        return STRESSID; }

"composit"|"COMPOSIT"              { yylval.ival = OutputInfo::Composit;     return STRESSID; }

"eletonod"|"ELASTONOD"             { yylval.ival = OutputInfo::ElemToNode;   return STRESSID; }

"nodtoele"|"NODTOELE"              { yylval.ival = OutputInfo::NodeToElem;   return STRESSID; }

"nodtonod"|"NODTONOD"              { yylval.ival = OutputInfo::NodeToNode;   return STRESSID; }

"raerotfl"|"RAEROTFL"              { yylval.ival = OutputInfo::AeroFlux;     return STRESSID; }

"heatflx"|"HEATFLX"                { yylval.ival = OutputInfo::HeatFlX;       return STRESSID; }

"grdtemp"|"GRDTEMP"                { yylval.ival = OutputInfo::GrdTemp;       return STRESSID; }

"gencoord"|"GENCOORD"              { yylval.ival = OutputInfo::ModeAlpha;    return STRESSID; }

"moderror"|"MODERROR"              { yylval.ival = OutputInfo::ModeError;    return STRESSID; }

"reaction"|"REACTION"              { yylval.ival = OutputInfo::Reactions;    return STRESSID; }

"heatreac"|"HEATREAC"              { yylval.ival = OutputInfo::HeatReactions; return STRESSID; }

"conpress"|"CONPRESS"              { yylval.ival = OutputInfo::ContactPressure; return STRESSID; }

"farfield"|"FARFIELD"              { yylval.ival = OutputInfo::Farfield;     return STRESSID; }

"kirchhoff"|"KIRCHHOFF"            { yylval.ival = OutputInfo::Kirchhoff;    return STRESSID; }

"geigpres"|"GEIGPRES"              { yylval.ival = OutputInfo::EigenPressure; return STRESSID; }

"geigslsh"|"GEIGSLSH"              { yylval.ival = OutputInfo::EigenSlosh; return STRESSID; }

"gfrmodes"|"GFRMODES"              { yylval.ival = OutputInfo::FreqRespModes; return STRESSID; }

"ghhmodes"|"GHHMODES"              { yylval.ival = OutputInfo::HelmholtzModes; return STRESSID; }

"sp1direc"|"SP1DIREC"              { yylval.ival = OutputInfo::StressPR1Direc; return STRESSID; }

"sp2direc"|"SP2DIREC"              { yylval.ival = OutputInfo::StressPR2Direc; return STRESSID; }

"sp3direc"|"SP3DIREC"              { yylval.ival = OutputInfo::StressPR3Direc; return STRESSID; }

"ep1direc"|"EP1DIREC"              { yylval.ival = OutputInfo::StrainPR1Direc; return STRESSID; }

"ep2direc"|"EP2DIREC"              { yylval.ival = OutputInfo::StrainPR2Direc; return STRESSID; }

"ep3direc"|"EP3DIREC"              { yylval.ival = OutputInfo::StrainPR3Direc; return STRESSID; }

"realimag"|"REALIMAG"              { yylval.ival = OutputInfo::realimag; return COMPLEXOUTTYPE; }

"modphase"|"MODPHASE"              { yylval.ival = OutputInfo::modulusphase; return COMPLEXOUTTYPE; }

"anim"[a-zA-Z]*|"ANIM"[a-zA-Z]*    { yylval.ival = OutputInfo::animate; return COMPLEXOUTTYPE; }

"spatial"|"Spatial"|"SPATIAL"      { yylval.ival = OutputInfo::spatial; return ANGULAROUTTYPE; }

"convected"|"Convected"|"CONVECTED" { yylval.ival = OutputInfo::convected; return ANGULAROUTTYPE; }

"total"|"Total"|"TOTAL"            { yylval.ival = OutputInfo::total; return ANGULAROUTTYPE; }

"euler"|"Euler"|"EULER"            { yylval.ival = OutputInfo::Euler; return ROTVECOUTTYPE; }

"complement"|"Complement"|"COMPLEMENT" { yylval.ival = OutputInfo::Complement; return ROTVECOUTTYPE; }

"ReducedEulerRodrigues"            { yylval.ival = OutputInfo::ReducedEulerRodrigues; return ROTVECOUTTYPE; }

"CayleyGibbsRodrigues"             { yylval.ival = OutputInfo::CayleyGibbsRodrigues; return ROTVECOUTTYPE; }

"WienerMilenkovic"                 { yylval.ival = OutputInfo::WienerMilenkovic; return ROTVECOUTTYPE; }

"BauchauTrainelli"                 { yylval.ival = OutputInfo::BauchauTrainelli; return ROTVECOUTTYPE; }

"rescal"[a-zA-Z]*|"RESCAL"[a-zA-Z]* { return RESCALING; }

"oeigenvalue"|"OEIGENVALUE"          { return EIGENVALUE; }

"xmatrix"|"XMATRIX"                { return XMATRIX; }

"qmatrix"|"QMATRIX"                { return QMATRIX; }

"rmatrix"|"RMATRIX"                { return RMATRIX; }

"conface"|"CONFACE"                                  { yylval.ival = 1; return TDENFORC; }
"normal_force_mag"|"NORMAL_FORCE_MAG"                { yylval.ival = 2; return TDENFORC; }
"normal_traction_mag"|"NORMAL_TRACTION_MAG"          { yylval.ival = 3; return TDENFORC; }
"tangential_force_mag"|"TANGENTIAL_FORCE_MAG"        { yylval.ival = 4; return TDENFORC; }
"tangential_traction_mag"|"TANGENTIAL_TRACTION_MAG"  { yylval.ival = 5; return TDENFORC; }
"cdirnorx"|"CDIRNORX"                                { yylval.ival = 6; return TDENFORC; }
"cdirnory"|"CDIRNORY"                                { yylval.ival = 7; return TDENFORC; }
"cdirnorz"|"CDIRNORZ"                                { yylval.ival = 8; return TDENFORC; }
"cdirtanx"|"CDIRTANX"                                { yylval.ival = 9; return TDENFORC; }
"cdirtany"|"CDIRTANY"                                { yylval.ival = 10; return TDENFORC; }
"cdirtanz"|"CDIRTANZ"                                { yylval.ival = 11; return TDENFORC; }
"slip_mag"|"SLIP_MAG"                                { yylval.ival = 12; return TDENFORC; }
"nodal_dissipation"|"NODAL_DISSIPATION"              { yylval.ival = 13; return TDENFORC; }
"contact_area"|"CONTACT_AREA"                        { yylval.ival = 14; return TDENFORC; }
"gap_cur"|"GAP_CUR"                                  { yylval.ival = 15; return TDENFORC; }
"gap_old"|"GAP_OLD"                                  { yylval.ival = 16; return TDENFORC; }

"elemdele"|"ELEMDELE"              { yylval.ival = OutputInfo::DeletedElements; return STRESSID; }

"damagesc"|"DAMAGESC"              { yylval.ival = OutputInfo::Damage; return STRESSID; }

"effpstrn"|"EFFPSTRN"              { yylval.ival = OutputInfo::EquivalentPlasticStrain; return STRESSID; }

"plstrnxx"|"PLSTRNXX"              { yylval.ival = OutputInfo::PlasticStrainXX; return STRESSID; }

"plstrnyy"|"PLSTRNYY"              { yylval.ival = OutputInfo::PlasticStrainYY; return STRESSID; }

"plstrnzz"|"PLSTRNZZ"              { yylval.ival = OutputInfo::PlasticStrainZZ; return STRESSID; }

"plstrnxy"|"PLSTRNXY"              { yylval.ival = OutputInfo::PlasticStrainXY; return STRESSID; }

"plstrnzx"|"PLSTRNZX"              { yylval.ival = OutputInfo::PlasticStrainXZ; return STRESSID; }

"plstrnxz"|"PLSTRNXZ"              { yylval.ival = OutputInfo::PlasticStrainXZ; return STRESSID; }

"plstrnyz"|"PLSTRNYZ"              { yylval.ival = OutputInfo::PlasticStrainYZ; return STRESSID; }

"bkstrsxx"|"BKSTRSXX"              { yylval.ival = OutputInfo::BackStressXX;    return STRESSID; }

"bkstrsyy"|"BKSTRSYY"              { yylval.ival = OutputInfo::BackStressYY;    return STRESSID; }

"bkstrszz"|"BKSTRSZZ"              { yylval.ival = OutputInfo::BackStressZZ;    return STRESSID; }

"bkstrsxy"|"BKSTRSXY"              { yylval.ival = OutputInfo::BackStressXY;    return STRESSID; }

"bkstrszx"|"BKSTRSZX"              { yylval.ival = OutputInfo::BackStressXZ;    return STRESSID; }

"bkstrsxz"|"BKSTRSXZ"              { yylval.ival = OutputInfo::BackStressXZ;    return STRESSID; }

"bkstrsyz"|"BKSTRSYZ"              { yylval.ival = OutputInfo::BackStressYZ;    return STRESSID; }

"dissipat"|"DISSIPAT"              { yylval.ival = OutputInfo::DissipatedEnergy;return STRESSID; }

"MATSPEC"|"matspec"                { return MATSPEC; }

"MATLAW"|"matlaw"                  { return MATSPEC; }

"MATUSAGE"|"matusage"              { return MATUSAGE; }

"bili"[a-zA-Z]*|"BILI"[a-zA-Z]*|"Bili"[a-zA-Z]*         { return BILINEARPLASTIC; }

"fini"[a-zA-Z]*|"FINI"[a-zA-Z]*|"Fini"[a-zA-Z]*         { return FINITESTRAINPLASTIC; }

"crush"[a-zA-Z]*|"CRUSH"[a-zA-Z]*|"Crush"[a-zA-Z]*         { return CRUSHABLEFOAM; }

"logs"[a-zA-Z]*|"LOGS"[a-zA-Z]*|"LogS"[a-zA-Z]*         { return LOGSTRAINPLASTIC; }


"linearplanestress"|"LINEARPLANESTRESS"|"LinearPlaneStress" { return LINPLSTRESS; }

"viscolinearplanestress"|"VISCOLINEARPLANESTRESS"|"ViscoLinearPlaneStress" { return VISCOLINPLSTRESS; }

"hyperelasticplanestress"|"HYPERELASTICPLANESTRESS"|"HyperElasticPlaneStress" { return SVKPLSTRESS; }

"viscohyperelasticplanestress"|"VISCOHYPERELASTICPLANESTRESS"|"ViscoHyperElasticPlaneStress" { return VISCOSVKPLSTRESS; }

"linear"|"LINEAR"|"Linear"|"linearelastic"|"LINEARELASTIC"|"LinearElastic" { return LINEARELASTIC; }

"isotropiclinearelastic"|"ISOTROPICLINEARELASTIC"|"IsotropicLinearElastic"  { return ISOTROPICLINEARELASTIC; }

"viscolinearelastic"|"VISCOLINEARELASTIC"|"ViscoLinearElastic"|"ViscoLinear"  { return VISCOLINEARELASTIC; }

"viscostve"[a-zA-Z]*|"VISCOSTVE"[a-zA-Z]*|"ViscoStVe"[a-zA-Z]*  { return VISCOSTVENANTKIRCHHOFF; }

"PlaneStressLinear"[a-zA-Z]*                            { return PLANESTRESSLINEAR; }

"PlaneStressStVenantKirchhoff"[a-zA-Z]*                 { return PLANESTRESSSTVENANTKIRCHHOFF; }

"PlaneStressNeoHookean"[a-zA-Z]*                        { return PLANESTRESSNEOHOOKEAN; }

"PlaneStressMooneyRivlin"[a-zA-Z]*                      { return PLANESTRESSMOONEYRIVLIN; }

"PlaneStressBilinearPlastic"[a-zA-Z]*                   { return PLANESTRESSBILINEARPLASTIC; }

"PlaneStressFiniteStrainPlastic"[a-zA-Z]*               { return PLANESTRESSFINITESTRAINPLASTIC; }

"PlaneStressViscoLinearElastic"[a-zA-Z]*                { return PLANESTRESSVISCOLINEARELASTIC; }

"PlaneStressViscoStVenantKirchhoff"[a-zA-Z]*            { return PLANESTRESSVISCOSTVENANTKIRCHHOFF; }

"PlaneStressViscoNeoHookean"[a-zA-Z]*                   { return PLANESTRESSVISCONEOHOOKEAN; }

"PlaneStressViscoMooneyRivlin"[a-zA-Z]*                 { return PLANESTRESSVISCOMOONEYRIVLIN; }

"fabricmap"|"FABRICMAP"|"FabricMap"                     { return FABRICMAP; }

"shellfabricmap"|"SHELLFABRICMAP"|"ShellFabricMap"      { return SHELLFABRICMAP; }

"fabricmat"|"FABRICMAT"|"FabricMat"                     { return FABRICMAT; }

"shellfabricmat"|"SHELLFABRICMAT"|"ShellFabricMat"      { return SHELLFABRICMAT; }

"viscofabricmap"|"VISCOFABRICMAP"|"ViscoFabricMap"                     { return VISCOFABRICMAP; }

"shellviscofabricmap"|"SHELLVISCOFABRICMAP"|"ShellViscoFabricMap"      { return SHELLVISCOFABRICMAP; }

"viscofabricmat"|"VISCOFABRICMAT"|"ViscoFabricMat"                     { return VISCOFABRICMAT; }

"shellviscofabricmat"|"SHELLVISCOFABRICMAT"|"ShellViscoFabricMat"      { return SHELLVISCOFABRICMAT; }

"stve"[a-zA-Z]*|"STVE"[a-zA-Z]*|"StVe"[a-zA-Z]*         { return STVENANTKIRCHHOFF; }

"tule"[a-zA-Z]*|"TULE"[a-zA-Z]*|"Tule"[a-zA-Z]*         { return TULERBUTCHER; }

"henc"[a-zA-Z]*|"HENC"[a-zA-Z]*|"Henc"[a-zA-Z]*         { return HENCKY; }

"neoh"[a-zA-Z]*|"NEOH"[a-zA-Z]*|"NeoH"[a-zA-Z]*         { return NEOHOOKEAN; }

"visconeoh"[a-zA-Z]*|"VISCONEOH"[a-zA-Z]*|"ViscoNeoH"[a-zA-Z]*         { return VISCONEOHOOKEAN; }

"moon"[a-zA-Z]*|"MOON"[a-zA-Z]*|"Moon"[a-zA-Z]*         { return MOONEYRIVLIN; }

"ogde"[a-zA-Z]*|"OGDE"[a-zA-Z]*|"Ogde"[a-zA-Z]*         { return OGDEN; }

"simoe"[a-zA-Z]*|"SIMOE"[a-zA-Z]*|"SimoE"[a-zA-Z]*      { return SIMOELASTIC; }

"simop"[a-zA-Z]*|"SIMOP"[a-zA-Z]*|"SimoP"[a-zA-Z]*      { return SIMOPLASTIC; }

"viscomoon"[a-zA-Z]*|"VISCOMOON"[a-zA-Z]*|"ViscoMoon"[a-zA-Z]*         { return VISCOMOONEYRIVLIN; }

"isotropiclinearelasticj2plastic"|"ISOTROPICLINEARELASTICJ2PLASTIC"|"IsotropicLinearElasticJ2Plastic" { return ISOTROPICLINEARELASTICJ2PLASTIC; }

"isotropiclinearelasticj2plasticplanestress"|"ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS"|"IsotropicLinearElasticJ2PlasticPlaneStress" { return ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS; }

"hypo"[a-zA-Z]*|"HYPO"[a-zA-Z]*|"Hypo"[a-zA-Z]*         { yylval.ival = 1; return OPTCTV; }
"elastoviscoplastic"|"ELASTOVISCOPLASTIC"|"evp"|"EVP"   { yylval.ival = 3; return OPTCTV; }
"j2"[a-zA-Z]*|"J2"[a-zA-Z]*                             { yylval.ival = 5; return OPTCTV; }
"kk"|"KK"|"kk1"|"KK1"                                   { yylval.ival = 6; return OPTCTV; }
"kkexp"|"KKEXP"|"kk2"|"KK2"                             { yylval.ival = 7; return OPTCTV; }
"kkexp2"|"KKEXP2"                                       { yylval.ival = 8; return OPTCTV; }

"read"|"READ"                      { return READ; }

"NODALCONTACT"|"nodalcontact"      { return NODALCONTACT; }

"BMPC"|"bmpc"                      { return BMPC; }
"DMPC"|"dmpc"                      { return DMPC; }
"CMPC"|"cmpc"                      { return CMPC; }

"FRIC"[a-zA-Z]*|"fric"[a-zA-Z]*    { return FRIC; }

"GAP"|"gap"                        { return GAP; }

"SURFTOPO"|"SURFACETOPO"           { return SURFACETOPOLOGY; }

"surftopo"|"surfacetopo"           { return SURFACETOPOLOGY; }

"REVERSE"[a-zA-Z_]*|"reverse"[a-zA-Z_]* { return REVERSENORMALS; }

"SURFACE_THICKNESS"|"surface_thickness" { return SHELLTHICKNESS; }

"MORTARTIED"|"mortartied"          { return MORTARTIED; }

"mortar_scal"[a-zA-Z]*             { return MORTARSCALING; }

"mortar_inte"[a-zA-Z_]*            { return MORTARINTEGRATIONRULE; }

"SEARCH_TOL"|"search_tol"          { return SEARCHTOL; }

"STD"|"std"                        { return STDMORTAR; }

"DUAL"|"dual"                      { return DUALMORTAR; }

"WETINTERFACE"|"wetinterface"      { return WETINTERFACE; }

"FSINTERFACE"|"fsinterface"        { return FSINTERFACE; }

"TIEDSURFACES"|"tiedsurfaces"      { return TIEDSURFACES; }

"CONTACTSURFACES"|"contactsurfaces" { return CONTACTSURFACES; }

"CONSTRAINEDSURFACES"|"constrainedsurfaces" { return CONSTRAINEDSURFACES; }

"POINTPOINT"|"pointpoint"|"PointPoint" { yylval.ival = BCond::PointPointDistance; return CSTYPE; }

"POINTLINE"|"pointline"|"PointLine"    { yylval.ival = BCond::PointLineDistance; return CSTYPE; }

"POINTPLANE"|"pointplane"|"PointPlane" { yylval.ival = BCond::PointPlaneDistance; return CSTYPE; }

"NO_SECONDARY"|"no_secondary"      { return NOSECONDARY; }
"NO_GHOSTING"|"no_ghosting"      { return NOGHOSTING; }
"SHELL_SIMPLE_LOFTING"|"shell_simple_lofting" { return SHELLSIMPLELOFTING; }
"NO_MULTIPLE_INTERACTIONS"|"no_multiple_interactions" { return NOMULTIPLEINTERACTIONS; }
"NO_NORMAL_SMOOTHING"|"no_normal_smoothing" { return NONORMALSMOOTHING; }
"NORMAL_SMOOTHING_DISTANCE"|"normal_smoothing_distance" { return NORMALSMOOTHINGDISTANCE; }
"RESOLUTION_METHOD"|"resolution_method" { return RESOLUTIONMETHOD; }
"OLD_DYNAMIC_SEARCH"|"old_dynamic_search" { return OLDDYNAMICSEARCH; }
"SHARP_NON_SHARP_ANGLE"|"sharp_non_sharp_angle" { return SHARPNONSHARPANGLE; }
"PARTITION_GAP"|"partition_gap" { return PARTITIONGAP; }
"DEFAULT_PENALTY"|"default_penalty" { return DEFAULTPENALTY; }
"GLOBAL_SEARCH_CULL"|"global_search_cull" { return GLOBALSEARCHCULL; }
"NO_WARPED_VOLUME"|"no_warped_volume" { return NOWARPEDVOLUME; }
"AUTO_TOL"|"auto_tol" { return AUTOTOL; }
"AGRESSIVE_TOLERANCES"|"agressive_tolerances" { return AGRESSIVETOLERANCES; }
"SKIP_PHYISCAL_FACES"|"skip_physical_faces" { return SKIPPHYSICALFACES; }

"TDEN"[a-zA-Z]*|"tden"[a-zA-Z]*    { return TDENFORCE; }

"PARA"[a-zA-Z]*|"para"[a-zA-Z]*    { return PARAMETERS; }

"ACME"[a-zA-Z]*|"acme"[a-zA-Z]*    { return ACMECNTL; }

"ANDES_CLR"|"andes_clr"            { return ANDESCLR; }
"ANDES_CQR"|"andes_cqr"            { return ANDESCQR; }
"ANDES_BETAB"|"andes_betab"        { return ANDESBETAB; }
"ANDES_ALPHA"|"andes_alpha"        { return ANDESALPHA; }
"ANDES_BETAM"|"andes_betam"        { return ANDESBETAM; }

"NLMEM_PTYPE"|"nlmem_ptype"        { return NLMEMPTYPE; }

"ffi_debug"                        { return FFIDEBUG; }

"NONINPC" { return NONINPC; }
"INPC" { return INPC; }
"GROU"[a-zA-Z]*|"grou"[a-zA-Z]* { return GROUP; }
"N"|"n" { yylval.ival = OutputInfo::Nodal; return GROUPTYPE; }
"A"|"a" { yylval.ival = OutputInfo::Attribute; return GROUPTYPE; }
"NG"|"ng" { yylval.ival = OutputInfo::NodeGroup; return GROUPTYPE; }
"RAND"[a-zA-Z]*|"rand"[a-zA-Z]* { return RANDOM; }
"E"|"e" { yylval.rprop = GeoSource::E; return RPROP; }
"KX"|"kx"|"Kx" { yylval.rprop = GeoSource::KX; return RPROP; }
"KY"|"ky"|"Ky" { yylval.rprop = GeoSource::KY; return RPROP; }
"KZ"|"kz"|"Kz" { yylval.rprop = GeoSource::KZ; return RPROP; }
"MEAN"|"mean" { yylval.ival = 1; return NDTYPE; }
"STDV"|"stdv" { yylval.ival = 2; return NDTYPE; }
"PDF"|"pdf" { yylval.ival = 3; return NDTYPE; }

"CONFOR"|"ConFor"|"confor"                { return REBUILDFORCE; }

"REBCON"|"RebCon"|"rebcon" 		  { return REBUILDCONSTRAINT; }

"SNSLOT"|"SnSlot"|"snslot"                { return SAMPNODESLOT; }

"REDSTIFF"|"RedStiff"|"redstiff"	  { return REDUCEDSTIFFNESS; }

"UDBASIS"|"UDBasis"|"udbasis"		  { return UDEIMBASIS; }

"FPODRB"|"FPODRb"|"fpodrb"                { return FORCEROB; }

"CPODRB"|"CPODrb"|"cpodrb" 	          { return CONSTRAINTROB; }

"FLTROW"|"FltRow"|"fltrow"		  { return FILTERROWS; }

"svdfrc"|"SvdFrc"|"svdfrc"		  { return SVDFORCESNAP; }

"svdcon"|"SvdCon"|"svdcon" 		  { return SVDCONSTRAINTSNAP; }

"DEIM"|"Deim"|"deim"                      { return DEIMINDICES; }
 
"UDEIM"|"UDeim"|"udeim"                   { return UDEIMINDICES; }

"ROBC"|"RobC"|"robc" 			  { return SVDTOKEN; } 

"SOLVE"|"solve"                           { return ROBCSOLVE; } /* deprecated */
"COMPRESS"|"Compress"|"compress"          { return ROBCSOLVE; }

"MNORMA"|"MNorma"|"mnorma"		  { return NORMALIZETOKEN; }

"mi"|"MI"                                 { yylval.ival = 1;  return FLAGMN; }
"im"|"IM"                                 { yylval.ival = 0;  return FLAGMN; }
"id"|"ID"                                 { yylval.ival = -1; return FLAGMN; }

"MPRO"[a-zA-Z]*|"mpro"[a-zA-Z]*           { return MPROJECT; }

"blcksz"|"BlckSz"|"BLCKSZ"                { return BLOCKSIZE; }

"RODC"|"RodC"|"rodc"                      { return CONVERSIONTOKEN; } /* reduced order data conversion token */
 
"CONFIL"|"ConFil"|"confil" 	 	  { return CONVFI; } /* reduced order data conversion files */

"ROMRES"|"RomRes"|"romres"                { return ROMRES; }

"SNAPFI"|"SnapFi"|"snapfi"		  { return SNAPFI; }
"VELSNAPFI"|"VelSnapFi"|"velsnapfi"       { return VELSNAPFI; }
"ACCSNAPFI"|"AccSnapFi"|"accsnapfi"       { return ACCSNAPFI; }
"DSVSNAPFI"|"DsvSnapFi"|"dsvsnapfi"       { return DSVSNAPFI; }
"MUVSNAPFI"|"MuvSnapFi"|"muvsnapfi"       { return MUVSNAPFI; }

"ROBFI"|"RobFi"|"robfi"		          { return ROBFI; }

"SNAPWEIGHT"|"SnapWeight"|"snapweight"    { return SNAPWEIGHT; } /* deprecated */
"SSCALI"|"SScali"|"sscali"                { return SNAPWEIGHT; }

"CONCFG"|"ConCfg"|"concfg"                { return CONWEPCFG; }

"POSCFG"|"PosCfg"|"poscfg"                { return SCALEPOSCOORDS; }

"NDSCFG"|"NdsCfg"|"ndscfg"		  { return NODEPOSCOORDS; }

"MESHSCALE"|"MeshScale"|"meshscale"       { return MESHSCALEFACTOR; }

"NORMFI"|"NormFi"|"normfi"		  { return VECTORNORM; }

"PODROB"|"PodRob"|"podrob"		  { return PODROB; }

"TRNVCT"|"TrnVct"|"trnvct"		  { return TRNVCT; } /* deprecated */
"TRNSOL"|"TrnSol"|"trnsol"                { return TRNVCT; }

"OFFSET"|"OffSet"|"offset"		  { return OFFSET; }

"DIMENS"|"Dimens"|"dimens"		  { return PODSIZEMAX; }

"ROMENR"|"RomEnr"|"romenr"                { return ROMENERGY; }

"ORTHOG"|"Orthog"|"orthog"		  { return ORTHOG; }

"RMSHC"|"Rmshc"|"rmshc" 		  { return SAMPLING; }

"SNAPPR"|"SnapPr"|"snappr"                { return SNAPSHOTPROJECT; }

"TOLOPT"|"TolOpt"|"tolopt" 		  { return TOLER; } /* deprecated */
"TOLERA"|"Tolera"|"tolera"                { return TOLER; }

"NUMSKP"|"NumSkp"|"numskp"		  { return SKIP; }
"SAMFRQ"|"SamFrq"|"samfrq"                { return SKIP; }
"RANSMP"|"RanSmp"|"ransmp"    	          { return RANDOMSAMPLE; }

"SOLVER"|"Solver"|"solver"                { return SOLVER; }

"NPMAX"|"npMax"|"npmax"                   { return NPMAX; }

"BSSPLH"|"bsSplh"|"bssplh"                { return BSSPLH; }

"PGSPLH"|"pgSplh"|"pgsplh"                { return PGSPLH; }

"lh"|"LH"                                 { yylval.ival = 0; return SPNNLSSOLVERTYPE; }
"nncgp"|"NNCGP"                           { yylval.ival = 1; return SPNNLSSOLVERTYPE; } /* deprecated */
"CGNNLS"|"cgnnls"                         { yylval.ival = 1; return SPNNLSSOLVERTYPE; }  
"gpfp"|"GPFP" 				  { yylval.ival = 2; return SPNNLSSOLVERTYPE; } /* deprecated */
"CGPFP"|"cgpfp"                           { yylval.ival = 2; return SPNNLSSOLVERTYPE; } 
"lars"|"LARS"				  { yylval.ival = 3; return SPNNLSSOLVERTYPE; } /* deprecated */
"CHLASSO"|"chlasso"                       { yylval.ival = 3; return SPNNLSSOLVERTYPE; }
"mp"|"MP"				  { yylval.ival = 4; return SPNNLSSOLVERTYPE; }
"omp"|"OMP"                               { yylval.ival = 5; return SPNNLSSOLVERTYPE; }
"splh"|"SPLH" 				  { yylval.ival = 6; return SPNNLSSOLVERTYPE; } /* deprecated */
"QRNNLS"|"qrnnls"                         { yylval.ival = 6; return SPNNLSSOLVERTYPE; } 
"cglars"|"CGLARS"                         { yylval.ival = 7; return SPNNLSSOLVERTYPE; } /* deprecated */
"CGLASSO"|"cglasso"                       { yylval.ival = 7; return SPNNLSSOLVERTYPE; }
"sppfp"|"SPPFP"				  { yylval.ival = 8; return SPNNLSSOLVERTYPE; } /* deprecated */
"QRPFP"|"qrpfp"                           { yylval.ival = 8; return SPNNLSSOLVERTYPE; } 

"MAXSZE"|"MaxSze"|"maxsze"                { return MAXSIZE; }

"REFSUBTRACT"|"refsubtract"|"RefSubtract" { return REFSUBTRACT; } 

"USE_NMF"|"use_nmf"                       { return USENMF; }
"USE_NMFC"|"use_nmfc"                     { return USENMFC; }
"USE_GREEDY"|"use_greedy"                 { return USEGREEDY; }
"USE_PQN"|"use_pqn"                       { return USEPQN; }

"DO_CLUSTERING"|"do_clustering"           { return DOCLUSTERING; }
"DO_ROWCLUSTERING"|"do_rowclustering"     { return DOROWCLUSTERING; }
"ANGLE"|"angle" 			  { return ANGLE; }
"CLUSTERSOLVER"|"clustersolver"           { return CLUSTERSOLVER; }
"CRANDOM"|"crandom"                       { yylval.ival = 0; return CLUSTERSOLVERTYPE; }
"KMEANS"|"kmeans"                         { yylval.ival = 1; return CLUSTERSOLVERTYPE; }
"SSC"|"ssc" 	           	          { yylval.ival = 2; return CLUSTERSOLVERTYPE; }

[A-Za-z][A-Za-z0-9_.]*  { yylval.strval = strdup(yytext); return FNAME; }

"\""[A-Za-z0-9_./]*"\"" { yytext[yyleng-1] = 0; yylval.strval = strdup(yytext+1);
                           return FNAME; }

"\n" { line_num++; /* Count line numbers */ return NewLine; }

{SPACES} /* Just discard spaces */

[+-]?{DIGIT}+ {yylval.ival = atoi(yytext) ; return IntConstant ; }

[+-]?{DIGIT}+"."{DIGIT}* |
[+-]?"."{DIGIT}+ |
[+-]?{DIGIT}+"."{DIGIT}*[eE][+-]?{DIGIT}+ |
[+-]?{DIGIT}+[eE][+-]?{DIGIT}+ |
[+-]?"."{DIGIT}+[eE][+-]?{DIGIT}+ {
       yylval.fval = atof(yytext) ; return DblConstant ; }

[+-]?{DIGIT}+"."{DIGIT}*[dD][+-]?{DIGIT}+ |
[+-]?{DIGIT}+[dD][+-]?{DIGIT}+ |
[+-]?"."{DIGIT}+[dD][+-]?{DIGIT}+ {
       int il;
       for(il = 0; yytext[il] != 'd' && yytext[il] != 'D';) ++il;
       yytext[il] = 'e'; yylval.fval = atof(yytext) ; return DblConstant ; }

"*" { while(yyinput() != '\n') {}; line_num++; }



%%
void
yyerror(const char  *)
{
 	fprintf(stderr," ** Syntax Error near line %d : %s\n", line_num,yytext);
}

int yywrap()
{
	//std::cout << "yywrap" << endl;
 	return 1;
}
