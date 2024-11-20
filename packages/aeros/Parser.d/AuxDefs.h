#ifndef _AUXDEFS_H_
#define _AUXDEFS_H_

#include <Element.d/Element.h>
#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/FabricMap.h>
#include <Element.d/NonLinearity.d/FabricMat.h>
#include <Element.d/Sommerfeld.d/SommerElement.h>
#include <Element.d/Sommerfeld.d/LineSommerBC.h>
#include <Element.d/Sommerfeld.d/Line2SommerBC.h>
#include <Element.d/Sommerfeld.d/TriangleSommerBC.h>
#include <Element.d/Sommerfeld.d/QuadSommerBC.h>
#include <Element.d/Sommerfeld.d/Triangle6SommerBC.h>
#include <Element.d/Sommerfeld.d/IsoParamQuadSommer.h>
#include <Element.d/Sommerfeld.d/IsoParamTriSommer.h>

#include <Parser.d/DecInit.h>

class Domain;
class GeoSource;
class Sfem;
template <class Scalar, class VecType> class SfemNonInpc;
template <class Scalar, class VecType> class SfemInpc;

typedef struct {
   int num;
   double xyz[3];
   int cp;
   int cd;
} NumedNode;
   
typedef struct {
   int num;
   int nd[500];
} NumList;

struct NumedElem {
  int num;
  Element *elem;
};

class BCList {
  public:
    BCond *d;
    int n;
    int maxbc;
    int loadsetid;

    BCList(int _loadsetid = 0);
    void add(BCond &);
    void add(int nd, int dof, double val) 
      { BCond bc; bc.nnum = nd; bc.dofnum = dof; bc.val = val; add(bc); }
};

class ComplexBCList {
  public:
    ComplexBCond *d;
    int n;
    int maxbc;

    ComplexBCList();
    void add(ComplexBCond &);
};

struct FrameData {
  int num;
  double d[9];
};

struct NodalFrameData {
  int id;
  double o[3];
  double d[9];
  int type;
};

struct LayerData {
  int lnum;
  int matid;
  double d[12];
};

struct ConstraintOptions {
  bool lagrangeMult;
  double penalty;
  int constraint_hess;
  double constraint_hess_eps;
};

struct StringList {
     int nval;
     const char *v[32];
};

template<typename T1, typename T2>
struct TrivialPair {
  T1 first;
  T2 second;
};

extern double fetiHTol;
extern int fetiHIterMax;

extern DecInit * decInit;

extern Domain *domain;
extern GeoSource *geoSource;

extern int line_num;
extern int numColumns;

int yylex(void);

void yyerror(const char*);

// dynamic weight for every type of elements
// element numbers are those assigned in AddElem.C
extern std::map<int, double> weightList;  
extern std::map<int, double> fieldWeightList;

#endif
