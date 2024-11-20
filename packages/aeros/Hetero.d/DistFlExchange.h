#ifndef _DISTFLEXCHANGE_H_
#define _DISTFLEXCHANGE_H_

#include <Element.d/Element.h>
#include <Utils.d/OutputInfo.h>
#include <Hetero.d/FlExchange.h>

#include <map>
#include <algorithm>

class CoordSet;
class State;
class Connectivity;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class VecType> class SysState;
class DistrGeomState;
class SurfaceEntity;
template <class Scalar> class GenSubDomain;
typedef GenSubDomain<double> SubDomain;

#define FL_NEGOT 10000

typedef std::map<int, InterpPoint> MatchMap;

class DofSetArray;

class DistFlExchanger {
  double *buffer, *buff;
  int bufferLen, buffLen;

  int nSender;             // dimension of sndTable
  int nbrReceivingFromMe;  // number of fluid mpi's w/matches in this mpi
  int *idSendTo;  	   // list of fluid mpi's to sendTo
  int *nbSendTo;	   // num of match data per fluid neighbor
  int *consOrigin;         // reverse table of idSendTo
  InterpPoint **sndTable;  // match data by local subdomain in mpi

  CoordSet **cs;           // nodes in this mpi process
  Elemset **eset;
  DofSetArray **cdsa;
  DofSetArray **dsa;

  double *localF;
  double aforce[4];
  double aflux;
  int rcvParity, sndParity;
  OutputInfo *oinfo;
  int algnum;
  int isCollocated;
  double alpha[2];
  double alph[2];
  double dt;
  double dtemp;
  DistrVector *tmpDisp;
  Vector *dsp;
  Vector *vel;
  Vector *acc; 
  Vector *pVel; 

  //PJSA (Nov.18,2010): FS Communication using Face Elements
  SurfaceEntity *surface;
  CoordSet *globalCoords;
  Connectivity *nodeToElem, *elemToSub;
  bool useFaceElem;
  SubDomain **sd;
  int **fnId2;

  bool wCracking;
  bool sentInitialCracking;
  Connectivity *faceElemToNode, *nodeToFaceElem;

public:

  DistFlExchanger(CoordSet **, Elemset **, DofSetArray **, 
                  DofSetArray **, OutputInfo *oinfo = 0);
  DistFlExchanger(CoordSet **, Elemset **, SurfaceEntity *, CoordSet *,
                  Connectivity *, Connectivity *, SubDomain **,
                  DofSetArray **, DofSetArray **, OutputInfo *oinfo, bool wCracking);
  ~DistFlExchanger();

  MatchMap* getMatchData();
  void negotiate();
  void thermoread(int);

  void sendDisplacements(SysState<DistrVector>&, double**, double**, int = -1, DistrGeomState* = 0);
  void sendTemperature(SysState<DistrVector>&);
  void sendStrucTemp(DistrVector&);
  
  double getFluidLoad(DistrVector&, int, double, double, int&, DistrGeomState* = 0);
  double getFluidFlux(DistrVector& flux, double time, double &bflux);
  void getStrucTemp(double*);

  void sendParam(int, double, double, int, int, double a[2]);
  void sendSubcyclingInfo(int sub);

  void sendTempParam(int algnum, double step, double totaltime,
                     int rstinc, double alphat[2]);

  void setBufferLength(int size)  { bufferLen = size; }
 
  void sendModeFreq(double *modFrq, int numFrq);
  void sendModeShapes(int numFrq, int nNodes, double (**)[6],
                      SysState<DistrVector>&, double factor = 1.0);

  void initSndParity(int pinit) { sndParity = pinit; }
  void initRcvParity(int pinit) { rcvParity = pinit; }
  void flipRcvParity() { if(rcvParity >= 0) rcvParity = 1-rcvParity; }
  void flipSndParity() { if(sndParity >= 0) sndParity = 1-sndParity; }

  void sendNoStructure();
  void sendNewStructure();

  int cmdCom(int);

  void sendEmbeddedWetSurface();

private:
  void transformVector(double *localF, Element *ele, int locSub);
  void transformVector(double *localF, FaceElement *ele, int locSub);
};

#define FLTOSTMT 1000
#define STTOFLMT 2000
#define STTOSTMT 4000
#define FLTOSTHEAT 5000
#define STTOFLHEAT 6000
#define STCMDMSG 8000
#define FLCMDMSG 9000
#define OPTPARMSG 8100
#define OPTRESMSG 9100
#define NBPRESSDATAMAX 7
#define FL_NEGOT 10000

#endif

