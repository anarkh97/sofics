#ifndef _FLEXCHANGE_H_
#define _FLEXCHANGE_H_

#include <Element.d/Element.h>
#include <Utils.d/OutputInfo.h>

class CoordSet;
class State;
class GeomState;
class SurfaceEntity;
class FaceElement;
class Connectivity;

// This is an inerpolation point. The element elemNum uses x and y to
// compute the interpolated displacements/velocities

struct InterpPoint {
    int subNumber;
    int elemNum;
    double xy[2];
    double gap[3];
    int *dofs;
};

class DofSetArray;

class FlExchanger {
 
     //KW (Jul.27,2010): FS Communication using Face Elements
     SurfaceEntity *surface;
     bool useFaceElem;

     double *buffer;
     double *buff;
     int bufferLen;
     int buffLen;

     int nbrSendingToMe;
     
     int nbrReceivingFromMe;
     int *idSendTo;
     int *nbSendTo;
     int *consOrigin; // reverse table of idSendTo
     InterpPoint **sndTable;

     CoordSet& cs;
     Elemset& eset;
     DofSetArray* dsa;

     double *localF;
     double aforce[4];
     double aflux;
     int rcvParity, sndParity;
     OutputInfo *oinfo;
     int algnum;
     int isCollocated;
     double alpha[2];
     double alphasv;
     double alph[2];
     double dt;
     double dtemp;
     Vector *tmpDisp;
     Vector *tmpVel;
     
     bool wCracking;
     bool sentInitialCracking;
     Connectivity *faceElemToNode, *nodeToFaceElem;

   public:
     //KW (Jul.27,2010): FS Communication using Face Elements
     FlExchanger(CoordSet&, Elemset&, SurfaceEntity*, DofSetArray *, OutputInfo *oinfo, bool wCracking);
     FlExchanger(CoordSet&, Elemset&, DofSetArray *, OutputInfo *oinfo = 0);
     ~FlExchanger();
     void read(int mynode, const char *filename);
     void negotiate();
     void thermoread(int &buffLen);

     void sendDisplacements(State& state, int tag = -1, GeomState* geomState = 0);
     void sendTemperature(State &state);
     void waitOnSend();
     double getFluidLoad(Vector& force, int tIndex, double time,
                         double alphaf, int& iscollocated, GeomState* geomState = 0);
     double getFluidLoadSensitivity(Vector& forceSen, int tIndex, double time,
                                    double alphaf, int& iscollocated, GeomState* geomState = 0);

     double getFluidFlux(Vector &flux, double time, double &fl);
     void sendStrucTemp(Vector &tempsent);
     void getStrucTemp(double *temprcvd);
     void sendHeatSource(Vector &heatsent);
     void getHeatSource(double *heatrcvd);
     
     void sendParam(int alg, double step, double totalTime,
                    int restartinc, int _isCollocated, double alphas[2], double alphasv);
     void sendSubcyclingInfo(int sub);

     void sendTempParam(int algnum, double step, double totaltime,
                        int rstinc, double alphat[2]);

     void sendModeFreq(double *modFrq, int numFrq);
     void sendModeShapes(int numFrq, int nNodes, double (**)[6],
                         State &st, double factor = 1.0, int numIDis6 = 0, BCond* iDis6=0);

     void sendEmbeddedWetSurface();
     void printreceiving();

      void initSndParity(int pinit) { sndParity = pinit; }
      void initRcvParity(int pinit) { rcvParity = pinit; }
      void flipSndParity() { if(sndParity >= 0) sndParity = 1-sndParity; }
      void flipRcvParity() { if(rcvParity >= 0) rcvParity = 1-rcvParity; }

      void sendNoStructure();
      void sendNewStructure(std::set<int> &newDeletedElements);

      void sendNumParam(int,int,double);
      void getNumParam(int&);
      void sendRelativeResidual(double);
      int cmdCom(int);
      int cmdComHeat(int);

   private:
      void transformVector(double *localF, Element *ele);
      void transformVector(double *localF, FaceElement *ele);
};

#define FLTOSTMT 1000
#define STTOFLMT 2000
#define STTOSTMT 4000
#define FLTOSTHEAT 5000
#define STTOFLHEAT 6000
#define STTOSTHEAT 7000
#define STNUMPAFL 7500
#define FLNUMPAST 7550
#define STRELRESFL 7600
#define STCMDMSG 8000
#define FLCMDMSG 9000
#define OPTPARMSG 8100
#define OPTRESMSG 9100
#define NBPRESSDATAMAX 7
#define FL_NEGOT 10000

#endif
