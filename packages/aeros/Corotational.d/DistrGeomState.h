#ifndef _DISTR_GEOM_STATE_H_
#define _DISTR_GEOM_STATE_H_

template <class Scalar> class GenDecDomain;
typedef GenDecDomain<double> DecDomain;
class GeomState;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;

#include <map>
#include <vector>
#include <cstddef>

class DistrGeomState {
   private:
     GeomState **gs;	// pointer to an array of GeomStates
     int numSub;

  public:
     std::map<std::pair<int,int>, double> *mu; // lagrange multipliers for the contact surfaces
     std::vector<double> *lambda;              // lagrange multipliers for all other constraints

     // Constructor
     DistrGeomState(DecDomain *domain);
     // Copy Constructor
     DistrGeomState(const DistrGeomState &);
     // Destructor
     ~DistrGeomState();

     void resize(DecDomain *domain, std::map<std::pair<int,int>,double> *mu = NULL);

     // return the ith subdomain's GeomState
     GeomState* operator[](int i) const { return gs[i]; }
     GeomState* getSubGeomState(int i) { return gs[i]; }

     // Update the GeomStates
     void update(DistrVector &v, int SO3param = 0);
     void update(DistrGeomState &ref, DistrVector &v, int S03param = 0);
     void explicitUpdate(GenDecDomain<double> *decDomain, DistrVector &v);
     void setVelocity(DistrVector &, int SO3param = 0);
     void setAcceleration(DistrVector &, int SO3param = 0);
     void setVelocityAndAcceleration(DistrVector &, DistrVector &);
     void setNodalTemperatures(DistrVector &);

// The following functions are necessary to implement NL dynamics and
// the arclength method
     int getTotalNumElemStates();
     void midpoint_step_update(DistrVector &veloc_n, DistrVector &accel_n, double &delta, DistrGeomState &ss,
                               double beta, double gamma, double alphaf, double alpham, bool zeroRot);
     void get_inc_displacement(DistrVector &inc_Vec, DistrGeomState &ss, bool zeroRot);
     void push_forward(DistrVector &f);
     void pull_back(DistrVector &f);
     void push_forward_trans(DistrVector &f);
     void pull_back_trans(DistrVector &f);
     void transform(DistrVector &f, int, bool unscaled = false);
     void get_tot_displacement(DistrVector &totVec, bool rescaled = true);
     void interp(double, DistrGeomState &, DistrGeomState &);
     void diff(DistrGeomState &unp, DistrVector &un);
     void print();

     DistrGeomState &operator=(DistrGeomState &unp);     
     void subCopy(int isub, DistrGeomState &unp);

     int getNumSub() const { return numSub; }
     bool getHaveRot();

     void getMultipliers(std::map<std::pair<int,int>,double> &mu);
     void setMultipliers(std::map<std::pair<int,int>,double> &mu);

  private:
     void subStep_update(int isub, DistrVector &veloc_n, DistrVector &accel_n,
                         double &delta, DistrGeomState &ss,
                         double beta, double gamma, double alphaf, double alpham, bool zeroRot);
     void subInc_get(int isub, DistrVector &inc_Vec, DistrGeomState &ss, bool zeroRot);
     void subTot_get(int isub, DistrVector &totVec, bool rescaled);
     void subInterp(int isub, double&, DistrGeomState &, DistrGeomState &);
     void subDiff(int isub, DistrGeomState &unp, DistrVector &un);
     void subUpdate(int isub, DistrVector &v, int SO3param);
     void subUpdateRef(int isub, DistrGeomState &ref, DistrVector &v, int SO3param);
     void subExplicitUpdate(int iSub,DistrVector &v, GenDecDomain<double> *decDomain);
     void subSetVelocity(int isub, DistrVector &v, int SO3param);
     void subSetAcceleration(int isub, DistrVector &a, int SO3param);
     void subSetVelocityAndAcceleration(int isub, DistrVector &v, DistrVector &a);
     void subSetNodalTemperatures(int isub, DistrVector &temps);
     void makeSubGeomStates(int isub, DecDomain *domain);
     void subCopyConstructor(int isub, const DistrGeomState &g2);
     void subPushForward(int isub, DistrVector &f);
     void subPullBack(int isub, DistrVector &f);
     void subTransform(int isub, DistrVector &f, int, bool);
     void subResize(int isub, DecDomain* domain);
     void subSetMultipliers(int isub, std::map<std::pair<int,int>,double> &mu);
};

#endif
