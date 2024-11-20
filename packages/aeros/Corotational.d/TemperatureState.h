#ifndef _TEMP_STATE_H_
#define _TEMP_STATE_H_

#include <Corotational.d/GeomState.h>

class DofSetArray;
class CoordSet;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
class ControlLawInfo;
class BCond;

class TemperatureState : public GeomState
{
   public:
     // Default Constructor
     TemperatureState();

     // Constructor
     TemperatureState(DofSetArray &dsa, DofSetArray &cdsa, CoordSet &cs);

     // Copy Constructor
     TemperatureState(const TemperatureState &);

     ~TemperatureState() { }

     void update(const Vector &, int = 0);
     void updatePrescribedDisplacement(BCond *dbc, int numDirichlet, 
                                       double delta);
     void updatePrescribedDisplacement(BCond *dbc, int numDirichlet,
                                       CoordSet &cs);
     void explicitUpdate(CoordSet &cs, const Vector &v);

     void midpoint_step_update(Vector &veloc_n, Vector &accel_n, double delta, GeomState &ss,
                               double beta, double gamma, double alphaf, double alpham,
                               bool);
     void get_inc_displacement(Vector &inc_Vec, GeomState &ss, bool);
     void get_tot_displacement(Vector &totVec, bool);

     void setVelocity(const Vector &, int = 0);
     void setAcceleration(const Vector &, int = 0);
     void setVelocityAndAcceleration(const Vector &, const Vector &);

     void push_forward(Vector &) {}
     void pull_back(Vector &) {}
     void transform(Vector &, int, bool = false) const {}
};

#endif
