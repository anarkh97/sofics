#ifndef _STATE_H_
#define _STATE_H_

#include <Utils.d/dofset.h>

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
extern Vector nillVec;
class CoordSet;

class State {
   DofSetArray *dsa;	// Constrained DOFSetArray
   DofSetArray *DSA;    // Entire dof Set Array
   double      *bcx;    // displacement prescribed values
   double      *vcx;    // velocity prescribed values
   Vector      &disp;	// displacement
   Vector      &veloc;	// velocity
   Vector      &accel;	// acceleration
   Vector      &prevVeloc;	// previous velocity
   CoordSet    *cs;
 public:
   State() : disp(nillVec), veloc(nillVec), accel(nillVec), prevVeloc(nillVec), cs(0)
         {} // Careful that you know what you are doing here

   State(const State &s) : dsa(s.dsa), DSA(s.DSA), bcx(s.bcx), vcx(s.vcx),
       disp(s.disp), veloc(s.veloc), accel(s.accel), 
       prevVeloc(s.prevVeloc), cs(s.cs)
     {}

   // Constructor
   State (DofSetArray *_dsa, Vector &d, Vector &v, Vector &a, Vector &pv) :
      disp(d), veloc(v), accel(a), prevVeloc(pv), cs(0)
    {dsa = _dsa; }

   // Constructor
   State (DofSetArray *_dsa, DofSetArray *_DSA, 
          double *_bcx, double *_vcx, Vector &d, Vector &v, Vector &a,
          Vector &pv, CoordSet *_cs = 0) :
      disp(d), veloc(v), accel(a), prevVeloc(pv), cs(_cs)
      { dsa = _dsa; DSA = _DSA; bcx = _bcx; vcx = _vcx; }

   // For heat problem
   State (DofSetArray *_dsa, DofSetArray *_DSA,
          double *_bcx, Vector &d, Vector &v, Vector &pv) :
      disp(d), veloc(v), accel(v), prevVeloc(pv), cs(0)   /*accel has to be here, but useless*/
      { dsa = _dsa; DSA = _DSA; bcx = _bcx; vcx = 0; }

   // For heat problem
   State (DofSetArray *_dsa, Vector &d, Vector &v, Vector &pv) :
      disp(d), veloc(v), accel(v), prevVeloc(pv), cs(0) { dsa = _dsa; }

   State (const State &s, Vector &n_d) : disp(n_d), veloc(s.veloc),
        accel(s.accel), prevVeloc(s.prevVeloc), cs(s.cs)
   {
     dsa = s.dsa;
     DSA = s.DSA;
     bcx = s.bcx;
     vcx = s.vcx;
   }

   State (const State &s, Vector &n_d, Vector &n_v) : disp(n_d), veloc(n_v),
        accel(s.accel), prevVeloc(s.prevVeloc), cs(s.cs)
   {
     dsa = s.dsa;
     DSA = s.DSA;
     bcx = s.bcx;
     vcx = s.vcx;
   }


   Vector &getDisp()  { return disp;  }
   Vector &getVeloc() { return veloc; }
   Vector &getAccel() { return accel; }
   Vector &getPrevVeloc() { return prevVeloc; }

   // get translational displacements and velocities
   void getDV(int node, double[3], double[3], bool noPrescribed = false);

   // get translational and rotational displacements and velocities
   void getDVRot(int node, double[6], double[6]);

   // get temperature and its first time derivative
   void getTemp(int node, double[1], double[1]);
};

#endif
