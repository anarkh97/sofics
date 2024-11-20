#ifndef _SYSSTATE_H_
#define _SYSSTATE_H_

template <class VecType>
class SysState {
   VecType &d_n, &v_n, &a_n, &v_n_p;
 public:

   SysState(VecType &v): d_n(v), v_n(v), a_n(v), v_n_p(v) { }

   SysState(VecType &d, VecType &v, VecType &a, VecType &vp):
     d_n(d), v_n(v), a_n(a), v_n_p(vp) { }

   SysState(VecType &d, VecType &v, VecType &vp):
     d_n(d), v_n(v), a_n(v), v_n_p(vp) { }

   SysState &operator=(const SysState &sys);

   VecType &getDisp()      { return d_n; }
   VecType &getVeloc()     { return v_n; }
   VecType &getPrevVeloc() { return v_n_p; }
   VecType &getAccel()     { return a_n; }

   VecType &getDispConst()      const { return d_n; }
   VecType &getVelocConst()     const { return v_n; }
   VecType &getPrevVelocConst() const { return v_n_p; }
   VecType &getAccelConst()     const { return a_n; }

};

#endif
