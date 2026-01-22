#ifndef _GEOM_STATE_H_
#define _GEOM_STATE_H_

#include<IoData.h>
#include<Vector3D.h>

//struct SpringElement {
//  double l; // natural length
//  double k; // stiffness
//  double m; // mass
//  double u0; // initial displacement
//};

struct GeomState {

  //! spring parameters
  double stiffness;
  double mass;

  //! spring state
  double disp;
  double velo;
  double accl;

  GeomState() = delete;
  GeomState(SpringData &spring); 
  ~GeomState() = default;

};

inline GeomState::GeomState(SpringData &spring)
{

  mass = spring.mass;
  stiffness = spring.stiffness;

  disp = spring.init_conds.disp;
  velo = spring.init_conds.velo;

}


#endif
