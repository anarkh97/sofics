/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <climits>
#include <cassert>
#include "UserDefinedState.h"
#include "Vector3D.h"
#include "Vector5D.h"

using namespace std;

//------------------------------------------------------------
// This is a template for users to fill. Not compiled w/ M2C.
// If multiple UserDefinedState need to be defined, they should
// have different names (e.g., StateCalculator1, StateCalculator2, etc.)
// Compilation script (an example):
// g++ -O3 -fPIC -I/path/to/folder/that/contains/UserDefinedState.h -c UserDefinedState.cpp; g++ -shared UserDefinedState.o -o UserDefinedState.so; rm UserDefinedState.o
//------------------------------------------------------------

class StateCalculator : public UserDefinedState{

public:
  void GetUserDefinedState(int i0, int j0, int k0, int imax, int jmax, int kmax,
                           Vec3D*** coords, Vec5D*** v, double*** id, vector<double***> phi);
};

//------------------------------------------------------------

void
StateCalculator::GetUserDefinedState(int i0, int j0, int k0, int imax, int jmax, int kmax,
                                     Vec3D*** coords, Vec5D*** v, double*** id, vector<double***> phi)
{
  // Note:
  // 1. Each processor core calls this functio with different inputs.
  // 2. Override the original values of v and id. 
  // 4. The state variables are each node are: rho,u,v,w,p.
  // 3. Do NOT change coords.

  for(int k=k0; k<kmax; k++)
    for(int j=j0; j<jmax; j++)
      for(int i=i0; i<imax; i++) {

        // left state
        if(coords[k][j][i][0] <= 990.0) {
          v[k][j][i][0] = 2e-6;
          v[k][j][i][1] = 0.0;
          v[k][j][i][2] = 0.0;
          v[k][j][i][3] = 0.0;
          v[k][j][i][4] = 2e5;
        }
        // right state
        else if(coords[k][j][i][0] > 990.0 and coords[k][j][i][0] <= 1000.0) {
          v[k][j][i][0] = 1e-6;
          v[k][j][i][1] = 0.0;
          v[k][j][i][2] = 0.0;
          v[k][j][i][3] = 0.0;
          v[k][j][i][4] = 1e5;
        }
        // state behind the structure
        else {
          v[k][j][i][0] = 1e-6;
          v[k][j][i][1] = 0.0;
          v[k][j][i][2] = 0.0;
          v[k][j][i][3] = 0.0;
          v[k][j][i][4] = 1e5;
        }

      }

}

//------------------------------------------------------------
// The class factory (Note: Do NOT change these functions except the word "StateCalculator".)
extern "C" UserDefinedState* Create() {
  return new StateCalculator;
}

extern "C" void Destroy(UserDefinedState* uds) {
  delete uds;
}

//------------------------------------------------------------

