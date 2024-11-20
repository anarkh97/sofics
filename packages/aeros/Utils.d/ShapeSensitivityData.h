#ifndef _SHAPESENSITIVITYDATA_H_
#define _SHAPESENSITIVITYDATA_H_

// Structure to contain Shape Sensitivity Data
// i.e. number of shape variables, number of nodes, 
// sensitivity of nodal coordinate wrt shape variables.

struct ShapeSensitivityData {
  int numVars;
  int numNodes;
  int *index;
  double (**sensitivities)[3];
  int *nodes;

  template<class VecType>
  void addMultY(int numY, BCond *Y, VecType& x, DofSetArray* dsa, int ndof=3) { //HB: x <- x + X.y 
    // The following implementation assumes that all the sensitivities have been given for the SAME nodes
    // to reduce the number of calls to DofSetArray::locate
    for(int j = 0; j < numNodes; ++j) // loop over nodes
      for(int k = 0; k < ndof; ++k) { // loop over dofs
         int dof = dsa->locate(nodes[j], 1 << k);
         if(dof >= 0)
           for(int i = 0; i<numY; ++i) // loop over selected sensitivities
             x[dof] += Y[i].val*sensitivities[Y[i].nnum][j][k];
      }
  }

  void addMultY(int numY, BCond *Y, BCond *iDis, int ndof=3) { // PJSA 7/14/2010
    for(int j = 0; j < numNodes; ++j) // loop over nodes
      for(int k = 0; k < ndof; ++k) { // loop over dofs
        iDis[j*ndof+k].nnum = j;
        iDis[j*ndof+k].dofnum = k;
        iDis[j*ndof+k].val = 0;
        for(int i = 0; i < numY; ++i) // loop over selected sensitivities
          iDis[j*ndof+k].val += Y[i].val*sensitivities[Y[i].nnum][j][k];
      }
  }


};

#endif
