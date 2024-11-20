#ifndef _MODES_H_
#define _MODES_H_

// Structure to contain Mode Data
// i.e. number of modes, number of nodes, 
// array of frequencies, and matrix of mode shapes.

struct ModeData {
  int numModes;
  int numNodes;
  double *frequencies;
  double (**modes)[6];
  int *nodes;

  template<class VecType>
  void addMultY(int numY, BCond *Y, VecType& x, DofSetArray* dsa, int ndof=6) { //HB: x <- x + X.y 
    for(int j = 0; j < numNodes; ++j) // loop over nodes
      for(int k = 0; k < ndof; ++k) { // loop over dofs
         int dof = dsa->locate(nodes[j], 1 << k);
         if(dof >= 0)
           for(int i = 0; i<numY; ++i) // loop over selected modes
             x[dof] += Y[i].val*modes[Y[i].nnum][j][k];
      }
  }

  void addMultY(int numY, BCond *Y, BCond *iDis, int ndof=6) { // PJSA 7/14/2010
    for(int j = 0; j < numNodes; ++j) // loop over nodes
      for(int k = 0; k < ndof; ++k) { // loop over dofs
        iDis[j*ndof+k].nnum = j;
        iDis[j*ndof+k].dofnum = k;
        iDis[j*ndof+k].val = 0;
        for(int i = 0; i < numY; ++i) // loop over selected modes
          iDis[j*ndof+k].val += Y[i].val*modes[Y[i].nnum][j][k];
      }
  }


};

#endif
