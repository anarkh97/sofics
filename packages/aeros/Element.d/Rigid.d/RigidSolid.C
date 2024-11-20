#ifdef USE_EIGEN3
#include <Element.d/Rigid.d/RigidSolid.h>
#include <Element.d/Rigid.d/RigidTwoNodeTruss.h>
#include <iostream>

RigidSolid::RigidSolid(int _nnodes, int* _nn)
 : SuperElement(true)
{
  if(_nnodes < 3) {
    std::cerr << " *** ERROR: minimum number of nodes is 3 for rigid solid. Exiting...\n";
    exit(-1);
  }
  nnodes = _nnodes;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
}

void
RigidSolid::buildFrame(CoordSet& cs)
{
  nSubElems = 3*(nnodes-2);
  subElems = new Element * [nSubElems];

  // Finding 3 nodes that make the largest triangle
  double LargestArea = 0.0;
  double x[3], y[3], z[3];
  int best[4];

  for (int i = 0; i < nnodes; ++i) {
    Node &Nd1 = cs.getNode(nn[i]);
    x[0] = Nd1.x; y[0] = Nd1.y; z[0] = Nd1.z;

    for (int j = 0; j < nnodes; ++j) {
      Node &Nd2 = cs.getNode(nn[j]);
      x[1] = Nd2.x; y[1] = Nd2.y; z[1] = Nd2.z;

      double dx_12 = x[1] - x[0];
      double dy_12 = y[1] - y[0];
      double dz_12 = z[1] - z[0];

      for (int k = 0; k < nnodes; ++k) {
        Node &Nd3 = cs.getNode(nn[k]);
        x[2] = Nd3.x; y[2] = Nd3.y; z[2] = Nd3.z;

        double dx_13 = x[2] - x[0];
        double dy_13 = y[2] - y[0];
        double dz_13 = z[2] - z[0];

        double diff1 = dy_12 * dz_13 - dy_13 * dz_12;
        double diff2 = dx_13 * dz_12 - dx_12 * dz_13;
        double diff3 = dx_12 * dy_13 - dx_13 * dy_12;

        double Area = sqrt( diff1 * diff1 + diff2 * diff2 + diff3 * diff3 );
        if (Area >= LargestArea) {
          LargestArea = Area;
          best[0] = i;
          best[1] = j;
          best[2] = k;
          best[3] = i;
        }
      }
    }
  }

  // Restraining the triangle
  int index = 0;
  for(int i = 0; i < 3; ++i)
    subElems[index++] = new RigidTwoNodeTruss(best+i);

  // Restraining the rest of the nodes
  for(int i = 0; i < nnodes; ++i) {
    if(i == best[0] || i == best[1] || i == best[2]) continue;
    for(int j = 0; j < 3; ++j) {
      int indices[2] = { i, best[j] };
      subElems[index++] = new RigidTwoNodeTruss(indices);
    }
  }

  SuperElement::buildFrame(cs);
}

int 
RigidSolid::getTopNumber() const
{ 
  switch(nnodes) {
    case 4 : return 123; // 4-node tetra
    case 6 : return 124; // 6-node penta
    case 8 : return 117; // 8-node hexa
    case 10 : return 125; // 10-node tetra
    case 15 : return 197; // 15-node penta
    case 20 : return 172; // 20-node hexa
    case 26 : return 192; // 26-node penta
    case 32 : return 191; // 32-node hexa
    default : return 101;
  }  
}

int 
RigidSolid::numTopNodes() const
{ 
  switch(nnodes) {
    case 4 : return 4; // 4-node tetra
    case 6 : return 6; // 6-node penta
    case 8 : return 8; // 8-node hexa
    case 10 : return 10; // 10-node tetra
    case 15 : return 15; // 15-node penta
    case 20 : return 20; // 20-node hexa
    case 26 : return 26; // 26-node penta
    case 32 : return 32; // 32-node hexa
    default : return 2;
  } 
}

bool
RigidSolid::isSafe() const
{
  return (nnodes >= 4);
}
#endif
