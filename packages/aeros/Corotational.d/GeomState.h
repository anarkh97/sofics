#ifndef _GEOM_STATE_H_
#define _GEOM_STATE_H_

#include <map>
#include <vector>

class DofSetArray;
class CoordSet;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
class ControlLawInfo;
class BCond;
class Node;
class Elemset;

class NodeState {
  public:
    double x, y, z;                     // x,y,z coordinates
    double v[6], a[6];                  // x,y,z velocities and accelerations
                                        // note: by convention angular velocities and accelerations are
                                        //       convected, except in the case of implicit ROM in which 
                                        //       case they are time derivatives of the rotation vector
    double theta[3];                    // rotation vector
    double R[3][3];                     // rotation tensor
    double temp;                        // temperature
    void operator=(const NodeState &);
    double diff(const Node &un, int dof);
    NodeState() { for(int i = 0; i < 6; ++i) v[i] = a[i] = 0;
                  for(int i = 0; i < 3; ++i) theta[i] = 0; }
};

class ElemState {
  public:
    int numInternalStates;
    double *internalStates;
    ElemState() : numInternalStates(0), internalStates(0) {}
    ~ElemState() { if(internalStates) delete [] internalStates; }
    void operator=(const ElemState &);
};

class GeomState {
  protected:
     std::vector<NodeState> ns; // node state (x,y,z position and rotation tensor)
     int numnodes, numnodesFixed; // number of nodes
     std::vector<std::vector<int> > loc; // dof location array
     double refCG[3];   // reference CG
     double gRot[3][3]; // Global Rotation Matrix
     const CoordSet *X0;
     std::vector<int> flag; // signifies if node is connected to element
     int numelems;
     std::vector<ElemState> es;
     std::map<int,int> emap;
     std::map<std::pair<int,int>,int> multiplier_nodes;
     bool haveRot;

   public:
     // Default Constructor
     GeomState();
     GeomState(const CoordSet &cs);

     // Constructor
     GeomState(DofSetArray &dsa, DofSetArray &cdsa, CoordSet &cs, Elemset *elems = 0, double *ndTemps = 0);

     // Copy Constructor
     GeomState(const GeomState &);

     virtual ~GeomState();

     void resize(DofSetArray &dsa, DofSetArray &cdsa, CoordSet &cs, Elemset *elems = 0);

     NodeState & operator[](int i)  { return ns[i]; }
     const NodeState & operator[](int i) const { return ns[i]; }
     void clearMultiplierNodes();
     void resizeLocAndFlag(DofSetArray &cdsa);
     int getNodeFlag(int i) const { return flag[i]; }

     double * getElemState(int glNum) { return (numelems > 0) ? es[emap[glNum]].internalStates : 0; }
     int getNumElemStates(int glNum) { return (numelems > 0) ? es[emap[glNum]].numInternalStates : 0; }
     int getTotalNumElemStates() const;

     // int getLocation(int inode, int dof) { return (loc[inode][dof]-1); }
     int numNodes() const { return numnodes; }
     int numNodesFixed() const { return numnodesFixed; }

     void getPositions(double *positions);
     void getRotations(double *rotations);
     void getVelocities(double *velocities);
     void getAccelerations(double *accelerations);
     void getElemStates(double *elemStates) const;

     void setPositions(double *positions);
     void setRotations(double *rotations);
     void setVelocities(double *velocities);
     void setAccelerations(double *accelerations);
     void setNodalTemperatures(double *ndTemps);
     void setElemStates(double *elemStates);

     void extract(double *p);

     virtual void update(const Vector &, int SO3param = 0);
     virtual void update(const Vector &, const std::vector<int> &, int SO3param = 0);
     virtual void explicitUpdate(CoordSet &cs, const Vector &v);
     virtual void explicitUpdate(CoordSet &cs, int numNodes, int* nodes, const Vector &v);
     virtual void update(GeomState &refState, const Vector &, int SO3param = 0);
     virtual void setVelocity(const Vector &, int SO3param = 0);
     virtual void setVelocity(int numNodes, int* nodes, const Vector &, int SO3param = 0);
     virtual void setAcceleration(const Vector &, int SO3param = 0);
     virtual void setAcceleration(int numNodes, int* nodes, const Vector &, int SO3param = 0);
     virtual void setVelocityAndAcceleration(const Vector &, const Vector &);
     virtual void updatePrescribedDisplacement(BCond *dbc, int numDirichlet, 
                                               double delta);
     virtual void updatePrescribedDisplacement(BCond *dbc, int numDirichlet,
                                               CoordSet &cs);
     void updatePrescribedDisplacement(double *userDefinedDisplacement,
                                       ControlLawInfo* claw,
                                       CoordSet &cs, double *userDefinedVel, double *userDefinedAcc);
     virtual void midpoint_step_update(Vector &veloc_n, Vector &accel_n, double delta, GeomState &ss,
                                       double beta, double gamma, double alphaf, double alpham,
                                       bool zeroRot);
     virtual void get_inc_displacement(Vector &inc_Vec, GeomState &ss, bool zeroRot);
     virtual void get_tot_displacement(Vector &totVec, bool rescaled = true);
     virtual void get_temperature(int numNodes, int* nodes, Vector &ndTemps, double Ta);
     virtual void push_forward(Vector &f);
     virtual void pull_back(Vector &f);

     virtual void transform(Vector &f, int flag, bool unscaled = false) const;
     virtual void transform(Vector &f, const std::vector<int> &, int flag, bool unscaled = false) const;
     void zeroRotDofs(Vector &vec);
     void interp(double, const GeomState &, const GeomState &);
     void diff(const GeomState &unp, Vector &un);
     void diff1(const GeomState &un, Vector &vD, int inode);
     
     GeomState &operator=(const GeomState &); // copy geometric data -- Assume similar geomState objects

     void print();
     void printNode(int nodeNumber);

     // these functions are used by the spring corotator
     void computeGlobalRotation();
     void getGlobalRot(double R[3][3]);
     void computeCG(double cg[3]);
     void computeRotMat(double angle[3], double mat[3][3]);
     void solve(double [3][3], double [3]);
     void computeRotGradAndJac(double cg [3], 
			       double  grad[3], double jac[3][3]);
     void rotate(double mat[3][3], double vec[3]);
     void setNewmarkParameters(double _beta, double _gamma, double _alpham, double _alphaf);

     void addMultiplierNode(const std::pair<int,int> &lmpc_id, double value);
     double getMultiplier(const std::pair<int,int> &lmpc_id);
     void getMultipliers(std::map<std::pair<int,int>,double> &mu);
     void setMultiplier(const std::pair<int,int> &lmpc_id, double mu);
     void setMultipliers(std::map<std::pair<int,int>,double> &mu);

     bool getHaveRot() const { return haveRot; }
     int getNumRotationDof(int inode) const;

     const CoordSet* getCoordSet() { return X0; }
     void transformCoords(double xScaleFactor, double yScaleFactor, double zScaleFactor);
     void setNewCoords(const Vector &X);
};

inline int
GeomState::getNumRotationDof(int i) const
{
  return int(loc[i][3] >= 0) + int(loc[i][4] >= 0) + int(loc[i][5] >= 0);
}

#endif
