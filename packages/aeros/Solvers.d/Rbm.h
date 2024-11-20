#ifndef _RBM_H_
#define _RBM_H_

#include <Utils.d/MyComplex.h>
#include <Math.d/matrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/Connectivity.h>
#include <Element.d/Element.h>
#include <Corotational.d/GeomState.h>
#include <Feti.d/DistrVectorSet.h>
#include <iostream>
#include <set>

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
typedef GenVector<DComplex> ComplexVector;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
class IntFullM;
template <class Scalar> class SubLMPCons;
template <class Scalar> class GenDistrVectorSet;
typedef GenDistrVectorSet<double> DistrVectorSet;

class Rbm 
{
  Vector  *grbm;
  int      ngrbm;       	// # of geometric rigid body modes
  Vector **allRbm;
  double   tolgrb;		// tolerance
  int     *numBC;               // number of BC per component
  int     *nRbmPerComp;		// number of rbm per component (array)
  int     *numDofPerComp;	// number of dof per component (array)
  int     *firstDofOfComp;	// first dof of a component    (array)
  double **xyzRot;              // xyz location about which rotation modes are calculated
  int      nComponents;		// number of components
  int      numUncon; 		// # of unconstrained dofs
  const ConstrainedDSA* c_dsa;
  const DofSetArray*   dsa;
  const compStruct*   comp;
  IntFullM *cornerModes;
  ComplexVector *cgrbm;  // complex grbm, temporary fix to avoid templating class

  // Nonlinear members to reuse memory instead of reallocating 
  // at each iteration
  FullM U;
  FullM *Rmat;
  FullM *Amat;
  int myMemory;

public:
  FullM R;
  FullM *Zstar;
  FullM *Rc;
  FullM *Zmpc;

	Rbm() { init(); }
	Rbm(Vector *zem, int numrbm, int numUncon, int myMemory = 0);
	Rbm(ComplexVector *zem, int numrbm, int numUncon, int myMemory = 0);
	Rbm(DofSetArray *dsa, ConstrainedDSA *c_dsa);
	Rbm(DofSetArray *dsa, ConstrainedDSA *c_dsa, const CoordSet &cs, double tolgrb,
	    compStruct &components, IntFullM *fm = 0);
	Rbm(DofSetArray *dsa, const ConstrainedDSA *c_dsa, const CoordSet &cs, double tolgrb,
	    const compStruct &components, int numMPC, ResizeArray<LMPCons *> &mpc, IntFullM *fm = 0);
	Rbm(const DofSetArray *_dsa, const ConstrainedDSA *_c_dsa, const CoordSet &cs,
	    double _tolgrb, double *centroid,
	    const std::vector<int> &cornerNodes, int numCRN, int numCRNdof, const std::vector<DofSet> &cornerDofs,
	    int numMPC, const std::vector<std::unique_ptr<SubLMPCons<double> > > &mpc);
	Rbm(const DofSetArray *_dsa, const ConstrainedDSA *_c_dsa, const CoordSet &cs,
	    double _tolgrb, double *centroid,
	    const std::vector<int> &cornerNodes, int numCRN, int numCRNdof,
	    const std::vector<DofSet> &cornerDofs, int numMPC,
	    const std::vector<std::unique_ptr<SubLMPCons<DComplex> > > &mpc)
	{ std::cerr << "Rbm(...) not implemented for complex LMPCs \n"; }
	~Rbm();
  
  void computeRbms(const CoordSet &cs, double *centroid, const std::vector<int> &cornerNodes,
                   int numCRN, int numCRNdof, const std::vector<DofSet> &cornerDofs,
                   int numMPC, const std::vector<std::unique_ptr<SubLMPCons<double> > >&mpc);
 
  void reBuildGeometricRbms(GeomState *gs);
  void initialize(int nComponents, int *numCol, int numUncon);
  void computeRbms(const CoordSet &cs);
  void computeRbms(const CoordSet& cs,  int numMPC, ResizeArray<LMPCons *> &mpc);
  void clean_up();

	int numRBM() const          { return ngrbm;               }
	int numRBM(int num) const   { return nRbmPerComp[num];    }
	int numComponents() const   { return nComponents;         }
	int firstDof(int num) const { return firstDofOfComp[num]; }
	int numDof(int num) const   { return numDofPerComp[num];  }
	int numDof() const          { return numUncon;            }
	void getxyzRot(int num, double* ans) {
		ans[0] = xyzRot[num][0]; ans[1] = xyzRot[num][1]; ans[2] = xyzRot[num][2];
	}
  Vector getGrbm(int i) { return grbm[i]; }
  void setGrbm(Vector *_grbm) { grbm = _grbm; }
  void print();
  template<class Scalar> void getRBMs(Scalar *, bool transpose = false); 
  void getRBMs(double *, std::set<int> &rbmFilters);
  void getRBMs(double *, int nc, int *dofs, int _numG=-1, int offset=0);
  void getRBMs(DComplex *, int nc, int *dofs, int _numG=-1, int offset=0);
  void getScaledRBMs(double *, int nc, int *dofs, double *scaling, int _numG=-1, int offset=0);
  void getRBMs(Vector* rigidBodyModes);
  void getRBMs(VectorSet& rigidBodyModes);
  void getRBMs(DistrVectorSet& rigidBodyModes) { std::cerr << "Rbm::getRBMs(DistrVectorSet& rigidBodyModes) is not implemented\n"; }

  void singularValueDecomposition(FullM &A, FullM &U, int ncol, int nrow,
                                  double max_value, int &numgrbm, int &rank);

 private:
  void init();
};


template<>
void
Rbm::getRBMs(double *rigidBodyModes, bool);

template<>
void
Rbm::getRBMs(DComplex *rigidBodyModes, bool);

#endif
