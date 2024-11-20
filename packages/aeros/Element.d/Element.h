#ifndef	_ELEMENT_H_
#define	_ELEMENT_H_

#include <Math.d/ComplexD.h>
#include <Utils.d/BlockAlloc.h>
#include <Utils.d/dofset.h>
#include <Utils.d/GlobalToLocalMap.h>
#include <Utils.d/MFTT.h>
#include <Utils.d/Conwep.d/BlastLoading.h>
#include <iostream>
#include <vector>
#include <cstddef>
#include <complex>
#include <set>
#include <map>
#include <Eigen/Dense>

class Corotator;
class State;
class PolygonSet;
class InterpPoint;
class NLMaterial;
class LMPCons;
class GeomState;
class NFrameData;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
typedef GenVector<DComplex> ComplexVector;
template <class Scalar> class GenFullM;
typedef GenFullM<double> FullM;
typedef GenFullM<DComplex> FullMC;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
typedef GenFullSquareMatrix<DComplex> FullSquareMatrixC;
typedef GlobalToLocalMap EleRenumMap;

// Boundary Condition Structure
struct BCond {
	int nnum;   // node number
	int dofnum; // dof number (0-6)
	double val; // value of bc
	enum BCType { Forces=0, Flux, Convection, Radiation, Hneu, Atdneu, Usdf, Actuators,
		Displacements, Temperatures, Hdir, Atddir, Usdd, Pdir, Hefrs,
		Idisplacements, Idisp6, Itemperatures, Ivelocities, Iaccelerations,
		Sensors, Undefined, Lmpc, PointPointDistance, PointLineDistance, PointPlaneDistance,
		Etemperatures } type;
	int loadsetid;
	enum MomentType { Axial=0, Rotational, Follower } mtype;
	void setData(int _nnum, int _dofnum, double _val, BCType _type = Undefined, int _loadsetid = -1,
	             MomentType _mtype = Axial) { nnum = _nnum; dofnum = _dofnum; val = _val; type = _type; loadsetid = _loadsetid; mtype = _mtype; }
};

// Complex Boundary Condition Structure
struct ComplexBCond {
	int    nnum;   // node number
	int    dofnum; // dof number
	double reval;  // real value of bc
	double imval;  // imaginary value of bc
	int loadsetid;
	void setData(int _nnum, int _dofnum, double _reval, double _imval, int _loadsetid = -1)
	{ nnum = _nnum; dofnum = _dofnum; reval = _reval; imval = _imval; loadsetid = _loadsetid; };
};

// Pressure Boundary Condition Structure
struct PressureBCond {
	union {
		int elnum;
		int surfid;
	};
	double val;
	int loadsetid;
	bool conwepswitch;
	MFTTData *mftt;
	BlastLoading::BlastData *conwep;
	double loadfactor;
	int face;
	void setData(int _elnum, double _val, int _loadsetid, bool _conwepswitch, MFTTData *_mftt = NULL, BlastLoading::BlastData *_conwep = NULL)
	{ elnum = _elnum; val = _val; loadsetid = _loadsetid; conwepswitch = _conwepswitch; mftt = _mftt; conwep = _conwep; loadfactor = 1; face = -1; }
};

struct NumExc {
	int gN, n1, n2, n3;
	NumExc(int g, int a1, int a2, int a3)
	{ gN = g; n1 = a1; n2 = a2; n3 = a3; }
};

typedef double EFrame[3][3];

// ****************************************************************
// * Class StructProp contains the defined structural properties  *
// ****************************************************************

// contains material and geometrical properties of elements

struct PMLProps {
	int PMLtype;
	double gamma, Rx, Ry, Rz, Sx, Sy, Sz;
};

struct FreeplayProps {
	double k, ll, ul, lz, dz, uz;
};

class StructProp {
public:
	union {
		double  A;      // Cross-sectional area
		double kx;
		double cx;      // x-component of discrete mass offset
	};
	union {
		double	E;      // Elastic modulus
		double d0;      // Initial stiffness of nonlin element
		double ky;
		double cy;      // y-component of discrete mass offset
	};
	union {
		double	nu; 	// Poisson's ratio
		double a;	// shear-t contribution ratio
		double kz;
		double lambda;  // damage control
		double omega;
		double cz;      // z-component of discrete mass offset
	};
	union {
		double  rho; 	// Mass density per unit volume
	};
	union {
		double  eh;	// Element thickness
		double  C1;	// Non-uniform torsion constant
		double b;       // shear-s contribution ratio
		double xo;      // Plastic parameter
		double phase;
		double c1;      // 1st parameter of an elementary prescribed motion function
	};
	union {
		double  Ixx;	// Cross-sectional moment of inertia about local x-axis
		double  ss;     // speed of sound
		double  c2;     // 2nd parameter of an elementary prescribed motion function
	};
	union {
		double  Iyy;	// Cross-sectional moment of inertia about local y-axis
		double  c3;     // 3rd parameter of an elementary prescribed motion function
	};
	union {
		double  Izz;	// Cross-sectional moment of inertia about local z-axis
		double  c4;     // 4th parameter of an elementary prescribed motion function
	};
	union {
		double c;       // Thermal convection coefficient
		double alphaY;  // Shear deflection constant associated to Iyy
		double sigmax;
		double sigE;    // Elastic limit
		double ft;      // Tensile strength
		double eps;     // Radiation: Emissivity of the body (for a black body, eps=1)
		double amplitude;
	};
	union {
		double k;       // Heat conduction coefficient
		double alphaZ;  // Shear deflection constant associated to Izz
		double v2;      // Fracture Energy
		double fc;      // Compressive strength
		double Ep;      // Hardening modulus
		double offset;
	};
	union {
		double Q;       // Specific heat coefficient
		double Mx;      // global x-component of applied moment
		double Fx;      // global x-component of applied force
		double F0;      // magnitude of an applied force
	};
	union {
		double W;       // Thermal expansion coefficient
		double My;      // global y-component of applied moment
		double Fy;      // global y-component of applied force
	};
	union {
		double P;       // Perimeter/circumference of thermal truss/beams
		double Mz;      // global z-component of applied moment
		double Fz;      // global z-component of applied force
	};
	union {
		double Ta;      // Ambient temperature
		double Tr;      // Temperature of the enclosure receiving the radiation
	};
	double sigma;   // Stefan's constant (5.6704e-8 in SI)
	union {
		double ymin;    // minimum height (< 0) for cross section of beam (local y-direction)
		double Ixy;     // product of inertia
	};
	union {
		double ymax;    // maximum height (> 0) for cross section of beam (local y-direction)
		double Iyz;     // product of inertia
	};
	union {
		double zmin;    // minimum height (< 0) for cross section of beam (local z-direction)
		double Ixz;     // product of inertia
	};
	double zmax;    // maximum height (> 0) for cross section of beam (local z-direction)

	double betaDamp; // Rayleigh stiffness damping coefficient
	double alphaDamp; // Rayleigh mass damping coefficient
	double etaDamp; // Structural damping coefficient
	int etaDampTable; // Structural damping coefficient table id
	double eta_mu, deta_mu, eta_E, deta_E, mu0, dmu, E0, dE ;
	int rubDampTable;

	double kappaHelm; // wave number for Helmholtz proplem
	double kappaHelmImag; // imaginary part of the wavenumber for
	// Helmholtz problem
	complex<double> soundSpeed;

	bool lagrangeMult; // whether or not to use lagrange multiplier for mpc type elements
	double penalty, initialPenalty; // penalty parameter for mpc type elements
	int funtype; // prescribed motion function type: 0 for sinusoidal, 1 for bounded ramp
	double B, C;
	int relop; // 0: equality (==), 1: inequality (<=)
	int constraint_hess;
	double constraint_hess_eps;
	FreeplayProps freeplay[3];
	enum PropType { Undefined=0, Fluid, Fabric, Thermal, Constraint } type;
	double k1, k2, k3;
	MFTTData *ymtt, *ctett;

	// Fabric Material Options
	int F_op; // Fabric Material Option
	double F_Uc; // Critical Stretch of the Fibrils
	double F_Uf; // Failure Stretch of the Yarn
	union {
		double F_h; // Initial Height of the Yarn
		double F_lam_y; // Lambda when E is zero (linear relationship)
	};
	union {
		double F_d; // Standard Deviation for the Fibril Inclination
		double F_Estd; // Standard Deviation in E
	};
	double F_dlambda; // Standard Deviation in Lambda
	int F_np; // Number of Points in Curve Fit for Lambda
	int F_Nf; // Number of Fibrils in a Yarn
	int Seed; // Seed for Random Number Generator

	PMLProps fp;

	/** the W and E coefficient might encode integer values when they're negative
	 * (see manual for this). Heavily templated Sower needs a temporary storage that's addressable.
	 * For some reason it wont let us do that any other way than this. I'm hoping to improve this. TG
	 */
	int __SOWER_TMP;

	StructProp() { E = 0.0; A = 0.0; nu = 0.0; rho = 1.0; eh = 0.0; Ixx = 0.0;
		Iyy = 0.0; Izz = 0.0; c = 0.0; k = 0.0; Q = 0.0; W = 0.0;
		P = 0.0; Ta = 0.0; sigma = 0.0;
		kappaHelm = 0.0; kappaHelmImag = 0.0; fp.PMLtype = 0;
		soundSpeed = 1.0; alphaDamp = 0.0; betaDamp = 0.0;
		etaDamp = 0.0; etaDampTable = -1;
		ymin = 0.0; ymax = 0.0;
		zmin = 0.0; zmax = 0.0;
		lagrangeMult = true; penalty = 0.0; initialPenalty = 0.0;
		B = 1.0; C = 0.0; relop = 0; type = Undefined; funtype = 0;
		k1 = 0; k2 = 0; k3 = 0; constraint_hess = 1; constraint_hess_eps = 0.0;
		freeplay[0].dz = freeplay[1].dz = freeplay[2].dz = 1; ymtt = NULL; ctett = NULL;
		eta_mu=deta_mu=eta_E=deta_E=mu0=dmu=E0=dE = 0.0;
		rubDampTable = -1;
	}

};

typedef std::map<int, StructProp, std::less<int>, block_allocator<std::pair<const int, StructProp> > > SPropContainer;

// ****************************************************************
// *                                                              *
// *     class Node: Keeps the coordinates of a node              *
// *                                                              *
// ****************************************************************

class Node {
public:
	// Constructors
	Node() {}
	Node(const double *xyz, int _cp=0, int _cd=0) { x = xyz[0]; y = xyz[1]; z = xyz[2]; cp = _cp, cd = _cd; }
	Node(const Node &node) = default;
	~Node() {};
	double distance2(const Node& node) const
	{ return (node.x-x)*(node.x-x)+(node.y-y)*(node.y-y)+(node.z-z)*(node.z-z); }
	double distance(const Node &node) const { return sqrt(distance2(node)); }

	// Coordinates
	double 	x;
	double 	y;
	double 	z;

	// Frames
	int cp;
	int cd;
};

// ****************************************************************
// *                        WARNING                               *
// *       Nodes in the CoordSet are from 0 and not from 1        *
// *       It is a good idea to subtract 1 from node input        *
// *       from the user at any time a node is read in rather     *
// *       then keeping it from 1 and then subtracting 1 here     *
// *       and there....                                          *
// ****************************************************************

class CoordSet {
	int nmax;
	int last;
	std::vector<Node *> nodes;
	BlockAlloc ba;

public:
	// Constructors
	CoordSet(int = 256);
	CoordSet(const CoordSet&) = delete;
	CoordSet(CoordSet &&) = default;

	// Destructor
	~CoordSet();

	// Assignment operator
	CoordSet & operator = (const CoordSet & other);

	/// \brief Obtain the highest node number plus 1.
	int size() const;
	/** \brief Insert a new node to the set.
	 *
	 * @param n Node index.
	 * @param xyz 3D coordinates of the node to add.
	 * @param cp Frame index.
	 * @param cd Frame index
	 */
	void  nodeadd(int n, const double *xyz, int cp = 0, int cd = 0);
	/// \brief Inset a new node.
	void  nodeadd(int n, Node &node);
	Node &getNode(int n);
	const Node &getNode(int n) const;
	void getCoordinates(const int *nn, int numNodes,
	                    double *xx, double *yy, double *zz) const;

	Node * operator[] (int i) const { return (i >= nmax) ? 0 : nodes[i]; }
	Node *& operator[] (int i);

	/// \brief Count the actual number of defined nodes, skipping over gaps in numbering.
	int nnz() const;

	std::pair<Eigen::Matrix<double,3,1>, int> computeSums() const;
	NFrameData * dofFrame(int i) const;
};


struct PrioInfo {
	int priority; // the level of priority (undefined in isReady == 0)
	bool isReady; // whether this element is ready to be inserted
};

class MultiFront;

/****************************************************************
 *                 Abstract Element class                       *
 *                                                              *
 * Class Element defines functions for finite elements.         *
 * Each element has a structural property and a pressure        *
 * associated with it. Each element defines it's own dofs, sets *
 * it's own node numbers, calculates stiffness, calculates mass *
 * von mises stress, internal force and displacements. Where    *
 * these functions are written for each type of element and if  *
 * a function is not defined for an element type, it will       *
 * default to the appropriate function in Element.C             *
 * i.e. an element with a zero mass matrix, will just return a  *
 * matrix of the appropriate size containing zeroes.            *
 ****************************************************************/

class Element {
public:
	enum Category { Structural=0, Acoustic, Thermal, Fluid, Undefined };
protected:
	StructProp *prop;	// structural properties for this element
	bool myProp;
	int glNum, subNum, stateOffset;

	// AN: attribute id associated with the element.
	int attributeId;

	void lumpMatrix(FullSquareMatrix&, std::vector<double>&) const;
public:
	Element() { prop = 0; myProp = false; };
	virtual ~Element() { if(myProp && prop) delete prop; }
	const StructProp * getProperty() const { return prop; }
	StructProp * getProperty() { return prop; }

	virtual Element *clone() { return nullptr; }
	virtual void renum(const int *)=0;
	virtual void renum(EleRenumMap& m)=0;

	virtual void setProp(StructProp *p, bool _myProp = false) {
		if(myProp && prop)  {
			delete prop;
			prop=0;
		}
		prop = p; myProp = _myProp;
	}

	// By default ignore any element pressure
	virtual void setPressure(PressureBCond *pbc) {}
	virtual PressureBCond* getPressure() { return NULL; }

	// By default ignore any element preload
	virtual void setPreLoad(std::vector<double> &load) {}
	virtual std::vector<double> getPreLoad() { return std::vector<double>(0); }

	virtual void setGlNum(int gn, int sn=0) { glNum = gn; subNum = sn; }
	int getGlNum() const { return glNum; }

	// By default, an element has no frame
	virtual void setFrame(EFrame *) {}
	virtual const EFrame *getFrame() const { return NULL; }
	virtual void buildFrame(CoordSet&) {}

	virtual void setOffset(double*) {}
	virtual void setCompositeData(int _type, int nlays, double *lData,
	                              double *coefs, double *frame);
	virtual double * setCompositeData2(int _type, int nlays, double *lData,
	                                   double *coefs, CoordSet &cs, double theta);
	virtual void getCFrame(CoordSet& cs, double cFrame[3][3]) const;

	virtual FullSquareMatrix stiffness(const CoordSet& cs,double *k,int flg=1) const;
	virtual void getStiffnessThicknessSensitivity(CoordSet& cs,FullSquareMatrix &dStiffdThick, int flg=1);
	virtual void getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs);
	virtual FullSquareMatrix massMatrix(const CoordSet& cs,double *m,int cmflg=1) const;
	virtual FullSquareMatrix imStiffness(CoordSet& cs,double *k,int flg=1);
	FullSquareMatrix massMatrix(const CoordSet& cs, double* m, double mratio) const;
	virtual FullSquareMatrixC stiffness(const CoordSet&, complex<double> *d) const;
	virtual FullSquareMatrixC massMatrix(const CoordSet&, complex<double> *d) const;
	virtual void aRubberStiffnessDerivs(CoordSet&, complex<double> *d, int n,
	                                    double omega);

	virtual FullSquareMatrix dampingMatrix(CoordSet& cs,double *m,int cmflg=1);

	virtual double getMass(const CoordSet&) const { return 0; }
	virtual double getMassThicknessSensitivity(CoordSet&) { return 0; }
	virtual double weight(CoordSet&, double *);
	virtual double getWeightThicknessSensitivity(CoordSet& cs, double *gravityAcceleration);
	virtual void getWeightNodalCoordinateSensitivity(Vector& dwdx, CoordSet& cs, double *gravityAcceleration);
	virtual double getDCmass(CoordSet &,Vector &, double&) { return 0; }

	virtual void getGravityForce(CoordSet&,double *gravity,Vector &force,
	                             int gravflg, GeomState *gs=0);

	virtual void getGravityForceThicknessSensitivity(CoordSet&,double *gravity, Vector &force,
	                                                 int gravflg, GeomState *gs=0);

	virtual void getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration,
	                                                       GenFullM<double> &dGfdx, int gravflg, GeomState *geomState = 0);

	virtual void   getThermalForce(CoordSet& cs, Vector &ndT, Vector &force,
	                               int glflag, GeomState *gs=0);
	virtual void   getThermalForceAdj(CoordSet& cs, Vector &ndT, Vector &force,
	                                  int glflag);

	virtual void   getIntrnForce(Vector &elForce, CoordSet& cs,
	                             double *elDisp, int Index, double *ndTemps);

	virtual void   getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
	                           Vector &elDisp, int strInd, int surface=0,
	                           double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum = 1);

	virtual void   getVonMises(ComplexVector &stress, Vector &weight, CoordSet &cs,
	                           ComplexVector &elDisp, int strInd, int surface=0,
	                           double *ndTemps=0, double ylayer=0.0, double zlayer=0.0, int avgnum = 1);

	virtual void   getVonMisesInt(CoordSet &,Vector &,double &,double &, int,
	                              double &,double &, double* dT=0 );

	virtual void getVonMisesThicknessSensitivity(Vector &dStdThick, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
	                                             double * = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);

	virtual void getVonMisesThicknessSensitivity(ComplexVector &dStdThick, ComplexVector &weight, CoordSet &cs, ComplexVector &elDisp, int strInd, int surface,
	                                             double * = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);

	virtual void getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *,
	                                                CoordSet &cs, Vector &elDisp, int strInd, int surface,
	                                                double * = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);

	virtual void getVonMisesDisplacementSensitivity(GenFullM<DComplex> &dStdDisp, ComplexVector &weight, GenFullM<DComplex> *,
	                                                CoordSet &cs, ComplexVector &elDisp, int strInd, int surface,
	                                                double * = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);

	virtual void getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
	                                                   double * = 0, int avgnum = 1, double ylayer = 0, double zlayer = 0);

	virtual void   getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
	                            Vector &elDisp, int strInd, int surface=0,
	                            double *ndTemps=0);

	virtual void   getAllStress(FullMC &stress, Vector &weight, CoordSet &cs,
	                            ComplexVector &elDisp, int strInd, int surface=0,
	                            double *ndTemps=0);

	virtual void   computeHeatFluxes(Vector& heatflux, CoordSet &cs,
	                                 Vector& elTemp, int hflInd);

	virtual void  trussHeatFluxes(double &trussflux, CoordSet &cs,
	                              Vector& elTemp, int hflInd) {}

	virtual void   computeSloshDisp(Vector& fluidDispSlosh, CoordSet &cs,
	                                Vector& elPotSlosh, int hflInd);

	virtual void   computeSloshDispAll(Vector& fluidDispSlosh, CoordSet &cs,
	                                   Vector& elPotSlosh);

	virtual void   computeDisp(CoordSet&, State&, const InterpPoint &,
	                           double*, GeomState *gs=0) {}

	virtual void   getFlLoad(CoordSet &, const InterpPoint &,
	                         double *, double *, GeomState *gs=0) {}

	virtual void   computeTemp(CoordSet&, State &, double[2],
	                           double*)  {}
	virtual void   getFlFlux(double[2], double *, double *) {}

	virtual void   markDofs(DofSetArray &) const = 0;
	virtual int*   dofs(DofSetArray &, int *p=0) const =0;
	virtual int    numDofs() const =0;

	virtual int    numNodes() const = 0;
	virtual int*   nodes(int * = 0) const = 0;
	virtual int*   allNodes(int *x = 0) const { return nodes(x); }

	virtual Corotator *getCorotator(CoordSet &, double *, int = 2, int = 2);

	virtual int getTopNumber() const;
	virtual int numTopNodes() const { return numNodes() - numInternalNodes(); }   // this is the number of nodes printed in the top file
	// can make it different to numNodes for elements that aren't
	// supported by xpost eg RigidSolid6Dof
	virtual void computePressureForce(CoordSet& cs,Vector& elPressureForce,
	                                  GeomState *gs = 0, int cflg = 0, double t = 0);
	virtual double * getMidPoint(CoordSet &)  { return 0; }
	/* toxa: use this midPoint instead */
	//virtual void getMidPoint(CoordSet &, double* result)  { assert(0); }
	virtual double * getCompositeData(int )   { return 0; }
	virtual double * getCompositeFrame()      { return 0; }

	virtual int      getCompositeLayer()      { return 0; }

	virtual int  dim() const { return 3; }

	virtual void addFaces(PolygonSet *pset);

	virtual void setMaterial(NLMaterial *);

	// AN: attach attribute id to the element.
	void setElementAttribute(int a) { attributeId = a; };
        int getElementAttribute() { return attributeId; };

	virtual int numInternalNodes() const { return 0; }

	virtual void setInternalNodes(int *) {}

	virtual bool isSafe() const { return true; }

	// from DEC
	// TOP/DOMDEC Element Functions
	virtual int facelist(PolygonSet &, int * = 0) { return 0; }

	// DEC : Routines for the decomposer
	// isStart indicates if an element is suitable to
	// be a germination center for a subdomain (bars are not)
	virtual bool isStart() { return true; }

	virtual bool isSpring() const { return false; }

	virtual bool isMass() { return false; }

	virtual bool hasRot() const { return false; }
	virtual PrioInfo examine(int sub, MultiFront *mf);
	virtual double weight() const;
	virtual double trueWeight() const;

	void getCG(CoordSet &cset, double &xcg, double &ycg, double &zcg);
	virtual int nDecFaces() const { return 0; }
	virtual int getDecFace(int iFace, int *fn) { return 0; }
	// END FROM DEC

	// this need to be defined for 6 node tri shell & 8 node quad shell
	virtual bool isRotMidSideNode(int iNode) { return false; }

	// the following function is used for element pressure with the face keyword
	// and differs from getDecFace in that (a) midside nodes are not ignored, and
	// (b) for 3D solid element all face normals are outward pointing.
	virtual int getFace(int iFace, int *fn);

	virtual bool hasDamping() final { return false; }
	bool isFluidElement();
	virtual bool isSommerElement() const { return false; }
	virtual bool isRadiationElement() { return false; }
	virtual bool isSloshingElement() { return false; }
	virtual bool isMpcElement() { return false; }
	virtual bool isFsiElement() { return false; }
	virtual bool isHEVFluidElement() { return false; }
	virtual int  fsiFluidNode() { return -1; }
	virtual int  fsiStrutNode() { return -1; }
	//virtual bool isRigidMpcElement(const DofSet & = DofSet::nullDofset, bool forAllNodes=false) { return false; }
	virtual bool isRigidElement() const { return false; }
	virtual int getNumMPCs() { return 0; }
	virtual LMPCons** getMPCs() { return 0; }

	virtual double helmCoef() { return prop ? prop->kappaHelm * prop->kappaHelm : 0; }
	virtual complex<double> helmCoefC() { return prop ?
	                                             complex<double>(prop->kappaHelm,prop->kappaHelmImag)*
	                                             complex<double>(prop->kappaHelm,prop->kappaHelmImag) : 0;
	}
	virtual double helmCoef(double omega) {
		return prop ? omega*omega*real(prop->soundSpeed)*real(prop->soundSpeed) : 0;
	}
	virtual complex<double> helmCoefC(double omega) {
		return prop ? (omega*omega)*prop->soundSpeed*prop->soundSpeed : 0;
	}

	virtual bool isShell() const { return false; }

	virtual bool isConstraintElement() { return (isRigidElement() || isMpcElement() || isFsiElement()); }
	virtual bool isConstraintElementIeq() { return (isMpcElement() && prop->relop != 0); }
	virtual bool isFreeplayElement() const { return false; }
	virtual bool isPhantomElement() { return (!(prop || isConstraintElement() || isSommerElement())); }

	virtual int getElementType() const = 0;
	void setElementType(int _type) { if( _type != getElementType() ) throw "Inconsistency"; }

	// friend class Domain;
	friend class SuperElement;

	// complex interface for the Element
	virtual bool isComplex() { return false; }
	virtual FullSquareMatrixC complexStiffness(CoordSet& cs, DComplex *k,int flg=1);
	virtual FullSquareMatrixC complexDampingMatrix(CoordSet& cs, DComplex* c, int flg=1);
	virtual FullSquareMatrixC complexMassMatrix(CoordSet& cs, DComplex* m, double mratio);

	virtual Category getCategory() const = 0;//{ return category; }
	void setCategory(Category _category) {
		if(_category != getCategory())
			throw "Bad category";
	}
	bool isDamped() { return (getCategory() != Thermal && !isSpring()) ? (prop && (prop->alphaDamp != 0.0 || prop->betaDamp != 0.0)) : false; }
	bool isSDamped() { return (getCategory() != Thermal && !isSpring()) ? (prop && (prop->etaDamp != 0.0 )) : false; }

	virtual int getMassType() const { return 1; }  // 0: lumped, 1: consistent, 2: both
	// notes: (a) if getMassType returns 0 then lumped gravity force will always be used for dynamics
	//        (b) is getMassType returns 1 then lumping is done using diagonal scaling if required (default)
	virtual void writeHistory(int) const {} // Only BelytschkoTsayShell overrides this method.
	virtual void readHistory(int) {}
	virtual int numStates() { return 0; }
	virtual void initStates(double *) {}
	virtual void setStateOffset(int _stateOffset) { stateOffset = _stateOffset; }

	virtual double computeStabilityTimeStep(FullSquareMatrix &K, FullSquareMatrix &M, CoordSet &cs,
	                                        GeomState *gs, double stable_tol, int stable_maxit);
};

// *****************************************************************
// *                        WARNING                                *
// *       The same remark as for node is valid for elements       *
// *****************************************************************

template <typename T>
class SetAccess;

class Elemset
{
protected:
	Element **elem;
	int emax;
	BlockAlloc ba;
	bool myData;
	int dampingFlag;
	std::vector<std::pair<int,int> > etypes;
public:
	Elemset(int = 256);
	virtual ~Elemset() { deleteElems(); }
	Elemset(Elemset &globalset, int numlocal, int *localToGlobalMap);
	int size() const { return emax; }
	int last() const;
	Element *operator[] (int i) const { return elem[i]; }
	Element *& operator[] (int i);
	void elemadd(int num, Element *);
	void elemadd(int num, int type, int nnodes, int *nodes);
	void mpcelemadd(int num, LMPCons *mpc, bool nlflag = false);
	void fsielemadd(int num, LMPCons *fsi);
	void setEmax(int max)  { emax = max; }
	void list();
	void deleteElems();
	void remove(int num) { elem[num] = 0; }
	void removeAll();
	void setMyData(bool _myData) { myData = _myData; }
	bool hasDamping();
	void collapseRigid6(std::set<int> &);
	void deleteElem(int i);
	void setWeights();

	SetAccess<Elemset> asSet() const;
};

template<>
class SetAccess<Elemset> {
public:
	SetAccess(const Elemset &set) : set(set) {}
	int size() const { return set.last(); }; //<! returns the number of members of the set
	int numNodes(int i) const { return (set[i] ? set[i]->numNodes() : 0); } //<! returns the number of targets for member i
	void nodes(int i, int *nd) const { if(set[i]) set[i]->nodes(nd); }; //<! copies into nd the targets for member i
private:
	const Elemset& set;
};

inline
SetAccess<Elemset>
Elemset::asSet() const {
	return SetAccess<Elemset>(*this);
}

class ElementFactory
{
public:
	ElementFactory() {}
	virtual ~ElementFactory() {}
	virtual Element* elemadd(int num, int type, int nnodes, int *nodes, BlockAlloc& ba);
};

#endif
