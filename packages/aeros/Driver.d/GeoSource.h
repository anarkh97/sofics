#ifndef _GEO_SOURCE_H_
#define _GEO_SOURCE_H_

#include <cassert>
#include <cstring>
#include <string>
#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <set>

#include <Element.d/Element.h>
#include <Utils.d/OutputInfo.h>
#include <Math.d/DistVector.h>
#include <Driver.d/Attrib.h>
#include <Driver.d/EFrameData.h>
#include <Driver.d/OffsetData.h>
#include <Control.d/ControlInfo.h>
#include <Utils.d/SparseConnectivityT.h>

#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif

using namespace NewVec;
#define PI 3.14159265358979

class Connectivity;
class Communicator;
template <class Scalar> class GenSubDomain;
typedef GenSubDomain<double> SubDomain;
class Domain;
class ControlInterface;
class MatrixTimer;
class LMPCons;
class Decomposition;
class NLMaterial;
class ControlInfo;
class CoordSet;
class OutputInfo;
class LayMat;
class Attrib;
class CoefData;
class LayInfo;
struct Group;
struct AttributeToElement;
class DistrGeomState;
class ControlLawInfo;
struct DispNode;

enum {SXX=0,SYY=1,SZZ=2,SXY= 3,SYZ= 4,SXZ= 5,VON=6,
	EXX=7,EYY=8,EZZ=9,EXY=10,EYZ=11,EXZ=12,STRAINVON=13,
	VONTOP=14,VONBOT=15,CONPRESS=16,DAMAGE=17,EQPLSTRN=18,
	BACKSXX=19,BACKSYY=20,BACKSZZ=21,BACKSXY=22,BACKSYZ=23,
	BACKSXZ=24,PLASTICEXX=25,PLASTICEYY=26,PLASTICEZZ=27,
	PLASTICEXY=28,PLASTICEYZ=29,PLASTICEXZ=30,AGGREGATEDVON=31};
enum {INX,INY,INZ,AXM,AYM,AZM};
enum {YOUNG,MDENS,THICK};
enum {PSTRESS1=0,PSTRESS2=1,PSTRESS3=2,
	PSTRAIN1=3,PSTRAIN2=4,PSTRAIN3=5,
	PSTRESS1DIREC=6,PSTRESS2DIREC=7,PSTRESS3DIREC=8,
	PSTRAIN1DIREC=9,PSTRAIN2DIREC=10,PSTRAIN3DIREC=11};
// Controls the output of heat fluxes or grad(Temp)
enum {HFLX=0, HFLY=1, HFLZ=2, GRTX=3, GRTY=4, GRTZ=5};
enum {SLDX=0, SLDY=1, SLDZ=2};

struct MatchData
{
	int glPos;
	int elemNum;
	double xi, eta;  // natural coordinates
	bool operator < (const MatchData &v) const { return glPos < v.glPos; }
};

typedef NLMaterial *(*MatLoader)(int n, double *params);

struct DoubleList {
	int nval;
	double v[32];
};


class GeoSource {
protected:
	struct OrderPair {
		int node, pos;
		double xyz[3];
		double &operator[] (int i) { return xyz[i]; }
		bool operator<(const OrderPair &op) const { return node < op.node; }
	};

	int curCluster;

	bool decJustCalled; // if true, no need for an external DECOMPOSITION FILE
	bool exitAfterDec; // if true no need to save Subdomain info for FEM in memory

	// input file names
	const char *conName;
	const char *geoName;
	const char *decName;
	const char *mapName;
	const char *matchName;

	// output file info
	ControlInfo *cinfo;    // contains nodeset, elemset & timer file name
	ResizeArray<OutputInfo> oinfo; // all output files except sensitivity outputs
	int numOutInfo;                // number of Output requests
	int outLimit; // maximum number of frequencies, eigenvectors or timesteps per file
	int numNodalOutput;   // number of disp output for single nodes
	std::vector<int> outputNodes;     // list of node for singular output
	std::vector<int> outNodeIndex;
	std::vector<int> headLen;  	// header description length for binary output

	int numClusters;
	int maxGlobNode;
	int maxClusNode;
	int numClusNodes;
	int numClusElems;
	int nGlobNodes;
	int numNodes, numInternalNodes;
	int numConstraintElementsIeq;
	bool lmpcflag; // true if all ieq constraints are linear, otherwise false.
	CoordSet nodes;
	Elemset elemSet;
	Elemset *packedEsetFluid;
	Elemset *packedEsetConstraintElementIeq;
	int nElem, nMpcElem;
	int nElemFluid;
	int nElemAllClusters;
	int *allNumClusElems;
	int phantomFlag;

	int *elemTypeNumNodesMap;
	std::map<int, int> glToPckElems;  // glElemNum  -> packedElem Num

	// Match Data
	int *numMatchData;
	MatchData **matchData;
	int *numGapVecs;
	double (**gapVec)[3];      // gap vectors at struct/fluid interface

	ControlLawInfo *claw;      // user defined control law info

	// Material data
	int numLayInfo;
	ResizeArray<LayInfo *> layInfo;
	int numCoefData;
	ResizeArray<CoefData *> coefData;
	int numLayMat;
	ResizeArray<LayMat *> layMat;

	std::map<std::string, MatLoader> userDefMat;
	std::map<int, NLMaterial *> materials;
	std::map<int, int> matUsage;

	int numProps;
	SPropContainer sProps;
	int na;  			// number of attributes
	int namax;
	std::map<int, Attrib> attrib;
	int maxattrib;
	std::map<int, int> optg;
	std::map<int, int> mortar_attrib;

	int numEframes;
	ResizeArray<EFrameData> efd;
	std::vector<OffsetData> offsets;

	int numNframes;
	ResizeArray<NFrameData> nfd;

	int numCframes;
	ResizeArray<double *> cframes;

	int numCSframes;
	ResizeArray<EFrameData> csfd;

	int prsflg;
	int constpflg, constqflg; // consistent pressure and gravity
	typedef std::vector<PressureBCond> ElemPressureContainer;
	ElemPressureContainer eleprs;
	typedef std::vector<std::pair<int,std::vector<double> > > ElemPreloadContainer;
	ElemPreloadContainer eleprl;

	// Connectivities
	std::unique_ptr<Connectivity> clusToSub;
	std::vector<int> subToClus;
	Connectivity *subToSub;
	Connectivity *subToNode;
	// TODO Make a consistent ownership. Ownership used to be passed along, yet can be delete
	// by contact problems.
	Connectivity *subToElem;
	Connectivity *unsortedSubToElem;

	int *subToCPU;
	std::shared_ptr<Connectivity> cpuToSub;
	Connectivity *cpuToCPU;

	double mratio; // consistent-lumped matrix ratio; 1==consistent, 0==lumped

public:
	SparseConnectivityType1 *subToNode_sparse;
	SparseConnectivityType2 *nodeToSub_sparse;
	int fixedEndM; // 0: don't include fixed end moments in gravity force vector for beams and shells
	int *gl2ClSubMap;

	// BC Data

	int numTextDirichlet;
	BCond *textDBC;
	int numDirichlet;           // number of dirichlet bc
	int numDirichletFluid;      // number of dirichlet bc in fluid
	BCond *dbc;                 // set of those dirichlet bc
	BCond *dbcFluid;            // set of those dirichlet bc in fluid

	int numTextNeuman;
	BCond *textNBC;
	int numNeuman;              // number of Neuman bc
	BCond *nbc;                 // set of Neuman bc
	int numNeumanModal;
	BCond *nbcModal;

	int numIDis;                // number of initial displacements
	BCond *iDis;                // set of those initial displacements

	int numIDisModal;
	BCond *iDisModal;

	int numIDis6;               // number of initial displacements (6 column)
	BCond *iDis6;               // set of those intitial displacements

	// PITA
	std::map<int, std::pair<int, int> > timeSliceOutputFiles;
	BCond *PitaIDis6;    // Array of initial seed displacement array
	int numPitaIDis6;    // # of initial seed initial displacement conditions
	int numTSPitaIDis6;  // # of initial seed displacement vectors
	BCond *PitaIVel6;    // Array of initial seed velocity array
	int numPitaIVel6;    // # of initial seed initial displacement conditions
	int numTSPitaIVel6;  // # of initial seed initial velocity conditions

	int numIVel;                // number of initial velocities
	BCond *iVel;                // set of those initial velocities
	int numIVelModal;
	BCond *iVelModal;

	int numComplexDirichlet;
	ComplexBCond *cdbc;

	int numComplexNeuman;
	ComplexBCond *cnbc;

	bool isShift;
	double shiftV;

	int numDampedModes;   // number of modes that have damping
	BCond *modalDamping;  // the value of damping for those modes with damping reuses BCond class

	Decomposition *optDec = nullptr, *optDecCopy = nullptr;

	std::map<int, Group> group;
	std::map<int, std::set<int> > nodeGroup;
	std::map<int, std::list<int> > surfaceGroup;

	std::map<int, AttributeToElement> atoe;

	int numSurfaceDirichlet;
	BCond *surface_dbc = nullptr;
	int numSurfaceNeuman;
	BCond *surface_nbc = nullptr;
	int numSurfacePressure;
	PressureBCond *surface_pres = nullptr;
	int numSurfaceConstraint;
	BCond *surface_cfe;

public:
	bool binaryInput, binaryOutput;
	bool binaryInputControlLeft;

	GeoSource(int iniSize = 16);
	virtual ~GeoSource();

	std::vector<int> localToGlobalElementsNumber;

	// Input Read Functions
	void readCpuToSub();
	int readRanges(BinFileHandler &, int &, int (*&r)[2]);
	void readMatchInfo(BinFileHandler &, int (*)[2], int, int, int *, int);
	template<class Scalar>
	GenSubDomain<Scalar> *getSubDomain(int glSub, Domain *, int locSub = -1);
	template<class Scalar>
	GenSubDomain<Scalar> *getSubDomain(Domain *, int, Connectivity *);
	void setElemTypeMap();
	void getClusterData(BinFileHandler &);
	void applyAuxData(int *, int *, int, int);
	template<class Scalar>
	void distributeBCs(GenSubDomain<Scalar> *&, int *, int *gl2clNodeMap = 0);
	int getBC(BCond *, int, int *, BCond *&, int *gl2clNodeMap = 0);
	void augmentBC(int, BCond *, BCond *&, int &);
	int getCPUMap(FILE *f, Connectivity *, int glNumSub, int nCpu);
	void createSingleCpuToSub(int numSub);
	int getSubCtrl(BCond *, int, BCond *&, int *, int *, int *&); //bin geo
	int getSubCtrl(BCond *, int, BCond *&, int, int *&); // text geo input
	template<class Scalar>
	void distributeCtrlLaw(GenSubDomain<Scalar> *&, int *, int *);  // binary geo input
	template<class Scalar>
	void distributeCtrlLaw(GenSubDomain<Scalar> *&, int); // text geo input
	template<class Scalar>
	void distributeOutputNodes(GenSubDomain<Scalar> *&, int *, int *);
	template<class Scalar>
	void distributeOutputNodesX(GenSubDomain<Scalar> *, Connectivity *nodeToSub);

	void setGeo(char *file) { geoName = file; }
	void setDecomp(char *file) { decName = file; }
	void setMatch(char *file) { matchName = file; }
	void setCpuMap(char *file) { mapName = file; }
	void setGlob(char *file) { conName = file; }
	const char *getGlob() const { return conName; }
	void setExitAfterDec(bool exit);
	void setNumLocSub(int);
	void deleteMatchArrays(int);
	void setTextBC();
	void computeGlobalNumElements();

	// duplicate files for the time parallel method
	void duplicateFilesForPita(int, const int*);

	// decomp read functions
	std::unique_ptr<Connectivity> getDecomposition();
	void getTextDecomp(bool sowering = false);
	void getBinaryDecomp();

	// Parser support Functions
	void setControl(char *checkfile, char *nodeset, char *elemset, char *bcondset = 0);

	// Parser support Functions - Geometry
	int  addNode(int nd, double xyz[3], int cp, int cd);
	int  addElem(int en, int type, int nn, int *nodeNumbers);
	int  addMat(int, const StructProp &);
	int  addLay(int, LayInfo *);
	int  addCoefInfo(int, CoefData &);
	CoefData* getCoefData(int i) { return (i >= 0 && i < numCoefData) ? coefData[i] : NULL; }
	int  addLayMat(int m, double *);
	int  setAttrib(int n, int a, int ca = -1, int cfrm = -1, double ctheta = 0.0);
	void setMortarAttrib(int n, int a);
	int  setFrame(int, double *);
	int  setNodalFrame(int, double *, double *, int);
	int  setCSFrame(int, double *);
	int  addCFrame(int, double *);
	void setElementPressure(PressureBCond&);
	void setElementPreLoad(int, double);
	void setElementPreLoad(int, double[3]);
	void setConsistentPFlag(int);
	void setConsistentQFlag(int, int=1);
	void addOffset(OffsetData &od) { offsets.push_back(od); }

	// Parser support Functions - Boundary Conditions
	int  setDirichlet(int, BCond *);
	void convertHEVDirToHelmDir(); // added to use HEFRS and Helmholtz per Charbel's request
	int  setDirichletFluid(int, BCond *);
	int  setNeuman(int, BCond *);
	int  setNeumanModal(int, BCond *);
	int  setIDis(int, BCond *);
	int  setIDisModal(int, BCond *);
	int  setIDis6(int, BCond *);
	int  setIVel(int, BCond *);
	int  setIVelModal(int, BCond *);
	int  addSurfaceDirichlet(int, BCond *);
	int  addSurfaceNeuman(int, BCond *);
	int  addSurfacePressure(int, PressureBCond *);
	int  addSurfaceConstraint(int, BCond *);

	// PITA
	int getLocalTimeSliceCount() { return timeSliceOutputFiles.size(); }
	std::pair<int, int> getTimeSliceOutputFileIndices(int timeSliceRank);
	// Set initial seed conditions
	int setPitaIDis6(int, BCond *, int);
	int setPitaIVel6(int, BCond *, int);

	// Parser support Function - Damping
	int setModalDamping(int, BCond *);

	// Parser support Functions - "New" Material Specifications
	void addMaterial(int i, NLMaterial *m) { materials[i] = m; }
	void addMaterial(int i, const char *, DoubleList &d);
	void loadMaterial(const char *, const char *);
	void setMatUsage(int i, int t) { matUsage[i] = t; }
	// Parser support Functions - Output Functions
	void addOutput(OutputInfo &);

	// Parser support Functions - User Defined Control Functions
	void setControlFile(char *_filename);
	void setControlRoutine(char *_routinename);
	int  setSensorLocations(int, BCond *);
	int  setActuatorLocations(int, BCond *);
	int  setUsddLocation(int, BCond *);
	int  setUsdfLocation(int, BCond *);

	void transformCoords();
	void setNewCoords(std::string nodeFile);
	void checkInputs();
	void setUpData(int topFlag);

	Elemset* getPackedEsetFluid() { return packedEsetFluid; }

	// PITA: Access to initial seed conditions
	int getNumPitaIDis6() const  { return numPitaIDis6; }     // # of initial seed initial displacement conditions
	int getNumPitaIVel6() const  { return numPitaIVel6; }     // # of initial seed initial velocity conditions
	int getNumTSPitaIDis6() const { return numTSPitaIDis6; }   // # of initial seed displacement vectors
	int getNumTSPitaIVel6() const { return numTSPitaIVel6; }   // # of initial seed velocity vectors
	int getUserProvidedSeedCount() const { return std::min(numTSPitaIDis6, numTSPitaIVel6); }
	const BCond &getPitaIDis6(int i) const { return PitaIDis6[i]; }  // Initial seed displacement for each time-slice
	BCond &getPitaIDis6(int i) { return PitaIDis6[i]; }
	const BCond &getPitaIVel6(int i) const { return PitaIVel6[i]; }  // Initial seed velocity for each time-slice
	BCond &getPitaIVel6(int i) { return PitaIVel6[i]; }

	ControlInfo *getCheckFileInfo()  { return cinfo; }
	int  getNumClusNodes()  { return numClusNodes; }
	int  getNodes(CoordSet &);
	CoordSet&  GetNodes(); // HB: return the node CoordSet of GeoSource
	int  getMaxNodeNum()  { return maxGlobNode; }
	int  getNumGlobNodes()  { return nGlobNodes; }
	void  setNumGlobNodes(int n) { nGlobNodes = n; }
	int  getElems(Elemset &, int = 0, int * = 0);
	int  getNumClusters() { return numClusters; }
	const char *getCpuMapFile() const { return mapName; }
	const char *getMatchFileName() const { return matchName; }
	int  numElem() { return nElem; }
	int  numMpcElem() { return nMpcElem; }
	int  numContactSurfaceElem() { return contactSurfElems.size(); }
	int  numElemFluid() { return nElemFluid; }
	int  numNode() { return numNodes; }
	void setNumNodes(int n) { numNodes = n; }
	int  totalNumNodes() { return numNodes + numInternalNodes; }
	int  internalNumNodes() { return numInternalNodes; }
	int  getNumConstraintElementsIeq() { return numConstraintElementsIeq; }
	Elemset& getPackedEsetConstraintElementIeq() { return *packedEsetConstraintElementIeq; }
	bool getLmpcFlag() { return lmpcflag; }
	//int  getPhantomFlag()  { return phantomFlag; }
	//int  glToPack(int i) { return glToPck[i]; }
	int  glToPackElem(int i) const;
	const Connectivity &getClusToSub() const { return *clusToSub; }
	int *getSubToClus()  { return subToClus.data(); }
	const Connectivity *getSubToSub() const { return subToSub; }
	const Connectivity *getSubToElem() const { return subToElem; }
	void setSubToElem(Connectivity *ste) { subToElem = ste; }
	Connectivity *getSubToNode()  { return subToNode; }

	LayInfo *getLayerInfo(int num)  { return layInfo[num]; }
	SPropContainer &getStructProps() { return sProps; }
	const std::map<int, NLMaterial *> &getMaterialLaws() const { return materials; }
	const std::map<int, int> &getMaterialLawMapping() const { return matUsage; }
	EFrameData *getEframes()  { return efd+0; }
	EFrameData *getCSframes() { return csfd+0; }
	NFrameData *getNFrames() { return nfd+0; }
	double **getCframes()  { return cframes+0; }
	LayInfo **getLayInfo() { return layInfo+0; }
	int getNumLayInfo() { return numLayInfo; }
	int getNumEframes() { return numEframes; }
	int getNumCSframes() { return numCSframes; }
	int getNumCframes() { return numCframes; }
	int getNumNframes() { return numNframes; }
	const ElemPressureContainer &getElementPressure() const { return eleprs; }
	int pressureFlag() { return prsflg; }
	int consistentPFlag() { return constpflg; }
	int consistentQFlag() { return constqflg; }
	std::map<int, Attrib> &getAttributes()  { return attrib; }
	std::map<int, int> &getMortarAttributes() { return mortar_attrib; }
	int getNumAttributes() { return na; }

	int getMaxAttributeNum() { return maxattrib+1; }

	int getNumDirichlet() { return numDirichlet; }
	int getNumDirichletFluid() { return numDirichletFluid; }
	int getNumNeuman() { return numNeuman; }
	int getNumNeumanModal() { return numNeumanModal; }
	int getNumIDisModal() { return numIDisModal; }
	int getDirichletBC(BCond *&);
	int getDirichletBCFluid(BCond *&);
	int getTextDirichletBC(BCond *&);
	int getNeumanBC(BCond *&);
	int getNeumanBCModal(BCond *&);
	int getTextNeumanBC(BCond *&);
	int getIDis(BCond *&);
	int getIDisModal(BCond *&);
	int getIDis6(BCond *&);
	int getIVel(BCond *&);
	int getIVelModal(BCond *&);
	int getITemp(BCond *&);
	int getNumProps() { return numProps; }
	void setNumProps(int n) { numProps = n; }

	int getSurfaceDirichletBC(BCond *&);
	int getNumSurfaceDirichlet() { return numSurfaceDirichlet; }
	int getSurfaceNeumanBC(BCond *&);
	int getNumSurfaceNeuman() { return numSurfaceNeuman; }
	int getSurfacePressure(PressureBCond *&);
	int getNumSurfacePressure() { return numSurfacePressure; }
	int getSurfaceConstraint(BCond *&);
	int getNumSurfaceConstraint() { return numSurfaceConstraint; }

	int getModalDamping(BCond *&);

	int *getNumMatchData()  { return numMatchData; }
	MatchData *getMatchData(int iSub)  { return matchData[iSub]; }
	double (*getGapVecs(int iSub))[3]  { return gapVec[iSub]; }

	int getNumOutInfo()  { return numOutInfo; }
	void setOutLimit(int _outLimit) { outLimit = _outLimit; }
	int getOutLimit() { return outLimit; }
	void setNumNodalOutput();
	int getNumNodalOutput() { return numNodalOutput; }
	OutputInfo *getOutputInfo()  { return oinfo+0; }
	bool elemOutput();
	bool energiesOutput();
	bool romExtForceOutput();
	bool noOutput(int x, int ndflag = 0);

	int *getSubToCPU()  { return subToCPU; }
//	void setCpuToSub(Connectivity *c) { cpuToSub=c; }
	std::shared_ptr<Connectivity> getCpuToSub() { return cpuToSub; }
	Connectivity *getCpuTOCPU()  { return cpuToCPU; }

	ControlLawInfo *getControlLaw()  { return claw; }
//  void makeGlobalClaw(ControlLawInfo *subClaw);
	ControlInterface *getUserSuppliedFunction();

	Elemset* getElemSet(void){return(&elemSet);}

	void simpleDecomposition(int numSubdomains, bool estFlag, bool weightOutFlag,
							 bool makeTrivial, bool fsglFlag);
	void modifyDecomposition(int maxEleKeep);

	// Output Functions
	template<int bound>
	void outputNodeVectors(int, double (*)[bound], int, double time = -1.0);
	template<int bound>
	void outputNodeVectors(int, DComplex (*)[bound], int, double time = -1.0);
	template<int bound>
	void outputNodeVectors6(int, double (*)[bound], int, double time = -1.0);
	template<int bound>
	void outputNodeVectors6(int, DComplex (*)[bound], int, double time = -1.0);
	template<int bound>
	void outputNodeVectors9(int, double (*)[bound], int, double time = -1.0);
	template<int bound>
	void outputNodeVectors4(int, double (*)[bound], int, double time = -1.0);
#ifdef USE_EIGEN3
	template<class Scalar>
	void outputSensitivityScalars(int, Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *, double time = 0.0,
								  Eigen::Matrix<double, Eigen::Dynamic, 1> *dwr = 0);
	template<class Scalar>
	void outputSensitivityVectors(int, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *,
								  double time = 0.0, Eigen::Matrix<double, Eigen::Dynamic, 1> *dwr = 0);
	template<class Scalar>
	void outputSensitivityDispVectors(int, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> **,
									  double time = 0.0, int numParams = 0, int numnodes = 0);
	template<class Scalar>
	void outputSensitivityAdjointStressVectors(int, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *, Scalar *,
											   double time, int numParams, std::vector<int>,
											   Eigen::Matrix<double,Eigen::Dynamic, 1> *dwr = 0);
	template<class Scalar>
	void outputSensitivityAdjointDispVectors(int, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> **, Scalar *,
											 double time, int numParams, std::vector<DispNode>,
											 Eigen::Matrix<double, Eigen::Dynamic, 1> *dwr = 0);
	template<class Scalar>
	void outputSensitivityDispVectors(int, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *,
									  double time = 0.0, int numnodes = 0);
#endif
	void outputNodeScalars(int, double *, int, double time = -1.0);
	void outputNodeScalars(int, DComplex *, int, double time = -1.0);
	void outputEnergy(int, double, double);
	void outputEnergy(int, double, DComplex);
	void outputEnergies(int, double, double, double, double, double, double, double);
	void outputEnergies(int, double, DComplex, DComplex, DComplex, DComplex, DComplex, DComplex);
	void outputEnergyPerAttribute(int, double, double *, int);
	void outputElemVectors(int, double *, int, double time = -1.0);
	void outputElemVectors(int, DComplex *, int, double time = -1.0);
	void outputElemStress(int, double *, int, const std::vector<size_t> &offsets, double = -1.0);
	void outputElemStress(int, DComplex *, int, const std::vector<size_t> &, double = -1.0);
	void openOutputFiles(int *outNodes = 0, int *outIndex = 0, int num = 0);
	void openSensorOutputFiles();
	void openOutputFilesForPita(int sliceRank);
	void closeOutputFiles();
	void closeOutputFilesForPita(int sliceRank);
	void createBinaryOutputFile(int, int, int iter = 0);
	void computeAndCacheHeaderLength(int);
	void outputHeader(int);
	void outputRange(int, int *, int, int, int , int iter = 0);

private:
	int getHeaderDescriptionAndLength(char *, int);
	void getHeaderDescription(char *, int);

public:
	void readGlobalBinaryData();
	void computeClusterInfo(int glSub, Connectivity *subToNode = NULL);

	void writeDistributedInputFiles(int nCluster, Domain*, int nCpu);
#ifdef SOWER_SURFS
	void readDistributedSurfs(int subNum);
#endif
	/** \brief Read aand build a list of subdomains.
	 *
	 * @tparam Scalar Scalar type for the subdomain
	 * @param glSubIndices List of global indices of the subdomains to construct
	 * @return A vector of pointers to the constructed subdoains.
	 */
	template<class Scalar>
	std::vector<GenSubDomain<Scalar> *> readDistributedInputFiles(gsl::span<const int> glSubIndices);

	// Output functions, implemented in Driver.d/BinaryOutput.C and Driver.d/BinaryOutputInclude.C
	void writeNodeScalarToFile(double *data, int numData, int glSub, int offset, int fileNumber, int iter,
							   int numRes, double time, int numComponents, int *glNodeNums);
	void writeNodeScalarToFile(DComplex *complexData, int numData, int glSub, int offset, int fileNumber, int iter,
							   int numRes, double time, int numComponents, int *glNodeNums);
	template<class Scalar, int dim>
	void writeNodeVectorToFile(SVec<Scalar, dim> &, int glSub, int offset, int fileNumber, int iter,
							   int numRes, double time, int numComponents, int startComponent, int *glNodeNums);
	void writeElemScalarToFile(double *data, int numData, int glSub, int offset, int fileNumber, int iter,
							   int numRes, double time, int totData, int *glElemNums);
	void writeElemScalarToFile(DComplex *complexData, int numData, int glSub, int offset, int fileNumber, int iter,
							   int numRes, double time, int totData, int *glElemNums);

private:
	int getHeaderNameBytes(int fileId) const;

	void getOutputFileName(char *result, int fileId, int clusterId, int iter);
	BinFileHandler* openBinaryOutputFile(int fileId, int clusterId, int iter, const char *flag);
	void outputHeader(BinFileHandler&, int, int);
	void writeArrayToBinFile(const double *data, int dataSize, int subId, int inDataOffset, int fileId,
							 int iterRank, int resultRank, double timeStamp, int inStateDataCount, int clusterItemCount);

public:
	// Shifting functions
	bool isShifted() { return isShift; }
	void initShift() {
	  if(!isShift) { isShift = true; shiftV = 0.0; }
	}
	void setShift(double w2) { isShift = true; shiftV = w2; }
	void setImpe(double f) { isShift = true; shiftV = 4.0*PI*PI*f*f; }
	void setOmega(double w) { isShift = true; shiftV = w*w; }
	void resetShift(double w2) { isShift = false; shiftV = w2; }
	double shiftVal() { return shiftV; }
	double freq() { return sqrt(shiftV)/(2.0*PI); }
	double omega() { return sqrt(shiftV); }
	double kappa() { if(numProps > 1) std::cerr << "Warning: assuming homogenous fluid (attr #1), k = " << sProps[0].kappaHelm << std::endl; return sProps[0].kappaHelm; }

	void setMRatio(double _mratio) { assert(_mratio >= 0.0 && _mratio <= 1.0); mratio = _mratio; }
	double getMRatio() const { return mratio; }

	// Housekeeping functions
	void cleanUp();
	void cleanAuxData();

	void makeDirectMPCs(int &numLMPC, ResizeArray<LMPCons *> &lmpc);
	int reduceMPCs(int numLMPC, ResizeArray<LMPCons *> &lmpc);
	bool checkLMPCs(int numLMPC, ResizeArray<LMPCons *> &lmpc);
	void transformLMPCs(int numLMPC, ResizeArray<LMPCons *> &lmpc);
	void addMpcElements(int numLMPC, ResizeArray<LMPCons *> &lmpc);
	void addMpcElementsIeq(int numLMPC, ResizeArray<LMPCons *> &lmpc);
	void addFsiElements(int numFSI, ResizeArray<LMPCons *> &fsi);
	Element* getElem(int topid) { return elemSet[topid]; }
	void UpdateContactSurfaceElements(DistrGeomState *, std::map<std::pair<int,int>,double> &);
	std::vector<int> contactSurfElems;
	void initializeParameters();
	void updateParameters();

	double global_average_E, global_average_nu, global_average_rhof;
	int num_arubber;
	void getARubberLambdaMu(double omega,
							complex<double> *lambda, complex<double> *mu);

	void makeEframe(int ele, int refnode, double *d);

	// Group stuff
	void setAttributeGroup(int a, int g);
	void setNodeGroup(int nn, int id);
	std::set<int> & getNodeGroup(int id) { return nodeGroup[id]; }
	std::map<int, std::set<int> > & getNodeGroups() { return nodeGroup; }
	void setSurfaceGroup(int sn, int id);

	// AN: required to access group number
	int getGroupNumber(int fileNum) { return oinfo[fileNum].groupNumber; }

	// AN: simple checks for id in node or attribute group
    bool isInNodeGroup(int);
	bool isInAttrGroup(int);

	// Sfem stuff
	enum Rprop { A, E, NU, RHO, T, KX, KY, KZ }; // sfem
	void setGroupRandomProperty(int g, Rprop prop_type, double mean, double std_dev);
	void printGroups();

	// POD-ROM sample nodes
	int sampleNodeCount() const { return sampleNode_.size(); }
	void sampleNodeAdd(int id) { sampleNode_.push_back(id); }

	typedef std::vector<int> SampleNodeList;
	SampleNodeList::const_iterator sampleNodeBegin() const { return sampleNode_.begin(); }
	SampleNodeList::const_iterator sampleNodeEnd()   const { return sampleNode_.end();   }

private:
	SampleNodeList sampleNode_;

public:
	// POD-ROM elementary lumping weights
	typedef std::map<int, double> ElementWeightMap;
	typedef std::vector<std::pair<int, int> > NodeDofPairVec;
	typedef std::vector<std::pair<int, int> > ElemDofPairVec;

	int elementLumpingWeightSize() const { return elementLumpingWeights_.size(); }
	int elementLumpingWeightLocalSize(int j=0) { return elementLumpingWeights_[j].size(); }
	ElementWeightMap::const_iterator elementLumpingWeightBegin(int j=0) const { return elementLumpingWeights_[j].begin(); }
	ElementWeightMap::const_iterator elementLumpingWeightEnd(int j=0)   const { return elementLumpingWeights_[j].end();   }

	NodeDofPairVec::const_iterator nodeDofSlotBegin() const { return nodeDofSlotPairVec_.begin(); }
	NodeDofPairVec::const_iterator nodeDofSlotEnd()   const { return nodeDofSlotPairVec_.end();   }

	std::vector<double>::const_iterator RedKVecBegin() const { return ReducedStiffVec.begin();}
	std::vector<double>::const_iterator RedKVecEnd() const { return ReducedStiffVec.end();}

	std::vector<double>::const_iterator UDEIMVecBegin() const { return UDEIMBasisVec.begin();}
	std::vector<double>::const_iterator UDEIMVecEnd() const { return UDEIMBasisVec.end();}

	std::vector<double>::const_iterator ROMLMPCVecBegin() const { return ROMLMPCVec.begin();}
	std::vector<double>::const_iterator ROMLMPCVecEnd() const { return ROMLMPCVec.end();}

	double * RedKData() { return ReducedStiffVec.data(); }
	double * UDEIMData() { return UDEIMBasisVec.data(); }

	ElemDofPairVec::const_iterator elemDofBegin() const { return elemDofPairVec_.begin(); }
	ElemDofPairVec::const_iterator elemDofEnd()   const { return elemDofPairVec_.end();   }

	void setElementLumpingWeight(int iele, double value);
	void pushBackStiffVec(double Kelem);
	void pushBackUDEIMVec(double Uelem);
	void pushBackROMLMPCVec(double value);
	void setSampleNodesAndSlots(int node, int dof);
	void setSampleElemsAndDOFs(int elem,int dof);
	void setLocalIndex(int j);
	int getLocalIndex() { return localIndex_; }

private:
	int localIndex_; // local bases
	std::vector<ElementWeightMap> elementLumpingWeights_;
	NodeDofPairVec nodeDofSlotPairVec_;
	ElemDofPairVec elemDofPairVec_;
	std::vector<double> ReducedStiffVec;
	std::vector<double> UDEIMBasisVec;
	std::vector<double> ROMLMPCVec;

protected:
	void closeOutputFileImpl(int fileIndex);
};


struct RandomProperty
{
	GeoSource::Rprop rprop;
	double mean;
	double std_dev;
	RandomProperty(GeoSource::Rprop rp, double m, double sd) { rprop = rp; mean = m; std_dev = sd; }
	RandomProperty(const RandomProperty &rp) { rprop = rp.rprop; mean = rp.mean; std_dev = rp.std_dev; }
};

struct Group
{
	std::vector<int> attributes;
	std::vector<RandomProperty> randomProperties;

	bool isAttributeInGroup(int a) {
		// AN: number of attributes in each group
		// are typically low; looping over the
		// vector.
		bool flag = false;
                for(int i=0; i<attributes.size(); i++) {
                        if(attributes[i] == a) {
                                flag = true;
                                break;
                        }
                }
                return flag;
        };

};

struct AttributeToElement
{
	std::vector<int> elems;
};

#endif
