// --------------------------------------------------------------------------------
// HB - 07/05/03
// --------------------------------------------------------------------------------
// NOTES: (1) WE ASSUME THAT BY DEFINITION, A MORTAR CONDITION
//            INVOLVES ONLY A PAIR OF SURFACE ENTITIES 
//            for example, if Surf1 is to be tied (or contact)
//            with Surf2 & Surf3, then we should define 2 Mortar
//            conditions: Surf1 <-> Surf2 AND  Surf1 <-> Surf3
//            note that Surf1 SHOULD retain the SAME "type" (i.e.
//            being the master or the salve side) in the 2 Mortar 
//            conditions 
// --------------------------------------------------------------------------------
// SPECIFIC ISSUES: 
//        (1) Case of surface entity made of different type of
//            face elements ; for example Quad4 & Tri3 face elements 
//            -> CURRENTLY NOT SUPORTED
//        (2) Case of surface entities involved in several Mortar
//            conditions (see above NOTES (1))
//            -> CURRENTLY PARTIALY SUPORTED
//               SHOULD BE OK ON SHARED MEMORY MACHINE 
//            Here the pb is about data access: currently the surface
//            entities data are NOT duplicate for each Mortar condition 
//            i.e. the MortarHandler object (one per Mortar condition)
//            "points" to the pair of surface entities, but do NOT OWN
//            HIS OWN COPIES of the surface entities involved in the
//            Mortar condition. 
//            Making the copy of the surface entities shouldn't be a big
//            pb, but the question is then: do we always want to duplicate
//            data even if not needed (i.e. on a shared memory machine/node)
//            but then it can be difficult to decide which data to copy.  
//        (3) Case of a surface entity that contains non connex 
//            parts -> the mesh of the surface entity can be split
//            into several INDEPENDANT (no COMMON nodes) submeshes. 
//            for example, if SurfI(=Surf1) is to be tied (or contact)
//            with SurfII, but that SurfII is made of 2 indenpendants
//            surfaces (SurfII = Surf2+Surf3). In a sense, we violate
//            our first assumption (see above NOTES (1)). 
//            -> CURRENTLY NOT WELL SUPORTED
//               ACTUALY, THE CURRENT CODE WILL "WORK" IN THIS CASE,
//               BUT ACCORDING TO ME THIS CASE IS NOT WELL TREATED !!
//            Note that it would be better if SurfI is the slave side. 
//            The difficulties are BOTH theoritical and practical:
//              * how to define the edges of the mesh to do the Mortar
//                shape fcts modfication (if we choose to do them) 
//                Here is a tricky case:
//                  |                              | 
//                  |          Surf1=SurfI         | 
//                 a+-----+----+-----+------+------+d
//                             b     c 
//                              e f 
//                ----+-----+---+ +-----+-------+------+-
//                    Surf2     | |      Surf3
//                              + |
//                              | +
//                         Surf4| |Surf5
//                              + |
//                              | +
//                              | |
//
//                Mortar conditions A: SurfI <-> SurfII(=Surf2+Surf3)
//                                  B: Surf4 <-> Surf5
//
//                -SurfI being the slave side, for the current code
//                 (treating SurfII as ONE surface) the nodes a & b
//                 will be the ends nodes. BUT if the Mortar cond. A
//                 is split (manualy / automaticaly ??) into 2 Mortar
//                 cond. A1: SurfI <-> Surf2
//                       A2: SurfI <-> Surf3, then nodes a & b will
//                 be the ends nodes for A1, and nodes c & d for A2.
//                 So basicaly, this will eliminate the 2 Lagranges
//                 multipliers associated with the nodes b & c (if we
//                 choose to do the end nodes shape fcts modification) 
//                 Note that if we do NOT elimite any of the end nodes 
//                 Lagrange multipliers (either in the case of 2 Mortar
//                 cond. A & B, or 3 Mortar cond, A1, A2 & B), then the 
//                 global CCt matrix is singular. Note also, that if we
//                 only do the end nodes modification in Mortar cond. B
//                 (and NOT in Mortar cond. A or A1 & A2), we will 
//                 elinimate this sigularity.
//                -what happen if SurfII is the slave side. NOW depending
//                 on how we search for the end nodes, we will probably
//                 find nodes e & f as end nodes, which in a sense is good!!
//                -So basicaly the question is: should we (a) first LOOK for
//                 non connex parts in a given surface entity (i.e. find
//                 the SurfII is made of Surf2 & Surf3) and then REDIFINE
//                 the Mortar cond. accordingly (SPLIT A into A1 & A2) ; 
//                 or (b) do NOTHING (i.e. treat directly Mortar cond A)
//                 BUT then we will either need to always do the end nodes
//                 modification (can be expensive, and how to deal with 
//                 quadratic face elements ??) or check for the singularity of
//                 CCt (algebraicly ??) and eliminate it ?? (but what does this
//                 means "eliminate it" -> eliminate the associated LMPC ??) 
//                 What about the contact, when CCt can change (and so its
//                 possible singularity) throughout the iterations ??
//              * If we go for the "automatic" spliting of SurfII into
//                Surf2 & Surf3, how to handle it in a distributed machine,
//                or with the Sandia data interface (data given at the subdomain
//                level) 
// --------------------------------------------------------------------------------
#ifndef _MORTARHANDLER_H_
#define _MORTARHANDLER_H_

// Locally define flags
#include <Mortar.d/MortarDriver.d/MortarHandlerDefs.h>

#include <cstdio>
#include <vector>
#include <map>

#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
#include <Utils.d/resize_array.h>
#include <Utils.d/MyComplex.h>
#include <Parser.d/AuxDefs.h>

class CoordSet;
class LMPCons;
class Connectivity;
class Elemset;
class BCond;

class FaceElement;
class MortarElement;
class NodalMortarShapeFct;
template <class MasterFaceElementType, class SolveFaceElement, class MortarElementType> class FFIPolygonTemplate;
typedef FFIPolygonTemplate<FaceElement, FaceElement, MortarElement> FFIPolygon;

class ContactSearch;
class ContactTDEnforcement;

template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
class ConstrainedDSA;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
template <class Scalar> class GenSubDomain;
typedef GenSubDomain<double> SubDomain;
class FSCommunicator;
template <class Scalar> class FSCommPattern;
template <class Scalar> class GenSubDOp;
typedef GenSubDOp<double> SubDOp;


class MortarHandler {

        int Id;
	int MasterEntityId;
	int SlaveEntityId ;
        int TDEnfNumIter;
        int FrictionModel; // TD_FRICTIONLESS=1, TD_CONSTANT_FRICTION=2, TD_TIED=3, TD_SPOT_WELD=4, TD_PRESSURE_DEPENDENT=5, TD_VELOCITY_DEPENDENT=6, 
        int DIST_ACME; // 0: sequential, 1: parallel with centralized input on host (cpu with id 0), 2: parallel with distributed input by subdomain
                       // NOTE: currently only dist_acme == 0 is supported for Mortar method (statics and implicit dynamics)
        int MortarIntegrationRule;
        int CtcMode;
        ConstraintOptions *ConstraintOptionsData;

        double NormalSearchTol;
        double TangSearchTol;
        double TDEnfConvTol;
        double FrictionCoef[4];
        double MortarScaling;

        bool NoSecondary;
        bool AveragedNodalNormals;

        // Currently, the following member data ONLY POINT to 
        // the (pair of) surface entities involved in the 
        // Mortar condition. (So currently, the MortarHandler object
        // should NOT delete these pointed surface entities: they
        // could be involved in other Mortar condition. see after) 
        // For the (near!!) future, it could be interesting 
        // (& necessary !!) to COPY the data: the MortarHandler 
        // object will then OWNS a copy of the pair of surface entities.
        // NECESSARY if a SAME surface entity is involved in MORE
        // than only one Mortar condition. 
        // 2 main reasons:
        //  1) be able to renumber the nodes in a Mortar condition
        //     (i.e a PAIR of surface entities) INDEPENDANTLY of 
        //     the other Mortar conditions which can potentialy
        //     involve the SAME surface entity
        //     (see NOTES at the header)
        //  2) more parallel: enable CONCURENT access to the data 
        //     of the (originaly) SAME surface entity 
        //     (see distributed memory)  
        // -> introduce a Copy/Own flag ??     
        SurfaceEntity* PtrMasterEntity; 
        SurfaceEntity* PtrSlaveEntity; 
        SurfaceEntity* PtrGlobalMasterEntity;
        SurfaceEntity* PtrGlobalSlaveEntity;

        // ACME output data (see ACME API guide)
        int nFFI;
        int nACMEFFIData;
        int* Slave_face_block_id;       
        int* Slave_face_index_in_block;
        int* Master_face_block_id;
        int* Master_face_index_in_block;
        int* Slave_face_procs;
        int* Master_face_procs;
        int* ACMEFFI_index;
        double* ACMEFFI_data;
    
        // Mortar data (FFIPolygon, etc, ...)
        std::vector<FFIPolygon*> CtcPolygons;         // Face-Face-Interaction polygons built from the ACME data
	std::vector<MortarElement*>      MortarEls;   // "Mortar els" built on the "active" slace face els
	std::vector<NodalMortarShapeFct> NodalMortars;// "global" mortar shape fct associated to each "active" slave node

        std::vector<int>  ActiveSlaveNodes;           // slave nodes that participate to the LMPCs
        FaceElemSet       ActiveSlaveElemSet;         // slave face els that actually interact with the master surface
        FaceElemSet       ActiveMasterElemSet;        // master face els that actually interact with the slave surface
        Connectivity*     ActiveSlaveNodeToElem;      // maps an "active" slave node to the face els it is connected to
        Connectivity*     ActiveMasterNodeToElem;     // maps an "active" master node to the face els it is connected to
        Connectivity*     SlaveFaceToFFIConnect;      // maps an "active" slave face els to the FFI polygons it supports
        std::map<int,int> ActiveSlaveFacesToMortarEl; // maps an "active" slave face el to the "mortar el" it supports

        int nMortarLMPCs;
        int gIdFirstLMPC;         
        int gIdLastLMPC;         
        
	void ComputeOneFFIMandN(int iFFI, CoordSet& cs, std::vector<FFIPolygon*> &CtcPolygons);
	void MakeOneNodalMortarLMPC(int i, std::vector<FFIPolygon*> &CtcPolygons, bool Dual=false);

        bool SelfContact;

  public:
        // Public data 
        // ~~~~~~~~~~~
        enum Mortar_Type{STD=0, DUAL} MortarType;
        enum Interaction_Type{TIED=0, CTC, FSI} InteractionType;
        enum Geom_Type{NON_MATCHING=0, EQUIVALENCED} GeomType;  

        // Constructors & destructor
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        MortarHandler();
        MortarHandler(int _Id);
        MortarHandler(int _MasterEntityId, int _SlaveEntityId);
        MortarHandler(int _Id, int _MasterEntityId, int _SlaveEntityId);
        
        MortarHandler(int _MasterEntityId, int _SlaveEntityId, double _NormalTol);
        MortarHandler(int _MasterEntityId, int _SlaveEntityId, double _NormalTol, double _TangTol);
        
        ~MortarHandler();

        // Initialization & clean/clear methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        void Initialize();
        void DeleteFFIData();

        // Set methods
        // ~~~~~~~~~~~
        void SetId(int);
        void SetMasterEntityId(int _MasterEntityId);
        void SetSlaveEntityId(int _SlaveEntityId);
	void SetSurfaceEntityId(int _MasterEntityId, int _SlaveEntityId);

	void SetPtrMasterEntity( SurfaceEntity* _PtrMasterEntity);
	void SetPtrSlaveEntity(  SurfaceEntity* _PtrSlaveEntity );
	void SetPtrSurfaceEntity(SurfaceEntity* _PtrMasterEntity, SurfaceEntity* _PtrSlaveEntity);

        void SetNormalSearchTol(double _NormalTol);
        void SetTangSearchTol(double _TangTol);
        void SetNormalAndTangSearchTol(double _NormalTol, double _TangTol);
        void SetSearchTol(double _NormalTol, double _TangTol);

        void SetMortarType(int _MortarType);
        void SetMortarType(MortarHandler::Mortar_Type _MortarType);
        void SetInteractionType(int _InteractionType);
        void SetInteractionType(MortarHandler::Interaction_Type _InteractionType);
        void SetGeomType(int _GeomType);
        void SetGeomType(MortarHandler::Geom_Type _GeomType);

        void SetTDEnfParams(int _TDEnfNumIter, double _TDEnfConvTol);
        void SetFrictionCoef(double _FrictionCoef); // TD_CONSTANT_FRICTION
        void SetFrictionCoef(double _FrictionCoef, double _ReferencePressure, double _OffsetPressure, double _PressureExponent); // TD_PRESSURE_DEPENDENT
        void SetFrictionCoef(double _StaticCoef, double _DynamicCoef, double _VelocityDecay); // TD_VELOCITY_DEPENDENT

        void SetNoSecondary(bool _NoSecondary);
        void SetAveragedNodalNormals(bool _AveragedNodalNormals);
        void SetSelfContact(bool _SelfContact);
        void SetDistAcme(int _DistAcme);
        void SetMortarScaling(double _MortarScaling);
        void SetMortarIntegrationRule(int _MortarIntegrationRule);
        void SetCtcMode(int _CtcMode);

        void SetConstraintOptions(ConstraintOptions& _ConstraintOptionsData);

        // Get methods
        // ~~~~~~~~~~~
        int GetId();
        int ID();
        int GetMasterEntityId();
        int GetSlaveEntityId() ;
        int GetnFFI();

        double GetNormalTol();
        double GetTangentialTol();
        int GetCtcMode();

        MortarHandler::Mortar_Type GetMortarType();
        MortarHandler::Interaction_Type GetInteractionType();
        MortarHandler::Geom_Type GetGeomType();
 
	SurfaceEntity* GetPtrMasterEntity();
        SurfaceEntity* GetPtrSlaveEntity() ;

	FaceElemSet* GetPtrMasterFaceElemSet();
	FaceElemSet* GetPtrSlaveFaceElemSet();
	
	//FaceElemSet& GetRefMasterFaceElemSet();
	//FaceElemSet& GetRefSlaveFaceElemSet();
        
	int GetnMortarLMPCs();
        int GetIdFirstMortarLMPC();
        int GetIdLastMortarLMPC();

        ConstraintOptions* GetConstraintOptions();

        // Print, display, ... methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        void Print();
#ifdef HB_ACME_FFI_DEBUG
	void PrintFFIPolyVertices(FILE* file, int& firstVertId);	
	void PrintFFIPolyTopo(FILE* file, int& EdgeOffset, int& VertOffset, int elCode=101);
#endif

	// FFI/NFI search methods
        // ~~~~~~~~~~~~~~~~~~~~~~
        void PerformACMEFFISearch();
        void CreateACMEFFIData();

        // Mortar methods
        // ~~~~~~~~~~~~~~
        void CreateFFIPolygon();
    
        void AddMortarLMPCs(ResizeArray<LMPCons*>*, int& numLMPC, int& numCtc, int nDofs=0, int* Dofs=0);
        void AddWetFSI(ResizeArray<LMPCons*>* FSIArray, int& numFSI);

  private:
        // new persistent variables to support new methods below
        ContactSearch *search_obj;
        ContactTDEnforcement *contact_obj;
        int nMasterNodes, nSlaveNodes;
        int num_entity_keys;
        double *mass; // nodal mass
        int *dofmap;
        int *node_global_ids;
        double *data;
        int *share;
        double *density;
        double *wavespeed;
        int  num_comm_partners;
        int* comm_proc_id;
        int* number_nodes_to_partner;
        int* comm_node;

  public:
        // New methods to allow more direct interfacing with ACME
        void build_search(bool tdenforceFlag = false, int numSub = 0, SubDomain **sd = 0);
        void set_search_data(int interaction_type); // interaction_type: 1 = NodeFace, 2 = NodeSurface, 3 = NodeNode, 3 = FaceFace, 4 = FaceCoverage, 5 = ElementElement
        void set_node_configuration(int config_type); // config_type: 1 = current, 2 = predicted
        void set_node_configuration(int config_type, int numSub, SubDomain **sd);
        void set_node_constraints(int numDBC, BCond *dbc);
        void set_node_constraints(int numSub, SubDomain **sd);
        void set_search_options();
        void perform_search(int search_algorithm, double dt_old = 0.0, double dt = 0.0); // search_algorithm: 1 = static-1, 2 = static-2, 3 = dynamic-2, 4 = augmented dynamic-2
        void get_interactions(int interaction_type); // interaction_type: 1 = NodeFace, 2 = NodeSurface, 3 = NodeNode, 4 = FaceFace, 5 = FaceCoverage, 6 = ElementElement
        void build_td_enforcement();
        void make_nodal_mass(SparseMatrix *M, ConstrainedDSA *c_dsa, SparseMatrix *Mcc);
        void make_nodal_mass(SubDOp *M, SubDomain **sd);
        void make_kinematic_partitioning(Elemset &packedEset, Connectivity *nodeToElem);
        void make_kinematic_partitioning(int numSub, SubDomain **sd);
        void compute_td_contact_force(double dt_old, double dt, Vector &f);
        void compute_td_contact_force(double dt_old, double dt, DistrVector &force); // multidomain version
        void get_plot_variable(int plot_var, double *all_data); // plot_var: 1 = CONFACE, 2 = NORMAL_FORCE_MAG, 3 = NORMAL_TRACTION_MAG, 4 = TANGENTIAL_FORCE_MAG, 5 = TANGENTIAL_TRACTION_MAG,
                                                                // 6 = CDIRNORX, 7 = CDIRNORY, 8 = CDIRNORZ, 9 = CDIRTANX, 10 = CDIRTANY, 11 = CDIRTANZ, 12 = SLIP_MAG, 13 = NODAL_DISSIPATION,
                                                                // 14 = CONTACT_AREA, 15 = GAP_CUR, 16 = GAP_OLD
        void get_plot_variable(int plot_var, double **plot_data, int numSub, gsl::span<SubDomain *> sd);
        void get_plot_variable(int plot_var, double **plot_data, int numSub, gsl::span<GenSubDomain<complex<double> > *>sd) { }
        void get_global_variable(int var, double &value); // var: 1 = FORCE_X, 2 = FORCE_Y, 3 = FORCE_Z, 4 = FORCE_NORM, 5 = DISSIPATION, 6 = CONSTRAINT_NORM, 7 = INC_FORCE_NORM
        int get_num_nodes() { return nMasterNodes + nSlaveNodes; }
        int * get_node_global_ids() { return node_global_ids; }
        void remove_gap(Vector &d);
        void make_share(int numSub, SubDomain **sd);
};

#endif
