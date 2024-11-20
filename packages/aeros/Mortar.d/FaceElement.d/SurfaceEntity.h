// ----------------------------------------------------------------
// HB - 06/30/03
// Modified by HB - 09/12/04 
//    * add NodeSet and related data 
//    * may be optimized for case when the face element set contains
//      only Tri3 or Quad4 face elements  (i.e node = vertex)                       
//      In such case, the vertex related data may be just pointed on
//      the node related data.
// ----------------------------------------------------------------
#ifndef _SURFACEENTITY_H_
#define _SURFACEENTITY_H_

// STL
#include <map>

// FEM headers
#include <Mortar.d/FaceElement.d/FaceElemSet.h>
#include <Mortar.d/FaceElement.d/FaceElement.h>
#include <Element.d/Element.h>

#ifdef SOWER_SURFS
#include <Utils.d/BinFileHandler.h>
#endif

#define HB_NODALNORMAL

class CoordSet;
class Connectivity;
class DistrGeomState;
template <class Scalar> class GenSubDomain;
typedef GenSubDomain<double> SubDomain;

class SurfaceEntity {
  private:	
	int Id;
	FaceElemSet ElemSet; 
        std::map<int,Node>* NodeCoordMap;

        int nNodes;
        int* gNodeIds;
        int* LlToGlNodeMap;
        std::map<int,int>* GlToLlNodeMap;

        CoordSet* NodeSet;
      
        int nVertices;
        int* gVertexIds;
        int* LlToGlVertexMap;
        std::map<int,int>* GlToLlVertexMap;
        int* LlVertexToLlNodeMap;

	bool LocalNumbering;

        int nTri3, nTri6, nQuad4, nQuad8; 
        Connectivity* ACMEBlocksMap;

#ifdef HB_NODALNORMAL
        double (*NdNormals)[3];
#endif
        bool ReverseNormals;
        bool IsShellFace;
        double ShellThickness;
        bool PreserveOrdering;

  public:
	// Constructors & destructor
	// ~~~~~~~~~~~~~~~~~~~~~~~~~
	SurfaceEntity();
	SurfaceEntity(int);
        SurfaceEntity(const SurfaceEntity&);

        ~SurfaceEntity();
	
        // Initialization & clean/clear methods
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        void Initialize();

        // Setup & update methods
        // ~~~~~~~~~~~~~~~~~~~~~~
	void MakeGlVertexIds();
	void MakeGlNodeIds();
        void MakeVertexMaps();
        void MakeNodeMaps();
        void MakeNodeSet(CoordSet&);
        void MakeNodeSet(std::map<int,Node>& NodeCoordMap);
        void ExtractNodeSet(CoordSet&);
        void MakeLlVertexToLlNodeMap();

        void SetUpData(CoordSet* cs=0);
        void SetUpVertexData();
        void SetUpNodeData();
        void UpdateNodeData(GeomState *);
        void UpdateNodeData(DistrGeomState *geomState, SubDomain **sd);
        void Renumber(std::map<int,int>& OldToNewNodeIds);
        void Renumber();

        void MakeACMEBlocksMap();
        int  FillACMEFaceBlocks(int* face_connectivity, std::map<int,int>& OldToNewNodeIds, bool vertexOnly=false);

#ifdef HB_NODALNORMAL
        double* ComputeNodalNormals(CoordSet& cs); // STILL EXPERIMENTAL ...
#endif
        void Reset(CoordSet* cs=0);

	// Set methods
	// ~~~~~~~~~~~
	void SetId(int);
	void AddFaceElement(int, int, int, int*);
	void AddFaceElement(FaceElement*);
        void RemoveFaceElement(int num);
        void SetNodeCoordMap(std::map<int,Node> *_map) { NodeCoordMap = _map; };
        void SetPtrNodeSet(CoordSet* ndSet);
        void SetReverseNormals(bool);
        void SetIsShellFace(bool);
        void SetShellThickness(double);
        void SetPreserveOrdering(bool);

	// Get/Accessor methods
	// ~~~~~~~~~~~~~~~~~~~~
	int GetId();
	int ID();
	int nFaceElements(); 	
	int GetnFaceElems();	
	int GetnNodes();	
	int GetnVertices();	
	FaceElemSet* GetPtrFaceElemSet();	
        CoordSet*    GetPtrNodeSet();
	FaceElemSet& GetFaceElemSet();	
        CoordSet&    GetNodeSet();
        bool IsRenumbered();
	int* GetPtrGlNodeIds();
	int* GetPtrGlVertexIds();
	int* GetPtrLlToGlNodeMap();
	int* GetPtrLlToGlVertexMap();

        std::map<int,int>* GetPtrGlToLlNodeMap();
        std::map<int,int>* GetPtrGlToLlVertexMap();

        int GetGlNodeId(int);
        Node& GetNode(int);

        int GetLlVertexInLlNode(int);
        int GetGlVertexId(int);
        Node& GetVertex(int);
 
        Connectivity* GetPtrACMEBlocksMap();
#ifdef HB_NODALNORMAL
        double* ViewNodalNormals() { return(reinterpret_cast<double*>(NdNormals)); }
#endif
        bool GetReverseNormals();
        bool GetIsShellFace();
        double GetShellThickness();

	// Print, display, ... methods
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~
	void PrintFaceElemSet();
	void Print();
#ifdef HB_NODALNORMAL
        void PrintNodalNormals();
#endif
        void PrintFaceNormal(CoordSet& cs);
#ifdef SOWER_SURFS
        void WriteSower(BinFileHandler& file);
#endif
};

#endif
