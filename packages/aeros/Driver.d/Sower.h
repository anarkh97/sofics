#ifndef _SOWER_H_
#define _SOWER_H_

/**
   FEM's SOWER
   DISTRIBUTED INPUT
   June 2004
**/

#include <string>
#include <sstream>
#define FILENAME_LENGTH 80
extern std::string clusterData_;
extern std::string subdomains_;
extern std::string decomposition_;
extern std::string connectivity_;
//#define SUBTOSUBINFILE
//#define SOWER_DEBUG 

#include <Utils.d/BinFileHandler.h>
#include <Utils.d/Connectivity.h>
#include <Element.d/Element.h>
#include <Element.d/ElemAccess.h>
#include <Element.d/Sommerfeld.d/LineSommerBC.h>
#include <Element.d/Sommerfeld.d/Line2SommerBC.h>
#include <Element.d/Sommerfeld.d/TriangleSommerBC.h>
#include <Element.d/Sommerfeld.d/QuadSommerBC.h>
#include <Element.d/Sommerfeld.d/Triangle6SommerBC.h>
#include <Element.d/Sommerfeld.d/IsoParamQuadSommer.h>
#include <Element.d/Sommerfeld.d/IsoParamTriSommer.h>
#include <iostream>
#include <list>
#include <map>
#include <vector>
#include <Driver.d/Access.h>
#include <Utils.d/CompositeInfo.h>
#include <Utils.d/MFTT.h>

//#define SOWER_DEBUG_SURFS
#ifdef SOWER_SURFS
#include <Mortar.d/FaceElement.d/SurfaceEntity.h>
template <class Type> class ResizeArray;
#endif

/* TypeTag defines tags for each category of data NODES, ELEMENTS ... */
enum TypeTag{END_TYPE=0, NODES_TYPE=1, ELEMENTS_TYPE=2, ATTRIBUTES_TYPE=3, FORCES_TYPE=4,
             MATERIALS_TYPE=5, DISPLACEMENTS_TYPE=6, EFRAMES_TYPE=7, COMPOSITEL_TYPE=8,
             CFRAMES_TYPE=9, CLUSTER_TYPE=10, CONV_TYPE=11, IDISP_TYPE=12, IDISP6_TYPE=13,
             IVEL_TYPE=14, ITEMP_TYPE=15, HDIR_TYPE=16, HNEU_TYPE=17,
             DNB_TYPE=18, SCAT_TYPE=19, ARB_TYPE=20, SENSOR_TYPE=21, ACTUATOR_TYPE=22,
             USDD_TYPE=23, USDF_TYPE=24, SURFACES=25, PRESSURE_TYPE=26, PRELOAD_TYPE=27, BOFFSET_TYPE=28,
             EFRAME_TYPE=29, DIMASS_TYPE=32, COMPOSITEC_TYPE=33, TETT_TYPE=34, YMTT_TYPE=35, LMPC_TYPE=36, HLMPC_TYPE=37,
             RAD_TYPE=38, WET_TYPE=39
};

typedef long long INT64BIT;

class ObjectOrdering 
{
  bool* subIsInClus;
  Connectivity* objToSub; /* has to be sorted */
 public:
  ObjectOrdering(Connectivity* c, bool* b) : subIsInClus(b), objToSub(c) { };
  bool operator()(int n1, int n2);
};

class Sower;

/*************** DATA STRUCTURES ************/

class DataStruct
{
	Connectivity subToData; /* cluster to data connectivity */
	Connectivity clusterToData;
	//TypeTag parentType;
public:
	DataStruct(Connectivity subToData) : subToData(std::move(subToData))
	{ }
	const Connectivity& getClusterToData(Connectivity * clusToSub) {
		if(clusterToData.csize() == 0)
			clusterToData = clusToSub->transcon(subToData); // PJSA
		return clusterToData;
		//return(clusToSub->transcon(subToData));
	}
	const Connectivity &getSubToData(void) { return(subToData); }

	virtual size_t writeObjectData(int index, BinFileHandler& file, int curObjID) { return 0; };

	size_t write(int clusNumber, Connectivity* clusToSub, int numSubdomains,
	             BinFileHandler& file, INT64BIT& curRangeSetLocation);

	virtual bool isVariableSize() { return false; };
};

template<typename DataType, typename IOObject>
class GenDataStruct : public DataStruct
{
 private:
  DataType data;
 public:
  GenDataStruct(DataType d, Connectivity ctd) : DataStruct(ctd), data(d)
    { }
  size_t writeObjectData(int index, BinFileHandler& file,  int curObjID)
    { return(IOObject::write(data,index,file,curObjID)); }
};


/* represents a range od data in the binary file @see RangeSet*/
class RangeEntry
{
 public:
  INT64BIT offset; /* offset in file where the entry starts */
  int beginNum;    /* number of first element in range */
  int endNum;      /* number of last element in range */
  RangeEntry(INT64BIT off, int n) : offset(off), beginNum(n) { };
  RangeEntry(BinFileHandler& file)
    {
      file.read(&offset, 1);
      file.read(&beginNum, 1);
      file.read(&endNum, 1);
    }
  size_t write(BinFileHandler& file)
    {
      file.write(&offset, 1);
      file.write(&beginNum, 1);
      file.write(&endNum, 1);
      return 0;
    }
};

/*
  RangeSet will hold the location of the data of a given subdomain for a given datatype
  using RangeEntries
  creation : call start and end when you encounter begin or end of a range of data for the subdomain
*/
class RangeSet
{
 private:
  std::list<RangeEntry*> rs;
  RangeEntry* current;
 public:
  RangeSet() { current = 0; };
  RangeSet(BinFileHandler& file) 
    {
      int size;
      file.read(&size, 1);
      for(int i=0; i<size; ++i)
        rs.push_back(new RangeEntry(file));
    }
  void print()
    {
      std::cerr << "Size : " << rs.size() << std::endl;
      for(std::list<RangeEntry*>::iterator it = rs.begin(); it!=rs.end(); ++it)
        std::cerr << "range from element " << (*it)->beginNum << " to " << (*it)->endNum << std::endl;
    }
  bool empty() { return(rs.empty()); } // true if range is empty
  void start(INT64BIT off, int n)
    {
      if(current == 0) // then start has been called twice in a row, which shall not happen !
	current = new RangeEntry(off,n);
      else
	std::cerr << " ** Sower.h Error start called twice in a row !" << std::endl;
    }
  void end(int n)
    {
      if(current) {
        current->endNum = n;
        rs.push_back(current);
        current=0;
      }
      else
	std::cerr << " ** Sower.h Error end called twice in a row or before start !" << std::endl;
    }
  size_t write(BinFileHandler& file)
    {
      int size = rs.size();
      file.write(&size,1);
      for(std::list<RangeEntry*>::iterator it = rs.begin();it != rs.end(); ++it)
	(*it)->write(file);
      return 0;
    }
  int size() { return(rs.size()); } /* number of entries in rangeset*/
  std::list<RangeEntry*>* getRangeEntries() { return(&rs); };
};

class TOCEntry
{
 public:
  TypeTag dataType;
  INT64BIT data;
  INT64BIT range;
  TOCEntry(BinFileHandler& file)
    {
      file.read(&dataType, 1);
      file.read(&data, 1);
      file.read(&range, 1);
#ifdef SOWER_DEBUG
      std::cerr << "Entry : type " << dataType << " data@ " << data << " range@ " << range << std::endl;   
#endif
    }
};

/*************** END OF DATA STRUCTURES ************/

/*
 * Sower class
 * stores all data parsed by fem that will be distributed to the clusters, its member functions create DataStruct objects for 
 * each categories (NODES, ELEMENTS, ATTRIBUTES ...) (one per category) with associated cluster to object connectivities that will help determinate 
 * which data goes to what cluster.
 */
class Sower 
{
private:
    std::map<TypeTag,DataStruct*> entries;
    int nCluster = 0;                                /* true number of cluster */
    int numSubdomains = 0;                           /* number of subdomains in mesh decomposition */
    std::unique_ptr<Connectivity> clusToSub; // used for writing
    std::unique_ptr<Connectivity> subToClus; // used for reading
    std::map<TypeTag,std::map<int,int>* > globalToLocal;   // mapsglobalnumber to local indexes
    std::map<TypeTag,std::map<int,RangeSet*>* > rangeSets; // rangesets for each type of data
    std::map<TypeTag,TOCEntry*> toc;                  // we read the table of content only one time
    bool tocRead = true;

    int nSurfaces = 0; //HB
    ResizeArray<SurfaceEntity*>* Surfaces; //HB

    void nomask( Connectivity* mn, int* maskno, int nelem ); /**< utility function for cpuGreedyBuild **/

public:
  /*
    This constructor will create the cluster to elements connectivity from which we will
    compute all the others
   */
  Sower(Connectivity* subToElem, Elemset& eset, int nClus, ResizeArray<SurfaceEntity*>* Surfs, Connectivity *cpuToSub); //HB
  Sower() = default; // PJSA: constructor for reading
  ~Sower();

  /*** ACCESS ***/
  const Connectivity* getSubToClus() const { return subToClus.get(); }
  
  /*** WRITING ***/
  void writeSubFile(BinFileHandler& f)
    {
      Connectivity * _subToClus = clusToSub->alloc_reverse();
      _subToClus->write(f);
      delete _subToClus;
    }

  void readSubFile(BinFileHandler& f)
    {
      subToClus = std::make_unique<Connectivity>(f);
    }

	std::unique_ptr<BinFileHandler> openBinaryFile(int sub)
	{
		if(!subToClus) {
			BinFileHandler fp(subdomains_.c_str(), "r");
			readSubFile(fp);
		}
		std::stringstream fName;
		fName << clusterData_;
		if ( (*subToClus)[sub].size() > 0 )
			fName << ( (*subToClus)[sub][0]+1 );
		// opening cluster file
		return std::make_unique<BinFileHandler>(fName.str().c_str(), "r");
	}

    int clusterIndex(int glSub) {
        if(!subToClus) {
            BinFileHandler fp(subdomains_.c_str(), "r");
            readSubFile(fp);
        }
        return (*subToClus)[glSub][0];
    }

  /*
   * addChildToParentData
   * adds data 'data' to the Sower object given a child to parent connectivity from which the cluster 
   * to child connectivity will be computed
   * @param thisType child type defined by enum @see TypeTag
   * @param parentType parent type defined by enum @see TypeTag
   * @param ndata unused
   * @param data considered data
   * @param childToParent connectivity
   * template parameters :
   * @param Datatype will ususally be a pointer to the object data we want to write : for example Elemset* for NODES
   * @param MapType will be a parent to child connectivity or implicit connectivity
   */
  template<typename IOObject, typename DataType, typename MapType>
    void addChildToParentData(TypeTag thisType, TypeTag parentType, int ndata, DataType data, MapType childToParent);

  /*
   * addParentToChildData
   * adds data 'data' to the Sower object given a parent to child connectivity from which the cluster 
   * to child connectivity will be computed
   * NOTE : the only difference with the previous function childToParent is that we will have to reverse the connectivity
   * @param thisType child type defined by enum @see TypeTag
   * @param parentType parent type defined by enum @see TypeTag
   * @param ndata unused
   * @param data considered data
   * @param childToParent connectivity
   * template parameters :
   * @param Datatype will ususally be a pointer to the object data we want to write : for example Elemset* for NODES
   * @param MapType will be a parent to child connectivity or implicit connectivity
   */
  template<typename IOObject, typename DataType, typename MapType>
    void addParentToChildData(TypeTag thisType, TypeTag parentType, int ndata, DataType data, MapType parentToChild); 

  /*
   * aka cpuGreedy build from the fluid's sower
   * greedy will compute the cluster to subdomain connectivity
   * @return cluster to subdomain connectivity
   * @param etoe subdomain to subdomain connectivity
   * @param ntoe node to subdomain connectivity
   * @param eton subdomain to node conncectivity
   * @param sizes weight of each subdomain\
   *        for now the array contains the number of elements in each subdomain
   */
  Connectivity* greedy(Connectivity *etoe,
                       Connectivity *ntoe, Connectivity *eton, long long *sizes, 
                       int decompNumSubdomains, int expectedNumSubdomains);

  /* writes the distributed data in as many files as clusters */
  void write();
  
  /*** READING ***/

  /* reads the data associated with the type of IOObject and returns it in an appropriate 
   * fem element for use with geoSource & friends */
  template<typename IOObject>
    typename IOObject::oType read(BinFileHandler& file, int subNum, int*& localToGlobal, bool rereadToc=false);

  /* returns the rangeSet for a particlar datatype and subdomain -- bufferised so we read rangesets once */
  RangeSet * getRangeSet(TypeTag datatype, int subdomain, BinFileHandler &file)
    {
      std::map<TypeTag, std::map<int,RangeSet*>* >::iterator it = rangeSets.find(datatype);
      std::map<int,RangeSet*>* rs;
      if(it == rangeSets.end()) { // range sets never read for this type of data
        // find location of rangeset in TOC
        std::map<TypeTag,TOCEntry*>::iterator ite = toc.find(datatype);
        if(ite == toc.end()) {
#ifdef SOWER_DEBUG
          std::cerr << "Warning : no entry in toc for this datatype" << datatype << std::endl;
#endif
          return 0;
        }
        file.seek(toc[datatype]->range);
	// now read the rangeset
        rs = new std::map<int, RangeSet*>();
        int numOfSubs;
        file.read(&numOfSubs, 1);
        for(int i=0; i<numOfSubs; ++i) {
          int subNumber;
          file.read(&subNumber, 1);
          (*rs)[subNumber] = new RangeSet(file);
        }
        rangeSets[datatype] = rs;
      }
      else {
        rs = (*it).second;
      }
#ifdef SOWER_DEBUG
      // show me the rangeset
      //for(std::map<int, RangeSet*>::iterator ite = rs->begin(); ite!=rs->end(); ++ite) {
      //  std::cerr<<"Sub "<<(*ite).first<<std::endl;
      //  (*ite).second->print();
      //}
#endif      
      std::map<int, RangeSet*>::iterator ite = rs->find(subdomain); 
      if(ite == rs->end()) {
#ifdef SOWER_DEBUG
        std::cerr << " ** WARNING getRangeSet :: subdomain " << subdomain << " has no data of type " << datatype << std::endl;
#endif
        return(0);
      }

      return((*ite).second);
    }

  /* reads table of content of the file */
  void readToc(BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
      std::cerr << " ** TABLE OF CONTENT : " << std::endl; 
#endif
      file.seek(0); // change if something comes before !!
      while(true) {
        TOCEntry* te = new TOCEntry(file);
	if(te->dataType == END_TYPE)
	  break;
	toc[te->dataType]=te;
      }
      if(toc.size()==0)
	std::cerr << " ** Warning : readToc : Table of content is empty !" << std::endl;
    }
  
  /* reads nobjs Scalar objects from BinFileHandler */
  template <class Scalar>
    size_t read(Scalar *p, int nobjs, BinFileHandler& file)
    {
      file.read(p,nobjs);
      return 0;
    }
  
  /* reads nobjs Scalar objects from BinFileHandler and transforms the global numbers to local ones using
   * the appropriate translation table according to the type identifier */
  template <class Scalar>
    size_t readNum(TypeTag type, Scalar *p, int nobjs, BinFileHandler& file)
    {
      file.read(p, nobjs);
      std::map<TypeTag, std::map<int, int>* >::iterator it = globalToLocal.find(type);
      if(it == globalToLocal.end()) {
        std::cerr << " ** Error : globalToLocal table not found for "<< type << std::endl;
        exit(1);
      }
      std::map<int,int>* table = (*it).second;
      for(int i = 0; i < nobjs; ++i) {
        std::map<int, int>::iterator fs = (*table).find(p[i]); 
	if(fs == table->end())
	  std::cerr << " ** WARNING (CRITICAL) : could'nt find translation for entry " << p[i] 
                    << " in table "<< type << std::endl;
	else
	  p[i] = (*table)[ p[i] ];
      }
      return 0;
    }

    /* prints the connectivities contained in the sower object */
    void printDebug() const;
    void printTables() const
    {
        for( const auto &ite : globalToLocal ) {
            std::cerr << "Table for type " << ite.first << std::endl;
            for( const auto &it : *ite.second ) {
                std::cerr << it.first << " -> " << it.second << std::endl;
            }
        }
    }

  //HB
  void writeSurfaces(BinFileHandler& file, INT64BIT& tocCurrentOffset, 
                     INT64BIT& nextEntryOffset, INT64BIT& curRangeSetLocation);
};

/**
   BEGINNING of IO objects for datatypes
**/

class ElemsetIO
{
  public:
  static size_t write(Elemset* eset, int index, BinFileHandler& file, int curObjID);
  typedef Elemset* oType;
  const static TypeTag dataType = ELEMENTS_TYPE;

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
std::cerr << "Sower.h, readData, ElemsetIO" << std::endl;
#endif
      int elType, nNodes, nodes[64];
      // Element * el;
      s->read(&elType, 1, file);
      s->read(&nNodes, 1, file);
      if(nNodes == 0) { std::cerr << "Sower.h found void element in readData, elType = " << elType << std::endl; return; }
      s->readNum<>(NODES_TYPE, nodes, nNodes, file);
      obj->elemadd(localIndex, elType, nNodes, nodes);
      double pressure;
      s->read(&pressure, 1, file);
      int preload_size;
      s->read(&preload_size, 1, file);
      std::vector<double> preload;
      preload.reserve(preload_size);
      s->read(&(preload[0]), preload_size, file);
      PressureBCond *pbc = new PressureBCond;
      pbc->setData(-1, pressure, 0, false); // XXX add support for cases
      (*obj)[localIndex]->setPressure(pbc);
      (*obj)[localIndex]->setPreLoad(preload);
    }

  static oType create(int size)
    {
      return(new Elemset(size));
    }
};

class AttribIO
{
 public:
  static size_t write(std::pair<int, std::map<int,Attrib>* >* attribp, int index, BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing attribute index = " << index << " curObjID = " << curObjID << std::endl;
#endif
      file.write(&curObjID, 1);
      std::map<int, Attrib>* attrib = attribp->second;
      file.write(&(*attrib)[index].nele, 1);
      file.write(&(*attrib)[index].attr, 1);
      return 0;
    }

  typedef std::pair<int, std::map<int,Attrib>* >* oType;
  const static TypeTag dataType = ATTRIBUTES_TYPE;

  static oType create(int size)
    {
      Attrib a;
      return(new std::pair<int, std::map<int,Attrib>* >(0, new std::map<int,Attrib>() ));
    }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
std::cerr << "Sower.h, readData, AttribIO" << std::endl;
#endif
      Attrib a;
      obj->first++;
      s->readNum(ELEMENTS_TYPE, &(a.nele), 1, file); 
      s->readNum(MATERIALS_TYPE, &(a.attr), 1, file);
      (*(obj->second))[localIndex]=a; 
    }
};

class NodesIO
{
 public:
  /* will write to file the node at index index's coordinates */
  static size_t write(CoordSet* nodes, int index, BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing node index = " << index << " curObjID = " << curObjID << std::endl;
#endif
      Node* n = (*nodes)[index];
      double xyz[3];
      xyz[0] = n->x;
      xyz[1] = n->y;
      xyz[2] = n->z;
      file.write(&curObjID, 1);	     
      file.write(xyz,3);
      return 0;
    }

  typedef CoordSet* oType;
  const static TypeTag dataType = NODES_TYPE;

  static oType create(int size)
    {
      return(new CoordSet(size));
    }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
std::cerr << "Sower.h, readData, NodesIO" << std::endl;
#endif
      double xyz[3];
      s->read(xyz, 3, file);
      obj->nodeadd(localIndex, xyz);
    }
};

/* PJSA: not used
class PressureIO
{
 public:
  static size_t write(Elemset* el, int index, BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing pressure index = " << index << " curObjID = " << curObjID << std::endl;
#endif
      double press = (*el)[index]->getPressure();
      file.write(&curObjID, 1);	     
      file.write(&press,1);
      return 0;
    }

  typedef Elemset* oType;
  const static TypeTag dataType = PRESSURE_TYPE;

  static oType create(int size) { return(new Elemset(size)); }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
      std::cerr << "Sower.h, readData, PressureIO" << std::endl;
#endif
      double press;
      s->read(&press, 1, file);
      (*obj)[localIndex]->setPressure(press);
    }
};

class PreLoadIO
{
 public:
  static size_t write(Elemset* el, int index, BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing preload index = " << index << " curObjID = " << curObjID << std::endl;
#endif
      double preload = (*el)[index]->getPreLoad();
      file.write(&curObjID, 1);	    
      file.write(&preload,1);
      return 0;
    }

  typedef Elemset* oType;
  const static TypeTag dataType = PRELOAD_TYPE;

  static oType create(int size) { return(new Elemset(size)); }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
      std::cerr << "Sower.h, readData, PreLoadIO" << std::endl;
#endif
      double preload;
      s->read(&preload, 1, file);
      (*obj)[localIndex]->setPreLoad(preload); 
    }
};
*/

class DMassIO
{
 public:
  /* will write to file the node at index index's coordinates */
  static size_t write(std::vector<DMassData*>* dmd, int index, BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing discrete mass index = " << index << " curObjID = " << curObjID << std::endl;
#endif
      file.write(&curObjID, 1);	    // node #
      file.write(& (*dmd)[index]->dof , 1);
      file.write(& (*dmd)[index]->diMass , 1);
      return 0;
    }

  typedef std::vector<DMassData*>* oType;
  const static TypeTag dataType = DIMASS_TYPE;

  static oType create(int size)
    {
      return new std::vector<DMassData*> ();
    }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
      std::cerr << "Sower.h, readData, DMassIO" << std::endl;
#endif
      DMassData* dmd = new DMassData();
      dmd->node = localIndex;
      s->read(& (dmd->dof) , 1 , file);
      s->read(& (dmd->diMass), 1, file);
      obj->push_back(dmd);
    }
};

class BoffsetIO
{
 public:
  /* will write to file the node at index index's coordinates */
  static size_t write(std::vector<OffsetData>* data, int index, BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing beam offset index = " << index << " curObjID = " << curObjID << std::endl;
#endif
      file.write(&curObjID, 1);
      int fl[2];
      fl[0] = (*data)[index].first; // must be global
      fl[1] = (*data)[index].last;  // must be global
      file.write(fl,2);
      file.write( ( (*data)[index].o ), 3);
      return 0;
    }

  typedef std::vector<OffsetData>* oType;
  const static TypeTag dataType = BOFFSET_TYPE;

  static oType create(int size)
    {
      return(new std::vector<OffsetData>());
    }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
      std::cerr << "Sower.h, readData, BoffsetIO" << std::endl;
#endif
      OffsetData od;
      s->readNum(ELEMENTS_TYPE, &(od.first), 1, file); 
      s->readNum(ELEMENTS_TYPE, &(od.last), 1, file); 
      double ofd[3];
      s->read(ofd, 3, file);
      od.o[0] = ofd[0];
      od.o[1] = ofd[1];
      od.o[2] = ofd[2];

      obj->push_back(od);

    }
};

class EFrameIO
{
   public:
  /* will write to file the node at index index's coordinates */
  static size_t write(std::pair<int,ResizeArray<EFrameData>* >* data, int index, BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing eframe index = " << index << " curObjID = " << curObjID << std::endl;
#endif
      file.write(&curObjID, 1);
      //int fl;                                // curObjID !
      //fl = (*(*data).second)[index].elnum;
      //file.write(&fl,1);
      double efr[9];
      efr[0] = (*(*data).second)[index].frame[0][0];
      efr[1] = (*(*data).second)[index].frame[0][1];
      efr[2] = (*(*data).second)[index].frame[0][2];
      efr[3] = (*(*data).second)[index].frame[1][0];
      efr[4] = (*(*data).second)[index].frame[1][1];
      efr[5] = (*(*data).second)[index].frame[1][2];
      efr[6] = (*(*data).second)[index].frame[2][0];
      efr[7] = (*(*data).second)[index].frame[2][1];
      efr[8] = (*(*data).second)[index].frame[2][2];
      file.write( efr, 9);
      return 0;
    }

  typedef std::vector<EFrameData>* oType;
  const static TypeTag dataType = EFRAME_TYPE;

  static oType create(int size)
    {
      return(new std::vector<EFrameData>());
    }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
      std::cerr << "Sower.h, readData, EFrameIO" << std::endl;
#endif
      EFrameData efd;
      //int fl;
      //s->read(&fl, 1, file);
      double efr[9];
      s->read(efr, 9, file);
      efd.elnum = localIndex;
      
      efd.frame[0][0] = efr[0];
      efd.frame[0][1] = efr[1];
      efd.frame[0][2] = efr[2];
      efd.frame[1][0] = efr[3];
      efd.frame[1][1] = efr[4];
      efd.frame[1][2] = efr[5];
      efd.frame[2][0] = efr[6];
      efd.frame[2][1] = efr[7];
      efd.frame[2][2] = efr[8];

      obj->push_back(efd);

    }
};

template<TypeTag DATATYPE>
class BCDataIO
{
 public:
  static size_t write(std::pair<int, BCond *>* data, int index, BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing boundary condition type " << DATATYPE << " index = " << index << " curObjID = " << curObjID << std::endl;
#endif
      if(index < data->first) {
	BCond bc = (data->second)[index];
	file.write(&curObjID, 1);
	file.write(&bc.nnum,1);
	file.write(&bc.dofnum,1);
	file.write(&bc.val,1);
      }
      else {
        std::cerr << "out of range" << std::endl;
      }
      return 0;
    }

  typedef std::list<BCond *>* oType;
  const static TypeTag dataType = DATATYPE;
  static oType create(int size)
    {
      return(new std::list<BCond *>());
    }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
std::cerr << "Sower.h, readData, BCDataIO" << std::endl;
#endif
      BCond* bc = new BCond();
      s->readNum(NODES_TYPE, &(bc->nnum), 1, file); // read node number on which BC applies and translate it
      s->read(&(bc->dofnum), 1, file);
      s->read(&(bc->val), 1, file);
      obj->push_back(bc);
    }
};

template<TypeTag DATATYPE>
class ComplexBCDataIO
{
 public:
  static size_t write(std::pair<int, ComplexBCond *>* data, int index, BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing complex dirichlet bc index = " << index << " curObjID = " << curObjID << std::endl;
#endif
      if(index < data->first) {
        ComplexBCond cbc = (data->second)[index];
	file.write(&curObjID, 1);
        file.write(&cbc.nnum,1);
        file.write(&cbc.dofnum,1);
        file.write(&cbc.reval,1);
        file.write(&cbc.imval,1);
      }
      else {
        std::cerr << "out of range" << std::endl;
      }
      return 0;
    }

  typedef std::list<ComplexBCond *>* oType;
  const static TypeTag dataType = DATATYPE;
  static oType create(int size)
    {
      return(new std::list<ComplexBCond *>());
    }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
      std::cerr << "Sower.h, readData, ComplexBCDataIO" << std::endl;
#endif
      ComplexBCond* cbc = new ComplexBCond();
      s->readNum(NODES_TYPE, &(cbc->nnum), 1, file); // read node number on which BC applies and translate it
      s->read(&(cbc->dofnum), 1, file);
      s->read(&(cbc->reval), 1, file);
      s->read(&(cbc->imval), 1, file);
      obj->push_back(cbc);
    }
};

template<TypeTag DATATYPE>
class SommerDataIO
{
 public:
  static size_t write(std::pair<int, SommerElement **>* data, int index, BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing sommerfeld element index = " << index << " vs curObjID = " << curObjID << std::endl;
#endif
      if(index < data->first) {
	SommerElement* so = (data->second)[index];
        int elType = so->getElementType();//not sure it works....
        int sFlag = so->sFlag;
        complex<double> ss = so->soundSpeed;
        int numN = so->numNodes();
        auto nodes = so->getNodes();
	file.write(&curObjID, 1);
        file.write(&elType, 1);
        file.write(&sFlag, 1);
        file.write(&ss, 1);
        file.write(&numN, 1);
        file.write(nodes, numN);
      }
      else {
        std::cerr << "out of range" << std::endl;
      }
      return 0;
    }

  typedef std::list<SommerElement *>* oType;
  const static TypeTag dataType = DATATYPE;
  static oType create(int size)
    {
      return(new std::list<SommerElement *>());
    }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
std::cerr << "Sower.h, readData, SommerDataIO" << std::endl;
#endif
      int etype;
      int nNodes;
      int sflag;
      complex<double> ss;
      s->read(&etype,1,file);
      s->read(&sflag,1,file);
      s->read(&ss,1,file);
      s->read(&nNodes,1,file);
      int* n = new int[nNodes];
      s->readNum(NODES_TYPE, n, nNodes, file);
      SommerElement *ele;

      switch(etype)
      {
        case 1:
          ele = new LineSommerBC(n[0], n[1]);
          ele->sFlag = sflag;
          ele->soundSpeed = ss;
          break;
        case 2:
          ele = new Line2SommerBC(n[0], n[1], n[2]);
          ele->sFlag = sflag;
          ele->soundSpeed = ss;
          break;
        case 3:
          ele = new TriangleSommerBC(n[0], n[1], n[2]);
          ele->sFlag = sflag;
          ele->soundSpeed = ss;
          break;
        case 4:
          ele = new QuadSommerBC(n[0], n[1], n[2], n[3]);
          ele->sFlag = sflag;
          ele->soundSpeed = ss;
          break;
        case 6:
          ele = new Triangle6SommerBC(n[0], n[1], n[2], n[3], n[4], n[5]);
          ele->sFlag = sflag;
          ele->soundSpeed = ss;
          break;
        case 10:
          ele = new IsoParamQuadSommer(nNodes,n);
          ele->sFlag = sflag;
          ele->soundSpeed = ss;
          break;
        case 11:
          ele = new IsoParamTriSommer(nNodes,n);
          ele->sFlag = sflag;
          ele->soundSpeed = ss;
          break;
        default:
          fprintf(stderr,"  Sower.h: etype = %d not implemented (yet ?)\n",etype);
          return;
      }
      obj->push_back(ele);
      delete[] n;
    }
};

class MatIO
{
 public:
  static size_t write(std::pair<int, SPropContainer* >* sProps, int index,BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing material index = " << index << " curObjID = " << curObjID << std::endl;
#endif
      if(index < sProps->first) {
	StructProp sp = (*(sProps->second))[index];
	double d[28+10];
	d[0] = sp.E;
	d[1] = sp.A;
	d[2] = sp.nu;
	d[3] = sp.rho;
	d[4] = sp.eh;
	d[5] = sp.Ixx;
	d[6] = sp.Iyy;
	d[7] = sp.Izz;
	d[8] = sp.c;
	d[9] = sp.k;
	d[10] = sp.Q;
	d[11] = sp.W;
	d[12] = sp.P;
	d[13] = sp.Ta;
	d[14] = sp.ymin;
	d[15] = sp.ymax;
	d[16] = sp.zmin;
	d[17] = sp.zmax;
	d[18] = sp.kappaHelm;
        // start new
        d[19] = double(sp.F_op);
        d[20] = sp.F_Uc;
        d[21] = sp.F_Uf;
        d[22] = sp.F_h;
        d[23] = sp.F_d;
        d[24] = sp.F_dlambda;
        d[25] = double(sp.F_np);
        d[26] = double(sp.F_Nf);
        d[27] = double(sp.Seed);
        // end new
        d[28] = double(sp.fp.PMLtype);
        d[29] = sp.fp.gamma;
        d[30] = sp.fp.Rx;
        d[31] = sp.fp.Ry;
        d[32] = sp.fp.Rz;
        d[33] = sp.fp.Sx;
        d[34] = sp.fp.Sy;
        d[35] = sp.fp.Sz;
        d[36] = real(sp.soundSpeed);
        d[37] = imag(sp.soundSpeed);
 
	file.write(&curObjID, 1);
	file.write(d, 28+10);
      }
      else {
        std::cerr << "out of range" << std::endl;
      }
      return 0;
    }

  typedef std::pair<int,std::map<int,StructProp>* >* oType;
  const static TypeTag dataType = MATERIALS_TYPE;

  static oType create(int size)
    {
      StructProp h;
      return(new std::pair<int, std::map<int,StructProp>* >(0, new std::map<int,StructProp>() ));
    }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
std::cerr << "Sower.h, readData, MatIO" << std::endl;
#endif
      StructProp sp;
      double d[28+10];
      s->read(d, 28+10, file);
      
      obj->first++;

      sp.E = d[0];
      sp.A = d[1];
      sp.nu = d[2] ;
      sp.rho = d[3];
      sp.eh = d[4];
      sp.Ixx = d[5];
      sp.Iyy = d[6];
      sp.Izz = d[7];
      sp.c = d[8];
      sp.k = d[9];
      sp.Q = d[10];
      sp.W = d[11];
      sp.P = d[12];
      sp.Ta = d[13];
      sp.ymin = d[14];
      sp.ymax = d[15];
      sp.zmin = d[16];
      sp.zmax = d[17];
      sp.kappaHelm = d[18];
      // start new
      sp.F_op = int(d[19]);
      sp.F_Uc = d[20];
      sp.F_Uf = d[21];
      sp.F_h = d[22];
      sp.F_d = d[23];
      sp.F_dlambda = d[24];
      sp.F_np = int(d[25]);
      sp.F_Nf = int(d[26]);
      sp.Seed = int(d[27]);
      // end new
      sp.fp.PMLtype = int(d[28]);
      sp.fp.gamma = d[29];
      sp.fp.Rx = d[30];
      sp.fp.Ry = d[31];
      sp.fp.Rz = d[32];
      sp.fp.Sx = d[33];
      sp.fp.Sy = d[34];
      sp.fp.Sz = d[35];
      sp.soundSpeed = complex<double>(d[36],d[37]);
      (*(obj->second))[localIndex]=sp; 
    }
};

template<TypeTag DATATYPE>
class MFTTDataIO
{
 public:
  static size_t write( std::pair<int, ResizeArray<MFTTData*>* >* sProps, int index,BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing MFTT index = " << index << " curObjID = " << curObjID << std::endl;
#endif
      MFTTData* sp = 0;
      for(int i = 0;i<sProps->first;++i)
	if((*(sProps->second))[i]->getID() == index)
	  {
	    sp = (*(sProps->second))[i];
	    break;
	  }
      if(sp!=0) {
	
	file.write(&curObjID, 1);
	int jav[2];
	jav[0] = sp->getID();
	jav[1] = sp->getNumPoints();

	file.write(jav,2);
	
	file.write((sp->time) ,sp->np);
	file.write((sp->value),sp->np);
      }
      else {
        std::cerr << "out of range for MFTTDataIO " << index << " " << sProps->first << std::endl;
      }
      return 0;
    }

  typedef std::pair<int, ResizeArray<MFTTData*>* >* oType;
  const static TypeTag dataType = DATATYPE;

  static oType create(int size)
    {
      return(new std::pair<int, ResizeArray<MFTTData*>* >(0, new ResizeArray<MFTTData*>(0,size) ));
    }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
std::cerr << "Sower.h, readData, MFTTDataIO" << std::endl;
#endif
      obj->first++;

      int jav[2];
      s->read(jav,2,file); // id, np
      
      double* ntime, *value;
      ntime = new double[jav[1]];
      value = new double[jav[1]];
      s->read(ntime,jav[1],file);
      s->read(value,jav[1],file);
      MFTTData* li = new MFTTData(jav[0]);
      for(int i = 0;i < jav[1];++i)
	li->add(ntime[i],value[i]);

      (*(obj->second))[localIndex]=li; 
      delete [] ntime;
      delete [] value;
    }
};


class LMPCIO
{
 public:
  static size_t write(std::pair<int, ResizeArray<LMPCons*>* >* mpc, int index, BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing lmpc index = " << index << " curObjID = " << curObjID << std::endl;
#endif
      if(index<mpc->first) {
	LMPCons* m = (*(mpc->second))[index];
	
	file.write(&curObjID, 1);
	int jav[4];
	jav[0] = m->lmpcnum;
	jav[1] = m->nterms;
        jav[2] = int(m->isComplex);
        jav[3] = m->type;
	file.write(jav,4);
	
        if(m->isComplex) {
	  DComplex ahce = m->getRhs<DComplex>();
	  double coeff[2];
	  coeff[0] = real(ahce);
	  coeff[1] = imag(ahce);
	  file.write(coeff,2);
        }
        else {
          double coeff = m->getRhs<double>();
          file.write(&coeff,1);
        }

	for(int j = 0;j<m->nterms; ++j)
	  {
	    int jav2[3];
	    jav2[0] = m->terms[j].nnum;
	    jav2[1] = m->terms[j].dofnum;
            jav2[2] = int(m->terms[j].isComplex);
	    file.write(jav2,3);
	    if(m->terms[j].isComplex) {
	      double coeff2[2];
	      coeff2[0] = real(m->terms[j].coef.c_value);
	      coeff2[1] = imag(m->terms[j].coef.c_value);
	      file.write(coeff2,2);
            }
            else {
              double coeff2 = m->terms[j].coef.r_value;
              file.write(&coeff2,1);
            }
	  }
      }
      else {
        std::cerr << "out of range for LPMCIO " << index << " " << mpc->first << std::endl;
      }
      return 0;
    }

  typedef std::pair<int, ResizeArray<LMPCons*>* >* oType;
  const static TypeTag dataType = LMPC_TYPE;

  static oType create(int size)
    {
      return(new std::pair<int, ResizeArray<LMPCons*>* >(0, new ResizeArray<LMPCons*>(0,size) ));
    }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
      std::cerr << "Sower.h, readData, LMPCIO" << std::endl;
#endif

      int jav[4];
      s->read(jav,4,file); // lmpcnum,nterms,isComplex,type

      LMPCons *lmp;
      if(jav[2]) { // isComplex
        double coeff[2];     // rhs
        s->read(coeff,2,file);
        lmp = new LMPCons(jav[0],coeff[0],coeff[1]);
      }
      else {
        double coeff; // rhs
        s->read(&coeff,1,file);
        lmp = new LMPCons(jav[0],coeff);
      }
      lmp->type = jav[3];

      obj->first++;

      for(int j = 0; j < jav[1]; ++j)
	{ // read term 
	  int nnum;
	  int dof;
          int isComplex;
	  s->readNum(NODES_TYPE, &nnum, 1, file);
	  s->read(&dof,1,file); 
          s->read(&isComplex,1,file);
          if(isComplex) {
	    double coeff2[2];     // c coeff
	    s->read(coeff2,2,file);
	    lmp->addterm(new LMPCTerm(nnum,dof,coeff2[0],coeff2[1]));
          }
          else {
            double coeff2; // r coeff
	    s->read(&coeff2,1,file);
            lmp->addterm(new LMPCTerm(nnum,dof,coeff2));
          }
        }

      (*(obj->second))[localIndex]=lmp; 
    }
};


class CompositeLIO
{
 public:
  static size_t write(std::pair<int, ResizeArray<LayInfo*>* >* sProps, int index,BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing compositeL index = " << index << " curObjID = " << curObjID << std::endl;
#endif
      if(index < sProps->first) {
	LayInfo* sp = (*(sProps->second))[index];
	
	file.write(&curObjID, 1);
	int jav[3];
	jav[0] = sp->type;
	jav[1] = sp->numLayers;
	jav[2] = sp->maxLayer;

	file.write(jav,3);
	
	//file.write(*(sp->data),9);
	//file.write(*(sp->grad),9);
      }
      else {
        std::cerr << "out of range" << std::endl;
      }
      return 0;
    }

  typedef std::pair<int,std::map<int,LayInfo*>* >* oType;
  const static TypeTag dataType = COMPOSITEL_TYPE;

  static oType create(int size)
    {
      return(new std::pair<int, std::map<int,LayInfo*>* >(0, new std::map<int,LayInfo*>() ));
    }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
std::cerr << "Sower.h, readData, CompositeLIO" << std::endl;
#endif
      obj->first++;

      int jav[3];
      s->read(jav,3,file);
      
      LayInfo* li = new LayInfo(jav[0]);
      li->numLayers = jav[1];
      li->maxLayer = jav[2];
      
      //s->read((li->data),9,file);
      //s->read((li->grad),9,file);
      (*(obj->second))[localIndex]=li; 
    }
};

class CompositeCIO
{
 public:
  static size_t write(std::pair<int, ResizeArray<CoefData*>* >* sProps, int index,BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing compositeC index = " << index << " curObjID = " << curObjID << std::endl;
#endif
      if(index < sProps->first) {
	CoefData* sp = (*(sProps->second))[index];

	file.write(&curObjID, 1);
	double c[36];
	for(int i = 0; i<6 ; ++i)
	  for(int j = 0; j<6 ; ++j)
	    {
	      c[i+j*6] = sp->c[i][j];
	    }
	file.write(c,36);
      }
      else {
        std::cerr << "out of range" << std::endl;
      }
      return 0;
    }

  typedef std::pair<int,std::map<int,CoefData*>* >* oType;
  const static TypeTag dataType = COMPOSITEC_TYPE;

  static oType create(int size)
    {
      return(new std::pair<int, std::map<int,CoefData*>* >(0, new std::map<int,CoefData*>() ));
    }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
std::cerr << "Sower.h, readData, CompositeCIO" << std::endl;
#endif
      obj->first++;
      CoefData* li = new CoefData();
      double c[36];
      s->read(c,36,file);
      for(int i = 0; i<6 ; ++i)
	for(int j = 0; j<6 ; ++j)
	  {
	    li->c[i][j] =  c[i+j*6];
	  }
      (*(obj->second))[localIndex]=li; 
    }
};

class CFramesIO
{
 public:
  static size_t write(std::pair<int, ResizeArray<double*>* >* sProps, int index,BinFileHandler& file, int curObjID)
    {
#ifdef SOWER_DEBUG
      std::cerr << "writing composite frame index = " << index << " curObjID = " << curObjID << std::endl;
#endif
      if(index < sProps->first) {
	double* sp = (*(sProps->second))[index];
	
	file.write(&curObjID, 1);
	file.write(sp,9);
      }
      else {
        std::cerr << "out of range" << std::endl;
      }
      return 0;
    }

  typedef std::pair<int,std::map<int,double*>* >* oType;
  const static TypeTag dataType = CFRAMES_TYPE;

  static oType create(int size)
    {
      return(new std::pair<int, std::map<int,double*>* >(0, new std::map<int,double*>() ));
    }

  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG
std::cerr << "Sower.h, readData, CompositeFramesIO" << std::endl;
#endif
      obj->first++;
      double* li = new double[9];
      s->read(li,9,file);
      (*(obj->second))[localIndex]=li; 
    }
};

#ifdef SOWER_SURFS
class SurfaceIO 
{
  public:
  static size_t write(SurfaceEntity* surf, int index, BinFileHandler& file, int curObjID);
  typedef SurfaceEntity* oType;
  const static TypeTag dataType = SURFACES;
                                                                                                              
  static void readData(oType obj, Sower* s, int localIndex, BinFileHandler& file)
    {
#ifdef SOWER_DEBUG_SURFS
    std::cerr << "Sower.h, readData, SurfaceIO" << std::endl;
#endif
     int surfId, nNodes, nElems, elNum, ndNum, elType, nodes[64];
     s->read(&surfId, 1, file); 
     obj->SetId(surfId);
     // read face els (using global node Ids)
     s->read(&nElems, 1, file); 
     for(int iel=0; iel<nElems; iel++) {
       s->read(&elNum, 1, file);
       s->read(&elType, 1, file);
       s->read(&nNodes, 1, file);
       s->read(nodes,nNodes,file);
       obj->AddFaceElement(elNum, elType, nNodes, nodes);
     }
     // read nodes 
     s->read(&nNodes, 1, file); 
     CoordSet* NodeSet = new CoordSet(nNodes);
     double xyz[3];
     for(int inode=0; inode<nNodes; inode++) {
       s->read(&ndNum,1,file);
       s->read(xyz,3,file);
       NodeSet->nodeadd(inode,xyz);  
     }
     obj->SetPtrNodeSet(NodeSet); //transfert ownership to the SurfaceEntity object

#ifdef SOWER_DEBUG_SURFS
     std::cerr << " Surface Id "<<surfId<<", nNodes = "<<nNodes<<", nElems = "<<nElems<< std::endl;
     obj->Print();
#endif
    }
                                                                                                              
  static oType create(int size)
    {
      return(new SurfaceEntity());
    }
};
#endif

/**
   END of IO objects for datatypes
**/


template<typename IOObject, typename DataType, typename MapType>
void 
Sower::addParentToChildData(TypeTag thisType, TypeTag parentType, int ndata, 
                            DataType data, MapType parentToChild)
{
  std::map<TypeTag,DataStruct*>::iterator it = entries.find(parentType);
  if(it != entries.end()) { // found cluster to parent connectivity
    // construct cluster to child with parent to child connectivity
    entries[thisType] = new GenDataStruct<DataType,IOObject>
                        (data, (*it).second->getSubToData().transcon(*parentToChild));
    return;
  }
  std::cerr << " *** ERROR : cannot find appropriate cluster to parent connectivity" << std::endl;
}

template<typename IOObject, typename DataType, typename MapType>
void 
Sower::addChildToParentData(TypeTag thisType, TypeTag parentType, int ndata, 
                            DataType data, MapType childToParent)
{
  std::map<TypeTag,DataStruct*>::iterator it = entries.find(parentType);
    if(it != entries.end()) { // found cluster to parent connectivity
      // construct cluster to child with parent to child connectivity
      Connectivity reChildToParent = childToParent->reverse();
      //Connectivity* reChildToParent2 = reChildToParent->reverse();
      entries[thisType] = new GenDataStruct<DataType, IOObject>
                         (data, (*it).second->getSubToData().transcon(reChildToParent));
      return;
    }
    std::cerr << " *** ERROR : cannot find appropriate cluster to parent connectivity" << std::endl;
}

/*
  warning reread toc when changing file
  @return 0 if no data for this subdomain is available
*/
#include <Utils.d/resize_array.h>

template<typename IOObject>
typename 
IOObject::oType Sower::read(BinFileHandler& file, int subNum, int*& localToGlobal, bool rereadToc)
{
  if(!tocRead || (rereadToc))
    readToc(file);
  // get the rangeSet
  RangeSet * rs = getRangeSet(IOObject::dataType, subNum, file);
  if(rs == 0)  // nothing for this subdomain in this datatype !
    return(0);
#ifdef SOWER_DEBUG
  else std::cerr  << "found some data of type " << IOObject::dataType << " in subdomain " << subNum << std::endl;
#endif
  std::list<RangeEntry*>* rangeEntries = rs->getRangeEntries();
  typename IOObject::oType obj = IOObject::create(rs->size());
  std::map<int,int>* table = new std::map<int, int>(); // global to local index table
  int localIndex = 0;
  ResizeArray<int> localToGlobalMap(0);  // PJSA
  for(std::list<RangeEntry*>::iterator it = rangeEntries->begin(); it != rangeEntries->end(); ++it) 
    { // for each range
      file.seek((*it)->offset);
      int globObjID;
      while(true) {
	file.read(&globObjID, 1); // element global number   
	(*table)[globObjID] = localIndex;
//      std::cerr << "localIndex = " << localIndex << ", globObjID = " << globObjID << std::endl;
	IOObject::readData(obj, this, localIndex, file);
        localToGlobalMap[localIndex] = globObjID;  // PJSA
        localIndex++;
	if(globObjID == (*it)->endNum) break; // TODO: reconsider this due to gaps in the element numbering
      }
    }
  TypeTag offsetagul = IOObject::dataType; // linker bug work around
  globalToLocal[offsetagul]=table;
  localToGlobal = localToGlobalMap.data(false);  // PJSA
  return(obj);
}

#endif
