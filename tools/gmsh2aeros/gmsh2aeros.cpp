//-----------------------------------------------------------
//  Discription: This routine converts mesh files in the GMSH
//               MSH format to AERO-S nodes & topology formats.
//        Usage: <path to the MSH file (input)> <path to the AEROS mesh (output)>
//       Author: Kevin Wang (created in Jan, 2019)
//-----------------------------------------------------------

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <bits/stdc++.h>
using namespace std;

static map<int,pair<int,int> > elemType; //Maps GMSH-MSH elem-type to a pair: #nodes and AERO-S elem-type.
static set<int> volTypes;

void ProcessMSH2(ifstream &input, ofstream &output);
void ProcessMSH4(ifstream &input, ofstream &output, string outputfilename);

int main(int argc, char* argv[])
{
  elemType[2] = pair<int,int>(3,3); //3-node triangle (SurfaceTOPO)
  elemType[3] = pair<int,int>(4,1); //4-node quadrangle (SurfaceTOPO)
  elemType[4] = pair<int,int>(4,23); volTypes.insert(4); //4-node tetrahedron
  elemType[5] = pair<int,int>(8,17); volTypes.insert(5); //8-node hexahedron
  elemType[6] = pair<int,int>(6,24); volTypes.insert(6); //6-node prism

  if(argc!=3) {
    fprintf(stderr,"input format: [binary] <path to MSH file (input)> <path to AERO-S mesh file (output)>\n");
    exit(-1);
  }

  ios_base::sync_with_stdio(false);
  cin.tie(NULL);

  ifstream input(argv[1], ios::in);
  ofstream output(argv[2], ios::out);
  string outputfilename(argv[2]);
  string text;

  // Read the header
  getline(input,text);
  double msh_version = 0;
  input >> msh_version;
  if (msh_version > 2.0 - 1.0e-8 && msh_version < 3.0 - 1.0e-8) {// MSH version 2
    cout << "Starting the MSH2 processor..." << endl;
    ProcessMSH2(input, output);
  } else if (msh_version > 4.0 - 1.0e-8 && msh_version < 5.0 - 1.0e-8) {// MSH version 4
    cout << "Starting the MSH4 processor..." << endl;
    ProcessMSH4(input, output, outputfilename);
  } else 
    cout << "WARNING: Input file has MSH format " << msh_version << ". The translater only supports Format 2 and Format 4" << endl;

  input.close();
  output.close();

  return 0;
}

void ProcessMSH4(ifstream &input, ofstream &output, string outputfilename)
{
  string text;
  getline(input,text);
  getline(input,text);

  // Read Entities
  int nPoints, nCurves, nSurfaces, nVolumes;
  multimap<int,int> entity2surface; //this is a multimap because one entity can belong to multiple surfaces
  multimap<int,int> entity2volume; //this is a multimap because one entity can belong to multiple volumes
  map<int,ofstream*> surfaceTopo;
  map<int,int> surfaceTopoCounter;
  getline(input,text);
  if(text.compare("$Entities")) {
    cerr << "ERROR: Expecting '$Entities', but got " << text << " instead. Aborting." << endl;
    exit(-1);
  }

  input >> nPoints >> nCurves >> nSurfaces >> nVolumes; getline(input,text); //read end-of-line

  for(int i=0; i<nPoints+nCurves; i++)
    getline(input,text); //read the "points" & "curves" --- we don't care about them (at least at the moment).

  double x0, x1, y0, y1, z0, z1;
  int id, nPhyTag, phyTag;
  for(int i=0; i<nSurfaces; i++) {
    input >> id >> x0 >> x1 >> y0 >> y1 >> z0 >> z1 >> nPhyTag;
    if(nPhyTag > 1)
      cout <<"Note: Detected multiple (" << nPhyTag << ") physical tags for surface entity " << id << ".\n";
    for(int j=0; j<nPhyTag; j++) {
      input >> phyTag;
      entity2surface.insert(pair<int,int>(id, phyTag));
      auto it = surfaceTopo.find(phyTag);
      if(it == surfaceTopo.end()) {
        surfaceTopo[phyTag] = new ofstream(outputfilename+".surf"+to_string(phyTag),ios::out);
        surfaceTopoCounter[phyTag] = 0;
        (*surfaceTopo[phyTag]) << "SURFACETOPO " << phyTag << "\n";
      }
    }
    getline(input,text);
  }

  // now, read volume entities
  map<int, vector<int> > volumeTopo; //stores the ids of elements in each (physical) volume group
  for(int i=0; i<nVolumes; i++) {
    input >> id >> x0 >> x1 >> y0 >> y1 >> z0 >> z1 >> nPhyTag;
    if(nPhyTag > 1)
      cout <<"Note: Detected multiple (" << nPhyTag << ") physical tags for volume entity " << id << ".\n";
    for(int j=0; j<nPhyTag; j++) {
      input >> phyTag;
      entity2volume.insert(pair<int,int>(id, phyTag));
      if(volumeTopo.find(phyTag) == volumeTopo.end())
        volumeTopo[phyTag] = {};
    }
    getline(input,text); 
  }

  getline(input,text);
  if(text.compare("$EndEntities")) {
    cerr << "ERROR: Expecting '$EndEntities', but got " << text << " instead. Aborting." << endl;
    exit(-1);
  }

  // Read nodes
  do getline(input,text);
  while (text.compare("$Nodes"));
  int nBlocks, nNodes, firstNode, lastNode;
  input >> nBlocks >> nNodes >> firstNode >> lastNode;
  double *xyz = new double [lastNode*3];
  bool    *used = new bool [lastNode];
  for(int i=0; i<lastNode; i++)
    used[i] = false;

  cout << "Number of nodes: " << nNodes << "." << endl;
  if(nNodes != lastNode)
    cerr << "Warning: Node numbers have gaps. Re-numbering nodes." << endl;

  output << "NODES\n";

  int entityDim, entityTag, parametric, nNodesInBlock;
 
  for(int iBlock=0; iBlock<nBlocks; iBlock++) {
    input >> entityDim >> entityTag >> parametric >> nNodesInBlock;
    int *loc2glob = new int[nNodesInBlock];
    for(int i=0; i<nNodesInBlock; i++) {
      input >> loc2glob[i];
      used[loc2glob[i]-1] = true;
    }
    for(int i=0; i<nNodesInBlock; i++) {
      input >> xyz[(loc2glob[i]-1)*3] >> xyz[(loc2glob[i]-1)*3+1] >> xyz[(loc2glob[i]-1)*3+2];
    }
    delete[] loc2glob;
  }

  getline(input,text);
  getline(input,text);
  if(text.compare("$EndNodes")) {
    cerr << "ERROR: Expecting '$EndNodes', but got " << text << " instead. Aborting.\n";
    exit(-1);
  }

  //if node numbers have gaps, we need to fix it.
  int *indexmap = new int[lastNode];
  int currentIndex = 0;
  for(int i=0; i<lastNode; i++)
    if(used[i])
      indexmap[i] = ++currentIndex;


  for(int i=0; i<lastNode; i++) {
    if(!used[i])
      continue;
    output.width(10);
    output << indexmap[i];
    output << "  " << scientific << setprecision(10) << xyz[i*3];
    output << "  " << scientific << setprecision(10) << xyz[i*3+1];
    output << "  " << scientific << setprecision(10) << xyz[i*3+2] << "\n";
  }
  delete[] xyz;

  // read and write elements (2D & 3D)  
  getline(input,text);
  if(text.compare("$Elements")) {
    cerr << "ERROR: Expecting '$Elements', but got " << text << " instead. Aborting.\n";
    exit(-1);
  }
  int nElems, firstElem, lastElem;
  input >> nBlocks >> nElems >> firstElem >> lastElem;
  if(nElems != lastElem || firstElem != 1) {
    cerr << "ERROR: Some issue with element numbering. Aborting.\n";
    exit(-1);
  }
  int currentElemType, nElemsInBlock, id0;

  int volumeCounter = 0;
  if (nVolumes > 0) 
    output << "TOPOLOGY\n";

  for(int iBlock=0; iBlock<nBlocks; iBlock++) {
    input >> entityDim >> entityTag >> currentElemType >> nElemsInBlock;
    if (entityDim == 2) { //found a 2d entity
      if (elemType.find(currentElemType) == elemType.end()) {
        cerr << "ERROR: Unknown/unsupported element type " << currentElemType << ". Aborting.\n";
        exit(-1); 
      }
      int nNodesInElem = elemType[currentElemType].first; //actually not used
      int aerosElemTag = elemType[currentElemType].second;
      if(entity2surface.find(entityTag) == entity2surface.end()) {
        cerr << "ERROR: Found elements that do not belong to a physical surface (entityTag = " << entityTag << "). Aborting.\n";
        exit(-1);
      }
      auto ret = entity2surface.equal_range(entityTag);
      for (int i=0; i<nElemsInBlock; i++) {
        input >> id0; /*not used*/
        getline(input,text);
        for(auto it=ret.first; it!=ret.second; it++) {
          phyTag = it->second;
          (*surfaceTopo[phyTag]) << std::setw(10) << ++surfaceTopoCounter[phyTag] << std::setw(6) << aerosElemTag;
        }
        istringstream is(text);
        int thisnode;
        while(is>>thisnode) {
          for(auto it=ret.first; it!=ret.second; it++) {
            phyTag = it->second;
            (*surfaceTopo[phyTag]) << "  " << std::setw(12) << indexmap[thisnode-1];
          }
        }
        for(auto it=ret.first; it!=ret.second; it++) {
          phyTag = it->second;
          (*surfaceTopo[phyTag]) << "\n";
        }
      }
    } else if (entityDim == 3) {
      if (elemType.find(currentElemType) == elemType.end()) {
        cerr << "ERROR: Unknown/unsupported element type " << currentElemType << ". Aborting.\n";
        exit(-1); 
      }
      if(entity2volume.find(entityTag) == entity2volume.end()) {
        cerr << "ERROR: Found elements that do not belong to a physical volume (entityTag = " << entityTag << "). Aborting.\n";
        exit(-1);
      }
      auto ret = entity2volume.equal_range(entityTag);
      for(auto it=ret.first; it!=ret.second; it++) {
        int phyTag = it->second;
        volumeTopo[phyTag].reserve(volumeTopo[phyTag].size() + nElemsInBlock);
      }
      for (int i=0; i<nElemsInBlock; i++) {
        input >> id0; /*not used*/
        getline(input,text);
        output << std::setw(10) << ++volumeCounter << std::setw(6) << elemType[currentElemType].second;
        istringstream is(text);
        int thisnode;
        while(is>>thisnode)
          output << "  " << std::setw(12) << indexmap[thisnode-1];
        output << "\n";

        for(auto it=ret.first; it!=ret.second; it++) {
          phyTag = it->second;
          volumeTopo[phyTag].push_back(volumeCounter);
        }
      }
    } else {
      cerr << "ERROR: Found EntityDim = " << entityDim << ". Currently the postprocessor only supports several types of 2D & 3D elements. Aborting\n";
      exit(-1);
    }
  }

  getline(input,text);
  if(text.compare("$EndElements")) {
    cerr << "ERROR: Expecting '$EndElements', but got " << text << " instead. Aborting.\n";
    exit(-1);
  }

  //Stop reading.
  cout << "Number of 3D elements: " << volumeCounter << " (in " << volumeTopo.size() << " physical volumes).\n";
  if(volumeTopo.size()>1) {
    for(auto it = volumeTopo.begin(); it != volumeTopo.end(); it++) {
      cout << "- Physical volume " << it->first << ": ";
      vector<int>& elems(it->second);
      sort(elems.begin(), elems.end());
      for(int i=0; i<elems.size(); i++) {
        if(i==0 || elems[i-1] != elems[i]-1)
          cout << "[" << elems[i] << ", ";
        if(i==elems.size()-1 || elems[i+1] != elems[i]+1)
          cout << elems[i] << "]    ";
      }
      cout << endl;
    }
  }

  for(auto it = surfaceTopoCounter.begin(); it != surfaceTopoCounter.end(); it++) {
    cout << "Number of 2D elements in SURFACETOPO " << it->first << ": " << it->second << endl;
    output << "*" << endl;
    output << "INCLUDE \"" << outputfilename << ".surf" << it->first << "\"" << endl;
  }
  for(auto it = surfaceTopo.begin(); it != surfaceTopo.end(); it++)
    it->second->close();

  cout << "Successful Completion!" << endl;
  return;
}

void ProcessMSH2(ifstream &input, ofstream &output)
{
  string text;
  getline(input,text);
  getline(input,text);
  
  // Read (and write) the nodes (Ignore "physical names")
  do getline(input,text);
  while (text.compare("$Nodes"));
  int nNodes = 0;
  input >> nNodes; getline(input,text); //read end-of-line
  cout << "Number of nodes: " << nNodes << "." << endl;
  output << "NODES" << "\n";
  for (int i=0; i<nNodes; i++) {
    getline(input,text);
    output << text << "\n";
  }
  getline(input,text);
  if(text.compare("$EndNodes")) {
    cerr << "ERROR: Expecting '$EndNodes', but got " << text << " instead. Aborting." << endl;
    exit(-1);
  }

  // Read the elements and write to the new file
  getline(input,text);
  if(text.compare("$Elements")) {
    cerr << "ERROR: Expecting '$Elements', but got " << text << " instead. Aborting." << endl;  
    exit(-1);
  }
  int nTotElems = 0;
  input >> nTotElems; getline(input,text); //read end-of-line
  cout << "Number of elements (of all types): " << nTotElems << "." << endl;

  map<int,pair<string,int> > groups; //Maps GMSH-MSH physical group to AERO-S topology strings, and number of elems
  groups[0] = pair<string,int> (string("TOPOLOGY"),0);  //added - for the volume mesh
  groups[1] = pair<string,int> (string("SURFACETOPO 1"),0);
  groups[2] = pair<string,int> (string("SURFACETOPO 2"),0);
  groups[3] = pair<string,int> (string("SURFACETOPO 3"),0);
  groups[4] = pair<string,int> (string("SURFACETOPO 4"),0);
  groups[5] = pair<string,int> (string("SURFACETOPO 5"),0);
  groups[6] = pair<string,int> (string("SURFACETOPO 6"),0);
  groups[7] = pair<string,int> (string("SURFACETOPO 7"),0);
  groups[8] = pair<string,int> (string("SURFACETOPO 8"),0);
  groups[9] = pair<string,int> (string("SURFACETOPO 9"),0);
  groups[10] = pair<string,int> (string("SURFACETOPO 10"),0);

  int ind, type, nTag, group, tag2;
  int nGroups = groups.size();
  int currentGroup = -1;
  for (int i=0; i<nTotElems; i++) {
    input >> ind >> type >> nTag >> group >> tag2;
  //  cout << ind << " " << type << " " << nTag << " " << group << " " << tag2 << endl;
    if (volTypes.find(type) != volTypes.end()) //this is a volume element
      group = 0; //overwrite the original group.
    getline(input,text); 

    if (group != currentGroup) {// write title line
      if (groups.find(group) == groups.end()) {
        cerr << "ERROR: Unknown/unsupported group ID " << group << " for element " << ind << ". Aborting." << endl;
        exit(-1);
      }
      if (currentGroup != -1)
        cout << "Written " << groups[currentGroup].second << " elements in " << groups[currentGroup].first << "." << endl;
      output << groups[group].first << "\n";
      cout << "Writing " << groups[group].first << " to file." << endl;
      if (groups[group].second != 0) 
        cerr << "WARNING: The same group '" << groups[group].first << "' has appeared before. You may need to manually merge different segments." << endl;
    }

    if (elemType.find(type) == elemType.end()) {
      cerr << "ERROR: Unknown/unsupported element type " << type << " for element " 
           << ind << ". Aborting." << endl;
      exit(-1); 
    }
    output << ++groups[group].second << "  " << elemType[type].second << "  " << text << "\n";
    currentGroup = group;
  }
  cout << "Written " << groups[currentGroup].second << " elements in " << groups[currentGroup].first << "." << endl;
  
  // finish up.
  getline(input,text);
  if (text.compare("$EndElements")) {
    cerr << "ERROR: Expecting '$EndElements', but got '" << text << "' instead. Aborting." << endl;
    exit(-1);
  }
  cout << "Done with elements." << endl;
  getline(input,text);
  if(!input.eof()) 
    cout << "WARNING: The input file contains additional information that is not translated into the output file." << endl; 
  cout << "Successful Completion!" << endl;

  return;
}
