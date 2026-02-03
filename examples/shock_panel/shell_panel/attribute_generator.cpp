#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<array>

using namespace std;

typedef array<double,3> Vec3D;

struct Element {
  int id;
  int type;
  array<int,4> nodes;

  Element(int id_, int type_, array<int,4> nodes_)
    : id(id_), type(type_), nodes(nodes_) { };
};

int main(int argc, char* argv[])
{

  //! design variables
  double tt = atof(argv[3]);
  double tm = atof(argv[4]);
  double tr = atof(argv[5]);

  /*
  cout << setw(16) << scientific << setprecision(5)
       << "----------------"
       << "tt = " << tt
       << "tm = " << tm
       << "tr = " << tr << "\n";
  */

  //! open I/O files
  ifstream mesh_file(argv[1], ios::in);
  ofstream out_file(argv[2], ios::out);

  if(!mesh_file) {
    cerr << "*** Error: Could not open mesh file \"" << string(argv[1])
         << "\".\n";
    exit(-1);
  }
  if(!out_file) {
    cerr << "*** Error: Could not open output file \"" << string(argv[2])
         << "\".\n";
    exit(-1);
  }

  //! create data containers
  vector<Vec3D> nodes;
  vector<Element> elements;

  //! read mesh file
  string line;

  getline(mesh_file, line);
  if(line.compare("NODES")!=0) {
    cerr << "*** Error: Expected NODES but got " << line << " instead.\n";
    exit(-1);
  }

  bool found_elements=false;
  while(getline(mesh_file, line)) {
    if(line.compare("TOPOLOGY")==0) {
      found_elements=true;
      break;
    }

    istringstream is(line);

    int id;
    Vec3D xyz;

    is >> id >> xyz[0] >> xyz[1] >> xyz[2];

    nodes.push_back(xyz);
  }

  if(!found_elements) {
    cerr << "*** Error: Could not find element topology in the mesh file.\n";
    exit(-1);
  }

  while(getline(mesh_file, line)) {
    if(line.compare("*") == 0) {
      break;
    }

    istringstream is(line);

    int id, eid;
    array<int,4> nids{-1};

    is >> id >> eid;

    int n_nodes = (eid==15) ? 3 : 4;

    for(int i=0; i<n_nodes; ++i) {
      is >> nids[i];
    }

    elements.push_back(Element(id,eid,nids));
  }

  //! compute element thickness and start writing element attributes
  double l2 = 25.0;
  for(int i=0; i<elements.size(); ++i) {
    const Element &ce = elements[i];

    //! compute element gauss-point. Using only one, i.e., centeroid.
    Vec3D gp{0};
    int n_nodes = (ce.type==15) ? 3 : 4;

    for(int n=0; n<n_nodes; ++n) {
      const Vec3D &xyz = nodes[ce.nodes[n]-1];
      gp[0] += xyz[0]/n_nodes;
      gp[1] += xyz[1]/n_nodes;
      gp[2] += xyz[2]/n_nodes;
    }

    //! compute element thickness
    double y_over_l2 = gp[1]/l2;
    double t = (gp[1] <= l2) ?
      (1-y_over_l2)*tr + y_over_l2*tm :
      (2-y_over_l2)*tm + (y_over_l2-1)*tt;

    //! output attribute
    out_file << fixed << setprecision(1)
             << ce.id << "  "
             << 0.0 << "  "
             << scientific << setprecision(1) << 2.20e11 << "  "
             << fixed << 0.3 << "  "
             << scientific << 7.6e-3 << "  "
             << fixed
             << 0.0 << "  "
             << 0.0 << "  "
             << scientific << setprecision(8) << t << "  "
             << fixed << setprecision(1)
             << 0.0 << "  "
             << 0.0 << "  "
             << 0.0 << "  "
             << 0.0 << "  "
             << 0.0 << "  "
             << 0.0 << "  "
             << 0.0 << endl;
  }

  mesh_file.close();
  out_file.close();

  return 0;

}
