//c headers
#include <climits>
#include <cunistd>

//stl headers
#include <set>
#include <string>
#include <vector>
#include <array>

//stream headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

//special functions
#include <iterator>
#include <algorithm> //std::max

using namespace std;

//#define no_argument 0
//#define required_argument 1
//#define optional_argument 2

const struct option command_opts[] = {
  {"help",                no_argument, nullptr, "h"},
  //provide the surface file when constraint is specific to particular surface.
  {"surface_topo",  optional_argument, nullptr, "t"},
  //provide the scalar aeros results file; like mass
  {"scalar_result", optional_argument, nullptr, "s"},
  //provide the 1D aeros results file; like stress, strain, etc.
  {"vector_result", optional_argument, nullptr, "v"},
  //provide the 3D aeros results file; like displacement
  {"matrix_result", optional_argument, nullptr, "m"},
  //provide the results file for a specific node
  {"node_result",   optional_argument, nullptr, "n"},
  //provide the dakota results file where function values will be written
  {"dakota_result", required_argument, nullptr, "d"},
  //end of options
  {nullptr,                    0, nullptr,   0}
};

int main(int argc, char* argv[])
{

  // relevant file names
  const char *topo_file = nullptr;
  const char *scalar_file = nullptr;
  const char *vector_file = nullptr;
  const char *matrix_file = nullptr;
  const char *node_file = nullptr;
  const char *dakota_file = nullptr;

  // whether calculations are surface specific or not
  bool surface_calculation = false;

  // parse the command line using gnu getopt_long
  int opt=0;
  while(opt != -1) {
    opt = getopt_long(argc, argv, "htsvmnd:", command_opts, nullptr);

    switch(opt) {
      case 't': {
        surface_calculation = true;
        if(optarg) {
          topo_file = optarg;
        } else {
          fprintf(stderr, "*** Error: Surface topology file not provided.\n");
          exit(EXIT_FAILURE);
        }
        break;
      }
      case 's': {
        if(optarg) {
          scalar_file = optarg;
        } else {
          fprintf(stderr, "*** Error: Aero-S output file with scalar value not provided.\n");
          exit(EXIT_FAILURE);
        }
        break;
      }
      case 'v': {
        if(optarg) {
          vector_file = optarg;
        } else {
          fprintf(stderr, "*** Error: Aero-S output file not provided.\n");
          exit(EXIT_FAILURE);
        }
        break;
      }
      case 'm': {
        if(optarg) {
          matrix_file = optarg;
        } else {
          fprintf(stderr, "*** Error: Aero-S output file value not provided.\n");
          exit(EXIT_FAILURE);
        }
        break;
      }
      case 'n': {
        if(optarg) {
          node_file = optarg;
        } else {
          fprintf(stderr, "*** Error: Aero-S output file with nodal history "
            "not provided.\n");
          exit(EXIT_FAILURE);
        }
        break;
      }
      case 'd': {
        if(optarg) {
          dakota_file = optarg;
        } else {
          fprintf(stderr, "*** Error: Path to dakota results file was not provided.\n");
          exit(EXIT_FAILURE);
        }
        break;
      }
      case '?': {
        //unrecognized option
        fprintf(stdout, "Usage: %s -d/--dakota_result <path-to-dakota-results-file>.\n", argv[0]);
        fprintf(stdout, "Options:\n");
        fprintf(stdout, "  -t or --surface_topo: ...\n");
        fprintf(stdout, "  -s or --scalar_result: ...\n");
        fprintf(stdout, "  -v or --vector_result: ...\n");
        fprintf(stdout, "  -m or --matrix_result: ...\n");
        fprintf(stdout, "  -n or --node_result: ...\n");
        exit(EXIT_FAILURE);
      }
      case 'h': {
        //help
        fprintf(stdout, "Usage: %s -d/--dakota_result <path-to-dakota-results-file>.\n", argv[0]);
        fprintf(stdout, "Options:\n");
        fprintf(stdout, "  -t or --surface_topo: ...\n");
        fprintf(stdout, "  -s or --scalar_result: ...\n");
        fprintf(stdout, "  -v or --vector_result: ...\n");
        fprintf(stdout, "  -m or --matrix_result: ...\n");
        fprintf(stdout, "  -n or --node_result: ...\n");
        exit(EXIT_SUCCESS);
      }
    }
  }

  // initialize dakota results file
  ofstream dakota_output(dakota_file, ios::out);
  if(!dakota_output) {
    fprintf(stderr, "*** Error: Could not open dakota results file %s.\n",
      dakota_file);
    exit(EXIT_FAILURE);
  }

  // common string to read lines
  string line;

  // read and write scalars
  if(scalar_file) {
    ifstream scalar_input(scalar_file, ios::in);

    if(!scalar_input) {
      fprintf(stderr, "*** Error: Could not open the Aero-S result file %s\n",
        scalar_file);
      exit(EXIT_FAILURE);
    }

    double scalar_data;
    getline(scalar_input, line);
    istringstream is(line);
    is >> scalar_data;

    dakota_output << "    "
                  << scientific << setprecision(6)
                  << scalar_data
                  //<< "  " << descriptor // TODO: get the dakota descriptor here
                  << endl;

    scalar_input.close();
  }

  // read and write nodal history value
  // NOTE: currently, we always write the maximum value in time.
  if(node_file) {

    //TODO
    fprintf(stderr, "*** Error: Currently, post-processing of nodal probe "
      "is not supported.\n");
    exit(EXIT_FAILURE);

    //ifstream node_input(node_file, ios::in);
   
    //if(!node_input) {
    //  fprintf(stderr, "*** Error: Could not open the Aero-S result file %s\n",
    //    node_file);
    //  exit(EXIT_FAILURE);
    //}

    
  }

  // read surface if requested
  set<int> surface_nodes;
  if(surface_calculation) {

    ifstream topo_input(topo_file, ios::in);

    if(!topo_input) {
      fprintf(stderr, "*** Error: Could not open the surface topology file %s.\n",
        topo_file);
      exit(EXIT_FAILURE);
    }

    topo_input.ignore(512, '\n'); // ignore header
    int r;
    for(r=0; r<INT_MAX; ++r) {
      getline(topo_input, line);
      istringstream iss(line);

      int id, type;
      iss >> id >> type;

      int num_nodes = (type == 3) ? 3 : 4;
      vector<int> data(num_nodes, -1);

      bool done = false;
      for(int i=0; i<num_nodes; ++i) {
        iss >> data[i];
        if(iss.fail()) {
          done = true;
          break;
        }
      }

      if(done) break;

      // only store unique node ids
      for(int i=0; i<num_nodes; ++i)
        surface_nodes.insert(data[i]-1);  
    }

    if(!surface_nodes.empty()) {
      fprintf(stderr, "*** Error: The surface topology contains either no "
        "nodes or duplicate nodes.\n");
      exit(EXIT_FAILURE);
    }
  }

  // read and write vector results
  // NOTE: Currently we always write the max value in space and time.
  if(vector_file) {
    ifstream vector_input(vector_file, ios::in);

    if(!vector_input) {
      fprintf(stderr, "*** Error: Could not open the Aero-S result file %s\n",
        vector_input);
      exit(EXIT_FAILURE);
    }

    // reading plastic strain
    int num_global_nodes = 0, num_topo_nodes = 0;
    vector_input.ignore(512, '\n'); // ignore header
    getline(vector_input, line); // number of elements/nodes
    istringstream is(line);
    is >> num_global_nodes;

    num_topo_nodes = (surface_nodes.empty()) ? num_global_nodes : surface_nodes.size();
    
    if(num_global_nodes < num_topo_nodes) {
      fprintf(stderr, "*** Error: Surface topology contains more nodes than the FE mesh "
        "used in Aero-S.\n");
      exit(EXIT_FAILURE);
    }

    bool done = false;
    double global_max = -DBL_MAX;
    while(getline(vector_input, line)) {
      //clear stream
      is.clear();

      // get time stamp
      double time;
      is.str(line);  
      is >> time;

      if(is.fail()) {
        done = true;
        break;
      }

      // start reading values
      // When a surface topology is provided, we only collect those values.
      // Else, we collect all values.
      vector<double> data(num_topo_nodes, -1);
      int idx=0;
      for(int i=0; i<num_topo_nodes; ++i) {
        getline(vector_input, line);
        if(surface_nodes.find(i) != surface_nodes.end()) {
          is.clear();
          is.str(line);
          is >> data[idx];
          ++idx;
        }
      }

      // compute the current max
      double current_max = *max_element(data.begin(), data.end());
      
      // update global max
      if(current_max >= global_max) global_max = current_max;

    }


    // write to results file
    dakota_output << "    "
                  << scientific << setprecision(6)
                  << global_max
                  //<< "  " << descriptor // TODO: get dakota descriptor here
                  << endl;
 
    vector_input.close();
  }

  // read and write matrix results
  // Note: Currently, we always record the maximum magnitude of the data
  // in both space and time.
  if(matrix_file) {
    ifstream matrix_input(vector_file, ios::in);

    if(!matrix_input) {
      fprintf(stderr, "*** Error: Could not open the Aero-S result file %s\n",
        matrix_input);
      exit(EXIT_FAILURE);
    }

    // reading plastic strain
    int num_global_nodes = 0, num_topo_nodes = 0;
    matrix_input.ignore(512, '\n'); // ignore header
    getline(matrix_input, line); // number of elements/nodes
    istringstream is(line);
    is >> num_global_nodes;

    num_topo_nodes = (surface_nodes.empty()) ? num_global_nodes : surface_nodes.size();
    
    if(num_global_nodes < num_topo_nodes) {
      fprintf(stderr, "*** Error: Surface topology contains more nodes than the FE mesh "
        "used in Aero-S.\n");
      exit(EXIT_FAILURE);
    }

    bool done = false;
    double global_max = -DBL_MAX;
    while(getline(matrix_input, line)) {
      //clear stream
      is.clear();

      // get time stamp
      double time;
      is.str(line);  
      is >> time;

      if(is.fail()) {
        done = true;
        break;
      }

      // start reading values
      // When a surface topology is provided, we only collect those values.
      // Else, we collect all values.
      vector<double> data(num_topo_nodes, -1);
      int idx=0;
      for(int i=0; i<num_topo_nodes; ++i) {
        getline(matrix_input, line);
        if(surface_nodes.find(i) != surface_nodes.end()) {
          is.clear();
          is.str(line);
          
          double sum=0; 
          for(int j=0; j<3; ++j) {
            double value;
            is >> value;
            sum += value*value;

          }

          data[idx] = sqrt(sum);
          ++idx;
        }
      }

      // compute the current max
      double current_max = *max_element(data.begin(), data.end());
      
      // update global max
      if(current_max >= global_max) global_max = current_max;

    }


    // write to results file
    dakota_output << "    "
                  << scientific << setprecision(6)
                  << global_max
                  //<< "  " << descriptor // TODO: get dakota descriptor here
                  << endl;
 
    matrix_input.close();
  }

  if(matrix_file) {

  }

  ofstream output(argv[4], ios::out);

  // find max
  auto max_ret = max_element(data.begin(), 
    data.end());
  
  double possible_max = *max_ret;

  // compare with previous max
  if(possible_max >= max_plastic_strain) {
    max_plastic_strain = possible_max;

/*
    // Debug: print the node Id for this max
    auto offset = distance(data.begin(), max_ret);
    auto node_iter = surface_nodes.begin();
    advance(node_iter, offset);
    max_id = *node_iter;
*/  
        

  return EXIT_SUCCESS;

}
