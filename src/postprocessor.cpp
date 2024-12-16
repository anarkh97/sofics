//c headers
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <climits>
#include <unistd.h>
#include <getopt.h>

//stl headers
#include <set>
#include <string>
#include <vector>

//stream headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

//special functions
#include <algorithm> //std::max_element

using namespace std;

//#define no_argument 0
//#define required_argument 1
//#define optional_argument 2

const struct option command_opts[] = {
  {"help",                 no_argument, nullptr, 'h'},
  //provide the surface file when constraint is specific to particular surface.
  {"surface_topo",   required_argument, nullptr, 't'},
  //provide the aeros scalar value output; like mass
  {"scalar_result",  required_argument, nullptr, 's'},
  //provide the aeros results file; like stress, strain, disp etc.
  {"aeros_result",   required_argument, nullptr, 'a'},
  //provide the aeros probe results file; like stress, strain, etc.
  {"scalar_probe",   required_argument, nullptr, 'p'},
  //provide the aeros probe results for vector variables; like displacement
  {"vector_probe",   required_argument, nullptr, 'v'},
  //provide the dakota results file where function values will be written
  {"dakota_result",  required_argument, nullptr, 'd'},
  //write flag
  {"write",                no_argument, nullptr, 'w'},
  //end of options
  {nullptr,                          0, nullptr,   0}
};

void print_commandline_message()
{
  
  fprintf(stdout, "Usage: sofics -d/--dakota_result <path-to-dakota-results-file>.\n");
  fprintf(stdout, "Options:\n");
  fprintf(stdout, "  -t or --surface_topo: Path to the specific surface topology when "
    "post-processing results on a surface.\n");
  fprintf(stdout, "  -s or --scalar_result: Path to the output file containing a scalar "
    "value. (e.g., mass).\n");
  fprintf(stdout, "  -a or --aeros_result: Path to the output file containing the time "
    "history of an Aero-S variable (e.g., gdisplac, vmstress, epstrain).\n");
  fprintf(stdout, "  -p or --scalar_probe: Path to the output file containing the time "
    "history of an Aero-S scalar variable (e.g., vmstress, epstrain) sampled at a probe.\n");
  fprintf(stdout, "  -v or --vector_probe: Path to the output file containing the time "
    "history of an Aero-S vector variable (e.g., gdisplac, gvelocit) sampled at a probe.\n");
  fprintf(stdout, "  -w or --write: Overwrites the existing Dakota results file.\n");

}

int main(int argc, char* argv[]) {

  if(argc == 1) { 
    print_commandline_message();
    exit(EXIT_SUCCESS);
  }

  // relevant file names
  const char *topo_file   = nullptr;
  const char *scalar_file = nullptr;
  const char *aeros_file  = nullptr;
  const char *probe_file  = nullptr;
  const char *vector_file = nullptr;
  const char *dakota_file = nullptr;

  // whether calculations are surface specific or not
  bool surface_calculation = false;

  // by default this utility appends to the results file.
  // using this flag forces a re-write of the dakota file.
  bool write_flag = false;

  // parse the command line using gnu getopt_long
  int opt=0;
  while(opt != -1) {
    opt = getopt_long(argc, argv, "hwt:s:a:p:v:d:", command_opts, nullptr);

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
          cout << optarg << endl;
          fprintf(stderr, "*** Error: Aero-S output file with scalar value not provided.\n");
          exit(EXIT_FAILURE);
        }
        break;
      }
      case 'a': {
        if(optarg) {
          aeros_file = optarg;
        } else {
          fprintf(stderr, "*** Error: Aero-S output file not provided.\n");
          exit(EXIT_FAILURE);
        }
        break;
      }
      case 'p': {
        if(optarg) {
          probe_file = optarg;
        } else {
          fprintf(stderr, "*** Error: Aero-S probe output file not provided.\n");
          exit(EXIT_FAILURE);
        }
        break;
      }
      case 'v': {
        if(optarg) {
          vector_file = optarg;
        } else {
          fprintf(stderr, "*** Error: Aero-S probe output file not provided.\n");
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
      case 'w': {
        write_flag = true;
        break;
      }
      case '?': {
        //unrecognized option
        print_commandline_message();
        exit(EXIT_FAILURE);
      }
      case 'h': {
        //help
        print_commandline_message();
        exit(EXIT_SUCCESS);
      }
    }
  }

  // initialize dakota results file
  ofstream dakota_output;
  if(write_flag)
    dakota_output.open(dakota_file, ios::out);
  else
    dakota_output.open(dakota_file, ios::out | ios::app);
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

    if(surface_nodes.empty()) {
      fprintf(stderr, "*** Error: The surface topology contains either no "
        "nodes or duplicate nodes.\n");
      exit(EXIT_FAILURE);
    }
  }

  // read and write vector results
  // NOTE: Currently we always write the max value in space and time.
  // If result variable is a vector (e.g. displacement, velocity, etc.)
  // then we write its max magnitude value.
  if(aeros_file) {
    ifstream aeros_input(aeros_file, ios::in);

    if(!aeros_input) {
      fprintf(stderr, "*** Error: Could not open the Aero-S result file %s\n",
        aeros_file);
      exit(EXIT_FAILURE);
    }

    // identify the dimensionality of the variable 
    string word;
    getline(aeros_input, line);
    istringstream is(line);
    is >> word;
    
    bool is_scalar = false;
    if(word.compare("Scalar") == 0) is_scalar = true;
    else if(word.compare("Vector") == 0) is_scalar = false;
    else {
      fprintf(stderr, "*** Error: Did not understand %s in the Aero-S input file.\n",
        word.c_str());
      exit(EXIT_FAILURE);
    }

    int num_global_nodes = 0, num_topo_nodes = 0;
    getline(aeros_input, line); // number of elements/nodes
    is.clear();
    is.str(line);
    is >> num_global_nodes;

    num_topo_nodes = (surface_nodes.empty()) ? num_global_nodes : surface_nodes.size();

    if(num_global_nodes < num_topo_nodes) {
      fprintf(stderr, "*** Error: Surface topology contains more nodes than the FE mesh "
        "used in Aero-S.\n");
      exit(EXIT_FAILURE);
    }

    bool done = false;
    double global_max = -DBL_MAX;
    while(getline(aeros_input, line)) {
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
      for(int i=0; i<num_global_nodes; ++i) {
        getline(aeros_input, line);
        is.clear();
        is.str(line);

        double value=0;
        if(is_scalar) is >> value;
        else {
          double sum=0, temp;
          for(int j=0; j<3; ++j) {
            is >> temp;
            sum += temp*temp;
          }
          value = sqrt(sum);
        }

        if(surface_nodes.empty()) data[i] = value;
        else if(surface_nodes.find(i+1) != surface_nodes.end()) {
          data[idx] = value;
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
 
    aeros_input.close();
  }

  // read and write probe results for a scalar variable
  if(probe_file) {
    ifstream probe_input(probe_file, ios::in);

    if(!probe_input) {
      fprintf(stderr, "*** Error: Could not open the Aero-S probe file %s\n",
        probe_file);
      exit(EXIT_FAILURE);
    }

    //ignore header
    probe_input.ignore(512, '\n');
    double global_max = -DBL_MAX;
    while(getline(probe_input, line)) {
      double time, value;
      istringstream iss(line);
      iss >> time >> value;
      
      if(value >= global_max) global_max = value; 
    }

    // write to results file
    dakota_output << "    "
                  << scientific << setprecision(6)
                  << global_max
                  //<< "  " << descriptor // TODO: get dakota descriptor here
                  << endl;

    probe_input.close();
  }

  // read and write probe results for a vector variable
  if(vector_file) {
    ifstream probe_input(vector_file, ios::in);

    if(!probe_input) {
      fprintf(stderr, "*** Error: Could not open the Aero-S probe file %s\n",
        vector_file);
      exit(EXIT_FAILURE);
    }

    //ignore header
    probe_input.ignore(512, '\n');
    double global_max = -DBL_MAX;
    while(getline(probe_input, line)) {
      double time, val1, val2, val3;
      istringstream iss(line);
      iss >> time >> val1 >> val2 >> val3;
      
      double value = sqrt(val1*val1 + val2*val2 + val3*val3);

      if(value >= global_max) global_max = value;
    }

    // write to results file
    dakota_output << "    "
                  << scientific << setprecision(6)
                  << global_max
                  //<< "  " << descriptor // TODO: get dakota descriptor here
                  << endl;

    probe_input.close();
  }

  dakota_output.close();

  return EXIT_SUCCESS;

}
