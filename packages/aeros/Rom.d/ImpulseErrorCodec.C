#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;

void print_help();
void print_syntax();

int main (int argc, char *argv[]) {

  if (argc == 1) {
    print_help();
    return EXIT_SUCCESS;
  }

  if (!(argc == 3 || 4)) {
    print_syntax();
    return EXIT_FAILURE;
  }

  // open input file streams
  ifstream truth_file (argv[1]);
  ifstream comp_file (argv[2]);

  string header_buffer;
  int num_nodes, length2;
  double time1, time2, tFinal;
  double a1, b1, c1, a2, b2, c2;
  double sumx, sumy, sumz, sumx2, sumy2, sumz2;
  double cum_normx, cum_normy, cum_normz, normalize_factorx, normalize_factory, normalize_factorz;
  double relative_errorx, relative_errory, relative_errorz;
  bool getTime1 = true, getTime2 = true;
  // check to see of both files were successfully opened
  if(truth_file.is_open() && comp_file.is_open()) {

    if( argc == 3 ) {
      std::cout << "calculate error up to time: ";
      std::cin >> tFinal;
      std::cout << std::endl;
    } else if (argc == 4) {
      std::string tFinalString(argv[3]);
 
      tFinal = atof(tFinalString.c_str()); 
    }

    // get header line and length of displacement vector 
    getline(truth_file, header_buffer);
    truth_file >> num_nodes;
    getline(comp_file, header_buffer);
    comp_file >> length2;

    // check to see if displacement vectors are same length
    if(num_nodes != length2) {
      cout << "incompatible files" << endl;
      return EXIT_FAILURE;
    }

    // initialize variables
    sumx = 0; sumy = 0; sumz = 0; sumx2 = 0; sumy2 = 0; sumz2 = 0;
    
    // begin Froebenius norm computation
    // first: loop over all timesteps
    int tcounter = 1;
    double told = 0.0;
    while(true) {

      if(getTime1) {
       truth_file >> time1;
       if(truth_file.eof() || time1 > tFinal) {
         break;
       }
     }

      if(getTime2) {
        comp_file >> time2;
        if(comp_file.eof() || time2 > tFinal) {
          break;
        }
      }

      printf("\r time stamp %d = %f",tcounter,time1);
      tcounter++;
      double dt = time1 - told;
      // second: loop over nodes
      for(int counter = 0; counter < num_nodes; counter++) {

        // if the timestamps are the same then read data from both files
        if(time1 == time2) {

          // third: read in all dofs
          truth_file >> a1; truth_file >> b1; truth_file >> c1;
          comp_file >> a2; comp_file >> b2; comp_file >> c2;
          getTime1 = true;
          getTime2 = true;
	  
          sumx += dt*(a1-a2);
          sumy += dt*(b1-b2);
          sumz += dt*(c1-c2);

          if(tcounter > 1) {
            sumx2 += a1*dt;
            sumy2 += b1*dt;
            sumz2 += c1*dt;
          } else {
            sumx2 += a1;
            sumy2 += b1;
            sumz2 += c1;
          }
        }
        else {

          if(time1 < time2) {
            truth_file >> a1; truth_file >> b1; truth_file >> c1;
            if(counter == 0) {
              std::cout << "\nskipping time step " << time1 << " in truthfile" << std::endl;
              getTime1 = true;
              getTime2 = false;
            }
          }

          if(time1 > time2) {
            comp_file >> a2; comp_file >> b2; comp_file >> c2;
            if(counter == 0) {
              std::cout << "\nskipping time step " << time2 << " in comparisonfile" << std::endl;
              getTime1 = false;
              getTime2 = true;
            }
          }
        }
      }

      told = time1;
      if(!truth_file)
        break;
    }

    relative_errorx = sumx/(sumx2);
    relative_errory = sumy/(sumy2);
    relative_errorz = sumz/(sumz2);

    cout << "*** absolute error: x ***  = " << sumx << endl;
    cout << "*** absolute error: y ***  = " << sumy << endl;
    cout << "*** absolute error: z ***  = " << sumz << endl;

    cout << "*** relative error: x ***  = " << relative_errorx*100 << "%" << endl;
    cout << "*** relative error: y ***  = " << relative_errory*100 << "%" << endl;
    cout << "*** relative error: z ***  = " << relative_errorz*100 << "%" << endl;

  }
}

void print_syntax() {
  std::printf("Syntax: relerr truthfile comparisonfile\n");
}

void print_help() {
  std::printf("Relative error computation executable\n");
  print_syntax();;
}
