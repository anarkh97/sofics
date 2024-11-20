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

  if (!(argc == 3 || 6)) {
    print_syntax();
    return EXIT_FAILURE;
  }

  // open input file streams
  ifstream truth_file (argv[1]);
  ifstream comp_file (argv[2]);

  string header_buffer;
  int num_nodes, length2, truthFlag, compFlag, i1, i2;
  double time1, time2, tFinal;
  double a1, b1, c1, d1, e1, f1, a2, b2, c2, d2, e2, f2;
  double sumx, sumy, sumz, sumx2, sumy2, sumz2;
  double sumrx, sumry, sumrz, sumrx2, sumry2, sumrz2;
  double cum_normx, cum_normy, cum_normz, normalize_factorx, normalize_factory, normalize_factorz;
  double cum_normrx, cum_normry, cum_normrz, normalize_factorrx, normalize_factorry, normalize_factorrz;
  double relative_errorx, relative_errory, relative_errorz;
  double relative_errorrx, relative_errorry, relative_errorrz;
  bool getTime1 = true, getTime2 = true;
  // check to see of both files were successfully opened
  if(truth_file.is_open() && comp_file.is_open()) {

    if( argc == 3 ) {
      std::cout << "calculate error up to time: ";
      std::cin >> tFinal;
      std::cout << std::endl;

      std::cout << "node numbers in truthfile? [0=no,1=yes]: "; 
      std::cin >> truthFlag;
      std::cout << std::endl;

      std::cout << "node numbers in comparisonfile? [0=no,1=yes]: ";                                
      std::cin >> compFlag;
      std::cout << std::endl;
    } else if( argc == 6 ) {
      std::string tFinalString(argv[3]);
      std::string truthFlagString(argv[4]);
      std::string compFlagString(argv[5]);
      
      tFinal    = atof(tFinalString.c_str());
      truthFlag = atoi(truthFlagString.c_str());
      compFlag  = atoi(compFlagString.c_str());
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
    sumrx = 0; sumry = 0; sumrz = 0; sumrx2 = 0; sumry2 = 0; sumrz2 = 0;
    
    // begin Froebenius norm computation
    // first: loop over all timesteps
    int tcounter = 1;
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

      // second: loop over nodes
      for(int counter = 0; counter < num_nodes; counter++) {

        // if the timestamps are the same then read data from both files
        if(time1 == time2) {

          // read node number if necessary
          if(truthFlag) truth_file >> i1;
          // third: read in all dofs
          truth_file >> a1; truth_file >> b1; truth_file >> c1;
          truth_file >> d1; truth_file >> e1; truth_file >> f1;
          if(compFlag) comp_file >> i2;
          comp_file >> a2; comp_file >> b2; comp_file >> c2;
          comp_file >> d2; comp_file >> e2; comp_file >> f2;
          getTime1 = true;
          getTime2 = true;
	  
          sumx += pow((a1-a2),2);      
          sumy += pow((b1-b2),2);
          sumz += pow((c1-c2),2);
          sumrx += pow((d1-d2),2);
          sumry += pow((e1-e2),2);
          sumrz += pow((f1-f2),2);

          sumx2 += pow(a1,2);
          sumy2 += pow(b1,2);
          sumz2 += pow(c1,2);
          sumrx2 += pow(d1,2);
          sumry2 += pow(e1,2);
          sumrz2 += pow(f1,2);
        }
        else {

          if(time1 < time2) {
            if(truthFlag) truth_file >> i1;
            truth_file >> a1; truth_file >> b1; truth_file >> c1;
            truth_file >> d1; truth_file >> e1; truth_file >> f1;
            if(counter == 0) {
              std::cout << "\nskipping time step " << time1 << " in truthfile" << std::endl;
              getTime1 = true;
              getTime2 = false;
            }
          }

          if(time1 > time2) {
            if(compFlag) comp_file >> i2;
            comp_file >> a2; comp_file >> b2; comp_file >> c2;
            comp_file >> d2; comp_file >> e2; comp_file >> f2;
            if(counter == 0) {
              std::cout << "\nskipping time step " << time2 << " in comparisonfile" << std::endl;
              getTime1 = false;
              getTime2 = true;
            }
          }
        }
      }

      if(!truth_file)
        break;
    }

    // square root of differences
    cum_normx = pow(sumx,0.5);
    cum_normy = pow(sumy,0.5);
    cum_normz = pow(sumz,0.5);
    cum_normrx = pow(sumrx,0.5);
    cum_normry = pow(sumry,0.5);
    cum_normrz = pow(sumrz,0.5);

    // square root of absolute
    normalize_factorx = pow(sumx2,0.5);
    normalize_factory = pow(sumy2,0.5);
    normalize_factorz = pow(sumz2,0.5);
    normalize_factorrx = pow(sumrx2,0.5);
    normalize_factorry = pow(sumry2,0.5);
    normalize_factorrz = pow(sumrz2,0.5);

    relative_errorx = cum_normx/(normalize_factorx);
    relative_errory = cum_normy/(normalize_factory);
    relative_errorz = cum_normz/(normalize_factorz);
    relative_errorrx = cum_normrx/(normalize_factorrx);
    relative_errorry = cum_normry/(normalize_factorry);
    relative_errorrz = cum_normrz/(normalize_factorrz);

    cout << "\n*** absolute error: x ***  = " << cum_normx << endl;
    cout << "*** absolute error: y ***  = " << cum_normy << endl;
    cout << "*** absolute error: z ***  = " << cum_normz << endl;
    cout << "*** absolute error: x rotation ***  = " << cum_normrx << endl;
    cout << "*** absolute error: y rotation ***  = " << cum_normry << endl;
    cout << "*** absolute error: z rotation ***  = " << cum_normrz << endl;

    cout << "*** relative error: x ***  = " << relative_errorx*100 << "%" << endl;
    cout << "*** relative error: y ***  = " << relative_errory*100 << "%" << endl;
    cout << "*** relative error: z ***  = " << relative_errorz*100 << "%" << endl;
    cout << "*** relative error: x rotation ***  = " << relative_errorrx*100 << "%" << endl;
    cout << "*** relative error: y rotation ***  = " << relative_errorry*100 << "%" << endl;
    cout << "*** relative error: z rotation ***  = " << relative_errorrz*100 << "%" << endl;

  }
}

void print_syntax() {
  std::printf("Syntax: relerr6 truthfile comparisonfile\n");
}

void print_help() {
  std::printf("Relative error computation executable for OUTPUT6 files\n");
  print_syntax();;
}
