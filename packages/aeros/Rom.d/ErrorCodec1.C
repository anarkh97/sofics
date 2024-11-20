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
  double time1=-1, time2=-1, time2b=-1, tFinal;
  double *a1=0, *a2=0, *b2=0;
  double sumx, sumx2, sumy, sumy2;
  double cum_normx, normalize_factorx;
  double relative_errorx;
  double maxabs = 0, maxrel = 0;
  int maxabs_node = -1, maxrel_node = -1;
  double maxabs_time = -1, maxrel_time = -1;
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

    // get header line and length of result vector 
    getline(truth_file, header_buffer);
    truth_file >> num_nodes;
    getline(comp_file, header_buffer);
    comp_file >> length2;

    // check to see if result vectors are same length
    if(num_nodes != length2) {
      cout << "incompatible files" << endl;
      return EXIT_FAILURE;
    }

    // initialize variables
    sumx = 0; sumx2 = 0;
    sumy = 0; sumy2 = 0;
    a1 = new double[num_nodes];
    a2 = new double[length2];
    b2 = new double[length2];
    
    // begin Froebenius norm computation
    // first: loop over all timesteps
    int tcounter = 1;
    while(true) {

      if(getTime1 = (time1 <= time2)) {
       truth_file >> time1;
       if(truth_file.eof() || time1 > tFinal) {
         std::cerr << "break #1 time1 = " << time1 << ", tFinal = " << tFinal << std::endl;
         break;
       }
       printf("\r time stamp %d = %f",tcounter,time1);
       tcounter++;
     }

      if(getTime2 = (time1 > time2)) {
        time2b = time2;
        comp_file >> time2;
        if(comp_file.eof() || time2 > tFinal) {
          std::cerr << "break #2 time2 = " << time2 << ", tFinal = " << tFinal << std::endl;
          break;
        }
      }

      // second: loop over nodes
      for(int counter = 0; counter < num_nodes; counter++) {

        // if the timestamps are the same then read data from both files
        if(time1 == time2) {

          // third: read in all dofs
          if(getTime1) truth_file >> a1[counter];
          if(getTime2) {
            b2[counter] = a2[counter];
            comp_file >> a2[counter];
          }
	  
          sumx += pow((a1[counter]-a2[counter]),2);
          sumx2 += pow(a1[counter],2);
          sumy += fabs(a1[counter]-a2[counter]);
          sumy2 += fabs(a1[counter]);

          maxabs = max(maxabs,fabs(a1[counter]-a2[counter]));
          if(maxabs == fabs(a1[counter]-a2[counter])) { maxabs_node = counter+1; maxabs_time = time1; }
          if(a1[counter] != 0) {
            maxrel = max(maxrel,fabs(a1[counter]-a2[counter])/fabs(a1[counter]));
            if(maxrel == fabs(a1[counter]-a2[counter])/fabs(a1[counter])) { maxrel_node = counter+1; maxrel_time = time1; }
          }

          time2b = time2;
        }
        else {

          if(time1 < time2) {
            if(getTime1) truth_file >> a1[counter];
            if(getTime2) {
              b2[counter] = a2[counter];
              comp_file >> b2[counter];
            }

            // interpolate: b2----c2----a2
            double c2 = b2[counter] + (time1-time2b)/(time2-time2b)*(a2[counter]-b2[counter]);

            sumx += pow((a1[counter]-c2),2);       
            sumx2 += pow(a1[counter],2);
            sumy += fabs(a1[counter]-c2);
            sumy2 += fabs(a1[counter]);

            maxabs = max(maxabs,fabs(a1[counter]-c2));
            if(maxabs == fabs(a1[counter]-c2)) { maxabs_node = counter+1; maxabs_time = time1; }
            if(a1[counter] != 0) { 
              maxrel = max(maxrel,fabs(a1[counter]-c2)/fabs(a1[counter]));
              if(maxrel == fabs(a1[counter]-c2)/fabs(a1[counter])) { maxrel_node = counter+1; maxrel_time = time1; }
            }
          }

          if(time1 > time2) {
            if(getTime1) truth_file >> a1[counter];
            if(getTime2) {
              b2[counter] = a2[counter];
              comp_file >> b2[counter];
            }
          }
        }
      }

      if(!truth_file)
        break;
    }

    // square root of differences
    cum_normx = pow(sumx,0.5);

    // square root of absolute
    normalize_factorx = pow(sumx2,0.5);

    relative_errorx = cum_normx/(normalize_factorx);

    cout << endl;
    cout << "*** absolute L2 error   = " << cum_normx << endl;
    cout << "*** relative L2 error   = " << relative_errorx*100 << "%" << endl;
    cout << "*** absolute L1 error   = " << sumy << endl;
    cout << "*** relative L1 error   = " << (sumy/sumy2)*100 << "%" << endl;
    cout << "*** absolute Linf error = " << maxabs << " at node = " << maxabs_node << ", time = " << maxabs_time << endl;
    cout << "*** relative Linf error = " << maxrel*100 << "%" << " at node = " << maxrel_node << ", time = " << maxrel_time << endl;

  }

  if(a1) delete [] a1;
  if(a2) delete [] a2;
  if(b2) delete [] b2;
}

void print_syntax() {
  std::printf("Syntax: relerr truthfile comparisonfile\n");
}

void print_help() {
  std::printf("Relative error computation executable\n");
  print_syntax();;
}
