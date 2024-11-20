#include "BasisBinaryFile.h"
#include "BasisInputFile.h"
#include "BasisOutputFile.h"
#include "RobCodec.h"

#include <cstdlib>
#include <cstdio>
#include <string>

void print_help();
void print_syntax();

int main(int argc, char *argv[]) {
  if (argc == 1) {
    print_help();
    return EXIT_SUCCESS;
  }

  if (argc == 2) {
    print_syntax();
    return EXIT_FAILURE;
  }

  const std::string input_file_name(argv[2]);
  const std::string output_file_name = input_file_name + ".out";
  
  const std::string mode(argv[1]);
  if (mode == "e") {
    convert_rob<Rom::BasisInputFile, Rom::BasisBinaryOutputFile>(input_file_name, output_file_name);
  } else if (mode == "d") {
    convert_rob<Rom::BasisBinaryInputFile, Rom::BasisOutputFile>(input_file_name, output_file_name);
  } else if (mode == "c") {
    std::vector<std::string> input_file_names(argc-3);
    for(int i=0; i<argc-3; i++) {
      input_file_names[i] = std::string(argv[2+i]);
    }
    std::string output_file_name = std::string(argv[argc-1]);
    convert_rob<Rom::BasisBinaryInputFile, Rom::BasisBinaryOutputFile>(input_file_names, output_file_name);
  } else if (mode == "t") {
    std::vector<std::string> input_file_names((argc-3)/2);
    std::vector<int> numvec((argc-3)/2);
    for(int i=0; i<argc-3; i+=2) {
      input_file_names[i/2] = std::string(argv[2+i]);
      numvec[i/2] = std::atoi(argv[2+i+1]);
    }
    std::string output_file_name = std::string(argv[argc-1]);
    convert_rob<Rom::BasisBinaryInputFile, Rom::BasisBinaryOutputFile>(input_file_names, output_file_name, numvec);
  } else if (mode == "i") {
    std::vector<std::string> input_file_names(argc-4);
    for(int i=0; i<argc-4; i++) {
      input_file_names[i] = std::string(argv[2+i]);
    }
    int stripe = std::atoi(argv[argc-2]);
    std::string output_file_name = std::string(argv[argc-1]);
    convert_rob<Rom::BasisBinaryInputFile, Rom::BasisBinaryOutputFile>(input_file_names, output_file_name, stripe);
  } else {
    print_syntax();
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

void print_syntax() {
  std::printf("Syntax: rob e|d inputfile\n");
  std::printf("      : rob c inputfile1 inputfile2 ... inputfileN outputfile\n");
  std::printf("      : rob t inputfile1 numvec1 inputfile2 numvec2 ... inputfileN numvec2 outputfile\n");
  std::printf("      : rob i inputfile1 inputfile2 ... inputfileN stripe outputfile\n");
}

void print_help() {
  std::printf("Reduced-order-basis (rob) file converter (binary/ascii)\n");
  print_syntax();
  std::printf("e: encode (ascii to binary), d: decode (binary to ascii),\n"
              "c: concatenate (binary to binary), t: truncate and concatenate (binary to binary)\n"
              "i: interleave (binary to binary)\n");
}
