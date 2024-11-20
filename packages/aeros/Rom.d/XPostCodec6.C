#include "BasisBinaryFile.h"
#include "BasisInputFile.h"
#include "BasisOutputFile.h"
#include "XPostInputFile.h"
#include "XPostOutputFile.h"
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
    convert_rob<Rom::XPostInputFile, Rom::BasisBinaryOutputFile>(input_file_name, output_file_name);
  } else if (mode == "d") {
    convert_rob<Rom::BasisBinaryInputFile, Rom::XPostOutputFile<6> >(input_file_name, output_file_name);
  } else {
    print_syntax();
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

void print_syntax() {
  std::printf("Syntax: xpo6 e|d inputfile\n");
}

void print_help() {
  std::printf("Reduced-order-basis (rob) file converter (binary/ascii)\n");
  print_syntax();
  std::printf("e: encode (ascii to binary), d: decode (binary to ascii)\n");
}
