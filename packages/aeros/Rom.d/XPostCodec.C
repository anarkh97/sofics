#include "BasisBinaryFile.h"
#include "BasisInputFile.h"
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

  if (argc > 2) {
    print_syntax();
    return EXIT_FAILURE;
  }

  const std::string input_file_name(argv[1]);
  const std::string output_file_name = input_file_name + ".xpost";
  convert_rob<Rom::BasisBinaryInputFile, Rom::XPostOutputFile<3> >(input_file_name, output_file_name);

  return EXIT_SUCCESS;
}

void print_syntax() {
  std::printf("Syntax: xpo inputfile\n");
}

void print_help() {
  std::printf("Reduced-order-basis (rob) file converter (binary to xpost)\n");
  print_syntax();
}
