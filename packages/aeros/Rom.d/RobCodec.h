#ifndef ROM_ROBCODEC_H
#define ROM_ROBCODEC_H

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include "NodeDof6Buffer.h"

template <typename InFile, typename OutFile>
int transfer_rob(InFile &input, OutFile &output, int numvec = std::numeric_limits<int>::max()) {
  Rom::NodeDof6Buffer state_buffer(input.nodeCount());
  int count = 0;
  while (input.validCurrentState() && count < numvec) {
    const double header = input.currentStateHeaderValue();
    input.currentStateBuffer(state_buffer);
    output.stateAdd(state_buffer, header);
    input.currentStateIndexInc();
    count++;
  }
  return count;
}

template <typename InFile, typename OutFile>
void convert_rob(const std::string &inFilename, const std::string &outFilename) {
  InFile input(inFilename);
  OutFile output(outFilename, input.nodeIdBegin(), input.nodeIdEnd(), false);
  transfer_rob(input, output);
}

template <typename InFile, typename OutFile>
void convert_rob(const std::vector<std::string> &inFilenames, const std::string &outFilename) {
  InFile firstInput(inFilenames.front());
  OutFile output(outFilename, firstInput.nodeIdBegin(), firstInput.nodeIdEnd(), false);
  std::string s = outFilename+".sources";
  std::ofstream sources(s.c_str());
  int fileIndex = 0;
  for(std::vector<std::string>::const_iterator it = inFilenames.begin(); it != inFilenames.end(); ++it) {
    InFile input(*it);
    if(!std::equal(firstInput.nodeIdBegin(), firstInput.nodeIdEnd(), input.nodeIdBegin())) {
      std::cerr << " *** ERROR: input files are not consistent\n";
      exit(-1);
    }
    int count = transfer_rob(input, output);
    for(int i = 0; i < count; ++i) sources << fileIndex+1 << " " << i+1 << std::endl;
    fileIndex++;
  }
  sources.close();
}

template <typename InFile, typename OutFile>
void convert_rob(const std::vector<std::string> &inFilenames, const std::string &outFilename, const std::vector<int> &numvec) {
  InFile firstInput(inFilenames.front());
  OutFile output(outFilename, firstInput.nodeIdBegin(), firstInput.nodeIdEnd(), false);
  std::vector<std::string>::const_iterator it = inFilenames.begin();
  std::vector<int>::const_iterator it2 = numvec.begin();
  for(; it != inFilenames.end(); ++it, ++it2) {
    InFile input(*it);
    if(!std::equal(firstInput.nodeIdBegin(), firstInput.nodeIdEnd(), input.nodeIdBegin())) {
      std::cerr << " *** ERROR: input files are not consistent\n";
      exit(-1);
    }
    transfer_rob(input, output, *it2);
  }
}

template <typename InFile, typename OutFile>
void convert_rob(const std::vector<std::string> &inFilenames, const std::string &outFilename, int stripe) {
  std::vector<InFile*> inputs;
  for(std::vector<std::string>::const_iterator it = inFilenames.begin(); it != inFilenames.end(); ++it) {
    inputs.push_back(new InFile(*it));
    if(it != inFilenames.begin()) {
      if(!std::equal(inputs.front()->nodeIdBegin(), inputs.front()->nodeIdEnd(), inputs.back()->nodeIdBegin())) {
        std::cerr << " *** ERROR: input files are not consistent\n";
        exit(-1);
      }
    }
  }
  OutFile output(outFilename, inputs.front()->nodeIdBegin(), inputs.front()->nodeIdEnd(), false);
  while(inputs.front()->validCurrentState()) {
    for(typename std::vector<InFile*>::iterator it = inputs.begin(); it != inputs.end(); ++it) {
      transfer_rob(*(*it), output, stripe);
    }
  }
  for(typename std::vector<InFile*>::iterator it = inputs.begin(); it != inputs.end(); ++it) delete *it;
}


#endif /* ROM_ROBCODEC_H */
