#include "FileNameInfo.h"

#include <Driver.d/GeoSource.h>
#include <Driver.d/Domain.h>

#include <sstream>

extern GeoSource *geoSource;
extern Domain *domain;

namespace Rom {

std::string
FileNameInfo::getPrefix(const std::string &str) {
  return str.substr(0, str.find('.'));
}

int
FileNameInfo::size(BasisId::Type type, BasisId::Level level)
{
  int ret = 0;
  if(domain->solInfo().svdPodRom) {
    if(level == 0) {
      if(type == 4)
        ret = domain->solInfo().accelPodRomFile.size();
      else if(type == 5)
        ret = domain->solInfo().velocPodRomFile.size();
      else
        ret = domain->solInfo().snapfiPodRom.size();
    }
    else if(level == 1) {
      ret = 1; // domain->solInfo().SVDoutput;
    } 
    else if(level == 2) {
      ret = domain->solInfo().robfi.size();
    }
  }
  else {
    if(level == 0) {
      if(type == 0)
        ret = domain->solInfo().statePodRomFile.size();
      if(type == 1)
        ret = 1; // domain->solInfo().residualPodRomFile;
      if(type == 2)
        ret = 1; // domain->solInfo().jacobianPodRomFile;
      if(type == 3)
        ret = 1; // domain->solInfo().forcePodRomFile;
      if(type == 4)
        ret = domain->solInfo().accelPodRomFile.size();
      if(type == 5)
        ret = domain->solInfo().velocPodRomFile.size();
      if(type == 6)
        ret = 1; // domain->solInfo().isvPodRomFile;
      if(type == 7)
        ret = domain->solInfo().dsvPodRomFile.size();
      if(type == 8)
        ret = 1;
      if(type == 9)
        ret = domain->solInfo().muvPodRomFile.size();
    }
    else if(level == 1) {
      ret = domain->solInfo().readInROBorModes.size();
    }
  }
  return ret;
}

FileNameInfo::FileNameInfo() :
  prefix_(getPrefix(geoSource->getCheckFileInfo()->checkfile))
{}

FileNameInfo::FileNameInfo(const std::string &prefix) :
  prefix_(prefix)
{}

std::string
FileNameInfo::basisFileName(const BasisId &id, int i) const {
  std::ostringstream builder;

  if(domain->solInfo().svdPodRom) {
    if(id.level() == 0) {
      if(id.type() == 4)
        builder << domain->solInfo().accelPodRomFile[i];
      else if(id.type() == 5)
        builder << domain->solInfo().velocPodRomFile[i];
      else if(id.type() == 7)
        builder << domain->solInfo().dsvPodRomFile[i];
      else if(id.type() == 9)
        builder << domain->solInfo().muvPodRomFile[i];
      else 
        builder << domain->solInfo().snapfiPodRom[i].c_str() ; 
    }
    else if(id.level() == 1) {
      builder << domain->solInfo().SVDoutput; 
    } 
    else if(id.level() == 2) {
      builder << domain->solInfo().robfi[i].c_str();
    }
  }
  else {
    if(id.level() == 0) {
      if(id.type() == 0) 
        builder << domain->solInfo().statePodRomFile[i];
      if(id.type() == 1) 
        builder << domain->solInfo().residualPodRomFile;
      if(id.type() == 2)
        builder << domain->solInfo().jacobianPodRomFile;
      if(id.type() == 3)
        builder << domain->solInfo().forcePodRomFile;
      if(id.type() == 4)
        builder << domain->solInfo().accelPodRomFile[i];
      if(id.type() == 5)
        builder << domain->solInfo().velocPodRomFile[i];
      if(id.type() == 6)
        builder << domain->solInfo().isvPodRomFile;
      if(id.type() == 7)
        builder << domain->solInfo().dsvPodRomFile[i];
      if(id.type() == 8)
        builder << domain->solInfo().constraintSnapshotFile;
      if(id.type() == 9)
        builder << domain->solInfo().muvPodRomFile[i];
    } else if(id.level() == 1) {
      if(id.type() == 7) {
        builder << domain->solInfo().readInDualROB[i];
      } else if(id.type() == 8) {
        builder << domain->solInfo().constraintPodRomFile;
      }
      
      else {
        builder << domain->solInfo().readInROBorModes[i];
      }
    }
    else if(id.level() == 2) {
      if(id.type() == 7) {
        builder << domain->solInfo().readInDualROB[i].c_str() << ".deim";
      } else {
        builder << domain->solInfo().SVDoutput;
      }
    }
  }

  std::string mystrg;
  mystrg = builder.str();

  return builder.str();
}

bool
BasisFileId::isBinary() const
{
  const char* modeFileName = name_.c_str();

  // Open file containing mode shapes and frequencies.
  FILE *f;
  if((f=fopen(modeFileName,"r"))==(FILE *) 0 ){
    std::cerr << " *** ERROR: Cannot open " << name_ << ", exiting...\n";
    exit(0);
  }
  fflush(f);

  // Read in number of modes and number of nodes
  char buf[80];
  char *str = fgets(buf, sizeof buf, f);
  bool b;
  if(strncmp("Vector", buf, 6) == 0) {
    b = false;
  }
  else {
    int numModes;
    int count = sscanf(buf, "%d", &numModes);
    b = (count != 1);
  }
  int numNodes;
  int count = fscanf(f, "%d", &numNodes);
  b = b || (count != 1);

  // If the file is not in one of the two valid ascii formats (see manual) then it is assumed to be binary
  return b;
}

} /* end namespace Rom */
