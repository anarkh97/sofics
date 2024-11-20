#include "PostProcessor.h"

#include <Driver.d/GeoSource.h>

namespace Pita {

PostProcessorRoot::PostProcessorRoot(GeoSource * gs, int localFileCount, const int * localFileId) :
  geoSource_(gs),
  fileStatus_()
{
  for (int i = 0; i < localFileCount; ++i) {
    fileStatus_.insert(fileStatus_.end(), std::make_pair(FileSetId(localFileId[i]), CLOSED));
  }
  this->geoSource()->duplicateFilesForPita(localFileCount, localFileId);
}

PostProcessorRoot::~PostProcessorRoot() {
  this->geoSource()->closeOutputFiles();
}

PostProcessorRoot::FileStatus
PostProcessorRoot::fileStatus(FileSetId fileSetId) const {
  std::map<FileSetId, FileStatus>::const_iterator it = fileStatus_.find(fileSetId);
  return it != fileStatus_.end() ? it->second : NO_FILE;
}

int
PostProcessorRoot::fileSetCount() const {
  return fileStatus_.size();
}

void
PostProcessorRoot::fileStatusIs(FileSetId fileSet, FileStatus status) {
  if (status == NO_FILE) {
    throw Fwk::RangeException("In PostProcessorRoot::fileStatusIs: Invalid status");
  }

  std::map<FileSetId, FileStatus>::iterator it = fileStatus_.find(fileSet);
  if (it == fileStatus_.end()) {
    throw Fwk::RangeException("In PostProcessorRoot::fileStatusIs: Invalid fileSetId");
  }

  if (it->second != status) {
    if (status == OPEN) {
      this->geoSource()->openOutputFilesForPita(fileSet.value());
    } else {
      this->geoSource()->closeOutputFilesForPita(fileSet.value());
    }
    it->second = status;
  }
}

} // end namespace Pita
