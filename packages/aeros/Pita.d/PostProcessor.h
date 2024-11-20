#ifndef PITA_POSTPROCESSOR_H
#define PITA_POSTPROCESSOR_H

#include "Fwk.h"

#include <map>

class GeoSource;

namespace Pita {

// class PostProcessorRoot
//
// - Handles the basic low-level file manipulation for Pita
// - Abstracts away GeoSource (global geoSource instance to be passed as a constructor parameter)
// 
// Constructor: Duplicate the underlying outputfiles
// Destructor: Close the output files
// Invariants: Track the open files, check that only "legal" files are opened

class PostProcessorRoot : public Fwk::PtrInterface<PostProcessorRoot> {
public:
  EXPORT_PTRINTERFACE_TYPES(PostProcessorRoot);

  enum FileStatus {
    NO_FILE = 0,
    CLOSED,
    OPEN
  };

private:
  class FileSet; // Dummy class  

public:
  class FileSetId : public Fwk::Ordinal<FileSet, int> {
  public:
    explicit FileSetId(int v = -1) :
      Fwk::Ordinal<FileSet, int>(v)
    {}
  }; 
  
  GeoSource * geoSource() const { return geoSource_; }
 
  FileStatus fileStatus(FileSetId fileSetId) const;
  int fileSetCount() const; 
  
  virtual void fileStatusIs(FileSetId fileSet, FileStatus status);

protected:
  PostProcessorRoot(GeoSource * geoSource,
                    int localFileCount,
                    const int * localFileId); // Pointer to array of size localFileCount
  
  virtual ~PostProcessorRoot();

private:
  GeoSource * geoSource_;
  std::map<FileSetId, FileStatus> fileStatus_;

  DISALLOW_COPY_AND_ASSIGN(PostProcessorRoot);
};

// class GenPostProcessor<IntegratorType>
//
// Provides a generic interface used to perform the output

template <typename IntegratorType>
class GenPostProcessor : public PostProcessorRoot {
public:
  EXPORT_PTRINTERFACE_TYPES(GenPostProcessor);

  typedef IntegratorType Integrator;

  virtual void outputNew(FileSetId fileSetId, const Integrator * integrator) = 0;

protected:
  GenPostProcessor(GeoSource * geoSource,
                   int localFileCount,
                   const int * localFileId) :
    PostProcessorRoot(geoSource, localFileCount, localFileId)
  {}
};

class DynamTimeIntegrator;
typedef GenPostProcessor<DynamTimeIntegrator> PostProcessor;

} // end namespace Pita

#endif /* PITA_POSTPROCESSOR_H */
