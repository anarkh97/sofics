#ifndef ROM_SNAPSHOTPROJECTIONDRIVER_H
#define ROM_SNAPSHOTPROJECTIONDRIVER_H

#include "DriverInterface.h"

#include "Problems.d/DynamDescr.h"

#include <Math.d/Vector.h>
#include "VecBasis.h"
#include <vector>
#include "FileNameInfo.h"

class Domain;
class Corotator;
class GeomState; 
class StaticTimers;
class FileNameInfo;
class MeshDesc;

template <typename Scalar> class GenFullSquareMatrix; 
template <typename Scalar> class GenVector;
typedef GenVector<double> Vector;

namespace Rom {

class SnapshotProjectionDriver : public SingleDomainDynamic, public DriverInterface {
public:
  virtual void solve(); // overriden
  int vectorSize() const;

  explicit SnapshotProjectionDriver(Domain *);
  ~SnapshotProjectionDriver();

protected:
  void preProcess();
  void postProcess();
  void compProjError();

  int elementCount() const;

  VecBasis podBasis_;
  VecBasis snapshots;
  VecBasis displac_;
  VecBasis *velocSnapshots;
  VecBasis *veloc_;
  VecBasis *accelSnapshots;
  VecBasis *accel_;
  std::vector<double> timeStamps_;
};

} // end namespace Rom

Rom::DriverInterface *snapshotProjectionDriverNew(Domain *);

#endif /* ROM_SNAPSHOTPROJECTIONDRIVER_H */
