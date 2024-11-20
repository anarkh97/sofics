#include "PitaNonLinDynam.h"

#include <Comm.d/Communicator.h>
#include <Corotational.d/Corotator.h>
#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>
#include <Math.d/SparseMatrix.h>
#include <Math.d/DBSparseMatrix.h>
#include <Timers.d/StaticTimers.h>
#include <Utils.d/dofset.h>
#include <Utils.d/DistHelper.h>

#include <sstream>

extern Communicator* structCom;

namespace Pita {

PitaNonLinDynamic::PitaNonLinDynamic(Domain *domain) :
  NonLinDynamic(domain)
{}

void
PitaNonLinDynamic::preProcess() {
  NonLinDynamic::preProcess();
  Kt = domain->constructDBSparseMatrix<double>();
  computeTimeInfo();
}

bool
PitaNonLinDynamic::getInitialAcceleration() const {
  return domain->solInfo().iacc_switch;
}

// Rebuild dynamic mass matrix and stiffness matrix (fine time-grid)
void
PitaNonLinDynamic::reBuildKonly() {
  times->rebuild -= getTime();

  Kt->zeroAll();

  Connectivity *allDofs = domain->getAllDOFs();
  for (int iele = 0; iele < domain->numElements(); ++iele) {
    Kt->add(kelArray[iele], (*allDofs)[iele]);
  }

  times->rebuild += getTime();
}

// Set rotational displacements equal to zero.
void
PitaNonLinDynamic::zeroRotDofs(VecType & vec) const {
  ::zeroRotDofs(*domain->getCDSA(), vec);
}

double
PitaNonLinDynamic::internalEnergy(const GeomState * configuration) const {
  double result = 0.0;
  for (int iele = 0; iele < domain->numElements(); ++iele) {
    result += allCorot[iele]->getElementEnergy(const_cast<GeomState &>(*configuration), domain->getNodes());
  }
  return result;
}

double
PitaNonLinDynamic::internalEnergy(const VecType & displacement) const {
  GeomState * configuration = const_cast<PitaNonLinDynamic *>(this)->createGeomState();
  configuration->update(const_cast<VecType &>(displacement));

  double result = internalEnergy(configuration);

  delete configuration;
  return result;
}

void
PitaNonLinDynamic::openResidualFile() {
  if (res != (FILE*) 0)
    fclose(res);
  
  int myCPU  = structCom->myID();
  
  std::stringstream s;
  s << "residuals." << myCPU;
  res = fopen(s.str().c_str(), "wt");
  
  if (res == (FILE *) 0) {
    filePrint(stderr, " *** ERROR: Cannot open residual file for CPU # %d\n", myCPU);
  }
}

void
PitaNonLinDynamic::pitaDynamOutput(int timeSliceRank, GeomState* geomState, VecType & velocity,
                                   VecType & vp, double time, int step, VecType & force, VecType & aeroF,
                                   VecType & acceleration) {
  // Note: No Aero
  times->output -= getTime();
  domain->pitaPostProcessing(timeSliceRank, geomState, force, aeroF,
                             time, step, velocity.data(), vcx,
                             allCorot, acceleration.data());
  times->output += getTime();
}

void
PitaNonLinDynamic::openOutputFiles(int sliceRank) {
  geoSource->openOutputFilesForPita(sliceRank);
}

void
PitaNonLinDynamic::closeOutputFiles() {
  geoSource->closeOutputFiles();
}

} /* end namespace Pita */
