#include "PitaNonLinDynam.h"
#include "DynamStateSet.h"

#include <Control.d/ControlInterface.h>
#include <Comm.d/Communicator.h>
#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>
#include <Driver.d/ControlLawInfo.h>
#include <Math.d/DBSparseMatrix.h>
#include <Timers.d/StaticTimers.h>
#include <Utils.d/dofset.h>
#include <Utils.d/DistHelper.h>

#include <iostream>

extern Communicator* structCom;

namespace Pita { namespace Old {

PitaNonLinDynamic::PitaNonLinDynamic(Domain *d) :
  NonLinDynamic(d),
  pitaTimers("NonLinear Pita"),
  defaultPostProcessor_(*this)
{ 
  this->preProcess();
  this->K = domain->constructDBSparseMatrix<double>();
  
  this->computeTimeInfo();
 
  kiter = d->solInfo().pitaMainIterMax;
  Jratio = d->solInfo().pitaTimeGridRatio;
  numTSonCPU = d->solInfo().pitaProcessWorkloadMax;

  coarseDt = this->getDt() * Jratio;
  coarseDelta = this->getDelta() * Jratio;

  numTS = int( ceil( ( this->getTotalTime() / this->getDt() ) / Jratio) );

  totalTime = numTS * coarseDt;

  basisImprovementMethod = d->solInfo().pitaGlobalBasisImprovement - 1;
  if (d->solInfo().pitaLocalBasisImprovement) {
    basisImprovementMethod = 2;
  }

  projTol = d->solInfo().pitaProjTol;
}

int PitaNonLinDynamic::getInitState(DynamState<double> & ds)
{
  // Dummy vectors: We do not need that information for PITA
  GenVector<double> dummy_acc(this->solVecInfo(), 0.0);
  GenVector<double> dummy_vp(this->solVecInfo(), 0.0);
  return NonLinDynamic::getInitState(ds.disp(), ds.vel(), dummy_acc, dummy_vp);
}

int PitaNonLinDynamic::getInitSeed(DynamState<double> & ds, int sliceRank)
{
  domain->initDispVelocOnTimeSlice(ds.disp(), ds.vel(), sliceRank);
  if (userSupFunc) // Get disp/velo when a user supplied control function has been specified
  {
    double sliceTime = domain->solInfo().getTimeStep() * domain->solInfo().pitaTimeGridRatio * sliceRank;
    double *ctrdisp = (double *) alloca(sizeof(double)*claw->numSensor);
    double *ctrvel  = (double *) alloca(sizeof(double)*claw->numSensor);
    double *ctracc  = (double *) alloca(sizeof(double)*claw->numSensor);
    GenVector<double> dummy_acc(this->solVecInfo(), 0.0);
    extractControlData(ds.disp(), ds.vel(), dummy_acc, ctrdisp, ctrvel, ctracc);
    userSupFunc->usd_disp(sliceTime, ctrdisp, ctrvel, ctracc);
  }  
  return 0; // Default value for int aeroAlg
}

// Rebuild dynamic mass matrix and stiffness matrix (fine time-grid)
void PitaNonLinDynamic::reBuildKonly()
{
  times->rebuild -= getTime();

  K->zeroAll();

  Connectivity *allDofs = domain->getAllDOFs();
  for (int iele = 0; iele < domain->numElements(); ++iele) {
    K->add(kelArray[iele], (*allDofs)[iele]);
  }

  times->rebuild += getTime();
}

// Set rotational displacements equal to zero.
void PitaNonLinDynamic::zeroRotDofs(Vector & vec) const
{
  ::zeroRotDofs(*domain->getCDSA(), vec);
}

double PitaNonLinDynamic::energyNorm(const Vector &disp, const Vector &velo)
{
  return sqrt(energyDot(disp, velo, disp, velo)); 
}

double PitaNonLinDynamic::energyDot(const Vector &disp1, const Vector &velo1, const Vector &disp2, const Vector &velo2)
{
  Vector Kdisp(disp1.size());
  Vector Mvelo(velo1.size());
  K->mult(disp1, Kdisp);
  M->mult(velo1, Mvelo);
  return (Mvelo * velo2) + (Kdisp * disp2); 
}

void PitaNonLinDynamic::openResidualFile()
{
  if (this->res != (FILE*) 0)
    fclose(this->res);
                                                                                                                                                   
  int myCPU  = structCom->myID();
  
  std::stringstream s;
  s << "residuals." << myCPU;
  this->res = fopen(s.str().c_str(), "wt");
                                                                                                                                                   
  if (this->res == (FILE *) 0)
    filePrint(stderr, " *** ERROR: Cannot open residual file for CPU # %d\n", myCPU);
}

// No Aero
void
PitaNonLinDynamic::pitaDynamOutput(int timeSliceRank, GeomState* geomState, Vector& velocity,
                                   Vector& vp, double time, int step, Vector& force, Vector &aeroF)
{
  times->output -= getTime();
  domain->pitaPostProcessing(timeSliceRank, geomState, force, aeroF, time, step + 1, velocity.data(), vcx, allCorot);
  times->output += getTime();
}

void
PitaNonLinDynamic::openOutputFiles(int sliceRank)
{
  geoSource->openOutputFilesForPita(sliceRank);
}

void
PitaNonLinDynamic::closeOutputFiles()
{
  geoSource->closeOutputFiles();
}

void 
PitaNonLinDynamic::printNLPitaTimerFile(int CPUid)
{
  std::string fileNameString(geoSource->getCheckFileInfo()->checkfile);
  std::stringstream s;
  s << ".pitaTiming." << CPUid;
  fileNameString.append(s.str());
  std::ofstream out(fileNameString.c_str());
  if (out.fail())
  {
    fprintf(stderr, "Failed to open %s\n", fileNameString.c_str());
    return;
  }
  out.precision(5);
  out << std::fixed << pitaTimers;
  out.close();
}

PitaNonLinDynamic::PitaPostProcessor::PitaPostProcessor(PitaNonLinDynamic & probDesc) :
  probDesc_(probDesc),
  sliceRank_(-1)
{
}

PitaNonLinDynamic::PitaPostProcessor::~PitaPostProcessor()
{
  probDesc_.closeOutputFiles();
}

void
PitaNonLinDynamic::PitaPostProcessor::sliceRank(int rank)
{
  probDesc_.closeOutputFiles();
  if (rank >= 0)
    probDesc_.openOutputFiles(rank);
  sliceRank_ = rank;
}

} /* end namespace Old */ } /* end namespace Pita */
