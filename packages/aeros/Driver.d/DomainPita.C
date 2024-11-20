#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>

/* Domain methods solely used for PITA method */

void Domain::initDispVelocOnTimeSlice(GenVector<double> & d_n, GenVector<double> & v_n, int sliceRank)
{
  // Default value if no condition is specified for any dof (not recommended)
  d_n = 0.0;
  v_n = 0.0;

  if (geoSource->getUserProvidedSeedCount() == 0)
  {
    // No initial seed was specified -- Should not happen
    fprintf(stderr, "Warning -- in Domain::initDispVelocOnTimeSlice(Vector&, Vector&, int) : No initial seeds specified\n");
    return;
  }
  
  int numDof = geoSource->getNumPitaIVel6();
 
  ConstrainedDSA *c_dsa = domain->getCDSA();
  int firstBCpos = sliceRank * numDof;
  int lastBCpos = firstBCpos + numDof;
  for (int i = firstBCpos; i < lastBCpos; ++i)
  {
    int dof = c_dsa->locate(geoSource->getPitaIVel6(i).nnum, 1 << geoSource->getPitaIVel6(i).dofnum);
    if (dof >= 0)
      v_n[dof] = geoSource->getPitaIVel6(i).val;
    dof = c_dsa->locate(geoSource->getPitaIDis6(i).nnum, 1 << geoSource->getPitaIDis6(i).dofnum);
    if (dof >= 0)
      d_n[dof] = geoSource->getPitaIDis6(i).val;
  }

  if (geoSource->getCheckFileInfo()->lastRestartFile)
    fprintf(stderr, "Warning -- in Domain::initDispVelocOnTimeSlice(Vector&, Vector&, int) : Restart not implemented for PITA...\n");
}

// Precondition: The files corresponding to the active time-slice have been opened
void Domain::pitaPostProcessing(int timeSliceRank, GeomState *geomState, Vector& force, Vector &aeroForce,
                                double time, int step, double* velocity, double *vcx, Corotator **allCorot,
                                double* acceleration, double *acx, GeomState *refState, Vector* reactions,
                                SparseMatrix *M, SparseMatrix *C)
{
  int numOutInfo = geoSource->getNumOutInfo();
  OutputInfo *oinfo = geoSource->getOutputInfo();
  for(int iInfo = 0; iInfo < numOutInfo; ++iInfo)
  {
    if (oinfo[iInfo].timeSliceRank == timeSliceRank)
      postProcessingImpl(iInfo, geomState, force, aeroForce, time, step, velocity, vcx, allCorot, acceleration,
                         acx, refState, reactions, M, C);
  }
}
