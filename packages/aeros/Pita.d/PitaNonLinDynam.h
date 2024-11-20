#ifndef PITA_PITANONLINDYNAM_H
#define PITA_PITANONLINDYNAM_H

#include <Problems.d/NonLinDynam.h>

namespace Pita {

// Extension of nonlinear dynamic problem descriptor
// Provides additional services and primitives required by PITA
class PitaNonLinDynamic : public NonLinDynamic {
public:
  typedef Vector VecType;
  typedef double ScalarType;

  virtual void preProcess(); // overriden

  // Solver information
  bool getInitialAcceleration() const;
 
  // Tangent stiffness matrix
  // Must call getStiffAndForce then reBuildKonly to update 
  const SparseMatrix * getStiffMatrix() const { return Kt; }
  void reBuildKonly();

  // Internal energy
  double internalEnergy(const GeomState * configuration) const;
  double internalEnergy(const VecType & displacement) const; // Helper function
  
  // Cancels out rotational dofs (ignored in velocity during time-integration)
  void zeroRotDofs(VecType & v) const;

  // Output (one file per MPI process)
  void openResidualFile();
  void pitaDynamOutput(int timeSliceRank, GeomState * geomState, VecType & velocity,
                       VecType & vp, double time, int step, VecType & force, VecType & aeroF, VecType & acceleration);
  void openOutputFiles(int sliceRank);
  void closeOutputFiles(); 

  // Constructor
  explicit PitaNonLinDynamic(Domain *);
 
private:
  SparseMatrix * Kt; // PITA explicitely requires the tangent stiffness matrix
};

} // end namespace Pita

#endif /* PITA_PITANONLINDYNAM_H */
