#include <Utils.d/NodeSpaceArray.h>
#include <Utils.d/pstress.h>
#include <cmath>

template<typename BaseMaterial>
NLMaterial *
BrittleFractureTB<BaseMaterial>::clone() const
{
  return new BrittleFractureTB<BaseMaterial>(*this);
}

template<typename BaseMaterial>
int
BrittleFractureTB<BaseMaterial>::getNumStates() const
{
  const int i = BaseMaterial::getNumStates();
  return i+1;
}

template<typename BaseMaterial>
void
BrittleFractureTB<BaseMaterial>::initStates(double *state)
{
  BaseMaterial::initStates(state);
  const int i = BaseMaterial::getNumStates();
  state[i] = 0;
}

template<typename BaseMaterial>
void
BrittleFractureTB<BaseMaterial>::getStress(Tensor *stress, Tensor &strain, double* state, double temp)
{
  // Note: this function is called for post-processing.
  double tol = std::numeric_limits<double>::epsilon();
  const int i = BaseMaterial::getNumStates();

  if(!state || (state[i] < Kf + tol)) {
    BaseMaterial::getStress(stress, strain, state, temp);
  }
  else {
    stress->setZero();
  }
}

template<typename BaseMaterial>
void 
BrittleFractureTB<BaseMaterial>::integrate(Tensor *_stress, Tensor *tm, Tensor &en, Tensor &enp,
                                           double *staten, double *statenp, double temp, Tensor *cache, double dt) const
{
  double tol = std::numeric_limits<double>::epsilon();
  const int i = BaseMaterial::getNumStates();

  if(!staten || staten[i] <  Kf + tol) {
    // compute elastic response
    BaseMaterial::integrate(_stress, tm, en, enp, staten, statenp, temp, cache, dt);
    if(!staten) return;

    // check failure criteria
    Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
    std::vector<double> pvec(3);
    double svec[6];
    // Tensor_d0s2_Ss12 has xx, xy, xz, yy, yz, zz while pstress takes in xx, yy, zz, xy, yz, xz
    svec[0] = (*stress)[0];
    svec[1] = (*stress)[3];
    svec[2] = (*stress)[5];
    svec[3] = (*stress)[1];
    svec[4] = (*stress)[4];
    svec[5] = (*stress)[2];
    pstress(svec, pvec.data());
    double pmax = std::max(std::max(pvec[0],pvec[1]), pvec[2]); 
    if(pmax > maxprs) {
      statenp[i] = staten[i] + dt*pow(pmax-maxprs, exponent);
      if(statenp[i] >= Kf + tol) {
        _stress->setZero();
        tm->setZero();
      }
    }
    else
      statenp[i] = staten[i];
  }
  else {
    statenp[i] = staten[i];
    _stress->setZero();
    tm->setZero();
  }
}

template<typename BaseMaterial>
void
BrittleFractureTB<BaseMaterial>::integrate(Tensor *_stress, Tensor &en, Tensor &enp,
                                           double *staten, double *statenp, double temp, Tensor *cache, double dt) const
{
  double tol = std::numeric_limits<double>::epsilon();
  const int i = BaseMaterial::getNumStates();

  if(!staten || staten[i] <  Kf + tol) {
    // compute elastic response
    BaseMaterial::integrate(_stress, en, enp, staten, statenp, temp, cache, dt);
    if(!staten) return;

    // check failure criteria
    Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
    std::vector<double> pvec(3);
    double svec[6];
    // Tensor_d0s2_Ss12 has xx, xy, xz, yy, yz, zz while pstress takes in xx, yy, zz, xy, yz, xz
    svec[0] = (*stress)[0];
    svec[1] = (*stress)[3];
    svec[2] = (*stress)[5];
    svec[3] = (*stress)[1];
    svec[4] = (*stress)[4];
    svec[5] = (*stress)[2];
    pstress(svec, pvec.data());
    double pmax = std::max(std::max(pvec[0],pvec[1]), pvec[2]); 
    if(pmax > maxprs) {
      statenp[i] = staten[i] + dt*pow(pmax-maxprs, exponent);
      if(statenp[i] >= Kf + tol) {
        _stress->setZero();
      }
    }
    else
      statenp[i] = staten[i];
  }
  else {
    statenp[i] = staten[i];
    _stress->setZero();
  }
}

template<typename BaseMaterial>
double
BrittleFractureTB<BaseMaterial>::getDamage(double *statenp) const
{
  double tol = std::numeric_limits<double>::epsilon();
  const int i = BaseMaterial::getNumStates();

  return (statenp && statenp[i] >= Kf + tol) ? 1.0 : 0.0;
}

template<typename BaseMaterial>
void
BrittleFractureTB<BaseMaterial>::print2(std::ostream &out) const
{
  out << " TulerButcher " << maxprs << " " << exponent << " " << Kf;
}
