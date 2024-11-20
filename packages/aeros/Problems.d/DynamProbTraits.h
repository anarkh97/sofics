#ifndef DYNAMPROBTRAITS_H
#define DYNAMPROBTRAITS_H

template <typename ProbType, typename VecType>
inline
void
handleForce(ProbType &, VecType &) {
  // Do nothing by default
}

template <typename ProbType, typename VecType>
inline
void
handleAcceleration(ProbType &, VecType &) {
  // Do nothing by default
}

template <typename ProbType, typename VecType>
inline
void
handleVelocity(ProbType &, VecType &) {
  // Do nothing by default
}

template <typename ProbType, typename VecType>
inline
void
handleDisplacement(ProbType &, VecType &) {
  // Do nothing by default
}

template <typename ProbType, typename VecType>
inline
void
handleResidual(ProbType &, VecType &) {
  // Do nothing by default
}

template <typename ProbType, typename VecType>
inline
void
handleJacobian(ProbType &, VecType &) {
  // Do nothing by default
}

#endif /* DYNAMPROBTRAITS_H */
